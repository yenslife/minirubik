/* bench_ida.c
 *
 * Host-side benchmark for:
 *   1. h = 0 + same-face pruning
 *   2. orientation PDB
 *   3. permutation PDB
 *   4. max(orientation PDB, permutation PDB)
 *
 * Put this file beside solver.c.
 */

#define main baseline_main
#include "../solver.c"
#undef main

#include <inttypes.h>
#include <limits.h>

/*
 * For the experiment we build all 9 HTM moves directly.
 *
 * perm_move[m][p]:
 *   permutation rank p --move m--> new permutation rank
 *
 * ori_move[m][o]:
 *   orientation rank o --move m--> new orientation rank
 */
static uint16_t bm_perm_move[MOVES][PERMUTATIONS];
static uint16_t bm_ori_move[MOVES][ORIENTATIONS];

/* Pattern databases */
static uint8_t bm_perm_pdb[PERMUTATIONS];
static uint8_t bm_ori_pdb[ORIENTATIONS];


/*
 * Build full 9-move transition tables from the baseline's
 * 3 quarter-turn transition functions.
 *
 * move numbering follows solver.c:
 *
 * 0 R
 * 1 R2
 * 2 R'
 * 3 B
 * 4 B2
 * 5 B'
 * 6 D
 * 7 D2
 * 8 D'
 */
static void bm_build_move_tables(void)
{
    uint16_t qperm[3][PERMUTATIONS];
    uint16_t qori[3][ORIENTATIONS];
    state_t state;

    /*
     * One quarter-turn transition for each permutation rank.
     */
    for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
        /*
         * Orientation = solved, because permutation transition
         * is independent of orientation.
         */
        unrank_state((uint32_t) rank * ORIENTATIONS, &state);

        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn(state, face);

            qperm[face][rank] =
                (uint16_t) (rank_state(&next) / ORIENTATIONS);
        }
    }

    /*
     * One quarter-turn transition for each orientation rank.
     */
    for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
        /*
         * Permutation = solved, because orientation transition
         * is independent of permutation.
         */
        unrank_state(rank, &state);

        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn(state, face);

            qori[face][rank] =
                (uint16_t) (rank_state(&next) % ORIENTATIONS);
        }
    }

    /*
     * Expand quarter-turn tables into all 9 HTM moves.
     *
     * turn = 0 -> face
     * turn = 1 -> face2
     * turn = 2 -> face'
     */
    for (uint8_t face = 0; face < 3; ++face) {
        for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
            uint16_t cur = rank;

            for (uint8_t turn = 0; turn < 3; ++turn) {
                cur = qperm[face][cur];

                bm_perm_move[face * 3U + turn][rank] = cur;
            }
        }

        for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
            uint16_t cur = rank;

            for (uint8_t turn = 0; turn < 3; ++turn) {
                cur = qori[face][cur];

                bm_ori_move[face * 3U + turn][rank] = cur;
            }
        }
    }
}


/*
 * Generic BFS for one abstract state space.
 *
 * dist[abstract_state] =
 *     exact shortest distance from solved abstract state.
 *
 * move_table is stored as:
 *
 *     move_table[move * nstates + state]
 */
static uint8_t bm_build_pdb(uint8_t *dist,
                            uint16_t nstates,
                            const uint16_t *move_table)
{
    uint16_t *queue =
        malloc((size_t) nstates * sizeof(*queue));

    if (!queue) {
        fputs("malloc failed\n", stderr);
        exit(1);
    }

    memset(dist, UINT8_MAX, nstates);

    uint16_t head = 0;
    uint16_t tail = 1;

    uint8_t max_distance = 0;

    /*
     * rank 0 is solved for both abstractions.
     */
    queue[0] = 0;
    dist[0] = 0;

    while (head < tail) {
        uint16_t here = queue[head++];
        uint8_t d = dist[here];

        if (d > max_distance)
            max_distance = d;

        for (uint8_t move = 0; move < MOVES; ++move) {
            uint16_t there =
                move_table[(size_t) move * nstates + here];

            if (dist[there] != UINT8_MAX)
                continue;

            /*
             * BFS:
             * first visit is the exact shortest distance.
             */
            dist[there] = (uint8_t) (d + 1U);
            queue[tail++] = there;
        }
    }

    /*
     * H2-style sanity check:
     * every abstract state must be reachable.
     */
    if (tail != nstates) {
        fprintf(stderr,
                "PDB incomplete: reached %u / %u states\n",
                tail,
                nstates);
        exit(1);
    }

    free(queue);

    return max_distance;
}


typedef enum {
    BM_H_ZERO,
    BM_H_ORIENTATION,
    BM_H_PERMUTATION,
    BM_H_MAX
} bm_heuristic_t;


static uint8_t bm_heuristic(bm_heuristic_t kind,
                            uint16_t p,
                            uint16_t o)
{
    switch (kind) {
    case BM_H_ZERO:
        return 0;

    case BM_H_ORIENTATION:
        return bm_ori_pdb[o];

    case BM_H_PERMUTATION:
        return bm_perm_pdb[p];

    case BM_H_MAX:
    default:
        return bm_ori_pdb[o] > bm_perm_pdb[p]
                   ? bm_ori_pdb[o]
                   : bm_perm_pdb[p];
    }
}


typedef struct {
    /*
     * Across all IDA* iterations.
     */
    uint64_t expanded_total;
    uint64_t generated_total;

    /*
     * Only current threshold iteration.
     */
    uint64_t expanded_iter;
    uint64_t generated_iter;

    int solution_depth;

    /*
     * Pocket Cube HTM diameter = 11.
     */
    uint8_t path[12];
} bm_stats_t;


#define BM_FOUND (-1)


/*
 * Host-side recursive IDA*.
 *
 * Recursion is okay here because this is only the Stage-2
 * native benchmark. The final RV32I implementation can replace
 * this with an explicit fixed-size stack.
 */
static int bm_dfs(uint16_t p,
                  uint16_t o,
                  int g,
                  int bound,
                  int prev_face,
                  bm_heuristic_t heuristic_kind,
                  bm_stats_t *stats)
{
    int h = bm_heuristic(heuristic_kind, p, o);
    int f = g + h;

    /*
     * IDA* heuristic pruning.
     */
    if (f > bound)
        return f;

    /*
     * rank 0 / rank 0 = complete solved state.
     */
    if (p == 0 && o == 0) {
        stats->solution_depth = g;
        return BM_FOUND;
    }

    /*
     * We count a node as expanded only when:
     *
     * 1. it passed f <= bound
     * 2. it is not the goal
     * 3. we are about to generate children
     */
    ++stats->expanded_total;
    ++stats->expanded_iter;

    int next_bound = INT_MAX;

    for (uint8_t move = 0; move < MOVES; ++move) {
        int face = move / 3;

        /*
         * Same-face pruning:
         *
         * if previous move was R/R2/R',
         * don't try R/R2/R' again.
         */
        if (face == prev_face)
            continue;

        ++stats->generated_total;
        ++stats->generated_iter;

        uint16_t next_p = bm_perm_move[move][p];
        uint16_t next_o = bm_ori_move[move][o];

        stats->path[g] = move;

        int result =
            bm_dfs(next_p,
                   next_o,
                   g + 1,
                   bound,
                   face,
                   heuristic_kind,
                   stats);

        if (result == BM_FOUND)
            return BM_FOUND;

        /*
         * Classical IDA*:
         *
         * next threshold =
         * minimum f-value exceeding current threshold.
         */
        if (result < next_bound)
            next_bound = result;
    }

    return next_bound;
}


static void bm_run(uint16_t p,
                   uint16_t o,
                   bm_heuristic_t heuristic_kind,
                   const char *name)
{
    bm_stats_t stats = {0};

    stats.solution_depth = -1;

    /*
     * Standard IDA* starts at h(root).
     */
    int bound =
        bm_heuristic(heuristic_kind, p, o);

    printf("\n[%s]\n", name);
    printf("h(root) = %d\n", bound);

    for (;;) {
        stats.expanded_iter = 0;
        stats.generated_iter = 0;

        int next =
            bm_dfs(p,
                   o,
                   0,
                   bound,
                   -1,
                   heuristic_kind,
                   &stats);

        printf("bound=%d"
               " expanded=%" PRIu64
               " generated=%" PRIu64 "\n",
               bound,
               stats.expanded_iter,
               stats.generated_iter);

        if (next == BM_FOUND)
            break;

        if (next == INT_MAX) {
            fputs("search failed\n", stderr);
            exit(1);
        }

        bound = next;
    }

    printf("solution depth = %d\n",
           stats.solution_depth);

    printf("total expanded = %" PRIu64 "\n",
           stats.expanded_total);

    printf("total generated = %" PRIu64 "\n",
           stats.generated_total);

    printf("solution =");

    for (int i = 0;
         i < stats.solution_depth;
         ++i) {
        printf(" %s", move_names[stats.path[i]]);
    }

    putchar('\n');
}


int main(int argc, char **argv)
{
    /*
     * Default to the assignment's distance-11 example.
     */
    const char *input =
        argc == 2
            ? argv[1]
            : "21345671111111";

    state_t state;

    if (!parse_state(input, &state)) {
        fprintf(stderr,
                "usage: %s PPPPPPPOOOOOOO\n",
                argv[0]);
        return 2;
    }

    /*
     * Host-only precomputation.
     */
    bm_build_move_tables();

    uint8_t orientation_diameter =
        bm_build_pdb(
            bm_ori_pdb,
            ORIENTATIONS,
            &bm_ori_move[0][0]);

    uint8_t permutation_diameter =
        bm_build_pdb(
            bm_perm_pdb,
            PERMUTATIONS,
            &bm_perm_move[0][0]);

    /*
     * Split baseline dense rank into p/o.
     */
    uint32_t rank =
        rank_state(&state);

    uint16_t p =
        (uint16_t) (rank / ORIENTATIONS);

    uint16_t o =
        (uint16_t) (rank % ORIENTATIONS);

    printf("state = %s\n", input);
    printf("permutation rank = %u\n", p);
    printf("orientation rank = %u\n", o);

    printf("orientation PDB max distance = %u\n",
           orientation_diameter);

    printf("permutation PDB max distance = %u\n",
           permutation_diameter);

    /*
     * Compare exactly the four Stage-2 candidates.
     */
    bm_run(
        p,
        o,
        BM_H_ZERO,
        "h=0 + same-face pruning");

    bm_run(
        p,
        o,
        BM_H_ORIENTATION,
        "orientation PDB");

    bm_run(
        p,
        o,
        BM_H_PERMUTATION,
        "permutation PDB");

    bm_run(
        p,
        o,
        BM_H_MAX,
        "max(orientation, permutation)");

    return 0;
}
