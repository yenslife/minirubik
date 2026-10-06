/* solver_heuristic.c
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
#include "solver.c"
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

#include <time.h>

typedef enum {
    BM_TRANS_3,
    BM_TRANS_6,
    BM_TRANS_9
} bm_transition_design_t;

/*
 * 3-table design:
 *   R, B, D
 *
 * 6-table design:
 *   R, R2, B, B2, D, D2
 *
 * 9-table design:
 *   R, R2, R', B, B2, B', D, D2, D'
 */
static uint16_t bm_perm_move3[3][PERMUTATIONS];
static uint16_t bm_ori_move3[3][ORIENTATIONS];

static uint16_t bm_perm_move6[6][PERMUTATIONS];
static uint16_t bm_ori_move6[6][ORIENTATIONS];


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
    state_t state;

    /*
     * First build one quarter-turn table for each face:
     *
     *   0 = R
     *   1 = B
     *   2 = D
     */
    for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
        unrank_state((uint32_t) rank * ORIENTATIONS, &state);

        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn(state, face);

            bm_perm_move3[face][rank] =
                (uint16_t) (rank_state(&next) / ORIENTATIONS);
        }
    }

    for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
        unrank_state(rank, &state);

        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn(state, face);

            bm_ori_move3[face][rank] =
                (uint16_t) (rank_state(&next) % ORIENTATIONS);
        }
    }

    /*
     * Derive the 6-table and 9-table designs from the
     * quarter-turn tables.
     */
    for (uint8_t face = 0; face < 3; ++face) {
        for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
            uint16_t q1 =
                bm_perm_move3[face][rank];

            uint16_t q2 =
                bm_perm_move3[face][q1];

            uint16_t q3 =
                bm_perm_move3[face][q2];

            /*
             * 6-table:
             * face * 2 + 0 = quarter turn
             * face * 2 + 1 = half turn
             */
            bm_perm_move6[face * 2U][rank] = q1;
            bm_perm_move6[face * 2U + 1U][rank] = q2;

            /*
             * 9-table:
             * face * 3 + 0 = R/B/D
             * face * 3 + 1 = R2/B2/D2
             * face * 3 + 2 = R'/B'/D'
             */
            bm_perm_move[face * 3U][rank] = q1;
            bm_perm_move[face * 3U + 1U][rank] = q2;
            bm_perm_move[face * 3U + 2U][rank] = q3;
        }

        for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
            uint16_t q1 =
                bm_ori_move3[face][rank];

            uint16_t q2 =
                bm_ori_move3[face][q1];

            uint16_t q3 =
                bm_ori_move3[face][q2];

            bm_ori_move6[face * 2U][rank] = q1;
            bm_ori_move6[face * 2U + 1U][rank] = q2;

            bm_ori_move[face * 3U][rank] = q1;
            bm_ori_move[face * 3U + 1U][rank] = q2;
            bm_ori_move[face * 3U + 2U][rank] = q3;
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

static const char *const bm_move_symbol[MOVES] = {
    "R",
    "R2",
    "Ri",
    "B",
    "B2",
    "Bi",
    "D",
    "D2",
    "Di",
};

static void bm_emit_u16_array(FILE *out,
                              const uint16_t *data,
                              size_t count)
{
    const size_t per_line = 8;

    for (size_t i = 0; i < count; i += per_line) {
        fprintf(out, "    .half ");

        size_t end = i + per_line;

        if (end > count)
            end = count;

        for (size_t j = i; j < end; ++j) {
            if (j != i)
                fprintf(out, ", ");

            fprintf(out, "%u", data[j]);
        }

        fputc('\n', out);
    }
}

static void bm_emit_u8_array(FILE *out,
                             const uint8_t *data,
                             size_t count)
{
    const size_t per_line = 16;

    for (size_t i = 0; i < count; i += per_line) {
        fprintf(out, "    .byte ");

        size_t end = i + per_line;

        if (end > count)
            end = count;

        for (size_t j = i; j < end; ++j) {
            if (j != i)
                fprintf(out, ", ");

            fprintf(out, "%u", data[j]);
        }

        fputc('\n', out);
    }
}

/*
 * Build exact distances for the complete 3,674,160-state space.
 *
 * Host benchmark only. This table will NOT be part of the final
 * RISC-V implementation.
 */
static uint8_t *bm_build_full_distance(uint8_t *diameter,
                                       uint32_t *distance11_count)
{
    uint8_t *dist = malloc(STATES);
    uint32_t *queue =
        malloc((size_t) STATES * sizeof(*queue));

    if (!dist || !queue) {
        free(dist);
        free(queue);
        fputs("full-distance allocation failed\n", stderr);
        exit(1);
    }

    memset(dist, UINT8_MAX, STATES);

    uint32_t head = 0;
    uint32_t tail = 1;

    queue[0] = 0;
    dist[0] = 0;
    *diameter = 0;

    while (head < tail) {
        uint32_t here = queue[head++];

        uint8_t d = dist[here];

        if (d > *diameter)
            *diameter = d;

        uint16_t p =
            (uint16_t) (here / ORIENTATIONS);

        uint16_t o =
            (uint16_t) (here % ORIENTATIONS);

        for (uint8_t move = 0; move < MOVES; ++move) {
            uint16_t next_p =
                bm_perm_move[move][p];

            uint16_t next_o =
                bm_ori_move[move][o];

            uint32_t there =
                (uint32_t) next_p * ORIENTATIONS
                + next_o;

            if (dist[there] != UINT8_MAX)
                continue;

            dist[there] = (uint8_t) (d + 1);
            queue[tail++] = there;
        }
    }

    if (tail != STATES) {
        fprintf(stderr,
                "full BFS incomplete: %u / %u states\n",
                tail,
                STATES);
        exit(1);
    }

    *distance11_count = 0;

    for (uint32_t rank = 0;
         rank < STATES;
         ++rank) {
        if (dist[rank] == 11)
            ++*distance11_count;
    }

    free(queue);

    return dist;
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

    /*
    * Stage-3 operation counters.
    */
    uint64_t transition_applications;
    uint64_t transition_table_loads;
    uint64_t pdb_table_loads;

    int solution_depth;

    /*
     * Pocket Cube HTM diameter = 11.
     */
    uint8_t path[12];
} bm_stats_t;

#define BM_MAX_DEPTH 11
#define BM_NO_FACE 3
#define BM_ORI_BITS 10
#define BM_ORI_MASK 0x3ffU
#define BM_FOUND (-1)

/*
 * Avoid move / 3 in the target-shaped hot loop.
 *
 * 0,1,2 -> R
 * 3,4,5 -> B
 * 6,7,8 -> D
 */
static const uint8_t bm_move_face[MOVES] = {
    0, 0, 0,
    1, 1, 1,
    2, 2, 2,
};

/*
 * Forward declarations used by the packed benchmark.
 */
static void bm_rank_to_string(uint32_t rank,
                              char output[15]);

static bm_stats_t bm_solve_iterative_transition_design(
    uint16_t p,
    uint16_t o,
    bm_transition_design_t design);

typedef uint32_t bm_packed_state_t;

static inline bm_packed_state_t bm_pack_state(uint16_t p, uint16_t o)
{
    return ((uint32_t) p << BM_ORI_BITS) | o;
}

static inline uint16_t bm_state_perm(bm_packed_state_t state)
{
    return (uint16_t) (state >> BM_ORI_BITS);
}

static inline uint16_t bm_state_ori(bm_packed_state_t state)
{
    return (uint16_t) (state & BM_ORI_MASK);
}


typedef struct {
    uint16_t p;
    uint16_t o;

    /*
     * Face used to reach this node.
     * 0 = R, 1 = B, 2 = D, 3 = root.
     */
    uint8_t prev_face;

    /*
     * Next move to try when we return to this frame.
     */
    uint8_t next_move;

    /*
     * Minimum f-value that exceeded the current threshold
     * in this subtree.
     */
    int min_excess;
} bm_frame_t;

typedef struct {
    bm_packed_state_t state;

    uint8_t prev_face;
    uint8_t next_move;
    uint8_t min_excess;
    uint8_t padding;
} bm_packed_frame_t;

static uint8_t bm_packed_heuristic(
    bm_packed_state_t state,
    bm_stats_t *stats)
{
    uint16_t p =
        bm_state_perm(state);

    uint16_t o =
        bm_state_ori(state);

    uint8_t hp =
        bm_perm_pdb[p];

    uint8_t ho =
        bm_ori_pdb[o];

    stats->pdb_table_loads += 2;

    return hp > ho ? hp : ho;
}

static bm_packed_state_t bm_packed_transition(
    bm_packed_state_t state,
    uint8_t move,
    bm_stats_t *stats)
{
    uint16_t p =
        bm_state_perm(state);

    uint16_t o =
        bm_state_ori(state);

    uint16_t next_p =
        bm_perm_move[move][p];

    uint16_t next_o =
        bm_ori_move[move][o];

    ++stats->transition_applications;
    stats->transition_table_loads += 2;

    return bm_pack_state(next_p, next_o);
}

static int bm_dfs_packed(
    bm_packed_state_t root,
    int bound,
    bm_stats_t *stats)
{
    bm_packed_frame_t stack[BM_MAX_DEPTH + 1];

    int depth = 0;

    int h =
        bm_packed_heuristic(
            root,
            stats);

    if (h > bound)
        return h;

    if (root == 0) {
        stats->solution_depth = 0;
        return BM_FOUND;
    }

    ++stats->expanded_total;
    ++stats->expanded_iter;

    stack[0] = (bm_packed_frame_t) {
        .state = root,
        .prev_face = BM_NO_FACE,
        .next_move = 0,
        .min_excess = UINT8_MAX,
        .padding = 0,
    };

    for (;;) {
        bm_packed_frame_t *frame =
            &stack[depth];

        int descended = 0;

        while (frame->next_move < MOVES) {
            uint8_t move =
                frame->next_move++;

            /*
             * No division by 3.
             */
            uint8_t face =
                bm_move_face[move];

            if (face == frame->prev_face)
                continue;

            ++stats->generated_total;
            ++stats->generated_iter;

            bm_packed_state_t next =
                bm_packed_transition(
                    frame->state,
                    move,
                    stats);

            int child_depth =
                depth + 1;

            stats->path[depth] =
                move;

            int child_h =
                bm_packed_heuristic(
                    next,
                    stats);

            int child_f =
                child_depth + child_h;

            if (child_f > bound) {
                if (child_f <
                    frame->min_excess) {
                    frame->min_excess =
                        (uint8_t) child_f;
                }

                continue;
            }

            /*
             * p == 0 && o == 0 is simply packed state == 0.
             */
            if (next == 0) {
                stats->solution_depth =
                    child_depth;

                return BM_FOUND;
            }

            if (child_depth >
                BM_MAX_DEPTH) {
                fputs(
                    "packed IDA* stack overflow\n",
                    stderr);

                exit(1);
            }

            ++stats->expanded_total;
            ++stats->expanded_iter;

            ++depth;

            stack[depth] =
                (bm_packed_frame_t) {
                    .state = next,
                    .prev_face = face,
                    .next_move = 0,
                    .min_excess = UINT8_MAX,
                    .padding = 0,
                };

            descended = 1;
            break;
        }

        if (descended)
            continue;

        int result =
            frame->min_excess;

        if (depth == 0)
            return result;

        --depth;

        if (result <
            stack[depth].min_excess) {
            stack[depth].min_excess =
                (uint8_t) result;
        }
    }
}

static bm_stats_t bm_solve_packed(
    uint16_t p,
    uint16_t o)
{
    bm_stats_t stats = {0};

    stats.solution_depth = -1;

    bm_packed_state_t root =
        bm_pack_state(p, o);

    /*
     * Initial IDA* threshold.
     */
    uint8_t hp =
        bm_perm_pdb[p];

    uint8_t ho =
        bm_ori_pdb[o];

    stats.pdb_table_loads += 2;

    int bound =
        hp > ho ? hp : ho;

    for (;;) {
        stats.expanded_iter = 0;
        stats.generated_iter = 0;

        int next =
            bm_dfs_packed(
                root,
                bound,
                &stats);

        if (next == BM_FOUND)
            return stats;

        if (next == UINT8_MAX) {
            fputs(
                "packed IDA* search failed\n",
                stderr);

            exit(1);
        }

        bound = next;
    }
}

static void bm_compare_packed_design(
    const uint8_t *full_dist,
    uint32_t d11_count)
{
    uint32_t count = 0;

    clock_t split_ticks = 0;
    clock_t packed_ticks = 0;

    for (uint32_t rank = 0;
         rank < STATES;
         ++rank) {

        if (full_dist[rank] != 11)
            continue;

        uint16_t p =
            (uint16_t)
            (rank / ORIENTATIONS);

        uint16_t o =
            (uint16_t)
            (rank % ORIENTATIONS);

        clock_t begin =
            clock();

        bm_stats_t split =
            bm_solve_iterative_transition_design(
                p,
                o,
                BM_TRANS_9);

        clock_t middle =
            clock();

        bm_stats_t packed =
            bm_solve_packed(
                p,
                o);

        clock_t end =
            clock();

        split_ticks +=
            middle - begin;

        packed_ticks +=
            end - middle;

        if (split.solution_depth !=
                packed.solution_depth ||
            split.expanded_total !=
                packed.expanded_total ||
            split.generated_total !=
                packed.generated_total ||
            split.transition_applications !=
                packed.transition_applications ||
            split.transition_table_loads !=
                packed.transition_table_loads ||
            split.pdb_table_loads !=
                packed.pdb_table_loads) {

            char state_string[15];

            bm_rank_to_string(
                rank,
                state_string);

            fprintf(
                stderr,
                "packed-state mismatch on %s\n",
                state_string);

            exit(1);
        }

        ++count;

        if (count % 100 == 0) {
            fprintf(
                stderr,
                "\rtested %u / %u",
                count,
                d11_count);

            fflush(stderr);
        }
    }

    fprintf(stderr, "\n");

    if (count != d11_count) {
        fprintf(
            stderr,
            "distance-11 count mismatch:"
            " %u != %u\n",
            count,
            d11_count);

        exit(1);
    }

    printf(
        "\n"
        "===== split vs packed state =====\n");

    printf(
        "distance-11 states = %u\n",
        count);

    printf(
        "split-state CPU time = %.3f s\n",
        (double) split_ticks /
            CLOCKS_PER_SEC);

    printf(
        "packed-state CPU time = %.3f s\n",
        (double) packed_ticks /
            CLOCKS_PER_SEC);

    printf(
        "packed frame size = %zu bytes\n",
        sizeof(bm_packed_frame_t));

    printf(
        "explicit stack size = %zu bytes\n",
        sizeof(bm_packed_frame_t)
            * (BM_MAX_DEPTH + 1));

    printf(
        "all search counters match = yes\n");
}

static void bm_apply_transition(
    bm_transition_design_t design,
    uint8_t move,
    uint16_t p,
    uint16_t o,
    uint16_t *next_p,
    uint16_t *next_o,
    bm_stats_t *stats)
{
    uint8_t face = (uint8_t) (move / 3U);
    uint8_t turn = (uint8_t) (move % 3U);

    /*
     * 3-table design
     *
     * R  -> 1 application
     * R2 -> 2 applications
     * R' -> 3 applications
     */
    if (design == BM_TRANS_3) {
        uint8_t applications =
            (uint8_t) (turn + 1U);

        uint16_t p2 = p;
        uint16_t o2 = o;

        for (uint8_t i = 0;
             i < applications;
             ++i) {
            p2 = bm_perm_move3[face][p2];
            o2 = bm_ori_move3[face][o2];

            ++stats->transition_applications;

            /*
             * One permutation load
             * + one orientation load.
             */
            stats->transition_table_loads += 2;
        }

        *next_p = p2;
        *next_o = o2;

        return;
    }

    /*
     * 6-table design
     *
     * R  -> one quarter-turn lookup
     * R2 -> one half-turn lookup
     * R' -> R2 + R, two applications
     */
    if (design == BM_TRANS_6) {
        uint16_t p2;
        uint16_t o2;

        if (turn == 0) {
            uint8_t index =
                (uint8_t) (face * 2U);

            p2 = bm_perm_move6[index][p];
            o2 = bm_ori_move6[index][o];

            ++stats->transition_applications;
            stats->transition_table_loads += 2;
        } else if (turn == 1) {
            uint8_t index =
                (uint8_t) (face * 2U + 1U);

            p2 = bm_perm_move6[index][p];
            o2 = bm_ori_move6[index][o];

            ++stats->transition_applications;
            stats->transition_table_loads += 2;
        } else {
            /*
             * inverse = half-turn followed by quarter-turn
             *
             * R' = R2 R
             */
            uint8_t half =
                (uint8_t) (face * 2U + 1U);

            uint8_t quarter =
                (uint8_t) (face * 2U);

            p2 = bm_perm_move6[half][p];
            o2 = bm_ori_move6[half][o];

            ++stats->transition_applications;
            stats->transition_table_loads += 2;

            p2 = bm_perm_move6[quarter][p2];
            o2 = bm_ori_move6[quarter][o2];

            ++stats->transition_applications;
            stats->transition_table_loads += 2;
        }

        *next_p = p2;
        *next_o = o2;

        return;
    }

    /*
     * 9-table design:
     *
     * every HTM move is one transition.
     */
    *next_p = bm_perm_move[move][p];
    *next_o = bm_ori_move[move][o];

    ++stats->transition_applications;
    stats->transition_table_loads += 2;
}

static uint8_t bm_combined_heuristic_counted(
    uint16_t p,
    uint16_t o,
    bm_stats_t *stats)
{
    uint8_t hp = bm_perm_pdb[p];
    uint8_t ho = bm_ori_pdb[o];

    stats->pdb_table_loads += 2;

    return hp > ho ? hp : ho;
}

static int bm_dfs_transition_design(
    uint16_t p,
    uint16_t o,
    int g,
    int bound,
    int prev_face,
    bm_transition_design_t design,
    bm_stats_t *stats)
{
    int h =
        bm_combined_heuristic_counted(
            p,
            o,
            stats);

    int f = g + h;

    if (f > bound)
        return f;

    if (p == 0 && o == 0) {
        stats->solution_depth = g;
        return BM_FOUND;
    }

    ++stats->expanded_total;
    ++stats->expanded_iter;

    int next_bound = INT_MAX;

    for (uint8_t move = 0;
         move < MOVES;
         ++move) {
        int face = move / 3;

        /*
         * Same-face pruning.
         */
        if (face == prev_face)
            continue;

        ++stats->generated_total;
        ++stats->generated_iter;

        uint16_t next_p;
        uint16_t next_o;

        bm_apply_transition(
            design,
            move,
            p,
            o,
            &next_p,
            &next_o,
            stats);

        stats->path[g] = move;

        int result =
            bm_dfs_transition_design(
                next_p,
                next_o,
                g + 1,
                bound,
                face,
                design,
                stats);

        if (result == BM_FOUND)
            return BM_FOUND;

        if (result < next_bound)
            next_bound = result;
    }

    return next_bound;
}

static bm_stats_t bm_solve_transition_design(
    uint16_t p,
    uint16_t o,
    bm_transition_design_t design)
{
    bm_stats_t stats = {0};

    stats.solution_depth = -1;

    /*
     * Root heuristic.
     */
    uint8_t hp = bm_perm_pdb[p];
    uint8_t ho = bm_ori_pdb[o];

    stats.pdb_table_loads += 2;

    int bound =
        hp > ho ? hp : ho;

    for (;;) {
        stats.expanded_iter = 0;
        stats.generated_iter = 0;

        int next =
            bm_dfs_transition_design(
                p,
                o,
                0,
                bound,
                -1,
                design,
                &stats);

        if (next == BM_FOUND)
            return stats;

        if (next == INT_MAX) {
            fputs("IDA* search failed\n", stderr);
            exit(1);
        }

        bound = next;
    }
}

static int bm_dfs_iterative_transition_design(
    uint16_t root_p,
    uint16_t root_o,
    int bound,
    bm_transition_design_t design,
    bm_stats_t *stats)
{
    bm_frame_t stack[BM_MAX_DEPTH + 1];

    int depth = 0;

    /*
     * Evaluate the root exactly as recursive bm_dfs() does.
     */
    int h =
        bm_combined_heuristic_counted(
            root_p,
            root_o,
            stats);

    int f = h;

    if (f > bound)
        return f;

    if (root_p == 0 && root_o == 0) {
        stats->solution_depth = 0;
        return BM_FOUND;
    }

    ++stats->expanded_total;
    ++stats->expanded_iter;

    stack[0] = (bm_frame_t) {
        .p = root_p,
        .o = root_o,
        .prev_face = BM_NO_FACE,
        .next_move = 0,
        .min_excess = INT_MAX,
    };

    for (;;) {
        bm_frame_t *frame =
            &stack[depth];

        int descended = 0;

        /*
         * Resume this frame from the next move that has
         * not yet been examined.
         */
        while (frame->next_move < MOVES) {
            uint8_t move =
                frame->next_move++;

            uint8_t face =
                (uint8_t) (move / 3U);

            /*
             * Same-face pruning.
             */
            if (face == frame->prev_face)
                continue;

            ++stats->generated_total;
            ++stats->generated_iter;

            uint16_t next_p;
            uint16_t next_o;

            bm_apply_transition(
                design,
                move,
                frame->p,
                frame->o,
                &next_p,
                &next_o,
                stats);

            int child_depth =
                depth + 1;

            /*
             * Store the move producing this child.
             */
            stats->path[depth] = move;

            /*
             * Evaluate the child before pushing it.
             *
             * This is equivalent to entering the recursive
             * bm_dfs() call and immediately evaluating f.
             */
            int child_h =
                bm_combined_heuristic_counted(
                    next_p,
                    next_o,
                    stats);

            int child_f =
                child_depth + child_h;

            if (child_f > bound) {
                if (child_f < frame->min_excess)
                    frame->min_excess = child_f;

                continue;
            }

            if (next_p == 0 && next_o == 0) {
                stats->solution_depth =
                    child_depth;

                return BM_FOUND;
            }

            /*
             * A non-solved state at depth 11 must have h >= 1,
             * and therefore cannot pass a bound <= 11.
             *
             * Keep this check anyway to guard the fixed stack.
             */
            if (child_depth > BM_MAX_DEPTH) {
                fputs("explicit IDA* stack overflow\n",
                      stderr);
                exit(1);
            }

            ++stats->expanded_total;
            ++stats->expanded_iter;

            ++depth;

            stack[depth] = (bm_frame_t) {
                .p = next_p,
                .o = next_o,
                .prev_face = face,
                .next_move = 0,
                .min_excess = INT_MAX,
            };

            descended = 1;
            break;
        }

        /*
         * We just pushed a child. Continue DFS from it.
         */
        if (descended)
            continue;

        /*
         * No moves remain in this frame:
         * equivalent to returning from recursive bm_dfs().
         */
        int result =
            frame->min_excess;

        /*
         * Returning from root means this entire IDA*
         * iteration has finished.
         */
        if (depth == 0)
            return result;

        --depth;

        /*
         * Propagate the minimum exceeded f-value to parent.
         */
        if (result < stack[depth].min_excess)
            stack[depth].min_excess = result;
    }
}

static bm_stats_t bm_solve_iterative_transition_design(
    uint16_t p,
    uint16_t o,
    bm_transition_design_t design)
{
    bm_stats_t stats = {0};

    stats.solution_depth = -1;

    /*
     * Initial IDA* threshold.
     */
    uint8_t hp = bm_perm_pdb[p];
    uint8_t ho = bm_ori_pdb[o];

    stats.pdb_table_loads += 2;

    int bound =
        hp > ho ? hp : ho;

    for (;;) {
        stats.expanded_iter = 0;
        stats.generated_iter = 0;

        int next =
            bm_dfs_iterative_transition_design(
                p,
                o,
                bound,
                design,
                &stats);

        if (next == BM_FOUND)
            return stats;

        if (next == INT_MAX) {
            fputs("iterative IDA* search failed\n",
                  stderr);
            exit(1);
        }

        bound = next;
    }
}

static void bm_compare_stack_design(
    const uint8_t *full_dist,
    uint32_t d11_count)
{
    uint32_t count = 0;

    uint64_t sum_expanded = 0;
    uint64_t sum_generated = 0;

    clock_t recursive_ticks = 0;
    clock_t iterative_ticks = 0;

    for (uint32_t rank = 0;
         rank < STATES;
         ++rank) {

        if (full_dist[rank] != 11)
            continue;

        uint16_t p =
            (uint16_t)
            (rank / ORIENTATIONS);

        uint16_t o =
            (uint16_t)
            (rank % ORIENTATIONS);

        /*
         * Use the selected Stage-3 candidate:
         *
         *   combined PDB
         *   + 9-move transitions
         */
        clock_t begin =
            clock();

        bm_stats_t recursive =
            bm_solve_transition_design(
                p,
                o,
                BM_TRANS_9);

        clock_t middle =
            clock();

        bm_stats_t iterative =
            bm_solve_iterative_transition_design(
                p,
                o,
                BM_TRANS_9);

        clock_t end =
            clock();

        recursive_ticks +=
            middle - begin;

        iterative_ticks +=
            end - middle;

        /*
         * The two implementations must produce identical
         * search behavior.
         */
        if (recursive.solution_depth !=
                iterative.solution_depth ||
            recursive.expanded_total !=
                iterative.expanded_total ||
            recursive.generated_total !=
                iterative.generated_total ||
            recursive.transition_applications !=
                iterative.transition_applications ||
            recursive.transition_table_loads !=
                iterative.transition_table_loads ||
            recursive.pdb_table_loads !=
                iterative.pdb_table_loads) {

            char state_string[15];

            bm_rank_to_string(
                rank,
                state_string);

            fprintf(stderr,
                    "\nstack mismatch on %s\n"
                    "recursive:"
                    " depth=%d"
                    " expanded=%" PRIu64
                    " generated=%" PRIu64
                    " transitions=%" PRIu64
                    " loads=%" PRIu64
                    " pdb=%" PRIu64 "\n"
                    "iterative:"
                    " depth=%d"
                    " expanded=%" PRIu64
                    " generated=%" PRIu64
                    " transitions=%" PRIu64
                    " loads=%" PRIu64
                    " pdb=%" PRIu64 "\n",
                    state_string,

                    recursive.solution_depth,
                    recursive.expanded_total,
                    recursive.generated_total,
                    recursive.transition_applications,
                    recursive.transition_table_loads,
                    recursive.pdb_table_loads,

                    iterative.solution_depth,
                    iterative.expanded_total,
                    iterative.generated_total,
                    iterative.transition_applications,
                    iterative.transition_table_loads,
                    iterative.pdb_table_loads);

            exit(1);
        }

        if (iterative.solution_depth != 11) {
            fprintf(stderr,
                    "unexpected solution depth %d\n",
                    iterative.solution_depth);
            exit(1);
        }

        sum_expanded +=
            iterative.expanded_total;

        sum_generated +=
            iterative.generated_total;

        ++count;

        if (count % 100 == 0) {
            fprintf(stderr,
                    "\rtested %u / %u",
                    count,
                    d11_count);

            fflush(stderr);
        }
    }

    fprintf(stderr, "\n");

    if (count != d11_count) {
        fprintf(stderr,
                "distance-11 count mismatch:"
                " %u != %u\n",
                count,
                d11_count);

        exit(1);
    }

    double recursive_seconds =
        (double) recursive_ticks /
        CLOCKS_PER_SEC;

    double iterative_seconds =
        (double) iterative_ticks /
        CLOCKS_PER_SEC;

    printf("\n"
           "===== recursive vs explicit stack =====\n");

    printf("distance-11 states = %u\n",
           count);

    printf("average expanded = %.2f\n",
           (double) sum_expanded / count);

    printf("average generated = %.2f\n",
           (double) sum_generated / count);

    printf("recursive CPU time = %.3f s\n",
           recursive_seconds);

    printf("explicit-stack CPU time = %.3f s\n",
           iterative_seconds);

    printf("all search counters match = yes\n");
}

static const char *bm_transition_design_name(
    bm_transition_design_t design)
{
    switch (design) {
    case BM_TRANS_3:
        return "3-move";
    case BM_TRANS_6:
        return "6-move";
    case BM_TRANS_9:
        return "9-move";
    }

    return "?";
}


static size_t bm_transition_static_bytes(
    bm_transition_design_t design)
{
    size_t moves;

    switch (design) {
    case BM_TRANS_3:
        moves = 3;
        break;

    case BM_TRANS_6:
        moves = 6;
        break;

    case BM_TRANS_9:
    default:
        moves = 9;
        break;
    }

    return moves
           * (PERMUTATIONS + ORIENTATIONS)
           * sizeof(uint16_t);
}

static void bm_compare_transition_design(
    const uint8_t *full_dist,
    uint32_t d11_count,
    bm_transition_design_t design)
{
    uint64_t sum_expanded = 0;
    uint64_t sum_generated = 0;

    uint64_t sum_apps = 0;
    uint64_t sum_transition_loads = 0;
    uint64_t sum_pdb_loads = 0;

    uint64_t worst_expanded = 0;
    uint64_t worst_apps = 0;
    uint64_t worst_transition_loads = 0;

    uint32_t worst_expanded_rank = 0;
    uint32_t worst_apps_rank = 0;

    uint32_t count = 0;

    clock_t begin = clock();

    for (uint32_t rank = 0;
         rank < STATES;
         ++rank) {
        if (full_dist[rank] != 11)
            continue;

        uint16_t p =
            (uint16_t)
            (rank / ORIENTATIONS);

        uint16_t o =
            (uint16_t)
            (rank % ORIENTATIONS);

        bm_stats_t stats =
            bm_solve_transition_design(
                p,
                o,
                design);

        if (stats.solution_depth != 11) {
            char state_string[15];

            bm_rank_to_string(
                rank,
                state_string);

            fprintf(stderr,
                    "%s: optimality failure on %s: %d\n",
                    bm_transition_design_name(design),
                    state_string,
                    stats.solution_depth);

            exit(1);
        }

        sum_expanded +=
            stats.expanded_total;

        sum_generated +=
            stats.generated_total;

        sum_apps +=
            stats.transition_applications;

        sum_transition_loads +=
            stats.transition_table_loads;

        sum_pdb_loads +=
            stats.pdb_table_loads;

        if (stats.expanded_total >
            worst_expanded) {
            worst_expanded =
                stats.expanded_total;

            worst_expanded_rank =
                rank;
        }

        if (stats.transition_applications >
            worst_apps) {
            worst_apps =
                stats.transition_applications;

            worst_transition_loads =
                stats.transition_table_loads;

            worst_apps_rank =
                rank;
        }

        ++count;
    }

    clock_t end = clock();

    if (count != d11_count) {
        fprintf(stderr,
                "count mismatch: %u != %u\n",
                count,
                d11_count);

        exit(1);
    }

    double seconds =
        (double) (end - begin)
        / CLOCKS_PER_SEC;

    char worst_expanded_state[15];
    char worst_apps_state[15];

    bm_rank_to_string(
        worst_expanded_rank,
        worst_expanded_state);

    bm_rank_to_string(
        worst_apps_rank,
        worst_apps_state);

    size_t transition_bytes =
        bm_transition_static_bytes(design);

    size_t pdb_bytes =
        PERMUTATIONS + ORIENTATIONS;

    size_t total_bytes =
        transition_bytes + pdb_bytes;

    printf("\n"
           "===== %s transition design =====\n",
           bm_transition_design_name(design));

    printf("distance-11 states = %u\n",
           count);

    printf("transition bytes = %zu"
           " (%.2f KiB)\n",
           transition_bytes,
           transition_bytes / 1024.0);

    printf("PDB bytes = %zu"
           " (%.2f KiB)\n",
           pdb_bytes,
           pdb_bytes / 1024.0);

    printf("total static bytes = %zu"
           " (%.2f KiB)\n",
           total_bytes,
           total_bytes / 1024.0);

    printf("average expanded = %.2f\n",
           (double) sum_expanded / count);

    printf("average generated = %.2f\n",
           (double) sum_generated / count);

    printf("average transition applications = %.2f\n",
           (double) sum_apps / count);

    printf("average transition table loads = %.2f\n",
           (double) sum_transition_loads / count);

    printf("average PDB table loads = %.2f\n",
           (double) sum_pdb_loads / count);

    printf("worst expanded = %" PRIu64
           " (%s)\n",
           worst_expanded,
           worst_expanded_state);

    printf("worst transition applications = %" PRIu64
           " (%s)\n",
           worst_apps,
           worst_apps_state);

    printf("transition loads for that case = %" PRIu64
           "\n",
           worst_transition_loads);

    printf("native CPU time = %.3f s\n",
           seconds);
}

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

static void bm_emit_rv32_tables(FILE *out)
{
    const size_t transition_bytes =
        MOVES *
        (PERMUTATIONS + ORIENTATIONS) *
        sizeof(uint16_t);

    const size_t pdb_bytes =
        PERMUTATIONS + ORIENTATIONS;

    const size_t move_face_bytes =
        MOVES;

    const size_t total_bytes =
        transition_bytes +
        pdb_bytes +
        move_face_bytes;

    fprintf(out,
            "# Generated by solver_heuristic.c\n"
            "# Do not edit manually.\n"
            "#\n"
            "# transition tables : %zu bytes\n"
            "# PDB tables        : %zu bytes\n"
            "# move-face table    : %zu bytes\n"
            "# total static data : %zu bytes (%.2f KiB)\n"
            "#\n",
            transition_bytes,
            pdb_bytes,
            move_face_bytes,
            total_bytes,
            total_bytes / 1024.0);

    /*
     * Ripes-compatible static-data section.
     *
     * These tables are immutable after assembly.
     */
    fprintf(out, ".data\n\n");

    /*
     * Each move has its own label.
     *
     * Layout:
     *
     * bm_perm_R[p]
     * bm_perm_R2[p]
     * ...
     */
    for (uint8_t move = 0; move < MOVES; ++move) {
        fprintf(out,
                "# permutation transition: %s\n"
                "bm_perm_%s:\n",
                bm_move_symbol[move],
                bm_move_symbol[move]);

        bm_emit_u16_array(
            out,
            bm_perm_move[move],
            PERMUTATIONS);

        fputc('\n', out);
    }

    /*
     * Orientation transition tables.
     */
    for (uint8_t move = 0; move < MOVES; ++move) {
        fprintf(out,
                "# orientation transition: %s\n"
                "bm_ori_%s:\n",
                bm_move_symbol[move],
                bm_move_symbol[move]);

        bm_emit_u16_array(
            out,
            bm_ori_move[move],
            ORIENTATIONS);

        fputc('\n', out);
    }

    /*
     * Pattern databases.
     */
    fprintf(out,
            "# permutation PDB\n"
            "bm_perm_pdb_data:\n");

    bm_emit_u8_array(
        out,
        bm_perm_pdb,
        PERMUTATIONS);

    fputc('\n', out);

    fprintf(out,
            "# orientation PDB\n"
            "bm_ori_pdb_data:\n");

    bm_emit_u8_array(
        out,
        bm_ori_pdb,
        ORIENTATIONS);

    fputc('\n', out);

    /*
     * move -> face
     *
     * 0,1,2 -> R
     * 3,4,5 -> B
     * 6,7,8 -> D
     */
    fprintf(out,
            "# move-to-face lookup\n"
            "bm_move_face_data:\n");

    bm_emit_u8_array(
        out,
        bm_move_face,
        MOVES);
}

static void bm_rank_to_string(uint32_t rank,
                              char output[15])
{
    state_t state;

    unrank_state(rank, &state);

    for (int i = 0; i < CUBIES; ++i)
        output[i] =
            (char) ('1' + state.p[i]);

    for (int i = 0; i < CUBIES; ++i)
        output[CUBIES + i] =
            (char) ('1' + state.o[i]);

    output[14] = '\0';
}

static void bm_emit_d11_packed(const uint8_t *full_dist)
{
    for (uint32_t rank = 0; rank < STATES; ++rank) {
        if (full_dist[rank] != 11)
            continue;

        uint16_t p =
            (uint16_t) (rank / ORIENTATIONS);

        uint16_t o =
            (uint16_t) (rank % ORIENTATIONS);

        uint32_t packed =
            ((uint32_t) p << BM_ORI_BITS) | o;

        char state_string[15];

        bm_rank_to_string(
            rank,
            state_string);

        /*
         * packed state rank p o
         */
        printf("%u\t%s\t%u\t%u\t%u\n",
               packed,
               state_string,
               rank,
               p,
               o);
    }
}

typedef struct {
    uint32_t rank;

    uint64_t expanded;
    uint64_t generated;
} bm_d11_result_t;

int main(int argc, char **argv)
{
    /*
     * Build host-side tables once.
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

    if (argc == 2 &&
        !strcmp(argv[1],
                "--emit-d11-packed")) {
    
        uint8_t full_diameter;
        uint32_t d11_count;
    
        uint8_t *full_dist =
            bm_build_full_distance(
                &full_diameter,
                &d11_count);
    
        if (full_diameter != 11 ||
            d11_count != 2644) {
            fprintf(stderr,
                    "unexpected full-state distribution\n");
    
            free(full_dist);
            return 1;
        }
    
        bm_emit_d11_packed(full_dist);
    
        free(full_dist);
        return 0;
    }

    /*
     * Generate static tables for the RV32I target.
     *
     * Nothing except assembly source is written to stdout.
     */
    if (argc == 2 &&
        !strcmp(argv[1],
                "--emit-rv32-tables")) {

        bm_emit_rv32_tables(stdout);
        return 0;
    }

    printf("orientation PDB max distance = %u\n",
           orientation_diameter);

    printf("permutation PDB max distance = %u\n",
           permutation_diameter);

    if (argc == 2
        && !strcmp(
            argv[1],
            "--compare-packed")) {
    
        uint8_t full_diameter;
        uint32_t d11_count;
    
        uint8_t *full_dist =
            bm_build_full_distance(
                &full_diameter,
                &d11_count);
    
        printf(
            "full-state diameter = %u\n",
            full_diameter);
    
        printf(
            "distance-11 states = %u\n",
            d11_count);
    
        bm_compare_packed_design(
            full_dist,
            d11_count);
    
        free(full_dist);
    
        return 0;
    } 
    if (argc == 2
        && !strcmp(argv[1], "--compare-stack")) {
    
        uint8_t full_diameter;
        uint32_t d11_count;
    
        uint8_t *full_dist =
            bm_build_full_distance(
                &full_diameter,
                &d11_count);
    
        printf("full-state diameter = %u\n",
               full_diameter);
    
        printf("distance-11 states = %u\n",
               d11_count);
    
        bm_compare_stack_design(
            full_dist,
            d11_count);
    
        free(full_dist);
    
        return 0;
    }


    /*
     * Full distance-11 benchmark.
     */
    if (argc == 2
    && !strcmp(argv[1], "--compare-transitions")) {

    uint8_t full_diameter;
    uint32_t d11_count;

    uint8_t *full_dist =
        bm_build_full_distance(
            &full_diameter,
            &d11_count);

    if (!full_dist) {
        fputs("failed to build full distance table\n",
              stderr);
        return 1;
    }

    printf("full-state diameter = %u\n",
           full_diameter);

    printf("distance-11 states = %u\n",
           d11_count);

    bm_compare_transition_design(
        full_dist,
        d11_count,
        BM_TRANS_3);

    bm_compare_transition_design(
        full_dist,
        d11_count,
        BM_TRANS_6);

    bm_compare_transition_design(
        full_dist,
        d11_count,
        BM_TRANS_9);

    free(full_dist);

    return 0;
}


    /*
     * Original single-state benchmark.
     */
    const char *input =
        argc == 2
            ? argv[1]
            : "21345671111111";

    state_t state;

    if (!parse_state(input, &state)) {
        fprintf(stderr,
                "usage:\n"
                "  %s PPPPPPPOOOOOOO\n"
                "  %s --all-d11\n",
                argv[0],
                argv[0]);

        return 2;
    }

    uint32_t rank =
        rank_state(&state);

    uint16_t p =
        (uint16_t)
        (rank / ORIENTATIONS);

    uint16_t o =
        (uint16_t)
        (rank % ORIENTATIONS);

    printf("state = %s\n", input);

    printf("permutation rank = %u\n", p);
    printf("orientation rank = %u\n", o);

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

