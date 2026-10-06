#!/usr/bin/env bash

set -u

RIPES="${RIPES:-../Ripes/build/Ripes.app/Contents/MacOS/Ripes}"
GENERATOR="${GENERATOR:-./solver_heuristic}"

TABLES="${TABLES:-generated_tables.s}"
BODY="${BODY:-solver_rv32i_body.s}"

JOBS="${JOBS:-4}"

RESULT_DIR="${RESULT_DIR:-benchmark_rv32_results}"
FINAL_RESULT="${FINAL_RESULT:-benchmark_rv32.tsv}"

SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"


# ------------------------------------------------------------
# Worker mode
# ------------------------------------------------------------

if [[ "${1:-}" == "--worker" ]]; then
    shift

    packed="$1"
    state="$2"
    rank="$3"
    p="$4"
    o="$5"

    mkdir -p "$RESULT_DIR"

    result_file="$RESULT_DIR/$state.tsv"

    # Resume support.
    if [[ -s "$result_file" ]]; then
        exit 0
    fi

    tmpdir="$(mktemp -d)"

    body="$tmpdir/body.s"
    source="$tmpdir/solver.s"
    output="$tmpdir/output.txt"

    # Replace only the real root_state .word.
    awk -v packed="$packed" '
        /^root_state:$/ {
            print
            getline
            print "    .word " packed
            next
        }

        { print }
    ' "$BODY" > "$body"

    cat "$TABLES" "$body" > "$source"

    "$RIPES" \
        --mode cli \
        --src "$source" \
        -t asm \
        --proc RV32_ISS \
        --iret \
        --cycles \
        > "$output" 2>&1

    status=$?

    if [[ "$status" -ne 0 ]]; then
        echo "Ripes failed for $state" >&2
        cat "$output" >&2
        rm -rf "$tmpdir"
        exit 1
    fi

    depth="$(
        awk -F'= ' '
            /^solution depth = / {
                print $2
                exit
            }
        ' "$output"
    )"

    instructions="$(
        awk '
            /^===== instructions retired$/ {
                getline
                print
                exit
            }
        ' "$output"
    )"

    cycles="$(
        awk '
            /^===== cycles$/ {
                getline
                print
                exit
            }
        ' "$output"
    )"

    if [[ -z "$depth" ||
          -z "$instructions" ||
          -z "$cycles" ]]; then
        echo "Unable to parse Ripes output for $state" >&2
        cat "$output" >&2
        rm -rf "$tmpdir"
        exit 1
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$state" \
        "$packed" \
        "$rank" \
        "$p" \
        "$o" \
        "$depth" \
        "$instructions" \
        "$cycles" \
        > "$result_file"

    rm -rf "$tmpdir"

    exit 0
fi


# ------------------------------------------------------------
# Controller
# ------------------------------------------------------------

for f in \
    "$RIPES" \
    "$GENERATOR" \
    "$TABLES" \
    "$BODY"
do
    if [[ ! -e "$f" ]]; then
        echo "missing: $f" >&2
        exit 1
    fi
done

mkdir -p "$RESULT_DIR"

state_file="$(mktemp)"
trap 'rm -f "$state_file"' EXIT

"$GENERATOR" \
    --emit-d11-packed \
    > "$state_file"

count="$(wc -l < "$state_file" | tr -d ' ')"

if [[ "$count" != "2644" ]]; then
    echo "expected 2644 distance-11 states, got $count" >&2
    exit 1
fi

echo "distance-11 states: $count"
echo "parallel workers:    $JOBS"
echo "results:             $RESULT_DIR/"
echo

# Each row:
#
# packed state rank p o

cat "$state_file" |
    xargs -n 5 -P "$JOBS" \
        "$SELF" --worker


# ------------------------------------------------------------
# Aggregate
# ------------------------------------------------------------

printf \
"state\tpacked\trank\tp\to\tdepth\tinstructions\tcycles\n" \
    > "$FINAL_RESULT"

cat "$RESULT_DIR"/*.tsv |
    sort -k1,1 \
    >> "$FINAL_RESULT"

actual="$(
    find "$RESULT_DIR" \
        -name '*.tsv' |
        wc -l |
        tr -d ' '
)"

if [[ "$actual" != "2644" ]]; then
    echo "only $actual / 2644 results completed" >&2
    exit 1
fi


# ------------------------------------------------------------
# Correctness
# ------------------------------------------------------------

bad_depth="$(
    awk -F '\t' '
        $6 != 11 {
            ++n
        }

        END {
            print n + 0
        }
    ' "$RESULT_DIR"/*.tsv
)"

over_budget="$(
    awk -F '\t' '
        $7 >= 50000000 {
            ++n
        }

        END {
            print n + 0
        }
    ' "$RESULT_DIR"/*.tsv
)"


# ------------------------------------------------------------
# Statistics
# ------------------------------------------------------------

awk -F '\t' '
BEGIN {
    max = -1
    sum = 0
    count = 0
}

{
    instructions = $7 + 0

    sum += instructions
    ++count

    if (instructions > max) {
        max = instructions

        worst_state = $1
        worst_packed = $2
        worst_rank = $3
        worst_p = $4
        worst_o = $5
        worst_depth = $6
        worst_cycles = $8
    }
}

END {
    printf "\n"
    printf "===== RV32_ISS exhaustive distance-11 benchmark =====\n"
    printf "count = %d\n", count
    printf "average instructions = %.2f\n", sum / count
    printf "\n"
    printf "===== worst case =====\n"
    printf "state = %s\n", worst_state
    printf "packed = %s\n", worst_packed
    printf "rank = %s\n", worst_rank
    printf "permutation rank = %s\n", worst_p
    printf "orientation rank = %s\n", worst_o
    printf "solution depth = %s\n", worst_depth
    printf "instructions retired = %d\n", max
    printf "cycles = %s\n", worst_cycles
    printf "budget utilization = %.2f%%\n",
           max / 50000000 * 100
}
' "$RESULT_DIR"/*.tsv

echo
echo "incorrect solution depths = $bad_depth"
echo "states >= 50M instructions = $over_budget"
echo
echo "full results: $FINAL_RESULT"
