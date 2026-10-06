# solver_rv32i_body.s
#
# RV32I Pocket Cube solver
#
# Design:
#   - IDA*
#   - h = max(permutation PDB, orientation PDB)
#   - same-face pruning
#   - 9 HTM transition tables
#   - packed state: (p << 10) | o
#   - explicit 12-frame stack
#
# generated_tables.s must appear before this file.

.data
.align 2

# ------------------------------------------------------------
# Test input
#
# 21345671111111:
#   permutation rank = 720
#   orientation rank = 0
#
# packed = (720 << 10) | 0 = 737280
#
# worst-case 54721631111111
# p = 3343, o = 0
# packed = (3343 << 10) | 0 = 3423232
# ------------------------------------------------------------

# root_state:
#     .word 737280

root_state:
	.word 3423232

# ------------------------------------------------------------
# Transition-table pointer arrays
#
# This avoids computing:
#
#   move * 5040
#   move * 729
#
# in the RV32I hot loop.
# ------------------------------------------------------------

perm_table_ptrs:
    .word bm_perm_R
    .word bm_perm_R2
    .word bm_perm_Ri
    .word bm_perm_B
    .word bm_perm_B2
    .word bm_perm_Bi
    .word bm_perm_D
    .word bm_perm_D2
    .word bm_perm_Di

ori_table_ptrs:
    .word bm_ori_R
    .word bm_ori_R2
    .word bm_ori_Ri
    .word bm_ori_B
    .word bm_ori_B2
    .word bm_ori_Bi
    .word bm_ori_D
    .word bm_ori_D2
    .word bm_ori_Di


# ------------------------------------------------------------
# Explicit DFS stack
#
# One frame = 8 bytes
#
# offset 0: uint32 packed state
# offset 4: uint8  previous face
# offset 5: uint8  next move
# offset 6: uint8  minimum exceeded f
# offset 7: unused
#
# 12 * 8 = 96 bytes
# ------------------------------------------------------------

stack_frames:
    .zero 96

solution_path:
    .zero 12


# ------------------------------------------------------------
# Output strings
# ------------------------------------------------------------

msg_depth:
    .string "solution depth = "

msg_solution:
    .string "solution = "

msg_fail:
    .string "search failed\n"

newline:
    .string "\n"

move_R:
    .string "R "

move_R2:
    .string "R2 "

move_Ri:
    .string "R' "

move_B:
    .string "B "

move_B2:
    .string "B2 "

move_Bi:
    .string "B' "

move_D:
    .string "D "

move_D2:
    .string "D2 "

move_Di:
    .string "D' "

move_name_ptrs:
    .word move_R
    .word move_R2
    .word move_Ri
    .word move_B
    .word move_B2
    .word move_Bi
    .word move_D
    .word move_D2
    .word move_Di


.text
.globl main

main:
    # --------------------------------------------------------
    # Register allocation
    #
    # s0  root packed state
    # s1  current IDA* threshold
    # s2  current DFS depth
    # s3  explicit stack base
    # s4  solution path base
    # s5  permutation table-pointer array
    # s6  orientation table-pointer array
    # s7  permutation PDB base
    # s8  orientation PDB base
    # s9  move -> face table
    # s10 move-name pointer array
    # s11 solution depth
    # --------------------------------------------------------

    la      s3, stack_frames
    la      s4, solution_path

    la      s5, perm_table_ptrs
    la      s6, ori_table_ptrs

    la      s7, bm_perm_pdb_data
    la      s8, bm_ori_pdb_data

    la      s9, bm_move_face_data
    la      s10, move_name_ptrs

    lw      s0, root_state


    # --------------------------------------------------------
    # Initial threshold:
    #
    # h(root) = max(
    #     perm_pdb[p],
    #     ori_pdb[o]
    # )
    #
    # packed state:
    #   p = state >> 10
    #   o = state & 0x3ff
    # --------------------------------------------------------

    srli    t0, s0, 10
    andi    t1, s0, 0x3ff

    add     t2, s7, t0
    lbu     t3, 0(t2)              # hp

    add     t2, s8, t1
    lbu     t4, 0(t2)              # ho

    mv      s1, t3

    bgeu    t3, t4, threshold_ready
    mv      s1, t4


threshold_ready:
    beq     s0, zero, solved_root


# ============================================================
# Start one IDA* threshold iteration
# ============================================================

start_iteration:
    li      s2, 0                  # depth = 0

    # root frame
    sw      s0, 0(s3)

    li      t0, 3                  # BM_NO_FACE
    sb      t0, 4(s3)

    sb      zero, 5(s3)            # next_move = 0

    li      t0, 255
    sb      t0, 6(s3)              # min_excess = INF


# ============================================================
# DFS main loop
# ============================================================

search_loop:
    # frame = stack + depth * 8

    slli    t0, s2, 3
    add     t0, s3, t0

    # next move

    lbu     t1, 5(t0)

    li      t2, 9
    bgeu    t1, t2, backtrack

    # frame->next_move++

    addi    t2, t1, 1
    sb      t2, 5(t0)


    # --------------------------------------------------------
    # face = move_face[move]
    # --------------------------------------------------------

    add     t2, s9, t1
    lbu     t3, 0(t2)              # face

    lbu     t4, 4(t0)              # previous face

    # same-face pruning

    beq     t3, t4, search_loop


    # --------------------------------------------------------
    # Load packed state and split p/o
    # --------------------------------------------------------

    lw      t4, 0(t0)

    srli    t5, t4, 10             # p
    andi    t6, t4, 0x3ff          # o


    # --------------------------------------------------------
    # permutation transition
    #
    # base = perm_table_ptrs[move]
    # next_p = base[p]
    # --------------------------------------------------------

    slli    a0, t1, 2

    add     a1, s5, a0
    lw      a1, 0(a1)

    slli    a2, t5, 1
    add     a1, a1, a2

    lhu     a3, 0(a1)              # next_p


    # --------------------------------------------------------
    # orientation transition
    # --------------------------------------------------------

    add     a1, s6, a0
    lw      a1, 0(a1)

    slli    a2, t6, 1
    add     a1, a1, a2

    lhu     a4, 0(a1)              # next_o


    # --------------------------------------------------------
    # Pack next state
    #
    # next = (next_p << 10) | next_o
    # --------------------------------------------------------

    slli    a5, a3, 10
    or      a5, a5, a4             # next packed state


    # child_depth = depth + 1

    addi    a6, s2, 1


    # --------------------------------------------------------
    # path[depth] = move
    # --------------------------------------------------------

    add     a7, s4, s2
    sb      t1, 0(a7)


    # --------------------------------------------------------
    # h(next)
    #
    # We already have next_p / next_o, so no need to unpack
    # the newly packed state again.
    # --------------------------------------------------------

    add     a0, s7, a3
    lbu     a0, 0(a0)              # hp

    add     a1, s8, a4
    lbu     a1, 0(a1)              # ho

    mv      a2, a0                 # h = hp

    bgeu    a0, a1, heuristic_ready
    mv      a2, a1                 # h = ho


heuristic_ready:
    # f = child_depth + h

    add     a2, a2, a6


    # --------------------------------------------------------
    # if f > bound:
    #     update frame.min_excess
    # --------------------------------------------------------

    bltu    s1, a2, prune_child


    # --------------------------------------------------------
    # Goal test
    #
    # packed state == 0 means p == 0 && o == 0
    # --------------------------------------------------------

    beq     a5, zero, found_solution


    # --------------------------------------------------------
    # Guard explicit stack
    # --------------------------------------------------------

    li      a0, 11
    bltu    a0, a6, search_failed


    # --------------------------------------------------------
    # Push child frame
    # --------------------------------------------------------

    mv      s2, a6

    slli    a0, s2, 3
    add     a0, s3, a0

    sw      a5, 0(a0)              # packed state
    sb      t3, 4(a0)              # previous face
    sb      zero, 5(a0)            # next_move = 0

    li      a1, 255
    sb      a1, 6(a0)              # min_excess = INF

    j       search_loop


# ============================================================
# Child exceeds IDA* threshold
# ============================================================

prune_child:
    # a2 = child f
    # t0 = current frame pointer

    lbu     a0, 6(t0)

    # if child_f >= current min_excess, nothing to do

    bgeu    a2, a0, search_loop

    sb      a2, 6(t0)

    j       search_loop


# ============================================================
# Finished all moves from current frame
# ============================================================

backtrack:
    # result = frame.min_excess

    lbu     t1, 6(t0)

    # root finished -> this threshold failed

    beq     s2, zero, next_threshold


    # pop

    addi    s2, s2, -1

    slli    t2, s2, 3
    add     t2, s3, t2             # parent frame

    lbu     t3, 6(t2)              # parent min_excess

    # parent.min_excess =
    #     min(parent.min_excess, result)

    bgeu    t1, t3, search_loop

    sb      t1, 6(t2)

    j       search_loop


# ============================================================
# Increase IDA* threshold
# ============================================================

next_threshold:
    li      t2, 255

    beq     t1, t2, search_failed

    mv      s1, t1

    j       start_iteration


# ============================================================
# Found solution
# ============================================================

found_solution:
    mv      s11, a6
    j       print_result


solved_root:
    li      s11, 0


# ============================================================
# Output
# ============================================================

print_result:
    # "solution depth = "

    la      a0, msg_depth
    li      a7, 4
    ecall

    # integer depth

    mv      a0, s11
    li      a7, 1
    ecall

    # newline

    la      a0, newline
    li      a7, 4
    ecall


    # "solution = "

    la      a0, msg_solution
    li      a7, 4
    ecall


    # print each move

    li      t0, 0

print_move_loop:
    bgeu    t0, s11, print_done

    add     t1, s4, t0
    lbu     t2, 0(t1)              # move index

    slli    t2, t2, 2
    add     t3, s10, t2

    lw      a0, 0(t3)

    li      a7, 4
    ecall

    addi    t0, t0, 1
    j       print_move_loop


print_done:
    la      a0, newline
    li      a7, 4
    ecall

    li      a7, 10
    ecall


# ============================================================
# Failure
# ============================================================

search_failed:
    la      a0, msg_fail
    li      a7, 4
    ecall

    li      a7, 10
    ecall
