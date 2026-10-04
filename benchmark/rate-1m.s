.text
main:
    li t0, 1000000
    li t1, 0

loop:
    addi t1, t1, 1
    addi t0, t0, -1
    bnez t0, loop

    li a7, 10
    ecall
