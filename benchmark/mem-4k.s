.text
main:
    li t0, 0x10000000
    li t1, 4096

loop:
    sb zero, 0(t0)
    addi t0, t0, 1
    addi t1, t1, -1
    bnez t1, loop

    li a7, 10
    ecall
