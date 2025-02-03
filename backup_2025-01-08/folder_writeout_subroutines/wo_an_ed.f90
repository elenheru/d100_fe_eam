subroutine wo_an_ed
    use compile_debug_options_mod, only : float_byte_size, integer_byte_size
    implicit none
    integer(integer_byte_size), parameter :: limit = 1699
    real(float_byte_size), parameter :: step_length = 3d-3
    integer(integer_byte_size) :: i
    real(float_byte_size) :: ed_ack4pot, ewo, r
    open(302, file = "electrondensity_acklandmendelev4.txt")
    do i = 1, limit
        r = i * step_length
        ewo = ed_ack4pot(r)
        write(302, 1302) r, ewo
    enddo
    close(302)
    print *, "Writing out Ackland-Mendelev fourth Fe electronic density in txt file, from 0 to ", limit * step_length
    1302 format(F8.5, 1x, ES18.10e3)
endsubroutine wo_an_ed
