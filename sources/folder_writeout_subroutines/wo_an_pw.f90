subroutine wo_an_pw
    use compile_debug_options_mod, only : float_byte_size, integer_byte_size
    integer(integer_byte_size), parameter :: limit = 1499
    real(float_byte_size), parameter :: step_length = 3d-3
    integer(integer_byte_size) :: i
    real(float_byte_size) :: pw_ack4pot, pwo, r
    open(301, file = "pairwise_acklandmendelev4.txt")
    do i = 1, limit
        r = i * step_length
        pwo = pw_ack4pot(r)
        write(301, 1301) r, pwo
    enddo
    close(301)
    print *, "Writing out Ackland-Mendelev fourth Fe pairwise potential in txt file, from 0 to ", limit * step_length
    1301 format(F8.5, 1x, ES18.10e3)
endsubroutine wo_an_pw
