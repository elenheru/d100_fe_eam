subroutine wo_an_mf
    use compile_debug_options_mod, only : float_byte_size, integer_byte_size
    integer(integer_byte_size), parameter :: limit = 799
    real(float_byte_size), parameter :: step_length = 3d-2
    integer(integer_byte_size) :: i
    real(float_byte_size) :: mf_ack4pot, pwo, rho
    open(303, file = "embedding_function_acklandmendelev4.txt")
    do i = 1, limit
        rho = i * step_length
        pwo = mf_ack4pot(rho)
        write(303, 1303) rho, pwo
    enddo
    close(303)
    print *, "Writing out Ackland-Mendelev fourth Fe embedding function in txt file, from 0 to ", limit * step_length
    1303 format(F10.5, 1x, ES18.10e3)
endsubroutine wo_an_mf
