subroutine wo_an_mf
    implicit none
    integer, parameter :: limit = 799
    real(8), parameter :: step_length = 3d-2
    integer :: i
    real(8) :: mf_ack4pot,pwo,rho
    open(303, file = "embedding_function_acklandmendelev4.txt")
    do i = 1, limit
        rho = i * step_length
        pwo = mf_ack4pot(rho)
        write(303,1303) rho, pwo
    enddo
    close(303)
    print *, "Writing out Ackland-Mendelev fourth Fe embedding function in txt file, from 0 to ", limit * step_length
    1303 format(F10.5,1x,ES18.10e3)
endsubroutine wo_an_mf
