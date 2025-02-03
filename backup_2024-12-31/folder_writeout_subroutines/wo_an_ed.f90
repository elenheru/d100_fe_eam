subroutine wo_an_ed
    implicit none
    integer, parameter :: limit = 1699
    real(8), parameter :: step_length = 3d-3
    integer :: i
    real(8) :: ed_ack4pot,ewo,r
    open(302, file = "electrondensity_acklandmendelev4.txt")
    do i = 1, limit
        r = i * step_length
        ewo = ed_ack4pot(r)
        write(302,1302) r, ewo
    enddo
    close(302)
    print *, "Writing out Ackland-Mendelev fourth Fe electronic density in txt file, from 0 to ", limit * step_length
    1302 format(F8.5,1x,ES18.10e3)
endsubroutine wo_an_ed
