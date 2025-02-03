subroutine wo_an_pw
    implicit none
    integer, parameter :: limit = 1499
    real(8), parameter :: step_length = 3d-3
    integer :: i
    real(8) :: pw_ack4pot,pwo,r
    open(301, file = "pairwise_acklandmendelev4.txt")
    do i = 1, limit
        r = i * step_length
        pwo = pw_ack4pot(r)
        write(301,1301) r, pwo
    enddo
    close(301)
    print *, "Writing out Ackland-Mendelev fourth Fe pairwise potential in txt file, from 0 to ", limit * step_length
    1301 format(F8.5,1x,ES18.10e3)
endsubroutine wo_an_pw
