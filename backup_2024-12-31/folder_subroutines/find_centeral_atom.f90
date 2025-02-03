subroutine find_centeral_atom(na)
    use positions_mod
    implicit none
    integer, intent(inout) :: na
    integer i,imin
    real(8) :: rmin
    imin = -1
    rmin = huge(rmin)

    do i =1,nat_mc
        if (norm2(R_at(1:3,i)) .gt. rmin) cycle
        rmin = norm2(R_at(1:3,i))
        imin = i
    enddo
    na = imin
    !print*, "centeral atom is ", na
endsubroutine find_centeral_atom
