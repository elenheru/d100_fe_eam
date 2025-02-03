subroutine find_far_atom(na)
    use positions_mod
    implicit none
    integer, intent(inout) :: na
    integer i,imax
    real(8) :: rmax
    imax = -1
    rmax = 0d0
    do i =1,nat_mc
        if ( norm2( R_at(1:3,i) ) .lt. rmax) cycle
        rmax = norm2( R_at(1:3,i) )
        imax = i
    enddo
    na = imax
endsubroutine find_far_atom
