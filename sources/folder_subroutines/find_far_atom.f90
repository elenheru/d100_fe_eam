subroutine find_far_atom(na)
    use positions_mod
    implicit none
    integer(integer_byte_size), intent(inout) :: na
    integer(integer_byte_size) i, imax
    real(float_byte_size) :: rmax
    imax = -1
    rmax = 0d0
    do i =1, last_relaxable
        if ( norm2( R_at(1:3, i) ) .lt. rmax) cycle
        rmax = norm2( R_at(1:3, i) )
        imax = i
    enddo
    na = imax
endsubroutine find_far_atom
