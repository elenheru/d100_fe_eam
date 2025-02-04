subroutine find_centeral_atom(na)
    use positions_mod
    implicit none
    integer(integer_byte_size), intent(inout) :: na
    integer(integer_byte_size) i, imin
    real(float_byte_size) :: rmin
    imin = -1
    rmin = huge(rmin)

    do i =1, last_relaxable
        if (norm2(R_at(1:3, i)) .gt. rmin) cycle
        rmin = norm2(R_at(1:3, i))
        imin = i
    enddo
    na = imin
    !print*, "centeral atom is ", na
endsubroutine find_centeral_atom
