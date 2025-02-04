subroutine shift_atom(na, magn)
    use positions_mod
    implicit none
    integer(integer_byte_size), intent(in) :: na
    real(float_byte_size), intent(in) :: magn
    real(float_byte_size) :: shift
    integer(integer_byte_size) i
    do i = 1, 3
        call random_number(shift)
        shift = (shift - 5d-1) * magn * 2d0
        R_at(i, na) = R_at(i, na) + shift
    enddo
endsubroutine shift_atom
