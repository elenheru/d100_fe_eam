subroutine shift_atom(na,magn)
    use positions_mod
    implicit none
    integer,intent(in) :: na
    real(8),intent(in) :: magn
    real(8) :: shift
    integer i
    do i = 1,3
        call random_number(shift)
        shift = (shift - 5d-1) * magn * 2d0
        R_at(i,na) = R_at(i,na) + shift
    enddo
endsubroutine shift_atom
