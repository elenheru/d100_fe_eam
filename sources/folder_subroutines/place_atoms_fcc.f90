subroutine place_atoms_fcc(lp)
    use positions_mod
    implicit none
    integer(integer_byte_size), parameter :: iw=10
    integer(integer_byte_size) i, j, k, na
    real(float_byte_size), intent(in) :: lp
    na=0
    STOP "not used in reconstruction, only in potential investigation"
    do i=-iw, iw
    do j=-iw, iw
    do k=-iw, iw
        na=na+1
        R_at(1, na)=i*lp
        R_at(2, na)=j*lp
        R_at(3, na)=k*lp
        na=na+1
        R_at(1, na)=i*lp
        R_at(2, na)=j*lp+lp*0.5
        R_at(3, na)=k*lp+lp*0.5
        na=na+1
        R_at(1, na)=i*lp+lp*0.5
        R_at(2, na)=j*lp
        R_at(3, na)=k*lp+lp*0.5
        na=na+1
        R_at(1, na)=i*lp+lp*0.5
        R_at(2, na)=j*lp+lp*0.5
        R_at(3, na)=k*lp
    enddo
    enddo
    enddo
    !print*, na, " atoms placed in FCC"
    last_wall = na
endsubroutine place_atoms_fcc
