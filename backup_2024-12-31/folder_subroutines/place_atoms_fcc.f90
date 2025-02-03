subroutine place_atoms_fcc(lp)
    use positions_mod
    implicit none
    integer,parameter :: iw=10
    integer i,j,k,na
    real(8), intent(in) :: lp
    na=0
    do i=-iw,iw
    do j=-iw,iw
    do k=-iw,iw
        na=na+1
        R_at(1,na)=i*lp
        R_at(2,na)=j*lp
        R_at(3,na)=k*lp
        na=na+1
        R_at(1,na)=i*lp
        R_at(2,na)=j*lp+lp*0.5
        R_at(3,na)=k*lp+lp*0.5
        na=na+1
        R_at(1,na)=i*lp+lp*0.5
        R_at(2,na)=j*lp
        R_at(3,na)=k*lp+lp*0.5
        na=na+1
        R_at(1,na)=i*lp+lp*0.5
        R_at(2,na)=j*lp+lp*0.5
        R_at(3,na)=k*lp
    enddo
    enddo
    enddo
    !print*, na, " atoms placed in FCC"
    nat = na
endsubroutine place_atoms_fcc
