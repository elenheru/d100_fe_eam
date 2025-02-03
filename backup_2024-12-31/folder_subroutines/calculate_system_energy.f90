subroutine calculate_system_energy(esys)
    use positions_mod
    use physical_parameters_mod
    use ackland4_parameters_mod, only : rphi_ack4pot
    implicit none
    real(8), intent(inout) :: esys
    integer i,j
    real(8) r,ec,eldc,pw_ack4pot,ed_ack4pot,mf_ack4pot,eat,cutoff
    esys = 0d0
    cutoff = rphi_ack4pot(15)
    do j=1,nat_mc
        eldc = 0d0 ; eat = 0d0
        do i=1,nat
            if(i .eq. j) cycle
            r=norm2(R_at(1:3,i) - R_at(1:3,j))
            if(r .gt. cutoff) cycle
            eat = eat + pw_ack4pot(r)
            eldc = eldc + ed_ack4pot(r)
        enddo
        esys = esys + eat + mf_ack4pot(eldc)
    enddo
endsubroutine calculate_system_energy
