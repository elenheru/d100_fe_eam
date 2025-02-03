subroutine calculate_energy_at(es,na)
    use positions_mod
    use physical_parameters_mod
    use ackland4_parameters_mod, only : rphi_ack4pot
    implicit none
    real(8), intent(inout) :: es
    integer, intent(in) :: na
    integer i
    real(8) r,ec,eldc,pw_ack4pot,ed_ack4pot,mf_ack4pot,cutoff
    eldc = 0d0 ; ec = 0d0
    cutoff = rphi_ack4pot(15)
    do i=1,nat
        if(i .eq. na) cycle
        r=norm2(R_at(1:3,i) - R_at(1:3,na))
        if(r .gt. cutoff) cycle
        ec = ec + pw_ack4pot(r)
        eldc = eldc + ed_ack4pot(r)
    enddo
    es = ec + mf_ack4pot(eldc)
endsubroutine calculate_energy_at
