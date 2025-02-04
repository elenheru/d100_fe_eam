subroutine calculate_pseudosystem_energy(esys, na)
    use positions_mod
    use physical_parameters_mod
    use ackland4_parameters_mod, only : rphi_ack4pot
    implicit none
    real(float_byte_size), intent(inout) :: esys
    integer(integer_byte_size), intent(in) :: na
    integer(integer_byte_size) i, j
    real(float_byte_size) r, ec, eldc, pw_ack4pot, ed_ack4pot, mf_ack4pot, eat, cutoff
    esys = 0d0
    cutoff = rphi_ack4pot(15)
    do j=1, last_relaxable
        if (norm2(R_at(1:3, na) - R_at(1:3, j)) .gt. cutoff) cycle
        eldc = 0d0 ; eat = 0d0
        do i=1, last_wall
            if(i .eq. j) cycle
            r=norm2(R_at(1:3, i) - R_at(1:3, j))
            if(r .gt. cutoff) cycle
            eat = eat + pw_ack4pot(r)
            eldc = eldc + ed_ack4pot(r)
        enddo
        esys = esys + eat + mf_ack4pot(eldc)
    enddo
endsubroutine calculate_pseudosystem_energy
