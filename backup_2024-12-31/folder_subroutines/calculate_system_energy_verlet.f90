subroutine calculate_system_energy_verlet(esys)
    use positions_mod
    use physical_parameters_mod
    use ackland4_parameters_mod, only : rphi_ack4pot
    use miscellaneous_parameters_mod, only : debug
    use verlet_list_mod
    implicit none
    real(8), intent(inout) :: esys
    integer i,j,iv,jv
    real(8) r,ec,eldc,pw_ack4pot,ed_ack4pot,mf_ack4pot,eat,cutoff,dist
    real(8) grad_num, grad_num_t1, grad_num_t2, grad_num_t3, grad_num_t4
    real(8), allocatable :: rho_at(:)


    if (debug) print*, "calculate_system_energy_verlet does not automaticly renew Verlet list"
    STOP "Not implemented"
    esys = 0d0
    cutoff = rphi_ack4pot(15)
    do j = 1,nat_mc
        eldc = 0d0 ; eat = 0d0
        do i = 1, verlet_list_length(j)
            iv = verlet_list(j,i)
            if(i .eq. j) STOP "own atom in verlet list"
            r = norm2(R_at(1:3,iv) - R_at(1:3,j))
            !if(r .gt. cutoff) cycle
            eat = eat + pw_ack4pot(r)
            eldc = eldc + ed_ack4pot(r)
        enddo
        esys = esys + eat + mf_ack4pot(eldc)
    enddo

! below is the copypast from test_verlet_list
! when finished, procedure shall calculate energy using verlet list, assuming it is already updated

! RHO_AT likely to be better moved in a module
! it is needed to know exactly when it shall be nullified

endsubroutine calculate_system_energy_verlet
