subroutine relax_fine(ncent)
    use positions_mod
    implicit none
    integer i,j
    integer, intent(in) :: ncent
    real(8) :: esystem
    STOP " not needed relax fine"
    do j=1,8
        call calculate_system_energy(esystem)
        print*, esystem, " eV is system energy, fine relax"
        do i=1,nat_mc
            if (i .eq. ncent) cycle
            !call relax_fine_quad_grad_rand(i,6d-3**(j))
        enddo
        call wo_xyz_snapshot
    enddo
    call calculate_system_energy(esystem)
    print*, esystem, " eV is system energy, fine relax."
endsubroutine relax_fine
