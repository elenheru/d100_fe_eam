subroutine relax_rude(ncent)
    use positions_mod
    implicit none
    integer i,j
    integer, intent(in) :: ncent
    real(8) :: esystem
    STOP " not needed relax rude"
    do j=1,26*7
        call calculate_system_energy(esystem)
        print*, esystem, " eV is system energy, rude relax"
        do i=1,nat_mc
            if (i .eq. ncent) cycle
            !call relax_rude_quad_grad_rand(i,6d-2**(j/3))
        enddo
        call wo_xyz_snapshot
    enddo
    call calculate_system_energy(esystem)
    print*, esystem, " eV is system energy, rude relax."
endsubroutine relax_rude
