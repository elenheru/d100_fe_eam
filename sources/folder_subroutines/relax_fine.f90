subroutine relax_fine(ncent)
    use positions_mod
    implicit none
    integer(integer_byte_size) i, j
    integer(integer_byte_size), intent(in) :: ncent
    real(float_byte_size) :: esystem
    STOP " not needed relax fine"
    do j=1, 8
        call calculate_system_energy(esystem)
        print*, esystem, " eV is system energy, fine relax"
        do i=1, last_relaxable
            if (i .eq. ncent) cycle
            !call relax_fine_quad_grad_rand(i, 6d-3**(j))
        enddo
        call wo_xyz_snapshot
    enddo
    call calculate_system_energy(esystem)
    print*, esystem, " eV is system energy, fine relax."
endsubroutine relax_fine
