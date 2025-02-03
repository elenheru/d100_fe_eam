subroutine test_field_work
    use positions_mod
    use potentials_ackland4fefe_linear_mod
    use verlet_list_mod
    use electronic_densities_allocated_mod
    use compile_debug_options_mod, only : atom_of_interest

    implicit none
    real(float_byte_size), parameter :: shift_magnitude = 1d-1
    integer(integer_byte_size) :: nprobe
    real(float_byte_size) :: float_random
    real(float_byte_size) :: energy_before_shift
    real(float_byte_size) :: energy_after_shift

    real(float_byte_size), dimension(3) :: force
    real(float_byte_size), dimension(3) :: shift

    print*, "Testing the fact that provided shift of atom is small "
    print*, "the dot product of force and shift (one atom) equals energy change"
    !I am not sure about sign of this product. I dont care much though


    call calculate_system_energy_verlet(energy_before_shift)
!    call calculate_system_energy(energy_before_shift)
    print *, "System energy before shift is: ", energy_before_shift

    call random_number(float_random)
    nprobe = floor(float_random * last_relaxable)
    nprobe = atom_of_interest
    !nprobe = 370
    print*, nprobe, &
        "Is atom probe number. Shift it slightly in random direction", &
        norm2(R_at(:, nprobe))

    call random_number(shift)
    shift = ( 2d0*shift - (/1d0, 1d0, 1d0/) ) * shift_magnitude
    shift = ( 2d0*(/0.3, 0.7, 0.02/) - (/1d0, 1d0, 1d0/) ) * shift_magnitude
    print 301, "shift is ", shift
    call inspection_system_energy_wrong_change
    R_at(:, nprobe) = R_at(:, nprobe) + shift
    call obtain_verlet_list
    !because of two outputs
    call calculate_system_energy_verlet(energy_after_shift)
!    call calculate_system_energy(energy_after_shift)
    print *, "System energy after shift is: ", energy_after_shift

    call calculate_force_atom_verlet(nprobe, force)
!    call calculate_force_atom_direct(nprobe, force)

    print 301, "Force is ", force
    print *, "Energy difference   is ", energy_after_shift - energy_before_shift
    print *, "Work of field force is ", dot_product(shift, force)

    call wo_xyz_snapshot(nprobe)
    R_at(:, nprobe) = R_at(:, nprobe) - shift
    call wo_xyz_snapshot(nprobe)
!    print *, "Time spent for direct energy calculation ", timers(2) - timers(1)
!    print *, "Time spent for verlet energy calculation ", timers(4) - timers(3)
301 format(SP, A, 3(1x, F16.12))
endsubroutine test_field_work
