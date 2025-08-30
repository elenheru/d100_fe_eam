subroutine test_field_work2
    use positions_mod
    use potentials_ackland4fefe_linear_mod
    use verlet_list_mod
    use electronic_densities_allocated_mod
    use compile_debug_options_mod, only : atom_of_interest

    implicit none
    real(float_byte_size), parameter :: shift_magnitude = 1e-2_float_byte_size
    integer(integer_byte_size) :: nprobe
    real(float_byte_size) :: float_random
    real(float_byte_size) :: energy_after_single_shift
    real(float_byte_size) :: energy_after_double_shift

    real(float_byte_size), dimension(3) :: force_1 = 0e0_float_byte_size
    real(float_byte_size), dimension(3) :: force_2 = 0e0_float_byte_size
    real(float_byte_size), dimension(3) :: force_average = 0e0_float_byte_size
    real(float_byte_size), dimension(3) :: shift

    print*, "Testing the fact that provided shift of atom is small "
    print*, "the dot product of force and shift (one atom) equals energy change"
    print*, "Shift is performed twice, because initial"
    print*, "atom position corresponds to strict minima."

    print*, "Modifying shift magnitude, beware the float cancellation."
    print*, "Cancellation may occur as far as system energies magnitudes"
    print*, "are too big compared to force work, and substraction is performed."

    !I am not sure about sign of this product. I dont care much though


    call random_number(float_random)

    nprobe = floor(float_random * last_relaxable)

    print*, nprobe, &
        "Is atom probe number. Shift it slightly in random direction", &
        norm2(R_at(:, nprobe))

    call random_number(shift)

    shift = &
        ( 2.0_float_byte_size * shift &
        - (/1e0_float_byte_size, 1e0_float_byte_size, 1e0_float_byte_size/) ) &
        * shift_magnitude

    print 301, "shift is ", shift

    R_at(:, nprobe) = R_at(:, nprobe) + shift

    call obtain_verlet_list

    call calculate_system_energy_verlet(energy_after_single_shift)

    print *, "System energy after 1 shift is: ", energy_after_single_shift

    call calculate_force_atom_verlet(nprobe, force_1)

    R_at(:, nprobe) = R_at(:, nprobe) + shift

    call obtain_verlet_list

    call calculate_system_energy_verlet(energy_after_double_shift)

    print *, "System energy after 2 shift is: ", energy_after_double_shift

    call calculate_force_atom_verlet(nprobe, force_2)

    force_average = 0.5 * (force_1 + force_2)
    print 301, "Force (average) is ", force_average
    print *, "Energy difference   is ", &
        energy_after_double_shift - energy_after_single_shift
    print *, "Work of field force is ", dot_product(shift, force_average)

    ! Division by two is correct, because force is not constant along shift.
    ! It is linearly increasing as shfit increases,
    ! provided minima is quadratic.
    ! Long story short - we are near to minima, so linear
    ! approximation is guaranteed to be wrong


    call wo_xyz_snapshot(nprobe)
    R_at(:, nprobe) = R_at(:, nprobe) - shift
    call wo_xyz_snapshot(nprobe)

301 format(SP, A, 3(1x, F16.12))

endsubroutine test_field_work2
