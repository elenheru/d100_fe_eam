subroutine test_field_work3
    use positions_mod
    use potentials_ackland4fefe_linear_mod
    use verlet_list_mod
    use electronic_densities_allocated_mod
    use compile_debug_options_mod, only : atom_of_interest
    use system_parameters_mod, only : a0bcc

    implicit none
    real(float_byte_size), parameter :: shift_magnitude = &
    a0bcc * 1e-3_float_byte_size
    integer(integer_byte_size) :: i_atom
    real(float_byte_size) :: float_random
    real(float_byte_size) :: energy_after_single_shift
    real(float_byte_size) :: energy_after_double_shift
    real(float_byte_size) :: work_of_force

    real(float_byte_size), dimension(:,:), allocatable :: shift
    real(float_byte_size), dimension(:,:), allocatable  :: force_1
    real(float_byte_size), dimension(:,:), allocatable  :: force_2
    real(float_byte_size), dimension(3) :: force_temporary
    real(float_byte_size), dimension(3) :: atom_shift
    real(float_byte_size), dimension(3), parameter :: &
        unit_vector = (/ &
        1e0_float_byte_size, &
        1e0_float_byte_size, &
        1e0_float_byte_size /)

    print*, " "
    print*, "Testing the fact that provided all shifts of atoms are small "
    print*, "the dot product of force and shift equals energy change"
    print*, "Shifts are performed twice, because initial"
    print*, "atomic positions correspond to strict minima."
    print*, "This test is for whole relaxable part of system"

    print*, " "
    print*, "!Modifying shift magnitude, beware the float cancellation."
    print*, "!Cancellation may occur as far as system energies magnitudes"
    print*, "!are too big compared to force work, and substraction is performed"

    print 302, " Shift magnitude is limited with", shift_magnitude, " angstroms"
    print*, " "


    !I am not sure about sign of this product. I dont care much though

    allocate(shift(3, last_relaxable))
    allocate(force_1(3, last_relaxable))
    allocate(force_2(3, last_relaxable))

    ! shifing all atoms from perfect lattice

    call random_number(shift)

    do i_atom = 1, last_relaxable

        shift(:, i_atom) = &
        ( 2.0_float_byte_size * shift(:, i_atom) - unit_vector ) &
        * shift_magnitude

        R_at(:, i_atom) = R_at(:, i_atom) + shift(:, i_atom)

    enddo

    ! print 301, "shift is ", shift


    call obtain_verlet_list

    call calculate_system_energy_verlet(energy_after_single_shift)

    print *, "System energy after 1 shift is: ", energy_after_single_shift

    do i_atom = 1, last_relaxable

        call calculate_force_atom_verlet(i_atom, force_temporary)

        force_1(:, i_atom) = force_temporary
        ! this may be too cautious, but i dont want to dive into array slice passing
    enddo

    do i_atom = 1, last_relaxable

        R_at(:, i_atom) = R_at(:, i_atom) + shift(:, i_atom)

    enddo

    call obtain_verlet_list

    call calculate_system_energy_verlet(energy_after_double_shift)

    print *, "System energy after 2 shift is: ", energy_after_double_shift

    do i_atom = 1, last_relaxable

        call calculate_force_atom_verlet(i_atom, force_temporary)

        force_2(:, i_atom) = force_temporary
        ! this may be too cautious, but i dont want to dive into array slice passing
    enddo

    print *, "Energy difference   is ", &
        energy_after_double_shift - energy_after_single_shift

    work_of_force = 0e0_float_byte_size

    do i_atom = 1, last_relaxable

        force_temporary = 0.5_float_byte_size &
        * (force_1(:, i_atom)  + force_2(:, i_atom) )

        atom_shift = shift(:, i_atom)

        work_of_force = work_of_force + dot_product(atom_shift, force_temporary)

    enddo

    print *, "Work of field force is ", work_of_force

    call wo_xyz_snapshot(atom_of_interest)

301 format(SP, A, 3(1x, F16.12))
302 format(A, 1x, F16.12, A)

endsubroutine test_field_work3
