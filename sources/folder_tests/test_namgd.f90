subroutine test_namgd
    use potentials_ackland4fefe_linear_mod
    use positions_mod
    use system_parameters_mod
    use electronic_densities_allocated_mod
    use wo_storage_mod, only: format_m3x3
    use compile_debug_options_mod, only : &
        debug_flag, debug_string_format


    real(float_byte_size), parameter :: shift_magnitude = 0.01 * a0bcc
    real(float_byte_size), parameter :: momentum_coefficient = 6e-1
    real(float_byte_size), parameter :: force_coefficient = 5e-1
    real(float_byte_size), parameter :: mass_reciprocal = 8e-2
    real(float_byte_size), parameter :: common_step_denominator = 0.9

    integer, parameter :: step_number = 50

    integer(integer_byte_size) i, nprobe

    real(float_byte_size) :: common_step_magnitude
    real(float_byte_size), dimension(3) :: shift
    real(float_byte_size), dimension(3) :: forestep

    real(float_byte_size), dimension(3) :: previous_position
    real(float_byte_size), dimension(3) :: current_position
    real(float_byte_size), dimension(3) :: forehead_position
    real(float_byte_size), dimension(3) :: next_best_position
    real(float_byte_size), dimension(3) :: step_calculated

    real(float_byte_size), dimension(3) :: force
    real(float_byte_size), dimension(3) :: force_prev


    print *, " "
    print *, "Test of Nesterov accelerated momentum gradient descend"
    print *, " "

    common_step_magnitude = 0.05 * a0bcc
    force = 0d0
    force_prev = 0d0

    call find_centeral_atom(nprobe)

    print *, "Taking atom # ", nprobe, " as probe"

    call random_number(shift)

    shift = ( 2e0*shift - (/1e0, 1e0, 1e0/) ) * shift_magnitude

    print 301, " Shift is ", shift

    R_at(:, nprobe) = R_at(:, nprobe) + shift

    previous_position = R_at(:, nprobe)
    current_position = R_at(:, nprobe)

    do i = 1, step_number
        forestep = momentum_coefficient &
        * (current_position - previous_position)

        forehead_position = current_position + forestep

        R_at(:, nprobe) = forehead_position

        call calculate_force_atom_verlet(nprobe, force)

        step_calculated = ( forestep &
        - force_coefficient * force * mass_reciprocal)

        step_calculated = step_calculated / norm2(step_calculated)

        next_best_position = current_position &
        + step_calculated * common_step_magnitude


        if ( norm2(force) .lt. norm2(force_prev) ) then

            R_at(:, nprobe) = next_best_position
            previous_position = current_position
            current_position = R_at(:, nprobe)
        else
            R_at(:, nprobe) = current_position
            common_step_magnitude = &
            common_step_magnitude * common_step_denominator
!            common_step_magnitude = &
!            common_step_magnitude * 1.1
        endif
        force_prev = force

        print 302, " Position ", R_at(:, nprobe)
        print 301, " Force ", - force

    enddo

301 format(A, 3(1x, F12.8))
302 format(A, 3(1x, F12.8), $)

endsubroutine test_namgd
