subroutine test_namgd
    use potentials_ackland4fefe_linear_mod
    use positions_mod
    use system_parameters_mod
    use electronic_densities_allocated_mod
    use wo_storage_mod, only: format_m3x3
    use compile_debug_options_mod, only : &
        debug_flag, debug_string_format


    real(float_byte_size), parameter :: shift_magnitude = 0.45 * a0bcc
    real(float_byte_size), parameter :: momentum_coefficient = 2e-1
    real(float_byte_size), parameter :: force_coefficient = 9e-1
    real(float_byte_size), parameter :: mass_reciprocal = 4e-2
    real(float_byte_size), parameter :: common_step_denominator_big = 0.85
    real(float_byte_size), parameter :: common_step_denominator = 1.2

    integer, parameter :: step_number = 160

    integer(integer_byte_size) i, nprobe

    real(float_byte_size) :: common_step_magnitude
    real(float_byte_size) :: common_step_coefficient
    real(float_byte_size) :: random_multiplicator
    real(float_byte_size), dimension(3) :: shift
    real(float_byte_size), dimension(3) :: forestep

    real(float_byte_size), dimension(3) :: previous_position
    real(float_byte_size), dimension(3) :: current_position
    real(float_byte_size), dimension(3) :: forehead_position
    real(float_byte_size), dimension(3) :: next_best_position
    real(float_byte_size), dimension(3) :: step_calculated

    real(float_byte_size), dimension(3) :: force_current
    real(float_byte_size), dimension(3) :: force
    real(float_byte_size), dimension(3) :: force_prev

    real(float_byte_size) :: field_half_work_1
    real(float_byte_size) :: field_half_work_2


    print *, " "
    print *, "Test of Nesterov accelerated momentum gradient descend"
    print *, "But this is not the test of procedure, but idea"
    print *, " "

    common_step_magnitude = 0.006 * a0bcc
    common_step_coefficient = 1
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

    call random_number(random_multiplicator)
    call calculate_force_atom_verlet(nprobe, force_prev)
    print 301, " Force initial is ", force_prev
    print 301, " -------------------------------- "

    do i = 1, step_number
        call renew_verlet_list(nprobe)
        forestep = momentum_coefficient &
        * (current_position - previous_position)

        forehead_position = current_position + forestep

    call calculate_force_atom_verlet(nprobe, force_current)

        R_at(:, nprobe) = forehead_position
        call renew_verlet_list(nprobe)

        call calculate_force_atom_verlet(nprobe, force)

        step_calculated = ( forestep &
        + force_coefficient * force * mass_reciprocal)

        step_calculated = common_step_magnitude * step_calculated &
        / norm2(step_calculated)

        next_best_position = current_position &
        + step_calculated * random_multiplicator

        R_at(:, nprobe) = next_best_position

        call renew_verlet_list(nprobe)
    call calculate_force_atom_verlet(nprobe, force)
    field_half_work_1 = dot_product(force_current, step_calculated) * 0.5
    field_half_work_2 = dot_product(force, step_calculated) * 0.5

        if ( (field_half_work_1 + field_half_work_2) .lt. 0 ) then
            print '(A,$)', 'A'
            R_at(:, nprobe) = current_position
            common_step_magnitude = &
            common_step_magnitude * common_step_denominator_big
            call random_number(random_multiplicator)
        else
            print '(A,$)', 'V'
            previous_position = current_position
            current_position = R_at(:, nprobe)
            !random_multiplicator = (random_multiplicator + 0.1) * 1.1
            common_step_magnitude = &
            common_step_magnitude * common_step_denominator
        endif

        print 303, " dist ", norm2(R_at(:, nprobe))
        print 303, " step ", norm2(step_calculated)
        print 302, " Position now ", R_at(:, nprobe)
        print 301, " Force forehead ", force
        !print*, "Work of field is ", field_half_work_1 + field_half_work_2

    enddo

301 format(A, 3(1x, F9.6)   )
302 format(A, 3(1x, F9.6), $)
303 format(A,  (1x, F9.6), $)

endsubroutine test_namgd
