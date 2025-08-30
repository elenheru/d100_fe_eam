subroutine calculate_force_atom_direct(number_atom, force)
    use electronic_densities_allocated_mod
    use potentials_ackland4fefe_linear_mod
    use positions_mod
    use system_parameters_mod
    use compile_debug_options_mod, only : debug_flag, debug_string_format

    real(float_byte_size), intent(inout), dimension(3) :: force
    integer(integer_byte_size), intent(in) :: number_atom

    real(float_byte_size), parameter :: shift_magnitude = 9e-3
    real(float_byte_size) :: distance
    real(float_byte_size), dimension(3) :: grad_num_t1, grad_num_t2, grad_num_t3, grad_num_t4, grad_num
    integer(integer_byte_size) i, j

    rho_at = 0e0

    grad_num = 0e0
    grad_num_t1 = 0e0
    grad_num_t2 = 0e0
    grad_num_t3 = 0e0
    grad_num_t4 = 0e0

    do i = 1, last_wall!_mc
        if (i .eq. number_atom) cycle
        distance = norm2( R_at(:, i) - R_at(:, number_atom) )
        grad_num_t1 = grad_num_t1 + pairwise_deriv1(distance) * ( R_at(:, i) - R_at(:, number_atom) ) / distance
    enddo

    do j = 1, last_wall
        if (j .eq. number_atom) cycle
        distance = norm2( R_at(:, j) - R_at(:, number_atom) )
        grad_num_t2 = grad_num_t2 + pairwise_deriv1(distance) * ( R_at(:, j) - R_at(:, number_atom) ) / distance
    enddo

    do j = 1, last_wall
        do i = 1, last_wall!_mc !this is the violation of my formula but I dont want it to be overcomplicated
            if (i .eq. j) cycle
            distance = norm2( R_at(:, j) - R_at(:, i) )
            rho_at(j) = rho_at(j) + edensity(distance)
        enddo
    enddo

    do j = 1, last_wall
        if (j .eq. number_atom) cycle
        distance = norm2( R_at(:, j) - R_at(:, number_atom) )
        grad_num_t3 = grad_num_t3 + edensity_deriv1(distance) * ( R_at(:, j) - R_at(:, number_atom) ) / distance
    enddo
    grad_num_t3 = grad_num_t3 * embeddin_deriv1(rho_at(number_atom))

    do i = 1, last_wall!_mc
        if (i .eq. number_atom) cycle
        distance = norm2( R_at(:, i) - R_at(:, number_atom) )
        grad_num_t4 = grad_num_t4 + &
        embeddin_deriv1(rho_at(i)) * edensity_deriv1(distance) * ( R_at(:, i) - R_at(:, number_atom) ) / distance
    enddo

    grad_num = - (grad_num_t1 + grad_num_t2 + grad_num_t3 + grad_num_t4)
    force = - grad_num

    if (debug_flag) then
        print 311, "Numerical g1", grad_num_t1
        print 311, "Numerical g2", grad_num_t2
        print 311, "Numerical g3", grad_num_t3
        print 311, "Numerical g4", grad_num_t4

        print"(A, 3(1x, es12.5), /, /)", "calculating gradient via analytical derivative formula: done", force
    endif

311 format(SP, A, 3(1x, F14.10))
endsubroutine calculate_force_atom_direct
