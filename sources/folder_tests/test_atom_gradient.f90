subroutine test_atom_gradient
    use potentials_ackland4fefe_linear_mod
    use positions_mod
    use system_parameters_mod
    use electronic_densities_allocated_mod
    use compile_debug_options_mod, only : &
        debug_flag, debug_string_format, timers


    real(float_byte_size), parameter :: shift_magnitude = 9d-3
    integer(integer_byte_size) i, j, nprobe
    real(float_byte_size) dist, rnum, d1v_ik, d1, rn, cosine_v
    real(float_byte_size), dimension(3) :: shift, move
    real(float_byte_size), dimension(3) :: grad_anl_t1, grad_anl_t2, grad_anl_t3, grad_anl_t4, grad_anl
    real(float_byte_size), dimension(3) :: grad_num_t1, grad_num_t2, grad_num_t3, grad_num_t4, grad_num
    real(float_byte_size) :: pw_ack4pot, ed_ack4pot, mf_ack4pot
    real(float_byte_size) :: d1pw_ack4pot, d1ed_ack4pot, d1mf_ack4pot

    print*, "Testing gradient formula for shift of a single atom "
    print*, "Also testing linearization in functions"

    call random_number(rn) ; nprobe = floor(rn * last_relaxable)
    print*, nprobe, "is atom probe number"

    print*, nprobe, "is probe atom"
    call random_number(shift)
    shift = ( 2d0*shift - (/1d0, 1d0, 1d0/) ) * shift_magnitude
    print 301, "shift is ", shift

    R_at(:, nprobe) = R_at(:, nprobe) + shift

    call system_clock(timers(1))

    print debug_string_format, "calculating gradient via analytical derivative formula ... "
    rho_at = 0d0
    grad_anl_t1 = 0d0
    grad_anl_t2 = 0d0
    grad_anl_t3 = 0d0
    grad_anl_t4 = 0d0

    do i = 1, last_wall!_mc
        if (i .eq. nprobe) cycle
        dist = norm2( R_at(:, i) - R_at(:, nprobe) )
        grad_anl_t1 = grad_anl_t1 + d1pw_ack4pot(dist) * ( R_at(:, i) - R_at(:, nprobe) ) / dist
    enddo

    do j = 1, last_wall
        if (j .eq. nprobe) cycle
        dist = norm2( R_at(:, j) - R_at(:, nprobe) )
        grad_anl_t2 = grad_anl_t2 + d1pw_ack4pot(dist) * ( R_at(:, j) - R_at(:, nprobe) ) / dist
    enddo

    do j = 1, last_wall
        do i = 1, last_wall!_mc !this is the violation of my formula but I dont want it to be overcomplicated
            if (i .eq. j) cycle
            dist = norm2( R_at(:, j) - R_at(:, i) )
            rho_at(j) = rho_at(j) + ed_ack4pot(dist)
        enddo
    enddo

    do j = 1, last_wall
        if (j .eq. nprobe) cycle
        dist = norm2( R_at(:, j) - R_at(:, nprobe) )
        grad_anl_t3 = grad_anl_t3 + d1ed_ack4pot(dist) * ( R_at(:, j) - R_at(:, nprobe) ) / dist
    enddo
    !print*, grad_anl_t3
    grad_anl_t3 = grad_anl_t3 * d1mf_ack4pot(rho_at(nprobe))

    do i = 1, last_wall!_mc
        if (i .eq. nprobe) cycle
        dist = norm2( R_at(:, i) - R_at(:, nprobe) )
        grad_anl_t4 = grad_anl_t4 + &
        d1mf_ack4pot(rho_at(i)) * d1ed_ack4pot(dist) * ( R_at(:, i) - R_at(:, nprobe) ) / dist
    enddo
    grad_anl = grad_anl_t1 + grad_anl_t2 + grad_anl_t3 + grad_anl_t4

    call system_clock(timers(2))
    print"(A, 3(1x, es12.5))", "Done ", grad_anl
    print 311, "Analytical g1", grad_anl_t1
    print 311, "Analytical g2", grad_anl_t2
    print 311, "Analytical g3", grad_anl_t3
    print 311, "Analytical g4", grad_anl_t4
    print *, "Time spent for spline function ", timers(2) - timers(1)

    cosine_v = dot_product(shift, grad_anl) / ( norm2(shift) * norm2(grad_anl) )
    print*, "cosine is ", cosine_v

    if (cosine_v .gt. -0.99) then
        print*, "inappropriate cosine", cosine_v
        print*, "number of atom, distance to probe, electronic density"
        do i = 1, last_relaxable
            dist = norm2( R_at(:, i) - R_at(:, nprobe) )
            if (dist .gt. 6.5) cycle
            print"(I7.2, SP, 2(1x, F15.10))", i, dist, rho_at(i)
        enddo
        do i = last_relaxable + 1, last_wall
            dist = norm2( R_at(:, i) - R_at(:, nprobe) )
            if (dist .gt. 6.5) cycle
            print"(I7.2, SP, 2(1x, F15.10)), A", i, dist, rho_at(i), " elastic"
        enddo
        print*, r_at(:, nprobe), "distance is ", norm2(r_at(:, nprobe))
    endif

    ! part
    call system_clock(timers(3))
    print debug_string_format, "calculating gradient via linearized functions derivative formula ... "
    rho_at = 0d0
    grad_num_t1 = 0d0
    grad_num_t2 = 0d0
    grad_num_t3 = 0d0
    grad_num_t4 = 0d0

    do i = 1, last_wall!_mc
        if (i .eq. nprobe) cycle
        dist = norm2( R_at(:, i) - R_at(:, nprobe) )
        grad_num_t1 = grad_num_t1 + pairwise_deriv1(dist) * ( R_at(:, i) - R_at(:, nprobe) ) / dist
    enddo

    do j = 1, last_wall
        if (j .eq. nprobe) cycle
        dist = norm2( R_at(:, j) - R_at(:, nprobe) )
        grad_num_t2 = grad_num_t2 + pairwise_deriv1(dist) * ( R_at(:, j) - R_at(:, nprobe) ) / dist
    enddo

    do j = 1, last_wall
        do i = 1, last_wall!_mc !this is the violation of my formula but I dont want it to be overcomplicated
            if (i .eq. j) cycle
            dist = norm2( R_at(:, j) - R_at(:, i) )
            rho_at(j) = rho_at(j) + edensity(dist)
        enddo
    enddo

    do j = 1, last_wall
        if (j .eq. nprobe) cycle
        dist = norm2( R_at(:, j) - R_at(:, nprobe) )
        grad_num_t3 = grad_num_t3 + edensity_deriv1(dist) * ( R_at(:, j) - R_at(:, nprobe) ) / dist
    enddo
    grad_num_t3 = grad_num_t3 * embeddin_deriv1(rho_at(nprobe))

    do i = 1, last_wall!_mc
        if (i .eq. nprobe) cycle
        dist = norm2( R_at(:, i) - R_at(:, nprobe) )
        grad_num_t4 = grad_num_t4 + &
        embeddin_deriv1(rho_at(i)) * edensity_deriv1(dist) * ( R_at(:, i) - R_at(:, nprobe) ) / dist
    enddo
    grad_num = grad_num_t1 + grad_num_t2 + grad_num_t3 + grad_num_t4

    call system_clock(timers(4))
    print"(A, 3(1x, es12.5), /, /)", " Done", grad_num
    print 311, "Numerical g1", grad_num_t1
    print 311, "Numerical g2", grad_num_t2
    print 311, "Numerical g3", grad_num_t3
    print 311, "Numerical g4", grad_num_t4
    print *, "Time spent for linear function ", timers(4) - timers(3)

    cosine_v = dot_product(shift, grad_num) / ( norm2(shift) * norm2(grad_num) )
    print*, "cosine is ", cosine_v

    if (cosine_v .gt. -0.99) then
        print*, "inappropriate cosine", cosine_v
        print*, "number of atom, distance to probe, electronic density"
        do i = 1, last_relaxable
            dist = norm2( R_at(:, i) - R_at(:, nprobe) )
            if (dist .gt. 6.5) cycle
            print"(I7.2, SP, 2(1x, F15.10))", i, dist, rho_at(i)
        enddo
        do i = last_relaxable + 1, last_wall
            dist = norm2( R_at(:, i) - R_at(:, nprobe) )
            if (dist .gt. 6.5) cycle
            print"(I7.2, SP, 2(1x, F15.10)), A", i, dist, rho_at(i), " elastic"
        enddo
        print*, r_at(:, nprobe), "distance is ", norm2(r_at(:, nprobe))
    endif


301 format(A, 3(1x, F12.8))
311 format(SP, A, 3(1x, F14.10))

endsubroutine test_atom_gradient
