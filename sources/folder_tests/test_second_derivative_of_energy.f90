subroutine test_second_derivative_of_energy
    use potentials_ackland4fefe_linear_mod
    use positions_mod
    use system_parameters_mod
    use electronic_densities_allocated_mod
    use wo_storage_mod, only: format_m3x3
    use compile_debug_options_mod, only : &
        debug_flag, debug_string_format

    real(float_byte_size), parameter :: shift_magnitude = 0.01 * a0bcc
    real(float_byte_size), parameter :: derivative_denominator = 1e-4 * a0bcc

    integer(integer_byte_size) i, j, nprobe

    real(float_byte_size), dimension(3) :: shift
    real(float_byte_size), dimension(3) :: default_position

    real(float_byte_size), dimension(3) :: grad_anl_t1
    real(float_byte_size), dimension(3) :: grad_anl_t2
    real(float_byte_size), dimension(3) :: grad_anl_t3
    real(float_byte_size), dimension(3) :: grad_anl_t4
    real(float_byte_size), dimension(3) :: grad_anl

    real(float_byte_size), dimension(3, 3) :: hess_anl

    real(float_byte_size), dimension(3) :: grad_num_xp
    real(float_byte_size), dimension(3) :: grad_num_xm
    real(float_byte_size), dimension(3) :: grad_num_yp
    real(float_byte_size), dimension(3) :: grad_num_ym
    real(float_byte_size), dimension(3) :: grad_num_zp
    real(float_byte_size), dimension(3) :: grad_num_zm
    real(float_byte_size), dimension(3) :: grad_num
    real(float_byte_size), dimension(3, 3) :: hess_num

    real(float_byte_size) :: dist, rnum, d1, cosine_v
    real(float_byte_size) :: d1v_ik
    real(float_byte_size) :: d2v_ik
    real(float_byte_size) :: random_float
    real(float_byte_size) :: pw_ack4pot
    real(float_byte_size) :: ed_ack4pot
    real(float_byte_size) :: mf_ack4pot
    real(float_byte_size) :: d1pw_ack4pot
    real(float_byte_size) :: d1ed_ack4pot
    real(float_byte_size) :: d1mf_ack4pot

    print*, ""
    print*, "Testing formula for second derivative of system energy"
    print*, "This test assumes that gradient is already tested"
    print*, "Comparison between numerical and anlytical calculation"
    print*, "Probe atom has minimal energy position at origin"
    print*, ""


    call find_centeral_atom(nprobe)

    print *, "Taking atom # ", nprobe, " as probe"

    call random_number(shift)

    shift = ( 2e0*shift - (/1e0, 1e0, 1e0/) ) * shift_magnitude

    print 301, " Shift is ", shift

    default_position = R_at(:, nprobe) + shift
    R_at(:, nprobe) = default_position

    ! call system_clock(timers(1))
    STOP "Not implemented"
    print debug_string_format, &
    "calculating hessian via analytical formula"

    rho_at = 0e0
    grad_anl_t1 = 0e0
    grad_anl_t2 = 0e0
    grad_anl_t3 = 0e0
    grad_anl_t4 = 0e0

!    do i = 1, last_wall!_mc
!        if (i .eq. nprobe) cycle
!        dist = norm2( R_at(:, i) - R_at(:, nprobe) )
!        grad_anl_t1 = grad_anl_t1 + d1pw_ack4pot(dist) * ( R_at(:, i) - R_at(:, nprobe) ) / dist
!    enddo

    do j = 1, last_wall
        if (j .eq. nprobe) cycle
        dist = norm2( R_at(:, j) - R_at(:, nprobe) )
!        grad_anl_t2 = grad_anl_t2 + d1pw_ack4pot(dist) * ( R_at(:, j) - R_at(:, nprobe) ) / dist
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

!        grad_anl_t3 = &
!        grad_anl_t3 &
!        + d1ed_ack4pot(dist) &
!        * ( R_at(:, j) - R_at(:, nprobe) ) / dist
    enddo

    grad_anl_t3 = grad_anl_t3 * d1mf_ack4pot(rho_at(nprobe))

    do i = 1, last_wall!_mc

        if (i .eq. nprobe) cycle

        dist = norm2( R_at(:, i) - R_at(:, nprobe) )
!
!        grad_anl_t4 = grad_anl_t4 &
!        + d1mf_ack4pot(rho_at(i)) * d1ed_ack4pot(dist) &
!        * ( R_at(:, i) - R_at(:, nprobe) ) / dist

    enddo

    !call system_clock(timers(2))
!    print"(A, 3(1x, es12.5))", "Done ", grad_anl
!    print 311, "Analytical g1", grad_anl_t1
!    print 311, "Analytical g2", grad_anl_t2
!    print 311, "Analytical g3", grad_anl_t3
!    print 311, "Analytical g4", grad_anl_t4
    !print *, "Time spent for spline function ", timers(2) - timers(1)
!
!    cosine_v = dot_product(shift, grad_anl) / ( norm2(shift) * norm2(grad_anl) )
!    print*, "cosine is ", cosine_v
!
!    if (cosine_v .gt. -0.99) then
!        print*, "inappropriate cosine", cosine_v
!        print*, "number of atom, distance to probe, electronic density"
!        do i = 1, last_relaxable
!            dist = norm2( R_at(:, i) - R_at(:, nprobe) )
!            if (dist .gt. 6.5) cycle
!            print"(I7.2, SP, 2(1x, F15.10))", i, dist, rho_at(i)
!        enddo
!        do i = last_relaxable + 1, last_wall
!            dist = norm2( R_at(:, i) - R_at(:, nprobe) )
!            if (dist .gt. 6.5) cycle
!            print"(I7.2, SP, 2(1x, F15.10)), A", i, dist, rho_at(i), " elastic"
!        enddo
!        print*, r_at(:, nprobe), "distance is ", norm2(r_at(:, nprobe))
!    endif

    ! part
    !call system_clock(timers(3))
!    print debug_string_format, &
!    "calculating gradient via linearized functions derivative formula ... "

    rho_at = 0e0

    grad_num_xp = 0e0
    grad_num_xm = 0e0
    grad_num_yp = 0e0
    grad_num_ym = 0e0
    grad_num_zp = 0e0
    grad_num_zm = 0e0


    do i = 1, last_wall!_mc
        if (i .eq. nprobe) cycle
        dist = norm2( R_at(:, i) - R_at(:, nprobe) )
!        grad_num_t1 = grad_num_t1 + pairwise_deriv1(dist) * ( R_at(:, i) - R_at(:, nprobe) ) / dist
    enddo

    do j = 1, last_wall
        if (j .eq. nprobe) cycle
        dist = norm2( R_at(:, j) - R_at(:, nprobe) )
!        grad_num_t2 = grad_num_t2 + pairwise_deriv1(dist) * ( R_at(:, j) - R_at(:, nprobe) ) / dist
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
!        grad_num_t3 = grad_num_t3 + edensity_deriv1(dist) * ( R_at(:, j) - R_at(:, nprobe) ) / dist
    enddo
!    grad_num_t3 = grad_num_t3 * embeddin_deriv1(rho_at(nprobe))

    do i = 1, last_wall!_mc
        if (i .eq. nprobe) cycle
        dist = norm2( R_at(:, i) - R_at(:, nprobe) )
!        grad_num_t4 = grad_num_t4 + &
!        embeddin_deriv1(rho_at(i)) * edensity_deriv1(dist) * ( R_at(:, i) - R_at(:, nprobe) ) / dist
    enddo

301 format(A, 3(1x, F12.8))
311 format(SP, A, 3(1x, F14.10))

endsubroutine test_second_derivative_of_energy
