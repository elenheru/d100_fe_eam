subroutine test_far_atom_gradient
    use potentials_ackland4fefe_linear_mod
    use positions_mod
    use system_parameters_mod
    implicit none
    real(float_byte_size), parameter :: shift_magnitude = 3d-3
    integer(integer_byte_size) i, j, nfar, ncent
    real(float_byte_size) dist, rnum, d1v_ik, d1
    real(float_byte_size), allocatable :: rho_at(:)
    real(float_byte_size), dimension(3) :: shift, move
    real(float_byte_size), dimension(3) :: grad_anl_t1, grad_anl_t2, grad_anl_t3, grad_anl_t4, grad_anl
    real(float_byte_size), dimension(3) :: grad_num_t1, grad_num_t2, grad_num_t3, grad_num_t4, grad_num
    real(float_byte_size) :: pw_ack4pot, ed_ack4pot, mf_ack4pot
    real(float_byte_size) :: d1pw_ack4pot, d1ed_ack4pot, d1mf_ack4pot

    print*, "testing gradient formula for shift of a single atom "
    call find_far_atom(nfar)

    print*, nfar, "is atom far"

    call random_number(shift)
    shift = ( 2d0*shift - (/1d0, 1d0, 1d0/) ) * shift_magnitude
    print 301, "shift is ", shift

    R_at(:, nfar) = R_at(:, nfar) + shift
    allocate(rho_at(last_wall))
    rho_at = 0d0
    grad_anl_t1 = 0d0
    grad_anl_t2 = 0d0
    grad_anl_t3 = 0d0
    grad_anl_t4 = 0d0
    print"(A)", "calculating gradient via analytical derivative formula"

    do i = 1, last_wall!_mc
        if (i .eq. nfar) cycle
        dist = norm2( R_at(:, i) - R_at(:, nfar) )
        grad_anl_t1 = grad_anl_t1 + d1pw_ack4pot(dist) * ( R_at(:, i) - R_at(:, nfar) ) / dist
    enddo

    do j = 1, last_wall
        if (j .eq. nfar) cycle
        dist = norm2( R_at(:, j) - R_at(:, nfar) )
        grad_anl_t2 = grad_anl_t2 + d1pw_ack4pot(dist) * ( R_at(:, j) - R_at(:, nfar) ) / dist
    enddo

    do j = 1, last_wall
        do i = 1, last_wall!_mc !this is the violation of my formula but I dont want it to be overcomplicated
            if (i .eq. j) cycle
            dist = norm2( R_at(:, j) - R_at(:, i) )
            rho_at(j) = rho_at(j) + ed_ack4pot(dist)
        enddo
    enddo

    do j = 1, last_wall
        if (j .eq. nfar) cycle
        dist = norm2( R_at(:, j) - R_at(:, nfar) )
        grad_anl_t3 = grad_anl_t3 + d1ed_ack4pot(dist) * ( R_at(:, j) - R_at(:, nfar) ) / dist
    enddo
    !print*, grad_anl_t3
    grad_anl_t3 = grad_anl_t3 * d1mf_ack4pot(rho_at(nfar))

    do i = 1, last_wall!_mc
        if (i .eq. nfar) cycle
        dist = norm2( R_at(:, i) - R_at(:, nfar) )
        grad_anl_t4 = grad_anl_t4 + &
        d1mf_ack4pot(rho_at(i)) * d1ed_ack4pot(dist) * ( R_at(:, i) - R_at(:, nfar) ) / dist
    enddo
    grad_anl = grad_anl_t1 + grad_anl_t2 + grad_anl_t3 + grad_anl_t4
    print 311, "g1", grad_anl_t1
    print 311, "g2", grad_anl_t2
    print 311, "g3", grad_anl_t3
    print 311, "g4", grad_anl_t4
    print *, "cosine is ", dot_product(shift, grad_anl) / ( norm2(shift) * norm2(grad_anl) )
    print"(A, 3(1x, es12.5), /, /)", "calculating gradient via analytical derivative formula: done", grad_anl

    print"(A)", "calculating gradient via numerical derivative formula"
    rho_at = 0d0
    grad_num_t1 = 0d0
    grad_num_t2 = 0d0
    grad_num_t3 = 0d0
    grad_num_t4 = 0d0
    do i = 1, last_wall!_mc
        if (i .eq. nfar) cycle
        dist = norm2( R_at(:, i) - R_at(:, nfar) )
        grad_num_t1 = grad_num_t1 + pairwise_deriv1(dist) * ( R_at(:, i) - R_at(:, nfar) ) / dist
    enddo

    do j = 1, last_wall
        if (j .eq. nfar) cycle
        dist = norm2( R_at(:, j) - R_at(:, nfar) )
        grad_num_t2 = grad_num_t2 + pairwise_deriv1(dist) * ( R_at(:, j) - R_at(:, nfar) ) / dist
    enddo

    do j = 1, last_wall
        do i = 1, last_wall!_mc !this is the violation of my formula but I dont want it to be overcomplicated
            if (i .eq. j) cycle
            dist = norm2( R_at(:, j) - R_at(:, i) )
            rho_at(j) = rho_at(j) + edensity(dist)
        enddo
    enddo

    do j = 1, last_wall
        if (j .eq. nfar) cycle
        dist = norm2( R_at(:, j) - R_at(:, nfar) )
        grad_num_t3 = grad_num_t3 + edensity_deriv1(dist) * ( R_at(:, j) - R_at(:, nfar) ) / dist
    enddo
    !print*, grad_num_t3
    grad_num_t3 = grad_num_t3 * embeddin_deriv1(rho_at(nfar))

    do i = 1, last_wall!_mc
        if (i .eq. nfar) cycle
        dist = norm2( R_at(:, i) - R_at(:, nfar) )
        grad_num_t4 = grad_num_t4 + &
        embeddin_deriv1(rho_at(i)) * edensity_deriv1(dist) * ( R_at(:, i) - R_at(:, nfar) ) / dist
    enddo
    grad_num = grad_num_t1 + grad_num_t2 + grad_num_t3 + grad_num_t4
    print 311, "g1", grad_num_t1
    print 311, "g2", grad_num_t2
    print 311, "g3", grad_num_t3
    print 311, "g4", grad_num_t4
    print *, "cosine is ", dot_product(shift, grad_num) / ( norm2(shift) * norm2(grad_num) )
    print"(A, 3(1x, es12.5), /, /)", "calculating gradient via numerical derivative formula: done", grad_num

    deallocate(rho_at)


301 format(A, 3(1x, F12.8))
311 format(SP, A, 3(1x, F14.10))

endsubroutine test_far_atom_gradient
