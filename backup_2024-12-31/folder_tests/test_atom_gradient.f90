subroutine test_atom_gradient
    use potentials_ackland4fefe_linear_mod
    use positions_mod
    use system_parameters_mod
    implicit none
    real(8), parameter :: shift_magnitude = 9d-3
    integer i,j,nprobe
    real(8) dist, rnum, d1v_ik, d1, rn, cosine_v
    real(8), allocatable :: rho_at(:)
    real(8), dimension(3) :: shift, move
    real(8), dimension(3) :: grad_anl_t1, grad_anl_t2, grad_anl_t3, grad_anl_t4, grad_anl
    real(8), dimension(3) :: grad_num_t1, grad_num_t2, grad_num_t3, grad_num_t4, grad_num
    real(8) :: pw_ack4pot, ed_ack4pot, mf_ack4pot
    real(8) :: d1pw_ack4pot, d1ed_ack4pot, d1mf_ack4pot

    print*, "testing gradient formula for shift of a single atom "
    call place_atoms_bcc_spheric

    call random_number(rn) ; nprobe = floor(rn * nat_mc)
    print*, nprobe, "is atom probe number"

    print*, nprobe, "is probe atom"
    call random_number(shift)
    shift = ( 2d0*shift - (/1d0,1d0,1d0/) ) * shift_magnitude
    print 301, "shift is ", shift

    R_at(:,nprobe) = R_at(:,nprobe) + shift
    allocate(rho_at(nat))
    rho_at = 0d0
    grad_anl_t1 = 0d0
    grad_anl_t2 = 0d0
    grad_anl_t3 = 0d0
    grad_anl_t4 = 0d0
    print"(A)", "calculating gradient via analytical derivative formula"

    do i = 1, nat!_mc
        if (i .eq. nprobe) cycle
        dist = norm2( R_at(:,i) - R_at(:,nprobe) )
        grad_anl_t1 = grad_anl_t1 + d1pw_ack4pot(dist) * ( R_at(:,i) - R_at(:,nprobe) ) / dist
    enddo

    do j = 1, nat
        if (j .eq. nprobe) cycle
        dist = norm2( R_at(:,j) - R_at(:,nprobe) )
        grad_anl_t2 = grad_anl_t2 + d1pw_ack4pot(dist) * ( R_at(:,j) - R_at(:,nprobe) ) / dist
    enddo

    do j = 1, nat
        do i = 1, nat!_mc !this is the violation of my formula but I dont want it to be overcomplicated
            if (i .eq. j) cycle
            dist = norm2( R_at(:,j) - R_at(:,i) )
            rho_at(j) = rho_at(j) + ed_ack4pot(dist)
        enddo
    enddo

    do j = 1, nat
        if (j .eq. nprobe) cycle
        dist = norm2( R_at(:,j) - R_at(:,nprobe) )
        grad_anl_t3 = grad_anl_t3 + d1ed_ack4pot(dist) * ( R_at(:,j) - R_at(:,nprobe) ) / dist
    enddo
    !print*, grad_anl_t3
    grad_anl_t3 = grad_anl_t3 * d1mf_ack4pot(rho_at(nprobe))

    do i = 1, nat!_mc
        if (i .eq. nprobe) cycle
        dist = norm2( R_at(:,i) - R_at(:,nprobe) )
        grad_anl_t4 = grad_anl_t4 + &
        d1mf_ack4pot(rho_at(i)) * d1ed_ack4pot(dist) * ( R_at(:,i) - R_at(:,nprobe) ) / dist
    enddo
    grad_anl = grad_anl_t1 + grad_anl_t2 + grad_anl_t3 + grad_anl_t4
    print 311, "g1", grad_anl_t1
    print 311, "g2", grad_anl_t2
    print 311, "g3", grad_anl_t3
    print 311, "g4", grad_anl_t4
    print"(A,3(1x,es12.5),/,/)", "calculating gradient via analytical derivative formula: done",grad_anl
    cosine_v = dot_product(shift,grad_anl) / ( norm2(shift) * norm2(grad_anl) )
    print*, "cosine is ", cosine_v

    if (cosine_v .gt. -0.99) then
        print*, "inappropriate cosine", cosine_v
        print*, "number of atom, distance to probe, electronic density"
        do i = 1, nat_mc
            dist = norm2( R_at(:,i) - R_at(:,nprobe) )
            if (dist .gt. 6.5) cycle
            print"(I7.2, SP, 2(1x, F15.10))", i, dist, rho_at(i)
        enddo
        do i = nat_mc + 1, nat
            dist = norm2( R_at(:,i) - R_at(:,nprobe) )
            if (dist .gt. 6.5) cycle
            print"(I7.2, SP, 2(1x, F15.10)), A", i, dist, rho_at(i), " elastic"
        enddo
        print*, r_at(:,nprobe), "distance is ", norm2(r_at(:,nprobe))
    endif
    deallocate(rho_at)
    print"(A)", "calculating gradient via numerical derivative formula"
    print"(A,/,/)", "calculating gradient via numerical derivative formula: done"

301 format(A,3(1x,F12.8))
311 format(SP, A,3(1x,F14.10))

endsubroutine test_atom_gradient
