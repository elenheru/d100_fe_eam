subroutine test_verlet_list
    use positions_mod
    use potentials_ackland4fefe_linear_mod
    use verlet_list_mod

    implicit none
    real(8), parameter :: shift_magnitude = 1d-3
    integer i, j, iv, jv, nprobe
    real(8) dist, rn, cosine_v
    real(8), allocatable :: rho_at(:)
    real(8), dimension(3) :: grad_num_t1, grad_num_t2, grad_num_t3, grad_num_t4, grad_num
    real(8), dimension(3) :: shift, move

    print*, "testing verlet list for atoms "

    print"(A)", "calculating force via numerical derivative formula, using verlet list"

    call random_number(rn) ; nprobe = floor(rn * nat_mc)
    print*, nprobe, "is atom probe number, which we are going to shift lightly in random direction"

    call random_number(shift)
    shift = ( 2d0*shift - (/1d0,1d0,1d0/) ) * shift_magnitude
    print 301, "shift is ", shift

    R_at(:,nprobe) = R_at(:,nprobe) + shift

    allocate(rho_at(nat))
    rho_at = 0d0
    grad_num_t1 = 0d0
    grad_num_t2 = 0d0
    grad_num_t3 = 0d0
    grad_num_t4 = 0d0


    do i = 1, verlet_list_length(nprobe)
        iv = verlet_list(nprobe,i)
        if (iv .eq. nprobe) cycle
        dist = norm2( R_at(:,iv) - R_at(:,nprobe) )
        grad_num_t1 = grad_num_t1 + pairwise_deriv1(dist) * ( R_at(:,iv) - R_at(:,nprobe) ) / dist
    enddo

    do j = 1, verlet_list_length(nprobe)
        jv = verlet_list(nprobe,j)
        if (jv .eq. nprobe) cycle
        dist = norm2( R_at(:,jv) - R_at(:,nprobe) )
        grad_num_t2 = grad_num_t2 + pairwise_deriv1(dist) * ( R_at(:,jv) - R_at(:,nprobe) ) / dist
    enddo

    do j = 1, verlet_list_length(nprobe)
        jv = verlet_list(nprobe,j)
        do i = 1, verlet_list_length(jv)
            iv = verlet_list(jv,i)
            if (iv .eq. jv) cycle
            dist = norm2( R_at(:,jv) - R_at(:,iv) )
            rho_at(jv) = rho_at(jv) + edensity(dist)
        enddo
    enddo
    jv = nprobe
    do i = 1, verlet_list_length(jv)
        iv = verlet_list(jv,i)
        if (iv .eq. jv) cycle
        dist = norm2( R_at(:,jv) - R_at(:,iv) )
        rho_at(jv) = rho_at(jv) + edensity(dist)
    enddo

    do j = 1, verlet_list_length(nprobe)
        jv = verlet_list(nprobe,j)
        if (jv .eq. nprobe) cycle
        dist = norm2( R_at(:,jv) - R_at(:,nprobe) )
        grad_num_t3 = grad_num_t3 + edensity_deriv1(dist) * ( R_at(:,jv) - R_at(:,nprobe) ) / dist
    enddo
    grad_num_t3 = grad_num_t3 * embeddin_deriv1(rho_at(nprobe))

    do i = 1, verlet_list_length(nprobe)
        iv = verlet_list(nprobe,i)
        if (iv .eq. nprobe) cycle
        dist = norm2( R_at(:,iv) - R_at(:,nprobe) )
        grad_num_t4 = grad_num_t4 + &
        embeddin_deriv1(rho_at(iv)) * edensity_deriv1(dist) * ( R_at(:,iv) - R_at(:,nprobe) ) / dist
    enddo
    grad_num = grad_num_t1 + grad_num_t2 + grad_num_t3 + grad_num_t4

    print 311, "force term 1", grad_num_t1
    print 311, "force term 2", grad_num_t2
    print 311, "force term 3", grad_num_t3
    print 311, "force term 4", grad_num_t4

    cosine_v = dot_product(shift,grad_num) / ( norm2(shift) * norm2(grad_num) )
    print *, "cosine is ", cosine_v
    print"(A,3(1x,es12.5),/,/)", "calculating gradient via numerical derivative formula: done",grad_num

    if (cosine_v .gt. -0.99) then
        print*, "inappropriate cosine", cosine_v
        print*, "electronic densities"
        print"(SP, 200(1x, F15.10, /))", rho_at(verlet_list(nprobe,1:verlet_list_length(nprobe)))
    endif

    deallocate(rho_at)
    R_at(:,nprobe) = R_at(:,nprobe) - shift

301 format(A,3(1x,F12.8))
311 format(SP, A,3(1x,F14.10))


endsubroutine test_verlet_list
