subroutine calculate_force_atom_verlet(number_atom, force)

    use compile_debug_options_mod, only : debug_flag, debug_string_format
    use potentials_ackland4fefe_linear_mod
    use positions_mod
    use system_parameters_mod
    use electronic_densities_allocated_mod
    use verlet_list_mod

    real(float_byte_size), intent(inout), dimension(3) :: force
    integer(integer_byte_size), intent(in) :: number_atom

    real(float_byte_size), dimension(3) :: grad_num_t1, grad_num_t2, grad_num_t3, grad_num_t4, grad_num
    real(float_byte_size) :: dist, rnum, d1v_ik, d1, rn, cosine_v
    integer(integer_byte_size) iv, i, jv, j

    rho_at = 0d0
    grad_num_t1 = 0d0
    grad_num_t2 = 0d0
    grad_num_t3 = 0d0
    grad_num_t4 = 0d0

    do iv = 1, verlet_list_length(number_atom)
        i = verlet_list(number_atom, iv)
        if (i .eq. number_atom) STOP "collision in atom indices 1"
        dist = norm2( R_at(:, i) - R_at(:, number_atom) )
        grad_num_t1 = grad_num_t1 + pairwise_deriv1(dist) * ( R_at(:, i) - R_at(:, number_atom) ) / dist
    enddo

    ! technically this is a reliable cycle, because sum has two parts
    ! these parts are always equal though, |r_i - r_j| = |r_j - r_i|
    ! but differentiation considers r_i and r_j as different variables
    ! this is where second cycle comes from

    do jv = 1, verlet_list_length(number_atom)
        j = verlet_list(number_atom, jv)
        if (j .eq. number_atom) STOP "collision in atom indices 2"
        dist = norm2( R_at(:, j) - R_at(:, number_atom) )
        grad_num_t2 = grad_num_t2 + pairwise_deriv1(dist) * ( R_at(:, j) - R_at(:, number_atom) ) / dist
    enddo

    do jv = 1, verlet_list_length(number_atom)
        j = verlet_list(number_atom, jv)
        do iv = 1, verlet_list_length(j)
            i = verlet_list(j, iv)
            if (i .eq. j) STOP "collision in atom indices 3"
            dist = norm2( R_at(:, j) - R_at(:, i) )
            rho_at(j) = rho_at(j) + edensity(dist)
            if (i .eq. number_atom) rho_at(i) = rho_at(i) + edensity(dist)
        enddo
    enddo

    ! The need in violation of formula arises from fact that we do not use
    ! periodic boundary conditions, while it is implied in energy formula given in article.
    ! To be consistent we have to split elastic area into two parts, outer and and inner
    ! outer is adjacent to vacuum, and inner is adjacent to atoms of main cell.
    ! Inner part is considered to be part of system in sence that its atoms embedding_function
    ! is required to be calculated for system energy, meanwhile for outer part
    ! such calculation is required only when the force is being calculated.
    ! Inner elastic atoms are not allowed to move - this is the difference
    ! from atoms in main computational cell.
    ! There exists a similar approach - to split main cell in two zones
    ! and then anchor atoms of its outer part.

    do jv = 1, verlet_list_length(number_atom)
        j = verlet_list(number_atom, jv)
        if (j .eq. number_atom) STOP "collision in atom indices 4"
        dist = norm2( R_at(:, j) - R_at(:, number_atom) )
        grad_num_t3 = grad_num_t3 + edensity_deriv1(dist) * ( R_at(:, j) - R_at(:, number_atom) ) / dist
    enddo
    grad_num_t3 = grad_num_t3 * embeddin_deriv1(rho_at(number_atom))

    do iv = 1, verlet_list_length(number_atom)
        i = verlet_list(number_atom, iv)
        if (i .eq. number_atom) STOP "collision in atom indices 5"
        dist = norm2( R_at(:, i) - R_at(:, number_atom) )
        grad_num_t4 = grad_num_t4 + &
        embeddin_deriv1(rho_at(i)) * edensity_deriv1(dist) * ( R_at(:, i) - R_at(:, number_atom) ) / dist
    enddo
    !there is a room for optimization - to store distances to verlet list atoms
    !this will save several negations mutiplications and squareroots
    grad_num = - grad_num_t1 - grad_num_t2 - grad_num_t3 - grad_num_t4
    force = grad_num

    if (debug_flag) then
        print 311, "Numerical g1", grad_num_t1
        print 311, "Numerical g2", grad_num_t2
        print 311, "Numerical g3", grad_num_t3
        print 311, "Numerical g4", grad_num_t4

        print"(A, 3(1x, es12.5), /, /)", &
        "calculating gradient via analytical derivative formula: done", force
    endif


301 format(A, 3(1x, F12.8))
311 format(SP, A, 3(1x, F14.10))


endsubroutine calculate_force_atom_verlet
