subroutine obtain_verlet_list
    use positions_mod
    use verlet_list_mod
    use compile_debug_options_mod, only : debug_flag
    implicit none
    integer(integer_byte_size) i, j
    real(float_byte_size) distance_squared
    real(float_byte_size) R_delta(3)

    if (debug_flag) print*, "building verlet list for atoms "

    verlet_list_length = 0
    verlet_list = 0

    do i = 1, last_wall - 1
        do j = 1 + i, last_wall
            !if (i .eq. j) cycle
            R_delta(1) = abs(R_at(1, i) - R_at(1, j))
            if (R_delta(1) .gt. verlet_cutoff) cycle

            R_delta(2) = abs(R_at(2, i) - R_at(2, j))
            if (R_delta(2) .gt. verlet_cutoff) cycle

            R_delta(3) = abs(R_at(3, i) - R_at(3, j))
            if (R_delta(3) .gt. verlet_cutoff) cycle

            distance_squared = dot_product(R_delta, R_delta)
            if (distance_squared .lt. verlet_cutoff ** 2) then
                verlet_list_length(i) = verlet_list_length(i) + 1
                verlet_list( i, verlet_list_length(i) ) = j

                verlet_list_length(j) = verlet_list_length(j) + 1
                verlet_list( j, verlet_list_length(j) ) = i
            endif
        enddo
    enddo

endsubroutine obtain_verlet_list
