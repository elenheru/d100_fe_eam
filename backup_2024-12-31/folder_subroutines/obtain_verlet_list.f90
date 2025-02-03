subroutine obtain_verlet_list
    use positions_mod
    use verlet_list_mod
    use miscellaneous_parameters_mod, only : debug
    implicit none
    integer i, j
    real(8) dist
    real(8) R_delta(3)

    if (debug) print*, "building verlet list for atoms "

    do i = 1, nat
        do j = 1 + i, nat

            R_delta(1) = abs(R_at(1,i) - R_at(1,j))
            if (R_delta(1) .gt. verlet_cutoff) cycle

            R_delta(2) = abs(R_at(2,i) - R_at(2,j))
            if (R_delta(2) .gt. verlet_cutoff) cycle

            R_delta(3) = abs(R_at(3,i) - R_at(3,j))
            if (R_delta(3) .gt. verlet_cutoff) cycle

            dist = norm2(R_delta)

            if (dist .lt. verlet_cutoff) then
                verlet_list_length(i) = verlet_list_length(i) + 1
                verlet_list( i, verlet_list_length(i) ) = j

                verlet_list_length(j) = verlet_list_length(j) + 1
                verlet_list( j, verlet_list_length(j) ) = i
            endif
        enddo
    enddo

endsubroutine obtain_verlet_list
