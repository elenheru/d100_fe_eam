subroutine obtain_verlet_list
    use positions_mod
    use verlet_list_mod
    use compile_debug_options_mod, only : debug_flag
    implicit none
    integer(integer_byte_size) :: i, j, jb
    integer(integer_byte_size) :: number_of_ajacents
    integer(integer_byte_size) :: ajacents(27 * max_atoms_in_box)
    real(float_byte_size) :: distance_squared = verlet_cutoff ** 2
    real(float_byte_size) :: R_delta(3) = 0

    if (debug_flag) print*, "building verlet list for atoms "

    verlet_list_length = 0
    verlet_list = 0
    call obtain_boxes_content
    do i = 1, last_wall
        call get_content_of_adjacent_boxes(i, ajacents, number_of_ajacents)
        do jb = 1, number_of_ajacents
            j = ajacents(jb)
            if (i .eq. j) cycle
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
            endif
        enddo
    enddo

    contains
    subroutine obtain_boxes_content
        integer(integer_byte_size) :: ix, iy, iz, n
        boxes_indices = 0
        how_much_atoms_in_box = 0
        which_atom_in_box = 0
        do n = 1, last_wall
            ix = nint((R_at(1,n) + x_offset) / verlet_cutoff)
            iy = nint((R_at(2,n) + y_offset) / verlet_cutoff)
            iz = nint((R_at(3,n) + z_offset) / verlet_cutoff)
            if (ix .lt. 1) ix = 1
            if (iy .lt. 1) iy = 1
            if (iz .lt. 1) iz = 1
            if (ix .gt. boxes_in_line) ix = boxes_in_line
            if (iy .gt. boxes_in_line) iy = boxes_in_line
            if (iz .gt. boxes_in_line) iz = boxes_in_line
            boxes_indices(:,n) = (/ix, iy, iz/)
            how_much_atoms_in_box(ix, iy, iz) = &
                how_much_atoms_in_box(ix, iy, iz) + 1
            which_atom_in_box(how_much_atoms_in_box(ix, iy, iz), ix, iy, iz) = n
        enddo
    endsubroutine obtain_boxes_content

    subroutine get_content_of_adjacent_boxes(n_atom, neighbors, neibors_number)
        integer, intent(in) :: n_atom
        integer, intent(inout) :: neighbors(27 * max_atoms_in_box)
        integer, intent(inout) :: neibors_number
        integer(integer_byte_size) ix, iy, iz
        integer(integer_byte_size) jx, jy, jz, jn

        neighbors = 0
        neibors_number = 0
        ix = boxes_indices(1, n_atom)
        iy = boxes_indices(2, n_atom)
        iz = boxes_indices(3, n_atom)
        do jx = ix - 1, ix + 1
        do jy = iy - 1, iy + 1
        do jz = iz - 1, iz + 1
            if (     (jx .lt. 1) &
                .or. (jy .lt. 1) &
                .or. (jz .lt. 1) &
                .or. (jx .gt. boxes_in_line) &
                .or. (jy .gt. boxes_in_line) &
                .or. (jz .gt. boxes_in_line))&
            then
                cycle
            endif
            jn = how_much_atoms_in_box(jx, jy, jz)
            if (jn .eq. 0) cycle
            neighbors(neibors_number + 1: neibors_number + jn) = &
                which_atom_in_box(1:jn, jx, jy, jz)
            neibors_number = neibors_number + jn
        enddo
        enddo
        enddo

    endsubroutine get_content_of_adjacent_boxes

endsubroutine obtain_verlet_list
























