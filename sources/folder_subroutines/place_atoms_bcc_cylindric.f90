subroutine place_atoms_bcc_cylindric(lp)
    use positions_mod
    use rotations_mod
    use system_parameters_mod

    real(float_byte_size), intent(in) :: lp
    logical :: zone_geometry_condition
    integer(integer_byte_size) :: iw, iz
    integer(integer_byte_size) :: i, j, k
    real(float_byte_size) :: atom_position(3)


    real(float_byte_size), dimension(3,3) :: matrix_100
    real(float_byte_size), dimension(3,3) :: matrix_110
    real(float_byte_size), dimension(3,3) :: matrix_111
    real(float_byte_size), dimension(3,3) :: matrix_to_use

    iw = nint( 1.41 * (&
        r_main_cell &
        + dr_casing_zone &
        + dr_elastic_zone &
        ) / lp )*3 + 10
    iz = nint( (&
        h_main_cell &
        + dh_casing_zone &
        + dh_elastic_zone) &
        / lp )*3 + 10

    matrix_100 = matrix_unit
    matrix_to_use = matrix_unit
    call build_rotation_matrix(matrix_to_use, "Y", 45.0_float_byte_size)
    matrix_110 = matrix_to_use
    call build_rotation_matrix(matrix_to_use, "X", degree_110_to_111)
    matrix_111 = matrix_to_use
    !matrix_to_use = matrix_110

    ! start relaxable placing
    matrix_to_use = matrix_111
    first_relaxable = 1
    last_relaxable = 1
    do i=-iw, iw
    do j=-iw, iw
    do k=-iz, iz
        atom_position(1) = i * lp
        atom_position(2) = j * lp
        atom_position(3) = k * lp
        atom_position = matmul(matrix_to_use, atom_position)
        zone_geometry_condition = &
            zone_geometry_condition_relaxable(atom_position)
        if ( zone_geometry_condition ) then
            R_at(:, last_relaxable) = atom_position(:)
            last_relaxable = last_relaxable + 1
        endif

        atom_position(1) = i * lp + lp * 0.5
        atom_position(2) = j * lp + lp * 0.5
        atom_position(3) = k * lp + lp * 0.5
        atom_position = matmul(matrix_to_use, atom_position)
        zone_geometry_condition = &
            zone_geometry_condition_relaxable(atom_position)
        if ( zone_geometry_condition ) then
            R_at(:, last_relaxable) = atom_position(:)
            last_relaxable = last_relaxable + 1
        endif
    enddo
    enddo
    enddo
    last_relaxable = last_relaxable - 1
    ! end relaxable placing

    ! start casing placing
    first_casing = last_relaxable + 1
    last_casing = first_casing
    do i=-iw, iw
    do j=-iw, iw
    do k=-iz, iz
        atom_position(1) = i * lp
        atom_position(2) = j * lp
        atom_position(3) = k * lp
        atom_position = matmul(matrix_to_use, atom_position)
        zone_geometry_condition = &
            zone_geometry_condition_casing(atom_position)

        if ( zone_geometry_condition ) then
            R_at(:, last_casing) = atom_position(:)
            last_casing = last_casing + 1
        endif

        atom_position(1) = i * lp + lp * 0.5
        atom_position(2) = j * lp + lp * 0.5
        atom_position(3) = k * lp + lp * 0.5
        atom_position = matmul(matrix_to_use, atom_position)
        zone_geometry_condition = &
            zone_geometry_condition_casing(atom_position)

        if ( zone_geometry_condition ) then
            R_at(:, last_casing) = atom_position(:)
            last_casing = last_casing + 1
        endif
    enddo
    enddo
    enddo
    last_casing = last_casing - 1
    ! end casing placing

    ! start elastic placing
    first_wall= last_casing + 1
    last_wall = first_wall
    do i=-iw, iw
    do j=-iw, iw
    do k=-iz, iz
        atom_position(1) = i * lp
        atom_position(2) = j * lp
        atom_position(3) = k * lp
        atom_position = matmul(matrix_to_use, atom_position)
        zone_geometry_condition = &
            zone_geometry_condition_wall(atom_position)

        if ( zone_geometry_condition ) then
            R_at(:, last_wall) = atom_position(:)
            last_wall = last_wall + 1
        endif

        atom_position(1) = i * lp + lp * 0.5
        atom_position(2) = j * lp + lp * 0.5
        atom_position(3) = k * lp + lp * 0.5
        atom_position = matmul(matrix_to_use, atom_position)
        zone_geometry_condition = &
            zone_geometry_condition_wall(atom_position)

        if ( zone_geometry_condition ) then
            R_at(:, last_wall) = atom_position(:)
            last_wall = last_wall + 1
        endif
    enddo
    enddo
    enddo
    last_wall = last_wall - 1
    ! end elastic placing

    print*, last_wall, " atoms placed in BCC cylindric cell." &
    // " Last atom in main cell is ", last_relaxable

    210 format(SP, 3(3(1x, F10.7), /))
    contains
        logical function zone_geometry_condition_relaxable(r)
            real(float_byte_size), intent(in) :: r(3)
            zone_geometry_condition_relaxable = (&
            ( norm2(r(1:2)) .lt. r_main_cell ) .and. &
            ( abs(r(3)) .lt. h_main_cell) )
        endfunction

        logical function zone_geometry_condition_casing(r)
            real(float_byte_size), intent(in) :: r(3)
            zone_geometry_condition_casing = ( (&
            ( norm2(r(1:2)) .lt. r_main_cell + dr_casing_zone ) .and. &
            ( abs(r(3)) .lt. h_main_cell + dh_casing_zone) )&
            .and. (.not. zone_geometry_condition_relaxable(r))  )
        endfunction

        logical function zone_geometry_condition_wall(r)
            real(float_byte_size), intent(in) :: r(3)
            zone_geometry_condition_wall = ( (&
            ( norm2(r(1:2)) .lt. &
            r_main_cell + dr_casing_zone + dr_elastic_zone) .and. &
            ( abs(r(3)) .lt. &
            h_main_cell + dh_casing_zone + dh_elastic_zone))&
            .and. ( (.not. zone_geometry_condition_relaxable(r))&
            .and. (.not. zone_geometry_condition_casing(r)) ) )
        endfunction
endsubroutine place_atoms_bcc_cylindric
