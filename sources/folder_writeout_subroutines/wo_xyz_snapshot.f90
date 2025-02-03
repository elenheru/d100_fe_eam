subroutine wo_xyz_snapshot(special_one)
    use wo_storage_mod
    use compile_debug_options_mod, only : debug_flag
    use positions_mod
    use system_parameters_mod, only : a0bcc
    integer i_atoms
    integer main_cell_counter
    integer elastic_zone_counter
    integer atom_counter
    integer(integer_byte_size), intent(in) :: special_one
    character (LEN = 29) filename_snapshot
    real(float_byte_size) :: wo_z_boundary_start
    real(float_byte_size) :: wo_z_boundary_end

    write(filename_snapshot, '(A, I5.5, A)'), &
    "dislocation_snaphot_", xyz_snapshot_counter, ".xyz"
    open (1405, file = filename_snapshot)

    wo_z_boundary_start = wo_z_boundary_start_111
    wo_z_boundary_end = wo_z_boundary_end_111

    if (show_single_z_plane) then

        atom_counter = 0
        do i_atoms = 1, last_wall
            if (r_at(3, i_atoms) .lt. wo_z_boundary_start) cycle
            if (r_at(3, i_atoms) .gt. wo_z_boundary_end) cycle
            atom_counter = atom_counter + 1
        enddo

        write(1405, *), atom_counter
        write(1405, *), "image of relaxation # ", xyz_snapshot_counter

        do i_atoms = first_relaxable, last_relaxable
            if (r_at(3, i_atoms) .lt. wo_z_boundary_start) cycle
            if (r_at(3, i_atoms) .gt. wo_z_boundary_end) cycle
            if (i_atoms .eq. special_one) then
                write(1405, 194) interesting_atom, r_at(1:3, i_atoms) &
                * wo_distance_multiplier
            else
                write(1405, 194) main_cell_atom, r_at(1:3, i_atoms) &
                * wo_distance_multiplier
            endif
        enddo
        do i_atoms = first_casing, last_casing
            if (r_at(3, i_atoms) .lt. wo_z_boundary_start) cycle
            if (r_at(3, i_atoms) .gt. wo_z_boundary_end) cycle
            write(1405, 194) casing_atom, r_at(1:3, i_atoms) &
            * wo_distance_multiplier
        enddo
        do i_atoms = first_wall, last_wall
            if (r_at(3, i_atoms) .lt. wo_z_boundary_start) cycle
            if (r_at(3, i_atoms) .gt. wo_z_boundary_end) cycle
            write(1405, 194) wall_atom, r_at(1:3, i_atoms) &
            * wo_distance_multiplier
        enddo

    elseif (only_centeral_nine_atoms) then
        write(1405, *), 9
        write(1405, *), "image of relaxation # ", xyz_snapshot_counter
        do i_atoms = first_relaxable, last_relaxable
            if (norm2(r_at(1:3, i_atoms)) .lt. a0bcc * sqrt(0.75 + 0.15)) then
                write(1405, 194) interesting_atom, r_at(1:3, i_atoms) &
                * wo_distance_multiplier
            endif
        enddo

    else

        write(1405, *), last_wall
        write(1405, *), "image of relaxation # ", xyz_snapshot_counter

        do i_atoms = first_relaxable, last_relaxable
            if (i_atoms .eq. special_one) then
                write(1405, 194) interesting_atom, r_at(1:3, i_atoms) &
                * wo_distance_multiplier
            else
                write(1405, 194) main_cell_atom, r_at(1:3, i_atoms) &
                * wo_distance_multiplier
            endif
        enddo
        do i_atoms = first_casing, last_casing
            write(1405, 194) casing_atom, r_at(1:3, i_atoms) &
            * wo_distance_multiplier
        enddo
        do i_atoms = first_wall, last_wall
            write(1405, 194) wall_atom, r_at(1:3, i_atoms) &
            * wo_distance_multiplier
        enddo
    endif

    if((mod(xyz_snapshot_counter, 25).eq.0) .and. debug_flag) then
        print'(A, I5.5, A)', &
        " Writing out XYZ picture # ", xyz_snapshot_counter, " called " // &
        filename_snapshot
    endif

    xyz_snapshot_counter = xyz_snapshot_counter+1

    close(1405)
194 format(A2, 3(3x, ES11.4))
endsubroutine wo_xyz_snapshot
