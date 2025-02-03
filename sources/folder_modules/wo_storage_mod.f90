module wo_storage_mod
    use compile_debug_options_mod, only : float_byte_size, integer_byte_size
    use system_parameters_mod, only : a0bcc
    save
    integer(integer_byte_size) :: xyz_snapshot_counter
    logical(1) :: show_single_z_plane = .false.
    logical(1) :: only_centeral_nine_atoms = .false.

    real(float_byte_size), parameter :: wo_z_boundary_start_100 = &
        -0.01 * a0bcc
    real(float_byte_size), parameter :: wo_z_boundary_end_100 = &
        0.99 * a0bcc

    real(float_byte_size), parameter :: wo_z_boundary_start_110 = &
        -0.01 * a0bcc * sqrt(2.0)
    real(float_byte_size), parameter :: wo_z_boundary_end_110 = &
        0.99 * a0bcc * sqrt(2.0)

    real(float_byte_size), parameter :: wo_z_boundary_start_111 = &
        -0.01 * a0bcc * sqrt(3.0)
    real(float_byte_size), parameter :: wo_z_boundary_end_111 = &
        0.99 * a0bcc * sqrt(3.0)

    real(float_byte_size), parameter :: wo_distance_multiplier = 8.0
    character(LEN = 3) :: main_cell_atom = "Fe "
    character(LEN = 3) :: wall_atom = "Al "
    character(LEN = 3) :: casing_atom = "Ca "
    character(LEN = 3) :: periodic_atom = "Pd "
    character(LEN = 3) :: interesting_atom = "Au "
    !character(LEN = 24) :: xyz_filename = "01234567890123456789_-_-"
endmodule wo_storage_mod
