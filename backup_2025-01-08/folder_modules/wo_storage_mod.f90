module wo_storage_mod
    use compile_debug_options_mod, only : float_byte_size, integer_byte_size
    use system_parameters_mod, only : a0bcc
    save
    integer(integer_byte_size) :: xyz_snapshot_counter
    logical(1) :: show_single_z_plane = .false.
    real(float_byte_size), parameter :: wo_z_boundary_start = -1d-1 * a0bcc
    real(float_byte_size), parameter :: wo_z_boundary_end = 6d-1 * a0bcc
    real(float_byte_size), parameter :: wo_distance_multiplier = 4.0
    character(LEN = 3) :: main_cell_atom = "Fe "
    character(LEN = 3) :: wall_atom = "Al "
    character(LEN = 3) :: casing_atom = "Ca "
    character(LEN = 3) :: periodic_atom = "Pd "
    character(LEN = 3) :: interesting_atom = "Au "
    !character(LEN = 24) :: xyz_filename = "01234567890123456789_-_-"
endmodule wo_storage_mod
