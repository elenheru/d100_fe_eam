module verlet_list_mod
    use compile_debug_options_mod
    use system_parameters_mod
    use positions_mod, only : max_at
    save
    real(float_byte_size), parameter :: verlet_cutoff = cutoff_pw + 0.1
    real(float_byte_size), parameter :: &
        reciprocal_verlet_cutoff = 1 / verlet_cutoff
    ! you may want to see potential parameters
    integer, parameter :: maximum_neibors = 200
    integer(integer_byte_size) :: verlet_list_length(max_at)
    integer(integer_byte_size) :: verlet_list(max_at, maximum_neibors)


    real(float_byte_size), parameter :: x_offset = &
        (r_main_cell + dr_casing_zone + dr_elastic_zone)
    real(float_byte_size), parameter :: y_offset = &
        (r_main_cell + dr_casing_zone + dr_elastic_zone)
    real(float_byte_size), parameter :: z_offset = &
        (h_main_cell + dh_casing_zone + dh_elastic_zone)
    integer, parameter :: max_atoms_in_box = &
        8 * 2 * nint((verlet_cutoff / a0bcc)**3)
    integer, parameter :: boxes_in_line = &
        2 * nint(max(&
            (r_main_cell + dr_casing_zone + dr_elastic_zone), &
            (h_main_cell + dh_casing_zone + dh_elastic_zone)) &
            / verlet_cutoff)
    integer(integer_byte_size) :: boxes_indices(3, max_at)
    ! boxes_indices(1:3, N) is (3,5,11) : the box where atom N is
    integer(integer_byte_size) :: &
        how_much_atoms_in_box(boxes_in_line, boxes_in_line, boxes_in_line)
    ! how_much_atoms_in_box(3,5,11) : how much atoms is in box 3 5 11
    integer(integer_byte_size) :: &
        which_atom_in_box(max_atoms_in_box, boxes_in_line, &
        boxes_in_line, boxes_in_line)
    ! which_atoms_in_box(:,3,5,11) : list of numbers of atoms in box 3 5 11
    integer(integer_byte_size) :: how_much_atoms_out_of_boxes
    integer(integer_byte_size) :: which_atoms_out_of_boxes(max_at)
endmodule verlet_list_mod
