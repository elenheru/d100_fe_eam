module periodicity_z_mod
    use compile_debug_options_mod, only : float_byte_size, integer_byte_size
    use positions_mod, only : max_at
    save
    integer(integer_byte_size), dimension(max_at) :: &
        paragon_atoms, acceptor_atoms
    integer(integer_byte_size) :: &
        paragon_atoms_total_number, &
        acceptor_atoms_total_number

endmodule periodicity_z_mod
