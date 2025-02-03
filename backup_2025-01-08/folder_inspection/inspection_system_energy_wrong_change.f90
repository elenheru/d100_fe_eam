subroutine inspection_system_energy_wrong_change
    use positions_mod
    use verlet_list_mod
    use compile_debug_options_mod, only : atom_of_interest

    print *, verlet_list_length(atom_of_interest)

endsubroutine inspection_system_energy_wrong_change
