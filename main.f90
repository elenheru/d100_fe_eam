! this uses 1/2 multiplier in pairwise term

program dislocation_reconstruction
    use compile_debug_options_mod, only : timers

    integer ir

    call random_seed()
    call system_clock(timers(1))
    call place_atoms_bcc
    call initiate_energetic_parameters
    call allocate_electronic_densities
    call obtain_verlet_list
!    call wo_an_pw
!    call wo_an_ed
!    call wo_an_mf
!    call test_lattice_parameter2
!    call test_lattice_parameter1
!    call test_lattice_parameter4
!    call test_lattice_parameter3
!    call test_lattice_parameter5
!    call test_atom_relax_bcc
!    call test_atom_relax_fcc
!    call test_vacancy_relax_bcc
!    call test_derivatives
!    call test_linear_approximation_of_potentials
!    call test_wo_values_la
!    call test_po_value_potential
!    call test_atom_gradient
!    call test_far_atom_gradient
!    call test_verlet_list
    call test_field_work
    call test_writeout_snapshot

!    call execute_command_line("sleep 0.4", wait = .true.)
endprogram dislocation_reconstruction
