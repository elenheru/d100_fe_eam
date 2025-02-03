subroutine place_atoms_bcc

    use system_parameters_mod, only : a0bcc
!    call place_atoms_bcc_spheric(a0bcc)
    call place_atoms_bcc_cylindric(a0bcc)

endsubroutine place_atoms_bcc
