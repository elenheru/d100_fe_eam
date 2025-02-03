subroutine test_vacancy_relax_bcc
    use positions_mod
    use system_parameters_mod
    implicit none
    integer(integer_byte_size) i, j, ncent
    real(float_byte_size) eat, dist, rnum

    call place_atoms_bcc_spheric
    call find_centeral_atom(ncent)
    print*, ncent, "is centeral atom"
    R_at(1, ncent) = (dr_elastic_zone + r_main_cell)*15d-1
    R_at(2, ncent) = (dr_elastic_zone + r_main_cell)*15d-1
    R_at(3, ncent) = (dr_elastic_zone + r_main_cell)*15d-1

    call wo_xyz_snapshot
    call relax_rude
    call wo_xyz_snapshot
    call relax_fine
    call wo_xyz_snapshot

endsubroutine test_vacancy_relax_bcc
