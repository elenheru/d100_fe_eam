subroutine test_different_relaxations_bcc
    use positions_mod
    use system_parameters_mod
    implicit none
    integer i,j,ncent
    real(8) eat,dist,rnum

    call place_atoms_bcc_spheric
    call find_centeral_atom(ncent)
    print*, ncent, "is centeral atom"
    R_at(1,ncent) = (dr_elastic_zone + r_main_cell)*15d-1
    R_at(2,ncent) = (dr_elastic_zone + r_main_cell)*15d-1
    R_at(3,ncent) = (dr_elastic_zone + r_main_cell)*15d-1

    call wo_xyz_snapshot
!    call relax_rude(ncent)
    call wo_xyz_snapshot

    call place_atoms_bcc_spheric
    call find_centeral_atom(ncent)
    print*, ncent, "is centeral atom"
    R_at(1,ncent) = (dr_elastic_zone + r_main_cell)*15d-1
    R_at(2,ncent) = (dr_elastic_zone + r_main_cell)*15d-1
    R_at(3,ncent) = (dr_elastic_zone + r_main_cell)*15d-1

    call wo_xyz_snapshot
!    call relax_fine(ncent)
    call wo_xyz_snapshot

    call place_atoms_bcc_spheric
    call find_centeral_atom(ncent)
    print*, ncent, "is centeral atom"
    R_at(1,ncent) = (dr_elastic_zone + r_main_cell)*15d-1
    R_at(2,ncent) = (dr_elastic_zone + r_main_cell)*15d-1
    R_at(3,ncent) = (dr_elastic_zone + r_main_cell)*15d-1

    call wo_xyz_snapshot
!    call relax_orthogonal_randomized(ncent)
    call wo_xyz_snapshot


endsubroutine test_different_relaxations_bcc
