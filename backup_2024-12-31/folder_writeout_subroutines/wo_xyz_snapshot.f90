subroutine      wo_xyz_snapshot !#104
    use positions_mod
    use system_parameters_mod
    use wo_storage_mod
    implicit none
    integer i_atoms
    character (LEN = 24) filename_trajectory
    xyz__WOcounter=xyz__WOcounter+1

    write(filename_trajectory,'(A,I5.5,A)'),"XYZ_trajectory_",xyz__WOcounter,".xyz"
    open (104, file = filename_trajectory)

    write(104,*),nat
    write(104,*),"IMAGE OF RELAXATION # ",xyz__wocounter

    do i_atoms=1,nat_mc
        write(104,194)  "Fe ", r_at(1:3,i_atoms)
    enddo
    do i_atoms=1+nat_mc,nat
        write(104,194)  " W ", r_at(1:3,i_atoms)
    enddo

    if(mod(xyz__wocounter,25).eq.0)&
    print'(A,I5.5,A)'," Writing out XYZ picture # ",xyz__WOcounter," . "
    close(104)
    194 format(A2,3(3x,ES11.4))
endsubroutine   wo_xyz_snapshot
