subroutine place_atoms_bcc_spheric
    use positions_mod
    use system_parameters_mod
    implicit none
    integer,parameter :: iw = nint( (dr_elastic_zone + r_main_cell) / a0bcc ) + 6
    integer i,j,k,na

    nat_mc=0
    do i=-iw,iw
    do j=-iw,iw
    do k=-iw,iw
        if ( (i**2 + j**2 + k**2)*a0bcc .gt. r_main_cell**2 ) cycle
        nat_mc=nat_mc+1
        R_at(1,nat_mc) = i * a0bcc
        R_at(2,nat_mc) = j * a0bcc
        R_at(3,nat_mc) = k * a0bcc
    enddo
    enddo
    enddo

    do i=-iw,iw
    do j=-iw,iw
    do k=-iw,iw
        if ( ((i+0.5)**2 + (j+0.5)**2 + (k+0.5)**2 )*a0bcc .gt. r_main_cell**2 ) cycle
        nat_mc=nat_mc+1
        R_at(1,nat_mc) = i * a0bcc + a0bcc * 0.5
        R_at(2,nat_mc) = j * a0bcc + a0bcc * 0.5
        R_at(3,nat_mc) = k * a0bcc + a0bcc * 0.5
    enddo
    enddo
    enddo

    nat = nat_mc

    do i=-iw,iw
    do j=-iw,iw
    do k=-iw,iw
        if ( (i**2 + j**2 + k**2)*a0bcc .lt. r_main_cell**2 ) cycle
        if ( (i**2 + j**2 + k**2)*a0bcc .gt. (r_main_cell + dr_elastic_zone)**2 ) cycle
        nat=nat+1
        R_at(1,nat) = i * a0bcc
        R_at(2,nat) = j * a0bcc
        R_at(3,nat) = k * a0bcc
    enddo
    enddo
    enddo

    do i=-iw,iw
    do j=-iw,iw
    do k=-iw,iw
        if ( ((i+0.5)**2 + (j+0.5)**2 + (k+0.5)**2 )*a0bcc .lt. r_main_cell**2 ) cycle
        if ( ((i+0.5)**2 + (j+0.5)**2 + (k+0.5)**2 )*a0bcc .gt. (r_main_cell + dr_elastic_zone)**2 ) cycle
        nat=nat+1
        R_at(1,nat) = i * a0bcc + a0bcc * 0.5
        R_at(2,nat) = j * a0bcc + a0bcc * 0.5
        R_at(3,nat) = k * a0bcc + a0bcc * 0.5
    enddo
    enddo
    enddo

    print*, nat, " atoms placed in BCC spheric cell. Last atom in main cell is ", nat_mc
endsubroutine place_atoms_bcc_spheric
