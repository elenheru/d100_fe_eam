subroutine place_atoms_bcc_spheric(lp)
    use positions_mod
    use system_parameters_mod

    real(float_byte_size), intent(in) :: lp
    integer(integer_byte_size) :: iw
    integer(integer_byte_size) :: i, j, k, na

    iw = nint( 1.41 * (dr_elastic_zone + r_main_cell) / lp ) + 6
    last_relaxable=0
    do i=-iw, iw
    do j=-iw, iw
    do k=-iw, iw
        if ( sqrt(i**2.0 + j**2 + k**2)*lp .gt. r_main_cell ) cycle
        last_relaxable=last_relaxable+1
        R_at(1, last_relaxable) = i * lp
        R_at(2, last_relaxable) = j * lp
        R_at(3, last_relaxable) = k * lp
    enddo
    enddo
    enddo

    do i=-iw, iw
    do j=-iw, iw
    do k=-iw, iw
        if ( sqrt((i+0.5)**2 + (j+0.5)**2 + (k+0.5)**2 )*lp .gt. r_main_cell ) cycle
        last_relaxable=last_relaxable+1
        R_at(1, last_relaxable) = i * lp + lp * 0.5
        R_at(2, last_relaxable) = j * lp + lp * 0.5
        R_at(3, last_relaxable) = k * lp + lp * 0.5
    enddo
    enddo
    enddo

    last_wall = last_relaxable

    do i=-iw, iw
    do j=-iw, iw
    do k=-iw, iw
        if ( sqrt(i**2.0 + j**2 + k**2)*lp .lt. r_main_cell ) cycle
        if ( sqrt(i**2.0 + j**2 + k**2)*lp .gt. (r_main_cell + dr_elastic_zone) ) cycle
        last_wall=last_wall+1
        R_at(1, last_wall) = i * lp
        R_at(2, last_wall) = j * lp
        R_at(3, last_wall) = k * lp
    enddo
    enddo
    enddo

    do i=-iw, iw
    do j=-iw, iw
    do k=-iw, iw
        if ( sqrt((i+0.5)**2 + (j+0.5)**2 + (k+0.5)**2 )*lp .lt. r_main_cell ) cycle
        if ( sqrt((i+0.5)**2 + (j+0.5)**2 + (k+0.5)**2 )*lp .gt. (r_main_cell + dr_elastic_zone) ) cycle
        last_wall=last_wall+1
        R_at(1, last_wall) = i * lp + lp * 0.5
        R_at(2, last_wall) = j * lp + lp * 0.5
        R_at(3, last_wall) = k * lp + lp * 0.5
    enddo
    enddo
    enddo
    print*, last_wall, " atoms placed in BCC spheric cell. Last atom in main cell is ", last_relaxable
    !print*, na, " atoms placed in BCC"
    !last_wall = na
    !last_relaxable = na

endsubroutine place_atoms_bcc_spheric
