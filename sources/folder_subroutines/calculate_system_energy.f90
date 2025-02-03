subroutine calculate_system_energy(e_system)
    use positions_mod
    use physical_parameters_mod
    use potentials_ackland4fefe_linear_mod
    use compile_debug_options_mod, only : debug_flag, debug_string_format
    use system_parameters_mod, only : cutoff_pw, cutoff_ed

    real(float_byte_size), intent(inout) :: e_system
    real(float_byte_size) r
    real(float_byte_size) electronic_density
    real(float_byte_size) e_atom
    real(float_byte_size) cutoff
    integer(integer_byte_size) i, j

    if (debug_flag) print debug_string_format, "Calculating system energy in a direct way ... "
    e_system = 0d0
    do j = 1, last_casing
        e_atom = 0d0
        electronic_density = 0d0
        do i = 1, last_wall
            if (i .eq. j) cycle
            r = norm2(R_at(1:3, i) - R_at(1:3, j))
            if(r .gt. cutoff_pw) cycle
            e_atom = e_atom + pairwise(r)
            if(r .gt. cutoff_ed) cycle
            electronic_density = electronic_density + edensity(r)
        enddo
        e_system = e_system + e_atom + embeddin(electronic_density)
    enddo
!    this is from investigation whether there should be 2 or 0.5 multiplier or not
!    do j = 1, last_relaxable
!        e_atom = 0d0
!        do i = 1, last_wall
!            if (i .eq. j) cycle
!            r=norm2(R_at(1:3, i) - R_at(1:3, j))
!            if(r .gt. cutoff_pw) cycle
!            e_atom = e_atom + pairwise(r) * 0.5
!        enddo
!        e_system = e_system + e_atom
!    enddo
!    do j = 1, last_relaxable
!        electronic_density = 0d0
!        do i = 1, last_wall
!            if (i .eq. j) cycle
!            r=norm2(R_at(1:3, i) - R_at(1:3, j))
!            if(r .gt. cutoff_ed) cycle
!            electronic_density = electronic_density + edensity(r)
!        enddo
!        e_system = e_system + embeddin(electronic_density)
!    enddo
    if (debug_flag) print *, "Done."
endsubroutine calculate_system_energy
