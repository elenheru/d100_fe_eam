subroutine calculate_system_energy_verlet(e_system)
    use positions_mod
    use physical_parameters_mod
    use potentials_ackland4fefe_linear_mod
    use compile_debug_options_mod, only : debug_flag, debug_string_format, atom_of_interest
    use verlet_list_mod
    use system_parameters_mod, only : cutoff_ed

    real(float_byte_size), intent(inout) :: e_system
    integer(integer_byte_size) i, j, iv
    real(float_byte_size) r
    real(float_byte_size) electronic_density
    real(float_byte_size) e_atom
        real(float_byte_size) pw_ack4pot, ed_ack4pot, mf_ack4pot
        ! cubic splines are better, or i have made a mistake in linearization

    if (debug_flag) then
        print debug_string_format, &
        "Calculate_system_energy_verlet " //&
        "(mind that it does not automaticly renew Verlet list) ... "
    endif
    e_system = 0

    do j = 1, last_relaxable
        electronic_density = 0
        e_atom = 0
        do i = 1, verlet_list_length(j)
            iv = verlet_list(j, i)

            !if(iv .eq. j) STOP "own atom in verlet list"
            r = norm2(R_at(1:3, iv) - R_at(1:3, j))
            if(r .lt. cutoff_ed) then
!                electronic_density = electronic_density + edensity(r)
                electronic_density = electronic_density + ed_ack4pot(r)
            endif
!            e_atom = e_atom + pairwise(r)
            e_atom = e_atom + pw_ack4pot(r)
            if (j .eq. atom_of_interest) then
                print *, iv , R_at(:, iv), "; ed , ea :", electronic_density, e_atom
            endif
            ! I just thought that G.J.Ackland made a mistake
            ! but no, line below does not help
            !if (iv .gt. last_relaxable) e_atom = e_atom + pw_ack4pot(r)
        enddo
!        e_system = e_system + e_atom + embeddin(electronic_density)
        e_system = e_system + e_atom + mf_ack4pot(electronic_density)
    enddo

    if (debug_flag)  print *, "Done."

endsubroutine calculate_system_energy_verlet
