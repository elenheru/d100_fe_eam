subroutine test_lattice_parameter5
    use positions_mod
    implicit none
    integer(integer_byte_size) i
    real(float_byte_size) latpar, lpmin, e_min, e_system, rnum

    e_min=1.0
    do i = 28557220, 28557360, 1
        call random_number(rnum)
        latpar = i*1d-7+1d-8*(rnum*2d0 -1d0)
        call place_atoms_bcc(latpar)
        call calculate_system_energy_verlet(e_system)
        print *, e_system, "eV is for latpar = ", latpar, " angstrem"
        if (e_system .lt. e_min) then
            e_min = e_system
            lpmin = latpar
        endif
    enddo
    call place_atoms_bcc(lpmin)
    call calculate_system_energy(e_system)
    print*, "Test 5 results:"
    print*, "expected lattice parameter    2.8557312"
    print*, "obtained lattice parameter ", lpmin
    print*, "respected verlet energy is ", e_min
    print*, "respected fair L energy is ", e_system
    print*, "Test 5 is over"
endsubroutine test_lattice_parameter5
