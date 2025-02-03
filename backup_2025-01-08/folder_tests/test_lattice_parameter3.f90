subroutine test_lattice_parameter3
    use positions_mod
    implicit none
    integer(integer_byte_size) i, ncent
    real(float_byte_size) latpar, lpmin, es, esmin, rnum
    es=1d0
    esmin=1d0
    do i = 36400000, 37200000, 15
        call random_number(rnum)
        latpar = i*1d-7+1d-8*(rnum*2d0 -1d0)
        call place_atoms_fcc(latpar)
        call find_centeral_atom(ncent)
        call calculate_energy_at(es, ncent)
        if (es .lt. esmin) then
            esmin = es
            lpmin = latpar
        endif
        !print *, es, "eV is for latpar = ", latpar, " angstrem"
    enddo
    print*, "Test 3 result:"
    print*, "expected lattice parameter 3.7083"
    print*, "obtained lattice parameter ", lpmin
    print*, "Test 3 is over"
endsubroutine test_lattice_parameter3
