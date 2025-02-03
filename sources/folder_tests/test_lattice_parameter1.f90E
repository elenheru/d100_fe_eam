subroutine test_lattice_parameter1
    use positions_mod
    implicit none
    integer(integer_byte_size) i, ncent
    real(float_byte_size) latpar, lpmin, es, esmin, rnum
    es=1d0
    esmin=1d0
    !ncent=9261
    do i = 28500000, 28600000, 1
        call random_number(rnum)
        latpar = i*1d-7+1d-8*(rnum*2d0 -1d0)
        call place_atoms_bcc(latpar)
        call find_centeral_atom(ncent)
        call calculate_energy_at(es, ncent)
        if (es .lt. esmin) then
            esmin = es
            lpmin = latpar
        endif
        !print *, es, "eV is for latpar = ", latpar, " angstrem"
    enddo
    print*, "Test 1 result:"
    print*, "expected lattice parameter 2.8557312"
    print*, "obtained lattice parameter ", lpmin
    print*, "Test 1 is over"
endsubroutine test_lattice_parameter1
