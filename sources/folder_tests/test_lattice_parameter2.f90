subroutine test_lattice_parameter2
    use positions_mod
    implicit none
    integer(integer_byte_size) i, ncent
    real(float_byte_size) latpar, lpmin, es, esmin, epair, rnum
    es=1d0
    esmin=1d0
    do i = 28500000, 28600000, 9
        call random_number(rnum)
        latpar = i*1d-7+1d-8*(rnum*2d0 -1d0)
!        call find_centeral_atom(ncent) ! why was it here at all?!
        call place_atoms_bcc(latpar)
        call find_centeral_atom(ncent)
        epair = 0d0
!        call calculate_energy_at(es, ncent+1)
!        epair = epair + es
        call calculate_energy_at(es, ncent-1)
        epair = epair + es
        if (epair .lt. esmin) then
            esmin = epair
            lpmin = latpar
        endif
        !print *, epair, "eV is for latpar = ", latpar, " angstrem"
    enddo
    print*, "Test 2 results:"
    print*, "expected lattice parameter    2.8557312"
    print*, "obtained lattice parameter ", lpmin
    print*, "Test 2 is over"
endsubroutine test_lattice_parameter2
