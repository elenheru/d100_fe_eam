subroutine test_lattice_parameter2
    use positions_mod
    implicit none
    integer i,ncent
    real(8) latpar, lpmin, es, esmin, epair,rnum
    es=1d0
    esmin=1d0
    do i = 28500000,28600000,1
        call random_number(rnum)
        latpar = i*1d-7+1d-8*(rnum*2d0 -1d0)
        call find_centeral_atom(ncent)
        call place_atoms_bcc(latpar)
        epair = 0d0
        call calculate_energy_at(es,ncent+1)
        epair = epair + es
        call calculate_energy_at(es,ncent-1)
        epair = epair + es
        if (epair .lt. esmin) then
            esmin = epair
            lpmin = latpar
        endif
        !print *, epair, "eV is for latpar = ", latpar, " angstrem"
    enddo
    print*,"Test 2 results:"
    print*, "expected lattice parameter 2.8557312"
    print*, "obtained lattice parameter ", lpmin
    print*,"Test 2 is over"
endsubroutine test_lattice_parameter2
