subroutine test_atom_relax_bcc
    use positions_mod
    implicit none
    integer(integer_byte_size) i, j, ncent
    real(float_byte_size) latpar, eat, dist, rnum

    dist=0
    latpar = 2.85573 !BCC ackland
    !latpar = 2.86000 !BCC ollson
    !latpar = 3.65216 !FCC ackland
    !latpar = 3.708303 !FCC olsson
    call place_atoms_bcc(latpar)
    call find_centeral_atom(ncent)


    do i=15000, 1, -1
        !call shift_atom(ncent, i*1.6d-6)
        call shift_atom(ncent, i*1.2d-5)
        print 151, "atom #", ncent, " shifted to ", R_at(1:3, ncent)
        do j=1, 6
            call random_number(rnum)
            call relax_rude_quad_grad_rand(ncent, 5d-1**(rnum + real(j/3, 8)**125d-3) )
            call calculate_energy_at(eat, ncent)
            print 152, "atom #", ncent, " relaxed to ", R_at(1:3, ncent), " energy= ", eat
            call EXECUTE_COMMAND_LINE('sleep 0.05', wait=.true.)
        enddo
        call calculate_energy_at(eat, ncent)
        if (norm2(R_at(1:3, ncent)) .gt. 0d0 ) dist = dist + log( norm2(R_at(1:3, ncent)) )
        if(mod(i, 200) .eq. 1)print 152, "atom #", ncent, " relaxed to ", R_at(1:3, ncent), " energy=", eat
        if(mod(i, 200) .eq. 1)print 153, "expected position :         ", 0d0, 0d0, 0d0, " geometric average dist =", exp(dist/i)
        print *, repeat("_", 150)
        call EXECUTE_COMMAND_LINE('sleep 5.0', wait=.true.)
    enddo
        print 153, "expected position :         ", 0d0, 0d0, 0d0, " geometric average dist =", exp(dist/i)
    151 format(A, I10.4, A, 3(F23.18, 1x))
    152 format(A, I10.4, A, 3(F23.18, 1x), A, ES24.14e3)
    153 format(A, 3(F18.14, 1x), A, ES24.14e3)
endsubroutine test_atom_relax_bcc
