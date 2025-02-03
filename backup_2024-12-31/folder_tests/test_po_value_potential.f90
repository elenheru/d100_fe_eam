subroutine test_po_value_potential
    use potentials_ackland4fefe_linear_mod
    implicit none
    integer i
    real(8) :: pw_ack4pot, ed_ack4pot, mf_ack4pot
    real(8) :: function_value_an
    real(8) :: function_value_la
    !real(8) :: argument = 4.5948
    real(8) :: argument = 9.3815032954688751E-006
    call initiate_energetic_parameters
    print *, "test of a single value"
    do i = 0,9,1
        function_value_an = mf_ack4pot(argument*(1+0.1*i))
        function_value_la = embeddin(argument*(1+0.1*i))
        print 314, " arg = ", argument*(1+0.1*i)
        print 314, " tru val = ", function_value_an
        print 314, " lin val = ", function_value_la
        print *, " "
    enddo
    314 format(A,ES11.4,$)
endsubroutine test_po_value_potential
