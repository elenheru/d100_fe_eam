subroutine test_linear_approximation_of_potentials
    use potentials_ackland4fefe_linear_mod
    implicit none
    integer(8), parameter :: rand_tries = 1000*1000*90!00
    real(float_byte_size), parameter :: argshift = 1d-6
    integer(8) :: i
    real(float_byte_size) :: argument, maxargument = 0d0
    real(float_byte_size) :: function_value_an
    real(float_byte_size) :: function_value_la
    real(float_byte_size) :: maxvalue
    real(float_byte_size) :: relative_error = 0d0, maxrelative_error = 0d0, averelative_error = 0d0
    real(float_byte_size) :: absolute_error = 0d0, maxabsolute_error = 0d0, aveabsolute_error = 0d0
    real(float_byte_size) :: pw_ack4pot, ed_ack4pot, mf_ack4pot
    real(float_byte_size) :: d1pw_ack4pot, d1ed_ack4pot, d1mf_ack4pot
    real(float_byte_size) :: pw_ack4pot_la, ed_ack4pot_la, mf_ack4pot_la
    real(float_byte_size) :: d1pw_ack4pot_la, d1ed_ack4pot_la, d1mf_ack4pot_la
    ! pairwise edensity embeddin pairwise_deriv1 edensity_deriv1 embeddin_deriv1
    print*, "testing linear approximation for functions, tries = ", rand_tries
    call random_seed()
    call initiate_energetic_parameters

    print*, "pairwise potential part"
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 6d0 + argshift
        function_value_an = pw_ack4pot(argument)
        function_value_la = pairwise(argument)
        absolute_error = function_value_la - function_value_an
        relative_error = 2d0 * absolute_error / ( function_value_la + function_value_an )
!        print 364, "arg = ", argument
!        print 364, "; absolute error = ", absolute_error
!        print 365, "; relative error =  ", relative_error*100
!        print *, "%"
        averelative_error = averelative_error + relative_error/rand_tries
        aveabsolute_error = aveabsolute_error + absolute_error/rand_tries
        if (abs(absolute_error) .gt. abs (maxabsolute_error)) then
            maxvalue = function_value_an
            maxargument = argument
            maxrelative_error = relative_error
            maxabsolute_error = absolute_error
        endif
    enddo
    print 364, "max arg = ", argument
    print 364, "max val = ", maxvalue
    print 364, "; max absolute error = ", maxabsolute_error
    print 364, "; max relative error =  ", maxrelative_error*100
    print 364, "; avg absolute error = ", aveabsolute_error
    print 364, "; avg relative error =  ", averelative_error*100
    print *, "%"

    maxrelative_error = 0d0 ; relative_error = 0d0 ; averelative_error = 0d0
    maxabsolute_error = 0d0 ; absolute_error = 0d0 ; aveabsolute_error = 0d0
    print*, "electronic density part"
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 4d0 + argshift
        function_value_an = ed_ack4pot(argument)
        function_value_la = edensity(argument)
        absolute_error = function_value_la - function_value_an
        relative_error = 2d0 * absolute_error / ( function_value_la + function_value_an )
!        if (( function_value_la + function_value_an ) .eq. 0d0) then
!            print*, function_value_la, function_value_an, argument
!        endif
        averelative_error = averelative_error + relative_error/rand_tries
        aveabsolute_error = aveabsolute_error + absolute_error/rand_tries
        if (abs(absolute_error) .gt. abs (maxabsolute_error)) then
            maxvalue = function_value_an
            maxargument = argument
            maxrelative_error = relative_error
            maxabsolute_error = absolute_error
        endif
    enddo
    print 364, "max arg = ", argument
    print 364, "max val = ", maxvalue
    print 364, "; max absolute error = ", maxabsolute_error
    print 364, "; max relative error =  ", maxrelative_error*100
    print 364, "; avg absolute error = ", aveabsolute_error
    print 364, "; avg relative error =  ", averelative_error*100
    print *, "%"

    maxrelative_error = 0d0 ; relative_error = 0d0 ; averelative_error = 0d0
    maxabsolute_error = 0d0 ; absolute_error = 0d0 ; aveabsolute_error = 0d0
    print*, "morphing function part"
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 25d0 + argshift
        function_value_an = mf_ack4pot(argument)
        function_value_la = embeddin(argument)
        absolute_error = function_value_la - function_value_an
        relative_error = 2d0 * absolute_error / ( function_value_la + function_value_an )
        averelative_error = averelative_error + relative_error/rand_tries
        aveabsolute_error = aveabsolute_error + absolute_error/rand_tries
        !if (abs(absolute_error) .gt. abs (maxabsolute_error)) then
        if (abs(relative_error) .gt. abs (maxrelative_error)) then
        !if (abs(function_value_an) .gt. abs (maxvalue)) then
            maxvalue = function_value_an
            maxargument = argument
            maxrelative_error = relative_error
            maxabsolute_error = absolute_error
        endif
    enddo
    print 364, "max arg = ", argument
    print 364, "max val = ", maxvalue
    print 364, "; max absolute error = ", maxabsolute_error
    print 364, "; max relative error =  ", maxrelative_error*100
    print 364, "; avg absolute error = ", aveabsolute_error
    print 364, "; avg relative error =  ", averelative_error*100
    print *, "%"
!    print*, "pairwise potential derivative part"
!
!    print*, "electronic density derivative part"
!
!    print*, "morphing function derivative part"


    361 format(A, I11.2, $)
    364 format(A, ES11.4, $)
    365 format(A, F14.9, $)
    367 format(/, A)
endsubroutine test_linear_approximation_of_potentials
