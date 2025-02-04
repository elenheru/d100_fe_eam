subroutine test_wo_values_la
    use potentials_ackland4fefe_linear_mod
    implicit none
    integer(8), parameter :: rand_tries = 1000*100!0!*90!00
    real(float_byte_size), parameter :: argshift = 1d-1
    integer(8) :: i
    real(float_byte_size) :: argument, maxargument = 0d0
    real(float_byte_size) :: function_value_an
    real(float_byte_size) :: function_value_la
    real(float_byte_size) :: relative_error, maxrelative_error = 0d0, averelative_error = 0d0
    real(float_byte_size) :: absolute_error, maxabsolute_error = 0d0, aveabsolute_error = 0d0
    real(float_byte_size) :: pw_ack4pot, ed_ack4pot, mf_ack4pot
    real(float_byte_size) :: d1pw_ack4pot, d1ed_ack4pot, d1mf_ack4pot
    real(float_byte_size) :: pw_ack4pot_la, ed_ack4pot_la, mf_ack4pot_la
    real(float_byte_size) :: d1pw_ack4pot_la, d1ed_ack4pot_la, d1mf_ack4pot_la
    ! pairwise edensity embeddin pairwise_deriv1 edensity_deriv1 embeddin_deriv1
    print*, "testing wo linear approximation for functions, tries = ", rand_tries
    !stop "test other functions too"
    call random_seed()
    call initiate_energetic_parameters
    open (990, file = 'pw_randoms_simple.txt')
    open (991, file = 'ed_randoms_simple.txt')
    open (992, file = 'mf_randoms_simple.txt')
    open (995, file = 'pw_randoms_tricky.txt')

    print*, "pairwise potential part"
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 5d0 + argshift
        function_value_an = pw_ack4pot(argument)
        function_value_la = pairwise(argument)
        absolute_error = function_value_la - function_value_an
        relative_error = 2d0 * absolute_error / ( function_value_la + function_value_an )
        averelative_error = averelative_error + abs(relative_error)/rand_tries
        aveabsolute_error = aveabsolute_error + abs(absolute_error)/rand_tries
        if (abs(absolute_error) .gt. abs (maxabsolute_error)) then
            maxargument = argument
            maxrelative_error = relative_error
            maxabsolute_error = absolute_error
        endif
        write(990, 369) argument, function_value_an, function_value_la
    enddo
    print 364, "max arg = ", maxargument
    print 364, "; max absolute error = ", maxabsolute_error
    print 364, "; max relative error =  ", maxrelative_error*100
    print 364, "; avg absolute error = ", aveabsolute_error
    print 364, "; avg relative error =  ", averelative_error*100
    print *, "%"

    print*, "electronic density part"
    maxrelative_error = 0d0 ; maxabsolute_error = 0d0
    averelative_error = 0d0 ; aveabsolute_error = 0d0
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 4d0 + argshift
        function_value_an = ed_ack4pot(argument)
        function_value_la = edensity(argument)
        absolute_error = function_value_la - function_value_an
        relative_error = 2d0 * absolute_error / ( function_value_la + function_value_an )
        averelative_error = averelative_error + abs(relative_error)/rand_tries
        aveabsolute_error = aveabsolute_error + abs(absolute_error)/rand_tries
        if (abs(absolute_error) .gt. abs (maxabsolute_error)) then
            maxargument = argument
            maxrelative_error = relative_error
            maxabsolute_error = absolute_error
        endif
        write(991, 369) argument, function_value_an, function_value_la
    enddo
    print 364, "max arg = ", maxargument
    print 364, "; max absolute error = ", maxabsolute_error
    print 364, "; max relative error =  ", maxrelative_error*100
    print 364, "; avg absolute error = ", aveabsolute_error
    print 364, "; avg relative error =  ", averelative_error*100
    print *, "%"

    print*, "morphing function part"
    maxrelative_error = 0d0 ; maxabsolute_error = 0d0
    averelative_error = 0d0 ; aveabsolute_error = 0d0
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 25d0 + argshift
        function_value_an = mf_ack4pot(argument)
        function_value_la = embeddin(argument)
        absolute_error = function_value_la - function_value_an
        relative_error = 2d0 * absolute_error / ( function_value_la + function_value_an )
        averelative_error = averelative_error + abs(relative_error)/rand_tries
        aveabsolute_error = aveabsolute_error + abs(absolute_error)/rand_tries
        if (abs(absolute_error) .gt. abs (maxabsolute_error)) then
            maxargument = argument
            maxrelative_error = relative_error
            maxabsolute_error = absolute_error
!            if (abs(relative_error) .gt. 1d-2) then
!                print*, argument, function_value_la, function_value_an, "@@@"
!            endif
        endif
        write(992, 369) argument, function_value_an, function_value_la
    enddo
    print 364, "max arg = ", maxargument
    print 364, "; max absolute error = ", maxabsolute_error
    print 364, "; max relative error =  ", maxrelative_error*100
    print 364, "; avg absolute error = ", aveabsolute_error
    print 364, "; avg relative error =  ", averelative_error*100
    print *, "%"

    print*, "pairwise potential derivative part"
    maxrelative_error = 0d0 ; maxabsolute_error = 0d0
    averelative_error = 0d0 ; aveabsolute_error = 0d0
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 5d0 + argshift
        function_value_an = d1pw_ack4pot(argument)
        function_value_la = pairwise_deriv1(argument)
        absolute_error = function_value_la - function_value_an
        relative_error = 2d0 * absolute_error / ( function_value_la + function_value_an )
        averelative_error = averelative_error + abs(relative_error)/rand_tries
        aveabsolute_error = aveabsolute_error + abs(absolute_error)/rand_tries
        if (abs(absolute_error) .gt. abs (maxabsolute_error)) then
            maxargument = argument
            maxrelative_error = relative_error
            maxabsolute_error = absolute_error
        endif
        write(990, 369) argument, function_value_an, function_value_la
    enddo
    print 364, "max arg = ", argument
    print 364, "; max absolute error = ", maxabsolute_error
    print 364, "; max relative error =  ", maxrelative_error*100
    print 364, "; avg absolute error = ", aveabsolute_error
    print 364, "; avg relative error =  ", averelative_error*100
    print *, "%"

    print*, "electronic density derivative part"
    maxrelative_error = 0d0 ; maxabsolute_error = 0d0
    averelative_error = 0d0 ; aveabsolute_error = 0d0
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 4d0 + argshift
        function_value_an = d1ed_ack4pot(argument)
        function_value_la = edensity_deriv1(argument)
        absolute_error = function_value_la - function_value_an
        relative_error = 2d0 * absolute_error / ( function_value_la + function_value_an )
        averelative_error = averelative_error + abs(relative_error)/rand_tries
        aveabsolute_error = aveabsolute_error + abs(absolute_error)/rand_tries
        if (abs(absolute_error) .gt. abs (maxabsolute_error)) then
            maxargument = argument
            maxrelative_error = relative_error
            maxabsolute_error = absolute_error
        endif
        write(991, 369) argument, function_value_an, function_value_la
    enddo
    print 364, "max arg = ", argument
    print 364, "; max absolute error = ", maxabsolute_error
    print 364, "; max relative error =  ", maxrelative_error*100
    print 364, "; avg absolute error = ", aveabsolute_error
    print 364, "; avg relative error =  ", averelative_error*100
    print *, "%"

    print*, "morphing function derivative part"
    maxrelative_error = 0d0 ; maxabsolute_error = 0d0
    averelative_error = 0d0 ; aveabsolute_error = 0d0
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 25d0 + argshift
        function_value_an = d1mf_ack4pot(argument)
        function_value_la = embeddin_deriv1(argument)
        absolute_error = function_value_la - function_value_an
        relative_error = 2d0 * absolute_error / ( function_value_la + function_value_an )
        averelative_error = averelative_error + abs(relative_error)/rand_tries
        aveabsolute_error = aveabsolute_error + abs(absolute_error)/rand_tries
        if (abs(relative_error) .gt. abs (maxrelative_error)) then
            maxargument = argument
            maxrelative_error = relative_error
            maxabsolute_error = absolute_error
!            if (abs(relative_error) .gt. 1d-2) then
!                print*, argument, function_value_la, function_value_an, "@@@"
!            endif
        endif
        write(992, 369) argument, function_value_an, function_value_la
    enddo
    print 364, "max arg = ", maxargument
    print 364, "; max absolute error = ", maxabsolute_error
    print 364, "; max relative error =  ", maxrelative_error*100
    print 364, "; avg absolute error = ", aveabsolute_error
    print 364, "; avg relative error =  ", averelative_error*100
    print *, "%"

    close(990); close(991); close(992)

    361 format(A, I11.2, $)
    364 format(A, ES11.4, $)
    365 format(A, F14.9, $)
    367 format(/, A)
    369 format(3(ES14.7, 1x))
endsubroutine test_wo_values_la
