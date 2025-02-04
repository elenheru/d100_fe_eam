subroutine test_derivatives
    use ackland4_parameters_mod
    implicit none
    integer(integer_byte_size), parameter :: rand_tries = 10
    real(float_byte_size), parameter :: argstep = 1d-4
    integer(integer_byte_size) :: i
    real(float_byte_size) :: argument
    real(float_byte_size) :: derivative_value_at_the_point
    real(float_byte_size) :: function_value_at_previous_point
    real(float_byte_size) :: function_value_at_next_point
    real(float_byte_size) :: derivative_value_approximation
    real(float_byte_size) :: relative_error
    real(float_byte_size) :: absolute_error
    real(float_byte_size) :: pw_ack4pot, ed_ack4pot, mf_ack4pot
    real(float_byte_size) :: d1pw_ack4pot, d1ed_ack4pot, d1mf_ack4pot
    real(float_byte_size) :: d2pw_ack4pot, d2ed_ack4pot, d2mf_ack4pot

    print*, "testing first and second derivatives via symmetric differential approximation"
    call random_seed()

    print*, "pairwise potential part"
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 2d0 + argstep
        derivative_value_at_the_point = d1pw_ack4pot(argument)
        function_value_at_previous_point = pw_ack4pot(argument - argstep)
        function_value_at_next_point = pw_ack4pot(argument + argstep)
        derivative_value_approximation = &
        (function_value_at_next_point - function_value_at_previous_point)/(2d0 * argstep)
        absolute_error = derivative_value_approximation - derivative_value_at_the_point
        relative_error = 2d0 * absolute_error / (derivative_value_approximation + derivative_value_at_the_point)
        print 354, "arg = ", argument
        print 354, "; step =  ", argstep
        print 354, "; absolute error = ", absolute_error
        print 354, "; relative error =  ", relative_error
        print *, "%"
    enddo
    print*, "electronic density part"
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 2d0 + argstep
        derivative_value_at_the_point = d1ed_ack4pot(argument)
        function_value_at_previous_point = ed_ack4pot(argument - argstep)
        function_value_at_next_point = ed_ack4pot(argument + argstep)
        derivative_value_approximation = &
        (function_value_at_next_point - function_value_at_previous_point)/(2d0 * argstep)
        absolute_error = derivative_value_approximation - derivative_value_at_the_point
        relative_error = 2d0 * absolute_error / (derivative_value_approximation + derivative_value_at_the_point)
        print 354, "arg = ", argument
        print 354, "; step =  ", argstep
        print 354, "; absolute error = ", absolute_error
        print 354, "; relative error =  ", relative_error
        print *, "%"
    enddo
    print*, "morphing function part"
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 2d0 + argstep
        derivative_value_at_the_point = d1mf_ack4pot(argument)
        function_value_at_previous_point = mf_ack4pot(argument - argstep)
        function_value_at_next_point = mf_ack4pot(argument + argstep)
        derivative_value_approximation = &
        (function_value_at_next_point - function_value_at_previous_point)/(2d0 * argstep)
        absolute_error = derivative_value_approximation - derivative_value_at_the_point
        relative_error = 2d0 * absolute_error / (derivative_value_approximation + derivative_value_at_the_point)
        print 354, "arg = ", argument
        print 354, "; step =  ", argstep
        print 354, "; absolute error = ", absolute_error
        print 354, "; relative error =  ", relative_error
        print *, "%"
    enddo
    print*, "pairwise potential derivative part"
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 2d0 + argstep
        derivative_value_at_the_point = d2pw_ack4pot(argument)
        function_value_at_previous_point = d1pw_ack4pot(argument - argstep)
        function_value_at_next_point = d1pw_ack4pot(argument + argstep)
        derivative_value_approximation = &
        (function_value_at_next_point - function_value_at_previous_point)/(2d0 * argstep)
        absolute_error = derivative_value_approximation - derivative_value_at_the_point
        relative_error = 2d0 * absolute_error / (derivative_value_approximation + derivative_value_at_the_point)
        print 354, "arg = ", argument
        print 354, "; step =  ", argstep
        print 354, "; absolute error = ", absolute_error
        print 354, "; relative error =  ", relative_error
        print *, "%"
    enddo
    print*, "electronic density derivative part"
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 2d0 + argstep
        derivative_value_at_the_point = d2ed_ack4pot(argument)
        function_value_at_previous_point = d1ed_ack4pot(argument - argstep)
        function_value_at_next_point = d1ed_ack4pot(argument + argstep)
        derivative_value_approximation = &
        (function_value_at_next_point - function_value_at_previous_point)/(2d0 * argstep)
        absolute_error = derivative_value_approximation - derivative_value_at_the_point
        relative_error = 2d0 * absolute_error / (derivative_value_approximation + derivative_value_at_the_point)
        print 354, "arg = ", argument
        print 354, "; step =  ", argstep
        print 354, "; absolute error = ", absolute_error
        print 354, "; relative error =  ", relative_error
        print *, "%"
    enddo
    print*, "morphing function derivative part"
    do i = 1, rand_tries
        call random_number(argument)
        argument = argument * 2d0 + argstep
        derivative_value_at_the_point = d2mf_ack4pot(argument)
        function_value_at_previous_point = d1mf_ack4pot(argument - argstep)
        function_value_at_next_point = d1mf_ack4pot(argument + argstep)
        derivative_value_approximation = &
        (function_value_at_next_point - function_value_at_previous_point)/(2d0 * argstep)
        absolute_error = derivative_value_approximation - derivative_value_at_the_point
        relative_error = 2d0 * absolute_error / (derivative_value_approximation + derivative_value_at_the_point)
        print 354, "arg = ", argument
        print 354, "; step =  ", argstep
        print 354, "; absolute error = ", absolute_error
        print 354, "; relative error =  ", relative_error
        print *, "%"
    enddo

    351 format(A, I11.2, $)
    354 format(A, ES11.4, $)
    357 format(/, A)
endsubroutine test_derivatives
