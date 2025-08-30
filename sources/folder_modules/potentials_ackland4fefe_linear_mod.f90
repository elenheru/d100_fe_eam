module      potentials_ackland4fefe_linear_mod
    use compile_debug_options_mod, only : float_byte_size, integer_byte_size
    use ackland4_approximation_arrays_mod
    contains

    elemental real(float_byte_size) function pairwise(r)
        implicit none
        real(float_byte_size), intent(in) :: r;real(float_byte_size) dist
        integer(integer_byte_size) intdistance
        pairwise = 0e0
        dist = r*pairwise_recistep
        intdistance = floor(dist)
        if (intdistance .gt. pot_steps) intdistance = pot_steps - 1
        if (intdistance .le. 1) intdistance = 1
        pairwise =  pairwise_pot(intdistance) + pairwise_ad1(intdistance) * &
                    (r - pairwise_step*intdistance)
    endfunction pairwise

    elemental real(float_byte_size) function edensity(r)
        implicit none
        real(float_byte_size), intent(in) :: r
        real(float_byte_size) dist
        integer(integer_byte_size) intdistance
        edensity = 0e0
        dist = r*elecdens_recistep
        intdistance = floor(dist)
        if (intdistance .gt. pot_steps) intdistance = pot_steps - 1
        if (intdistance .le. 0        ) intdistance = 1
        edensity =  elecdens_pot(intdistance) + elecdens_ad1(intdistance) * &
                    (r - elecdens_step*intdistance)
    endfunction edensity

    elemental real(float_byte_size) function embeddin(rho)
        use ackland4_approximation_arrays_mod
        implicit none
        real(float_byte_size), intent(in) :: rho
        integer(integer_byte_size) rho_knot
        embeddin = 0e0
        rho_knot = floor(rho * embefunc_recistep)
        if ( rho_knot .le. 0 ) rho_knot = 1
            if ( rho_knot .le. 0 ) rho_knot = 1
        if ( rho_knot .ge. pot_steps) rho_knot = pot_steps - 1
        embeddin =  embefunc_pot(rho_knot) + embefunc_ad1(rho_knot) * &
                    (rho - embefunc_step * rho_knot)
    endfunction embeddin

    elemental real(float_byte_size) function pairwise_deriv1(r)
        use ackland4_approximation_arrays_mod
        implicit none
        real(float_byte_size), intent(in) :: r
        real(float_byte_size) dist
        integer(integer_byte_size) intdistance
        pairwise_deriv1= 0e0
        dist = r*pairwise_recistep
        intdistance = floor(dist)
        if (intdistance .ge. pot_steps) return
        if (intdistance .le. 0        ) intdistance = 1
        pairwise_deriv1 =  pairwise_dr1(intdistance) + pairwise_ad2(intdistance) * &
                    (r - pairwise_step*intdistance)
    endfunction pairwise_deriv1

    elemental real(float_byte_size) function edensity_deriv1(r)
        use ackland4_approximation_arrays_mod
        implicit none
        real(float_byte_size), intent(in) :: r;real(float_byte_size) dist
        integer(integer_byte_size) intdistance
        edensity_deriv1 = 0e0
        dist = r*elecdens_recistep
        intdistance = floor(dist)
        if (intdistance .ge. pot_steps) return
        if (intdistance .le. 0) intdistance = 1
        edensity_deriv1 =  elecdens_dr1(intdistance) + elecdens_ad2(intdistance) * &
                    (r - elecdens_step*intdistance)
    endfunction edensity_deriv1

    elemental real(float_byte_size) function embeddin_deriv1(rho)
        use ackland4_approximation_arrays_mod
        implicit none
        real(float_byte_size), intent(in) :: rho
        integer(integer_byte_size) rho_knot
        embeddin_deriv1 = 0e0
        if(rho .lt. 0e0) return
        rho_knot = floor(rho*embefunc_recistep)
        if ( rho_knot .le. 0 ) rho_knot = 1
            ! dangerous errors for small
            !if ( rho_knot .le. 0 ) rho_knot = 1
        if ( rho_knot .ge. pot_steps) rho_knot = pot_steps - 1
        embeddin_deriv1 =  embefunc_dr1(rho_knot) + embefunc_ad2(rho_knot) * &
                    (rho - embefunc_step * rho_knot)
    endfunction embeddin_deriv1

endmodule   potentials_ackland4fefe_linear_mod
