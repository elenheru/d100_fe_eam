module      potentials_ackland4fefe_linear_mod
    contains
        !pairwise potential using  initiated massive
    elemental real(8) function pairwise(r)
        use ackland4_approximation_arrays_mod
        implicit none
        real(8),intent(in)::r;real(8) dist
        integer intdistance
        pairwise = 0d0
        dist = r*pairwise_recistep
        intdistance = floor(dist)
        if (intdistance .gt. pot_steps) intdistance = pot_steps - 1 !probably need fair function near cutoff
        if (intdistance .le. 1        ) intdistance = 1
        pairwise =  pairwise_pot(intdistance) + pairwise_ad1(intdistance) * &
                    (r - pairwise_step*intdistance)
        !    pairwise =  pairwise_pot(intdistance) + pairwise_dr1(intdistance) * &
        !                (r - pairwise_step*intdistance)
        !significant error when too close to the next node, because of non-linearity
    endfunction pairwise
        !electronic density using  initiated massive
    elemental real(8) function edensity(r)
        use ackland4_approximation_arrays_mod
        implicit none
        real(8),intent(in)::r;real(8) dist
        integer intdistance
        edensity = 0d0
        dist = r*elecdens_recistep
        intdistance = floor(dist)
        if (intdistance .gt. pot_steps) intdistance = pot_steps - 1 !probably need fair function near cutoff
        if (intdistance .le. 0        ) intdistance = 1
        edensity =  elecdens_pot(intdistance) + elecdens_ad1(intdistance) * &
                    (r - elecdens_step*intdistance)!? check
    endfunction edensity
        !embedding function using  initiated massive
    elemental real(8) function embeddin(rho)
        use ackland4_approximation_arrays_mod !of course argument is not really a distance
        implicit none
        real(8),intent(in)::rho;real(8) dist
        integer intdistance
        embeddin = 0d0
        dist = rho*embefunc_recistep
        intdistance = floor(dist)
        if ( intdistance .le. 0 ) intdistance = 1
        if ( intdistance .ge. pot_steps) intdistance = pot_steps - 1
        embeddin =  embefunc_pot(intdistance) + embefunc_ad1(intdistance) * &
                    (rho - embefunc_step*intdistance)
    endfunction embeddin
    !pairwise potential using  initiated massive
    elemental real(8) function pairwise_deriv1(r)
        use ackland4_approximation_arrays_mod
        implicit none
        real(8),intent(in)::r;real(8) dist
        integer intdistance
        pairwise_deriv1= 0d0
        dist = r*pairwise_recistep
        intdistance = floor(dist)
        if (intdistance .ge. pot_steps) return
        if (intdistance .le. 0        ) intdistance = 1
        pairwise_deriv1 =  pairwise_dr1(intdistance) + pairwise_ad2(intdistance) * &
                    (r - pairwise_step*intdistance)
    endfunction pairwise_deriv1
    !electronic density using  initiated massive
    elemental real(8) function edensity_deriv1(r)
        use ackland4_approximation_arrays_mod
        implicit none
        real(8),intent(in)::r;real(8) dist
        integer intdistance
        edensity_deriv1 = 0d0      !;RETURN!to check if III Newton's law is correct here
        dist = r*elecdens_recistep !        surprisingly it is incorrect.
        intdistance = floor(dist)
        if (intdistance .ge. pot_steps) return
        if (intdistance .le. 0        ) intdistance = 1
        edensity_deriv1 =  elecdens_dr1(intdistance) + elecdens_ad2(intdistance) * &
                    (r - elecdens_step*intdistance)
    endfunction edensity_deriv1
    !embedding function using  initiated massive
    elemental real(8) function embeddin_deriv1(rho)
        use ackland4_approximation_arrays_mod !of course argument is not really a distance
        implicit none
        real(8),intent(in)::rho; real(8) dist
        integer intdistance
        embeddin_deriv1 = 0d0
        if(rho .lt. 0d0) return
        dist = rho*embefunc_recistep
        intdistance = floor(dist)
        if ( intdistance .le. 0 ) intdistance = 1
        if ( intdistance .ge. pot_steps) intdistance = pot_steps - 1
        embeddin_deriv1 =  embefunc_dr1(intdistance) + embefunc_ad2(intdistance) * &
                    (rho - embefunc_step*intdistance)
    endfunction embeddin_deriv1

endmodule   potentials_ackland4fefe_linear_mod
