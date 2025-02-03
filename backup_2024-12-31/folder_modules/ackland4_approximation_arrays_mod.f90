module ackland4_approximation_arrays_mod
    implicit none
    save
    integer, parameter :: pot_steps = 100000!*100!from 00 A to 10 A. cutoff otherwise
    !for function   itself we only need its actual value
    !for first  derivative we both need its actual value and it's average value on step interval
    !for second derivative we only need it's average value on step interval
    real(8), parameter :: pairwise_step = (7d0/pot_steps), pairwise_recistep = 1d0/pairwise_step
    real(8), dimension(pot_steps) :: pairwise_pot, pairwise_dr1, pairwise_dr2
    real(8), dimension(pot_steps) ::               pairwise_ad1, pairwise_ad2
    !pairwise potential, argument step(reciprocal), first derivative, second derivative
    real(8), parameter :: elecdens_step = (7d0/pot_steps), elecdens_recistep = 1d0/elecdens_step
    real(8), dimension(pot_steps) :: elecdens_pot, elecdens_dr1, elecdens_dr2
    real(8), dimension(pot_steps) ::               elecdens_ad1, elecdens_ad2
    !electronic density, argument step(reciprocal), first derivative, second derivative
    real(8), parameter :: embefunc_step = (20d1/pot_steps), embefunc_recistep = 1d0/embefunc_step
    real(8), dimension(pot_steps) :: embefunc_pot, embefunc_dr1, embefunc_dr2
    real(8), dimension(pot_steps) ::               embefunc_ad1, embefunc_ad2
    !embedding function, argument step(reciprocal), first derivative, second derivative
endmodule ackland4_approximation_arrays_mod
