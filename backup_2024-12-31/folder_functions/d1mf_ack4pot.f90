real(8) function d1mf_ack4pot(rho)
    use ackland4_parameters_mod
    implicit none
    real(8), intent(in) :: rho
    !mf = embedding function
    d1mf_ack4pot = (-5d-1 / sqrt(rho) ) + 2 * aphibi_ack4pot * rho
endfunction d1mf_ack4pot
