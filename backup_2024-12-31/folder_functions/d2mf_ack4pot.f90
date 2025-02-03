real(8) function d2mf_ack4pot(rho)
    use ackland4_parameters_mod
    implicit none
    real(8), intent(in) :: rho
    !mf = embedding function
    d2mf_ack4pot = ( 25d-2 / (rho*sqrt(rho)) ) + 2 * aphibi_ack4pot
endfunction d2mf_ack4pot
