real(float_byte_size) function d2mf_ack4pot(rho)
    use ackland4_parameters_mod
    implicit none
    real(float_byte_size), intent(in) :: rho
    !mf = embedding function
    d2mf_ack4pot = ( 25e-2 / (rho*sqrt(rho)) ) + 2 * aphibi_ack4pot
endfunction d2mf_ack4pot
