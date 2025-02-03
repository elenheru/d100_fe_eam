real(float_byte_size) function d1mf_ack4pot(rho)
    use ackland4_parameters_mod
    implicit none
    real(float_byte_size), intent(in) :: rho
    !mf = embedding function
    d1mf_ack4pot = (-5d-1 / sqrt(rho) ) + 2 * aphibi_ack4pot * rho
endfunction d1mf_ack4pot
