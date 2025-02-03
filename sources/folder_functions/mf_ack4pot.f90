real(float_byte_size) function mf_ack4pot(rho)
    use ackland4_parameters_mod
    implicit none
    real(float_byte_size), intent(in) :: rho
    !mf = embedding function
    mf_ack4pot = -sqrt(rho) + aphibi_ack4pot * rho * rho
endfunction mf_ack4pot
