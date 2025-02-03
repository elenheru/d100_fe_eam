real(8) function mf_ack4pot(rho)
    use ackland4_parameters_mod
    implicit none
    real(8), intent(in) :: rho
    !mf = embedding function
    mf_ack4pot = -sqrt(rho) + aphibi_ack4pot * rho * rho
endfunction mf_ack4pot
