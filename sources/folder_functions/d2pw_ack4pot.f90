real(float_byte_size) function d2pw_ack4pot(r)
    use ackland4_parameters_mod
    implicit none
    real(float_byte_size), intent(in) :: r
    real(float_byte_size), parameter :: rb = 0.529177211, qe = sqrt(14.399764), zfe = 26e0
    real(float_byte_size) :: intermediate, rs
    !pw = pairwise term
    integer(integer_byte_size) i
    intermediate = 0e0
    do i=1, 15
        intermediate = intermediate + &
        aphi_ack4pot(i) * max( sign( 1e0_float_byte_size, rphi_ack4pot(i)-r ), 0e0 ) * 6e0 * (rphi_ack4pot(i)-r)
    enddo
    d2pw_ack4pot = 5e-1*intermediate
endfunction d2pw_ack4pot
