real(float_byte_size) function ed_ack4pot(r)
    use ackland4_parameters_mod
    implicit none
    real(float_byte_size), intent(in) :: r
    real(float_byte_size) :: intermediate
    !ed = electronic density
    integer(integer_byte_size) i
    intermediate = 0
    do i=1, 3
        intermediate = intermediate + &
        apsi_ack4pot(i) * max( sign( 1d0, rpsi_ack4pot(i)-r ), 0d0 ) * (rpsi_ack4pot(i)-r)**3
    enddo
    ed_ack4pot = intermediate
endfunction ed_ack4pot
