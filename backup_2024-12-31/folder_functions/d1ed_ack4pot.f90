real(8) function d1ed_ack4pot(r)
    use ackland4_parameters_mod
    implicit none
    real(8), intent(in) :: r
    real(8) :: intermediate
    !ed = electronic density
    integer i
    intermediate = 0
    do i=1,3
        intermediate = intermediate + &
        apsi_ack4pot(i) * max( sign( 1d0, rpsi_ack4pot(i)-r ), 0d0 ) * -3d0 * (rpsi_ack4pot(i)-r)**2
    enddo
    d1ed_ack4pot = intermediate
endfunction d1ed_ack4pot
