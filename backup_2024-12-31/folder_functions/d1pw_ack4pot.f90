real(8) function d1pw_ack4pot(r)
    use ackland4_parameters_mod
    implicit none
    real(8), intent(in) :: r
    real(8), parameter  :: rb = 0.529177211, qe = sqrt(14.399764), zfe = 26d0
    real(8) :: intermediate, rs
    !pw = pairwise term
    integer i
    intermediate = 0d0
    do i=1,15
        intermediate = intermediate + &
        aphi_ack4pot(i) * max( sign( 1d0, rphi_ack4pot(i)-r ), 0d0 ) * -3d0 * (rphi_ack4pot(i)-r)**2
    enddo
    d1pw_ack4pot = 5d-1*intermediate
endfunction d1pw_ack4pot
