real(float_byte_size) function pw_ack4pot(r)
    use ackland4_parameters_mod
    implicit none
    real(float_byte_size), intent(in) :: r
    real(float_byte_size), parameter :: rb = 0.529177211, qe = sqrt(14.399764), zfe = 26d0
    real(float_byte_size) :: intermediate, rs
    !pw = pairwise term
    integer(integer_byte_size) i

    intermediate = 0d0
    do i=1, 15
        intermediate = intermediate + &
        aphi_ack4pot(i) * max( sign( 1d0, rphi_ack4pot(i)-r ), 0d0 ) * (rphi_ack4pot(i)-r)**3
    enddo
    pw_ack4pot = 5d-1*intermediate
    return
    ! according to article, V(r) must be calculated as follows, but it does not work
    if ( r .lt. rphi_ack4pot(1) ) then
        rs = 0.88534 * rb / ( sqrt(2d0) * zfe**(1d0/3d0) )
        pw_ack4pot = zfe**2 * qe**2 * varphi_ack4pot(r/rs) / r
    elseif (( r .ge. rphi_ack4pot(1) ) .and. ( r .lt. rphi_ack4pot(2) )) then
        pw_ack4pot = exp( ((bbi_ack4pot(3)*r + bbi_ack4pot(2))*r + bbi_ack4pot(1))*r + bbi_ack4pot(0) )
    else
        intermediate = 0d0
        do i=2, 15
            intermediate = intermediate + &
            aphi_ack4pot(i) * max( sign( 1d0, rphi_ack4pot(i)-r ), 0d0 ) * (rphi_ack4pot(i)-r)**3
        enddo
        pw_ack4pot = intermediate
    endif
endfunction pw_ack4pot
