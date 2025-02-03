module ackland4_parameters_mod
    use compile_debug_options_mod, only : float_byte_size, integer_byte_size
    save
    real(float_byte_size), parameter, dimension(2) :: rsm_ack4pot = (/ 0.9, 1.95/)
    real(float_byte_size), parameter, dimension(0:3) :: bbi_ack4pot = &
        (/ 14.996917289290, -20.533174190155, 14.002591780752, -3.6473736591143/)
    real(float_byte_size), parameter, dimension(1:15) :: rphi_ack4pot = &
        (/ 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 3.3, 3.7, 4.2, 4.7, 5.3, 6.0/)
    real(float_byte_size), parameter, dimension(1:15) :: aphi_ack4pot = &
        (/ 195.92322853994, 17.516698453315, 1.4926525164290      &
        , 6.4129476125197, -6.8157461860553, 9.6582581963600      &
        , -5.3419002764419, 1.7996558048346, -1.4788966636288      &
        , 1.8530435283665, -0.64164344859316, 0.24463630025168     &
        , -0.057721650527383, 0.023358616514826, -0.0097064921265079 /)
    real(float_byte_size), parameter, dimension(1:3) :: rpsi_ack4pot = (/ 2.4, 3.2, 4.2/)
    real(float_byte_size), parameter, dimension(1:3) :: apsi_ack4pot = &
        (/  11.686859407970, -0.014710740098830, 0.47193527075943 /)
    real(float_byte_size), parameter :: aphibi_ack4pot = -0.00034906178363530

    contains
        real(float_byte_size) function varphi_ack4pot(x)
            implicit none
            real(float_byte_size), intent(in) :: x
            varphi_ack4pot = 0.1818 * exp( -3.2000 * x) + 0.5099 * exp(-0.9423 * x) &
                           + 0.2802 * exp( -0.4029 * x) + 0.02817* exp(-0.2016 * x)
        endfunction
        real(float_byte_size) function d1varphi_ack4pot(x)
            implicit none
            real(float_byte_size), intent(in) :: x
            d1varphi_ack4pot = -3.2000 * 0.1818 * exp( -3.2000 * x) &
                            -0.9423 * 0.5099 * exp(-0.9423 * x) &
                            -0.4029 * 0.2802 * exp( -0.4029 * x) &
                            -0.2016 * 0.02817* exp(-0.2016 * x)
        endfunction

endmodule ackland4_parameters_mod
