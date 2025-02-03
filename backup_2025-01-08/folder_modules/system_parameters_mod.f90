module system_parameters_mod
    use compile_debug_options_mod, only : float_byte_size, integer_byte_size
    use ackland4_parameters_mod, only : rphi_ack4pot, rpsi_ack4pot
    save
    real(float_byte_size), parameter :: a0bcc = 2.85573, a0fcc = 3.65216
    real(float_byte_size), parameter :: r_main_cell = 75d0 ! system is a cylinder with radius
    real(float_byte_size), parameter :: dr_casing_zone = 7.5d0 ! r_main_cell + dr_casing_zone + dr_elastic_zone
    real(float_byte_size), parameter :: dr_elastic_zone = 7.5d0
    real(float_byte_size), parameter :: h_main_cell = 9d0 ! system is a cylinder with height
    real(float_byte_size), parameter :: dh_periodic_zone = 7.5d0 ! h_main_cell + dh_periodic_zone
    real(float_byte_size), parameter :: dh_casing_zone = 7.5d0 ! + dh_casing_zone
    real(float_byte_size), parameter :: dh_elastic_zone = 7.5d0 ! + dh_elastic_zone
    real(float_byte_size), parameter :: cutoff_pw = rphi_ack4pot(15)
    real(float_byte_size), parameter :: cutoff_ed = rpsi_ack4pot(3)

endmodule system_parameters_mod
