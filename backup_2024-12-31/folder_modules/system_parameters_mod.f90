module system_parameters_mod
    use ackland4_parameters_mod, only : rphi_ack4pot, rpsi_ack4pot
    save
    real(8), parameter :: a0bcc = 2.85573, a0fcc = 3.65216
    real(8), parameter :: r_main_cell = 16d0 ! system is a sphere with radius =
    real(8), parameter :: dr_elastic_zone = 8.5d0 ! r_main_cell + dr_elastic_zone
    real(8), parameter :: cutoff_pw = rphi_ack4pot(15)
    real(8), parameter :: cutoff_ed = rpsi_ack4pot(3)

endmodule system_parameters_mod
