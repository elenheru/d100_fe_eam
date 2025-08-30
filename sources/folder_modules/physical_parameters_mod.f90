module physical_parameters_mod
    use compile_debug_options_mod, only : float_byte_size, integer_byte_size
    save
    real(float_byte_size), parameter :: iron_mass_amu = 55.9349363
    real(float_byte_size), parameter :: ev_in_si = 1.602176634e-19
    real(float_byte_size), parameter :: angstrem_in_si = 1e-10
    real(float_byte_size), parameter :: amu_in_kg = 1.660539e-27
    real(float_byte_size), parameter :: timestep_sec = 1000*1000*1e-18 !1d-12 is appr time of one atomic oscillation
endmodule physical_parameters_mod
