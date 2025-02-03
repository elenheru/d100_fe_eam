module physical_parameters_mod
    save

    real(8), parameter :: iron_mass_amu = 55.9349363
    real(8), parameter :: ev_in_si = 1.602176634d-19
    real(8), parameter :: angstrem_in_si = 1d-10
    real(8), parameter :: amu_in_kg = 1.660539d-27
    real(8), parameter :: timestep_sec = 1000*1000*1d-18 !1d-12 is appr time of one atomic oscillation
endmodule physical_parameters_mod
