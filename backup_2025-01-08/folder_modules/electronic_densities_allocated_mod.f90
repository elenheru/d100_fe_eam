module electronic_densities_allocated_mod
    use compile_debug_options_mod, only : float_byte_size
    save
    real(float_byte_size), allocatable :: rho_at(:)
endmodule electronic_densities_allocated_mod


