module verlet_list_mod
    use compile_debug_options_mod
    use system_parameters_mod, only : cutoff_pw
    use positions_mod, only : max_at
    save
    real(float_byte_size), parameter :: verlet_cutoff = cutoff_pw + 0.1 ! see potential parameters
    integer(integer_byte_size), parameter :: maximum_neibors = 200
    integer(integer_byte_size) :: verlet_list_length(max_at)
    integer(integer_byte_size) :: verlet_list(max_at, maximum_neibors)
endmodule verlet_list_mod
