module compile_debug_options_mod
    save
    integer, parameter :: float_byte_size = 8
    integer, parameter :: integer_byte_size = 4
    integer, parameter :: atom_of_interest = 370
    logical, parameter :: debug_flag = .not. .true.
    character(LEN = *), parameter :: debug_string_format = "('DBG: ', A, 1x, $)"
    integer(integer_byte_size), dimension(10)  :: timers
endmodule compile_debug_options_mod
