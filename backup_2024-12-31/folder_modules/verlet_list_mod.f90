module verlet_list_mod
    use positions_mod, only : max_at
    save
    real(8), parameter :: verlet_cutoff = 6.1d0 ! see potential parameters
    integer, parameter :: maximum_neibors = 200
    integer :: verlet_list_length(max_at), verlet_list(max_at, maximum_neibors)
endmodule verlet_list_mod
