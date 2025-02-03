subroutine allocate_electronic_densities
    use compile_debug_options_mod, only : debug_flag
    use electronic_densities_allocated_mod
    use positions_mod, only : last_wall

    if ( debug_flag ) print *, &
        "Allocating electronic densities array (known as RHO_i)"
    if ( allocated(rho_at) ) then
        STOP "Trying to allocate array of rho's but it is already allocated"
    else
        allocate(rho_at(last_wall))
    endif
endsubroutine allocate_electronic_densities
