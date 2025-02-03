subroutine test_verlet_list
    use positions_mod
    use potentials_ackland4fefe_linear_mod
    use verlet_list_mod
    use electronic_densities_allocated_mod
    use compile_debug_options_mod, only : timers

    implicit none
    real(float_byte_size), parameter :: shift_magnitude = 1d-3
    integer(integer_byte_size) i, j, iv, jv, nprobe
    real(float_byte_size) :: dist, float_random, cosine_v
    real(float_byte_size) :: energy_before_shift
    real(float_byte_size) :: energy_after_shift
    real(float_byte_size) :: energy_before_shift_verlet
    real(float_byte_size) :: energy_after_shift_verlet

    real(float_byte_size), dimension(3) :: grad_num_t1, grad_num_t2, grad_num_t3, grad_num_t4, grad_num
    real(float_byte_size), dimension(3) :: shift, move

    print*, "Testing verlet list for atoms "
    print*, "test procedure includes calculation of system energy in two ways"
    print*, "without verlet list and with it"

    !print"(A)", "calculating force via numerical derivative formula, using verlet list"

    call calculate_system_energy(energy_before_shift)
    print *, "System energy (direct calculation) before shift is: ", energy_before_shift
    call calculate_system_energy_verlet(energy_before_shift_verlet)
    print *, "System energy (verlet calculation) before shift is: ", energy_before_shift_verlet

    call random_number(float_random)
    nprobe = floor(float_random * last_relaxable)
    print*, nprobe, "is atom probe number, which we are going to shift lightly in random direction"

    call random_number(shift)
    shift = ( 2d0*shift - (/1d0, 1d0, 1d0/) ) * shift_magnitude
    print 301, "shift is ", shift

    R_at(:, nprobe) = R_at(:, nprobe) + shift

    call system_clock(timers(1))
    call calculate_system_energy(energy_after_shift)
    call system_clock(timers(2))
    print *, "System energy (direct calculation) after shift is: ", energy_after_shift
    call system_clock(timers(3))
    call calculate_system_energy_verlet(energy_after_shift_verlet)
    call system_clock(timers(4))
    print *, "System energy (verlet calculation) after shift is: ", energy_after_shift_verlet

    R_at(:, nprobe) = R_at(:, nprobe) - shift
    print *, "Time spent for direct energy calculation ", timers(2) - timers(1)
    print *, "Time spent for verlet energy calculation ", timers(4) - timers(3)
301 format(A, 3(1x, F12.8))

endsubroutine test_verlet_list
