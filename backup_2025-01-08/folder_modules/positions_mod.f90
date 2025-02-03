module positions_mod
    use compile_debug_options_mod, only : float_byte_size, integer_byte_size
    save
    integer(integer_byte_size), parameter :: max_at = 170000
    real(float_byte_size), dimension(3, max_at) :: R_at
    !integer(integer_byte_size) :: last_wall, last_relaxable !elastic close far and main cell
    integer(integer_byte_size) :: first_relaxable, last_relaxable ! main_cell movable atoms
    integer(integer_byte_size) :: first_periodic, last_periodic ! main_cell anchored atoms
    integer(integer_byte_size) :: first_casing, last_casing ! main_cell anchored atoms
    integer(integer_byte_size) :: first_wall, last_wall ! wall_
    ! Z order from inside:
    ! Relaxable - Periodic - Casing - Wall
    ! R order from inside
    ! Relaxable - Casing - Wall
    ! This partition scheme seems to be inevitable, as far we reject mere periodicity
    ! If you simply increase thickness of wall,
    ! then you either incorrectly calculate system energy
    ! or have to add a force for zones interface.
endmodule positions_mod
