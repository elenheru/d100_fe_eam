module rotations_mod
    use compile_debug_options_mod
    save
    real(float_byte_size), parameter, dimension(3,3) :: &
        matrix_unit = reshape( [real(float_byte_size) :: &
            1.0, 0.0, 0.0, &
            0.0, 1.0, 0.0, &
            0.0, 0.0, 1.0 ], (/3, 3/))
    real(float_byte_size), parameter :: &
            pi = asin(1.0_float_byte_size) * 2.0_float_byte_size
    real(float_byte_size), parameter :: &
        degree_110_to_111 = &
        ( asin(sqrt(1.0_float_byte_size/3.0_float_byte_size)) &
        * 180.0_float_byte_size / pi)
!    real(float_byte_size), parameter, dimension(3,3) :: &
!        rm_t_pnp_b_pp0 =
    ! rotation from default orientation to
    ! b || [1 1 0], t || (1 -1 1)

    contains

    subroutine build_rotation_matrix(matrix, axis, degrees)
        character(LEN=1), intent(in) :: axis
        real(float_byte_size), dimension(3,3), intent (inout) :: matrix
        real(float_byte_size), intent (in) :: degrees
        real(float_byte_size), dimension(3,3) :: matrix_requested
        real(float_byte_size) :: radians

        radians = pi * degrees / 180.0_float_byte_size
        matrix_requested = 0
        ! use invalid axis for identity matrix generation
        ! dont put junk into valid axes
        if (    (axis .eq. "x") .or. (axis .eq. "X")) then
            matrix_requested(1,1) = 1.0_float_byte_size
            matrix_requested(2,2) = cos(radians)
            matrix_requested(3,3) = cos(radians)

            matrix_requested(2,3) = sin(radians)
            matrix_requested(3,2) =-sin(radians)

        elseif ((axis .eq. "y") .or. (axis .eq. "Y")) then
            matrix_requested(1,1) = cos(radians)
            matrix_requested(2,2) = 1.0_float_byte_size
            matrix_requested(3,3) = cos(radians)

            matrix_requested(1,3) =-sin(radians)
            matrix_requested(3,1) = sin(radians)

        elseif ((axis .eq. "z") .or. (axis .eq. "Z")) then
            matrix_requested(1,1) = cos(radians)
            matrix_requested(2,2) = cos(radians)
            matrix_requested(3,3) = 1.0_float_byte_size

            matrix_requested(1,2) = sin(radians)
            matrix_requested(2,1) =-sin(radians)

        else
            matrix_requested(1,1) = 1.0_float_byte_size
            matrix_requested(2,2) = 1.0_float_byte_size
            matrix_requested(3,3) = 1.0_float_byte_size
            matrix = matrix_requested
            return
        endif

        matrix = matmul(matrix_requested, matrix)
        !matrix = matmul(matrix, matrix_requested)
    endsubroutine build_rotation_matrix

endmodule rotations_mod
