subroutine relax_fine_quad_grad_rand(na, h)
    use positions_mod
    implicit none
    integer(integer_byte_size), intent(inout) :: na
    real(float_byte_size), intent(inout) :: h
    real(float_byte_size), parameter :: leapover_shrink = 4d0
    integer(integer_byte_size) i
    real(float_byte_size) :: e000, ep00, em00, e0p0, e0m0, e00p, e00m, esss, enorm, eg, eran
    real(float_byte_size) :: xadd, yadd, zadd, addnorm, Rs(3), Rmin(3), R_ran(3)

    call random_number(R_ran)
    R_ran = R_ran * 1d-12 * (4d-1**(R_ran * 25d0))
    R_at(1:3, na) = R_at(1:3, na) + R_ran(1:3)

    Rs(1:3) = R_at(1:3, na)
    Rmin(1:3) = R_at(1:3, na)

    call calculate_pseudosystem_energy(e000, na)
    esss = e000

    ! x addition
    r_at(1, na) = r_at(1, na) + h
    call calculate_pseudosystem_energy(ep00, na)
    if (esss .gt. ep00) then
        esss = ep00
        Rmin(1:3) = R_at(1:3, na)
    endif

    r_at(1, na) = r_at(1, na) - h * 2d0
    call calculate_pseudosystem_energy(em00, na)
    if (esss .gt. em00) then
        esss = em00
        Rmin(1:3) = R_at(1:3, na)
    endif

    r_at(1, na) = r_at(1, na) + h

    ! y addition
    r_at(2, na) = r_at(2, na) + h
    call calculate_pseudosystem_energy(e0p0, na)
    if (esss .gt. e0p0) then
        esss = e0p0
        Rmin(1:3) = R_at(1:3, na)
    endif

    r_at(2, na) = r_at(2, na) - h * 2d0
    call calculate_pseudosystem_energy(e0m0, na)
    if (esss .gt. e0m0) then
        esss = e0m0
        Rmin(1:3) = R_at(1:3, na)
    endif

    r_at(2, na) = r_at(2, na) + h

    ! z addition
    r_at(3, na) = r_at(3, na) + h
    call calculate_pseudosystem_energy(e00p, na)
    if (esss .gt. e00p) then
        esss = e00p
        Rmin(1:3) = R_at(1:3, na)
    endif

    r_at(3, na) = r_at(3, na) - h * 2d0
    call calculate_pseudosystem_energy(e00m, na)
    if (esss .gt. e00m) then
        esss = e00m
        Rmin(1:3) = R_at(1:3, na)
    endif

    r_at(3, na) = r_at(3, na) + h


    xadd = ( h * 5d-1 * (ep00 - em00) / (2d0 * e000 - ep00 - em00) )
    yadd = ( h * 5d-1 * (e0p0 - e0m0) / (2d0 * e000 - e0p0 - e0m0) )
    zadd = ( h * 5d-1 * (e00p - e00m) / (2d0 * e000 - e00p - e00m) )

    if ( isnan(xadd) .or. isnan(xadd - xadd) ) then
!        print*, "X: can not relax quartic ", h, e000, ep00, em00
        R_at(1:3, na) = Rmin(1:3)
        return
    endif
    if ( isnan(yadd) .or. isnan(yadd - yadd) ) then
!        print*, "Y: can not relax quartic ", h, e000, e0p0, e0m0
        R_at(1:3, na) = Rmin(1:3)
        return
    endif
    if ( isnan(zadd) .or. isnan(zadd - zadd) ) then
!        print*, "Z: can not relax quartic ", h, e000, e00p, e00m
        R_at(1:3, na) = Rmin(1:3)
        return
    endif

    r_at(1, na) = r_at(1, na) + xadd
    r_at(2, na) = r_at(2, na) + yadd
    r_at(3, na) = r_at(3, na) + zadd

    call calculate_pseudosystem_energy(e000, na)

    if(e000 .gt. esss) then
!        print*, "minimization failed, minimum is not in quadratic extremum, but in octahedron vertice or center"
!        print*, "we'll try gradient"
        R_at(1:3, na) = Rmin(1:3)
        !RETURN !2.56e-10 !2.38e-6 !2.44e-6
    else
        return
    endif

    !trying gradient simple step
    enorm = sqrt( (ep00-em00)**2 + (e0p0 - e0m0)**2 + (e00p - e00m)**2 )
    if( isnan((1d0/enorm ) - (2d0/enorm )) ) then
        call random_number(R_ran)
        R_ran = R_ran * 2d-14
        R_at(1:3, na) = Rmin(1:3)!+R_ran(1:3)
        return
    endif
    xadd = -5d-1 * (ep00-em00) / enorm / leapover_shrink
    yadd = -5d-1 * (e0p0-e0m0) / enorm / leapover_shrink
    zadd = -5d-1 * (e00p-e00m) / enorm / leapover_shrink !arbitrary shrink for step, to avoid overleap

    r_at(1, na) = rs(1) + xadd
    r_at(2, na) = rs(2) + yadd
    r_at(3, na) = rs(3) + zadd

    call calculate_pseudosystem_energy(ep00, na)

    if(ep00 .gt. esss) then
!        print*, "gradien minimization is not successful also. we'll try quadratic via gradient"
        R_at(1:3, na) = Rmin(1:3)
        RETURN !1.62e-10! 9.53e-7 !9.32e-7
    else
        return
    endif

    !trying quadratic descent via gradient
    r_at(1, na) = rs(1) - xadd
    r_at(2, na) = rs(2) - yadd
    r_at(3, na) = rs(3) - zadd

    call calculate_pseudosystem_energy(em00, na)
    addnorm = sqrt(xadd*xadd + yadd*yadd + zadd*zadd)
    if( ( isnan( 1d0/(2d0 * e000 - ep00 - em00) - 2d0/(2d0 * e000 - ep00 - em00) ) ) &
        .or. ( isnan((1d0/addnorm ) - (2d0/addnorm )) ) ) then
        R_at(1:3, na) = Rmin(1:3)
        return
    endif

    xadd = ( xadd * 5d-1 * (ep00 - em00) / (2d0 * e000 - ep00 - em00) ) / addnorm
    yadd = ( yadd * 5d-1 * (ep00 - em00) / (2d0 * e000 - ep00 - em00) ) / addnorm
    zadd = ( zadd * 5d-1 * (ep00 - em00) / (2d0 * e000 - ep00 - em00) ) / addnorm
    r_at(1, na) = rs(1) + xadd
    r_at(2, na) = rs(2) + yadd
    r_at(3, na) = rs(3) + zadd
    call calculate_pseudosystem_energy(eg, na)

    if(eg .gt. esss) then
        !print*, "quadratic via gradient minimization is not successful also. we'll try several random steps"
        R_at(1:3, na) = Rmin(1:3)
        !RETURN !1.66e-10 !9.61e-7 ! 9.20e-7
    else
        return
    endif

    !trying random steps
    do i=1, 10
        call random_number(R_ran)
        R_ran = R_ran * 2d-1**(-i-7)
        r_at(1:3, na) = r_at(1:3, na) + R_ran(1:3)
        call calculate_pseudosystem_energy(eran, na)

        if(eran .gt. esss) then
            R_at(1:3, na) = Rmin(1:3)
        else
            return !1.364e-10 !9.38e-7! 9.33e-7
        endif
    enddo
    !print*, "random step minimization is never successful also"

endsubroutine relax_fine_quad_grad_rand
