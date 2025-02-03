subroutine      initiate_energetic_parameters
    use ackland4_approximation_arrays_mod
    use potentials_ackland4fefe_linear_mod
    integer(integer_byte_size) i
    real(float_byte_size) r, rho
    real(float_byte_size) :: pw_ack4pot, ed_ack4pot, mf_ack4pot
    real(float_byte_size) :: d1pw_ack4pot, d1ed_ack4pot, d1mf_ack4pot
    real(float_byte_size) :: d2pw_ack4pot, d2ed_ack4pot, d2mf_ack4pot
    open (90, file = 'pw.txt');    open (91, file = 'ed.txt');    open (92, file = 'mf.txt')
    do i = 1, pot_steps
        r   = i * pairwise_step
        pairwise_pot(i) = pw_ack4pot(r)
        pairwise_dr1(i) = d1pw_ack4pot(r)
        pairwise_dr2(i) = d2pw_ack4pot(r)
        if(i.ge.2) then! we need average values of derivatives to avoid gap due to higher parts
            pairwise_ad1(i-1)=(pairwise_pot(i)-pairwise_pot(i-1))*pairwise_recistep
            pairwise_ad2(i-1)=(pairwise_dr1(i)-pairwise_dr1(i-1))*pairwise_recistep
        endif
                write(90, *) r, pairwise_pot(i), pw_ack4pot(r)
        r   = i * elecdens_step
        elecdens_pot(i) = ed_ack4pot(r)
        elecdens_dr1(i) = d1ed_ack4pot(r)
        elecdens_dr2(i) = d2ed_ack4pot(r)
        if(i.ge.2) then! we need average values of derivatives to avoid gap due to higher parts
            elecdens_ad1(i-1)=(elecdens_pot(i)-elecdens_pot(i-1))*elecdens_recistep
            elecdens_ad2(i-1)=(elecdens_dr1(i)-elecdens_dr1(i-1))*elecdens_recistep
        endif
                write(91, *) r, elecdens_pot(i), ed_ack4pot(r)

        rho = i * embefunc_step !according to formula it must be rho
        embefunc_pot(i) = mf_ack4pot(rho)
        embefunc_dr1(i) = d1mf_ack4pot(rho)
        embefunc_dr2(i) = d2mf_ack4pot(rho)
        if(i.ge.2) then! we need average values of derivatives to avoid gap due to higher parts
            embefunc_ad1(i-1)=(embefunc_pot(i)-embefunc_pot(i-1))*embefunc_recistep
            embefunc_ad2(i-1)=(embefunc_dr1(i)-embefunc_dr1(i-1))*embefunc_recistep
        endif
                write(92, *) r, embefunc_pot(i), mf_ack4pot(rho)
    enddo
    pairwise_ad1(pot_steps) = pairwise_ad1(pot_steps-1)
    pairwise_ad2(pot_steps) = pairwise_ad2(pot_steps-1)
    elecdens_ad1(pot_steps) = elecdens_ad1(pot_steps-1)
    elecdens_ad2(pot_steps) = elecdens_ad2(pot_steps-1)
    embefunc_ad1(pot_steps) = embefunc_ad1(pot_steps-1)
    embefunc_ad2(pot_steps) = embefunc_ad2(pot_steps-1)
    close(90);close(91);close(92)
    print*, "Energetic parameters (lin) are obtained from analytical form. OK"
endsubroutine initiate_energetic_parameters
