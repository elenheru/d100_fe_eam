subroutine      set_seed_of_random_numbers(truly_random)
    logical, intent(in) :: truly_random
    integer(4) date_values(8), feed_input(100)
    integer(4), allocatable :: seed(:)
    integer(4) :: n
    if (truly_random) then
        call date_and_time(VALUES=date_values)
        feed_input = (                                  &
                    (date_values(5)*60+date_values(6)  &
                    ) *60 + date_values(7)             &
                    ) *1000 + date_values(8)
        call random_seed(put=feed_input)
        print *, "Random number generator got date_and_time as seed. OK"
    else
        call random_seed(size = n)
        allocate(seed(n))
        seed = 19940305
        print *, "Random number generator got 19940305 as seed. OK"
        call random_seed(get=seed)
    endif
endsubroutine   set_seed_of_random_numbers
