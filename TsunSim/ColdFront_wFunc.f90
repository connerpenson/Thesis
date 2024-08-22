program coldfront_wfunc

    implicit none

    real :: nhours(8)

    nhours = [6., 12., 18., 24., 30., 36., 42., 48.]
    print *, cold_front_temp(12., 24., 20., 960., nhours)

    contains

    pure elemental real function cold_front_temp(temp_A, temp_M, c, dx, dt) result(res)
        real, intent(in) :: temp_A, temp_M, c, dx, dt
        res = temp_M - (temp_M - temp_A)/dx * c *dt
    end function cold_front_temp

end program coldfront_wfunc
