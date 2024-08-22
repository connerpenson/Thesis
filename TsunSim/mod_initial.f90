module mod_initial

    use iso_fortran_env, only: int32, real32

    implicit none

    contains

    subroutine set_gaussian(x, icenter, decay)
        real(real32), intent(in out)    ::  x(:)
        integer(int32), intent(in)      ::  icenter
        real(real32), intent(in)        ::  decay

        integer(int32)                  ::  i

        space_loop: do concurrent (i = 1:size(x))
            x(i) = exp(-decay * (i-icenter)**2)
        end do space_loop
    end subroutine set_gaussian
end module mod_initial
