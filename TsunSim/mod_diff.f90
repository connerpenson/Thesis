module mod_diff

    use iso_fortran_env, only: int32, real32

    implicit none

    contains

    function diff(x) result(dx)
        real(real32), intent(in)    ::  x(:)
        real(real32)                ::  dx(size(x))
        integer(int32)             ::  n

        n = size(x)
        dx(1) = x(1) - x(n)
        dx(2:n) = x(2:n) - x(1:n-1)
    end function diff

end module mod_diff
