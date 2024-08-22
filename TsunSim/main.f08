program TsunamiSimulator

    use iso_fortran_env, only: int32, real32
    use mod_diff, only: diff
    use mod_initial, only: set_gaussian

    implicit none

    integer(int32) :: i, n
    integer(int32), parameter :: grid_size = 100
    integer(int32), parameter :: num_time_steps = 100
    integer(int32), parameter :: icenter = 25


    real(real32), parameter :: dt = 1, dx = 1, c = 1
    real(real32), parameter :: decay = 0.02

    real(real32) :: h(grid_size), dh(grid_size)

    !Input checks
    if (grid_size <= 0) stop "Grid size must be greater than zero!"
    if (dt <= 0) stop "Time step must be greater than zero!"
    if (dx <= 0) stop "X-interval size must be greater than zero!"
    if (c <= 0) stop "Flow speed must be greater than zero!"

    !Setting initial wave heights
    call set_gaussian(h, icenter, decay)

    print *, 0, h

    time_loop: do n = 1, num_time_steps
        h = h - c * diff(h) / dx * dt
        print *, n, h
    end do time_loop

end program TsunamiSimulator
