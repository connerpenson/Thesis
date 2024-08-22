program coldfront

    implicit none

    real, parameter :: temp_A = 12., temp_M = 24., dist = 960., speed = 20.
    real, parameter :: elapsed_time = 24.

    real :: new_temp_M, temp_diff, grad, rate

    temp_diff = temp_A - temp_M
    grad = temp_diff/dist
    rate = grad * speed

    new_temp_M = temp_M + rate * elapsed_time

    print *, "Temperature after", elapsed_time, "hours is", new_temp_M, "degrees"

end program coldfront
