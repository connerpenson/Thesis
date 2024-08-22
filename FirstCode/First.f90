
program dip_moment
    implicit none
    real, dimension(3) :: mom_array
    real :: r, WF_gr, Wf_ex, energy_val
    integer :: energy, i, j
    character (len = 51) :: WF_file
    character (len = 3) :: r_str
    character (len = 256) :: cwd
    real, dimension(2000) :: temp
    real, dimension(3) :: int_array, sum_array, temp_array
    real, dimension(500):: A1_dat, B1_dat, B2_dat
    real, dimension(500):: energy_array

call getcwd(cwd)

do energy = 0,499
    

    sum_array = 0


        do j = 0, 3
        r = 1.8 + .4*j

            call getcwd(cwd)


            write(r_str, "(F3.1)") r

            call dip_mag(trim(cwd) // "/bound_states/Dipole Data - CH_RVaried/A1_" // r_str, energy, mom_array(1))
            call dip_mag(trim(cwd) // "/bound_states/Dipole Data - CH_RVaried/B1_" // r_str, energy, mom_array(2))
            call dip_mag(trim(cwd) // "/bound_states/Dipole Data - CH_RVaried/B2_" // r_str, energy, mom_array(3))

            print "(3e20.12)", mom_array


                WF_file = trim(trim(cwd) // "/bound_states/wf000.dat")
                call WF_mag(WF_file, r, WF_gr)

                WF_file = trim(trim(cwd) // "/bound_states/wf001.dat")
                call WF_mag(WF_file, r, WF_ex)

                print *, WF_gr, Wf_ex

                int_array(1) = mom_array(1) * WF_gr * WF_ex
                int_array(2) = mom_array(2) * WF_gr * WF_ex
                int_array(3) = mom_array(3) * WF_gr * WF_ex

                temp_array = int_array * .4
                sum_array = sum_array + temp_array

        end do


    energy_val = .735E-3 + .735E-3*energy

    A1_dat(energy+1) = sum_array(1)
    B1_dat(energy+1) = sum_array(2)
    B2_dat(energy+1) = sum_array(3)

    energy_array(energy+1) = energy_val


end do

open(2, file = trim(trim(cwd) // "/bound_states/A1_dat.dat") )
do i=1,500
  write(2,'(2e20.12)') energy_array(i), A1_dat(i)
end do
close(2)

open(3, file = trim(trim(cwd) // "/bound_states/B1_dat.dat") )
do i=1,500
  write(3,'(4e20.12)') energy_array(i), B1_dat(i)
end do
close(3)

open(4, file = trim(trim(cwd) // "/bound_states/B2_dat.dat") )
do i=1,500
  write(4,'(4e20.12)') energy_array(i), B2_dat(i)
end do
close(4)

print "(3e20.12)", A1_dat
print "(3e20.12)", B1_dat
print "(3e20.12)", B2_dat
contains

subroutine dip_mag(file_path, energy, dip_mom)

    character (len = 100), intent(in) :: file_path

    real :: x_real, x_imag, z_real, z_imag, y_real, y_imag, norm_x, norm_y, norm_z
    integer, intent(in) :: energy
    real, dimension(2 + 8*energy) :: temp
    real, intent(out) :: dip_mom

    open (1, file = file_path, status= 'old')
    read (1, *) temp, x_real, x_imag, z_real, z_imag, y_real, y_imag
    close(1)

    norm_x = sqrt(x_real**2 + x_imag**2)
    norm_y = sqrt(y_real**2 + y_imag**2)
    norm_z = sqrt(z_real**2 + z_imag**2)

    dip_mom = sqrt(norm_x**2 + norm_y**2 + norm_z**2)

    return

end subroutine dip_mag

subroutine WF_mag(file_path, r, WF)

    character (len = 51), intent(in) :: file_path
    real :: diff, var, r
    integer :: counter
    real, intent(out) :: WF

    open(1, file = file_path, status = "old")
    read(1, *) temp
    do counter = 1, 499
        var = temp(4*counter +1)
        diff = var - r
        if (abs(diff) < .006) then
        WF = temp(4*counter + 2)
        end if
    end do
    close(1)

end subroutine

end program
