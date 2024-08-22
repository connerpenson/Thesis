
program dip_moment
    implicit none
    complex, dimension(3,3) :: mom_array
    real :: r, WF_n, Wf_i, energy_val
    integer :: energy, i, j, chan
    character (len = 57) :: WF_file
    character (len = 3) :: r_str
    character (len = 1) :: chan_str
    character (len = 256) :: cwd
    real, dimension(2000) :: temp
    complex, dimension(3,3) :: int_array, sum_array, temp_array
    complex, dimension(500,3):: A1_dat, B1_dat, B2_dat
    real, dimension(500):: energy_array

 do chan = 0, 4, 1
        print *, chan

  do energy = 0, 499, 1

    !set Riemann sum matrix to 0 before all loops
    sum_array = 0

        !integration variable loop
        !done in this weird way to avoid non-integer loop variable
        do j = 0, 3
        r = 1.8 + .4*j

            !getting cd for relative file paths
            call getcwd(cwd)

            !getting string for file paths
            write(r_str, "(F3.1)") r

            !subroutine returns dipole moment for specified symmetry at a specified energy
            call dip_mag(trim(cwd) // "/DipoleData/A1_" // r_str, energy, mom_array(1:3,1))
            call dip_mag(trim(cwd) // "/DipoleData/B1_" // r_str, energy, mom_array(1:3,2))
            call dip_mag(trim(cwd) // "/DipoleData/B2_" // r_str, energy, mom_array(1:3,3))

            write(chan_str, "(i1)") chan

                !These function calls get the value of the wavefunction at the specified r
                WF_file = trim(trim(cwd) // "/wf000_neut.dat")
                call WF_mag(WF_file, r, WF_n)

                WF_file = trim(trim(cwd) // "/wf00" // chan_str // "_ion.dat")
                call WF_mag(WF_file, r, WF_i)

                !regarding that the wavefunctions for both neutral and ion are fully real
                !integral must be done for xyz for all symmetries

                !A1
                !x
                int_array(1,1) = mom_array(1,1) * WF_n * WF_i
                !y
                int_array(1,2) = mom_array(1,2) * WF_n * WF_i
                !z
                int_array(1,3) = mom_array(1,3) * WF_n * WF_i

                !B1
                !x
                int_array(2,1) = mom_array(2,1) * WF_n * WF_i
                !y
                int_array(2,2) = mom_array(2,2) * WF_n * WF_i
                !z
                int_array(2,3) = mom_array(2,3) * WF_n * WF_i

                !B2
                !x
                int_array(3,1) = mom_array(3,1) * WF_n * WF_i
                !y
                int_array(3,2) = mom_array(3,2) * WF_n * WF_i
                !z
                int_array(3,3) = mom_array(3,3) * WF_n * WF_i

                !Multiply by dr
                temp_array = int_array * .4
                !Add across r
                sum_array = sum_array + temp_array
       
        !r do
        end do

    energy_val = .735E-3 + .735E-3*energy

    !Now printing out the x,y,z for each symmetry
    !x
    A1_dat(energy+1,1) = sum_array(1,1)
    !y
    A1_dat(energy+1,2) = sum_array(1,2)
    !z
    A1_dat(energy+1,3) = sum_array(1,3)

    !x
    B1_dat(energy+1,1) = sum_array(2,1)
    !y
    B1_dat(energy+1,2) = sum_array(2,2)
    !z
    B1_dat(energy+1,3) = sum_array(2,3)

    !x
    B2_dat(energy+1,1) = sum_array(3,1)
    !y
    B2_dat(energy+1,2) = sum_array(3,2)
    !z
    B2_dat(energy+1,3) = sum_array(3,3)

    energy_array(energy+1) = energy_val

  !energy do
  end do

!Writing all calculated data into .dat files
!format: energy real imaginary
print *, chan_str

open(2, file = trim(trim(cwd) // "/N0_I" // chan_str // "/A1x_dat.dat") )
do i=1,500
  write(2,'(3e20.12)') energy_array(i), real(A1_dat(i,1)), aimag(A1_dat(i,1))
end do
close(2)

print *, "printed1"

open(2, file = trim(trim(cwd) // "/N0_I" // chan_str // "/A1y_dat.dat") )
do i=1,500
  write(2,'(3e20.12)') energy_array(i), real(A1_dat(i,2)), aimag(A1_dat(i,2))
end do
close(2)

print *, "printed2"

open(2, file = trim(trim(cwd) // "/N0_I" // chan_str // "/A1z_dat.dat") )
do i=1,500
  write(2,'(3e20.12)') energy_array(i), real(A1_dat(i,3)), aimag(A1_dat(i,3))
end do
close(2)

print *, "printed3"

open(2, file = trim(trim(cwd) // "/N0_I" // chan_str // "/B1x_dat.dat") )
do i=1,500
  write(2,'(3e20.12)') energy_array(i), real(B1_dat(i,1)), aimag(B1_dat(i,1))
end do
close(2)

print *, "printed4"

open(2, file = trim(trim(cwd) // "/N0_I" // chan_str // "/B1y_dat.dat") )
do i=1,500
  write(2,'(3e20.12)') energy_array(i), real(B1_dat(i,2)), aimag(B1_dat(i,2))
end do
close(2)

print *, "printed5"

open(2, file = trim(trim(cwd) // "/N0_I" // chan_str // "/B1z_dat.dat") )
do i=1,500
  write(2,'(3e20.12)') energy_array(i), real(B1_dat(i,3)), aimag(B1_dat(i,3))
end do
close(2)

print *, "printed6"

open(2, file = trim(trim(cwd) // "/N0_I" // chan_str // "/B2x_dat.dat") )
do i=1,500
  write(2,'(3e20.12)') energy_array(i), real(B2_dat(i,1)), aimag(B2_dat(i,1))
end do
close(2)

print *, "printed7"

open(2, file = trim(trim(cwd) // "/N0_I" // chan_str // "/B2y_dat.dat") )
do i=1,500
  write(2,'(3e20.12)') energy_array(i), real(B2_dat(i,2)), aimag(B2_dat(i,2))
end do
close(2)

print *, "printed8"

open(2, file = trim(trim(cwd) // "/N0_I" // chan_str // "/B2z_dat.dat") )
do i=1,500
  write(2,'(3e20.12)') energy_array(i), real(B2_dat(i,3)), aimag(B2_dat(i,3))
end do
close(2)

print *, "printed9"

!channel do
end do 

contains

subroutine dip_mag(file_path, energy, dip_mom)

    character (len = *), intent(in) :: file_path

    real :: x_real, x_imag, z_real, z_imag, y_real, y_imag
    complex :: x, y, z
    integer, intent(in) :: energy
    real, dimension(2 + 8*energy) :: temp

    !dip_mom now an array with dip_mom(1,2,3) being x y and z comp
    complex, dimension(3), intent(out) :: dip_mom

    open (1, file = file_path, status= 'old')
    read (1, *) temp, x_real, x_imag, z_real, z_imag, y_real, y_imag
    close(1)

    x = cmplx(x_real, x_imag)
    y = cmplx(y_real, y_imag)
    z = cmplx(z_real, z_imag)

    dip_mom(1) = x
    dip_mom(2) = y
    dip_mom(3) = z

    return

end subroutine dip_mag

subroutine WF_mag(file_path, r, WF)

    !Routine pulls wf value from specified file for the given r value
    character (len = *), intent(in) :: file_path
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
