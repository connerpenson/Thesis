module read_lib

use constants

implicit none

contains

subroutine read_smat(file_path, num_channels, smat)

character (len = 256), intent(in)				::	file_path
integer, intent(in)						::	num_channels

integer								::	i
real(8), dimension(2, num_channels)				::	tmp
character (len = 256)						::	cwd

complex(8), dimension(num_channels, num_channels), intent(out)	::	smat

call getcwd(cwd)

open(1, file = trim(file_path))
do i = 1, num_channels
        read(1,*) tmp(1:2, 1:i)
        smat(1:i, i) = tmp(1, 1:i) + tmp(2, 1:i)*ci
        smat(i, 1:i-1) = smat(1:i-1,i)
end do
close(1)

open(1, file = trim(cwd) // "/smat.dat")
do i = 1, num_channels

        write(1,'(E12.4)') sum(abs(smat(i,1:num_channels))**2)

end do
close(1)

end subroutine

subroutine read_neut_WF(r_interval, interval_num, neut_WF_array)

real(8), intent(in)				::	r_interval
integer, intent(in)				::	interval_num

character (len = 64)				::	file_path

integer						::	j
character (len = 3)				::	r_str
character (len = 256)				::	cwd
real(8)						::	r

real(8), intent(out), dimension(interval_num)	::	neut_WF_array

call getcwd(cwd)

do j = 1, interval_num

        r = 1.8 + r_interval*(j-1)

        write(r_str, "(F3.1)") r
	file_path = trim(trim(cwd) // "/wf000_neut.dat")
        open(2, file = file_path, status = "old")
        call WF_mag(2, r, neut_WF_array(j))

end do

end subroutine

subroutine read_ion_WF(r_interval, interval_num, vib_num, ion_WF_array)

real(8), intent(in)				::	r_interval
integer, intent(in)				::	interval_num, vib_num

character (len = 64)				::	file_path

character (len = 256)				::	cwd
integer						::	j, chan_i
character (len = 3)				::	r_str, chan_str
real(8)						::	r

real(8), intent(out), dimension(interval_num, vib_num)	::	ion_WF_array

call getcwd(cwd)

do chan_i = 0, 4
        do j = 1, interval_num

                r = 1.8 + r_interval*(j-1)

                write(chan_str, "(i1)") chan_i

                file_path = trim(trim(cwd) // "/wf00" // chan_str // "_ion.dat")
                open(3, file = file_path, status = "old")
                call WF_mag(3, r, ion_WF_array((j),(chan_i+1)))
        end do
end do

end subroutine

subroutine WF_mag(file_num, r, WF)

    !Routine pulls wf value from specified file for the given r value
    real(8)                  :: diff, var, r, res
    integer                  :: counter, file_num
    real(8), intent(out)     :: WF
    real(8), dimension(4000) :: WF_data

    read(file_num, *) WF_data
    diff = 10
    counter = 0

    res = WF_data(5)-WF_data(1)

    do while (abs(diff) > res)
        var = WF_data(4*counter + 1)
        diff = var - r
        WF = WF_data(4*counter + 2)
        counter = counter + 1
    end do
    close(file_num)

end subroutine

subroutine read_dip_data(interval_num, r_min, r_interval, energy_i, dipole_data)

integer, intent(in)						::	interval_num, energy_i
real(8), intent(in)						::	r_min, r_interval

complex(8), dimension(3,3,9,interval_num+1), intent(out)	::	dipole_data

integer								::	sym_i, ang_i, ang_chan, r_i
real(8)								::	r
character (len = 3)						::	r_str, sym_str

do sym_i = 1,3

        if (sym_i == 1) then
        ang_chan = 9
        else if (sym_i == 2 .OR. sym_i == 3) then
        ang_chan = 6
        end if

        do ang_i = 1, ang_chan
                do r_i = 1, interval_num + 1

                        r = r_min + r_interval*(r_i - 1)
                        write(r_str, "(F3.1)") r
                        sym_str = trim(sym_strings(sym_i))

                        call dip_mag(energy_i, ang_i, r_str, sym_str, dipole_data(:, sym_i, ang_i, r_i))

                end do
        end do
end do

end subroutine

subroutine dip_mag(energy_i, ang_i, r_str, sym_str, dip_mom)

    real(8)							::	x_real, x_imag, z_real, z_imag, y_real, y_imag
    complex(8)							::	x, y, z
    integer, intent(in)						::	energy_i, ang_i
    real(8), dimension(2  + 8*(energy_i-1) + 4000*(ang_i-1))	::	temp1
    character (len = 256)					::	cwd

    character(len = 3), intent(in)				::	r_str 
    character(len = 2), intent(in)				::	sym_str

    complex(8), dimension(3), intent(out)			::	dip_mom

    call getcwd(cwd)

    open(1, file = trim(cwd) // "/DipoleData/" // sym_str // "_" // r_str, status = "old")
    read (1, *) temp1, x_real, x_imag, z_real, z_imag, y_real, y_imag
    close(1)


    x = cmplx(x_real, x_imag, 8)
    y = cmplx(y_real, y_imag, 8)
    z = cmplx(z_real, z_imag, 8)

    dip_mom(1) = x
    dip_mom(2) = y
    dip_mom(3) = z

    return

end subroutine dip_mag

end module read_lib
