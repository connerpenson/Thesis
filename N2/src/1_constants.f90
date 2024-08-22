module constants

implicit none

	real(8), parameter			::	alpha = 0.0072973525693d0
	real(8), parameter			::	pi = 3.14159265358979d0
	real(8), parameter			::	ev_au = .036749405469679d0
	real(8), parameter			::	au_ang = .529177d0
	real(8), parameter			::	omega_0 = 10.64d0*ev_au
	real(8), parameter			::	const = (4d0/3d0)*(pi**2)*alpha
	real(8), parameter			::	r_max = 2.3d0
	real(8), parameter			::	r_min = 1.9d0
    real(8), parameter			::	r_max_lin = 2.3d0
	real(8), parameter			::	r_min_lin = 1.9d0 
	real(8), parameter			::	gamma = .35*ev_au
	real(8), parameter          ::  sigma = .25*ev_au
    real(8), dimension(0:2), parameter     ::  delta_Es = (/14.70976, 14.37118, 14.0439/)
    real(8), dimension(0:2), parameter		::	weights = (/.50, .25, .25/)
    real(8), parameter          ::  unit_tol = .5
    real(8), parameter          ::  elec_en_range = .1*ev_au
	real(8), dimension(3), parameter :: dip_en_list = (/4.13805E-2,3.25605E-2,2.3935E-2/) !au
	integer, dimension(3)		::	num_waves = (/54,42,42/)

	character (len = 2), dimension(3)	::	sym_strings = (/"B1", "B2", "B3"/)
	character (len = 1), dimension(3)	::	comp_strings = (/"x", "y", "z"/)

	integer, parameter			::	num_targ = 3
	integer, parameter			::	interval_num = 9
	integer, parameter			::	vib_num = 5
	integer, parameter			::	l_max = 4
	integer					    ::	num_channels = (num_targ*vib_num*(l_max+1)**2)
	integer, parameter			::	num_conv_points = 6000
	integer, parameter          ::  lin_size = 9
	integer, parameter			::	neut_vib_num = 3


	complex(8), parameter			::	ci = (0d0,1d0)

	real(8), parameter			::	r_interval = (r_max - r_min)/interval_num
	real(8), parameter          ::  lin_interval = (r_max_lin - r_min_lin)/lin_size

contains

! subroutine determine_num_channels(l_max, num_targ, vib_num, num_channels)

! integer, intent(in)	::	l_max, num_targ, vib_num
! integer, intent(out)	::	num_channels

! num_channels = (num_targ*vib_num*(l_max+1)**2)
! print *, num_channels

! end subroutine

subroutine get_geom_list(r_interval, r_min, interval_num, geom_list)

integer, intent(in)				::	interval_num
real(8), intent(in)				::	r_min, r_interval

integer						::	r_i
real(8)						::	r
character(len = 4)				::	r_str

character(len = 4), dimension(0:interval_num), intent(out)	::	geom_list

print *, r_interval, r_max, r_min

do r_i = 0, interval_num
	r = r_min + (r_i * r_interval)
	write(r_str,"(F4.2)") r
	geom_list(r_i) = r_str
end do

end subroutine

subroutine get_quant_nums(quant_nums)

    integer ::  n,v,l,m,i,j 
    integer ::  quant_nums(num_channels,4)

            open(2, file = "./Smat/3states_lmax4/lmax4emax3.channel", status = "old")
			do i = 1, num_channels
				read(2, *) j, n, v, l, m
                quant_nums(i,:) = (/n, v, l, m/)
            end do
            close(2)

end subroutine get_quant_nums

subroutine read_neut_WF(neut_WF_array)

	character (len = 65)				::	file_path

	integer						::	j, vib_i
	real(8)						::	r
    character (len = 1)         ::  vib_str

	real(8), intent(out), dimension(lin_size, 0:2)	::	neut_WF_array

do j = 1, lin_size
    r = r_min_lin + lin_interval*(j-1)
    do vib_i = 0, neut_vib_num-1
        write(vib_str, "(i1)") vib_i
        file_path = "./Neutral_WF/wf00" // vib_str // "_neut.dat"
        open(vib_i+2, file = trim(file_path), status = "old")
        call WF_mag(vib_i+2, r, neut_WF_array(j,vib_i))
    end do
end do

end subroutine

subroutine read_ion_WF(ion_WF_array)

	character (len = 286)				::	file_path

	integer						::	state, j, chan_i
	character (len = 1)			::	chan_str, state_str
	real(8)						::	r

	real(8), intent(out), dimension(lin_size, vib_num, num_targ)	::	ion_WF_array


open(50, file = "./wf_data.dat")
do state = 1, num_targ
	do chan_i = 0, vib_num-1
	        do j = 1, lin_size
	                r = (r_min_lin + lin_interval*(j-1))

	                write(chan_str, "(i1)") chan_i
                    write(state_str, "(i1)") state

	                file_path = "./Ion_WF/Target_wf" // chan_str // "_N2+_ElecState00" // state_str // ".dat"
	                open(3, file = file_path, status = "old")
	                call WF_mag(3, r, ion_WF_array(j,(chan_i+1),state))
			write(50, *) r, chan_i, state, ion_WF_array(j,(chan_i+1),state)
	        end do
	end do
end do
close(50)

end subroutine

subroutine WF_mag(file_num, r, WF)

	real(8)                  :: r
	integer                  :: file_num, a, eof
	real(8), intent(out)     :: WF
	real(8), allocatable     :: WF_data(:), r_mat(:), sub_mat(:)
	!real(8)                  :: tmp(100)


    allocate(WF_data(1000), r_mat(1000), sub_mat(1000))
    r_mat = 100
    sub_mat = 100
    WF_data = 100
    do a = 1, 1000
        read(file_num, *, iostat=eof) r_mat(a), WF_data(a)
        if (eof<0) then
            exit
        end if
    end do
    sub_mat = abs(r_mat - r)
    !print *, "r:", r
    !print *, "loc:", minloc(sub_mat, 1)
    !print *, sub_mat(minloc(sub_mat, 1))
    
    WF = WF_data(minloc(sub_mat, 1))
    !print *, "WF:", WF

	close(file_num)
	deallocate(WF_data, r_mat, sub_mat)
end subroutine

end module constants
