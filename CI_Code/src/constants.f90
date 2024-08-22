module constants

implicit none

	real(8), parameter			::	alpha = 0.0072973525693d0
	real(8), parameter			::	diss_en = 14.70976
	real(8), parameter			::	pi = 3.14159265358979d0
	real(8), parameter			::	ev_au = .036749405469679d0
	real(8), parameter			::	omega_0 = 10.64d0*ev_au
	real(8), parameter			::	const = (8d0/3d0)*pi*alpha
	real(8), parameter			::	extra_energy = .1*ev_au
	real(8), parameter			::	r_max = 3.0d0
	real(8), parameter			::	r_min = 1.8d0
	real(8), parameter			::	gamma = .1*ev_au
	real(8), parameter          ::  sigma = .1*ev_au
    real(8), dimension(0:2), parameter      ::  delta_Es = ( (/0.000, .33858, .66586/) ) * ev_au
    real(8), dimension(0:2), parameter		::	weights = (/.50, .25, .25/)
    real(8), parameter          ::  unit_tol = .5
    real(8), parameter          ::  elec_en_range = .05*ev_au

	character (len = 2), dimension(3)	::	sym_strings = (/"A1", "A2", "B1"/)
	character (len = 1), dimension(3)	::	comp_strings = (/"x", "y", "z"/)

	integer, parameter			::	num_targ = 3
	integer, parameter			::	energy_step_num = 500000
	integer, parameter			::	interval_num = 12
	integer, parameter			::	vib_num = 5
	integer, parameter			::	l_max = 4
	integer					    ::	num_channels, preft_num_channels
	integer, parameter			::	num_conv_points = 6000

	real(8), parameter			::  interp_min = 1.5, interp_max = 3.3
	integer						::	num_interp_points = 19
	real(8), allocatable  		::	interp_geoms(:)


	complex(8), parameter			::	ci = (0d0,1d0)

	real(8)			::	int_r_interval
	real(8), parameter ::	r_interval = (r_max - r_min)/interval_num

contains

subroutine get_interp_grid

	integer :: i, count
	real(8) :: r

	allocate(interp_geoms(num_interp_points))

	int_r_interval = (interp_max - interp_min) / (num_interp_points-1)

	do i = 1, num_interp_points
		interp_geoms(i) = interp_min + int_r_interval * (i-1)
	end do

end subroutine

subroutine determine_num_channels(l_max, num_targ, vib_num, num_channels)

integer, intent(in)	::	l_max, num_targ, vib_num
integer, intent(out)	::	num_channels

preft_num_channels = num_targ*(l_max+1)**2
num_channels = preft_num_channels*vib_num

end subroutine

subroutine get_geom_list(r_interval, r_min, interval_num, geom_list)

integer, intent(in)				::	interval_num
real(8), intent(in)				::	r_min, r_interval

integer						::	r_i
real(8)						::	r
character(len = 3)				::	r_str

character(len = 3), dimension(0:interval_num), intent(out)	::	geom_list

do r_i = 0, interval_num
	r = r_min + (r_i * r_interval)
	write(r_str,"(F3.1)") r
	geom_list(r_i) = r_str
end do

end subroutine
end module constants
