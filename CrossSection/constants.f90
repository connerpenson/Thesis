module constants

implicit none

	real(8), parameter			::	alpha = 0.0072973525693
	real(8), parameter			::	pi = 3.14159265358979
	real(8), parameter			::	ev_au = .036749405469679
	real(8), parameter			::	omega_0 = 10.64*ev_au

	character (len = 2), dimension(3)	::	sym_strings = (/"A1", "B1", "B2"/)
	character (len = 2), dimension(3)	::	comp_strings = (/"x", "y", "z"/)

	complex(8), parameter			::	ci = (0d0,1d0)
	
end module constants
