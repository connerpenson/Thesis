module math_ops

implicit none

contains

subroutine find_cderiv(y_fin, y_in, dx, deriv)

complex(8), intent(in)  :: y_fin, y_in
real(8), intent(in)	:: dx
complex(8), intent(out) :: deriv

real(8)			:: real_part, imag_part

real_part = (real(y_fin,8) - real(y_in,8))/dx
imag_part = (dimag(y_fin) - dimag(y_in))/dx

deriv = complex(real_part, imag_part)

end subroutine

subroutine find_rderiv(y_fin, y_in, dx, deriv)

real(8), intent(in)	::	y_fin, y_in
real(8), intent(in)	::	dx
real(8), intent(out)	::	deriv

deriv = (y_fin - y_in)/dx

end subroutine
end module math_ops
