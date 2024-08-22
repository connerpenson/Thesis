module math_ops

use constants

implicit none

contains

subroutine find_cderiv(y_fin, y_in, dx, deriv)

	complex(8),	intent(in) 	:: y_fin, y_in
	real(8),	intent(in)	:: dx
	complex(8),	intent(out)	:: deriv

	real(8)				:: real_part, imag_part

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

!hi

subroutine cauchy_conv(size_old, TCS_array, conv_TCS_array)

integer, dimension(2), intent(in)				::	size_old
real(8), dimension(size_old(1), size_old(2)), intent(in)	::	TCS_array

real(8), dimension(num_conv_points, 2), intent(out)	::	conv_TCS_array

real(8)			::	E_in, E_fin, conv_E_step, E_curr, TCS_E_step, int_sum
integer			::	i, j

real(8), allocatable			::	my_list(:,:), cauchy_list(:)

E_in = 10.0d0
E_fin = 12.0d0
TCS_E_step = TCS_array(2,2) - TCS_array(1,2)
conv_E_step = (E_fin-E_in)/num_conv_points

allocate(my_list(size_old(1),2), cauchy_list(size_old(1)))

print *, "Using ", num_conv_points, " points for convolution grid"

do i = 1, num_conv_points
int_sum = 0
E_curr = E_in + (conv_E_step*(i-1))
    if (E_curr >= TCS_array(1,2) .AND. E_curr <= TCS_array(size_old(1),2)) then
        do j = 1, size_old(1)
            cauchy_list(j) = (gamma/pi)*(1/((TCS_array(j,2)-E_curr)**2 + gamma**2))
            int_sum = int_sum + TCS_array(j,1)*cauchy_list(j)*TCS_E_step
        end do
        conv_TCS_array(i,1) = int_sum
    else
        conv_TCS_array(i,1) = 0
    end if
conv_TCS_array(i,2) = E_curr
end do
deallocate(my_list, cauchy_list)
end subroutine

subroutine gauss_conv(size_old, TCS_array, conv_TCS_array)

integer, dimension(2), intent(in)				::	size_old
real(8), dimension(size_old(1), size_old(2)), intent(in)	::	TCS_array

real(8), dimension(num_conv_points, 2), intent(out)	::	conv_TCS_array

real(8)			::	E_in, E_fin, conv_E_step, E_curr, TCS_E_step, int_sum
integer			::	i, j

real(8), allocatable			::	my_list(:,:), cauchy_list(:)

E_in = 10.0d0
E_fin = 12.0d0
TCS_E_step = TCS_array(2,2) - TCS_array(1,2)
conv_E_step = (E_fin-E_in)/num_conv_points

allocate(my_list(size_old(1),2), cauchy_list(size_old(1)))

print *, "Using ", num_conv_points, " points for convolution grid"

do i = 1, num_conv_points
int_sum = 0
E_curr = E_in + (conv_E_step*(i-1))
    if (E_curr >= TCS_array(1,2) .AND. E_curr <= TCS_array(size_old(1),2)) then
        do j = 1, size_old(1)
            cauchy_list(j) = (1/(sigma*sqrt(2d0*pi)))*exp(-.5d0*((TCS_array(j,2)-E_curr)/sigma)**2)
            int_sum = int_sum + TCS_array(j,1)*cauchy_list(j)*TCS_E_step
        end do
        conv_TCS_array(i,1) = int_sum
    else
        conv_TCS_array(i,1) = 0
    end if
conv_TCS_array(i,2) = E_curr
end do
deallocate(my_list, cauchy_list)
end subroutine

function linearize(dips) result(lin_dips)

    real(8), intent(in)    ::  dips(:)
    real(8)                ::  lin_dips(0:lin_size)
    integer             ::  i, i_Q
    real(8)             ::   a

    lin_dips(0) = dips(1)
    a = (dips(size(dips))-dips(1))/(r_max-r_min)
    do i = 1, lin_size
        lin_dips(i) = lin_dips(i-1) + a * lin_interval
    end do

end function linearize

subroutine frame_transform(data, n, v, n_prime, v_prime, ion_WF_array, ft_data)

    complex(8), intent(in)      ::  data(:)
    real(8), intent(in)         ::  ion_WF_array(:,:,:)
    integer, intent(in)         ::  n, v, n_prime, v_prime

    integer                     ::  r_i
    complex                     ::  int_sum

    complex(8), intent(out)     ::  ft_data

    int_sum = 0
    do r_i = 1, lin_size
        int_sum = int_sum + data(r_i) * ion_WF_array(r_i, v+1, n) * ion_WF_array(r_i, v_prime+1, n_prime) * lin_interval
    end do
    ft_data = int_sum

end subroutine

end module math_ops
