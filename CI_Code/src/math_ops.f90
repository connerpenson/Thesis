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


subroutine cauchy_conv(E_min, E_max, size_old, TCS_array, conv_TCS_array)

integer, dimension(2), intent(in)				::	size_old
real(8), dimension(:, :), intent(in)	::	TCS_array
real(8), intent(in)                     ::  E_min, E_max

real(8), dimension(num_conv_points, 2), intent(out)	::	conv_TCS_array

real(8)			::	E_in, E_fin, conv_E_step, E_curr, TCS_E_step, int_sum
integer			::	i, j

real(8), allocatable			::	cauchy_list(:)

print *, size_old(1), size_old(2)

E_in = minval(TCS_array(:,2))
E_fin = maxval(TCS_array(:,2))
TCS_E_step = TCS_array(2,2) - TCS_array(1,2)
conv_E_step = (E_max-E_min)/num_conv_points

allocate(cauchy_list(size_old(1)))

print *, "Using ", num_conv_points, " points for convolution grid"

do i = 1, num_conv_points
int_sum = 0
E_curr = E_min + (conv_E_step*(i-1))
    if (E_curr >= E_in .AND. E_curr <= E_fin) then
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
deallocate(cauchy_list)
end subroutine

subroutine gauss_conv(E_min, E_max, size_old, TCS_array, conv_TCS_array)

integer, dimension(2), intent(in)				::	size_old
real(8), dimension(:, :), intent(in)	::	TCS_array
real(8), intent(in)                     ::  E_min, E_max

real(8), dimension(num_conv_points, 2), intent(out)	::	conv_TCS_array

real(8)			::	E_in, E_fin, conv_E_step, E_curr, TCS_E_step, int_sum
integer			::	i, j

real(8), allocatable			:: cauchy_list(:)

E_in = minval(TCS_array(:,2))
E_fin = maxval(TCS_array(:,2))
TCS_E_step = TCS_array(2,2) - TCS_array(1,2)
conv_E_step = (E_max-E_min)/num_conv_points

open(2, file = "SPES_missing_peak.dat")
do i = 1, size_old(1)

    write(2, *) TCS_array(i,1), TCS_array(i,2)
 
end do
close(2)

allocate(cauchy_list(size_old(1)))

print *, "Using ", num_conv_points, " points for convolution grid"


do i = 1, num_conv_points
int_sum = 0
E_curr = E_min + (conv_E_step*(i-1))
    if (E_curr >= E_in .AND. E_curr <= E_fin) then
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

deallocate(cauchy_list)
end subroutine

! subroutine prune_dips(dipole_data, reduced_dipole_data)

!     complex(8), intent(in)  ::  dipole_data(:,:,:)

!     integer                 ::  chan_i, coord_i

!     complex(8), intent(out) ::  reduced_dipole_data(:,:,:)

!     do chan_i = 1, num_channels
!         do coord_i = 1, 3

            


!         end do
!     end do


end module math_ops
