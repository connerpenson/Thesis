 program dip_moment

use constants
use read_lib
use math_ops

    implicit none

    real(8)      ::   row_sum, r, energy_begin, energy_end, energy_step, E, omega, diff, energy, E_min, E_max, temp1, temp2
    integer      ::   state_range, a, i, num_open, vib_counter, e_i
    integer      ::   vib_i, r_i, n, j, l, m, k, INFO, eof, index
    integer      ::   num_closed, coord_i, ch_vib, counter
    complex(8)   ::   x_comp, x_comp_integrated, y_comp, y_comp_integrated, z_comp, z_comp_integrated
    complex(8)	 ::   deriv_1, deriv_2, old_deriv
    integer, dimension(3)		:: rot_ex
    integer, dimension(2)		:: conv_size_array

    character (len = 256)		::   smat_file
    character (len = 3), dimension(0:interval_num)	::	geom_list
	character (len = 1)			::	vib_str
    real(8), dimension(2)		::	state_avg

    complex(8), allocatable		::	dipole_data(:,:,:), smat(:,:), interp_dips(:,:,:), slopes(:,:), beg_slopes(:,:), end_slopes(:,:)
    complex(8), allocatable		::	dipole_phys(:,:), dipole_open(:,:), dipole_closed(:,:), smat_cc(:,:), smat_co(:,:), D(:,:)
    complex(8), allocatable		::	part_vect(:), dip_vect(:,:,:)
    integer,    allocatable		::	IPIV(:), j_list(:), quant_nums(:,:)
    real(8),    allocatable		::	energy_array(:), neut_WF_array(:,:), ion_WF_array(:,:,:), channel_energies(:,:), Z_val(:,:,:)
    real(8),	allocatable		::	TCS_array(:,:,:), conv_TCS_array(:,:,:,:), conv_SPES_array(:,:,:,:), conv_SPES_array_ground(:,:,:,:), &
									& SPES_array(:,:,:), SPES_array_ground(:,:,:), unitarity_check(:), energies(:), temp_energies(:), Z_norms(:)
	integer, allocatable		::	ground_channels(:,:), exc_channels(:,:), max_chans(:,:,:)
	real(8), allocatable  :: nu(:,:,:), beta(:)

    real(8), dimension(100)     ::  temp

call determine_num_channels(l_max, num_targ, vib_num, num_channels)
call get_geom_list(r_interval, r_min, interval_num, geom_list)

allocate(ground_channels(5,3*num_channels), exc_channels(5,3*num_channels))
allocate(quant_nums(num_channels,4), unitarity_check(num_channels), energies(num_targ*vib_num), temp_energies(num_targ*vib_num) )
!remake x

print *, "----------------------------------------------------------------------------------------------------"
print *, "READING S-MATRIX IN..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(smat(num_channels, num_channels))

smat_file = "./Smat/3states_lmax4/vibronic_DR.Smat"

call read_smat(smat_file, num_channels, smat)

print *, "----------------------------------------------------------------------------------------------------"
print *, "EXTRACTING DIPOLE DATA..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(dipole_data(3,num_channels,interval_num+1))
dipole_data = 0
call new_dip_data(dipole_data, quant_nums)

print *, "----------------------------------------------------------------------------------------------------"
print *, "UNITARITY CHECK..."
print *, "----------------------------------------------------------------------------------------------------"

open(1, file = "./unitarity.dat")
write(1,*) "Row", "Sum of Squares"
unitarity_check = 0
do i = 1, num_channels
    do j = 1, num_channels
        unitarity_check(i) = unitarity_check(i) + abs(smat(i,j))**2
    end do
    if (unitarity_check(i) < 10) then
        !smat(i,:) = 0
        !dipole_data(:,i,:) = 0
        write(1,*) i, unitarity_check(i), quant_nums(i,1), quant_nums(i,2), quant_nums(i,3), quant_nums(i,4)
    end if
end do
close(1)
!do i = 1, num_channels
!    diff = abs(1.00d0 - unitarity_check(i))
!    if (diff > unit_tol) then
!        print *, "ERROR!"
!        print *, "The ", i, " row of the S-matrix is not unitary."
!        print *, "Sum of squares:", unitarity_check(i)
!        ERROR STOP
!    else
!        if (i == num_channels) then
!            print *, "Unitary!"
!        end if
!    end if
!end do

print *, "----------------------------------------------------------------------------------------------------"
print *, "DOCTORING DIPOLE DATA..."
print *, "----------------------------------------------------------------------------------------------------"

open(2, file = "./dips_real.dat")
do i = 1, num_channels
	do r_i = 1, interval_num + 1
		r = r_min + (r_interval*(r_i-1))
		if (real(dipole_data(1,i,r_i)) /= 0) then
			write(2, *) r, real(dipole_data(1,i,r_i)), i, quant_nums(i,:)
		else if (real(dipole_data(2,i,r_i)) /= 0) then
			write(2, *) r, real(dipole_data(2,i,r_i)), i, quant_nums(i,:)
		else if (real(dipole_data(3,i,r_i)) /= 0) then
			write(2, *) r, real(dipole_data(3,i,r_i)), i, quant_nums(i,:)
		else
			 write(2,*) r, real(dipole_data(3,i,r_i)), i, quant_nums(i,:)
		end if
	end do
write(2,*)
end do
close(2)

open(2, file = "./dips_imag.dat")
do i = 1, num_channels
	do r_i = 1, interval_num + 1
		r = r_min + (r_interval*(r_i-1))
		if (aimag(dipole_data(1,i,r_i)) /= 0) then
			write(2, *) r, aimag(dipole_data(1,i,r_i)), i, quant_nums(i,:)
		else if (aimag(dipole_data(2,i,r_i)) /= 0) then
			write(2, *) r, aimag(dipole_data(2,i,r_i)), i, quant_nums(i,:)
		else if (aimag(dipole_data(3,i,r_i)) /= 0) then
			write(2, *) r, aimag(dipole_data(3,i,r_i)), i, quant_nums(i,:)
		else
			write(2,*) r, aimag(dipole_data(3,i,r_i)), i, quant_nums(i,:)
		end if
	end do
write(2,*)
end do
close(2)


allocate(j_list(vib_num))

do j = 1, num_channels

	if (quant_nums(j,2) == 0) then

		rot_ex = (/quant_nums(j,1), quant_nums(j,3), quant_nums(j,4)/)
		print *, rot_ex
		vib_counter = 1
		do i = 1, num_channels
			if (rot_ex(1) == quant_nums(i,1) .AND. rot_ex(2) == quant_nums(i,3) .AND. rot_ex(3) == quant_nums(i,4)) then
				j_list(vib_counter) = i
				vib_counter = vib_counter + 1
				!print *, j_list
			else
				cycle
			end if
		end do

		do r_i = 1, interval_num + 1 
			do coord_i = 1,3
					if (r_i > 1) then

						call find_cderiv(dipole_data(coord_i, j, r_i), dipole_data(coord_i, j, r_i-1), r_interval, deriv_1)
						call find_cderiv(-dipole_data(coord_i, j, r_i), dipole_data(coord_i, j, r_i-1), r_interval, deriv_2)

							if (r_i == 2) then
								if (abs(deriv_1) >= abs(deriv_2)) then
									do i = 1, vib_num
										dipole_data(coord_i, j_list(i), r_i) = -dipole_data(coord_i, j_list(i), r_i)
									end do
								end if

							else
								if (abs(deriv_1 - old_deriv) >= abs(deriv_2 - old_deriv)) then
									do i = 1, vib_num
										dipole_data(coord_i, j_list(i), r_i) = -dipole_data(coord_i, j_list(i), r_i)
									end do
								end if
							end if

							call find_cderiv(dipole_data(coord_i, j, r_i), dipole_data(coord_i, j, r_i-1), r_interval, old_deriv)
						
					end if
			end do
		end do
	end if
end do

open(2, file = "./dips_real_exc.dat")
counter = 0
do i = 1, num_channels
	if (quant_nums(i,1) >= 2) then
		do a = 1,3
			counter = counter + 1
			exc_channels(:,counter) = (/quant_nums(i,:), a/)
			write(2,*) "#counter,n,l,m =", counter, quant_nums(i,1), quant_nums(i,3), quant_nums(i,4), a
			do r_i = 1, interval_num + 1
				r = r_min + (r_interval*(r_i-1))
				write(2,*) r, real(dipole_data(a,i,r_i))
			end do
			write(2,*)
		end do
		write(2,*)
	end if
end do
close(2)

open(2, file = "./dips_real_ground.dat")
counter = 0
do i = 1, num_channels
	if (quant_nums(i,1) == 1) then
		do a = 1,3
			counter = counter + 1
			ground_channels(:,counter) = (/quant_nums(i,:), a/)
			write(2,*) "#counter,n,l,m =", counter, quant_nums(i,1), quant_nums(i,3), quant_nums(i,4), a
			do r_i = 1, interval_num + 1
				r = r_min + (r_interval*(r_i-1))
				write(2,*) r, real(dipole_data(a,i,r_i))
			end do
			write(2,*)
		end do
		write(2,*)
	end if
end do
close(2)

open(2, file = "./dips_imag_exc.dat")
counter = 0
do i = 1, num_channels
	if (quant_nums(i,1) >= 2) then
		do a = 1,3
			counter = counter + 1
			write(2,*) "#counter,n,l,m =", counter, quant_nums(i,1), quant_nums(i,3), quant_nums(i,4), a
			do r_i = 1, interval_num + 1
				r = r_min + (r_interval*(r_i-1))
				write(2,*) r, aimag(dipole_data(a,i,r_i))
			end do
			write(2,*)
		end do
		write(2,*)
	end if
end do
close(2)

open(2, file = "./dips_imag_ground.dat")
counter = 0
do i = 1, num_channels
	if (quant_nums(i,1) == 1) then
		do a = 1,3
			counter = counter + 1
			write(2,*) "#counter,n,l,m =", counter, quant_nums(i,1), quant_nums(i,3), quant_nums(i,4), a
			do r_i = 1, interval_num + 1
				r = r_min + (r_interval*(r_i-1))
				write(2,*) r, aimag(dipole_data(a,i,r_i))
			end do
			write(2,*)
		end do
		write(2,*)
	end if
end do
close(2)


open(2, file = "./dips_abs_ground.dat")
counter = 0
do i = 1, num_channels
	if (quant_nums(i,1) == 1 .and. quant_nums(i,2) == 0) then
		do a = 1,3
			counter = counter + 1
			write(2,*) "#counter,n,l,m =", counter, quant_nums(i,1), quant_nums(i,3), quant_nums(i,4), a
			do r_i = 1, interval_num + 1
				r = r_min + (r_interval*(r_i-1))
				write(2,*) r, sqrt(real(dipole_data(a,i,r_i))**2 + aimag(dipole_data(a,i,r_i))**2)
			end do
			write(2,*)
		end do
		write(2,*)
	end if
end do
close(2)

open(2, file = "./dips_abs_exc.dat")
counter = 0
do i = 1, num_channels
	if (quant_nums(i,1) >= 2 .and. quant_nums(i,2) == 0) then
		do a = 1,3
			counter = counter + 1
			write(2,*) "#counter,n,l,m =", counter, quant_nums(i,1), quant_nums(i,3), quant_nums(i,4), a
			do r_i = 1, interval_num + 1
				r = r_min + (r_interval*(r_i-1))
				write(2,*) r, sqrt(real(dipole_data(a,i,r_i))**2 + aimag(dipole_data(a,i,r_i))**2)
			end do
			write(2,*)
		end do
		write(2,*)
	end if
end do
close(2)

print *, "----------------------------------------------------------------------------------------------------"
print *, "GET INTERPOLATED DIPOLE DATA..."
print *, "----------------------------------------------------------------------------------------------------"

call get_interp_grid

allocate(interp_dips(3, num_channels, num_interp_points), slopes(3,num_channels), beg_slopes(3,num_channels), end_slopes(3,num_channels))

print *, interp_geoms

do r_i = 1, num_interp_points
	if (r_i >= 4 .and. r_i <= 16) then
		interp_dips(:,:,r_i) = dipole_data(:,:,r_i-3)
	else 
		interp_dips(:,:,r_i) = 0
	end if
end do

do i = 1, num_channels
	do j = 1, 3
		slopes(j,i) = (dipole_data(j,i,11) - dipole_data(j,i,3)) / (interp_geoms(14) - interp_geoms(6))
	end do
end do

do i = 1, num_channels
	do j = 1, 3
		beg_slopes(j,i) = (dipole_data(j,i,2) - dipole_data(j,i,1)) / (interp_geoms(5) - interp_geoms(4))
	end do
end do

do i = 1, num_channels
	do j = 1, 3
		end_slopes(j,i) = (dipole_data(j,i,13) - dipole_data(j,i,12)) / (interp_geoms(16) - interp_geoms(15))
	end do
end do

do i = 1, num_channels
	if ( (quant_nums(i,1) == 3 .and. quant_nums(i,3) == 1 .and. quant_nums(i,4) == 1) .or. &
		& (quant_nums(i,1) == 3 .and. quant_nums(i,3) == 2 .and. quant_nums(i,4) == 1) .or. &
		& (quant_nums(i,1) == 2 .and. quant_nums(i,3) == 2 .and. quant_nums(i,4) == -1) .or. &
		& (quant_nums(i,1) == 1 .and. quant_nums(i,3) == 2 .and. quant_nums(i,4) == 2) .or. &
		& (quant_nums(i,1) == 2 .and. quant_nums(i,3) == 1 .and. quant_nums(i,4) == -1) ) then 
		do j = 1, 3
			do r_i = 1, num_interp_points
				r = interp_geoms(r_i)
				interp_dips(j,i,r_i) = dipole_data(j,i,3) + slopes(j,i) * (r - (r_min + r_interval*3))
			end do
		end do
	else
		do j = 1, 3
			do r_i = 1, 3
				r = interp_geoms(r_i)
				interp_dips(j,i,r_i) = dipole_data(j,i,1) + beg_slopes(j,i) * (r - interp_geoms(4))
			end do
			do r_i = 17, 19
				r = interp_geoms(r_i)
				interp_dips(j,i,r_i) = dipole_data(j,i,13) + end_slopes(j,i) * (r - interp_geoms(16))
			end do
		end do
	end if
end do


! open(333, file = "./Interpolated Dipoles/real_ground")
! open(444, file = "./Interpolated Dipoles/imag_ground")

! eof = 0
! counter = 0
! print *, size(ground_channels,1)
! do while (eof == 0)
! 	counter = counter + 1
! 	do i = 1, num_channels
! 		if (all(quant_nums(i,:)==ground_channels(1:4,counter))) then
! 			print *, ground_channels(:,counter)
! 			index = i
! 			Exit
! 		end if
! 	end do
! 	do i = 1, num_interp_points
! 		read(333,*,iostat=eof) r, temp1
! 		!print *, r, temp1
! 		read(444,*,iostat=eof) r, temp2
! 		interp_dips(ground_channels(5,counter),index,i) = complex(temp1, temp2)
! 	end do
! end do
! close(333)
! close(444)

! open(333, file = "./Interpolated Dipoles/real_exc")
! open(444, file = "./Interpolated Dipoles/imag_exc")

! eof = 0
! counter = 0
! do while (eof == 0)
! 	counter = counter + 1
! 	do i = 1, num_channels
! 		print *, quant_nums(i,:)
! 		!print *, exc_channels(1:4,counter)
! 		if (all(quant_nums(i,:)==exc_channels(1:4,counter))) then
! 			index = i
! 			Exit
! 		end if
! 	end do
! 	do i = 1, num_interp_points
! 		read(333,*,iostat=eof) r, temp1
! 		read(444,*,iostat=eof) r, temp2
! 		interp_dips(exc_channels(5,counter),index,i) = complex(temp1, temp2)
! 	end do
! end do
! close(333)
! close(444)

open(2, file = "./interp_dips_real_exc.dat")
counter = 0
do i = 1, num_channels
	if (quant_nums(i,1) >= 2) then
		do a = 1,3
			counter = counter + 1
			write(2,*) "#counter,n,l,m =", counter, quant_nums(i,1), quant_nums(i,3), quant_nums(i,4), a
			do r_i = 1, num_interp_points
				r = interp_geoms(r_i)
				write(2,*) r, real(interp_dips(a,i,r_i))
			end do
			write(2,*)
		end do
		write(2,*)
	end if
end do
close(2)

open(2, file = "./interp_dips_real_ground.dat")
counter = 0
do i = 1, num_channels
	if (quant_nums(i,1) == 1) then
		do a = 1,3
			counter = counter + 1
			write(2,*) "#counter,n,l,m =", counter, quant_nums(i,1), quant_nums(i,3), quant_nums(i,4), a
			do r_i = 1, num_interp_points
				r = interp_geoms(r_i)
				write(2,*) r, real(interp_dips(a,i,r_i))
			end do
			write(2,*)
		end do
		write(2,*)
	end if
end do
close(2)

open(2, file = "./interp_dips_imag_exc.dat")
counter = 0
do i = 1, num_channels
	if (quant_nums(i,1) >= 2) then
		do a = 1,3
			counter = counter + 1
			write(2,*) "#counter,n,l,m =", counter, quant_nums(i,1), quant_nums(i,3), quant_nums(i,4), a
			do r_i = 1, num_interp_points
				r = interp_geoms(r_i)
				write(2,*) r, aimag(interp_dips(a,i,r_i))
			end do
			write(2,*)
		end do
		write(2,*)
	end if
end do
close(2)

open(2, file = "./interp_dips_imag_ground.dat")
counter = 0
do i = 1, num_channels
	if (quant_nums(i,1) == 1) then
		do a = 1,3
			counter = counter + 1
			write(2,*) "#counter,n,l,m =", counter, quant_nums(i,1), quant_nums(i,3), quant_nums(i,4), a
			do r_i = 1, num_interp_points
				r = interp_geoms(r_i)
				write(2,*) r, aimag(interp_dips(a,i,r_i))
			end do
			write(2,*)
		end do
		write(2,*)
	end if
end do
close(2)

print *, "----------------------------------------------------------------------------------------------------"
print *, "EXTRACTING WAVEFUNCTION VALUES..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(neut_WF_array(num_interp_points, 0:2), ion_WF_array(num_interp_points, vib_num, num_targ))
call read_neut_WF(neut_WF_array)
call read_ion_WF(ion_WF_array)

open(333, file = "neutral_wavefunctions.dat")
do i = 0, 2
	write(333,*) "#v = ", i
	do j = 1, num_interp_points
		r = r_min + r_interval*(j-1)
		write(333,*) r, neut_WF_array(j,i)
	end do
	write(333,*)
end do

open(333, file = "ion_wavefunctions.dat")
do k = 1, num_targ
	do i = 1, vib_num
		write(333,*) "#n = ", k, "v = ", i
		do j = 1, num_interp_points
			r = r_min + r_interval*(j-1)
			write(333,*) r, ion_WF_array(j,i,k)
		end do
		write(333,*)
	end do
end do

print *, "----------------------------------------------------------------------------------------------------"
print *, "APPROXIMATING INTEGRAL..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(dip_vect(3, num_channels, 0:2))

do ch_vib = 0,2
	do j = 1, num_channels

			n = quant_nums(j,1)
			vib_i = quant_nums(j,2)

			x_comp_integrated = 0
			y_comp_integrated = 0
			z_comp_integrated = 0

			do r_i = 1, num_interp_points

				x_comp = interp_dips(1, j, r_i) * neut_WF_array(r_i, ch_vib) * ion_WF_array(r_i, vib_i+1, n) * int_r_interval
				x_comp_integrated = x_comp_integrated + x_comp

				y_comp = interp_dips(2, j, r_i) * neut_WF_array(r_i, ch_vib) * ion_WF_array(r_i, vib_i+1, n) * int_r_interval
				y_comp_integrated = y_comp_integrated + y_comp

				z_comp = interp_dips(3, j, r_i) * neut_WF_array(r_i, ch_vib) * ion_WF_array(r_i, vib_i+1, n) * int_r_interval
				! if (n == 1 .AND. l == 0 .AND. m == 0) then
				! 	print *, j
				! 	print *, "r_i =", r_i
				! 	print *, "dip =", dipole_data(3, j, r_i)
				! 	print *, "z_comp =", z_comp
				! 	print *, "ion WF =", ion_WF_array(r_i, vib_i+1, n)
				! end if
				z_comp_integrated = z_comp_integrated + z_comp

			!sum over r do
			end do

			dip_vect(1,j,ch_vib) = x_comp_integrated
			dip_vect(2,j,ch_vib) = y_comp_integrated
			dip_vect(3,j,ch_vib) = z_comp_integrated

	!sym do
	end do
end do

do ch_vib = 0,2
	write(vib_str, "(I1)") ch_vib
	open(2, file = "./integrated_dips_v0.dat")
	do i = 1, num_channels
		write(2, *) i, dip_vect(1,i,0), dip_vect(2,i,0), dip_vect(3,i,0)
	end do
	close(2)

	open(2, file = "./dips_integrated_imag_ground_"// vib_str //".dat")
	do i = 1, num_channels
		n = quant_nums(i,1)
		l = quant_nums(i,3)
		m = quant_nums(i,4)
			if (n == 1) then
				do k = 0,4
					do j = 1, num_channels
						if (quant_nums(j,1) == n .AND. quant_nums(j,2) == k .AND. quant_nums(j,3) == l .AND. quant_nums(j,4) == m) then
							if (aimag(dip_vect(1,i,0)) /= 0) then
								write(2, *) k, aimag(dip_vect(1,j,ch_vib))
							else if (aimag(dip_vect(2,i,ch_vib)) /= 0) then
								write(2, *) k, aimag(dip_vect(2,j,ch_vib))
							else if (aimag(dip_vect(3,i,0)) /= 0) then
								write(2, *) k, aimag(dip_vect(3,j,ch_vib))
							else
								write(2,*) k, aimag(dip_vect(3,j,ch_vib))
							end if
						end if
					end do
				end do
			write(2,*)
			end if
	end do
	close(2)


	open(2, file = "./dips_integrated_imag_exc_"// vib_str //".dat")
	do i = 1, num_channels
		n = quant_nums(i,1)
		l = quant_nums(i,3)
		m = quant_nums(i,4)
			if (n == 2 .OR. n==3) then
				do k = 0,4
					do j = 1, num_channels
						if (quant_nums(j,1) == n .AND. quant_nums(j,2) == k .AND. quant_nums(j,3) == l .AND. quant_nums(j,4) == m) then
							if (aimag(dip_vect(1,i,0)) /= 0) then
								write(2, *) k, aimag(dip_vect(1,j,ch_vib))
							else if (aimag(dip_vect(2,i,0)) /= 0) then
								write(2, *) k, aimag(dip_vect(2,j,ch_vib))
							else if (aimag(dip_vect(3,i,0)) /= 0) then
								write(2, *) k, aimag(dip_vect(3,j,ch_vib))
							else
								write(2,*) k, aimag(dip_vect(3,j,ch_vib))
							end if
						end if
					end do
				end do
			write(2,*)
			end if
	end do
	close(2)


	open(2, file = "./dips_integrated_real_ground_"// vib_str //".dat")
	do i = 1, num_channels
		n = quant_nums(i,1)
		l = quant_nums(i,3)
		m = quant_nums(i,4)
			if (n == 1) then
				do k = 0,4
	!				if (k == 0) then
	!					write(2,*) n, l, m
	!				end if
					do j = 1, num_channels
						if (quant_nums(j,1) == n .AND. quant_nums(j,2) == k .AND. quant_nums(j,3) == l .AND. quant_nums(j,4) == m) then
							if (real(dip_vect(1,i,0)) /= 0) then
								write(2, *) k, real(dip_vect(1,j,ch_vib))
							else if (real(dip_vect(2,i,0)) /= 0) then
								write(2, *) k, real(dip_vect(2,j,ch_vib))
							else if (real(dip_vect(3,i,0)) /= 0) then
								write(2, *) k, real(dip_vect(3,j,ch_vib))
							else
								write(2,*) k, real(dip_vect(3,j,ch_vib))
							end if
						end if
					end do
				end do
			write(2,*)
			end if
	end do
	close(2)


	open(2, file = "./dips_integrated_real_exc_"// vib_str //".dat")
	do i = 1, num_channels
		n = quant_nums(i,1)
		l = quant_nums(i,3)
		m = quant_nums(i,4)
			if (n == 2 .OR. n == 3) then
				do k = 0,4
					do j = 1, num_channels
						if (quant_nums(j,1) == n .AND. quant_nums(j,2) == k .AND. quant_nums(j,3) == l .AND. quant_nums(j,4) == m) then
							if (real(dip_vect(1,i,0)) /= 0) then
								write(2, *) k, real(dip_vect(1,j,ch_vib))
							else if (real(dip_vect(2,i,0)) /= 0) then
								write(2, *) k, real(dip_vect(2,j,ch_vib))
							else if (real(dip_vect(3,i,0)) /= 0) then
								write(2, *) k, real(dip_vect(3,j,ch_vib))
							else
								write(2,*) k, real(dip_vect(3,j,ch_vib))
							end if
						end if
					end do
				end do
			write(2,*)
			end if
	end do
	close(2)
end do



print *, "----------------------------------------------------------------------------------------------------"
print *, "CHANNEL ELIMINATION..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(energy_array(num_channels),channel_energies(num_targ, 0:vib_num-1), max_chans(energy_step_num, 0:2, 3), nu(num_channels, energy_step_num, 0:2))
open(1, file = "./vibronic_channel_energies.txt")

print *, "Reading channel energies..."

do i = 0,(num_targ*vib_num)-1
	read(1,*) n, vib_i, channel_energies(n,vib_i)
	channel_energies(n,vib_i) = ev_au*channel_energies(n,vib_i)
end do
close(1)

open(1, file = "./energies.dat")
do i = 1, num_channels

		n = quant_nums(i,1)
		vib_i = quant_nums(i,2)
		energy_array(i) = channel_energies(n,vib_i)

		write(1,*) n, vib_i, energy_array(i)/ev_au + diss_en
end do
close(1)

print *, ""
print *, "		GROUND STATE ENERGIES"
print *, ""
print *, "	   ", "n	       ", "v	         ", "Threshold"
print *, "_____________________________________________________________________"
print *, ""
do a = 1, num_channels
if (quant_nums(a, 1) == 1 .AND. quant_nums(a,3) == 0 .AND. quant_nums(a,4) == 0) then
	print *, quant_nums(a,1), quant_nums(a,2), "	", energy_array(a)/ev_au, "(eV)"
end if
end do
print*, ""
print *, "		EXCITED STATE ENERGIES"
print *, ""
print *, "	   ", "n	       ", "v	         ", "Threshold"
print *, "_____________________________________________________________________"
print *, ""
do a = 1, num_channels
if (quant_nums(a, 1) /= 1 .AND. quant_nums(a,3) == 0 .AND. quant_nums(a,4) == 0) then
	print *, quant_nums(a,1), quant_nums(a,2), "	", energy_array(a)/ev_au, "(eV)"
end if
end do

allocate(TCS_array(energy_step_num, 2, 0:2), SPES_array(energy_step_num, 2, 0:2), SPES_array_ground(energy_step_num, 2, 0:2), Z_val(energy_step_num, 0:2, 3))
SPES_array = 0d0
nu = 0d0
do ch_vib = 0,2
do i = 1, energy_step_num

	
	energy_begin = channel_energies(1,0) 
	energy_end = channel_energies(num_targ, vib_num-1)
	energy_step = (energy_end - energy_begin)/energy_step_num

	!print *, energy_begin
	E = energy_begin + i*energy_step

	num_open = 0
        Do a=1,num_channels
          If(energy_array(a)>E)Exit
        Enddo
        num_open=a-1
		!print *, num_open, delta_Es(ch_vib)

	if (num_open /= num_channels) then

		num_closed = num_channels - num_open
		allocate(dipole_open(3,num_open),dipole_closed(3,num_closed),IPIV(num_channels))
		allocate(smat_cc(num_closed,num_closed),smat_co(num_closed,num_open),beta(num_channels), Z_norms(num_closed), D(num_closed, num_open))
		dipole_closed = dip_vect(:,num_open+1:num_channels, ch_vib)
		dipole_open = dip_vect(:,1:num_open, ch_vib)
		smat_cc = Conjg(Transpose(smat(num_open+1:num_channels,num_open+1:num_channels)))
		smat_co = Conjg(Transpose(smat(1:num_open,num_open+1:num_channels)))

		beta = 0d0
		do a = 1, num_closed
			beta(a+num_open) = pi/sqrt(2*(energy_array(a+num_open)-E))
		end do

		do j = 1, num_closed
			!print *, energy_array(j+num_open), E + delta_Es(ch_vib)
			nu(j+num_open,i,ch_vib) = 1/sqrt(2*(energy_array(j+num_open)-E))
		end do

		allocate(dipole_phys(3,num_open))
		dipole_phys = 0d0

		do a = 1, num_closed
			smat_cc(a,a) = smat_cc(a,a) - exp(2d0*ci*beta(a+num_open))
		end do

		call ZGESV(num_closed,num_open,smat_cc,num_closed,IPIV(1:num_closed),smat_co,num_closed,INFO)

		!do coord_i = 1,3
			do k = 1, num_open
				do j = 1, num_closed
					D(j, k) = smat_co(j,k) * nu(j+num_open, i, ch_vib)**1.5
				end do
			end do
		!end do

		Z_norms = 0
		do j = 1, num_closed
			!do coord_i = 1, 3
				Z_norms(j) = Z_norms(j) + sum(abs(D(j,:))**2)
			!end do
		end do
 
		!print *, maxloc(Z_norms, 1) + num_open
		do j = 1, 3
			max_chans(i, ch_vib, j) = maxloc(Z_norms, 1) + num_open
			Z_val(i, ch_vib, j) = Z_norms(maxloc(Z_norms, 1))
			Z_norms(maxloc(Z_norms, 1)) = 0
		end do

 		do a = 1,3
			dipole_phys(a,:) = dipole_open(a,:) - MatMul(dipole_closed(a,:),smat_co)
		end do
		deallocate(dipole_open, dipole_closed, smat_cc, smat_co, IPIV, beta, Z_norms, D)

		if (i == 10) then
			print *, ""
			print *, "HALF WAY SUM"
			print *, "__________________________"
			print *, sum(abs(dipole_phys(1,:))**2)+sum(abs(dipole_phys(2,:))**2)+sum(abs(dipole_phys(3,:))**2)
			print *, ""
		end if

	else if (num_open == num_channels) then

	num_open = num_channels
	allocate(dipole_phys(num_open, 3))
	dipole_phys = dip_vect(:,:,ch_vib)

	end if

	omega = omega_0 + ( E - channel_energies(1,0) )
	!omega = omega_0 - channel_energies(0)

	TCS_array(i,1,ch_vib) = const*omega*(sum(abs(dipole_phys(1,:))**2)+sum(abs(dipole_phys(2,:))**2)+sum(abs(dipole_phys(3,:))**2))
	TCS_array(i,2,ch_vib) = ( E - delta_Es(ch_vib) ) / ev_au + diss_en
	SPES_array(i,2,ch_vib) = TCS_array(i,2,ch_vib)

	do n = 1, num_targ
        do vib_i = 0, vib_num-1
            if ((abs(E - channel_energies(n,vib_i)) < elec_en_range) .AND. (E - channel_energies(n,vib_i) >= 0)) then
                do j = 1, num_channels
                    if (quant_nums(j,1) == n .AND. quant_nums(j,2) == vib_i) then
                        SPES_array(i,1,ch_vib) = SPES_array(i,1,ch_vib) + const*omega*(abs(dipole_phys(1,j))**2+abs(dipole_phys(2,j))**2+abs(dipole_phys(3,j))**2)
                    end if
                end do
            end if
        end do
    end do


	deallocate(dipole_phys)

end do
end do
close(1)

open(1, file = "./TCS.dat")
do i = 1, energy_step_num

	write(1,*) TCS_array(i,2,0), TCS_array(i,1,0), TCS_array(i,2,1), TCS_array(i,1,1), TCS_array(i,2,2), TCS_array(i,1,2)

end do
close(1)

open(1, file = "./SPES.dat")
do i = 1, energy_step_num

	write(1,*) SPES_array(i,2,0), SPES_array(i,1,0), SPES_array(i,2,1), SPES_array(i,1,1), SPES_array(i,2,2), SPES_array(i,1,2)

end do
close(1)

E_min = minval(TCS_array(:,2,:))
E_max = maxval(TCS_array(:,2,:))

allocate(conv_TCS_array(num_conv_points,2,0:2,2), conv_SPES_array(num_conv_points,2,0:2,2), conv_SPES_array_ground(num_conv_points,2,0:2,2))
do ch_vib = 0,2
    conv_size_array(1) = size(TCS_array,1)
    conv_size_array(2) = size(TCS_array,2)
    call cauchy_conv(E_min, E_max, conv_size_array, TCS_array(:,:,ch_vib), conv_TCS_array(:,:,ch_vib,1))
    call gauss_conv(E_min, E_max, conv_size_array, TCS_array(:,:,ch_vib), conv_TCS_array(:,:,ch_vib,2))
    call cauchy_conv(E_min, E_max, conv_size_array, SPES_array(:,:,ch_vib), conv_SPES_array(:,:,ch_vib,1))
    call gauss_conv(E_min, E_max, conv_size_array, SPES_array(:,:,ch_vib), conv_SPES_array(:,:,ch_vib,2))
end do

open(1, file = "./conv_data_v_lorentz.dat")
do i = 1, num_conv_points
        write(1,*) conv_TCS_array(i,2,0,1), conv_TCS_array(i,1,0,1), conv_TCS_array(i,1,1,1), conv_TCS_array(i,1,2,1)
end do
close(1)

open(1, file = "./conv_data_v_gauss.dat")
do i = 1, num_conv_points
        write(1,*) conv_TCS_array(i,2,0,2), conv_TCS_array(i,1,0,2), conv_TCS_array(i,1,1,2), conv_TCS_array(i,1,2,2)
end do
close(1)

open(1, file = "./conv_cross_section_summed.dat")
do i = 1, num_conv_points

	write(1,*) conv_TCS_array(i,2,0,2), weights(0)*conv_TCS_array(i,1,0,2)+weights(1)*conv_TCS_array(i,1,1,2)+weights(2)*conv_TCS_array(i,1,2,2)

end do
close(1)


!open(1, file = "./SPES.dat")
!do i = 1, num_conv_points

!	write(1,*) new_conv_TCS_array(i,2,0), weights(0)*new_conv_TCS_array(i,1,0)+weights(1)*new_conv_TCS_array(i,1,1)+weights(2)*new_conv_TCS_array(i,1,2)

!end do
!close(1)

open(1, file = "./SPES_data_v_lorentz.dat")
do i = 1, num_conv_points
        write(1,*) conv_SPES_array(i,2,0,1), conv_SPES_array(i,1,0,1), conv_SPES_array(i,1,1,1), conv_SPES_array(i,1,2,1)
end do
close(1)

open(1, file = "./SPES_data_v_gauss.dat")
do i = 1, num_conv_points
        write(1,*) conv_SPES_array(i,2,0,2), conv_SPES_array(i,1,0,2), conv_SPES_array(i,1,1,2), conv_SPES_array(i,1,2,2)
end do
close(1)

open(27, file = "Resonances_v0.dat")
do e_i = 1, energy_step_num
  E = energy_begin + (e_i-1)*energy_step
  !print *, e_i, max_chans(e_i,1)
  write(27,'(F11.8,3(4I4,2F10.3))') E/ev_au + diss_en, quant_nums(max_chans(e_i,0,1),1), quant_nums(max_chans(e_i,0,1),2), quant_nums(max_chans(e_i,0,1),3), quant_nums(max_chans(e_i,0,1),4), Z_val(e_i,0,1), nu(max_chans(e_i,0,1), e_i, 0),&
				& quant_nums(max_chans(e_i,0,2),1), quant_nums(max_chans(e_i,0,2),2), quant_nums(max_chans(e_i,0,2),3), quant_nums(max_chans(e_i,0,2),4), Z_val(e_i,0,2), nu(max_chans(e_i,0,2), e_i, 0), &
				& quant_nums(max_chans(e_i,0,3),1), quant_nums(max_chans(e_i,0,3),2), quant_nums(max_chans(e_i,0,3),3), quant_nums(max_chans(e_i,0,3),4), Z_val(e_i,0,3), nu(max_chans(e_i,0,3), e_i, 0)
  write(27,*)
end do
close(27)


state_range = 0
do i = 1, energy_step_num
	if (TCS_array(i,2,0) <= -2.835) then
		state_range = state_range +1
	end if
end do

!state_avg(1) = sum(TCS_array(1:state_range,1,0))/state_range
!state_avg(2) = sum(TCS_array(state_range:energy_step_num,1,0))/(energy_step_num - state_range)

!print *, ""
!print *, "	AVERAGE TCS"
!print *, "_________________________"
!print *, ""
!print *, "GS:	", state_avg(1)
!print *, "ES:	", state_avg(2)

end program
