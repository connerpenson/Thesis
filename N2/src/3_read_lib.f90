module read_lib

use constants
use math_ops, only: frame_transform

implicit none

contains

subroutine read_smat(file_path, num_channels, smat)

    integer, intent(in)                         ::  num_channels
	character (len = 256), intent(in)			::	file_path

	integer								::	i
	real(8), dimension(2, num_channels)				::	tmp

	complex(8), dimension(num_channels, num_channels), intent(out)	::	smat

open(1, file = trim(file_path))
do i = 1, num_channels
        read(1,*) tmp(1:2, 1:i)
        smat(1:i, i) = tmp(1, 1:i) + tmp(2, 1:i)*ci
        smat(i, 1:i-1) = smat(1:i-1,i)
end do
close(1)

end subroutine

subroutine read_dip_data(dipole_data, num_channels, quant_nums)

	use constants, only: dip_en_list, num_waves

    integer, intent(in)                 ::  num_channels
	complex(8), dimension(3, num_channels,interval_num+1), intent(out)	::	dipole_data

	integer								::	i, sym_i, ang_i, ang_chan, r_i, block_size, channel_order, i_Q
	integer								::	n_josh, l_josh, m_josh, v, n, a, l, m, skip_num, j
	real(8)								::	energy
	character(1)	::	interval_str
	character(45)				::	file_name
	character (len = 2)						    ::	sym_str
	complex(8), allocatable, dimension(:,:)					::	mom_array
	integer, allocatable, dimension(:,:)					::	order_list

	real(8), dimension(20)						::	temp2
	real(8), dimension(interval_num+1)			::	exc_levels
	integer, dimension(num_channels,4), intent(out)				::	quant_nums


call get_quant_nums(quant_nums)

do r_i = 1, interval_num+1

	if (r_i /= 10) then
		write(interval_str, '(I1)') r_i
	else
		write(interval_str, '(I1)') 1
	end if

	if (r_i /= 10) then
		open(1, file = "./Dipole_Data_CI/new/0" // interval_str // "/B1u.pwDipoles.formatted", status = "old")
	else 
		open(1, file = "./Dipole_Data_CI/new/" // interval_str // "0/B1u.pwDipoles.formatted", status = "old")
	end if

	do a = 1, 16
		read(1,*)
	end do

	read(1,*) temp2(1:7), exc_levels(r_i)
	close(1)
	exc_levels(r_i) = .5*exc_levels(r_i)
	print *, exc_levels(r_i)
end do

!dipole_data = 0

do sym_i = 1,3

	skip_num = num_waves(sym_i)

    if (sym_i == 1) then
        ang_chan = 9
    else if (sym_i == 2 .OR. sym_i == 3) then
        ang_chan = 12
    end if

	allocate(order_list(ang_chan,3), mom_array(ang_chan,3))

	sym_str = sym_strings(sym_i)

		do r_i = 1, interval_num+1

			if (r_i /= 10) then
				write(interval_str, '(I1)') r_i
			else
				write(interval_str, '(I1)') 1
			end if
		
			if (r_i /= 10) then
				open(1, file = "./Dipole_Data_CI/new/0" // interval_str // "/" // sym_str // "u.pwDipoles.formatted", status = "old")
			else 
				open(1, file = "./Dipole_Data_CI/new/" // interval_str // "0/" // sym_str // "u.pwDipoles.formatted", status = "old")
			end if

			do a = 1,2
				read(1,*)
			end do

			read(1, *) temp2(1:5), block_size
			!print *, block_size

			do a = 1,7
				read(1,*)
			end do

			do ang_i = 1, ang_chan
				read(1,*) temp2(1), order_list(ang_i,1), temp2(1:3), order_list(ang_i,2), order_list(ang_i,3)
			end do

			energy = exc_levels(r_i) + .00367 !dip_en_list(r_i)

			call dip_mag(order_list, skip_num, energy, block_size, mom_array)

			do j = 1, ang_chan

				if (r_i == interval_num .or. r_i == interval_num + 1) then
					if (order_list(j,1) == 1) then
						order_list(j,1) = 11
					end if
					if (order_list(j,1) == 2) then
						order_list(j,1) = 22
					end if
					if (order_list(j,1) == 3) then
						order_list(j,1) = 33
					end if

					if (order_list(j,1) == 11) then
						order_list(j,1) = 2
					end if
					if (order_list(j,1) == 22) then
						order_list(j,1) = 3
					end if
					if (order_list(j,1) == 33) then
						order_list(j,1) = 1
					end if
				end if

				n = order_list(j,1)
				l = order_list(j,2)
				m = order_list(j,3)

				print *, sym_strings(sym_i), j, n, l, m

				open(10, file = "./Smat/3states_lmax4/vibronic_DR.channel", status = "old")
				read(10, *)
				read(10, *)

				do i = 1, num_channels
					n_josh = quant_nums(i,1)
					l_josh = quant_nums(i,3)
					m_josh = quant_nums(i,4)

					if (n == 1) then
						if (n == n_josh .AND. l == l_josh .AND. m == m_josh) then
							dipole_data(:,i,r_i) = mom_array(j,:)
						end if
					else
						if (n_josh /= 1 .AND. l == l_josh .AND. m == m_josh) then
							dipole_data(:,i,r_i) = mom_array(j,:)
						end if
					end if
				end do
				close(10)
				close(1)
			end do
		end do
	deallocate(order_list, mom_array)
end do


open(50, file = "./dipole_data.dat")
    do a = 1, num_channels
		write(50, *) quant_nums(a,:)
		do i_Q = 1, interval_num + 1
        	write(50, *) dipole_data(:,a,i_Q)
    	end do
		write(50,*)
	end do
close(50)


end subroutine

subroutine dip_mag(order_list, skip_num, energy, block_size, dip_mom)

	integer, allocatable, dimension(:,:), intent(in)	:: order_list
	real(8)							::	res, step, en_1, en_2, cur_energy, des_energy, exc_en
	real(8)							::	x_real, x_imag, z_real, z_imag, y_real, y_imag
	complex(8)						::	x, y, z
	integer							::	ang_i, block_size, a, order, skip_order
	real(8), dimension(0:10000000)				::	temp1
	integer, dimension(2)			::	array_dims

	real(8), intent(in)					::	energy
	integer, intent(in)					::	skip_num

	complex(8), allocatable, dimension(:,:), intent(out)			::	dip_mom
	real(8), dimension(100)					::	temp

	array_dims = shape(order_list)

	allocate(dip_mom(array_dims(1),3))

	do while (order /= skip_num)
		read(1,*) order
	end do
	
	read(1,*)
	read(1,*)
	read(1,*)

	read(1, *) temp(1), en_1
	read(1, *) temp(1), en_2
	!print *, en_2, en_1
	step = en_2 - en_1

	do ang_i = 1, array_dims(1)

		order = 0
		do while (order /= ang_i)
			read(1, *) order
			!print *, order
			if (order == ang_i) then
				res = 10
				do while (abs(res) > step)
					read(1, *) temp(1), cur_energy, y_real, y_imag, z_real, z_imag, x_real, x_imag
					res = energy - cur_energy
					!print *, res
				end do
			else
				do a = 1, block_size
					read(1,*) skip_order
					if (skip_order == ang_i) then
						exit
					end if
				end do
			end if
		end do

		x = cmplx(x_real, x_imag, 8)
		y = cmplx(y_real, y_imag, 8)
		z = cmplx(z_real, z_imag, 8)

		!print *, x, y, z

		dip_mom(ang_i, 1) = x
		dip_mom(ang_i, 2) = y
		dip_mom(ang_i, 3) = z

	end do

!print *, x_real, x_imag, y_real, y_imag, z_real, z_ima

end subroutine dip_mag

! subroutine smat_ft(trans_smat)

!     integer                 ::  a, l, l_prime, m, m_prime, n, n_prime, v, v_prime, i, j, eof, r_i, pre_ft_num_channels, vib_i, vib_i_prime
!     integer, dimension(3)   ::  nums, prime_nums
!     integer, dimension(2, vib_num, vib_num)   ::  inds
!     integer, allocatable    ::  quant_nums(:,:), quant_mat(:,:,:)
!     complex                 ::  temp
!     real(8)                 ::  skip, temp_r, temp_i, int_sum

!     real(8),    allocatable ::  ion_WF_array(:,:,:)

!     complex(8), allocatable ::  smat(:,:,:)
! 	complex(8), intent(out), allocatable	::	trans_smat(:,:)
!     character(len = 50)     ::  smat_file
!     character(len = 25)     ::  channel_file
!     character(len = 10)     ::  read_string

!     pre_ft_num_channels = num_channels/vib_num

!     allocate(smat(num_channels, num_channels, lin_size), quant_nums(pre_ft_num_channels, 3), quant_mat(num_channels, num_channels, 8))
!     allocate(ion_WF_array(lin_size, vib_num, num_targ))
!     allocate(trans_smat(num_channels, num_channels))

!     call read_ion_WF(ion_WF_array)

!     smat_file = "./Smat/Re_Smatr_interpol_new_grid_continuity_Q.dat"
!     channel_file = "./Smat/smat_ft.channel"

!     open(2, file = channel_file)
!     read(2,*)
!     do i = 1, pre_ft_num_channels
!         read(2,*) n, l, m
!         quant_nums(i,:) = (/n,l,m/)
!        !print *, quant_nums(i,:)
!     end do

! 	do vib_i = 1, vib_num
! 		do vib_i_prime = 1, vib_num
! 			do i = 1, pre_ft_num_channels
! 				do j = 1, pre_ft_num_channels
! 					quant_mat(i + (vib_i -1)*pre_ft_num_channels, j + (vib_i_prime -1)*pre_ft_num_channels, 1)   = quant_nums(i, 1)
! 					quant_mat(i + (vib_i -1)*pre_ft_num_channels, j + (vib_i_prime -1)*pre_ft_num_channels, 2)   = vib_i - 1
! 					quant_mat(i + (vib_i -1)*pre_ft_num_channels, j + (vib_i_prime -1)*pre_ft_num_channels, 3:4) = quant_nums(i, 2:3)
! 					quant_mat(i + (vib_i -1)*pre_ft_num_channels, j + (vib_i_prime -1)*pre_ft_num_channels, 5)   = quant_nums(j, 1)
! 					quant_mat(i + (vib_i -1)*pre_ft_num_channels, j + (vib_i_prime -1)*pre_ft_num_channels, 6)   = vib_i_prime - 1
! 					quant_mat(i + (vib_i -1)*pre_ft_num_channels, j + (vib_i_prime -1)*pre_ft_num_channels, 7:8) = quant_nums(j, 2:3)
! 				end do
! 			end do
! 		end do
! 	end do
!     close(2)


!     open(2, file = smat_file)
!     read(2,*)
!     read(2,*)
!     read(2,*)
!     do
!         read(2,*, IOSTAT = eof) n, l, m, n_prime, l_prime, m_prime
!         if (eof < 0) then
!             exit
!         end if
!         !print *, n, l, m, n_prime, l_prime, m_prime
!         call find_smat_index(quant_mat, (/n,l,m/), (/n_prime,l_prime,m_prime/), inds)
! 		if (all(inds == 0)) then
! 			cycle
! 		end if
!         do j = 1, lin_size
!             read(2,*) skip, temp_r, temp_i
! 			!print *, temp_r, temp_i, inds
!             temp = complex(temp_r,temp_i)
! 			do vib_i = 1, vib_num
! 				do vib_i_prime = 1, vib_num
! 					!print *, "inds:", inds(1,vib_i,vib_i_prime), inds(2,vib_i,vib_i_prime)
!             		smat(inds(1,vib_i,vib_i_prime), inds(2,vib_i,vib_i_prime), j) = temp
! 				end do
! 			end do
!         end do
! 		close(2)
!     end do
 

!     open(2, file = "./smat_r1.dat")
!     do i = 1, pre_ft_num_channels
!         do j = 1, i
!             write (2,*) quant_mat(i,j,1), quant_mat(i,j,2), quant_mat(i,j,3), quant_mat(i,j,4), &
!             & quant_mat(i,j,5), quant_mat(i,j,6), smat(i,j,1)
!         end do
!     end do
	
! 	open(2, file = "quant_mat.dat")
! 	do a = 1, num_channels
! 		write(2,*) quant_mat(1,a,:)
! 	end do
! 	close(2)

! 	do i = 1, num_channels
! 		do j = 1, num_channels

! 			! Get the frame transform quantum numbers from quant mat
! 			n		= quant_mat(i,j,1)
! 			v 		= quant_mat(i,j,2)
! 			n_prime = quant_mat(i,j,5)
! 			v_prime = quant_mat(i,j,6) 

! 			!print *, i, j
! 			!print *, "look here", quant_mat(i,j,1), quant_mat(i,j,2)

! 			!approximate the integral on the data in smat
! 			call frame_transform(smat(i,j,:), n, v, n_prime, v_prime, ion_WF_array, trans_smat(i,j))
! 		end do
! 	end do

! 	open(2, file = "trans_smat.dat")
! 	do i = 1, num_channels

! 			write(2, *) quant_mat(17,i,:), trans_smat(17,i)

! 	end do
! 	close(2)

! 	! open(2, file = "trans_smat_entry.dat")
! 	! call find_smat_index(quant_mat, (/1,4,-4/), (/1,4,-4/), inds)
! 	! write(2, *) trans_smat(inds(1),inds(2),1,1)
! 	! close(2)

!     end subroutine smat_ft

    subroutine find_smat_index(quant_mat, nums, prime_nums, inds)

        integer, intent(in)                 ::  nums(3), prime_nums(3)
        integer, dimension(3), intent(in)   ::  quant_mat(:,:,:)
        integer                             ::  a, b, v, v_prime
        integer, dimension(2, 0:vib_num-1, 0:vib_num-1), intent(out)  ::  inds

		inds = 0

        do a = 1, num_channels
            do b = 1, num_channels
                if ((quant_mat(a,b,1) == nums(1)) .AND. (all(quant_mat(a,b,3:4) == nums(2:3))) .AND. &
					(quant_mat(a,b,5) == prime_nums(1)) .AND. (all(quant_mat(a,b,7:8) == prime_nums(2:3)))) then
                
						v = quant_mat(a, b, 2)
						v_prime = quant_mat(a, b, 6)

						inds(1, v, v_prime) = a
                    	inds(2, v, v_prime) = b
                    
                end if
            end do
        end do

    end subroutine

	subroutine read_smat_ft(S_ft)

		integer :: i, j, k, ind1, ind2
		integer, dimension(8) :: nums
		real(8), dimension(2) :: S_temp
		integer, dimension(num_channels,4) :: quant_nums

		complex(8), dimension(num_channels,num_channels), intent(out) :: S_ft 

		call get_quant_nums(quant_nums)

		open(27, file = "FT_Smat.dat")
		do i = 1, num_channels
			do j = 1, num_channels
				read(27,*)
				read(27,*) nums(:)
				read(27,*) S_temp(1), S_temp(2)

				do k = 1, num_channels
					if (quant_nums(k,1) == nums(1) .and. quant_nums(k,2) == nums(2) .and. quant_nums(k,3) == nums(3) .and. quant_nums(k,4) == nums(4)) then
						ind1 = k
						exit
					end if
				end do
				do k = 1, num_channels
					if (quant_nums(k,1) == nums(5) .and. quant_nums(k,2) == nums(6) .and. quant_nums(k,3) == nums(7) .and. quant_nums(k,4) == nums(8)) then
						ind2 = k
						exit
					end if
				end do
				
				S_ft(ind1,ind2) = complex(S_temp(1), S_temp(2))

			end do
		end do
		close(27)

	end subroutine

end module read_lib
