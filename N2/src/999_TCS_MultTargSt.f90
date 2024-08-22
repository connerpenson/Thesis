program dip_moment

use constants
use read_lib
use math_ops
use print_ops

    implicit none

    real(8)      ::   row_sum, r, energy_begin, energy_end, energy_step, E, omega, diff, energy
    integer      ::   state_range, energy_step_num, a, i, num_open, vib_counter
    integer      ::   vib_i, r_i, n, j, l, m, k, INFO
    integer      ::   num_closed, coord_i, neut_vib
    complex(8)   ::   x_comp, x_comp_integrated, y_comp, y_comp_integrated, z_comp, z_comp_integrated
    complex(8)	 ::   deriv_1, deriv_2, old_deriv
    integer, dimension(3)		:: rot_ex
    integer, dimension(2)		:: conv_size_array

    character (len = 256)		::   smat_file
    character (len = 3), dimension(0:interval_num)	::	geom_list
    real(8), dimension(2)		::	state_avg

    complex(8), allocatable		::	dipole_data(:,:,:), smat(:,:)
    complex(8), allocatable		::	dipole_phys(:,:), dipole_open(:,:), dipole_closed(:,:), smat_cc(:,:), smat_co(:,:)
    complex(8), allocatable		::	beta(:), part_vect(:), dip_vect(:,:,:)
    integer,    allocatable		::	IPIV(:), j_list(:), quant_nums(:,:)
    real(8),    allocatable		::	energy_array(:), neut_WF_array(:,:), ion_WF_array(:,:,:), channel_energies(:,:)
    real(8),	allocatable		::	TCS_array(:,:,:), conv_TCS_array(:,:,:,:), conv_SPES_array(:,:,:,:), conv_SPES_array_ground(:,:,:,:), SPES_array(:,:,:), SPES_array_ground(:,:,:), unitarity_check(:), energies(:), temp_energies(:)

    real(8), dimension(100)     ::  temp

!call determine_num_channels(l_max, num_targ, vib_num, num_channels)
call get_geom_list(r_interval, r_min, interval_num, geom_list)
allocate(quant_nums(num_channels,4), unitarity_check(num_channels), energies(num_targ*vib_num), temp_energies(num_targ*vib_num) )

!remake x

print *, "----------------------------------------------------------------------------------------------------"
print *, "READING S-MATRIX IN..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(smat(num_channels, num_channels))
!if (l_max == 4) then
	!smat_file = "./Smat/3states_lmax4/vibronic_DR.Smat"
!else if (l_max == 2) then
	smat_file = "./Smat/3states_lmax2/vibronic_VE_2states.Smat"
!end if
call read_smat_ft(smat)

print *, "----------------------------------------------------------------------------------------------------"
print *, "UNITARITY CHECK..."
print *, "----------------------------------------------------------------------------------------------------"

open(1, file = "./smat.dat")
do i = 1, num_channels
	row_sum = 0
	do j = 1, num_channels
		row_sum = row_sum + abs(smat(i,j))**2
	end do
	write(1,*) quant_nums(i, 1), quant_nums(i, 2), quant_nums(i, 3), quant_nums(i, 4), row_sum
end do
close(1)

open(1, file = "./unitarity.dat")
write(1,*) "Row", "Sum of Squares"
unitarity_check = 0
do i = 1, num_channels
    do j = 1, num_channels
        unitarity_check(i) = unitarity_check(i) + abs(smat(i,j))**2
    end do
    write(1,*) i, unitarity_check(i)
end do
close(1)
do i = 1, num_channels
    diff = abs(1.00d0 - unitarity_check(i))
    if (diff > unit_tol) then
        print *, "ERROR!"
        print *, "The ", i, " row of the S-matrix is not unitary."
        print *, "Sum of squares:", unitarity_check(i)
    else
        if (i == num_channels) then
            print *, "Unitary!"
        end if
    end if
end do


print *, "----------------------------------------------------------------------------------------------------"
print *, "EXTRACTING WAVEFUNCTION VALUES..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(neut_WF_array(lin_size, 0:2), ion_WF_array(lin_size, vib_num, num_targ))
call read_neut_WF(neut_WF_array)
call read_ion_WF(ion_WF_array)

print *, "----------------------------------------------------------------------------------------------------"
print *, "EXTRACTING DIPOLE DATA..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(dipole_data(3,num_channels,interval_num+1))
call read_dip_data(dipole_data, num_channels, quant_nums)

print *, "----------------------------------------------------------------------------------------------------"
print *, "DOCTORING DIPOLE DATA..."
print *, "----------------------------------------------------------------------------------------------------"

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
								if (abs(deriv_1) <= abs(deriv_2)) then
									cycle
								else
									do i = 1, vib_num
										dipole_data(coord_i, j_list(i), r_i) = -dipole_data(coord_i, j_list(i), r_i)
									end do
								end if

							else
								if (abs(deriv_1 - old_deriv) <= abs(deriv_2 - old_deriv)) then
									cycle
								else
									do i = 1, vib_num
										dipole_data(coord_i, j_list(i), r_i) = -dipole_data(coord_i, j_list(i), r_i)
									end do
								end if
							end if

							call find_cderiv(dipole_data(coord_i, j, r_i), dipole_data(coord_i, j, r_i-1), r_interval, old_deriv)

					else
						cycle
					end if
			end do
		end do
	else
		cycle
	end if
end do

call write_dips(dipole_data)

print *, "----------------------------------------------------------------------------------------------------"
print *, "APPROXIMATING INTEGRAL..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(dip_vect(3, num_channels, 0:2))

do neut_vib = 0, neut_vib_num -1
do j = 1, num_channels

		n = quant_nums(j,1)
		vib_i = quant_nums(j,2)

		x_comp_integrated = 0
		y_comp_integrated = 0
		z_comp_integrated = 0

		do r_i = 1, lin_size

			x_comp = dipole_data(1, j, r_i) * neut_WF_array(r_i, neut_vib) * ion_WF_array(r_i, vib_i+1, n) * lin_interval
			x_comp_integrated = x_comp_integrated + x_comp

			y_comp = dipole_data(2, j, r_i) * neut_WF_array(r_i, neut_vib) * ion_WF_array(r_i, vib_i+1, n) * lin_interval
			y_comp_integrated = y_comp_integrated + y_comp

			z_comp = dipole_data(3, j, r_i) * neut_WF_array(r_i, neut_vib) * ion_WF_array(r_i, vib_i+1, n) * lin_interval
			!if (n == 1 .AND. l == 0 .AND. m == 0) then
				!print *, j
				!print *, "r_i =", r_i
				!print *, "dip =", dipole_data(3, j, r_i)
				!print *, "z_comp =", z_comp
				!print *, "ion WF =", ion_WF_array(r_i, vib_i+1, n)
			!end if
			z_comp_integrated = z_comp_integrated + z_comp

		!sum over r do
		end do

		dip_vect(1,j,neut_vib) = x_comp_integrated
		dip_vect(2,j,neut_vib) = y_comp_integrated
		dip_vect(3,j,neut_vib) = z_comp_integrated

!sym do
end do
end do

open(2, file = "./integrated_dips_v0.dat")
do i = 1, num_channels
	write(2, *) i, dip_vect(1,i,0), dip_vect(2,i,0), dip_vect(3,i,0)
end do
close(2)

stop

print *, "----------------------------------------------------------------------------------------------------"
print *, "CHANNEL ELIMINATION..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(energy_array(num_channels),channel_energies(num_targ, 0:vib_num-1))
open(1, file = "./vibronic_channel_energies.txt")

print *, "Reading channel energies..."

do i = 0,(num_targ*vib_num)-1
	read(1,*) n, vib_i, channel_energies(n,vib_i)
	channel_energies(n,vib_i) = ev_au*channel_energies(n,vib_i)
end do
close(1)

energy_begin = channel_energies(1,0)
energy_end = channel_energies(num_targ, vib_num-1)
energy_step_num = 30000
energy_step = (energy_end - energy_begin)/energy_step_num

open(1, file = "./energies.dat")
do i = 1, num_channels

		n = quant_nums(i,1)
		vib_i = quant_nums(i,2)
		energy_array(i) = channel_energies(n,vib_i)

		write(1,*) energy_array(i)/ev_au
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

allocate(TCS_array(energy_step_num, 2, 0:2), SPES_array(energy_step_num, 2, 0:2), SPES_array_ground(energy_step_num, 2, 0:2))
SPES_array = 0
do neut_vib = 0, neut_vib_num - 1
do i = 1, energy_step_num

	E = energy_begin + i*energy_step

	num_open = 0
        Do a=1,num_channels
          If(energy_array(a)>E)Exit
        Enddo
        num_open=a-1

	if (num_open /= num_channels) then

		num_closed = num_channels - num_open
		allocate(dipole_open(3,num_open),dipole_closed(3,num_closed),IPIV(num_channels))
		allocate(smat_cc(num_closed,num_closed),smat_co(num_closed,num_open),beta(num_channels))
		dipole_closed = dip_vect(:,num_open+1:num_channels, neut_vib)
		dipole_open = dip_vect(:,1:num_open, neut_vib)
		smat_cc = Conjg(Transpose(smat(num_open+1:num_channels,num_open+1:num_channels)))
		smat_co = Conjg(Transpose(smat(1:num_open,num_open+1:num_channels)))

		beta = 0d0
		do a = 1, num_closed
			beta(a+num_open) = pi/sqrt(2*(energy_array(a+num_open)-E))
		end do

		allocate(dipole_phys(3,num_open))
		dipole_phys = 0d0

		do a = 1, num_closed
			smat_cc(a,a) = smat_cc(a,a) - exp(2d0*ci*beta(a+num_open))
		end do

		call ZGESV(num_closed,num_open,smat_cc,num_closed,IPIV(1:num_closed),smat_co,num_closed,INFO)


 		do a = 1,3
			dipole_phys(a,:) = dipole_open(a,:) - MatMul(dipole_closed(a,:),smat_co)
		end do
		deallocate(dipole_open, dipole_closed, smat_cc, smat_co, IPIV, beta)

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
	dipole_phys = dip_vect(:,:,neut_vib)

	end if

	omega = omega_0 + E
	!omega = omega_0 - channel_energies(0)

	TCS_array(i,1,neut_vib) = const*omega*(sum(abs(dipole_phys(1,:))**2)+sum(abs(dipole_phys(2,:))**2)+sum(abs(dipole_phys(3,:))**2))
	TCS_array(i,2,neut_vib) = E/ev_au + delta_Es(neut_vib)
	SPES_array(i,2,neut_vib) = TCS_array(i,2,neut_vib)

	do n = 1, num_targ
        do vib_i = 0, vib_num-1
            if ((abs(E - channel_energies(n,vib_i)) < elec_en_range) .AND. (E - channel_energies(n,vib_i) >= 0)) then
                do j = 1, num_channels
                    if (quant_nums(j,1) == n .AND. quant_nums(j,2) == vib_i) then
                        SPES_array(i,1,neut_vib) = SPES_array(i,1,neut_vib) + const*omega*(abs(dipole_phys(1,j))**2+abs(dipole_phys(2,j))**2+abs(dipole_phys(3,j))**2)
                    end if
                end do
            end if
        end do
    end do

	SPES_array_ground(i,2,neut_vib) = TCS_array(i,2,neut_vib)

    do vib_i = 0, vib_num-1
        if ((abs(E - channel_energies(1,vib_i)) < elec_en_range) .AND. (E - channel_energies(1,vib_i) >= 0)) then
            do j = 1, num_channels
                if (quant_nums(j,1) == 1 .AND. quant_nums(j,2) == vib_i) then
                    SPES_array_ground(i,1,neut_vib) = SPES_array_ground(i,1,neut_vib) + const*omega*(abs(dipole_phys(1,j))**2+abs(dipole_phys(2,j))**2+abs(dipole_phys(3,j))**2)
                end if
            end do
        end if
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

allocate(conv_TCS_array(num_conv_points,2,0:2,2), conv_SPES_array(num_conv_points,2,0:2,2), conv_SPES_array_ground(num_conv_points,2,0:2,2))
do neut_vib = 0, neut_vib_num - 1
    conv_size_array(1) = size(TCS_array,1)
    conv_size_array(2) = size(TCS_array,2)
    call cauchy_conv(conv_size_array, TCS_array(:,:,neut_vib), conv_TCS_array(:,:,neut_vib,1))
    call gauss_conv(conv_size_array, TCS_array(:,:,neut_vib), conv_TCS_array(:,:,neut_vib,2))
    call cauchy_conv(conv_size_array, SPES_array(:,:,neut_vib), conv_SPES_array(:,:,neut_vib,1))
    call gauss_conv(conv_size_array, SPES_array(:,:,neut_vib), conv_SPES_array(:,:,neut_vib,2))
    call cauchy_conv(conv_size_array, SPES_array_ground(:,:,neut_vib), conv_SPES_array_ground(:,:,neut_vib,1))
    call gauss_conv(conv_size_array, SPES_array_ground(:,:,neut_vib), conv_SPES_array_ground(:,:,neut_vib,2))
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

open(1, file = "./SPES_data_v_lorentz.dat")
do i = 1, num_conv_points
        write(1,*) conv_SPES_array_ground(i,2,0,1), conv_SPES_array_ground(i,1,0,1), conv_SPES_array_ground(i,1,1,1), conv_SPES_array_ground(i,1,2,1)
end do
close(1)

open(1, file = "./SPES_data_v_gauss.dat")
do i = 1, num_conv_points
        write(1,*) conv_SPES_array_ground(i,2,0,2), conv_SPES_array_ground(i,1,0,2), conv_SPES_array_ground(i,1,1,2), conv_SPES_array_ground(i,1,2,2)
end do
close(1)

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
