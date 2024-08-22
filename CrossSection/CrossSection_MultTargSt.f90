program dip_moment

use constants
use read_lib
use math_ops

    implicit none

    complex(8),  dimension(3,5,9,3)       ::   dip_vect
    complex(8),  dimension(125,3)	  ::   arr_dip_vect
    real(8),     dimension(125)           ::   energy_array
    integer,     dimension(4)             ::   par_list
    real(8),     dimension(0:4)           ::   channel_energies

    real(8)      ::   r_max, r_min, r_interval, r, energy_begin, energy_end, energy_step
    real(8)      ::   E, const, omega
    integer      ::   energy_step_num, degeneracy, l_max, a, i, energy_i, sym_i, ang_i, num_open
    integer	 ::   ang_chan, vib_num, vib_i, r_i, interval_num, l, m, INFO
    integer      ::   num_channels, num_closed, channel_indicator
    complex(8)   ::   x_comp, x_comp_integrated, y_comp, y_comp_integrated, z_comp, z_comp_integrated
    complex(8)	 ::   deriv_1, deriv_2, old_deriv

    character (len = 256)		::   smat_file
    character (len = 3)			::   r_str
    character (len = 2)			::   sym_str
    character (len = 256)		::   cwd

    complex(8), allocatable		::	dipole_data(:,:,:,:), smat(:,:)
    complex(8), allocatable		::	dipole_phys(:,:), dipole_open(:,:), dipole_closed(:,:), smat_cc(:,:), smat_co(:,:)
    complex(8), allocatable		::	beta(:), part_vect(:)
    integer,    allocatable		::	IPIV(:)
    real(8),    allocatable		::	cross_section(:,:), neut_WF_array(:), ion_WF_array(:,:)

call getcwd(cwd)

vib_num = 5
interval_num = 6

const = (4./3.)*(pi**2)*alpha

arr_dip_vect = 0

l_max = 4
energy_i = 50

r_max = 3.
r_min = 1.8
r_interval = (r_max - r_min)/interval_num

degeneracy = (l_max + 1)**2
num_channels = degeneracy*vib_num

print *, "----------------------------------------------------------------------------------------------------"
print *, "READING S-MATRIX IN..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(smat(num_channels, num_channels))
smat_file = trim(cwd) // "/Smat/vibronic_DR.Smat"
call read_smat(smat_file, num_channels, smat)

print *, "----------------------------------------------------------------------------------------------------"
print *, "EXTRACTING WAVEFUNCTION VALUES..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(neut_WF_array(interval_num), ion_WF_array(interval_num, vib_num))
call read_neut_WF(r_interval, interval_num, neut_WF_array)
call read_ion_WF(r_interval, interval_num, vib_num, ion_WF_array)

print *, "----------------------------------------------------------------------------------------------------"
print *, "EXTRACTING DIPOLE DATA..."
print *, "----------------------------------------------------------------------------------------------------"

allocate(dipole_data(3,3,9,interval_num+1))
call read_dip_data(interval_num, r_min, r_interval, energy_i, dipole_data)

print *, "----------------------------------------------------------------------------------------------------"
print *, "DOCTORING DIPOLE DATA..."
print *, "----------------------------------------------------------------------------------------------------"


do sym_i = 1,3

	if (sym_i == 1) then
	ang_chan = 9
	else if (sym_i == 2 .OR. sym_i == 3) then
	ang_chan = 6
	end if	
	
	do ang_i = 1, ang_chan
		do r_i = 1, interval_num + 1
			do i = 1,3
				r = r_min + r_interval*(r_i - 1)
				write(r_str, "(F3.1)") r

				if (r_i > 1) then
					
					print *, r_str, "sym_i:", sym_i, "coord:", i	
					print *, real(dipole_data(i,sym_i,ang_i,r_i)), aimag(dipole_data(i,sym_i,ang_i,r_i))

					call find_cderiv(dipole_data(i, sym_i, ang_i, r_i), dipole_data(i, sym_i, ang_i, r_i-1), r_interval, deriv_1)
					call find_cderiv(-dipole_data(i, sym_i, ang_i, r_i), dipole_data(i, sym_i, ang_i, r_i-1), r_interval, deriv_2)

					if (r_i == 2) then
						if (abs(deriv_1) <= abs(deriv_2)) then
							cycle
						else 
							dipole_data(i, sym_i, ang_i, r_i) = -dipole_data(i, sym_i, ang_i, r_i)
						end if 
						
					else 
						if (abs(deriv_1 - old_deriv) <= abs(deriv_2 - old_deriv)) then
							cycle
						else 
							dipole_data(i, sym_i, ang_i, r_i) = -dipole_data(i, sym_i, ang_i, r_i)
						end if
					end if

					call find_cderiv(dipole_data(i, sym_i, ang_i, r_i), dipole_data(i, sym_i, ang_i, r_i-1), r_interval, old_deriv)

				else
					cycle
				end if
			end do	
		end do
	end do
end do

print *, "----------------------------------------------------------------------------------------------------"
print *, "APPROXIMATING INTEGRAL..."
print *, "----------------------------------------------------------------------------------------------------"

do sym_i = 1, 3

	print *, "Beginning ", sym_strings(sym_i), " symmetry..."

	if (sym_i == 1) then
	ang_chan = 9
	else if (sym_i == 2 .OR. sym_i == 3) then
	ang_chan = 6
	end if	

	do ang_i = 1, ang_chan 

		print *, "Beginning ", ang_i, " angular channel..."
			
		    do vib_i = 1, vib_num

			x_comp_integrated = 0
			y_comp_integrated = 0
			z_comp_integrated = 0

			do r_i = 1, interval_num+1
 
				x_comp = dipole_data(1, sym_i, ang_i, r_i) * neut_WF_array(r_i) * ion_WF_array(r_i, vib_i) * r_interval
 				x_comp_integrated = x_comp_integrated + x_comp

				y_comp = dipole_data(2, sym_i, ang_i, r_i) * neut_WF_array(r_i) * ion_WF_array(r_i, vib_i) * r_interval
 				y_comp_integrated = y_comp_integrated + y_comp

				z_comp = dipole_data(3, sym_i, ang_i, r_i) * neut_WF_array(r_i) * ion_WF_array(r_i, vib_i) * r_interval
 				z_comp_integrated = z_comp_integrated + z_comp

			!sum over r do
			end do

			dip_vect(1, vib_i, ang_i, sym_i) = x_comp_integrated
			dip_vect(2, vib_i, ang_i, sym_i) = y_comp_integrated
			dip_vect(3, vib_i, ang_i, sym_i) = z_comp_integrated

			print *, "v = ", vib_i, "ang = ", ang_i, "sym = ", sym_i

		!channels do
		end do 
	 	
	!angular do
	end do	   

!sym do
end do

open(2, file = trim(cwd) // "/dipolevsR.dat")
do r_i = 1, interval_num+1

	ang_i = 3
	sym_str = sym_strings(2)
	r = r_min + r_interval*(r_i - 1)
	write(r_str, "(F3.1)") r
	write(2, *) r, real(dipole_data(2,2,ang_i,r_i)), aimag(dipole_data(2,2,ang_i,r_i))

end do
close(2)

print *, "----------------------------------------------------------------------------------------------------"
print *, "REARRANGING DIPOLE MATRICES..."
print *, "----------------------------------------------------------------------------------------------------"

open(1, file = trim(cwd) // "/Smat/vibronic_DR.channel")

read(1,"(I2,A)") a
read(1,*)

do i = 1, num_channels
	read(1,*) par_list(1:4)

	vib_i = par_list(2)+1
	l = par_list(3)
	m = par_list(4)

	call get_index(l, m, sym_i, ang_i)

	if (sym_i /= 0) then
		arr_dip_vect(i, 1) = dip_vect(1, vib_i, ang_i, sym_i) 
		arr_dip_vect(i, 2) = dip_vect(2, vib_i, ang_i, sym_i)
		arr_dip_vect(i, 3) = dip_vect(3, vib_i, ang_i, sym_i)

	else 
		arr_dip_vect(i, 1:3) = 0

	end if
end do

print *, "Arrange Complete."
close(1)

print *, "----------------------------------------------------------------------------------------------------"
print *, "CHANNEL ELIMINATION..."
print *, "----------------------------------------------------------------------------------------------------"

open(1, file = trim(cwd) // "/vibronic_channel_energies.txt")

print *, "Reading channel energies..."

do i = 0,4
	read(1,*) channel_energies(i)
	channel_energies(i) = ev_au*channel_energies(i)
end do
close(1)

energy_begin = channel_energies(0)
energy_end = channel_energies(4)
energy_step_num = 50000
allocate(cross_section(energy_step_num,3))
energy_step = (energy_end - energy_begin)/energy_step_num

open(1, file = trim(cwd) // "/energies.dat")
do i = 1, num_channels
		
		channel_indicator = floor((i-1)/real(degeneracy))
		energy_array(i) = channel_energies(channel_indicator)

		write(1,*) energy_array(i)
end do
close(1)

do i = 1, energy_step_num

	E = energy_begin + i*energy_step
	
	num_open = 0
        Do a=1,num_channels
          If(energy_array(a)>E)Exit
        Enddo
        num_open=a-1

	if (num_open /= num_channels) then
		
		num_closed = num_channels - num_open
		allocate(dipole_open(num_open,3),dipole_closed(num_closed,3),IPIV(num_channels))
		allocate(smat_cc(num_closed,num_closed),smat_co(num_closed,num_open),beta(num_channels))
		dipole_closed = arr_dip_vect(num_open+1:num_channels,:)
		dipole_open = arr_dip_vect(1:num_open,:)
		smat_cc = Conjg(Transpose(smat(num_open+1:num_channels,num_open+1:num_channels)))
		smat_co = Conjg(Transpose(smat(1:num_open,num_open+1:num_channels)))
	
		beta = 0.
		do a = 1, num_closed
			beta(a+num_open) = pi/sqrt(2.*(energy_array(a+num_open)-E)) 
		end do

		allocate(dipole_phys(num_open,3))
		dipole_phys = 0.		
			
		do a = 1, num_closed
			smat_cc(a,a) = smat_cc(a,a) - exp(2.*ci*beta(a+num_open))
		end do
	
		call ZGESV(num_closed,num_open,smat_cc,num_closed,IPIV(1:num_closed),smat_co,num_closed,INFO)

 		do a = 1,3
			dipole_phys(:,a) = dipole_open(:,a) - MatMul(dipole_closed(:,a),smat_co)
		end do
		deallocate(dipole_open, dipole_closed, smat_cc, smat_co, IPIV, beta)
	else

	num_open = num_channels
	allocate(dipole_phys(num_open, 3))
	dipole_phys = arr_dip_vect

	end if

	omega = omega_0 + E - channel_energies(0)

	do a = 1,3
		allocate(part_vect(num_open))
		part_vect = dipole_phys(:,a)
		cross_section(i,a) = const*omega*sum(abs(part_vect)**2)
		deallocate(part_vect)
	end do

	deallocate(dipole_phys) 
 
end do
close(1)

open(1, file = trim(cwd) // "/cross_section.dat")
do i = 1, energy_step_num

	write(1,*) energy_begin + i*energy_step, cross_section(i,1), cross_section(i,2), cross_section(i,3) 

end do
close(1)

contains

subroutine get_index(l, m, sym, ind)

   integer                           ::   l, m
   integer, intent(out)              ::   sym, ind
   integer, dimension(3,2,9)         ::   quant_num_array = 0
   integer, dimension(0:4,-4:4,2)    ::   rev_num_array = 0

quant_num_array(1,1,1:9) = (/0,1,2,2,3,3,4,4,4/)
quant_num_array(1,2,1:9) = (/0,0,0,2,0,2,0,4,2/)

quant_num_array(2,1,1:9) = (/1,2,3,3,4,4,0,0,0/)
quant_num_array(2,2,1:9) = (/1,1,1,3,1,3,0,0,0/)

quant_num_array(3,1,1:9) = (/1,2,3,3,4,4,0,0,0/)
quant_num_array(3,2,1:9) = (/-1,-1,-1,-3,-3,-1,0,0,0/)

rev_num_array(0:4,0,1) = 1
rev_num_array(0:4,0,2) = (/1,2,3,5,7/)

rev_num_array(0:4,2,1) = 1
rev_num_array(0:4,2,2) = (/0,0,4,6,9/)

rev_num_array(0:4,4,1) = 1
rev_num_array(0:4,4,2) = (/0,0,0,0,8/)

rev_num_array(0:4,1,1) = 2
rev_num_array(0:4,1,2) = (/0,1,2,3,5/)

rev_num_array(0:4,3,1) = 2
rev_num_array(0:4,3,2) = (/0,0,0,4,6/)

rev_num_array(0:4,-1,1) = 3
rev_num_array(0:4,-1,2) = (/0,1,2,3,6/)

rev_num_array(0:4,-3,1) = 3
rev_num_array(0:4,-3,2) = (/0,0,0,4,5/)

sym = rev_num_array(l, m, 1)
ind = rev_num_array(l, m, 2)

end subroutine

end program
