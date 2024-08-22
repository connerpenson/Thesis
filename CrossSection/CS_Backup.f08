program dip_moment

    implicit none

    complex(8),  dimension(3)             ::   mom_array
    !complex(8), dimension(501,3,5,9,3)   ::   write_mat
    complex(8),  dimension(3,5,9,3)       ::   dip_vect
    complex(8),  dimension(125,3)	  ::   arr_dip_vect
    real(8),     dimension(500)           ::   energy_array
    real(8),     dimension(4)             ::   neut_WF_arr
    real(8),     dimension(4,5)           ::   ion_WF_arr
    real(8),     dimension(2)             ::   ex_array
    complex(8),  dimension(125,125)       ::   s_mat
    real(8),     dimension(2, 125)        ::   tmp
    !real(8),     dimension(3)             ::   dip_mags
    integer,     dimension(4)             ::   par_list

    real(8)      ::   r
    integer      ::   a, i, count_sum, energy_i, j, chan_i, sym_i, ang_i, ang_chan, chan_num, comp_i, vib_i, r_i, interval_num, l, m
    complex(8)   ::   data_val, x_comp, x_comp_integrated, y_comp, y_comp_integrated, z_comp, z_comp_integrated, ci

    character (len = 59)               ::   WF_file_neut, Wf_file_ion, file_str
    character (len = 3)                ::   r_str
    character (len = 2)                ::   sym_str
    character (len = 1)                ::   chan_str, ang_str, comp_str
    character (len = 2), dimension(3)  ::   sym_strings, comp_strings
    character (len = 256)              ::   cwd


call getcwd(cwd)

ex_array = shape(ion_WF_arr)
chan_num = ex_array(2)
interval_num = ex_array(1)

arr_dip_vect(:,:) = (0,0)

ci = (0d0,1d0)

sym_strings = (/"A1", "B1", "B2"/)
comp_strings = (/"x", "y", "z"/)


print *, "----------------------------------------------------------------------------------------------------"
print *, "READING S-MATRIX IN..."
print *, "----------------------------------------------------------------------------------------------------"

open(1, file = trim(cwd) // "/Smat/vibronic_DR.Smat")	
do i = 1,125
	read(1,*) tmp(1:2, i)
	s_mat(1:i, i) = tmp(1, 1:i) + tmp(2, 1:i)*ci
	s_mat(1:i, i) = s_mat(i, 1:i)
end do
close(1)

!do energy_i = 0, 499
!	energy_array(energy_i+1) = .000735 + .000735*energy_i 
!end do

energy_i = 50 !selecting 1ev over ionization

print *, "----------------------------------------------------------------------------------------------------"
print *, "EXTRACTING WAVEFUNCTION VALUES..."
print *, "----------------------------------------------------------------------------------------------------"

do j = 0, 3

	r = 1.8 + .4*j
		
	write(r_str, "(F3.1)") r

	WF_file_neut = trim(trim(cwd) // "/wf000_neut.dat")
	open(2, file = WF_file_neut, status = "old")
	call WF_mag(2, r, neut_WF_arr(j+1))

end do

print *, "Neutral Wavefunction:"
print *, "r = 1.8, ", neut_WF_arr(1)
print *, "r = 2.2, ", neut_WF_arr(2)
print *, "r = 2.6, ", neut_WF_arr(3)
print *, "r = 3.0, ", neut_WF_arr(4)


do chan_i = 0, 4

	do j = 0,3

		r = 1.8 + .4*j

		write(chan_str, "(i1)") chan_i

		WF_file_ion = trim(trim(cwd) // "/wf00" // chan_str // "_ion.dat")
		open(3, file = WF_file_ion, status = "old")
		call WF_mag(3, r, ion_WF_arr((j+1),(chan_i+1)))
	
	end do

	print *, "Ionic Wavefunction:", chan_i
	print *, "r = 1.8, ", ion_WF_arr(1,chan_i+1)
	print *, "r = 2.2, ", ion_WF_arr(2,chan_i+1)
	print *, "r = 2.6, ", ion_WF_arr(3,chan_i+1)
	print *, "r = 3.0, ", ion_WF_arr(4,chan_i+1)
	
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

		!do energy_i = 0, 499
		energy_i = 50
		
		    sym_str = trim(sym_strings(sym_i))
		    open(1, file = trim(cwd) // "/DipoleData/" // sym_str // "_" // r_str, status = "old")
		    call dip_mag(energy_i, ang_i, mom_array)
			
		    do vib_i = 1, chan_num
	
			x_comp_integrated = 0
			y_comp_integrated = 0
			z_comp_integrated = 0

			do r_i = 1, interval_num
 
				x_comp = mom_array(1) * neut_WF_arr(r_i) * ion_WF_arr(r_i, vib_i) * .4
 				x_comp_integrated = x_comp_integrated + x_comp

				y_comp = mom_array(2) * neut_WF_arr(r_i) * ion_WF_arr(r_i, vib_i) * .4
 				y_comp_integrated = y_comp_integrated + y_comp

				z_comp = mom_array(3) * neut_WF_arr(r_i) * ion_WF_arr(r_i, vib_i) * .4
 				z_comp_integrated = z_comp_integrated + z_comp

			!sum over r do
			end do

		   ! write_mat(energy_i + 1, 1, vib_i, ang_i, sym_i) = x_comp_integrated
		  !  write_mat(energy_i + 1, 2, vib_i, ang_i, sym_i) = y_comp_integrated
		 !   write_mat(energy_i + 1, 3, vib_i, ang_i, sym_i) = z_comp_integrated
			
		     dip_vect(1, vib_i, ang_i, sym_i) = x_comp_integrated
		     dip_vect(2, vib_i, ang_i, sym_i) = y_comp_integrated
		     dip_vect(3, vib_i, ang_i, sym_i) = z_comp_integrated

		     print *, "v = ", vib_i, "ang = ", ang_i, "sym = ", sym_i

		    !channels do
		    end do 
	 	
		!energy do
		!end do
	
	!angular do
	end do	   

!sym do
end do

!open(1, file = trim(cwd) // "/dipole_mags.txt", status = "old")
!do sym_i = 1, 3
!
!	if (sym_i == 1) then
!	ang_chan = 9
!	else if (sym_i == 2 .OR. sym_i == 3) then
!	ang_chan = 6
!	end if	
!
!	do ang_i = 1, ang_chan 
!   			
!		    do sum_i = 1, chan_num
!
!			dip_mags(1) = sqrt(real(dip_vect(1, sum_i, ang_i, sym_i))**2 + aimag(dip_vect(1, sum_i, ang_i, sym_i))**2)
!			dip_mags(2) = sqrt(real(dip_vect(2, sum_i, ang_i, sym_i))**2 + aimag(dip_vect(2, sum_i, ang_i, sym_i))**2)
!			dip_mags(3) = sqrt(real(dip_vect(3, sum_i, ang_i, sym_i))**2 + aimag(dip_vect(3, sum_i, ang_i, sym_i))**2)
!
!		
!			write(1, *) "L = ", quant_num_array(sym_i, 1, ang_i), "m = ", quant_num_array(sym_i, 1, ang_i), "Vib = ", sum_i
!			write(1, *) dip_mags(1), dip_mags(2), dip_mags(3)
!			write(1, *) "---------------------------------------------------------------------"			
!
!		     end do
!	end do
!end do

!close(1) 			
	
!print *, "----------------------------------------------------------------------------------------------------"
!print *, "WRITING INTO DATA FILES..."
!print *, "----------------------------------------------------------------------------------------------------"

!do chan_i = 0, 4
!
!	do sym_i = 1,3
!				
!		if (sym_i == 1) then
!		ang_chan = 9
!		else if (sym_i == 2 .OR. sym_i == 3) then
!		ang_chan = 6
!		end if	
!	
!		do ang_i = 1, ang_chan
!			
!			do comp_i = 1, 3
!		
!				write(chan_str, "(i1)") chan_i
!				write(ang_str, "(i1)") ang_i
!				sym_str = trim(sym_strings(sym_i))
!				comp_str = trim(comp_strings(comp_i))
!							
!				file_str = trim(cwd) // "/N0_I" // chan_str // "/" // sym_str // "/" // ang_str // "/" // comp_str // ".dat"
!				open(1, file = file_str, status = "old")
!				
!				do energy_i = 1,500
!					
!					data_val = write_mat(energy_i, comp_i, chan_i+1, ang_i, sym_i)
!					write(1, *) energy_array(energy_i), real(data_val), aimag(data_val)
!				
!				!energy	
!				end do
!					close(1)
!
			!components
!			end do
		!angulular
!		end do
	!symmetry
!	print *, sym_str, " symmetry complete."
!	end do
!channel
!print *, chan_i, " channel complete."
!end do

print *, "----------------------------------------------------------------------------------------------------"
print *, "REARRANGING DIPOLE MATRICES..."
print *, "----------------------------------------------------------------------------------------------------"

open(1, file = trim(cwd) // "/Smat/vibronic_DR.channel")

read(1,"(I2,A)") a
read(1,*)

do i = 1,125
	read(1,*) par_list(1:4)

	vib_i = par_list(2)
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

close(1)


contains

subroutine dip_mag(energy_i, ang_i, dip_mom)

    real(8)                                              :: x_real, x_imag, z_real, z_imag, y_real, y_imag
    complex(8)                                          :: x, y, z
    integer, intent(in)                                  :: energy_i, ang_i
    real(8), dimension(2  + 8*energy_i + 4000*(ang_i-1)) :: temp1

    !dip_mom now an array with dip_mom(1,2,3) being x y and z comp
    complex(8), dimension(3), intent(out)               :: dip_mom

    read (1, *) temp1, x_real, x_imag, z_real, z_imag, y_real, y_imag
    close(1)

    !print *, x_real, x_imag, z_real, z_imag, y_real, y_imag

    x = cmplx(x_real, x_imag)
    y = cmplx(y_real, y_imag)
    z = cmplx(z_real, z_imag)

    dip_mom(1) = x
    dip_mom(2) = y
    dip_mom(3) = z

    return

end subroutine dip_mag

subroutine WF_mag(file_num, r, WF)

    !Routine pulls wf value from specified file for the given r value
    real(8)                  :: diff, var, r, res
    integer                  :: counter, file_num
    real(8), intent(out)     :: WF
    real(8), dimension(4000) :: WF_data 

    read(file_num, *) WF_data
    diff = 10
    counter = 0
    
    res = WF_data(5)-WF_data(1) 

    do while (abs(diff) > res)
        var = WF_data(4*counter + 1)
        diff = var - r
        WF = WF_data(4*counter + 2)
	counter = counter + 1
    end do
    close(file_num)

end subroutine

subroutine get_index(l, m, sym, ind)

   integer                           ::   l, m
   integer, intent(out)              ::   sym, ind
   integer, dimension(3,2,9)         ::   quant_num_array = 0
   integer, dimension(0:4,-3:4,2)    ::   rev_num_array = 0

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
