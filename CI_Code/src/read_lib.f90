module read_lib

	use constants
	
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
	
	subroutine read_neut_WF(neut_WF_array)

		use constants, only: interp_geoms, num_interp_points
	
		character (len = 65)				::	file_path
	
		integer						::	j, vib_i
		real(8)						::	r
		character (len = 1)         ::  vib_str
	
		real(8), intent(out), dimension(num_interp_points, 0:2)	::	neut_WF_array
	
	do j = 1, num_interp_points
		r = interp_geoms(j)
		do vib_i = 0, 2
			write(vib_str, "(i1)") vib_i
			file_path = "./Neutral_WF/wf00" // vib_str // "_neut.dat"
			open(vib_i+2, file = trim(file_path), status = "old")
			call WF_mag(vib_i+2, r, neut_WF_array(j,vib_i))
		end do
	end do
	
	end subroutine
	
	!make x


	subroutine read_ion_WF(ion_WF_array)

		use constants, only: interp_geoms, num_interp_points
	
		character (len = 286)				::	file_path
	
		integer						::	state, j, chan_i
		character (len = 1)			::	chan_str, state_str
		real(8)						::	r
	
		real(8), intent(out), dimension(num_interp_points, vib_num, num_targ)	::	ion_WF_array
	

	do state = 1, num_targ
		do chan_i = 0, 4
			do j = 1, num_interp_points

				r = interp_geoms(j)

				write(chan_str, "(i1)") chan_i
				write(state_str, "(i1)") state

				file_path = "./Ion_WF/Target_wf" // chan_str // "_CH+_ElecState" // state_str // ".dat"
				open(3, file = file_path, status = "old")
				call WF_mag(3, r, ion_WF_array(j,(chan_i+1),state))
			end do
		end do
	end do

	
	end subroutine
	
	subroutine WF_mag(file_num, r, WF)
	
		real(8)                  :: diff, var, r, res
		integer                  :: counter, file_num
		real(8), intent(out)     :: WF
		real(8), allocatable     :: WF_data(:)
	
	diff = 10
	counter = 0
	print *, file_num
	!	if (file_num == 2) then 
	!		allocate(WF_data(4000))
	!		read(file_num, *) WF_data
	!		res = WF_data(5)-WF_data(1
	
	!		do while (abs(diff) > res)
	!	    		var = WF_data(4*counter + 1)
	!	    		diff = var - r
	!	    		WF = WF_data(4*counter + 2)
	!	    		counter = counter + 1
	!		end do
	!	else
			allocate(WF_data(6000))
			read(file_num, *) WF_data
			res = WF_data(7)-WF_data(1)
	
			do while (abs(diff) > res)
					var = WF_data(6*counter + 1)
					diff = var - r
					WF = WF_data(6*counter + 2)
					counter = counter + 1
			end do
	 !   	end if
		close(file_num)
		deallocate(WF_data)
	end subroutine
	
	! subroutine read_dip_data(dipole_data, num_channels, quant_nums)
	
	!     integer, intent(in)                 ::  num_channels
	! 	complex(8), dimension(3, num_channels,interval_num+1), intent(out)	::	dipole_data
	
	! 	integer								::	i, sym_i, ang_i, ang_chan, r_i, block_size, channel_order
	! 	integer								::	n_josh, l_josh, m_josh, v, n, a, l, m, skip_num
	! 	real(8)								::	special_en
	! 	character (len = 3), dimension(0:interval_num)			::	geom_list
	! 	character(:), allocatable					::	file_name
	! 	character (len = 2)						    ::	sym_str
	! 	complex(8), dimension(3)					::	mom_array
	
	! 	real(8), dimension(20)						::	temp2
	! 	real(8), dimension(interval_num+1)			::	exc_levels
	! 	integer, dimension(num_channels,4), intent(out)				::	quant_nums
	
	
	! call get_geom_list(r_interval, r_min, interval_num, geom_list)
	
	! do r_i = 1, interval_num+1
	! 	open(1, file = "./Dipole_Data_CI/A1_" // geom_list(r_i-1) // ".formatted", status = "old")
	
	! 	do a = 1, 20
	! 		read(1,*)
	! 	end do
	
	! 	read(1,*) temp2(1:7), exc_levels(r_i)
	! 	close(1)
	! 	exc_levels(r_i) = .5*exc_levels(r_i)
	! end do
	
	! dipole_data = 0
	
	! do sym_i = 1,3
	
	!         if (sym_i == 1) then
	!         	ang_chan = 21
	!         else if (sym_i == 2 .OR. sym_i == 3) then
	!         	ang_chan = 19
	!         end if
	
	! 	sym_str = sym_strings(sym_i)
	
	! 	do ang_i = 1, ang_chan
	! 		do r_i = 1, interval_num+1
	
	! 			file_name = "./Dipole_Data_CI/" // sym_str // "_" // geom_list(r_i-1) // ".formatted"
	! 			open(1, file = file_name, status = "old")
	
	! 			do a = 1,2
	! 				read(1,*)
	! 			end do
	
	! 			read(1, *) temp2(1:5), block_size
	
	! 			do a = 1,(6+ang_i)
	! 				read(1,*)
	! 			end do
	
	! 			read(1,*) channel_order, n, temp2(1:3), l, m
	
	! 			!No triplet states, the QM excited state is the first excited singlet
	! 			if (n==2) then
	!                 n = 4
	!             else if (n==3) then
	!                 n = 5
	!             end if
	
	! 			print *, sym_strings(sym_i), channel_order, n, l, m
	! 			file_name = "./Dipole_Data_CI/" // sym_str // "_" // geom_list(r_i-1) // ".formatted"
	
	! 			if (sym_str == "A1") then
	! 				skip_num = 115
	! 			else
	! 				skip_num = 106
	! 			end if
	
	! 			if (geom_list(r_i-1) == "1.8") then
	! 				special_en = 0.03762036637931
	! 			else if (geom_list(r_i-1) == "2.0") then
	! 				special_en = 0.017323669738407
	!             else if (geom_list(r_i-1) == "2.2") then
	! 				special_en = 0.007323669738407
	!             else
	! 				special_en = extra_energy
	! 			end if
	
	! 			call dip_mag(ang_i, file_name, skip_num, special_en, block_size, exc_levels(r_i), mom_array)
	
	! 			if (l_max == 4) then
	! 				open(2, file = "./Smat/3states_lmax4/vibronic_VE_2states.channel", status = "old")
	! 			else if (l_max == 2) then
	! 				open(2, file = "./Smat/5state_lmax2/vibronic_DR.channel", status = "old")
	! 			end if
	! 			read(2, *)
	! 			read(2, *)
	! 			do i = 1, num_channels
	! 				read(2, *) n_josh, v, l_josh, m_josh
	! 				if (r_i == 1) then
	! 					quant_nums(i,:) = (/n_josh, v, l_josh, m_josh/)
	! 				end if
	! 				if (n_josh == 2 .OR. n_josh == 3) then
	!                     dipole_data(:,i,r_i) = 0
	! 				end if
	
	! 				if (n == 1) then
	! 					if (n == n_josh .AND. l == l_josh .AND.  m == m_josh) then
	! 						dipole_data(:,i,r_i) = mom_array
	! 					end if
	! 				else if (n==4 .OR. n==5) then
	! 					if ((n_josh /= 2 .OR. n_josh /=3)  .AND. l == l_josh .AND.  m == m_josh) then
	! 						dipole_data(:,i,r_i) = mom_array
	! 					end if
	! 				end if
	! 			end do
	! 			close(2)
	! 			close(1)
	! 		end do
	! 	end do
	! end do
	
	
	! 			open(50, file = "./dipole_data.dat")
	! 			do a = 1, num_channels
	! 				write(50, *) quant_nums(a,:)
	! 				write(50, *) dipole_data(:,a,1)
	! 			end do
	! 			close(50)
	

	! end subroutine
	
	subroutine new_dip_data(dipole_data, quant_nums)
	
		complex(8), intent(out)			::	dipole_data(:,:,:)
		integer, intent(out)			::	quant_nums(:,:)
	
		character(56)					::	file_string
		real(8)							::	x_real, x_imag, y_real, y_imag, z_real, z_imag, energy, r
		integer							::	sym_i, r_i, n_i, i, j, n, l, m, n_josh, l_josh, m_josh, num_rows, v
		complex(8)						::	x,y,z
		complex(8), dimension(3)		::	mom_array
		character(23)					::	temp

		character(4)					::	r_str
	
		do r_i = 1, interval_num +1
			
			r = r_min + (r_i-1)*r_interval
			write(r_str, "(F4.2)") r
			print *, "R str: ", r_str

			do sym_i = 1, 3
	
				file_string = "./New_Rmat_data/R-" // r_str // "/photoionization-amplitudes-" // sym_strings(sym_i) // ".dat"
				print *, file_string
				!call FindStr(file_string, 2, "  Energy=   5.0000000000000003E-002")
				Open(2 ,File=file_string, status = "old")
				read(2,*)
				read(2,*)
				read(2,*)
				read(2,*)
			

	
				do n_i = 1, num_targ

					if (n_i /= 1) read(2,*)
					read(2,"(A23,I1)") temp, n
	
					if (sym_i == 1) then
						if (n == 1) then
							num_rows = 12
						else
							num_rows = 10
						end if
					else if (sym_i == 2) then
						if (n == 1) then
							num_rows = 5
						else
							num_rows = 10
						end if
					else if (sym_i == 3) then
						if (n == 1) then
							num_rows = 10
						else if (n == 2) then
							num_rows = 5
						else if (n == 3) then
							num_rows = 12
						end if
					end if
	
					do j = 1, num_rows
						read(2,*) energy, l, m, x_real, x_imag, y_real, y_imag, z_real, z_imag 
						x = complex(x_real,x_imag)
						y = complex(y_real,y_imag)
						z = complex(z_real,z_imag)
						mom_array = (/x,y,z/)
	
						if (l_max == 4) then
							open(3, file = "./Smat/3states_lmax4/vibronic_DR.channel", status = "old")
						else if (l_max == 2) then
							open(3, file = "./Smat/5state_lmax2/vibronic_DR.channel", status = "old")
						end if
						read(3, *)
						read(3, *)
						do i = 1, num_channels
							read(3, *) n_josh, v, l_josh, m_josh
							if (r_i == 1) then
								quant_nums(i,:) = (/n_josh, v, l_josh, m_josh/)
							end if
							if (n == n_josh .AND. l == l_josh .AND.  m == m_josh) then
								dipole_data(:,i,r_i) = mom_array
							end if
						end do
						close(3)
					end do
				end do
				close(2)
			end do

			if (r_i == 12) then
				open(333, file = "rewrite_new_data.dat")
				do i = 1, num_channels
					write(333,*) quant_nums(i, 1), quant_nums(i, 2), quant_nums(i, 3), quant_nums(i, 4), dipole_data(1,i,r_i), dipole_data(2,i,r_i), dipole_data(3,i,r_i)
				end do
				close(333)
			end if


		end do
	end subroutine new_dip_data
	
	subroutine dip_mag(ang_i, file_name, skip_num, special_en, block_size, exc_en, dip_mom)
	
		real(8)							::	res, step, en_1, en_2, cur_energy, des_energy, exc_en
		real(8)							::	x_real, x_imag, z_real, z_imag, y_real, y_imag
		complex(8)						::	x, y, z
		integer							::	ang_i, block_size, a
		real(8), dimension(0:10000000)				::	temp1
	
		real(8), intent(in)					::	special_en
		integer, intent(in)					::	skip_num
		character(:), allocatable, intent(in)				::	file_name
	
		complex(8), dimension(3), intent(out)			::	dip_mom
		real(8), dimension(100)					::	temp
	
	open(3, file = trim(file_name), status = "old")
		do a = 1, skip_num
			read(3,*)
		end do
	read(3, *) temp1(0:(block_size*8)*(ang_i-1)), en_1
	read(3, *) temp(1), en_2
	step = en_2 - en_1
	res = 10
	
	des_energy = special_en + exc_en
	
	print *, des_energy
	do while (abs(res) > step)
		read(3, *) temp(1), cur_energy, y_real, y_imag, z_real, z_imag, x_real, x_imag
		res = des_energy - cur_energy
	end do
	close(3)
	
	!print *, x_real, x_imag, y_real, y_imag, z_real, z_imag
	
	x = cmplx(x_real, x_imag, 8)
	y = cmplx(y_real, y_imag, 8)
	z = cmplx(z_real, z_imag, 8)
	
	dip_mom(1) = x
	dip_mom(2) = y
	dip_mom(3) = z
	
	end subroutine dip_mag
	
	Subroutine FindStr(Name_of_in,num_file,str_to_find)
		Integer num_file,ifile
		Character*35 str_to_find, blockname
		Character*56 Name_of_in
		Open(num_file,File=Name_of_in, status = "old")
		do ifile = 1,1000
		  read(num_file,'(A35)') blockname
		  print *, blockname
		  If (blockname.eq.str_to_find) Return
		End do
		Close(num_file)
		Print *,'The string ',str_to_find,' is not found in the file ',Name_of_in
	  End Subroutine FindStr
	
	end module read_lib
	