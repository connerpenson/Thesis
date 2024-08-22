program EI_Exc_Calc

    use global_params
    use stdlib_sorting, only: sort_index
    use iso_fortran_env, only: int64, iostat_end

    implicit none

    integer ::  i, j, n, v, l, m, n_prime, v_prime, l_prime, m_prime, i_Q, elec_i, vib_i, counter, ie, ei, ei_p, vi, vi_p, i_Q_
    real(8) ::  E_tot, E_elec, E, energy_begin, energy_end, energy_step

    real(8), dimension(dim_lm, dim_lm, num_Q)      :: K_Q, Delta_Q, fixed_Delta_Q, trimmed_Delta_Q
    real(8), dimension(dim_lm, dim_lm, num_Q_)     :: lin_Delta_Q, final_Delta_Q
    real(8), dimension(dim_lm, num_Q)      :: eigenphases
    complex*16, dimension(dim_lm, dim_lm, num_Q_)   :: S_Q
    integer, dimension(100,num_gamma_Q,3,num_Q)    :: D2h_Dinf
    integer, dimension(num_Q, num_gamma_Q)         :: md
    real(8), dimension(0:vib_num-1, num_Q_, n_max)      :: ion_wf
    complex*16, dimension(num_tot, num_tot)        :: S_ft
    complex*16, dimension(num_tot, num_tot)        :: replaced_S_ft
    complex*16, dimension(num_tot, num_tot)        :: sorted_S_ft
    integer, dimension(num_tot, 4)                 :: ft_channels, sorted_ft_channels
    real(8), dimension(n_max, n_max, 0:vib_num, 0:vib_num)               :: eics_matrix
    !real(8), dimension(1, 2:2, 0:vib_num, 0:vib_num)  :: diff_mat
    real(8), dimension(num_scatE)                  :: energies
    real(8), dimension(num_tot)                    :: energy_array
    integer(kind=int64), dimension(num_tot) ::  index
    integer, dimension(energy_step_num) ::  max_chans

    good_elem_list(1,:) = (/1,1,0,1,1,0/) 
    good_elem_list(2,:) = (/1,1,0,4,0,0/)
    good_elem_list(3,:) = (/2,2,1,2,2,1/)
    good_elem_list(4,:) = (/3,2,-1,3,2,-1/)
    good_elem_list(5,:) = (/4,0,0,1,1,0/)

    call get_geoms
    call get_lin_geoms
    call get_channel_ens
    call make_electron_energies

    call get_KQ(K_Q, md, D2h_Dinf, energies)

    open(333, file = "kmat.dat")
    counter = 0
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        if (sum(abs(K_Q(i,j,:))) > 1e-3) then
          counter = counter +1
          write(333, *) "#Partial Wave:", counter, n, l, m, n_prime, l_prime, m_prime
          do i_Q = 1, num_Q
            write(333, *) geoms(i_Q), K_Q(i,j,i_Q)
          end do
          write(333, *)
        end if
      end do
    end do
    close(333)

    do i_Q = 1, num_Q
      call K_delta(K_Q(:,:,i_Q), Delta_Q(:,:,i_Q), eigenphases(:,i_Q))
    end do

    call sign_fixer(Delta_Q, fixed_Delta_Q)
    call trim_Delta(fixed_Delta_Q, trimmed_Delta_Q)

    do i = 1, dim_lm
      do j = 1, dim_lm

        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)

        if ( ( (n == 2 .and. l == 2 .and. m == 1) .and. (n == 2 .and. l == 2 .and. m == 1) ) .or. &
        & ( (n_prime == 3 .and. l_prime == 2 .and. m_prime == -1) .and. (n_prime == 3 .and. l_prime == 2 .and. m_prime == -1) ) ) then

          fixed_Delta_Q(i,j,num_Q) = fixed_Delta_Q(i,j,num_Q-1)

        end if
      end do
    end do

    call linearize_Delta(trimmed_Delta_Q(:,:,:), lin_Delta_Q(:,:,:))
    call add_nice_elems(fixed_Delta_Q, lin_Delta_Q, final_Delta_Q)

    open(333, file = "delta_interp_points.dat")
    counter = 0
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        if (sum(abs(Delta_Q(i,j,:))) > 1e-3) then
          counter = counter +1
          write(333, *) "#Partial Wave:", counter, n, l, m, n_prime, l_prime, m_prime
          write(333, *) geoms(lin_start), Delta_Q(i,j,lin_start)
          write(333, *) geoms(lin_end), Delta_Q(i,j,lin_end)
          write(333, *)
        end if
      end do
    end do
    close(333)

    open(333, file = "lin_delta.dat")
    counter = 0
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        if (sum(abs(final_Delta_Q(i,j,:))) > 0) then
          counter = counter +1
          write(333, *) "#Partial Wave:", counter, n, l, m, n_prime, l_prime, m_prime
          do i_Q_ = 1, num_Q_
            write(333, *) lin_geoms(i_Q_), final_Delta_Q(i,j,i_Q_)
          end do
          write(333, *)
        end if
      end do
    end do
    close(333)

    open(333, file = "eigenphases.dat")
    counter = 0
    do i=1, dim_lm
      call i2lm(i, n, l, m)
      counter = counter +1
      write(333, *) "#Partial Wave:", counter, n, l, m
      do i_Q = 1, num_Q
        write(333, *) geoms(i_Q), eigenphases(i,i_Q)
      end do
      write(333, *)
    end do
    close(333)

    
    open(333, file = "fixed_delta.dat")
    counter = 0
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        if (sum(abs(fixed_Delta_Q(i,j,:))) > 5) then
          counter = counter +1
          write(333, *) "#Partial Wave:", counter, n, l, m, n_prime, l_prime, m_prime
          do i_Q = 1, num_Q
            write(333, *) geoms(i_Q), fixed_Delta_Q(i,j,i_Q)
          end do
          write(333, *)
        end if
      end do
    end do
    close(333)

    !call check_unitarity(fixed_S_Q(:,:,:,i_energy))

    open(27, file = "check_like_elements.dat")
    do i = 1, dim_lm
      call i2lm(i, n, l, m)
      do j = 1, dim_lm
        call i2lm(j, n_prime, l_prime, m_prime)

        if ((n == 1 .and. l == 4 .and. (m == 4 .or. m == -4)) .and. (n_prime == 2 .or. n_prime == 3 )) then
          if (K_Q(i,j,4) /= 0) then
            write(27, *) m, n_prime, l_prime, m_prime, K_Q(i,j,4)
          end if
        end if

      end do
    end do
    close(27)

    call get_wf(ion_wf)
            
    open(27, file = "WF_vals.dat")
    do vib_i = 0, vib_num-1
        do elec_i = 1, n_max
            write(27,*) "#e=",elec_i,"v=", vib_i
            do i_Q_ = 1, num_Q_
                write(27,*) lin_geoms(i_Q_), ion_wf(vib_i, i_Q_, elec_i)
            end do
            write(27,*)
        end do
        write(27,*)
    end do
    close(27)

    do i_Q_ = 1, num_Q_
      call Delta_S(final_Delta_Q(:,:,i_Q_), S_Q(:,:,i_Q_))
    end do
 
    open(333, file = "smat.dat")
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        if (sum(abs(S_Q(i,j,:))) > .1) then
          write(333, *) "#Partial Wave:", n, l, m, n_prime, l_prime, m_prime
          do i_Q_ = 1, num_Q_
            write(333, *) lin_geoms(i_Q_), real(S_Q(i,j,i_Q_)), aimag(S_Q(i,j,i_Q_))
          end do
          write(333, *)
        end if
      end do
    end do
    close(333)

    call frame_transform(S_Q, ion_wf, S_ft, ft_channels, energy_array)

    open(333, file = "ft_channels.dat")
    do i = 1, num_tot
      write(333,*) i, ft_channels(i,:)
    end do
    close(333)

    call replace_like_vals(S_ft(:,:), ft_channels, replaced_S_ft)

    call sort_index(energy_array, index)

    do i = 1, num_tot
      do j = 1, num_tot
        sorted_S_ft(i,j) = replaced_S_ft(index(i),index(j))
      end do
    end do

    do i = 1, num_tot
      sorted_ft_channels(i, :) = ft_channels(index(i), :)
    end do

    call channel_elimination(sorted_S_ft, energy_array, sorted_ft_channels, eics_matrix)

    
    !---------------------------------------------------------------------------------------------------------------
    !                                          WRITING TO FILES
    !---------------------------------------------------------------------------------------------------------------


    open(27, file = "channel_energies.dat")
    do ei = 1, n_max
      do vi = 0, vib_num-1

        write(27,*) "#ei,vi", ei, vi
        write(27,*) chan_ens(ei,vi)
        write(27,*)

      end do
    end do

    open(27, file = "check_like_elements_replaced.dat")
    do i = 1, num_tot
      n = ft_channels(i,1)
      v = ft_channels(i,2)
      l = ft_channels(i,3)
      m = ft_channels(i,4)
      do j = 1, num_tot
        n_prime = ft_channels(j,1)
        v_prime = ft_channels(j,2)
        l_prime = ft_channels(j,3)
        m_prime = ft_channels(j,4)

        if ((n == 1 .and. l == 1 .and. (m == 1 .or. m == -1)) .and. (n_prime == 2 .or. n_prime == 3 )) then
          if (replaced_S_ft(i,j) /= 0) then
            write(27, *) m, n_prime, v_prime, l_prime, m_prime, replaced_S_ft(i,j)
          end if
        end if

      end do
    end do
    close(27)

    ! open(333, file = "fixed_kmat.dat")
    ! do i=1, dim_lm
    !   do j=1, dim_lm
    !     call i2lm(i, n, l, m)
    !     call i2lm(j, n_prime, l_prime, m_prime)
    !     if (fixed_K_Q(i,j,1,i_energy) /= 0) then
    !       write(333, *) "#Partial Wave:", n, l, m, n_prime, l_prime, m_prime
    !       do i_Q = 1, num_Q
    !         write(333, *) geoms(i_Q), fixed_K_Q(i,j,i_Q,i_energy)
    !       end do
    !       write(333, *)
    !     end if
    !   end do
    ! end do
    ! close(333)


end program

subroutine add_nice_elems(fixed_Delta_Q, lin_Delta_Q, final_Delta_Q)

  use global_params, only: dim_lm, num_Q, num_Q_, good_elem_list

  real(8), dimension(dim_lm, dim_lm, num_Q), intent(in)      :: fixed_Delta_Q
  real(8), dimension(dim_lm, dim_lm, num_Q_), intent(in)      :: lin_Delta_Q
  real(8), dimension(dim_lm, dim_lm, num_Q_), intent(out)     :: final_Delta_Q

  integer ::  i, j, num_elems, ind, ind_p

  real(8), dimension(num_Q)      ::  in_elem
  real(8), dimension(num_Q_)     ::  fin_elem

  num_elems = size(good_elem_list, 1)
  print *, num_elems

  final_Delta_Q = lin_Delta_Q
  do i = 1, num_elems

    call lm2i(good_elem_list(i,1), good_elem_list(i,2), good_elem_list(i,3), ind)
    call lm2i(good_elem_list(i,4), good_elem_list(i,5), good_elem_list(i,6), ind_p)
    print *, good_elem_list(i,1), good_elem_list(i,2), good_elem_list(i,3), ind, good_elem_list(i,4), good_elem_list(i,5), good_elem_list(i,6), ind_p
    in_elem = fixed_Delta_Q(ind, ind_p, :)

    call change_grid(in_elem, fin_elem)

    final_Delta_Q(ind, ind_p, :) = fin_elem

  end do

end subroutine

subroutine change_grid(in_elem, fin_elem)

  use global_params, only:  num_Q, num_Q_, geoms, lin_geoms

  real(8), dimension(num_Q), intent(in)       ::  in_elem
  real(8), dimension(num_Q_), intent(out)     ::  fin_elem

  real(8), dimension(num_Q)                   ::  cm

  integer   ::    i, j
  real(8)   ::    SPL

  real(8), dimension(num_Q) ::  sub_mat

  call spline(num_Q, geoms, in_elem, cm)
  do i = 1, num_Q_
    fin_elem(i) = SPL(num_Q, geoms, in_elem, cm, lin_geoms(i))
  end do


end subroutine

subroutine trim_Delta(Delta_Q, trimmed_Delta_Q)

  use global_params, only: dim_lm, num_Q, trim_value, lin_start, lin_end, deriv_trim, geoms

  real(8), dimension(dim_lm, dim_lm, num_Q), intent(in)       :: Delta_Q
  real(8), dimension(dim_lm, dim_lm, num_Q), intent(out)      :: trimmed_Delta_Q

  real(8) ::  deriv, avg

  integer ::  i, j

  trimmed_Delta_Q = Delta_Q
  do i = 1, dim_lm
    do j = 1, i

      if (sum(abs(Delta_Q(i,j,:))) <= trim_value) then
        trimmed_Delta_Q(i,j,:) = 0
        trimmed_Delta_Q(j,i,:) = 0
      end if

      call find_rderiv(trimmed_Delta_Q(i, j, lin_end), trimmed_Delta_Q(i, j, lin_start), (geoms(lin_end) - geoms(lin_start)), deriv)

      if (abs(deriv) > deriv_trim) then

        avg = sum(trimmed_Delta_Q(i,j,:))/size(trimmed_Delta_Q(i,j,:))

        trimmed_Delta_Q(i,j,:) = avg
        trimmed_Delta_Q(j,i,:) = avg

      end if

    end do
  end do

  

end subroutine

subroutine channel_elimination(S_ft, energy_array, sorted_ft_channels, eics_matrix)

  use global_params, only: num_tot, chan_ens, energy_step_num, pi, ci, n_max, vib_num, elec_energies
  use iso_fortran_env, only: int64, iostat_end

  complex*16, dimension(num_tot, num_tot), intent(in)        :: S_ft
  real(8), dimension(num_tot), intent(in)                    :: energy_array
  real(8), dimension(n_max, n_max, 0:vib_num, 0:vib_num), intent(out) ::  eics_matrix
  integer ::  max_chan
  integer, dimension(num_tot, 4), intent(in)                 :: sorted_ft_channels
  real(8), dimension(num_tot) :: nu


  complex*16, allocatable       :: phys_S_ft(:,:)

  real(8) ::  energy_begin, energy_end, energy_step, E, eics_sum, E_elec

  integer ::  en_i, ei, ei_p, vi, vi_p, num_open, num_closed, j, INFO, i, n, v, n_p, v_p, k

  complex*16, allocatable		::	smat_cc(:,:), smat_co(:,:), smat_oc(:,:), smat_oo(:,:)
  complex*16, allocatable		::	beta(:)
  integer,    allocatable		::	IPIV(:)
  

  energy_begin = chan_ens(1,0) !Origin is ground state of N2+
  energy_end = chan_ens(n_max, vib_num-1)
  energy_step = (energy_end - energy_begin)/energy_step_num
  print *, energy_begin, energy_end, energy_step

  open(27, file = "eics_vs_energy_fromground.dat")
  open(28, file = "eics_vs_energy_toground.dat")
  open(333, file = "Resonances_exc.dat")

  eics_sum = 0
  do en_i = 1, energy_step_num

    if (en_i == int(energy_step_num/4)) then
      print *, "Quarter of the way done."
    else if (en_i == int(energy_step_num/2)) then
      print *, "Half of the way done."
    else if (en_i == int(.75*energy_step_num)) then
      print *, "Three quarters of the way done."
    end if

    E = energy_begin + en_i*energy_step

    num_open = 0
    Do j=1,num_tot
      If(energy_array(j)>E)Exit
    Enddo
    num_open=j-1
    num_closed = num_tot - num_open
    
    allocate(IPIV(num_tot), phys_S_ft(num_open,num_open))
    allocate(smat_cc(num_closed,num_closed),smat_co(num_closed,num_open),smat_oc(num_open,num_closed), &
              & smat_oo(num_open,num_open),beta(num_tot))

    if (num_open /= num_tot) then

      smat_oo = S_ft(1:num_open,1:num_open)
      smat_cc = S_ft(num_open+1:num_tot,num_open+1:num_tot)
      smat_oc = S_ft(1:num_open,num_open+1:num_tot)
      smat_co = S_ft(num_open+1:num_tot,1:num_open)

      beta = 0d0
      do j = 1, num_closed
        beta(j+num_open) = pi/sqrt(2*(energy_array(j+num_open)-E))
      end do

      nu = pi * beta   

      call Z_analysis(en_i, sorted_ft_channels, nu, smat_co, num_closed, num_open)

      do j = 1, num_closed
      smat_cc(j,j) = smat_cc(j,j) - exp(-2d0*ci*beta(j+num_open))
      end do

      call ZGESV(num_closed,num_open,smat_cc,num_closed,IPIV(1:num_closed),smat_co,num_closed,INFO)

      phys_S_ft(:,:) =  MatMul(smat_oc,smat_co)

    else if (num_open == num_tot) then

      phys_S_ft(:,:) = S_ft

    end if

    print *, en_i
    open(29, file = "Sum_Check.dat")
      do ei = 1, n_max
        do ei_p = 1, n_max
          do vi = 0, vib_num-1
            do vi_p = 0, vib_num-1
  
              eics_sum = 0
              do i = 1, num_open
                n = sorted_ft_channels(i,1)
                v = sorted_ft_channels(i,2)
   
                do j = 1, num_open
                  n_p = sorted_ft_channels(j,1)
                  v_p = sorted_ft_channels(j,2)
  
                  if ( (n == ei) .and. (n_p == ei_p) .and. (v == vi) .and. (v_p == vi_p) ) then
              
                      eics_sum = eics_sum + ( real(phys_S_ft(i,j))**2 + aimag(phys_S_ft(i,j))**2 )
  
                     if (phys_S_ft(i,j) /= 0) then
                      write(29, *) n, v, sorted_ft_channels(i,3), sorted_ft_channels(i,4), n_p, v_p, sorted_ft_channels(j,3), sorted_ft_channels(j,4), phys_S_ft(i,j)
                     end if
                     
                  end if
  
                end do
              end do

              if (E - chan_ens(ei_p,vi_p) /= 0) then
                eics_matrix(ei, ei_p, vi, vi_p) =  eics_sum  
              end if

            end do
          end do
        end do
      end do
      close(29)

      write(27,*) E, (pi / (2*E) ) * eics_matrix(1,1,0,0), (pi / (2*E) ) * eics_matrix(1,1,0,1), (pi / (2*E) ) * eics_matrix(1,1,0,2), (pi / (2*E) ) * eics_matrix(1,1,0,3)
      write(28,*) elec_energies(2,0,en_i), (pi / (2*elec_energies(2,0,en_i)) ) * eics_matrix(2,1,0,0), elec_energies(2,1,en_i), (pi / (2*elec_energies(2,1,en_i)) ) * eics_matrix(2,1,1,0), &
                  & elec_energies(2,2,en_i), (pi / (2*elec_energies(2,2,en_i)) ) * eics_matrix(2,1,2,0), elec_energies(2,3,en_i), (pi / (2*elec_energies(2,3,en_i)) ) * eics_matrix(2,1,3,0)

    deallocate(phys_S_ft, smat_cc, smat_co, smat_oc, smat_oo, IPIV, beta)

  end do
  close(27)
  close(28)
  close(333)

end subroutine

subroutine Z_analysis(en_i, sorted_ft_channels, nu, smat_co, num_closed, num_open)

  use global_params, only: num_tot, energy_step_num, n_max, vib_num, chan_ens

  real(8), intent(in) ::  nu(num_tot)
  integer, intent(in) ::  num_closed, num_open, en_i
  complex*16, dimension(num_closed, num_open)  ::   smat_co
  integer, dimension(num_tot, 4), intent(in)                 :: sorted_ft_channels

  integer ::  j,k
  integer, dimension(3) ::  max_chans
  real(8), dimension(num_closed)  :: Z_norms
  real(8) ::  energy_begin, energy_end, energy_step, E
  complex*16, dimension(num_closed, num_open)  :: D

  energy_begin = chan_ens(1,0) !Origin is ground state of N2+
  energy_end = chan_ens(n_max, vib_num-1)
  energy_step = (energy_end - energy_begin)/energy_step_num

  E = energy_begin + en_i*energy_step

  do k = 1, num_open
    do j = 1, num_closed
      D(j, k) = smat_co(j,k) * nu(j+num_open)**1.5
    end do
  end do

  Z_norms = 0
  do j = 1, num_closed
    !do coord_i = 1, 3
      Z_norms(j) = Z_norms(j) + sum(abs(D(j,:))**2)
    !end do
  end do

  do j = 1, 3
    max_chans(j) = maxloc(Z_norms, 1) 
    Z_norms(maxloc(Z_norms, 1)) = 0
  end do

  
  !print *, ei, max_chans(ei,1)
  write(333,'(F11.8,3(4I4,2F10.3))') E, sorted_ft_channels(max_chans(1)+num_open,1), sorted_ft_channels(max_chans(1)+num_open,2), sorted_ft_channels(max_chans(1)+num_open,3), sorted_ft_channels(max_chans(1)+num_open,4), Z_norms(max_chans(1)), nu(max_chans(1)),&
        & sorted_ft_channels(max_chans(2)+num_open,1), sorted_ft_channels(max_chans(2)+num_open,2), sorted_ft_channels(max_chans(2)+num_open,3), sorted_ft_channels(max_chans(2)+num_open,4), Z_norms(max_chans(2)), nu(max_chans(2)), &
        & sorted_ft_channels(max_chans(3)+num_open,1), sorted_ft_channels(max_chans(3)+num_open,2), sorted_ft_channels(max_chans(3)+num_open,3), sorted_ft_channels(max_chans(3)+num_open,4), Z_norms(max_chans(3)), nu(max_chans(3))
  write(333,*)

end subroutine

subroutine replace_like_vals(S_ft, ft_channels, replaced_S_ft)

  use global_params, only: num_tot, num_Q, vib_num, l_max

  complex*16, dimension(num_tot, num_tot), intent(in)         :: S_ft
  integer, dimension(num_tot, 4), intent(in)                    :: ft_channels
  complex*16, dimension(num_tot, num_tot), intent(out)       :: replaced_S_ft

  integer ::  i, j, k, n, v, l, m, n_p, v_p, l_p, m_p, n_r, v_r, l_r, m_r, n_p_r, v_p_r, l_p_r, m_p_r, j_r, i_r

  
    replaced_S_ft(:,:) = S_ft
    do i = 1, num_tot

      if (i == int(num_tot/4)) then
        print *, "Quarter of the way done."
      else if (i == int(num_tot/2)) then
        print *, "Half of the way done."
      else if (i == int(.75*num_tot)) then
        print *, "Three quarters of the way done."
      end if

      n = ft_channels(i,1)
      v = ft_channels(i,2)
      l = ft_channels(i,3)
      m = ft_channels(i,4)
      if (n /= 1) then; cycle; end if
      
      do j = 1, num_tot
        n_p = ft_channels(j,1)
        v_p = ft_channels(j,2)
        l_p = ft_channels(j,3)
        m_p = ft_channels(j,4)
        if (n_p /= 2 .and. n_p /= 3) then; cycle; end if

        if (n_p == 2) then

          if ( m == 0 ) then
            n_r = n; v_r = v; l_r = l; m_r = m; n_p_r = 3; v_p_r = v_p; l_p_r = l_p; m_p_r = -m_p
            do j_r = 1, num_tot
              if (ft_channels(j_r,1) == n_p_r .and. ft_channels(j_r,2) == v_p_r .and. ft_channels(j_r,3) == l_p_r .and. ft_channels(j_r,4) == m_p_r) then
                replaced_S_ft(i,j_r) = replaced_S_ft(i,j)
                replaced_S_ft(j_r,i) = replaced_S_ft(i,j)
              end if
            end do

          else if ( m /= 0 ) then
            n_r = n; v_r = v; l_r = l; m_r =-m; n_p_r = 3; v_p_r = v_p; l_p_r = l_p; m_p_r = m_p
            do i_r = 1, num_tot
              do j_r = 1, num_tot
                if (ft_channels(i_r,1) == n_r .and. ft_channels(i_r,2) == v_r .and. ft_channels(i_r,3) == l_r .and. ft_channels(i_r,4) == m_r) then
                  if (ft_channels(j_r,1) == n_p_r .and. ft_channels(j_r,2) == v_p_r .and. ft_channels(j_r,3) == l_p_r .and. ft_channels(j_r,4) == m_p_r) then
                    replaced_S_ft(i_r,j_r) = replaced_S_ft(i,j)
                    replaced_S_ft(j_r,i_r) = replaced_S_ft(i,j)
                  end if
                end if
              end do
            end do
          end if

        end if

      end do
    end do

end subroutine

subroutine check_unitarity(S_Q)

  use global_params, only: dim_lm, num_Q

  complex*16, dimension(dim_lm, dim_lm, num_Q)  ::  S_Q

  real(8) ::  sum

  integer ::  i_Q, i, j

  do i_Q = 1, num_Q
    do i = 1, dim_lm
      sum = 0
      do j = 1, dim_lm
        sum = sum + (real(S_Q(i,j,i_Q))**2 + aimag(S_Q(i,j,i_Q))**2)
      end do
      print *, i_Q, i, sum
    end do
  end do

end subroutine

! subroutine calculate_eics(S_ft, num_open, sorted_ft_channels, eics_matrix)

!   use global_params, only: num_tot, n_max, vib_num, l_max
!   use iso_fortran_env, only: int64, iostat_end 

!   integer, intent(in)                                       :: num_open
!   complex*16, dimension(num_open, num_open), intent(in)       :: S_ft
!   integer, dimension(num_tot, 4), intent(in)                :: sorted_ft_channels
!   real(8), dimension(n_max, n_max, 0:vib_num, 0:vib_num), intent(out)             :: eics_matrix

!   real(8) ::  eics_sum
!   integer ::  i, j, ei, ei_p, vi, vi_p
!   integer, dimension((l_max+1)**2 * vib_num)  ::  n_trim
!   integer, dimension((l_max+1)**2 * vib_num, 4)  ::  n_trimmed_chans
!   integer, dimension((l_max+1)**2)  ::  v_trim
!   integer, dimension((l_max+1)**2)  ::  in_indices, fin_indices

!   !do i_energy = 1, num_scatE
!     do ei = 1, n_max
!       do vi = 0, vib_num-1

!         n_trim = findloc(sorted_ft_channels(:,1), ei, 1)
!         do i = 1, (l_max+1)**2 * vib_num
!           n_trimmed_chans(i,:) = sorted_ft_channels(n_trim(i),:)
!         end do
!         v_trim = findloc(n_trimmed_chans(:,2), vi, 1)
!         do i = 1, (l_max+1)**2 
!           in_indices(i) = n_trim(v_trim(i))
!         end do

!         do ei_p = 1, n_max
!           do vi_p = 0, vib_num-1

!             n_trim = findloc(sorted_ft_channels(:,1), ei_p, 1)
!             do i = 1, (l_max+1)**2 * vib_num
!               n_trimmed_chans(i,:) = sorted_ft_channels(n_trim(i),:)
!             end do
!             v_trim = findloc(n_trimmed_chans(:,2), vi_p, 1)
!             do i = 1, (l_max+1)**2 
!               fin_indices(i) = n_trim(v_trim(i))
!             end do

!             eics_sum = 0
!             do i = 1, (l_max+1)**2
!               do j = 1, (l_max+1)**2
!                 eics_sum = eics_sum + ( real(S_ft(in_indices(i),fin_indices(j)))**2 + aimag(S_ft(in_indices(i),fin_indices(j)))**2 )
!               end do
!             end do

!             ! do i = 1, num_tot
!             !   n = ft_channels(index(i),1)
!             !   v = ft_channels(index(i),2)
 
!             !   do j = 1, num_tot
!             !     n_p = ft_channels(index(j),1)
!             !     v_p = ft_channels(index(j),2)

!             !     if ( (n == ei) .and. (n_p == ei_p) .and. (v == vi) .and. (v_p == vi_p) ) then
!             !        eics_sum = eics_sum + ( real(S_ft(index(i),index(j)))**2 + aimag(S_ft(index(i),index(j)))**2 )
!             !        if (S_ft(index(i),index(j)) /= 0) then
!             !         write(27, *) n, v, ft_channels(index(i),3), ft_channels(index(i),4), n_p, v_p, ft_channels(index(j),3), ft_channels(index(j),4), S_ft(index(i),index(j))
!             !        end if
!             !     end if

!             !   end do
!             ! end do
!             eics_matrix(ei, ei_p, vi, vi_p) = eics_sum

!           end do
!         end do
!       end do
!     end do
!   !end do
  

! end subroutine

subroutine frame_transform(S_Q, ion_wf, S_ft, ft_channels, energy_array)

  use global_params, only : vib_num, n_max, num_Q_, dim_lm, num_tot, lin_geoms, all_open_i, chan_ens, lin_r_step
  
  
  complex*16, dimension(dim_lm, dim_lm, num_Q_), intent(in) ::  S_Q
  real(8), dimension(0:vib_num-1, num_Q_, n_max), intent(in)	::	ion_wf
  
  complex*16, dimension(num_tot, num_tot), intent(out) ::  S_ft
  integer, dimension(num_tot, 4), intent(out) ::  ft_channels
  real(8), dimension(num_tot), intent(out)    ::  energy_array
  
  complex*16 :: S_sum
  integer :: i, j, iv, i_p, j_p, iv_p, i_Q_, n, l, m, n_p, l_p, m_p, ie 
  
  ie = maxval(all_open_i)
  open(27, file = "S_ft.dat")
  !do ie = 1, num_scatE
    do i = 1, dim_lm
      do iv = 0, vib_num-1

        j = (i-1)*vib_num + (iv+1)

        call i2lm(i, n, l, m)
        ft_channels(j, :) = (/n, iv, l, m/)
        energy_array(j) = chan_ens(n, iv)

        do i_p = 1, dim_lm
          call i2lm(i_p, n_p, l_p, m_p)
            do iv_p = 0, vib_num-1

              j_p = (i_p-1)*vib_num + (iv_p+1)

              S_sum = 0
              do i_Q_ = 1,num_Q_
                S_sum = S_sum + S_Q(i, i_p, i_Q_) * ion_wf(iv, i_Q_, n) * ion_wf(iv_p, i_Q_, n_p) * lin_r_step
              end do
              S_ft(j, j_p) = S_sum

             
                if (S_ft(j, j_p) /= 0) then
                  if (iv_p == 0) then
                    write(27,*) "#e=", n,"v=", iv, "l=", l, "m=", m, "#e_p=", n_p, "l_p=", l_p, "m_p=", m_p
                  end if
                  write(27,*) "v_p=", iv_p, real(S_ft(j, j_p)), aimag(S_ft(j, j_p))
                end if
              

            end do
            if (S_ft(j, j_p) /= 0) then
              write(27,*)
            end if
        end do

      end do
    end do
  !end do
  close(27)
  
  end subroutine

subroutine get_wf(ion_wf)

  use global_params, only: n_max, num_Q_, vib_num, lin_geoms

  real(8), dimension(0:vib_num-1, num_Q_, n_max), intent(out)     :: ion_wf

  character(len = 44) ::  wf_name

  integer ::  elec_i, i_Q_, vib_i

  do elec_i = 0, n_max-1
    do vib_i = 0, vib_num-1
      call make_wf_name(elec_i, vib_i, wf_name)
      open(333, file = wf_name, status = "old")
      do i_Q_ = 1, num_Q_
          call get_wf_vals(wf_name, lin_geoms(i_Q_), ion_wf(vib_i, i_Q_, elec_i+1))
      end do
      close(333)
    end do
  end do

end subroutine


subroutine make_wf_name(elec_num, vib_i, wf_name)

  use global_params, only : wf_dir, state_names

  integer, intent(in)     ::  elec_num, vib_i
  character (len = 2)     ::  cex2
  character(len = 44), intent(out) ::  wf_name
  integer                 ::  cnj, unj, dnj

  cnj = vib_i/100
  dnj = (vib_i - cnj*100)/10
  unj = vib_i - cnj*100 - dnj*10
  ! cex2 = char(cnj+48)//char(dnj+48)//char(unj+48)
  cex2 = char(dnj+48)//char(unj+48)

  print *, cex2

  wf_name = wf_dir // state_names(elec_num) // "/Target_wf0" // cex2 // "_N2__ElecState001.dat"

end subroutine

subroutine get_wf_vals(file_path, r, WF)

  use global_params, only : vib_num

  real(8)                  :: r
  integer                  :: a, io, vib_i
  character(len = 44), intent(in) ::  file_path
  real(8), intent(out)     :: WF
  real(8)    :: WF_data(1000), r_mat(1000), sub_mat(1000), x1, x2, y1, y2

 
  r_mat = 100
  sub_mat = 100
  WF_data = 0
  open(333, file = file_path, status = "old")
  do a = 1, 1000
      read(333, *, iostat=io) r_mat(a), WF_data(a)
      if (io /= 0) exit
      !print *, WF_data(:,a)
  end do
  close(333)
  
  
  !Getting norm, but imaginary parts are tiny so its basically just the real part
  !do a = 1, 1000
      !do vib_i = 0, vib_num
          !print 
          !WF_mags(vib_i, a) = sqrt(WF_data( (vib_i * 2)+1, a)**2 + WF_data( ((vib_i+1) * 2), a)**2)
          !print *, WF_mags(:,a)
      !end do
  !end do


  sub_mat = abs(r_mat - r)
  !print *, "r:", r
  !print *, "loc:", minloc(sub_mat, 1)
  !print *, sub_mat(minloc(sub_mat, 1))
  
  if ( (r_mat(minloc(sub_mat, 1)) - r) > 0 ) then
    x1 = r_mat(minloc(sub_mat, 1)-1)
    x2 = r_mat(minloc(sub_mat, 1))
    y1 = WF_data(minloc(sub_mat, 1)-1)
    y2 = WF_data(minloc(sub_mat, 1))
    WF = y1 + (r-x1)*( (y2-y1) / (x2-x1) )
    
  else if ((r_mat(minloc(sub_mat, 1)) - r) < 0) then
    x1 = r_mat(minloc(sub_mat, 1))
    x2 = r_mat(minloc(sub_mat, 1)+1)
    y1 = WF_data(minloc(sub_mat, 1))
    y2 = WF_data(minloc(sub_mat, 1)+1)
    WF = y1 + (r-x1)*( (y2-y1) / (x2-x1) )

  else if ((r_mat(minloc(sub_mat, 1)) - r) == 0) then
    WF = WF_data(minloc(sub_mat, 1))

  end if
  !print *, "WF:", WF


end subroutine

subroutine sign_fixer(Delta_Q, fixed_Delta_Q)

  use global_params, only: num_Q, dim_lm, geoms, all_open_i

  real(8), intent(in)     :: Delta_Q(dim_lm, dim_lm, num_Q)
  real(8), intent(out)    :: fixed_Delta_Q(dim_lm, dim_lm, num_Q)

  real(8) ::  old_deriv, deriv_1, deriv_2
  integer ::  i, j, i_Q

  fixed_Delta_Q = Delta_Q
  old_deriv = 0
    do i = 1, dim_lm
      do j = 1, dim_lm
        do i_Q = 1, num_Q
          if (i_Q > 1) then
            call find_rderiv(fixed_Delta_Q(i, j, i_Q), fixed_Delta_Q(i, j, i_Q - 1), geoms(i_Q) - geoms(i_Q-1), deriv_1)
            call find_rderiv(-fixed_Delta_Q(i, j, i_Q), fixed_Delta_Q(i, j, i_Q - 1), geoms(i_Q) - geoms(i_Q-1), deriv_2)

                if (i_Q == 2) then
                  if (abs(deriv_2) <= abs(deriv_1)) then
                    fixed_Delta_Q(i, j, i_Q) = -fixed_Delta_Q(i, j, i_Q)
                    !fixed_Delta_Q(j, i, i_Q, :) = -fixed_Delta_Q(j, i, i_Q, :)
                  end if
                else
                  if (abs(deriv_2 - old_deriv) <= abs(deriv_1 - old_deriv)) then
                    fixed_Delta_Q(i, j, i_Q) = -fixed_Delta_Q(i, j, i_Q)
                    !fixed_Delta_Q(j, i, i_Q, :) = -fixed_Delta_Q(j, i, i_Q, :)
                  end if
                end if

                call find_rderiv(fixed_Delta_Q(i, j, i_Q), fixed_Delta_Q(i, j, i_Q - 1), geoms(i_Q) - geoms(i_Q-1), old_deriv)
            end if
        end do
      end do
    end do


end subroutine

subroutine make_abs_states(md, D2h_Dinf, fixed_D2h_Dinf)

  use global_params, only: num_gamma_Q, num_Q, name_of_file, geoms

  integer, dimension(100, num_gamma_Q, 3, num_Q), intent(out)    ::  fixed_D2h_Dinf
  integer, dimension(100, num_gamma_Q, 3, num_Q), intent(in)     ::  D2h_Dinf
  integer, dimension(num_Q,num_gamma_Q), intent(in)      ::  md

  integer ::  i, i_name, i_Q, abs_states(100, num_gamma_Q, 3)
  logical, dimension(num_Q) ::  switch_mask

  print *, "called abs states"

  !Making 1Ag, 6th geom the reference. Seems like switch only depends on geom
  do i_name = 1, num_gamma_Q
    abs_states(1:md(1,i_name), i_name, :) = D2h_Dinf(1:md(1,i_name), i_name, :, 6)
  end do

  !Now we make a mask which tells which geoms flipped
  switch_mask = .false.
  do i_Q = 1, num_Q
    do i = 1, md(i_Q, 1) !just 1Ag sym
      !print *, D2h_Dinf(i, 1, 1, i_Q)
      if (D2h_Dinf(i, 1, 1, i_Q) /= 1) then !Only look at excited states, if n /= 1
        !print *, D2h_Dinf(i, 1, 3, i_Q), abs_states(i, 1, 3)
        if ( (D2h_Dinf(i, 1, 3, i_Q) /= abs_states(i, 1, 3)) .and. (D2h_Dinf(i, 1, 1, i_Q) == abs_states(i, 1, 1))) then !Compare l and m vals and see if they are the same
          switch_mask(i_Q) = .true.
          print *, "Syms flipped,", i_Q
          exit
        end if
      end if
    end do
  end do

  !Flip all the geoms back
  fixed_D2h_Dinf = D2h_Dinf
  do i_name = 1, num_gamma_Q !loops over syms
    do i_Q = 1, num_Q         !loops over geom
      do i = 1, md(i_Q, i_name) !loops over channel

        if (D2h_Dinf(i, i_name, 1, i_Q) /= 1) then !specify exc states
          if (switch_mask(i_Q)) then                !should it be flipped?

            if (D2h_Dinf(i, i_name, 1, i_Q) == 2) then  !change n=2 -> n=3 in fixed
              fixed_D2h_Dinf(i, i_name, 1, i_Q) = 3
              !print *, "2 -> 3", i_Q
            else if (D2h_Dinf(i, i_name, 1, i_Q) == 3) then !change n=3 -> n=2 in fixed
              fixed_D2h_Dinf(i, i_name, 1, i_Q) = 2
              !print *, "3 -> 2", i_Q
            end if
          
          end if
        end if

      end do
    end do
  end do

  do i_name = 1, num_gamma_Q
    open(27, file="CHECKFIXED_" // name_of_file(i_name))
    write(27,*) "ABS STATES:"
    do i = 1, md(6, i_name)
      write(27,*) abs_states(i, i_name, :)
    end do
    write(27,*)
    do i_Q = 1, num_Q
      write(27,*) geoms(i_Q)
      do i = 1, md(i_Q, i_name)
        write(27,*) fixed_D2h_Dinf(i,i_name,:,i_Q)
      end do
      write(27,*)
    end do
    close(27)
  end do

end subroutine

subroutine get_KQ(K_Q, md, D2h_Dinf, energies)

    use global_params, only: num_Q, num_gamma_Q, ordered_geoms, num_scatE, n_max, dim_lm, skip_en, num_scatE, all_open_i, name_of_file

    real(8), intent(out)   :: K_Q(dim_lm, dim_lm, num_Q)
    
    real(8), dimension(n_max, num_Q)            :: thresh_ens
    real(8), allocatable, dimension(:)          :: temp_arr
    real(8), allocatable, dimension(:,:)        :: k_mtrx
    integer, dimension(100,num_gamma_Q,3,num_Q), intent(out)    :: D2h_Dinf
    real(8), dimension(num_scatE), intent(out)  :: energies

    integer                                   ::  i, i_name, j, i_rec, records(5,2000), itmp(1:4), i_Q, dims_extended(num_gamma_Q), &
                                                  & num_recs, n_num, a1, b, c, skip, block_size, ie, my_count, i2, j2
    integer, dimension(num_Q, num_gamma_Q), intent(out)      ::  md
    character(len=2)                          ::  cex2
    real(8), dimension(2000)                  ::  ener_ch_lm(2000)
    real(8)                                   ::  E
    
    geom_Q1: Do i_Q=1,num_Q

       i_rec=0
       do i_name = 1, num_gamma_Q
         !print *, name_of_file(i_name)   
         
         open(unit = 7,&
         file= 'K-mats/'//cex2(i_Q)//'/matrices/'//name_of_file(i_name),&
         action='read')
       
         do j = 1, 2;  read(7,*) ; enddo !- ignore useless lines
         read(7,*) itmp(1:4)
         dims_extended(i_name) = itmp(4)
        
         !D_infty,v symmetric str.
         ! md(1:8)=(/6, 1, 3, 3, 3, 3, 3, 3 /) 
   
   
         if(i_name==1)md(i_Q,i_name)=15
         if(i_name==2)md(i_Q,i_name)=10
         if(i_name==3)md(i_Q,i_name)=10
         if(i_name==4)md(i_Q,i_name)=15
         if(i_name==5)md(i_Q,i_name)=10
         if(i_name==6)md(i_Q,i_name)=15
         if(i_name==7)md(i_Q,i_name)=10
         if(i_name==8)md(i_Q,i_name)=15
         ! if(i_name==9)md(i_Q,i_name)=12
         ! if(i_name==10)md(i_Q,i_name)=7
         ! if(i_name==11)md(i_Q,i_name)=9
         ! if(i_name==12)md(i_Q,i_name)=9
         ! if(i_name==13)md(i_Q,i_name)=7
         ! if(i_name==14)md(i_Q,i_name)=12
         ! if(i_name==15)md(i_Q,i_name)=7
         ! if(i_name==16)md(i_Q,i_name)=12
   
         do j = 1, itmp(1)+1;  read(7,*) ; enddo !- ignore useless lines
   
         do j = 1, md(i_Q, i_name)
           i_rec=i_rec+1
           records(1,i_rec)=i_name
           read(7,*)records(2,i_rec),records(3,i_rec),records(4,i_rec),records(5,i_rec),ener_ch_lm(i_rec)
           !Write(6,*)'iQ=',i_Q,records(1:5,i_rec),ener_ch_lm(i_rec)
           !print *, records(2,i_rec),records(3,i_rec),records(4,i_rec) 
           
            do i = 1, md(i_Q, i_name)
              if (records(3,i_rec) == 5 .or. records(3,i_rec) == 7) then 
                records(3,i_rec) = 4
              end if
            end do
           

           if (i_Q == 1 .or. i_Q == 2 .or. i_Q == 3 .or. i_Q == 4) then
            do i = 1, md(i_Q, i_name)
              if (records(3,i_rec) == 4) then 
                records(3,i_rec) = 2
              else if (records(3,i_rec) == 2) then 
                records(3,i_rec) = 4
              end if
            end do
           end if
           thresh_ens(records(3,i_rec), i_Q) = ener_ch_lm(i_rec)
           D2h_Dinf(records(2,i_rec),i_name,1:3,i_Q) = records(3:5,i_rec) ! i_el, l,lambda
         enddo
         ! Quartet crossing makes quantemol call B state 5 instead of 4 in 15 and 16

         Close(7)
   
         ! if (i_Q == num_Q .and. i_name == num_gamma_Q) then
         !   print *, "HERE"
         !   print *, D2h_Dinf(4,3,:,i_Q-1)
         !   print *, D2h_Dinf(4,3,:,i_Q)
         !   call lm2i(D2h_Dinf(4,3,1,i_Q-1), D2h_Dinf(4,3,2,i_Q-1), D2h_Dinf(4,3,3,i_Q-1),i)
         !   call lm2i(D2h_Dinf(4,3,1,i_Q), D2h_Dinf(4,3,2,i_Q), D2h_Dinf(4,3,3,i_Q),i2) 
         !   print *, i, i2
         ! end if
   
       enddo
    end do geom_Q1

    do i = 1, md(1,1)
    print *, D2h_Dinf(i,1, :,2)
    end do
    call make_abs_states(md, D2h_Dinf, D2h_Dinf)
    do i = 1, md(1,1)
    print *, D2h_Dinf(i,1, :,6)
    end do

    geom_Q2: Do i_Q=1,num_Q
       num_recs=i_rec
    
       name_of_file(1:num_gamma_Q) = (/ '1Ag.kmat ', '1Au.kmat ', '1B1g.kmat', '1B1u.kmat', &
                                           '1B2g.kmat', '1B2u.kmat', '1B3g.kmat', '1B3u.kmat'/)!, &
                                                         !'3Ag.kmat ', '3Au.kmat ', '3B1g.kmat', '3B1u.kmat', &
                                                                     !'3B2g.kmat', '3B2u.kmat', '3B3g.kmat', '3B3u.kmat' /)
   
       file_read:  do i_name = 1, num_gamma_Q
         ! print *, md(1, i_name)
         !- Matrix allocation
         !n_num=md(i_Q,i_name)*(md(i_Q,i_name)+1)/2 ! it works only in dimension does not change with energy (md*md+1)/2 ! +1 for the diagonal and md*md/2 off-diag
         n_num=records(3,num_recs)*dim_lm*(records(3,num_recs)*dim_lm+1)/2 
         allocate( temp_arr(n_num), k_mtrx(md(i_Q,i_name),md(i_Q,i_name)) )
        
         ! print *, name_of_file(i_name)
      
         open(unit = 7, file= 'K-mats/'//cex2(i_Q)//'/matrices/'//name_of_file(i_name),&
          action='read')
         do j = 1, 5;  read(7,*) ; enddo !- ignore useless lines
         !- loop over energy
       
         skip = 0
         energy : Do ie = 1, num_scatE
         !print *, skip
       
             read(7,*,end=1) a1, b, c, E

              !print *, "a, b, c, E :", a1, b, c, E, i_Q ! /2.d0 * au_eV
   
             if (mod(c,4)==0) then
               block_size = c/4
             else
               block_size = ceiling(real(c/4.))
             end if
   
             !print *, E, thresh_ens(n_max, i_Q)
             
             if (E < maxval(thresh_ens(:, i_Q))) then     !  skipping down to energies where the excited states are open
              do i = 1, block_size; read(7,*); end do
              cycle
             else
               skip = skip+1
               if (skip >= skip_en) then
                 do i = 1, block_size; read(7,*); end do   !Skip first blocks because they have bad vals
                 read(7,*,end=1) a1, b, c, E 
                 !print *, "a, b, c, E :", a1, b, c, E, i_Q
                 energies(ie) = E
                 !print *, E   !Get the vals for the first real block
                 !print *, i_Q, E
               else 
                 do i = 1, block_size; read(7,*); end do
                 cycle
               end if
             end if
   
             !print *, a1, b, c, E
              

         !=====================================
         !- read from file and save in 1D array
         !=====================================
             do j = 1, c ! n_num
               read(7,'(D20.13)',advance="NO") temp_arr(j)
               if(mod(j,4)==0) Then
                 read(7,*) !- to go to the next line
               endif
               ! print *, "temp : ", j, temp_arr(j)
             enddo
          
         !=================================
         !- constructing K-matrix from temp
         !=================================
             my_count = 0
             ! print *, md(i_Q,i_name)
             do i = 1, a1  
               do j = 1, i
                 my_count = my_count + 1
                 k_mtrx(i,j) = temp_arr(my_count)
                 k_mtrx(j,i) = k_mtrx(i,j)
             ! write(*,'(2i4, 3d20.13)') i,j,temp_arr(my_count),real(k_mtrx(i,j))
             !If(i_name==1)Then
               !print * ,i, j, k_mtrx(i,j) , "--> test : OK"
             !EndIf
               enddo
             enddo
   
             !If(i_Q==1 .and. i_name==1) Then
               !print *, k_mtrx
             !Endif
            
   
             ! do i=1, a1
   
             !     !Crossing over ppoint means quantemol say ground is n=3 and exc are n=1,2
             !     if (i_Q == num_Q .or. i_Q == num_Q-1) then
             !       if (D2h_Dinf(i,i_name,1,i_Q) == 1) then
             !         D2h_Dinf(i,i_name,1,i_Q) = 11
             !       end if
             !       if (D2h_Dinf(i,i_name,1,i_Q) == 2) then
             !         D2h_Dinf(i,i_name,1,i_Q) = 22
             !       end if
             !       if (D2h_Dinf(i,i_name,1,i_Q) == 3) then
             !         D2h_Dinf(i,i_name,1,i_Q) = 33
             !       end if
         
             !       if (D2h_Dinf(i,i_name,1,i_Q) == 11) then
             !         D2h_Dinf(i,i_name,1,i_Q) = 2
             !       end if
             !       if (D2h_Dinf(i,i_name,1,i_Q) == 22) then
             !         D2h_Dinf(i,i_name,1,i_Q) = 3
             !       end if
             !       if (D2h_Dinf(i,i_name,1,i_Q) == 33) then
             !         D2h_Dinf(i,i_name,1,i_Q) = 1
             !       end if
   
             !     end if
   
             ! end do

             do i=1, a1
               do j=1, a1

                call lm2i(D2h_Dinf(i,i_name,1,i_Q), D2h_Dinf(i,i_name,2,i_Q), D2h_Dinf(i,i_name,3,i_Q),i2)
                call lm2i(D2h_Dinf(j,i_name,1,i_Q), D2h_Dinf(j,i_name,2,i_Q), D2h_Dinf(j,i_name,3,i_Q),j2)
   
                If(i2<1 .or. i2>dim_lm .or. j2<1 .or. j2>dim_lm )Stop 'Problem with indexes!'
                K_Q(i2,j2,i_Q) = k_mtrx(i,j)

                ! if (i_Q == 3) then
                !   print *, K_Q(i2,j2,i_Q,ie)
                ! end if

               enddo
             enddo
             Exit
          enddo energy
              
          1     deallocate( temp_arr, k_mtrx)
                close(7)

          enddo file_read
              
      Enddo geom_Q2



end subroutine

Subroutine i2lm(i,n,l,m)
  Use global_params
  Integer i,n,l,m,j
  j=0
  Do n = 1, n_max; Do m=-l_max,l_max; Do l=abs(m),l_max;
    j=j+1 
    !print *, n, m, l, j
    if(i==j)return
  Enddo;enddo;enddo
  Print *,'in Subroutine i2lm i is too large, i=',i;stop
End Subroutine i2lm

Subroutine lm2i(n,l,m,j)
  Use global_params
  Integer np,lp,mp,n,l,m,j
  j=0
  Do np=1,n; Do mp=-l_max,l_max;Do lp=abs(mp),l_max;
    j=j+1 
    !print *, np, mp, lp, j
    if(np==n .and. mp==m .and. lp==l)return
  Enddo;enddo;enddo
 
  Print *,'in Subroutine lm2i i,j, or m is out of range l,m,j,i=',l,m,j;stop
End Subroutine lm2i

subroutine linearize_Delta(Delta, lin_Delta)

  use global_params, only: dim_lm, num_Q, lin_start, lin_end, geoms, lin_geoms, num_Q_

  real(8), dimension(dim_lm, dim_lm, num_Q), intent(in)    ::  Delta
  real(8), dimension(dim_lm, dim_lm, num_Q_), intent(out)   ::  lin_Delta

  real(8), dimension(dim_lm, dim_lm)  ::  a
  integer ::  i, j, i_Q, i_Q_


  do i = 1, dim_lm
    do j = 1, dim_lm
      a(i,j) = ( Delta(i,j,lin_end) - Delta(i,j,lin_start) ) / ( geoms(lin_end) - geoms(lin_start) )
    end do
  end do

  do i = 1, dim_lm
    do j = 1, dim_lm
      do i_Q_ = 1, num_Q_
      lin_Delta(i,j,i_Q_) = Delta(i,j,lin_start) + ( a(i,j) * (lin_geoms(i_Q_) - geoms(lin_start)) )
      end do
    end do
  end do

end subroutine

subroutine find_rderiv(y_fin, y_in, dx, deriv)

	real(8), intent(in)	::	y_fin, y_in
	real(8), intent(in)	::	dx
	real(8), intent(out)	::	deriv

  deriv = (y_fin - y_in)/dx

end subroutine

subroutine K_to_S(md, D2h_Dinf, K_Q, S_Q)

  use global_params, only: num_Q, num_gamma_Q, dim_lm

  real(8), dimension(dim_lm, dim_lm, num_Q), intent(in)   ::  K_Q
  integer, dimension(num_Q, num_gamma_Q), intent(in)      ::  md
  integer, dimension(100, num_gamma_Q, 3, num_Q), intent(in)     ::  D2h_Dinf

  complex*16, dimension(dim_lm, dim_lm, num_Q), intent(out)  ::  S_Q

  integer ::  i_Q, i_name, i, j, a, ind, ind_p

  integer, allocatable  ::  nums(:,:)
  real(8), allocatable  ::  k_mtrx(:,:)
  complex*16, allocatable  :: s_mtrx(:,:)

    do i_Q = 1, num_Q
      do i_name = 1, num_gamma_Q

        !print *, i_Q, i_name
        allocate( nums(md(i_Q,i_name),3), k_mtrx(1:md(i_Q,i_name),1:md(i_Q,i_name)), s_mtrx(1:md(i_Q,i_name),1:md(i_Q,i_name)) )

        do a = 1, md(i_Q, i_name)
          nums(a,:) = D2h_Dinf(a, i_name, :, i_Q)
        end do

        do i = 1, md(i_Q,i_name)
          do j = 1, md(i_Q,i_name)

            call lm2i(nums(i,1), nums(i,2), nums(i,3), ind)
            call lm2i(nums(j,1), nums(j,2), nums(j,3), ind_p)

            !print *, nums(i,:)
            !print *, nums(j,:)
            k_mtrx(i,j) = K_Q(ind, ind_p, i_Q)
            !print *, ind, ind_p, ie

          end do
        end do


        s_mtrx=(0.d0,0.d0)
        ! if (i_Q == 10 .and. (i_name == 5 .or. i_name == 4)) then
        !   print *, k_mtrx
        ! end if
        call K_S(k_mtrx(1:md(i_Q,i_name),1:md(i_Q,i_name)), s_mtrx(1:md(i_Q,i_name),1:md(i_Q,i_name)), md(i_Q,i_name))
        
        do i = 1, md(i_Q,i_name)
          do j = 1, md(i_Q,i_name)

            call lm2i(nums(i,1), nums(i,2), nums(i,3), ind)
            call lm2i(nums(j,1), nums(j,2), nums(j,3), ind_p)

            S_Q(ind, ind_p, i_Q) = s_mtrx(i,j)

          end do
        end do
        
        deallocate(k_mtrx, s_mtrx, nums)

      end do
    end do



end subroutine

Subroutine K_S(K,S,N)
  Integer i,INFO,N,IPIV(N)
  Complex*16 S(N,N),K2(N,N),K1(N,N),ci
  Real*8 K(N,N)
  ci=(0.d0,1.d0)
  K2=ci*K
  K1=-1.d0*K2
  Do i=1,N
    K2(i,i)=K2(i,i)+1.d0
    K1(i,i)=K1(i,i)+1.d0
  Enddo
!  Inversion and Multiplication
  K1=Transpose(K1); K2=Transpose(K2)
  CALL ZGESV(N,N,K1,N,IPIV,K2,N,INFO)
  S=Transpose(K2)
  !print *, S
End Subroutine K_S


Subroutine K_delta(Kmat,delta_mat,eigenvalues)

  use global_params, only: dim_lm

  Integer i,j,k
  Real*8 Kmat(dim_lm,dim_lm),delta_mat(dim_lm,dim_lm),eigenvectors(dim_lm,dim_lm),eigenvalues(dim_lm)

  call diagonalization(dim_lm,Kmat,eigenvalues,eigenvectors)
  eigenvalues=atan(eigenvalues)
  delta_mat=0.d0

  Do i=1,dim_lm;Do j=1,dim_lm
    Do k=1,dim_lm
      delta_mat(i,j)=delta_mat(i,j)+eigenvectors(i,k)*eigenvalues(k)*eigenvectors(j,k)
    Enddo
  Enddo; Enddo

End Subroutine K_delta

subroutine Delta_S(Delta, S)

  use global_params, only: dim_lm, ci


  real(8), dimension(dim_lm, dim_lm), intent(in)        :: Delta
  complex*16, dimension(dim_lm, dim_lm), intent(out)    :: S

  real(8), dimension(dim_lm, dim_lm)                    :: eigenvectors
  real(8), dimension(dim_lm)                            :: eigenvalues
  integer                                               :: i,j,k

  call diagonalization(dim_lm,Delta,eigenvalues,eigenvectors)

  eigenvalues = exp(2*ci*eigenvalues)

  Do i=1,dim_lm
    Do j=1,dim_lm
      Do k=1,dim_lm
        S(i,j)=S(i,j)+eigenvectors(i,k)*eigenvalues(k)*eigenvectors(j,k)
      Enddo
    Enddo
  Enddo

end subroutine


Subroutine diagonalization(N,A,eigenvalues,eigenvectors)
  Integer N,info
  Real*8 A(N,N),eigenvalues(N),eigenvectors(N,N),work(3*N)
  eigenvectors=A
  Call dsyev('V','U',N,eigenvectors,N,eigenvalues,work,3*N,info)
End Subroutine diagonalization


function cex2(nj)
  character*2 cex2
  integer cnj,dnj,unj,nj
  cnj = nj/100
  dnj = (nj - cnj*100)/10
  unj = nj - cnj*100 - dnj*10
!  cex2 = char(cnj+48)//char(dnj+48)//char(unj+48)
  cex2 = char(dnj+48)//char(unj+48)
end function cex2

!    Standart cubic spline procedure
SUBROUTINE SPLINE(N,X,Y,CM)
  Integer i,n
  Real*8 X(1:n),Y(1:n),CM(1:n),A,C,ALPHA(1:n),BETA(1:n),GAMMA(1:n),B(1:n)
  CM(1)=0.
  CM(N)=0.
  DO I=3,N
    A=X(I-1)-X(I-2)
    C=X(I)-X(I-1)
    ALPHA(I-2)=(A+C)/3.
    BETA(I-2)=C/6.
    GAMMA(I-2)=BETA(I-2)
    B(I-2)=(Y(I)-Y(I-1))/C-(Y(I-1)-Y(I-2))/A
  Enddo
  CALL TRIDIA(ALPHA,BETA,GAMMA,B,CM(2:n),N-2)
END

SUBROUTINE TRIDIA(ALPHA,BETA,GAMMA,B,X,N)
  Integer i,n,j
  Real*8 ALPHA(1:n+2),BETA(1:n+2),GAMMA(1:n+2),B(1:n+2),X(1:n+1),RAP
  DO I=2,N
    RAP=BETA(I-1)/ALPHA(I-1)
    ALPHA(I)=ALPHA(I)-RAP*GAMMA(I-1)
    B(I)=B(I)-RAP*B(I-1)
  Enddo
  X(N)=B(N)/ALPHA(N)
  DO J=2,N
    I=N-J+1
    X(I)=(B(I)-GAMMA(I)*X(I+1))/ALPHA(I)
  Enddo
END

real(8) FUNCTION SPL(N,X,Y,M,T)
  Integer i,k,n
  Real*8 X(1:n),Y(1:n),M(1:n),G,E,F,T
  IF(T.LE.X(1)) GO TO 30
    IF(T.GE.X(N)) GO TO 40
  K=2
10   IF(T.LE.X(K)) GO TO 20
  K=K+1
  GO TO 10
20   E=X(K)-X(K-1)
  F=X(K)-T
  G=T-X(K-1)
  SPL=(M(K-1)*F*F*F+M(K)*G*G*G+(6.*Y(K)-M(K)*E*E)*G+(6.*Y(K-1)- M(K-1)*E*E)*F)/(6.*E)
  RETURN
30   E=X(2)-X(1)
  SPL=((Y(2)-Y(1))/E-M(2)*E/6.)*(T-X(1))+Y(1)
  RETURN
40   E=X(N)-X(N-1)
  SPL=((Y(N)-Y(N-1))/E+M(N-1)*E/6.)*(T-X(N))+Y(N)
END