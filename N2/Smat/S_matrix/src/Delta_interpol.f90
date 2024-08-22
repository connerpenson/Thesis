! Module global_params
!   complex(kind=8), parameter :: ci = (0.d0, 1.d0)
!   integer, parameter :: num_gamma_Q = 8,l_max = 4, n_max = 3, vib_num = 9, num_Q_ = 10, neut_vib_num = 3, new_param = 3
!   real*8 :: ryd,au_eV,au_cm,au_to_sec
!   real(8) :: pi,ev_au
!   integer :: num_Q, num_v ,num_ev,num_scatE,neut_vib
!   integer, allocatable, dimension(:,:) ::  qn_ev
!   real(kind=8), allocatable, dimension(:) ::  geom,energies_ev,E_el_Rmat,energies_v
!   Complex*16, allocatable, dimension(:,:,:,:) :: S_Q,S_Q_check
!   real(8), parameter, dimension(num_Q_) ::  geom_ = (/1.9, 1.944, 1.989, 2.0333, 2.0778, 2.122, 2.1667, 2.211, 2.255, 2.3/) ! in atomic units
!   integer, parameter :: dim_lm = (l_max+1)**2 * n_max
!   real(8), allocatable, dimension(:,:) :: chan_ens
!   integer, parameter  ::  num_tot = dim_lm*vib_num
  

! End Module global_params
!---------------------------------------------------
program extrapolate_SQ
  Use global_params
  implicit none

  real(8) ::  omega_0
  
  Ryd=219474.6313710d0
  au_eV=27.2113834d0
  au_cm=0.5291772083d-8
  au_to_sec=2.4188843d-17 
  pi=atan(1.d0)*4.d0
  ev_au = .036749405469679d0


  num_Q=10 ! number of calculated geometries
  
  call Get_Extrapolated_SQ_matrix
  
end program extrapolate_SQ

Subroutine Get_Extrapolated_SQ_matrix
    Use global_params
    Use stdlib_sorting, only: sort_index
    Use iso_fortran_env, only: int64, iostat_end
    implicit none
    integer lwork,info,ie,i,j,i_Q,i_Q_,j2,i2,n,l,m,n_p,l_p,m_p,&
        num_ener(100),ordered_geom(100), a1, b, c, i_name, k, my_count, n_num,&
        itmp(1:4), dims_extended(num_gamma_Q), i_rec, records(5,2000), num_recs, file_num,&
        i_QQ, n_prime, l_prime, m_prime, ie_p,v,vib_i,v_p,e_i, error
    real*8, allocatable, dimension(:,:,:) :: eigenphases, a
    real*8, allocatable, dimension(:) :: tmp_eigenphases, rwork
    real*8, allocatable, dimension(:,:) :: tmp_Delta, tmp_Delta_, U, k_mtrx, eigenphases_matrix, M_print
    real(kind=8) :: E, tmp, ener_ch_lm(2000), old_deriv, deriv_1, deriv_2, temp_en, en, en_p
    complex(kind=8), allocatable, dimension(:,:) :: s_mtrx, s_mtrx_
    real(kind=8), allocatable, dimension(:,:,:,:) :: K_Q, Delta_Q
    complex(kind=8), allocatable, dimension(:,:,:,:) :: Delta_Q_interpol, S_Q_interpol, Delta_Q_extrapol, S_Q_extrapol
    complex*16, allocatable, dimension(:,:) :: tmp_eigenphases_matrix, tmp_eigenphases_matrix_
    integer iQ_eq, block_size, skip_en
    real(kind=8), allocatable, dimension(:) :: temp_arr
    real(kind=8), allocatable, dimension(:,:) :: thresh_ens
    integer, allocatable, dimension(:,:,:,:) :: D2h_Dinf
    integer, allocatable, dimension(:,:) :: md
    character(len=15), dimension(num_gamma_Q) :: name_of_file
    character*2 cex2
    integer skip, eof, neut_iv
    real(8), allocatable  ::  ion_WF_array(:,:,:)
    real(8), allocatable  ::  neut_WF_array(:,:)
    complex(8), allocatable, dimension(:,:) :: S_ft, sorted_S_ft
    complex(8), allocatable, dimension(:,:,:) :: D_ft, sorted_D_ft
    integer, allocatable, dimension(:,:) :: ft_channels, sorted_ft_channels
    complex(8), allocatable, dimension(:,:,:) :: dipole_data, interp_dipole_data
    real(8), dimension(100) :: temp
    real(8), allocatable, dimension(:)  ::  energy_matrix, Z_norms
    integer, dimension(4)  :: num_tmp
    integer(int64), dimension(num_tot) ::  index
    real(8),  allocatable, dimension(:,:,:) ::  TCS_array, SPES_array, conv_TCS_array
    integer,    allocatable		::	IPIV(:), max_chans(:,:)
    complex(8), allocatable		::	dipole_phys(:,:), dipole_open(:,:), dipole_closed(:,:), smat_cc(:,:), smat_co(:,:)
    complex(8), allocatable		::	beta(:)
    real(8)      ::   energy_begin, energy_end, energy_step, omega, diff, energy, const, elec_en_range, alpha, omega_0, en_shift, r
    integer :: num_closed, num_open
    real(8), dimension(0:2)    ::  delta_Es 
    logical ::  trim_cond = .true., lin_approx = .true.
    real(8), allocatable  :: nu(:,:,:)
    real(8) ::  trim_thresh = 1
 
    omega_0 = 15.577501268901065d0*ev_au
    alpha = 0.0072973525693d0
    const = (4d0/3d0)*(pi**2)*alpha
    elec_en_range = .1*ev_au
    
    delta_Es = (/.0000, .2873222056, .5614759886/)*ev_au 
    
    iQ_eq = 3
    skip_en = 5
    
    ! print *, dim_lm
  
    num_scatE = 1000

    num_ener=num_scatE ! number of energies for each i_Q 
    allocate(E_el_Rmat(num_scatE),S_Q(dim_lm,dim_lm,num_Q,num_scatE),geom(num_Q),D2h_Dinf(1000,num_gamma_Q,3,num_Q))
    allocate(md(num_Q,num_gamma_Q),S_Q_check(dim_lm,dim_lm,num_Q,num_scatE))
    allocate(K_Q(dim_lm,dim_lm,num_Q,num_scatE), a(dim_lm,dim_lm,num_scatE))
    allocate(U(dim_lm,dim_lm), M_print(dim_lm,dim_lm))
    allocate(tmp_eigenphases(1:dim_lm))
    allocate(eigenphases(dim_lm,num_Q,num_scatE),Delta_Q(dim_lm,dim_lm,num_Q,num_scatE),&
    tmp_Delta(dim_lm,dim_lm),tmp_Delta_(dim_lm,dim_lm))
    allocate(eigenphases_matrix(dim_lm,dim_lm))
    allocate(Delta_Q_interpol(dim_lm, dim_lm, num_Q, num_scatE),S_Q_interpol(dim_lm, dim_lm, num_Q, num_scatE),&
    Delta_Q_extrapol(dim_lm, dim_lm, num_Q_, num_scatE),S_Q_extrapol(dim_lm, dim_lm, num_Q_, num_scatE))
    allocate(tmp_eigenphases_matrix(dim_lm,dim_lm), tmp_eigenphases_matrix_(dim_lm,dim_lm))
    allocate(thresh_ens(100, num_Q)) ! Hard coded 100 says we'll never have an n higher than 100

! Getting the calculated geometries and ordering them
  
    Do i_Q=1,num_Q
      open(unit = 7, file='K-mats/'//cex2(i_Q)//'/Q.txt', action='read')
      read(7,*)geom(i_Q)
      ! print *,'geometries are',geom(i_Q)
    Enddo
  
    !  ordering geometries
    Do i=1,num_Q
      ordered_geom(i)=i
    Enddo
    Do i=1,num_Q
      Do j=i+1,num_Q
        If(geom(i)>geom(j))Then
          tmp=geom(i)
          geom(i)=geom(j)
          geom(j)=tmp
          k=ordered_geom(i)
          ordered_geom(i)=ordered_geom(j)
          ordered_geom(j)=k
        Endif
      Enddo
    Enddo

! For each channel, establish its dimensions and read the kmat file to write it into the global matrix K_Q (until ~ line 221)

 geom_Q: Do i_Q=1,num_Q
 
 ! print *, i_Q
 
    name_of_file(1:num_gamma_Q) = (/ '1Ag.channel ', '1Au.channel ', '1B1g.channel', '1B1u.channel', &
                                        '1B2g.channel', '1B2u.channel', '1B3g.channel', '1B3u.channel'/)!, &
                                                      !'3Ag.channel ', '3Au.channel ', '3B1g.channel', '3B1u.channel', &
                                                                  !'3B2g.channel', '3B2u.channel', '3B3g.channel', '3B3u.channel' /)
    i_rec=0
    do i_name = 1, num_gamma_Q
      ! print *, name_of_file(i_name)   
     
      
      open(unit = 7,&
      file= 'K-mats/'//cex2(ordered_geom(i_Q))//'/matrices/'//name_of_file(i_name),&
      action='read')
    

      do j = 1, 2;  read(7,*) ; enddo !- ignore useless lines
      read(7,*)itmp(1:4)
      dims_extended(i_name)=itmp(4)
     
      !D_infty,v symmetric str.
      ! md(1:8)=(/6, 1, 3, 3, 3, 3, 3, 3 /) 


      if(i_name==1)md(i_Q,i_name)=12
      if(i_name==2)md(i_Q,i_name)=7
      if(i_name==3)md(i_Q,i_name)=9
      if(i_name==4)md(i_Q,i_name)=9
      if(i_name==5)md(i_Q,i_name)=7
      if(i_name==6)md(i_Q,i_name)=12
      if(i_name==7)md(i_Q,i_name)=7
      if(i_name==8)md(i_Q,i_name)=12
      ! if(i_name==9)md(i_Q,i_name)=12
      ! if(i_name==10)md(i_Q,i_name)=7
      ! if(i_name==11)md(i_Q,i_name)=9
      ! if(i_name==12)md(i_Q,i_name)=9
      ! if(i_name==13)md(i_Q,i_name)=7
      ! if(i_name==14)md(i_Q,i_name)=12
      ! if(i_name==15)md(i_Q,i_name)=7
      ! if(i_name==16)md(i_Q,i_name)=12

      do j = 1, itmp(1)+1;  read(7,*) ; enddo !- ignore useless lines

      do j = 1, dims_extended(i_name)
        i_rec=i_rec+1
        records(1,i_rec)=i_name
        read(7,*)records(2,i_rec),records(3,i_rec),records(4,i_rec),records(5,i_rec),ener_ch_lm(i_rec)
        !Write(6,*)'iQ=',i_Q,records(1:5,i_rec),ener_ch_lm(i_rec)
        thresh_ens(records(3,i_rec), i_Q) = ener_ch_lm(i_rec) ! Getting the threshold energies for each n, geom
        D2h_Dinf(records(2,i_rec),i_name,1:3,i_Q) = records(3:5,i_rec) ! i_el, l,lambda
      enddo
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
      allocate( temp_arr(n_num), k_mtrx(md(i_Q,i_name),md(i_Q,i_name)))
      allocate(s_mtrx(md(i_Q,i_name),md(i_Q,i_name)), s_mtrx_(md(i_Q,i_name),md(i_Q,i_name)))
     
      ! print *, name_of_file(i_name)
   
      open(unit = 7, file= 'K-mats/'//cex2(ordered_geom(i_Q))//'/matrices/'//name_of_file(i_name),&
 action='read')
      do j = 1, 5;  read(7,*) ; enddo !- ignore useless lines
      !- loop over energy
    
      skip = 0
      
      energy : Do ie = 1, num_scatE
    
          read(7,*,end=1) a1, b, c, E
          ! print *, "a, b, c, E :", a1, b, c, E ! /2.d0 * au_eV

          if (mod(c,4)==0) then
            block_size = c/4
          else
            block_size = ceiling(real(c/4.))
          end if

          !print *, E, thresh_ens(n_max, i_Q)
          
          if (E < thresh_ens(n_max, i_Q)) then     !  skipping down to energies where the excited states are open
            do i = 1, block_size; read(7,*); end do
            cycle
          else 
            skip = skip+1
            if (skip == skip_en) then
              do i = 1, block_size; read(7,*); end do   !Skip first block because it has bad vals
              read(7,*,end=1) a1, b, c, E 
              !print *, E   !Get the vals for the first real block
              ie_p = ie
              print *, i_Q, E
            else 
              do i = 1, block_size; read(7,*); end do
              cycle
            end if
          end if

          !print *, a1, b, c, E
	    
	  !If(a1.ne.md(i_Q,i_name)) exit
          If(i_Q==1 .and. i_name==1) Then
            E_el_Rmat(ie)=E/2.d0
          Endif
	       
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
          my_count =0
          ! print *, md(i_Q,i_name)
          do i=1, a1  
            do j=1, i
              my_count = my_count + 1
              k_mtrx(i,j) = temp_arr(my_count)
              k_mtrx(j,i) = k_mtrx(i,j)
	      ! write(*,'(2i4, 3d20.13)') i,j,temp_arr(my_count),real(k_mtrx(i,j))
	      !If(i_name==1)Then
	        !print*,i,j,real(k_mtrx(i,j)),'--> test : OK'
	      !EndIf
            enddo
          enddo

          !If(i_Q==1 .and. i_name==1) Then
            !print *, k_mtrx
          !Endif
          
          
          s_mtrx=(0.d0,0.d0)
          call K_S(k_mtrx(1:md(i_Q,i_name),1:md(i_Q,i_name)),s_mtrx(1:md(i_Q,i_name),1:md(i_Q,i_name)),md(i_Q,i_name))
          

        if (ie == ie_p) then

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
              K_Q(i2,j2,i_Q,1) = k_mtrx(i,j)
		          S_Q(i2,j2,i_Q,1) = s_mtrx(i,j)

              !Need to acconut for possible flipping of degenerate states
              if (D2h_Dinf(i,i_name,1,i_Q) == 2 .or. D2h_Dinf(j,i_name,1,i_Q) == 2) then
                if (D2h_Dinf(i,i_name,1,i_Q) == 2) then
                  call lm2i(3, D2h_Dinf(i,i_name,2,i_Q), D2h_Dinf(i,i_name,3,i_Q),i2)
                end if

                if (D2h_Dinf(j,i_name,1,i_Q) == 2) then
                  call lm2i(3, D2h_Dinf(j,i_name,2,i_Q), D2h_Dinf(j,i_name,3,i_Q),j2)
                end if
                K_Q(i2,j2,i_Q,1) = k_mtrx(i,j)
                S_Q(i2,j2,i_Q,1) = s_mtrx(i,j)
              end if

              if (D2h_Dinf(i,i_name,1,i_Q) == 3 .or. D2h_Dinf(j,i_name,1,i_Q) == 3) then
                if (D2h_Dinf(i,i_name,1,i_Q) == 3) then
                  call lm2i(2, D2h_Dinf(i,i_name,2,i_Q), D2h_Dinf(i,i_name,3,i_Q),i2)
                end if

                if (D2h_Dinf(j,i_name,1,i_Q) == 3) then
                  call lm2i(2, D2h_Dinf(j,i_name,2,i_Q), D2h_Dinf(j,i_name,3,i_Q),j2)
                end if
                K_Q(i2,j2,i_Q,1) = k_mtrx(i,j)
                S_Q(i2,j2,i_Q,1) = s_mtrx(i,j)
              end if

            enddo
          enddo
        end if
 	  
      enddo energy
      
1     num_ener(i_Q)=min(num_ener(i_Q),ie-1)
      ! Print *,'num_ener=',num_ener(i_Q)
      deallocate( temp_arr, k_mtrx,s_mtrx,s_mtrx_)
      close(7)
    enddo file_read
    
    Enddo geom_Q


    open(333, file = "smat.dat")
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        !write(333, *) "Partial Wave:", n, l, m, n_prime, l_prime, m_prime
        do i_QQ = 1, num_Q
          write(333, *) geom_(i_QQ), real(S_Q(i,j,i_QQ,1))
        end do
        write(333, *)
      end do
    end do
    close(333)

    open(333, file = "kmat.dat")
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        !write(333, *) "Partial Wave:", n, l, m, n_prime, l_prime, m_prime
        do i_QQ = 1, num_Q
          write(333, *) geom_(i_QQ), K_Q(i,j,i_QQ,1)
        end do
        write(333, *)
      end do
    end do
    close(333)

    old_deriv = 0
    do i = 1, dim_lm
      do j = 1, dim_lm
        do i_Q = 1, num_Q
          if (i_Q > 1) then

            call find_rderiv(real(S_Q(i, j, i_Q, 1)), real(S_Q(i, j, i_Q - 1, 1)), geom(i_Q) - geom(i_Q-1), deriv_1)
              call find_rderiv(-real(S_Q(i, j, i_Q, 1)), real(S_Q(i, j, i_Q - 1, 1)), geom(i_Q) - geom(i_Q-1), deriv_2)

                if (i_Q == 2) then
                  if (abs(deriv_1) <= abs(deriv_2)) then
                    cycle
                  else
                    S_Q(i, j, i_Q, 1) = -S_Q(i, j, i_Q, 1)
                  end if
                else
                  if (abs(deriv_1 - old_deriv) <= abs(deriv_2 - old_deriv)) then
                    cycle
                  else
                    S_Q(i, j, i_Q, 1) = -S_Q(i, j, i_Q, 1)
                  end if
                end if

                call find_rderiv(real(S_Q(i, j, i_Q, 1)), real(S_Q(i, j, i_Q - 1, 1)), geom(i_Q) - geom(i_Q-1), old_deriv)

            else
              cycle
            end if
        end do
      end do
    end do

    
    open(333, file = "smat_fixed.dat")
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        !write(333, *) "Partial Wave:", n, l, m, n_prime, l_prime, m_prime
        do i_QQ = 1, num_Q
          write(333, *) geom_(i_QQ), real(S_Q(i,j,i_QQ,1))
        end do
        write(333, *)
      end do
    end do
    close(333)

      

! Get Delta from the K matrix and calculate the eigenphases matrix
    
    num_scatE=minval(num_ener(1:num_Q))
 
!   get eigenphases matrix
!   for num_Q geometries    
    lwork = 3*dim_lm
    !do ie=1,num_scatE
        do i_Q=1,num_Q
           allocate(rwork(lwork))
           rwork=0.d0
           U=0.d0
           tmp_eigenphases=0.d0
           do i=1, dim_lm
             do j=1, dim_lm
               U(i,j) = K_Q(i,j,i_Q,1)
             enddo
           enddo
          ! diagonalizing U
           Call dsyev('V','U',dim_lm,U(:,:),dim_lm,tmp_eigenphases(1:dim_lm),rwork,lwork,info)
           deallocate(rwork)
           eigenphases(1:dim_lm,i_Q,1)=atan(tmp_eigenphases(1:dim_lm))
           eigenphases_matrix=0.d0
           do i=1,dim_lm
              eigenphases_matrix(i,i)=eigenphases(i,i_Q,1)
           enddo
          ! calculating Delta
           Delta_Q(:,:,i_Q,1)=matmul(U(:,:),matmul(eigenphases_matrix,transpose(U(:,:))))
          !  if(i_Q==1) then
          !    if(ie==num_scatE/2) then
          !      M_print = matmul(U(:,:),transpose(U(:,:)))
          !      do i = 1,dim_lm
          !        print *, M_print(i,i)
          !      enddo
          !    endif
          !  endif
           

          tmp_Delta=0.d0
          tmp_Delta(1:dim_lm,1:dim_lm)=Delta_Q(1:dim_lm,1:dim_lm,i_Q,1)

          allocate(rwork(lwork))
          rwork=0.d0
          tmp_eigenphases=0.d0
          Call dsyev('V','U',dim_lm,tmp_Delta(:,:),dim_lm,tmp_eigenphases(1:dim_lm),rwork,lwork,info)
          deallocate(rwork)
          
          tmp_eigenphases_matrix=(0.d0,0.d0)
          do i=1,dim_lm
            tmp_eigenphases_matrix(i,i)=exp(2*ci*tmp_eigenphases(i))
            ! print *, 'ok'
          enddo
          
          s_mtrx=matmul(tmp_Delta,matmul(tmp_eigenphases_matrix,transpose(tmp_Delta)))
          S_Q_check(1:dim_lm,1:dim_lm,i_Q,1)=s_mtrx(1:dim_lm,1:dim_lm)
        enddo
    !enddo

    open(333, file = "delta_unfixed.dat")
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        !write(333, *) "Partial Wave:", n, l, m, n_prime, l_prime, m_prime
        do i_QQ = 1, num_Q
          write(333, *) geom_(i_QQ), Delta_Q(i,j,i_QQ,1)
        end do
        write(333, *)
      end do
    end do
    close(333)


    old_deriv = 0
    do i = 1, dim_lm
      do j = 1, dim_lm
        do i_Q = 1, num_Q
          if (i_Q > 1) then

            call find_rderiv(Delta_Q(i, j, i_Q, 1), Delta_Q(i, j, i_Q - 1, 1), geom(i_Q) - geom(i_Q-1), deriv_1)
              call find_rderiv(-Delta_Q(i, j, i_Q, 1), Delta_Q(i, j, i_Q - 1, 1), geom(i_Q) - geom(i_Q-1), deriv_2)

                if (i_Q == 2) then
                  if (abs(deriv_1) <= abs(deriv_2)) then
                    cycle
                  else
                    Delta_Q(i, j, i_Q, 1) = -Delta_Q(i, j, i_Q, 1)
                  end if
                else
                  if (abs(deriv_1 - old_deriv) <= abs(deriv_2 - old_deriv)) then
                    cycle
                  else
                    Delta_Q(i, j, i_Q, 1) = -Delta_Q(i, j, i_Q, 1)
                  end if
                end if

                call find_rderiv(Delta_Q(i, j, i_Q, 1), Delta_Q(i, j, i_Q - 1, 1), geom(i_Q) - geom(i_Q-1), old_deriv)

            else
              cycle
            end if
        end do
      end do
    end do


    open(333, file = "delta_fixed.dat")
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        !write(333, *) "Partial Wave:", n, l, m, n_prime, l_prime, m_prime
        do i_QQ = 1, num_Q
          write(333, *) geom_(i_QQ), Delta_Q(i,j,i_QQ,1)
        end do
        write(333, *)
      end do
    end do
    close(333)
    
    
    print *, 'Finished calculating eigenphases matrix'

! Get the linear coefficients for the Delta elements


       do i=1,dim_lm
          do j=1,dim_lm
             a(i,j,1)=(Delta_Q(i,j,5,1)-Delta_Q(i,j,1,1))/(geom(5)-geom(1))
          enddo
       enddo


!   Phase matrix extrapolation and build S-matrix

    print*,'phase matrix extrapolation and S-matrix building'
    do i_Q=1,num_Q

          tmp_Delta=(0.d0,0.d0)
          tmp_Delta(1:dim_lm,1:dim_lm)=Delta_Q(1:dim_lm,1:dim_lm,5,1)+a(1:dim_lm,1:dim_lm,1)*(geom(i_Q)-geom(5))
 

          Delta_Q_interpol(1:dim_lm,1:dim_lm,i_Q,1)=tmp_Delta(1:dim_lm,1:dim_lm)
 

          allocate(rwork(lwork))
          rwork=0.d0;tmp_eigenphases=0.d0
          Call dsyev('V','U',dim_lm,tmp_Delta,dim_lm,tmp_eigenphases(1:dim_lm),rwork,lwork,info)
          deallocate(rwork)
          tmp_eigenphases_matrix=(0.d0,0.d0)
          do i=1,dim_lm
            tmp_eigenphases_matrix(i,i)=exp(2*ci*tmp_eigenphases(i))
            ! print *, 'ok'
          enddo
          s_mtrx=matmul(tmp_Delta,matmul(tmp_eigenphases_matrix,transpose(tmp_Delta)))
          S_Q_interpol(1:dim_lm,1:dim_lm,i_Q,1)=s_mtrx(1:dim_lm,1:dim_lm)

    enddo

    open(333, file = "interp_smat_real.dat")
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        !write(333, *) "Partial Wave:", n, l, m, n_prime, l_prime, m_prime
        do i_QQ = 1, num_Q
          write(333, *) geom_(i_QQ), real(S_Q_interpol(i,j,i_QQ,1))
        end do
        write(333, *)
      end do
    end do
    close(333)

    open(333, file = "interp_smat_imag.dat")
    do i=1, dim_lm
      do j=1, dim_lm
        call i2lm(i, n, l, m)
        call i2lm(j, n_prime, l_prime, m_prime)
        !write(333, *) "Partial Wave:", n, l, m, n_prime, l_prime, m_prime
        do i_QQ = 1, num_Q
          write(333, *) geom_(i_QQ), aimag(S_Q_interpol(i,j,i_QQ,1))
        end do
        write(333, *)
      end do
    end do
    close(333)
    
!   Phase matrix extrapolation and build S-matrix for new geometries


    print*,'phase matrix extrapolation and S-matrix building for new geometries'
    do i_Q_=1,num_Q_
       
          
          tmp_Delta_=(0.d0,0.d0)
          do i=1,dim_lm
            do j=1,dim_lm
              tmp_Delta_(i,j)=Delta_Q(i,j,5,1)+a(i,j,1)*(geom_(i_Q_)-geom(5))
            enddo
          enddo
          !if (i_Q_ == 1 .AND. ie == num_scatE/2) then
            !print *, 'tmp',tmp_Delta_(7,1:dim_lm), 'Delta',Delta(7,1:dim_lm,3,ie), 'a',a(7,1:dim_lm,ie), 'diff',geom_(i_Q_)-geom(3)
          !endif
          Delta_Q_extrapol(1:dim_lm,1:dim_lm,i_Q_,1)=tmp_Delta_(1:dim_lm,1:dim_lm)

 

          allocate(rwork(lwork))
          rwork=0.d0;tmp_eigenphases=0.d0
          Call dsyev('V','U',dim_lm,tmp_Delta_,dim_lm,tmp_eigenphases,rwork,lwork,info)
          deallocate(rwork)
          tmp_eigenphases_matrix_=(0.d0,0.d0)
          do i=1,dim_lm
            tmp_eigenphases_matrix_(i,i)=exp(2*ci*tmp_eigenphases(i))
          enddo
          s_mtrx_=matmul(matmul(tmp_Delta_,tmp_eigenphases_matrix_),transpose(tmp_Delta_))
          !if (i_Q_ == 1 .AND. ie == num_scatE/2) then
            !print *, s_mtrx_(7,1:dim_lm)
          !endif
          S_Q_extrapol(1:dim_lm,1:dim_lm,i_Q_,1)=s_mtrx_(1:dim_lm,1:dim_lm)
       
    enddo


! Read in the dipole elements from Q-mol
    print *, "reading in the dipole elements"

    allocate(dipole_data(3,dim_lm,num_Q))

    open(27, file = "dipoles.dat")
    do j = 1, dim_lm
      read(27,*) n,l,m
      call lm2i(n,l,m,i)
      do i_Q_ = 1, num_Q
        read(27,*) temp(1:6)
        dipole_data(1,i,i_Q_) = complex(temp(1),temp(2))
        dipole_data(2,i,i_Q_) = complex(temp(3),temp(4))
        dipole_data(3,i,i_Q_) = complex(temp(5),temp(6))
      end do
      read(27,*)
    end do
    close(27)

    print *, "reading in the linear dipole elements"

    allocate(interp_dipole_data(3,dim_lm,num_Q))

    open(27, file = "interp_dipoles.dat")
    do j = 1, dim_lm
      read(27,*) n,l,m
      call lm2i(n,l,m,i)
      do i_Q_ = 1, num_Q
        read(27,*) temp(1:6)
        interp_dipole_data(1,i,i_Q_) = complex(temp(1),temp(2))
        interp_dipole_data(2,i,i_Q_) = complex(temp(3),temp(4))
        interp_dipole_data(3,i,i_Q_) = complex(temp(5),temp(6))
      end do
      read(27,*)
    end do
    close(27)

    if (trim_cond .eqv. .true.) then
      do j = 1, dim_lm
        if (sum(abs(real(dipole_data(:,j,:)))) < trim_thresh .and. sum(abs(aimag(dipole_data(:,j,:)))) < trim_thresh) then
          dipole_data(:,j,:) = 0
        end if
        if (sum(abs(real(interp_dipole_data(:,j,:)))) < trim_thresh .and. sum(abs(aimag(interp_dipole_data(:,j,:)))) < trim_thresh) then
          interp_dipole_data(:,j,:) = 0
        end if
      end do
    end if

    Open(27,file='dips_check.dat')
  do j = 1, 3
    do i = 1,dim_lm
        write(27,*)
        call i2lm(i,n,l,m)
        write(27,*) n,l,m
        do i_Q = 1, num_Q
          Write(27,*) geom_(i_Q), real(dipole_data(j,i,i_Q))
        end do
    end do
  end do
  do j = 1, 3
    do i = 1,dim_lm
        write(27,*)
        call i2lm(i,n,l,m)
        !write(27,*) n,l,m
        do i_Q = 1, num_Q
          Write(27,*) geom_(i_Q), aimag(dipole_data(j,i,i_Q))
        end do
    end do
  end do
    Close(27)

    Open(27,file='lin_dips_check.dat')
    do j = 1, 3
      do i = 1,dim_lm
          write(27,*)
          call i2lm(i,n,l,m)
          !write(27,*) n,l,m
          do i_Q = 1, num_Q
            Write(27,*) geom_(i_Q), real(interp_dipole_data(j,i,i_Q))
          end do
      end do
    end do
    do j = 1, 3
      do i = 1,dim_lm
          write(27,*)
          call i2lm(i,n,l,m)
          !write(27,*) n,l,m
          do i_Q = 1, num_Q
            Write(27,*) geom_(i_Q), aimag(interp_dipole_data(j,i,i_Q))
          end do
      end do
    end do
      Close(27)



    print*, "reading in channel energies"

    allocate(chan_ens(n_max, 0:vib_num-1), energy_matrix(num_tot))
    do neut_iv = 0, neut_vib_num-1
      open(27, file = "vibronic_channel_energies.txt")
      do 
          read(27,*, iostat = error) n, v, temp_en
          Select Case(error)
          case(0)
            if (v < vib_num) then
              if (n == 2 .or. n == 3) then
                chan_ens(n,v) = temp_en * ev_au - delta_Es(neut_iv) + ion_shift * ev_au
              else
                chan_ens(n,v) = temp_en * ev_au - delta_Es(neut_iv)
              end if
            end if
          case(iostat_end)
            exit
          case Default
            print *, "Error in reading channel energies"
            Stop
          end Select
      end do
      close(27)
    end do

    print *, "correcting channel energies"
    en_shift = chan_ens(1,0) - omega_0
    chan_ens = chan_ens - en_shift


! Frame transforming the S-matrix
    print*,'frame transforming the extrapolated S-matrix and dipole vectors'

    allocate(ion_WF_array(num_Q, 0:vib_num-1, n_max), neut_WF_array(num_Q,0:neut_vib_num-1), S_ft(num_tot,num_tot),&
    &sorted_S_ft(num_tot,num_tot), D_ft(3, num_tot,0:neut_vib_num-1), sorted_D_ft(3, num_tot,0:neut_vib_num-1), ft_channels(num_tot,4), sorted_ft_channels(num_tot,4))
    call read_ion_WF(geom, ion_WF_array)
    call read_neut_WF(geom, neut_WF_array)

    if (lin_approx .eqv. .true.) then
      call frame_transform(S_Q_interpol(:,:,:,1), dipole_data, geom, neut_WF_array, ion_WF_array, S_ft, D_ft, ft_channels, energy_matrix)
    else
      call frame_transform(S_Q(:,:,:,1), dipole_data, geom, neut_WF_array, ion_WF_array, S_ft, D_ft, ft_channels, energy_matrix)
    end if


    Open(27,file='FT_Smat.dat')
    do i = 1,num_tot
      do j = 1,num_tot
        write(27,*)
        write(27,*) ft_channels(i,:), ft_channels(j,:)
        Write(27,'(100e12.4)') real(S_ft(i,j)), aimag(S_ft(i,j))
      end do
    end do
    Close(27)

    Open(27,file='FT_dips_v0.dat')
    do i = 1,num_tot
        write(27,*)
        write(27,*) ft_channels(i,:)
        Write(27,*) real(D_ft(1,i,0)), aimag(D_ft(1,i,0)), real(D_ft(2,i,0)), aimag(D_ft(2,i,0)), real(D_ft(3,i,0)), aimag(D_ft(3,i,0))
    end do
    Close(27)

    Open(27,file='lmax4emax3.channel')
    do i = 1,num_tot
        write(27,*) i, ft_channels(i,:)
    end do
    Close(27)

! Need to sort all channels by energy in the ft smatrix and dipole moments

    print *, "sorting by energy"
    call sort_index(energy_matrix, index)

    do neut_iv = 0, neut_vib_num-1
      do i = 1, num_tot
        sorted_D_ft(:,i,neut_iv) = D_ft(:,index(i),neut_iv)
        sorted_ft_channels(i,:) = ft_channels(index(i),:)
      end do
    end do

    do i = 1, num_tot
      do j = 1, num_tot
        sorted_S_ft(i,j) = S_ft(index(i),index(j))
      end do
    end do

    Open(27,file='sorted_FT_dips_v0.dat')
    do i = 1,num_tot
        write(27,*)
        write(27,*) sorted_ft_channels(i,:)
        Write(27,*) real(sorted_D_ft(1,i,0)), aimag(sorted_D_ft(1,i,0)), real(sorted_D_ft(2,i,0)), aimag(sorted_D_ft(2,i,0)), real(sorted_D_ft(3,i,0)), aimag(sorted_D_ft(3,i,0))
    end do
    Close(27)

   
    open(27, file = "index.dat")
    do i = 1, num_tot
      write(27,*) index(i)
    end do
    close(27)

    print *, "performing channel elimination"
    print *, ""
    print *, "		GROUND STATE ENERGIES"
    print *, ""
    print *, "	   ", "n	       ", "v	         ", "Threshold"
    print *, "_____________________________________________________________________"
    print *, ""
    do i = 1, num_tot
    if (ft_channels(index(i), 1) == 1 .AND. ft_channels(index(i),3) == 0 .AND. ft_channels(index(i),4) == 0) then
      print *, ft_channels(index(i),1), ft_channels(index(i),2), "	", energy_matrix(i)/ev_au, "(eV)"
    end if
    end do
    print*, ""
    print *, "		EXCITED STATE ENERGIES"
    print *, ""
    print *, "	   ", "n	       ", "v	         ", "Threshold"
    print *, "_____________________________________________________________________"
    print *, ""
    do i = 1, num_tot
    if (ft_channels(index(i), 1) /= 1 .AND. ft_channels(index(i),3) == 0 .AND. ft_channels(index(i),4) == 0) then
      print *, ft_channels(index(i),1), ft_channels(index(i),2), "	", energy_matrix(i)/ev_au, "(eV)"
    end if
    end do

    allocate(TCS_array(energy_step_num, 2,0:neut_vib_num-1), SPES_array(energy_step_num, 2,0:neut_vib_num-1), conv_TCS_array(num_conv_points,2,0:neut_vib_num-1), &
              & max_chans(energy_step_num, 0:neut_vib_num-1), nu(num_tot,energy_step_num, 0:neut_vib_num))
    SPES_array = 0

  do i = 1, num_tot
    if (ft_channels(index(i),1) == 1) then
      sorted_D_ft(:,index(i),:) = dip_mult * sorted_D_ft(:,index(i),:)
    end if
  end do
  nu = 0d0
  do neut_iv = 0, neut_vib_num-1

    energy_begin = chan_ens(1,0) - delta_Es(neut_iv)
    energy_end = chan_ens(n_max, vib_num-1) - delta_Es(neut_iv)
    energy_step = (energy_end - energy_begin)/energy_step_num

    
    do e_i = 1, energy_step_num

      E = energy_begin + (e_i-1)*energy_step

      num_open = 0
            Do j=1,num_tot
              If(energy_matrix(j)-delta_Es(neut_iv)>E)Exit
            Enddo
            num_open=j-1

      if (num_open /= num_tot) then

        num_closed = num_tot - num_open
        allocate(dipole_open(3,num_open),dipole_closed(3,num_closed),IPIV(num_tot))
        allocate(smat_cc(num_closed,num_closed),smat_co(num_closed,num_open),beta(num_tot), Z_norms(num_closed))
        dipole_closed = sorted_D_ft(:,num_open+1:num_tot,neut_iv)
        dipole_open = sorted_D_ft(:,1:num_open,neut_iv)
        smat_cc = Conjg(Transpose(sorted_S_ft(num_open+1:num_tot,num_open+1:num_tot)))
        smat_co = Conjg(Transpose(sorted_S_ft(1:num_open,num_open+1:num_tot)))

        beta = 0d0
        do j = 1, num_closed
          beta(j+num_open) = pi/sqrt(2*((energy_matrix(j+num_open)-delta_Es(neut_iv))-E))
        end do

        
        do j = 1, num_closed
          nu(j+num_open,e_i,neut_iv) = 1/sqrt(2*(energy_matrix(j+num_open)-delta_Es(neut_iv)-E))
        end do

        allocate(dipole_phys(3,num_open))
        dipole_phys = 0d0

        do j = 1, num_closed
          smat_cc(j,j) = smat_cc(j,j) - exp(2d0*ci*beta(j+num_open))
        end do

        call ZGESV(num_closed,num_open,smat_cc,num_closed,IPIV(1:num_closed),smat_co,num_closed,INFO)


        ! do j = 1,3
        !   res_mat(j,:) = MatMul(dipole_closed(j,:),smat_co)
        !   do i = 1, num_open
        !     !Want to pull resonances and possible baseline contributions from only lower v
        !     if (ft_channels(index(i),2) >0) then
        !       dipole_phys(j,i) = dipole_open(j,i)
        !     else
        !       print *, ft_channels(index(i),2)
        !       dipole_phys(j,i) = dipole_open(j,i) - res_mat(j,i)
        !     end if
        !   end do
        ! end do

        do i = 1, num_closed
           Z_norms(i) = sum(abs(dipole_closed(:,i))**2) * (sum(abs(smat_co(i,:))**2) + sum(abs(smat_co(i,:))**2) + sum(abs(smat_co(i,:))**2))
        end do

        max_chans(e_i, neut_iv) = maxloc(Z_norms, 1) + num_open

        do j = 1,3
          dipole_phys(j,:) = dipole_open(j,:) - MatMul(dipole_closed(j,:),smat_co)
        end do

        deallocate(dipole_open, dipole_closed, smat_cc, smat_co, IPIV, beta, Z_norms)

      else if (num_open == num_tot) then

      num_open = num_tot
      allocate(dipole_phys(num_open, 3))
      dipole_phys = sorted_D_ft(:,:,neut_iv)

      end if

      !omega = omega_0 - channel_energies(0)

      TCS_array(e_i,1,neut_iv) = const*E*(sum(abs(dipole_phys(1,:))**2)+sum(abs(dipole_phys(2,:))**2)+sum(abs(dipole_phys(3,:))**2))
      TCS_array(e_i,2,neut_iv) = E
      SPES_array(e_i,2,neut_iv) = TCS_array(e_i,2,neut_iv)

      ! do n = 1, n_max
      !       do vib_i = 0, vib_num-1
      !           if ((abs(E - chan_ens(n,vib_i)) < elec_en_range) .AND. (E - chan_ens(n,vib_i) >= 0)) then
      !               do j = 1, num_open
      !                   if (ft_channels(index(j),1) == n .AND. ft_channels(index(j),2) == vib_i) then
      !                       SPES_array(e_i,1,neut_iv) = SPES_array(e_i,1,neut_iv) + const*omega*(abs(dipole_phys(1,j))**2+abs(dipole_phys(2,j))**2+abs(dipole_phys(3,j))**2)
      !                   end if
      !               end do
      !           end if
      !       end do
      !   end do

      deallocate(dipole_phys)

    end do
    call gauss_conv(TCS_array(:,:,neut_iv), conv_TCS_array(:,:,neut_iv))
  end do
    close(1)


! Writing all the results into output files
  
  
if (lin_approx .eqv. .true.) then
  open(27, file = "TCS_wavelength_lin.dat")
  do e_i = 1, energy_step_num
    write(27,*) 12398.0/(TCS_array(e_i,2,0)/ev_au), TCS_array(e_i,1,0)*(.529**2)*100, 12398.0/(TCS_array(e_i,2,1)/ev_au), TCS_array(e_i,1,1)*(.529**2)*100, 12398.0/(TCS_array(e_i,2,2)/ev_au), TCS_array(e_i,1,2)*(.529**2)*100
  end do
  close(27)

  open(27, file = "TCS_photonenergy_lin.dat")
  do e_i = 1, energy_step_num
    write(27,*) (TCS_array(e_i,2,0)/ev_au), TCS_array(e_i,1,0)*(.529**2)*100, (TCS_array(e_i,2,1)/ev_au), TCS_array(e_i,1,1)*(.529**2)*100, (TCS_array(e_i,2,2)/ev_au), TCS_array(e_i,1,2)*(.529**2)*100
  end do
  close(27)

  open(27, file = "conv_TCS_wavelength_lin.dat")
  do e_i = 1, num_conv_points
    write(27,*) 12398.0/(conv_TCS_array(e_i,2,0)/ev_au), conv_TCS_array(e_i,1,0)*(.529**2)*100, 12398.0/(conv_TCS_array(e_i,2,1)/ev_au), conv_TCS_array(e_i,1,1)*(.529**2)*100, 12398.0/(conv_TCS_array(e_i,2,2)/ev_au), conv_TCS_array(e_i,1,2)*(.529**2)*100
  end do
  close(27)

  open(27, file = "conv_TCS_photonenergy_lin.dat")
  do e_i = 1, num_conv_points
    write(27,*) (conv_TCS_array(e_i,2,0)/ev_au), conv_TCS_array(e_i,1,0)*(.529**2)*100, (conv_TCS_array(e_i,2,1)/ev_au), conv_TCS_array(e_i,1,1)*(.529**2)*100, (conv_TCS_array(e_i,2,2)/ev_au), conv_TCS_array(e_i,1,2)*(.529**2)*100
  end do
  close(27)

else
  open(27, file = "TCS_wavelength.dat")
  do e_i = 1, energy_step_num
    write(27,*) 12398.0/(TCS_array(e_i,2,0)/ev_au), TCS_array(e_i,1,0)*(.529**2)*100, 12398.0/(TCS_array(e_i,2,1)/ev_au), TCS_array(e_i,1,1)*(.529**2)*100, 12398.0/(TCS_array(e_i,2,2)/ev_au), TCS_array(e_i,1,2)*(.529**2)*100
  end do
  close(27)

  open(27, file = "TCS_photonenergy.dat")
  do e_i = 1, energy_step_num
    write(27,*) (TCS_array(e_i,2,0)/ev_au), TCS_array(e_i,1,0)*(.529**2)*100, (TCS_array(e_i,2,1)/ev_au), TCS_array(e_i,1,1)*(.529**2)*100, (TCS_array(e_i,2,2)/ev_au), TCS_array(e_i,1,2)*(.529**2)*100
  end do
  close(27)

  open(27, file = "conv_TCS_wavelength.dat")
  do e_i = 1, num_conv_points
    write(27,*) 12398.0/(conv_TCS_array(e_i,2,0)/ev_au), conv_TCS_array(e_i,1,0)*(.529**2)*100, 12398.0/(conv_TCS_array(e_i,2,1)/ev_au), conv_TCS_array(e_i,1,1)*(.529**2)*100, 12398.0/(conv_TCS_array(e_i,2,2)/ev_au), conv_TCS_array(e_i,1,2)*(.529**2)*100
  end do
  close(27)

  open(27, file = "conv_TCS_photonenergy.dat")
  do e_i = 1, num_conv_points
    write(27,*) (conv_TCS_array(e_i,2,0)/ev_au), conv_TCS_array(e_i,1,0)*(.529**2)*100, (conv_TCS_array(e_i,2,1)/ev_au), conv_TCS_array(e_i,1,1)*(.529**2)*100, (conv_TCS_array(e_i,2,2)/ev_au), conv_TCS_array(e_i,1,2)*(.529**2)*100
  end do
  close(27)
end if


  Open(27,file='Re_Smatr_continuity_Q.dat')
  ie=num_scatE/2
   Write(27,'(A10,f15.12)')'#E in eV',E_el_Rmat(ie)*au_eV
   do i2=1,dim_lm;do j2=i2,dim_lm
      call i2lm(i2 ,n,l,m)
      call i2lm(j2 ,n_p,l_p,m_p)
        Write(27,'(A12,100i4)')'# n,l,m,n_p,l_p,m_p',n,l,m,n_p,l_p,m_p
        Do i_Q=1,num_Q
           Write(27,'(100e12.4)')geom(i_Q),real(S_Q(i2,j2,i_Q,ie)) 
        Enddo
        Write(27,'(A12,100i4)')
   enddo;enddo
  Close(27)

  open(27, file = "Resonances_v0.dat")
  do e_i = 1, energy_step_num
    E = chan_ens(1,0) + (e_i-1)*energy_step
    write(27,*) E/ev_au, ft_channels(index(max_chans(e_i,1)),1), ft_channels(index(max_chans(e_i,1)),2), ft_channels(index(max_chans(e_i,1)),3), ft_channels(index(max_chans(e_i,1)),4), nu(max_chans(e_i,1),e_i,1)
    write(27,*)
  end do
  close(27)
  
  Open(27,file='Re_Smatr_check.dat')
  ie=num_scatE/2
   Write(27,'(A10,f15.12)')'#E in eV',E_el_Rmat(ie)*au_eV
   do i2=1,dim_lm;do j2=i2,dim_lm
      call i2lm(i2 ,n,l,m)
      call i2lm(j2 ,n_p,l_p,m_p)
        Write(27,'(A12,100i4)')'# n,l,m,n_p,l_p,m_p',n,l,m,n_p,l_p,m_p
        Do i_Q=1,num_Q
           Write(27,'(100e12.4)')geom(i_Q),real(S_Q_check(i2,j2,i_Q,ie)) 
        Enddo
        Write(27,'(A12,100i4)')
   enddo;enddo
  Close(27)
  
  Open(27,file='Re_Smatr_interpol.dat')
  ie=num_scatE/2
   Write(27,'(A10,f15.12)')'#E in eV',E_el_Rmat(ie)*au_eV
   do i2=1,dim_lm;do j2=i2,dim_lm
      call i2lm(i2 ,n,l,m)
      call i2lm(j2 ,n_p,l_p,m_p)
        Write(27,'(A12,100i4)')'# n,l,m,n_p,l_p,m_p',n,l,m,n_p,l_p,m_p
        Do i_Q=1,num_Q
           Write(27,'(100e12.4)')geom(i_Q),real(S_Q_interpol(i2,j2,i_Q,ie)) 
        Enddo
        Write(27,'(A12,100i4)')
   enddo;enddo
  Close(27)
  
  Open(27,file='Re_Smatr_extrapol.dat')
   Write(27,'(A10,f15.12)')'#E in eV',E_el_Rmat(1)*au_eV
   do i2=1,dim_lm;do j2=i2,dim_lm
      call i2lm(i2 ,n,l,m)
      call i2lm(j2 ,n_p,l_p,m_p)
        Write(27,'(A12,100i4)')'# n,l,m,n_p,l_p,m_p',n,l,m,n_p,l_p,m_p
        Do i_Q_=1,num_Q_
           Write(27,'(100e12.4)')geom_(i_Q_),real(S_Q_extrapol(i2,j2,i_Q_,1)), aimag(S_Q_extrapol(i2,j2,i_Q_,1))
        Enddo
        Write(27,'(A12,100i4)')
   enddo;enddo
  Close(27)

  Open(27,file='Smatr_extrapol.dat')
  ie=num_scatE/2
   Write(27,'(A9,f15.12)')'# E in eV',E_el_Rmat(ie)*au_eV
   Write(27,'(A19)')'# n,l,m,n_p,l_p,m_p'
   Write(27,'(A20,100i4)')'Number of channels: ', dim_lm ** 2 !n_max = 1
   Write(27,'(A12,100i4)')
   do i2=1,dim_lm;do j2=i2,dim_lm
      call i2lm(i2 ,n,l,m)
      call i2lm(j2 ,n_p,l_p,m_p)
        Write(27,'(100i4)')n,l,m,n_p,l_p,m_p
        Do i_Q_=1,num_Q_
           Write(27,*)geom_(i_Q_),S_Q_extrapol(i2,j2,i_Q_,ie)
        Enddo
        Write(27,'(A12,100i4)')
   enddo;enddo
  Close(27)
  
  Open(27,file='Delta_Q.dat')
  ie=num_scatE/2
   Write(27,'(A10,f15.12)')'#E in eV',E_el_Rmat(ie)*au_eV
   do i2=1,dim_lm;do j2=i2,dim_lm
      call i2lm(i2 ,n,l,m)
      call i2lm(j2 ,n_p,l_p,m_p)
        Write(27,'(A12,100i4)')'# n,l,m,n_p,l_p,m_p',n,l,m,n_p,l_p,m_p
        Do i_Q=1,num_Q
           Write(27,'(100e12.4)')geom(i_Q),Delta_Q(i2,j2,i_Q,ie)
        Enddo
        Write(27,'(A12,100i4)')
   enddo;enddo
  Close(27)
  
  Open(27,file='Delta_Q_interpol.dat')
  ie=num_scatE/2
   Write(27,'(A10,f15.12)')'#E in eV',E_el_Rmat(ie)*au_eV
   do i2=1,dim_lm;do j2=i2,dim_lm
      call i2lm(i2 ,n,l,m)
      call i2lm(j2 ,n_p,l_p,m_p)
        !Write(27,'(A12,100i4)')'# n,l,m,n_p,l_p,m_p',n,l,m,n_p,l_p,m_p
        Do i_Q=1,num_Q
           Write(27,'(100e12.4)')geom(i_Q),real(Delta_Q_interpol(i2,j2,i_Q,1))
        Enddo
        Write(27,'(A12,100i4)')
   enddo;enddo
  Close(27)
  
  Open(27,file='Delta_Q_extrapol.dat')
  ie=num_scatE/2
   !Write(27,'(A10,f15.12)')'#E in eV',E_el_Rmat(ie)*au_eV
   do i2=1,dim_lm;do j2=i2,dim_lm
      call i2lm(i2 ,n,l,m)
      call i2lm(j2 ,n_p,l_p,m_p)
        Write(27,'(A12,100i4)')'# n,l,m,n_p,l_p,m_p',n,l,m,n_p,l_p,m_p
        Do i_Q_=1,num_Q_
           Write(27,'(100e12.4)')geom_(i_Q_),Delta_Q_extrapol(i2,j2,i_Q_,1)
        Enddo
        Write(27,'(A12,100i4)')
   enddo;enddo
  Close(27)

  open(27, file = 'ion_WF.dat')
  do i = 0, vib_num-1
    do j = 1, n_max
      write(27,*)
      !write(27,*) i, j
      do i_Q_ = 1,num_Q_
        write(27,*) geom_(i_Q_), ion_WF_array(i_Q_,i,j)
      end do
    end do
  end do
  close(27)

  open(2, file = "./dips_imag_exc.dat")
  do i = 1, dim_lm
    call i2lm(i,n,l,m)
      if (n >= 2) then
          write(2, *) i,n,l,m
          do i_Q_ = 1, num_Q_
              r = geom_(i_Q_)
              if (aimag(dipole_data(1,i,i_Q_)) /= 0) then
                  write(2, *) r, aimag(dipole_data(1,i,i_Q_))
              else if (aimag(dipole_data(2,i,i_Q_)) /= 0) then
                  write(2, *) r, aimag(dipole_data(2,i,i_Q_))
              else if (aimag(dipole_data(3,i,i_Q_)) /= 0) then
                  write(2, *) r, aimag(dipole_data(3,i,i_Q_))
              else
                  write(2,*) r, aimag(dipole_data(3,i,i_Q_))
              end if
          end do
          write(2,*)
      end if
  end do
  close(2)
  open(2, file = "./dips_real_exc.dat")
  do i = 1, dim_lm
    call i2lm(i,n,l,m)
      if (n >= 2) then
          write(2, *) i,n,l,m
          do i_Q_ = 1, num_Q_
              r = geom_(i_Q_)
              if (real(dipole_data(1,i,i_Q_)) /= 0) then
                  write(2, *) r, real(dipole_data(1,i,i_Q_))
              else if (real(dipole_data(2,i,i_Q_)) /= 0) then
                  write(2, *) r, real(dipole_data(2,i,i_Q_))
              else if (real(dipole_data(3,i,i_Q_)) /= 0) then
                  write(2, *) r, real(dipole_data(3,i,i_Q_))
              else
                  write(2,*) r, real(dipole_data(3,i,i_Q_))
              end if
          end do
          write(2,*)
      end if
  end do
  close(2)
  open(2, file = "./dips_imag_ground.dat")
  do i = 1, dim_lm
    call i2lm(i,n,l,m)
      if (n == 1) then
          write(2, *) i,n,l,m
          do i_Q_ = 1, num_Q_
              r = geom_(i_Q_)
              if (aimag(dipole_data(1,i,i_Q_)) /= 0) then
                  write(2, *) r, aimag(dipole_data(1,i,i_Q_))
              else if (aimag(dipole_data(2,i,i_Q_)) /= 0) then
                  write(2, *) r, aimag(dipole_data(2,i,i_Q_))
              else if (aimag(dipole_data(3,i,i_Q_)) /= 0) then
                  write(2, *) r, aimag(dipole_data(3,i,i_Q_))
              else
                  write(2,*) r, aimag(dipole_data(3,i,i_Q_))
              end if
          end do
          write(2,*)
      end if
  end do
  close(2)
  open(2, file = "./dips_real_ground.dat")
  do i = 1, dim_lm
    call i2lm(i,n,l,m)
      if (n == 1) then
          write(2, *) i,n,l,m
          do i_Q_ = 1, num_Q_
              r = geom_(i_Q_)
              if (real(dipole_data(1,i,i_Q_)) /= 0) then
                  write(2, *) r, real(dipole_data(1,i,i_Q_))
              else if (real(dipole_data(2,i,i_Q_)) /= 0) then
                  write(2, *) r, real(dipole_data(2,i,i_Q_))
              else if (real(dipole_data(3,i,i_Q_)) /= 0) then
                  write(2, *) r, real(dipole_data(3,i,i_Q_))
              else
                  write(2,*) r, real(dipole_data(3,i,i_Q_))
              end if
          end do
          write(2,*)
      end if
  end do
  close(2)


  ! open(2, file = "superunitpart.dat")
  ! call lm2i(4,-3,i)
  ! call i2lm(i, l, m)
  ! ie = num_scatE/2
  !     write(2, *) "Partial Wave:", l, m, l, m
  !     write(2, *) "Geom  ", "Smat  ", "Delta  ", "eigenphase_matrix  "
  !     do i_Q_ = 1, num_Q_
  !     write(2,*) geom_(i_Q_), S_Q_extrapol(i,i,i_Q_,ie), Delta_Q_extrapol(i,i,i_Q_,ie), tmp_eigenphases_matrix_(i,i)
  !     end do
  !     write(2, *) 
  !     write(2, *) "Kmat"
  !     do i_Q = 1, num_Q
  !       write(2,*) geom_(i_Q), K_Q(i,i,i_Q,ie)
  !     end do
  ! close(2)

  
deallocate(E_el_Rmat,S_Q,geom,D2h_Dinf,md,S_Q_check,K_Q,&
a,U,tmp_eigenphases,eigenphases,Delta_Q,M_print,&
tmp_Delta,tmp_Delta_,eigenphases_matrix,Delta_Q_interpol,&
S_Q_interpol,Delta_Q_extrapol,S_Q_extrapol,&
tmp_eigenphases_matrix,tmp_eigenphases_matrix_)

End subroutine Get_Extrapolated_SQ_matrix
    
function cex2(nj)
  character*2 cex2
  integer cnj,dnj,unj,nj
  cnj = nj/100
  dnj = (nj - cnj*100)/10
  unj = nj - cnj*100 - dnj*10
!  cex2 = char(cnj+48)//char(dnj+48)//char(unj+48)
  cex2 = char(dnj+48)//char(unj+48)
end function cex2
!----------------------------------------------------
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
!----------------------------------------------------
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
! ****************************************************************
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
 End Subroutine K_S
 ! ***************************************************************
 subroutine find_rderiv(y_fin, y_in, dx, deriv)

	real(8), intent(in)	::	y_fin, y_in
	real(8), intent(in)	::	dx
	real(8), intent(out)	::	deriv

deriv = (y_fin - y_in)/dx

end subroutine
 ! ***************************************************************
subroutine read_ion_WF(geom, ion_WF_array)

  use global_params, only : vib_num, n_max, num_Q

	character (len = 286)				::	file_path

	integer						::	state, i_Q_, chan_i
	character (len = 1)			::	chan_str, state_str
	real(8)						::	r
  real(8), dimension(num_Q), intent(in)  ::  geom

	real(8), dimension(num_Q,0:vib_num-1,n_max), intent(out)	::	ion_WF_array


do state = 1, n_max
	do chan_i = 0, vib_num-1
    print *, "n:", state, "v:", chan_i
	        do i_Q_ = 1, num_Q
	                r = geom(i_Q_)

	                write(chan_str, "(i1)") chan_i
                    write(state_str, "(i1)") state

	                file_path = "./Ion_WF/Target_wf" // chan_str // "_N2+_ElecState00" // state_str // ".dat"
	                open(3, file = file_path, status = "old")
	                call WF_mag(3, r, ion_WF_array(i_Q_,chan_i,state))

	        end do
	end do
end do


end subroutine
 ! ***************************************************************
subroutine WF_mag(file_num, r, WF)

	real(8)                  :: r
	integer                  :: file_num, a, eof
	real(8), intent(out)     :: WF
	real(8), allocatable     :: WF_data(:), r_mat(:), sub_mat(:)
	!real(8)                  :: tmp(100)


    allocate(WF_data(1000), r_mat(1000), sub_mat(1000))
    r_mat = 100
    sub_mat = 100
    WF_data = 100
    do a = 1, 1000
        read(file_num, *, iostat=eof) r_mat(a), WF_data(a)
        if (eof<0) then
            exit
        end if
    end do
    sub_mat = abs(r_mat - r)
    !print *, "r:", r
    !print *, "loc:", minloc(sub_mat, 1)
    !print *, sub_mat(minloc(sub_mat, 1))
    
    WF = WF_data(minloc(sub_mat, 1))
    !print *, "WF:", WF

	close(file_num)
	deallocate(WF_data, r_mat, sub_mat)
end subroutine
! ***************************************************************
subroutine frame_transform(S_Q_, dipole_data, geom, neut_WF_array, ion_WF_array, S_ft, D_ft, ft_channels, energy_matrix)

use global_params, only : vib_num, n_max, num_Q, dim_lm, chan_ens, num_tot, neut_vib_num


complex(8), dimension(dim_lm,dim_lm,num_Q), intent(in) ::  S_Q_
complex(8), dimension(3,dim_lm,num_Q) :: dipole_data
real(8), dimension(num_Q), intent(in)  ::  geom
real(8), dimension(num_Q,0:vib_num-1,n_max), intent(in)	::	ion_WF_array
real(8), dimension(num_Q, 0:neut_vib_num-1), intent(in)	::	neut_WF_array

complex(8), dimension(num_tot,num_tot), intent(out) ::  S_ft
complex(8), dimension(3,num_tot,0:neut_vib_num-1), intent(out) ::  D_ft
integer, dimension(num_tot,4), intent(out) ::  ft_channels
real(8), dimension(num_tot) :: energy_matrix

complex(8) :: S_sum, D_sum
integer :: i, j, iv, i_p, j_p, iv_p, i_Q_,n,l,m,n_p,l_p,m_p,i_coord,vib_i, neut_iv
real(8) :: step

step = geom(2) - geom(1)

do i = 1, dim_lm
  do iv = 0, vib_num-1
    j = (i-1)*vib_num + iv+1
    call i2lm(i,n,l,m)
    energy_matrix(j) = chan_ens(n,iv)
    ft_channels(j,:) = (/n,iv,l,m/)
    do i_p = 1, dim_lm
        do iv_p = 0, vib_num-1
          call i2lm(i_p,n_p,l_p,m_p)
          j_p = (i_p-1)*vib_num + iv_p+1
          S_sum = 0
          do i_Q_ = 1, num_Q
            S_sum = S_sum + S_Q_(i,i_p,i_Q_) * ion_WF_array(i_Q_,iv,n) * ion_WF_array(i_Q_, iv_p, n_p) * step
          end do
          S_ft(j,j_p) = S_sum
        end do
    end do
    do neut_iv = 0,neut_vib_num-1
      do i_coord = 1,3
        D_sum = 0
        do i_Q_ = 1, num_Q
          D_sum = D_sum + dipole_data(i_coord,i,i_Q_) * ion_WF_array(i_Q_, iv, n) * neut_WF_array(i_Q_,neut_iv) * step
        end do
        D_ft(i_coord,j,neut_iv) = D_sum
      end do
    end do
  end do
end do

end subroutine

subroutine read_neut_WF(geom_, neut_WF_array)

  use global_params, only : num_Q, neut_vib_num

	character (len = 65)				::	file_path

  real(8), dimension(num_Q), intent(in)  ::  geom_

	integer						::	i_Q_, vib_i
	real(8)						::	r
    character (len = 1)         ::  vib_str

	real(8), intent(out), dimension(num_Q,0:neut_vib_num-1)	::	neut_WF_array

do i_Q_ = 1, num_Q
    r = geom_(i_Q_)
    do vib_i = 0, neut_vib_num-1
        write(vib_str, "(i1)") vib_i
        file_path = "./Neutral_WF/wf00" // vib_str // "_neut.dat"
        open(3, file = trim(file_path), status = "old")
        call WF_mag(3, r, neut_WF_array(i_Q_,vib_i)) !Replaced vib_i w/ 0 for one state
    end do
end do

open(27, file = "neut_WF_vals.dat")
do vib_i = 0, neut_vib_num-1 
  do i_Q_ = 1, num_Q
  write(27,*) neut_WF_array(i_Q_,vib_i)
  end do
  write(27,*)
end do
close(27)

end subroutine

subroutine gauss_conv(TCS_array, conv_TCS_array)

  use global_params, only: num_conv_points, energy_step_num, sigma, pi

  real(8), dimension(energy_step_num,2), intent(in)	::	TCS_array
  
  real(8), dimension(num_conv_points, 2), intent(out)	::	conv_TCS_array
  
  real(8)			::	E_in, E_fin, conv_E_step, E_curr, TCS_E_step, int_sum
  integer			::	i, j
  
  real(8), allocatable			::	my_list(:,:), cauchy_list(:)

  
  allocate(my_list(energy_step_num,2), cauchy_list(energy_step_num))

  E_in = TCS_array(1,2)
  E_fin = TCS_array(energy_step_num,2)
  TCS_E_step = TCS_array(2,2) - TCS_array(1,2)
  conv_E_step = (E_fin-E_in)/num_conv_points
  
  print *, "Using ", num_conv_points, " points for convolution grid"
  
  do i = 1, num_conv_points
  int_sum = 0
  E_curr = E_in + (conv_E_step*(i-1))
      if (E_curr >= TCS_array(1,2) .AND. E_curr <= TCS_array(energy_step_num,2)) then
          do j = 1, energy_step_num
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