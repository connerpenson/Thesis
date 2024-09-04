Module global_params

    integer, parameter :: num_gamma_Q = 8, l_max = 4, n_max = 4, num_Q = 19, num_Q_ = 500, skip_en = 4, num_scatE = 85, energy_step_num = 500
    integer, parameter :: dim_lm    = (l_max+1)**2 * n_max
    integer, dimension(n_max), parameter    ::  	vib_nums = (/5,5,5,5/) !specify num v for each elec state
    integer, parameter                      ::      vib_max = maxval(vib_nums,1)
    integer, parameter :: num_tot   = (l_max+1)**2 * sum(vib_nums)
    real(8), parameter :: pi = 3.14159265359
    real(8), dimension(num_Q)     ::   geoms
    integer, dimension(num_Q)     ::   ordered_geoms
    integer, dimension(num_Q)     ::   all_open_i
    real(8), dimension(n_max,0:vib_max-1) ::  chan_ens
    character(len=1), dimension(1:n_max) :: state_names
    real(8), dimension(energy_step_num) ::  total_energies
    complex*16, parameter   ::  ci = (0.0,1.0)
    integer, parameter      ::  lin_start = 8, lin_end = 13
    real(8), parameter      ::  lin_initial_r = 1.7, lin_final_r = 3.3, trim_value = 5, deriv_trim = 1
    real(8)                 ::  lin_r_step = (lin_final_r - lin_initial_r)/num_Q_
    real(8), dimension(num_Q_)  ::  lin_geoms
    integer, dimension(5,6)    ::  good_elem_list
    character(len=15), dimension(num_gamma_Q) ::  name_of_file(1:num_gamma_Q) = (/ '1Ag.channel ', '1Au.channel ', '1B1g.channel', '1B1u.channel', &
                                                                                '1B2g.channel', '1B2u.channel', '1B3g.channel', '1B3u.channel'/)!, &
                                                                                !'3Ag.channel ', '3Au.channel ', '3B1g.channel', '3B1u.channel', &
                                                                                !'3B2g.channel', '3B2u.channel', '3B3g.channel', '3B3u.channel' /)

    character(len=9)   ::  wf_dir = "./Ion WF/"
    
    contains

    subroutine get_lin_geoms

        implicit none

        integer ::  i_Q_

        do i_Q_ = 1, num_Q_
            lin_geoms(i_Q_) = lin_initial_r + lin_r_step*(i_Q_ - 1)
        end do

    end subroutine

    subroutine get_geoms

        implicit none

        integer ::  i_Q, i, j, k
        real(8)     ::  tmp
        character(len=2)    ::  cex2
        integer ::  cnj, dnj, unj

        Do i_Q=1,num_Q
            cnj = i_Q/100
            dnj = (i_Q - cnj*100)/10
            unj = i_Q - cnj*100 - dnj*10
            !  cex2 = char(cnj+48)//char(dnj+48)//char(unj+48)
            cex2 = char(dnj+48)//char(unj+48)
            open(unit = 7, file='K-mats/'//cex2//'/Q.txt', action='read')
            read(7,*) geoms(i_Q)
            !print *,'geometries are',geom(i_Q)
        Enddo
        
          ! ordering geometries
        Do i=1,num_Q
            ordered_geoms(i) = i
        Enddo
        Do i=1,num_Q
            Do j=i+1,num_Q
                If (geoms(i) > geoms(j)) Then
                    tmp = geoms(i)
                    geoms(i) = geoms(j)
                    geoms(j) = tmp
                    k = ordered_geoms(i)
                    ordered_geoms(i) = ordered_geoms(j)
                    ordered_geoms(j) = k
                Endif
            Enddo
        Enddo

    end subroutine

    subroutine get_channel_ens

        implicit none

        integer :: ie, iv, a
        real(8) ::  tmp

        state_names = [character(len=1) :: "X", "A", "A", "B"]

        do ie = 1, n_max
            ! cnj = ie/100
            ! dnj = (ie - cnj*100)/10
            ! unj = ie - cnj*100 - dnj*10
            ! cex2 = char(cnj+48)//char(dnj+48)//char(unj+48)
            ! cex2 = char(dnj+48)//char(unj+48)
            open(333, file = "./Energies/" // state_names(ie) // "/Target_N2_N2__ElecState001_E.dat", status = "old")
            read(333,*)
            read(333,*)
            read(333,*)
            do iv = 0, vib_nums(ie)-1
                read(333,*) a, tmp, chan_ens(ie, iv)
            end do
            close(333)
        end do
        chan_ens = chan_ens - chan_ens(1,0)

        do ie = 1, n_max
            do iv = 0, vib_nums(ie)-1
                open(333, file = "chan_ens_au.dat")
                write(333, *) chan_ens(ie,iv)
            end do
        end do
        close(333)

    end subroutine

End Module global_params