program Photo_Trans_Calc

    use global_params, only : e_max, v_max, dips_dir

    implicit NONE
    
    character (len = 18)				::	bsp_file_path
    
    real(8), allocatable  ::  geoms(:)
    
    integer						::	i_q, vib_i, elec_i, io, num_q, v_p, v_pp, e_pp, e_p
    real(8)						::	r
    character (len = 1)         ::  vib_str, elec_str

    real(8), dimension(0:v_max, 0:v_max, 0:0, e_max)     ::  einst_coeffs, franck_condon
    real(8), dimension(0:e_max, 0:v_max)    ::  vib_ens    
    
    real(8), allocatable	::	WF_array(:,:,:) ! (geom,vib,elec)
    real(8), allocatable    ::  dip_data(:,:)     ! (geom,elec)


    !make geometry array
    !how many points?
    num_q = 0
    open(333, file = dips_dir // "N2+_2_transdip.dat", status = "old")
    do
        read(333, *, iostat = io)
        
        if (io /= 0) exit
        num_q = num_q + 1
        
    end do
    close(333)

    print *, num_q

    !make allocations with num_q
    allocate(geoms(num_q), WF_array(num_q, 0:v_max, 0:e_max), dip_data(num_q, e_max))

    !get geoms
    open(333, file = dips_dir // "N2+_2_transdip.dat", status = "old")
    do i_q = 1, num_q
        read(333, *) geoms(i_q)
    end do
    close(333)

    do elec_i = 0, e_max
        call make_bsp_name(elec_i, bsp_file_path)
        do i_q = 1, num_q
            call get_vals_bsp_wf(bsp_file_path, geoms(i_q), WF_array(i_q,:,elec_i))
        end do
    end do

    
    open(27, file = "bsp_WF_vals.dat")
    do vib_i = 0, v_max-1 
        do elec_i = 0, e_max
            write(27,*) "#e=",elec_i,"v=", vib_i
            do i_q = 1, num_q
                write(27,*) geoms(i_q), WF_array(i_q,vib_i,elec_i)
            end do
            write(27,*)
        end do
        write(27,*)
    end do
    close(27)

    call read_dips(num_q, geoms, dip_data)

    call get_vib_ens(vib_ens)

    call calc_einst_coeffs(num_q, geoms, WF_array, dip_data, vib_ens, einst_coeffs)

    call calc_franck_condon(num_q, geoms, WF_array, franck_condon)

    open(333, file = "EC_BSP_XtoA.dat")
    do v_p = 0, v_max
        do v_pp = 0, v_max
            write(333, *) "vp = ", v_p, "v_pp = ", v_pp, "D = ", einst_coeffs(v_p,v_pp,0,1)
            write(333, *)
        end do
    end do
    close(333)

    open(333, file = "EC_BSP_XtoB.dat")
    do v_p = 0, v_max
        do v_pp = 0, v_max
            write(333, *) "vp = ", v_p, "v_pp = ", v_pp, "D = ", einst_coeffs(v_p,v_pp,0,2)
            write(333, *)
        end do
    end do
    close(333)

    open(333, file = "FC_BSP_XtoA.dat")
    do v_p = 0, v_max
        do v_pp = 0, v_max
            write(333, *) "vp = ", v_p, "v_pp = ", v_pp, "D = ", franck_condon(v_p,v_pp,0,1)
            write(333, *)
        end do
    end do
    close(333)

    open(333, file = "FC_BSP_XtoB.dat")
    do v_p = 0, v_max
        do v_pp = 0, v_max
            write(333, *) "vp = ", v_p, "v_pp = ", v_pp, "D = ", franck_condon(v_p,v_pp,0,2)
            write(333, *)
        end do
    end do
    close(333)
    
    open(333, file = "vib_ens.dat")
    do v_p = 0, v_max
        do v_pp = 0, v_max
            do e_p = 0, e_max
                do e_pp = 0, e_max
                    write(333, *) "vp = ", v_p, "vpp = ", v_pp, "ep = ", e_p, "epp = ", e_pp, "trans_freq = ", ( vib_ens(e_pp, v_pp) - vib_ens(e_p, v_p) )
                    write(333, *)
                end do 
            end do
        end do
    end do
    close(333)

end program

subroutine get_vib_ens(vib_ens)

    use global_params, only : v_max, e_max, au_to_cm

    real(8), dimension(0:e_max, 0:v_max), intent(out)    ::  vib_ens

    integer     ::      elec_i, vib_i
    character(len=24)   ::  file_path
    real(8)     ::      energy

    do elec_i = 0, e_max
        call make_ens_name(elec_i, file_path)
        open(333, file = file_path, status = "old")
        read(333,*)
        do vib_i = 0, v_max
            read(333,*) energy !reads in au, formula needs cm^-1
            print *, elec_i, vib_i, energy, energy * au_to_cm
            vib_ens(elec_i,vib_i) = energy * au_to_cm
        end do
    end do

end subroutine

subroutine make_ens_name(elec_num, ens_name)

    use global_params, only : ens_dir

    integer, intent(in)     ::  elec_num
    character (len = 1)     ::  elec_str
    character(len = 24), intent(out) ::  ens_name

    write(elec_str, "(i1)") elec_num
    ens_name = ens_dir // elec_str // "_bsp.dat"


end subroutine
      
subroutine make_bsp_name(vib_num, bsp_name)

    use global_params, only : bsp_wf_dir

    integer, intent(in)     ::  vib_num
    character (len = 1)     ::  vib_str
    character(len = 18), intent(out) ::  bsp_name

    write(vib_str, "(i1)") vib_num

    bsp_name = bsp_wf_dir // "wf" // vib_str // ".dat"

end subroutine

subroutine get_vals_bsp_wf(file_path, r, WF)

    use global_params, only : v_max

    real(8)                  :: r
    integer                  :: a, io, vib_i
    character(len = 18), intent(in) ::  file_path
    real(8), dimension(0:v_max), intent(out)     :: WF
    real(8), dimension(0:v_max, 1000)     :: WF_mags
    real(8), allocatable     :: WF_data(:,:), r_mat(:), sub_mat(:)
    !real(8)                  :: tmp(100)

    allocate(WF_data( ((v_max+1)*2), 1000), r_mat(1000), sub_mat(1000))
    r_mat = 100
    sub_mat = 100
    WF_data = 0
    open(333, file = file_path, status = "old")
    do a = 1, 1000
        read(333, *, iostat=io) r_mat(a), WF_data(:,a)
        if (io /= 0) exit
        !print *, WF_data(:,a)
    end do
    close(333)
    
    
    !Getting norm, but imaginary parts are tiny so its basically just the real part
    !do a = 1, 1000
        !do vib_i = 0, v_max
            !print 
            !WF_mags(vib_i, a) = sqrt(WF_data( (vib_i * 2)+1, a)**2 + WF_data( ((vib_i+1) * 2), a)**2)
            !print *, WF_mags(:,a)
        !end do
    !end do

    !Just taking the real part
    do a = 1, 1000
        do vib_i = 0, v_max
            WF_mags(vib_i, a) = WF_data( (vib_i * 2)+1, a)
        end do
        !print *, WF_mags(:, a)
    end do

    sub_mat = abs(r_mat - r)
    !print *, "r:", r
    !print *, "loc:", minloc(sub_mat, 1)
    !print *, sub_mat(minloc(sub_mat, 1))
    
    WF = WF_mags(:, minloc(sub_mat, 1))
    !print *, "WF:", WF

    deallocate(WF_data, r_mat, sub_mat)

end subroutine

subroutine read_dips(num_q, geoms, dip_data)

    use global_params, only : e_max


    integer, intent(in) ::  num_q
    real(8), dimension(num_q), intent(in) ::  geoms
    real(8), dimension(num_q, e_max), intent(out)   ::  dip_data

    integer ::  i_q, elec_i
    character(len = 25) ::  dip_file_path
    real(8) ::  tmp

    do elec_i = 1, e_max
        call make_dip_name(elec_i, dip_file_path)
        open(333, file = dip_file_path, status = "old")
        do i_q = 1, num_q
            read(333, *) tmp, dip_data(i_q, elec_i)
        end do
    end do

    open(27, file = "dips.dat")
    do elec_i = 1, e_max
        write(27,*) "e=",elec_i
        do i_q = 1, num_q
            write(27,*) geoms(i_q), dip_data(i_q,elec_i)
        end do
        write(27,*)
    end do
    close(27)

end subroutine

subroutine make_dip_name(elec_num, dip_name)

    use global_params, only : dips_dir

    integer, intent(in)     ::  elec_num 
    character (len = 1)     ::  elec_str
    character(len = 25), intent(out) ::  dip_name

    write(elec_str, "(i1)") elec_num

    dip_name = dips_dir // "N2+_" // elec_str // "_transdip.dat"

end subroutine

subroutine calc_einst_coeffs(num_q, geoms, WF_array, dip_data, vib_ens, einst_coeffs)

    use global_params, only : e_max, v_max, const

    real(8), dimension(0:v_max, 0:v_max, 0:0, e_max), intent(out)  ::  einst_coeffs

    integer, intent(in)     ::  num_q
    real(8), dimension(num_q), intent(in) ::  geoms
    real(8), dimension(num_q, 0:v_max, 0:e_max), intent(in)    ::  WF_array
    real(8), dimension(0:e_max, 0:v_max), intent(in)           ::  vib_ens
    real(8), dimension(num_q, e_max), intent(in)                  ::  dip_data

    integer     ::  v_p, v_pp, e_p, e_pp

    e_p = 0 !Calculating coeffs from the ground state
 
    do e_pp = 1, e_max  !To first and second exc
        do v_p = 0, v_max
            do v_pp = 0, v_max
                call integrate_EC(num_q, geoms, WF_array(:,v_p,e_p), WF_array(:,v_pp,e_pp), dip_data(:,e_pp), einst_coeffs(v_p,v_pp,e_p,e_pp))
                einst_coeffs(v_p,v_pp,e_p,e_pp) = const * (vib_ens(e_pp,v_pp) - vib_ens(e_p,v_p))**(3) * einst_coeffs(v_p,v_pp,e_p,e_pp)**2
                print *, vib_ens(e_pp,v_pp) - vib_ens(e_p,v_p)
            end do
        end do
    end do
  

    
end subroutine

subroutine calc_franck_condon(num_q, geoms, WF_array, franck_condon)

    use global_params, only : e_max, v_max

    real(8), dimension(0:v_max, 0:v_max, 0:0, e_max), intent(out)  ::  franck_condon

    integer, intent(in)     ::  num_q
    real(8), dimension(num_q), intent(in) ::  geoms
    real(8), dimension(num_q, 0:v_max, 0:e_max), intent(in)    ::  WF_array

    integer     ::  v_p, v_pp, e_p, e_pp

    e_p = 0 !Calculating coeffs from the ground state
    do e_pp = 1, e_max  !To first and second exc
        do v_p = 0, v_max
            do v_pp = 0, v_max
                call integrate_FC(num_q, geoms, WF_array(:,v_p,e_p), WF_array(:,v_pp,e_pp), franck_condon(v_p,v_pp,e_p,e_pp))
                franck_condon(v_p,v_pp,e_p,e_pp) = franck_condon(v_p,v_pp,e_p,e_pp)**2
            end do
        end do
    end do

    
end subroutine

subroutine integrate_EC(num_q, geoms, WF_array_p, WF_array_pp, dip_data, int_val)

    use global_params, only : e_max, v_max

    integer, intent(in) ::  num_q
    real(8), dimension(num_q), intent(in) ::  geoms
    real(8), dimension(num_q), intent(in) ::  WF_array_p, WF_array_pp
    real(8), dimension(num_q), intent(in) ::  dip_data

    real(8), intent(out)    ::  int_val

    real(8)    ::  step
    integer ::  i_q

    int_val = 0
    do i_q = 2, num_q
        step = geoms(i_q) - geoms(i_q-1)
        !print *, geoms(i_q), WF_array_p(i_q), WF_array_pp(i_q), dip_data(i_q), step
        int_val = int_val + ( WF_array_p(i_q) * WF_array_pp(i_q) * dip_data(i_q) * step )
    end do
    !print *, "------------------------------", int_val

end subroutine


subroutine integrate_FC(num_q, geoms, WF_array_p, WF_array_pp, int_val)

    use global_params, only : e_max, v_max

    integer, intent(in) ::  num_q
    real(8), dimension(num_q), intent(in) ::  geoms
    real(8), dimension(num_q), intent(in) ::  WF_array_p, WF_array_pp

    real(8), intent(out)    ::  int_val

    real(8)    ::  step
    integer ::  i_q

    int_val = 0
    do i_q = 2, num_q
        step = geoms(i_q) - geoms(i_q-1)
        int_val = int_val + ( WF_array_p(i_q) * WF_array_pp(i_q) * step )
    end do

end subroutine