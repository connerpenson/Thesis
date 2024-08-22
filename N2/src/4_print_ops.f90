module print_ops

    use constants, only: r_min_lin, r_min, r_interval, interval_num, num_channels, get_quant_nums, lin_interval, lin_size, au_ang
    use math_ops, only: linearize
    implicit none

    contains

    subroutine write_dips(dipole_data)

        complex(8), intent(in)  ::  dipole_data(:,:,:)
        real(8)                 ::  lin_dips_real(3,num_channels,0:lin_size), lin_dips_imag(3,num_channels,0:lin_size)
        real(8)                 ::  lin_dips(3,num_channels,0:lin_size)

        real(8)     ::  r

        integer     ::  quant_nums(num_channels,4)
        integer     ::  r_i, i, j

        call get_quant_nums(quant_nums)
        open(2, file = "./dips_real_exc.dat")
        do i = 1, num_channels
                if (quant_nums(i,1) >= 2 .and. quant_nums(i,2) == 0) then
                    !write(2, *) i, quant_nums(i,:)
                    do r_i = 1, interval_num + 1
                        r = r_min + (r_interval*(r_i-1))
                        if (real(dipole_data(1,i,r_i)) /= 0) then
                            write(2, *) r, real(dipole_data(1,i,r_i))
                        else if (real(dipole_data(2,i,r_i)) /= 0) then
                            write(2, *) r, real(dipole_data(2,i,r_i))
                        else if (real(dipole_data(3,i,r_i)) /= 0) then
                            write(2, *) r, real(dipole_data(3,i,r_i))
                        else
                            write(2,*) r, real(dipole_data(3,i,r_i))
                        end if
                    end do
                    write(2,*)
                end if
        end do
        close(2)

        open(2, file = "./interp_dipoles.dat")
        do i = 1, num_channels
            do j = 1,3
                lin_dips_real(j,i,:) = linearize(real(dipole_data(j,i,:)))
                lin_dips_imag(j,i,:) = linearize(aimag(dipole_data(j,i,:)))
            end do

            if (quant_nums(i,2) == 0) then
                write(2, *) quant_nums(i,1), quant_nums(i,3:4)
                do r_i = 0, lin_size
                    write(2, *) lin_dips_real(1,i,r_i),lin_dips_imag(1,i,r_i),lin_dips_real(2,i,r_i),lin_dips_imag(2,i,r_i),lin_dips_real(3,i,r_i),lin_dips_imag(3,i,r_i)
                end do
                write(2,*)
            end if
        end do
        close(2)

        open(2, file = "./dipoles.dat")
        do i = 1, num_channels
            if (quant_nums(i,2) == 0) then
                write(2, *) quant_nums(i,1), quant_nums(i,3:4)
                do r_i = 1, interval_num+1
                    write(2, *) real(dipole_data(1,i,r_i)),aimag(dipole_data(1,i,r_i)),real(dipole_data(2,i,r_i)),aimag(dipole_data(2,i,r_i)),real(dipole_data(3,i,r_i)),aimag(dipole_data(3,i,r_i))
                end do
                write(2,*)
            end if
        end do
        close(2)

        open(2, file = "./dips_real_ground_lin.dat")
        do i = 1, num_channels
            do j = 1,3
                lin_dips(j,i,:) = linearize(real(dipole_data(j,i,:)))
            end do
            if (quant_nums(i,1) == 1 .and. quant_nums(i,2) == 0) then
                do r_i = 1, lin_size
                    r = r_min_lin + (lin_interval*(r_i-1))
                    if (lin_dips(1,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(1,i,r_i)
                    else if (lin_dips(2,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(2,i,r_i)
                    else if (lin_dips(3,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(3,i,r_i)
                    else
                        write(2,*) r, lin_dips(3,i,r_i)
                    end if
                end do
                write(2,*)
            end if
        end do
        close(2)

        open(2, file = "./dips_real_exc_lin.dat")
        do i = 1, num_channels
            do j = 1,3
                lin_dips(j,i,:) = linearize(real(dipole_data(j,i,:)))
            end do
            if (quant_nums(i,1) == 2 .and. quant_nums(i,2) == 0) then
                do r_i = 1, lin_size
                    r = r_min_lin + (lin_interval*(r_i-1))
                    if (lin_dips(1,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(1,i,r_i)
                    else if (lin_dips(2,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(2,i,r_i)
                    else if (lin_dips(3,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(3,i,r_i)
                    else
                        write(2,*) r, lin_dips(3,i,r_i)
                    end if
                end do
                write(2,*)
            end if
        end do
        close(2)

        open(2, file = "./dips_imag_ground_lin.dat")
        do i = 1, num_channels
            do j = 1,3
                lin_dips(j,i,:) = linearize(aimag(dipole_data(j,i,:)))
            end do
            if (quant_nums(i,1) == 1 .and. quant_nums(i,2) == 0) then
                do r_i = 1, lin_size
                    r = r_min_lin + (lin_interval*(r_i-1))
                    if (lin_dips(1,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(1,i,r_i)
                    else if (lin_dips(2,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(2,i,r_i)
                    else if (lin_dips(3,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(3,i,r_i)
                    else
                        write(2,*) r, lin_dips(3,i,r_i)
                    end if
                end do
                write(2,*)
            end if
        end do
        close(2)

        open(2, file = "./dips_imag_exc_lin.dat")
        do i = 1, num_channels
            do j = 1,3
                lin_dips(j,i,:) = linearize(aimag(dipole_data(j,i,:)))
            end do
            if (quant_nums(i,1) == 2 .and. quant_nums(i,2) == 0) then
                do r_i = 1, lin_size
                    r = r_min_lin + (lin_interval*(r_i-1))
                    if (lin_dips(1,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(1,i,r_i)
                    else if (lin_dips(2,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(2,i,r_i)
                    else if (lin_dips(3,i,r_i) /= 0) then
                        write(2, *) r, lin_dips(3,i,r_i)
                    else
                        write(2,*) r, lin_dips(3,i,r_i)
                    end if
                end do
                write(2,*)
            end if
        end do
        close(2)

        open(2, file = "./dips_imag_exc.dat")
        do i = 1, num_channels
            if (quant_nums(i,1) >= 2 .and. quant_nums(i,2) == 0) then
                !write(2, *) quant_nums(i,:)
                do r_i = 1, interval_num + 1
                    r = r_min + (r_interval*(r_i-1))
                    if (aimag(dipole_data(1,i,r_i)) /= 0) then
                        write(2, *) r, aimag(dipole_data(1,i,r_i))
                    else if (aimag(dipole_data(2,i,r_i)) /= 0) then
                        write(2, *) r, aimag(dipole_data(2,i,r_i))
                    else if (aimag(dipole_data(3,i,r_i)) /= 0) then
                        write(2, *) r, aimag(dipole_data(3,i,r_i))
                    else
                        write(2,*) r, aimag(dipole_data(3,i,r_i))
                    end if
                end do
                write(2,*)
            end if
        end do
        close(2)


        open(2, file = "./dips_real_ground.dat")
        do i = 1, num_channels
            if (quant_nums(i,1) == 1 .and. quant_nums(i,2) == 0) then
                !write(2, *) quant_nums(i,:)
                do r_i = 1, interval_num + 1
                    r = r_min + (r_interval*(r_i-1))
                    if (real(dipole_data(1,i,r_i)) /= 0) then
                        write(2, *) r, real(dipole_data(1,i,r_i))
                    else if (real(dipole_data(2,i,r_i)) /= 0) then
                        write(2, *) r, real(dipole_data(2,i,r_i))
                    else if (real(dipole_data(3,i,r_i)) /= 0) then
                        write(2, *) r, real(dipole_data(3,i,r_i))
                    else
                         write(2,*) r, real(dipole_data(3,i,r_i))
                    end if
                end do
                write(2,*)
            end if
        end do
        close(2)

        open(2, file = "./dips_imag_ground.dat")
        do i = 1, num_channels
            if (quant_nums(i,1) == 1 .and. quant_nums(i,2) == 0) then
                !write(2, *) quant_nums(i,:)
                do r_i = 1, interval_num + 1
                    r = r_min + (r_interval*(r_i-1))
                    if (aimag(dipole_data(1,i,r_i)) /= 0) then
                        write(2, *) r, aimag(dipole_data(1,i,r_i))
                    else if (aimag(dipole_data(2,i,r_i)) /= 0) then
                        write(2, *) r, aimag(dipole_data(2,i,r_i))
                    else if (aimag(dipole_data(3,i,r_i)) /= 0) then
                        write(2, *) r, aimag(dipole_data(3,i,r_i))
                    else
                        write(2,*) r, aimag(dipole_data(3,i,r_i))
                    end if
                end do
                write(2,*)
            end if
        end do
        close(2)

    end subroutine write_dips

end module print_ops
