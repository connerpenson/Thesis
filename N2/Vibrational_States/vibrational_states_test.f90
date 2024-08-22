!************************************************************************************************!
         program main
         implicit none
         integer i_elec,num_elec,num_Q,num_v,i_Q,iv_index,tmp_index_of_resonances
         real*8, allocatable :: geom(:)
         complex*16, allocatable :: wf_on_geom(:,:,:),energies_v(:,:)
         Character*2 molecule

         num_v=9;num_Q=20;num_elec=1
         allocate(wf_on_geom(num_v,num_Q,num_elec),geom(num_Q),energies_v(num_v,num_elec))

         do i_Q=1,num_Q
            geom(i_Q)=1.8+(i_Q-1.)*(3-1.8)/(num_Q-1)
         enddo

         do i_elec=1,num_elec
         print*,'********************************************************************************'
         Call Program_of_LAC(molecule,num_Q,num_v,geom,energies_v(:,i_elec),&
               wf_on_geom(:,:,i_elec),i_elec,iv_index,tmp_index_of_resonances)
         enddo

         print*,'bound (0) or resonant states caculations (1)',tmp_index_of_resonances

         deallocate(wf_on_geom,geom,energies_v)
         end program main
!*************************************************************************************************!
         Real*8 function delt(i,j)
         integer i,j
         delt=0.
         If(i.eq.j)delt=1.
         end function delt
!*************************************************************************!
         include "vib_eigen_states.f90"
