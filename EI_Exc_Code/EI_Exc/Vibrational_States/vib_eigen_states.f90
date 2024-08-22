 Module MainMod
   SAVE
   Integer npmax,nsmax,ns,np,ij_cor,index_of_resonances
   Real*8 ryd,w,Pi,bet,au_amu,au_to_sec,au_eV
   Character*10 Name_of_in
   Real*8, Allocatable :: rs0(:),rs1(:),rs2(:),rs3(:),x(:),E_bv(:,:,:)
   Real*8, Allocatable ::  dR_all(:,:),R_all(:,:)
   Real*8 U_infinity_2,U_infinity_Elec1
   Complex*16, Allocatable :: VL(:,:),wf_all(:,:)
   integer i_elec
   real*8 wread,kmax
  Character*3 molecule
 End Module MainMod
! --------------------------------------
!             MAIN PROGRAM_of_LACx
! --------------------------------------
 Subroutine Program_of_LAC(molecule_Q,num_Q,num_v,geom,energies_v,&
    wf_geom,tmp_i_elec,iv_index,tmp_index_of_resonances)
   Use MainMod
   implicit none
   Integer j_max,j,nvj(100),jmin,jmax,jpas,i,i_change,num_v,num_Q,tmp_i_elec,iv_index
   Integer tmp_index_of_resonances
   Real*8 bv,tmp,dx,Rmin,Rmax,U_infinity,geom(num_Q)

   Character*2 atome
   Complex*16 energies_v(num_v), wf_geom(num_v,num_Q)
   Character*2 molecule_Q
   Character*3 cex
   Complex*16, parameter :: ci = (0.d0, 1.d0)

   npmax=5000 ; nsmax=1; j_max=1!81   ! This can be changed by the user

   Allocate(E_bv(j_max,npmax*nsmax,5))
   wf_geom=(0.d0,0.d0)
   energies_v=(0.d0,0.d0)
   E_bv=0.d0

   Ryd=219474.6313710
   pi=atan(1.d0)*4.
   au_amu=5.4857990946d-4
   au_to_sec=2.4188843e-17
   au_eV=27.2113834d0

   U_infinity=0.d0
   Name_of_in='BSCuFGR.in'
   Call FindStr(Name_of_in,5,'#initia')
   read(5,*)atome,Rmin,Rmax,ns,wread,U_infinity,index_of_resonances,molecule
   read(5,*)jmin,jmax,jpas
   Read(5,*)i_change,bet
   Close(5)

   molecule_Q=molecule
   i_elec=tmp_i_elec
   w=wread/au_amu

   print*,'molecule is ',molecule
   print*,'mass=',w,'in a.u.'

   Do ij_cor=jmin,jmax,jpas
     Print *,'J=',ij_cor
     Call BoundWithChGrid(nvj,num_v,num_Q,geom,wf_geom)
   Enddo

!  The main ouput file
   Open(22,Status='REPLACE',File='Target/Target_'//atome//'_'//molecule//'_ElecState'//&
         cex(i_elec)//'_E.dat')
   Do ij_cor=jmin,jmax,jpas
     Write(22,'(A)')''
     Write(22,'(A3,I4.1)')'#J= ',ij_cor
     Write(22,'(A90)')'#iv, ib, E_v (bound and diss. states for E_v<0 and E_v>0 resp. in eV), Bv, Gamma '
     Do i=1,nvj(ij_cor+1)
      Write(22,'(i6.2,20e20.10)')i-1,E_bv(ij_cor+1,i,1),E_bv(ij_cor+1,i,2)!*au_eV
      !Write(22,'(i6.2,20e20.10)')i-1,E_bv(ij_cor+1,i,1),(-U_infinity_Elec1+E_bv(ij_cor+1,i,2))*au_eV,&
      !E_bv(ij_cor+1,i,3)
       if(-U_infinity_2+E_bv(ij_cor+1,i,2)<=0)iv_index=i
     Enddo
   Enddo
   Close(22)

   Open(22,Status='REPLACE',File='Target/Target_'//atome//'_'//molecule//'_ElecState'//&
          cex(i_elec)//'_Transition_Frequencies.dat')
   Write(22,'(A10,2e20.12)')molecule,w,w*au_amu
   Do ij_cor=jmin,jmax,jpas
       Do i=2,nvj(ij_cor+1)
           Write(22,'(i6.2,A6,i6.2,e20.10)')i-1,'-->',i-2,(E_bv(ij_cor+1,i,2)-E_bv(ij_cor+1,i-1,2))*Ryd
       Enddo
   Enddo
   Close(22)

   Open(22,Status='REPLACE',File='Target/Target_'//atome//'_'//molecule//'_ElecState'//&
          cex(i_elec)//'_Beta_vs_E.dat')
   Do ij_cor=jmin,jmax,jpas
       Do i=1,nvj(ij_cor+1)
           !if(E_bv(ij_cor+1,i,2)<0.d0)
           Write(22,'(i6.2,e20.10)')i-1,E_bv(ij_cor+1,i,2)
       Enddo
   Enddo
   Close(22)

   ! data for the excitation program
    do i=1,num_v
     energies_v(i)=complex(E_bv(jmin+1,i,2),-E_bv(jmin+1,i,4)/2.0)   !modified by mehdi Oct. 20th 2019
   enddo
   tmp_index_of_resonances=index_of_resonances

   Deallocate(E_bv)
 End SUBROUTINE Program_of_LAC
! --------------------------------------------------------------------
!             Solution of SE for one partial wave
! --------------------------------------------------------------------
 Subroutine BoundWithChGrid(nvj,num_v,num_Q,geom,wf_geom)
   Use MainMod
   Integer i,k,iv,is,iwf,j,nvj(100),iff,num_v,num_Q
   Real*8 Rmin,Rmax,U_infinity,F(npmax),nl,nr,U_adiab(npmax,nsmax),geom(num_Q)
   Real*8, Allocatable ::  T(:,:),v_diag(:,:),T5(:,:),S(:,:,:),v(:,:,:),v_non_modified(:,:,:)
   Character*2 atome
   Character*3 cex
   Complex*16 wf_geom(num_v,num_Q),psi_on_geom(num_v,num_Q,nsmax)

   Allocate(rs0(npmax),rs1(npmax),rs2(npmax),rs3(npmax),x(npmax),S(npmax,nsmax,nsmax),v(npmax,nsmax,nsmax),&
      v_non_modified(npmax,nsmax,nsmax))

   Call FindStr(Name_of_in,5,'#initia')
   read(5,*)atome,Rmin,Rmax,ns,wread,U_infinity,index_of_resonances!,molecule
   Close(5)

   Call Read_V_and_ChangeGrid(v,v_non_modified,Rmin,Rmax, U_infinity)
   Call S_matrix(v_non_modified,S,U_adiab)

   Allocate(v_diag(np*ns,2),VL(np*ns,np*ns),T5(np*ns,np*ns),T(np,np))
   Call KinOperator_in_x_basis_Odd(np,x(1:np),rs1(1:np),T,w)
   v_diag=0.
   Call FindStr(Name_of_in,5,'#initwf')
   read(5,*)iv,k,iwf
   Close(5)
   If(index_of_resonances==0)Then
     print*,'Bound states calculations'
     Call Sum_and_Diag(v,T,v_diag(1:np*ns,1),T5)
     VL=T5
   ElseIF(index_of_resonances==1)Then
     print*,'resonance calculations'
     F=0.
     Call OptPot(np,rs0,F(1:np),Name_of_in)
!    Ouput of optical potential
    Open(21,Status='REPLACE',File='CAP'//'_ElecState'//cex(i_elec)//'.dat')
    Do i=1,np
      Write(21,'(10e20.12)')rs0(i),(max(F(i)+(-U_infinity_Elec1+U_infinity_2),1e-20))*au_eV ! in eV
    Enddo
    Close(21)
     Call Sum_and_Diag_C(v,T,F(1:np),v_diag,iwf)
   Else
     Stop 'Choose the right resonance index'
   Endif

   Call FindStr(Name_of_in,5,'#FC_Fac')
   read(5,*)iff
   Close(5)
   If(iff==1)Then
     Open(unit=20,Status='REPLACE',file='T5.bin',form='unformatted')
     write(20)w
     write(20)ns,np
     write(20)(x(i),rs0(i),rs1(i),i=1,np)
     write(20)(v_diag(i,1),i = 1,np*ns)
     write(20)((T5(i,j),i = 1,np*ns),j=1,np*ns)
     Close(20)
   Endif

    Call Bv_sort(v_diag,nvj(ij_cor+1))

   If(iwf==1)Call GoodWFunction(S,num_v)
   print*,'Wave functions plot done...'

   ! get Fourier Transform of the WF
   If(iwf==1.and.index_of_resonances==1)call GoodWFunction_FT(num_v)
   If(iwf==1.and.index_of_resonances==1)print*,'FT wave functions plot done...'

   ! data for the excitation program
   call GoodWFunction2(num_Q,num_v,geom,psi_on_geom)
   wf_geom(:,:)=psi_on_geom(:,:,1)
   print*,'wave functions for the excitation program done...'

   Deallocate(v,T5,rs0,rs1,rs2,rs3,x,S,v_diag,VL,T,v_non_modified)
   If(iwf==1)Deallocate(wf_all,R_all,dR_all)

 End Subroutine BoundWithChGrid
! --------------------------------------------------------------------
!             Calculation of a change of variable
! --------------------------------------------------------------------
 Subroutine Read_V_and_ChangeGrid(Umo_MO,v_non_modified,Rmin,Rmax,U_infinity)
   Use MainMod
   Integer i_change,i_non_ad, is,js,id,ncolumnx,ncolumnU,i,nd2(nsmax),&
           n_add,ii,j,N_sup, ndmax,nd,i0,nd0,i1,i2,Na
   Parameter(ndmax=35000)
   Real*8, Dimension (:,:), Allocatable ::Umo_tmp,Umo_tmp2,xtmp2
   Real*8, Dimension (:), Allocatable ::U_sup,cm,xtmp,xtmp3
   Real*8 Umo_MO(npmax,nsmax,nsmax),alpha,Ua,xm,dx,t,t3,r,v_non_modified(npmax,nsmax,nsmax),&
           U_coupl(npmax,nsmax,nsmax),dr_tmp,Rmin_tmp,Rmax_tmp,&
           spl,U_infinity,tmp,tmp2,Rmin,Rmax,E_inf,E_infinity
   Character*8 file_with_curve
   Character*3 cex
   Allocate (Umo_tmp(0:ndmax,nsmax),Umo_tmp2(ndmax,nsmax),&
          U_sup(ndmax),xtmp2(ndmax,nsmax),cm(ndmax), xtmp(ndmax),xtmp3(ndmax))
   Call FindStr(Name_of_in,5,'#initia')
   Read(5,*)
   Read(5,*)
   Read(5,*)i_change,bet,Na
   Close(5)
   Call FindStr(Name_of_in,5,'#E_infi')
   read(5,*)E_inf
   Close(5)

   Call FindStr(Name_of_in,5,'#Pot_di')
   xtmp2=0.
   Umo_tmp=0.

   Do is=1,ns
     Read(5,*)file_with_curve,ncolumnx, ncolumnU
   !print*,is,file_with_curve//cex(i_elec)//'.dat',i_elec
     Open(11,File=file_with_curve//cex(i_elec)//'.dat')
     Do id = 1,ndmax
       read(11,*,end=99) (xtmp2(id,is),i=1,ncolumnx),(Umo_tmp(id,is),i=ncolumnx+1,ncolumnU)
     Enddo
99   nd2(is)=id-1
     Close (11)
   Enddo
   Close(5)

    print *, file_with_curve//cex(i_elec)//'.dat'
   Rmin_tmp=xtmp2(1,1)
   Rmax_tmp=xtmp2(1,1)

   !Umo_tmp(1:nd2(1),1)=Umo_tmp(1:nd2(1),1)+0.02*exp(-(xtmp2(1:nd2(1),1)-7)**2/10.)

   U_infinity_2=Umo_tmp(nd2(1),1)
   if(i_elec==1)U_infinity_Elec1=Umo_tmp(nd2(1),1)

   Do is=1,ns
     Do i=1,nd2(is)
       Rmin_tmp=min(Rmin_tmp,xtmp2(i,is))
       Rmax_tmp=max(Rmax_tmp,xtmp2(i,is))
     Enddo
   Enddo

   dr_tmp=0.1
   nd=int((Rmax_tmp-Rmin_tmp)/dr_tmp)+1
   Do i=1,nd
     xtmp(i)=Rmin_tmp+(i-1)*dr_tmp
   Enddo


   print*,'U_infinity',U_infinity
   print*,'U_infinity_2',U_infinity_2
   print*,'U_infinity_Elec1=',U_infinity_Elec1
   Do is=1,ns
     call spline(nd2(is),xtmp2(1:nd2(is),is),Umo_tmp(1:nd2(is),is),cm(1:nd2(is)))
     Do i=1,nd
       r=xtmp(i)
       If (r.lt.xtmp2(1,is))Then
         Umo_tmp2(i,is)=Umo_tmp(1,is)
       Elseif(r.gt.xtmp2(nd2(is),is))Then
         Umo_tmp2(i,is)=Umo_tmp(nd2(is),is)
       Else
         Umo_tmp2(i,is)=spl(nd2(is),xtmp2(1:nd2(is),is),Umo_tmp(1:nd2(is),is),cm(1:nd2(is)),r)
       Endif
     Enddo
   Enddo
   Do is=1,ns
     Do i=1,nd
       Umo_tmp2(i,is)=Umo_tmp2(i,is)-U_infinity
     Enddo
   Enddo


   Umo_tmp=0.
   Do is=1,ns
     Do i=1,nd
       Umo_tmp(i,is)=Umo_tmp2(i,is)+ij_cor*(ij_cor+1.)/(2.*w*xtmp(i)*xtmp(i))
     Enddo
   Enddo


   Print *, 'i_elec',i_elec,Umo_tmp(nd,1),nd

   E_inf=E_inf+E_infinity(nd,Umo_tmp(1:nd,1:ns),xtmp(1:nd))
   Print *,'E_inf=',E_inf
   Rmin=max(Rmin_tmp,Rmin)
   Call U_sup_calc(nd,xtmp(1:nd),N_sup,xtmp3(1:nd),Umo_tmp(1:nd,1:ns),U_sup(1:nd),E_inf)
   If(i_change==0)Then
     Call ChangeGrid0(Umo_tmp,xtmp,nd,ndmax, Rmin,Rmax,E_inf)
   Elseif(i_change==2)Then
     Call ChangeGrid2(xtmp,Rmin,Rmax,U_sup,N_sup,xtmp3,Na)
   Else
     Stop 'Choose the correct coordinate transformation'
   Endif

   call spline(np,rs0(1:np),x(1:np),cm(1:np))
   xtmp3=0.
   Do i=1,nd
     xtmp3(i)=spl(np,rs0(1:np),x(1:np),cm(1:np),xtmp(i))
   Enddo
   Umo_MO=0.
   Do is=1,ns
     Do i=2,nd
       If(Umo_tmp2(i,is).ne.Umo_tmp2(i-1,is))Exit
     Enddo
     i0=i-1
     nd0=nd-i0+1
     call spline(nd0,xtmp3(i0:nd),Umo_tmp2(i0:nd,is),cm(i0:nd))
     Do i=1,np
       r=x(i)
       tmp=ij_cor*(ij_cor+1.)/(2.*w*rs0(i)*rs0(i))
       If(r.lt.xtmp3(i0))Then
         Umo_MO(i,is,is)=Umo_tmp2(1,is)
       Elseif(r.gt.xtmp3(nd))Then
         Umo_MO(i,is,is)=Umo_tmp2(nd,is)+tmp
       Else
         Umo_MO(i,is,is)=spl(nd0,xtmp3(i0:nd), Umo_tmp2(i0:nd,is),cm(i0:nd),r)+tmp
       Endif
     Enddo
   Enddo
   DeAllocate (Umo_tmp,Umo_tmp2,U_sup,xtmp2,cm,xtmp,xtmp3)

   Call FS_on_R(U_coupl)
   Umo_MO(1:np,1:ns,1:ns)=Umo_MO(1:np,1:ns,1:ns)+U_coupl(1:np,1:ns,1:ns)

!  Output of potentials and functions R,J,J',J''
   Open (9,Status='REPLACE',file='Target/'//molecule//'_ElecState'//cex(i_elec)//'_in_calc.dat')
   write(9,'(A60,e20.12,A8)')'  #energies calculated with respect to the asymptotic limit :',U_infinity_Elec1*au_eV,&
          '  in eV'
   Do ii=1,np-1   !mehdi comments: see line 296 such as np=2*k+1 coming from the Jacobian transformation of the variable grid
     write(9,'(17(e20.12))')rs0(ii),((-U_infinity_Elec1+Umo_MO(ii,i,i))*au_eV, i=1,ns) !modified in Oct 24th 2019
   Enddo

   Close(9)
!    Open(13,Status='REPLACE',file='v_coupl.dat')
!    Do ii=1,np
!      write(13,'(17(e20.12))')rs0(ii),((Umo_MO(ii,i,j), i=1,ns),j=1,ns)
!    Enddo
!    Close(13)
!    Open(12,Status='REPLACE',file='rs.dat')
!    Do i=1,np
!      write(12,'(9(f20.12))') x(i),rs0(i),rs1(i),rs2(i),rs3(i)
!    Enddo
!    Close(12)

   v_non_modified=Umo_MO
!  Calculation of the modified potential (see thesis)
   do j = 1,ns ; do i = 1,np
     t=1./(rs1(i)); t3=t*t*t
     Umo_MO(i,j,j)=Umo_MO(i,j,j)+1./(2.*w)*(7./4.*t3*t*rs2(i)*rs2(i)-1./2.*t3*rs3(i) )
   enddo; enddo

 End Subroutine Read_V_and_ChangeGrid
! --------------------------------------------------------------------
!             Calculation of a change of variable N0 : constant step
! --------------------------------------------------------------------
 Subroutine ChangeGrid0(Umo_tmp,xtmp,nd,ndmax,Rmin,Rmax,E_inf)
   Use MainMod
   Integer nd,ndmax,i,is,k
   Real*8 Rmin,Rmax,Umo_tmp(ndmax,nsmax),xtmp(ndmax),U_sup(nd),Umax,Umin,E_inf
   Umax=Umo_tmp(nd,1)
   Umin=Umo_tmp(nd,1)
   Do is=1,ns
     Umax=min(Umo_tmp(nd,is),Umax) !It needs to be changed
     Do i=1,nd
       If((xtmp(i).ge.Rmin).and.(xtmp(i).le.Rmax))Umin=min(Umo_tmp(i,is),Umin)
     Enddo
   Enddo
   Umax=E_inf-Umin+Umax
   rs1(1)=bet*pi/(sqrt(2*w*Umax))
   k = (Rmax-Rmin)/rs1(1)+1
   k=(k+1)/2
   print *, Umax
   np=2*k+1
   rs1(2:np)=rs1(1)
   Print *,'Number of grid poins, np=',np
   If (np.gt.npmax)Stop 'Too many points'

   Do i=1,np
     x(i)=i
     rs0(i)=Rmin+(i-1)*rs1(1)
     rs2(i)=0.
     rs3(i)=0.
   Enddo
 End Subroutine ChangeGrid0
! --------------------------------------------------------------------
!             Calculation of a change of variable N2 :
!             variable step adapted to the potential.
! --------------------------------------------------------------------
 Subroutine ChangeGrid2(xtmp,Rmin,Rmax,U_sup,N_sup,x_sup,Na)
   Use MainMod
   Integer i,Na,k,N_sup
   Real*8 alf,r_tmp,tmp,Rmin,Rmax,xtmp(1:N_sup),U_sup(1:N_sup),cm(3*N_sup),spl,x_sup(1:N_sup)
   alf=dsqrt(2*w)/pi
   call spline(N_sup,x_sup(1:N_sup),U_sup(1:N_sup),cm(1:N_sup))
   tmp=dabs(spl(N_sup,x_sup(1:N_sup),U_sup(1:N_sup),cm(1:N_sup),Rmin))
   r_tmp=Rmin
   k=0
1  Do While(r_tmp<Rmax)
     tmp=dabs(spl(N_sup,x_sup(1:N_sup),U_sup(1:N_sup),cm(1:N_sup),r_tmp))
     If (r_tmp.ge.x_sup(N_sup))Then
       rs1(k+1)=rs1(k)
     Else
       rs1(k+1)=bet/(alf*dsqrt(tmp))
     Endif
     r_tmp=r_tmp+rs1(k+1)
     k=k+1
   Enddo
   If(2*(k/2).eq.k)k=k-1
   If(k.gt.npmax)Then
     Print *, 'Attention: the grid is very large, np=',k
     Stop
   Endif
   np=k+Na
   rs0(1)=Rmin
   Print *,'Number of grid poins, np=',np
   Call Der_Odd(np,rs0,rs1,rs2,rs3,Na)
   Do i=1,np
     x(i)=i
   Enddo
 End Subroutine ChangeGrid2
! --------------------------------------------------------------------
!             Calculation of the enveloping potential.
! --------------------------------------------------------------------
 Subroutine U_sup_calc(nd,R,N,R_sup,Umo_tmp,U_sup,E_inf)
   Use MainMod
   Integer nd,i,is,j,N
   Real*8 R(1:nd),U_sup(1:nd),Umo_tmp(1:nd,1:ns),R_sup(1:nd),Umo_tmp2(nd,ns),tmp,E_inf,x_loc,dx,U_sup2(nd),C,Del

   Do is=1,ns
     Do i=1,nd
       Umo_tmp2(i,is)=Umo_tmp(i,is)-Umo_tmp(nd,is)
     Enddo
   Enddo
   U_sup(nd)=0.
   Do i=nd-1,1,-1
     tmp=Umo_tmp2(i,1)
     Do is=1,ns
       tmp=min(tmp,Umo_tmp2(i,is))
     Enddo
     U_sup(i)=min(U_sup(i+1),tmp)
   Enddo
   U_sup(1:nd)=U_sup(1:nd)-E_inf


   U_sup2=0.
   dx=0.3
   j=1
   N=(R(nd)-R(1))/dx-1
   If(N>nd)Stop 'The number of points for U_sup is too big.'
   Do i=1,N
     R_sup(i)=R(1)+(i-1)*dx
 1   If(R_sup(i).ge.R(j+1))Then
       j=j+1
       Goto 1
     Endif
     U_sup2(i)=U_sup(j)
     If(((R(j)+R(j+1))/2.).gt.R_sup(i))U_sup2(i)=U_sup(j+1)
   Enddo

   U_sup(1:N)=0.
   Del=1.d0
   C = 1.d0/(Del*sqrt(2.d0*asin(1.d0)))
   Do i=1,N
     Do j=1,N
       x_loc=abs(R_sup(i)-R_sup(j))
       If(x_loc<4.*Del)U_sup(i)=U_sup(i)+U_sup2(j)*exp(-((x_loc/Del)**2))*dx
     Enddo
     Do j=1,N
       x_loc=abs(R_sup(i)-(R_sup(1)-j*dx))
       If(x_loc<4.*Del) U_sup(i)=U_sup(i)+U_sup2(1)*exp(-((x_loc/Del)**2))*dx
       x_loc=abs(R_sup(i)-(R_sup(N)+j*dx))
       If(x_loc<4.*Del) U_sup(i)=U_sup(i)+U_sup2(N)*exp(-((x_loc/Del)**2))*dx
     Enddo
   Enddo

   U_sup(1:N)=C*U_sup(1:N)
!  Output of the supporting potential
!   Open (10,Status='REPLACE',file='U_sup.dat')
!   Do i=1,N
!     write(10,'(9(e20.12))') R_sup(i),U_sup(i)
!   Enddo
!   Close(10)

 End Subroutine U_sup_calc
! --------------------------------------------------------------------
!             Calculation of the asymptotic energy needed
!             to be represented by the mapping at large distances
! --------------------------------------------------------------------
 Real*8 Function E_infinity(nd,Umo_tmp,xtmp)
   Use MainMod
   implicit none
   Integer nd,is,i,Ml(1)
   Real*8 Umo_tmp(1:nd,1:ns),xtmp(1:nd),tmp1,tmp2

   tmp1=Umo_tmp(nd,1)
   tmp2=Umo_tmp(nd,1)
   Do is=1,ns
     Ml=Minloc(Umo_tmp(1:nd,is))
     tmp1=min(Umo_tmp(nd,is),tmp1)
!    If(index_of_resonances==1)Then
     Do i=nd-1,Ml(1),-1
       If((Umo_tmp(i,is).gt.Umo_tmp(i-1,is)).and. (Umo_tmp(i,is).gt.Umo_tmp(i+1,is)))tmp2=max(tmp2,Umo_tmp(i,is))
     Enddo
!    Endif
   Enddo
   E_infinity=max(tmp2-tmp1,0.d0)

 End Function E_infinity
! --------------------------------------------------------------------
!             This procedure applies a smooth periodic
!             boundary condition for the calculated mapping function
! --------------------------------------------------------------------
 Subroutine Der_Odd(N,rs0,rs1,rs2,rs3,Na)
   Integer N,i,Na
   Real*8 rs0(N),rs1(N),rs2(N),rs3(N),dx,r_tmp
   Complex*16 r_p(N),ak(N),r(N),tmp(N),tmp2(N)
   dx=1.
   r_tmp=rs0(1)
   rs1(N-Na+1:N)=rs1(1)
   Call GTrans(N,rs1,Na)
   r=rs1
   Call FT(N,r,r_p,-1)
   i=N/Na
   r_p((N-1)/2-i:(N-1)/2+i)=0.

   ak=1.
   tmp=r_p*ak/N
   Call FT(N,tmp,tmp2,1)
   rs1=tmp2

   Call Initak(ak,N,dx,1)
   tmp=r_p*ak/N
   Call FT(N,tmp,tmp2,1)
   rs2=tmp2

   Call Initak(ak,N,dx,2)
   tmp=r_p*ak/N
   Call FT(N,tmp,tmp2,1)
   rs3=tmp2

   Call Initak(ak,N,dx,-1)
   tmp=r_p*ak/N
   Call FT(N,tmp,tmp2,1)
   Do i=1,N
     rs0(i)=tmp2(i)-tmp2(1)+r_p(N)/N*(i-1)+r_tmp
   Enddo
 End Subroutine Der_Odd
! --------------------------------------------
 Subroutine GTrans(N,r,Na)
   Integer N,Na,i,j
   Real*8 r(N),rt(-N:2*N),pi,norm
   pi = 4.d0*datan(1.d0)
   norm=4./(Na*sqrt(pi))
   rt(1:N)=r(1:N)
   rt(-N+1:0)=r(1:N)
   rt(N+1:2*N)=r(1:N)
   r=0.
   Do i=1,N
     Do j=i-Na,i+Na
       r(i)=r(i)+rt(j)*exp(-(abs(4.*(j-i)/Na))**2)
     Enddo
   Enddo
   r=r*norm
 End Subroutine GTrans
! --------------------------------------------------------------------
!      To be adapted for a given application
!      This procedure defines non-diagonal couplings between potential curves
! --------------------------------------------------------------------
 Subroutine FS_on_R(U_coupl)
   Use MainMod
   Integer i
   Real*8 c12,c13,c23,FS,U_coupl(npmax,nsmax,nsmax)
   U_coupl=0.
   Select Case (ns)
   Case(1)
     U_coupl=0.
   Case(2)
     Call FindStr(Name_of_in,5,'#coupl2')
     read(5,*)FS
     Close(5)
     U_coupl(:,1,1)=-FS/3.
     U_coupl(:,1,2)=sqrt(2.)*FS/3.
     U_coupl(:,2,1)=U_coupl(:,1,2)
     U_coupl(:,2,2)=0.
   Case(3)
     Call FindStr(Name_of_in,5,'#coupl3')
     read(5,*)c12,c13,c23
     Close(5)
     U_coupl(:,1,2)=c12/2.
     U_coupl(:,2,1)=c12/2.
     U_coupl(:,1,3)=c13
     U_coupl(:,3,1)=c13
     Do i=1,np
       U_coupl(i,2,3)=c23*Sqrt(ij_cor*(ij_cor+1.))/(w*rs0(i)*rs0(i))
       U_coupl(i,3,2)=U_coupl(i,2,3)
!      States 1 and 2 are \Pi states, therefore:
       U_coupl(i,1,1)=-1./(2.*w*rs0(i)*rs0(i))
       U_coupl(i,2,2)=-1./(2.*w*rs0(i)*rs0(i))
     Enddo
   End Select
 End Subroutine FS_on_R
! --------------------------------------------------------------------
!      Discrete Fourier Transform
! --------------------------------------------------------------------
 Subroutine FT(N,f,f_t,isign)
   Integer N,isign,k,j
   Complex*16 a,f(N),f_t(N),ci
   ci=(0.d0,1.d0)
   a = 8.*atan(1.d0)*ci/N
   If(isign.lt.0)a=-1.*a
   f_t=0.
   Do j=1,N
     Do k=1,N
       f_t(j)=f_t(j)+f(k)*exp(a*k*j)
     Enddo
   Enddo
 End Subroutine FT
! --------------------------------------------------------------------
!     Auxiliary procedure needed for the Fourier differentiation
!                       of the order "iorder"
! --------------------------------------------------------------------
 Subroutine Initak(ak,n,delta,iorder)
   Integer n,j,iord4,iorder,iord2
   Real*8 delta,pi,sign2,anorm
   Complex*16 ak(n),c

   pi = 4.*atan(1.d0)
   anorm = 2.*pi/(n*delta)
   c=(0.,1.)

   sign2 = (-1)**iorder
   do j = 1,n/2
     ak(j) = (c*j*anorm)**iorder
   enddo
   do j = n/2+1,n-1
     ak(j) = sign2*(c*(n-j)*anorm)**iorder
   enddo
   ak(n) = 0.
 End Subroutine Initak
! -----------------------------------------------------------
!     The matrix of the kinetic energy (-1/(2m)*d^2/dR^2)
! -----------------------------------------------------------
 Subroutine KinOperator_in_x_basis_Odd(np,x,rs1,T,w)
!  Only for odd number of grid points (np)
   Integer i,j,np
   Real*8 x(1:np),rs1(1:np),w,T(np,np),L,coef,coefdiag,s,ss,pi
   T=0.
   pi = 2.*asin(1.d0)
   L=(x(2)-x(1))*np
   coef=pi*pi/(w*L*L)
   coefdiag=coef*(np*np-1.)/6.
   Do i=1,np
     T(i,i)=coefdiag/(rs1(i)*rs1(i))
     Do j=i+1,np
       s=dcos((i-j)*pi/np)/(dsin((i-j)*pi/np))**2
       ss=1./(rs1(i)*rs1(i))+1./(rs1(j)*rs1(j))
       T(i,j)=(-1)**(i-j)*coef/2.*ss*s
       T(j,i)=T(i,j)
     Enddo
   Enddo

  !   L=rs1(1)*(x(2)-x(1))*np !added by Mehdi 18/12/2019
  !
  !    do i=1,np
  !      do j=1,np
  !         if(i==j)then
  !          T(i,j)=pi**2/2.*1./(2*w*L**2)*(2*(np+1)**2/3.+1./3.-1./sin(pi*i/(np+1))**2)
  !
  !          else
  !          T(i,j)=(-1)**(i-j)*pi**2/2.*1./(2*w*L**2)*(1./sin(pi*(i-j)/(2*(np+1)))**2-1./sin(pi*(i+j)/(2*(np+1)))**2)
  !         endif
  !      enddo
  !   enddo

 End Subroutine KinOperator_in_x_basis_Odd
! -----------------------------------------------------------
!     The construction and diagonalization of the matrix H=T+V.
!     H is hermitian.
! -----------------------------------------------------------
 Subroutine Sum_and_Diag(v,T,v_diag,T5)
   Use MainMod
   Integer i,j,Iwork(5*np*ns),info,neig,ifail(np*ns),l,ia,ib,m,lwork
   Real*8 T(np,np),v(npmax,nsmax,nsmax),v_diag(np*ns),vll,vu,work(16*np*ns),abstol,delt
   Real*8 T4(np*ns,np*ns),T5(np*ns,np*ns)

   T4=0.

   Do i=1,ns
     Do l=1,np
       ia=(i-1)*np+l
       Do j=1,ns
         Do m=1,np
           ib=(j-1)*np+m
           T4(ib,ia)=T(l,m)*delt(i,j)+v(l,i,j)*delt(l,m)
         Enddo
       Enddo
     Enddo
   Enddo

!  diagonalisation
   lwork=8*ns*np
   abstol = 1.d-12
   neig = 0
   info=1
   call dsyevx('v','i','u',np*ns,T4,np*ns,vll,vu,1,np*ns,&
          abstol,neig,v_diag,T5,np*ns,work, lwork,iwork,ifail,info)
   print*,'diagonalization done...'

 End  Subroutine Sum_and_Diag
! -----------------------------------------------------------
!     The construction and diagonalization of the matrix H=T+V.
!     H is NOT  hermitian.
! -----------------------------------------------------------
 Subroutine Sum_and_Diag_C(v,T,F,v_diag,in_v)
   Use MainMod
   Integer i,j,info,l,ia,ib,m,in_v,ilo,ihi,   lwork,is,iv
   Real*8 T(np,np),v(npmax,nsmax,nsmax),nl,nr,delt,ABNRM,F(np),RCONDE(np*ns),RCONDV(np*ns),&
      RWORK(2*np*ns),SCALE(np*ns),v_diag(np*ns,2)
   Complex*16 E_v(np*ns),ci
   Character*1 jobvl,jobvr
   Complex*16, Dimension (:,:), Allocatable ::CT,VR
   Complex*16, Dimension (:), Allocatable ::WORK
   Allocate (CT(np*ns,np*ns),VR(np*ns,np*ns),WORK(2*np*ns*np*ns))

   CT=0.
   ci=(0.,1.)
   Do i=1,ns
     Do l=1,np
       ia=(i-1)*np+l
       Do j=1,ns
         Do m=1,np
           ib=(j-1)*np+m
           CT(ib,ia)=T(l,m)*delt(i,j)+v(l,i,j)*delt(l,m)
         Enddo
       Enddo
       CT(ia,ia)=CT(ia,ia)-ci*F(l)
     Enddo
   Enddo
!  diagonalisation
   lwork=2*ns*np*ns*np
   If(in_v.eq.1)Then
     jobvl='V'
     jobvr='V'
   Else
     jobvl='N'
     jobvr='N'
   Endif

   Call zgeevx('N',jobvl,jobvr,'N',ns*np,CT,np*ns,E_v,VL,np*ns,VR,np*ns,ilo,ihi,SCALE,ABNRM, RCONDE, RCONDV,WORK,lwork,RWORK,info)
   If(info.ne.0)Print *,'Attention, info=',info
   print*,'diagonalization done...'

   Do i=1,np*ns
     v_diag(i,1)=real(E_v(i))
     v_diag(i,2)=dimag(E_v(i))
   Enddo
   DeAllocate(WORK,VR,CT)
 End Subroutine Sum_and_Diag_C
! -----------------------------------------------------------
!     Optical potential. See works by
!     Vibok et al. (J. Phys. Chem. 96, 8712 (1992)) for details and parameters
! -----------------------------------------------------------
 Subroutine OptPot(np,rs0,F,Name_of_in)
   Integer np,i,CAP
   Real*8 rs0(np),L,A5,p,r0,F(np)
   Character*10 Name_of_in
   Call FindStr(Name_of_in,5,'#OptPot')
   read(5,*)CAP,L,A5
   Close(5)
   p=13.22
   r0=rs0(np)-L
   F=0.
   Do i=1,np
     If((rs0(i).ge.r0).and.(rs0(i).le.r0+L))then
         If(CAP==0)then
            F(i)=A5*p*exp(-2.*L/(rs0(i)-r0))
         ElseIF(CAP==1)then
            F(i)=A5*3./2.*((rs0(i)-r0)/L)**2
         EndIf
     Endif
   Enddo
 End Subroutine OptPot

! -----------------------------------------------------------
!     Here the block of output starts
! -----------------------------------------------------------

! -----------------------------------------------------------
!     The procedure calculates and writes wavefunctions selected
!     by the user in the input file 'BSCuFGR.in'
! -----------------------------------------------------------
 Subroutine GoodWFunction(S,num_v)
   Use MainMod
   Integer N_g,iv,i,ii,js,is,j,k,N_g_max,num_v
   Parameter(N_g_max=2000)
   Real*8 psi(np,ns),dx,q,q_tmp,r0,r1,cm1(np),E(1:np*ns),cm(np),spl,norme
   Real*8 tmp_r0,S(npmax,nsmax,nsmax),cmS(np,nsmax,nsmax)
   Complex*16 psi_g(nsmax), psi_psi(nsmax)
   Character*3 cex
   Character*20 blockname
   tmp_r0=0.

   allocate(wf_all(N_g_max,num_v),R_all(N_g_max,num_v),dR_all(N_g_max,num_v))
   wf_all=(0.,0.);R_all=0.;dR_all=0.

   Open(5,File=Name_of_in)
   do ii = 1,1000
     read(5,'(A7)',end=99)blockname
     If (blockname.eq.'#initwf')Then
       read(5,*)iv,N_g
       print*,'i_elec',i_elec,'iv',iv,' from ', Name_of_in,' file'

      if(iv+1<=(num_v))Open(7,Status='REPLACE',File='Target/Target_'//'wf'//cex(iv)//'_'//molecule//&
             '_ElecState'//cex(i_elec)//'.dat')
       dx=x(2)-x(1)
       call spline(np,x(1:np),rs0(1:np),cm(1:np))
       call spline(np,x(1:np),rs1(1:np),cm1(1:np))
       Do i=1,ns
         Do j=1,ns
           call spline(np,x(1:np),S(1:np,i,j),cmS(1:np,i,j))
         Enddo
       Enddo
       Do k=1,N_g
         q=x(1)+(k-1)*(x(np)-x(1))/(N_g-1)
         norme=0.
         Do is=1,ns
           psi_g(is)=0.
           Do j=1,np
             q_tmp=pi/dx*(q-x(j))
             If (dabs(q_tmp).le.1.e-7)then
               psi_g(is)=psi_g(is)+VL(j+(is-1)*np,iv+1)
             Else
               psi_g(is)=psi_g(is)+ VL(j+(is-1)*np,iv+1)*dsin(q_tmp)/(q_tmp)
             Endif
           Enddo
         Enddo
!        This block can be uncommented if one needs to
!        have wavefunctions in the adiabatic representation
!         Do is=1,ns
!           psi_psi(is)=0.
!           Do js=1,ns
!             psi_psi(is)= psi_psi(is)+ psi_g(js)*spl(np,x(1),S(1,is,js),cmS(1,is,js),q)
!           Enddo
!         Enddo
!         psi_g=psi_psi


         r1=spl(np,x(1:np),rs1(1:np),cm1(1:np),q)
         Do is=1,ns
           psi_g(is)=psi_g(is)/dsqrt(r1)
           norme=norme+abs(psi_g(is))**2
         Enddo
         r0=spl(np,x(1:np),rs0(1:np),cm(1:np),q)

!        save data for Fourier transform
         if(iv+1<=(num_v))then
            wf_all(k,iv+1)=psi_g(1)   ! added by mehdi Oct. 20th  2019
	    R_all(k,iv+1)=r0
!           build dR grid
            if(k>1)dR_all(k,iv+1)=R_all(k,iv+1)-R_all(k-1,iv+1)
         endif

         if(iv+1<=(num_v))Write(7,'(20(e20.12))')r0,(psi_g(is),is=1,ns),norme,R_all(k,iv+1),dR_all(k,iv+1)

       Enddo

       if(iv+1<=(num_v))dR_all(1,iv+1)=dR_all(2,iv+1)
       Close(7)
     Endif
   enddo
99 close(5)

 End Subroutine GoodWFunction
!**************************************************************************************************!
Subroutine GoodWFunction_FT(num_v)
   Use MainMod
   Integer N_g,iv,i,Nk,ik,ii,num_v
   Real*8 k
   Complex*16 FT_wf,ci
   Character*3 cex
   Character*20 blockname

   Nk=2000
   ci=(0.d0,1.d0)

   Open(5,File=Name_of_in)
   do ii = 1,1000
     read(5,'(A7)',end=99)blockname
     If (blockname.eq.'#initwf')Then
         read(5,*)iv,N_g

         if(iv+1<=(num_v))then

!        k-space grid
         kmax=pi/minval(dR_all(1:N_g,iv+1))
         print*,'i_elec=',i_elec,' iv= ',iv,' kmax=',kmax

	 Open(7,Status='REPLACE',File='Target/Target_'//'FT_wf'//cex(iv)//'_'//molecule//&
           '_ElecState'//cex(i_elec)//'.dat')
         do ik=1,Nk
             ! k-space grid
             k=-kmax+2*kmax*(ik-1)/(Nk-1)
             FT_wf=(0.d0,0.d0)
            ! The Normalized Fourier Transform WF
!           Fourier transfrom of psi(x) (== \bar{\psi}=int exp(-ikr)psi(r)dr/sqrt(2pi))
             do i=1,N_g
               FT_wf=FT_wf+exp(-ci*R_all(i,iv+1)*k)*wf_all(i,iv+1)*dR_all(i,iv+1)/sqrt(2*pi)
             enddo
	    Write(7,'(17e20.12)')k,abs(FT_wf)**2,FT_wf
	 enddo
         Close(7)
         endif
      Endif
   enddo
99 close(5)
 End Subroutine GoodWFunction_FT
! -----------------------------------------------------------
!     Another procedure that calculates and writes wavefunctions selected
!     by the user in the list of parameters
! -----------------------------------------------------------
 Subroutine GoodWFunction2(num_Q,num_v,geom,psi)
   Use MainMod
   Integer N_g,iv,i,js,is,j,k,num_v,num_Q
   Real*8 dx,q,q_tmp,r1,cm1(np),E(1:np*ns),cm(np),spl,norme,geom(num_Q)
   Complex*16 psi_g(nsmax), psi(num_v,num_Q,nsmax)

   do iv = 1,num_v
    dx=x(2)-x(1)
     call spline(np,rs0(1:np),x(1:np),cm(1:np))
     call spline(np,x(1:np),rs1(1:np),cm1(1:np))
     Do k=1,num_Q
       q=spl(np,rs0(1:np),x(1:np),cm(1:np),geom(k))
       norme=0.d0
       psi_g=(0.d0,0.d0)
       Do is=1,ns
         Do j=1,np
           q_tmp=pi/dx*(q-x(j))
           If (abs(q_tmp).le.1.e-7)then
             psi_g(is)=psi_g(is)+VL(j+(is-1)*np,iv) ! modified by mehdi 19th august 2018 iv+1--> iv
           Else
             psi_g(is)=psi_g(is)+ VL(j+(is-1)*np,iv)*sin(q_tmp)/(q_tmp) ! modified by mehdi 19th august 2018 iv+1--> iv
           Endif
         Enddo
       Enddo
       r1=spl(np,x(1:np),rs1(1:np),cm1(1:np),q)
       Do is=1,ns
         psi_g(is)=psi_g(is)/sqrt(r1)
         norme=norme+abs(psi_g(is))**2
       Enddo
       psi(iv,k,1:ns)=psi_g(1:ns)
     Enddo
   enddo

 End Subroutine GoodWFunction2

! -----------------------------------------------------------
!     This procedure calculates the adiabatic potential curves
!     by a simple diagonalization of the total potential for
!     each internuclear distance.
!     The procedure can be used for the representation of wavefunctions
!     in the adiabatic basis (if this basis is preferable)
! -----------------------------------------------------------
 Subroutine S_matrix(v,S,U_adiab)
   Use MainMod
   Integer lwork,neig,info,i,is,js,ind,i0,k,Iwork(5*ns),ifail(ns)
   Real*8 v(npmax,nsmax,nsmax),S(npmax,nsmax,nsmax),U_adiab(npmax,nsmax),vtmp(ns,ns),v_diag(ns),&
     vll,vu,work(8*ns),abstol,stmp(ns,ns),dS(np,ns,ns),t1(np,ns,ns),t2(np,ns,ns),gam,tmp,cm(np),spl
   lwork=8*ns; abstol = 1.d-12;  neig = 0
   Do i=1,np
     Do is=1,ns
       vtmp(is,is)=v(i,is,is)
       Do js=is+1,ns
         vtmp(is,js)=v(i,is,js)
         vtmp(js,is)=vtmp(is,js)
       Enddo
     Enddo
     call dsyevx('v','i','u',ns,vtmp,ns,vll,vu,1,ns,abstol,neig,v_diag,stmp,ns,work,lwork,iwork,ifail,info)
     Do is=1,ns
       U_adiab(i,is)=v_diag(is)
       If (stmp(is,is).ge.0.)Then
         ind=1
       Else
         ind=-1
       Endif
       Do js=1,ns
         S(i,is,js)=ind*stmp(js,is)
       Enddo
     Enddo
   Enddo

!  Ouput of adiabatic potential curves
!    Open(11, Status='REPLACE',file='v_adiab.dat')
!    Do i=1,np
!      Write(11,'(10(e20.12))')rs0(i),(U_adiab(i,is),is=1,ns)
!    Enddo
!    Close(11)
 End Subroutine S_matrix

! -----------------------------------------------------------
!     Calculation of different parameters using calculated eigen-energies
!     and eige-wavefunctions.
! -----------------------------------------------------------
 Subroutine Bv_sort(v_diag,nvj)
   Use MainMod
   Integer j_max,ib,i_loc,i,j,ind(np*ns),k,l,i_vac,nvj,comp_of_interest,i_sort
   Real*8 v_diag(np*ns,2),bv,tmp,dx,v_cut,R_cut,Bv_cut,ge_cut,au_to_k,au_to_Ghz,au_ps
   Complex*16 vl_t(np*ns)

!  Calculation of different parameters
   au_ps=2.4189d-05
   au_to_Ghz=30.*Ryd
   dx=x(2)-x(1)
   Do ib=1,np*ns
!    Calculation of rotational constants B_v=<i|1/(2mu*R^2)|i>
     bv=0.
     Do i=1,np*ns
       i_loc=Mod(i-1,np)+1
       bv=bv+(( abs(VL(i,ib)) /rs0(i_loc))**2)*dx
     Enddo
     bv=bv/(2.*w)
     E_bv(ij_cor+1,ib,1)=ib
     E_bv(ij_cor+1,ib,2)=v_diag(ib,1)                   ! Eigenenergy in hartree, get the real part
     E_bv(ij_cor+1,ib,3)=bv*Ryd                         ! B_v in cm-1
     !E_bv(ij_cor+1,ib,4)=au_to_Ghz*abs(2.*v_diag(ib,2)) ! Widhs in GHZ
     E_bv(ij_cor+1,ib,4)=abs(2.0*v_diag(ib,2)) ! Widhs in a.u. modified by mehdi Oct. é0th 2019

!    Population of a particular component of wf
     bv=0.
     comp_of_interest=1
     Do i=np*(comp_of_interest-1)+1,np*comp_of_interest
       bv=bv+abs(VL(i,ib))**2*dx
     Enddo
     E_bv(ij_cor+1,ib,5)=bv
   Enddo

!  Sorting according to the increasing energy
   Do i=1,np*ns
     ind(i)=i
   Enddo
   Do k=1,np*ns
     Do i=k+1,np*ns
       If(E_bv(ij_cor+1,i,2).lt.E_bv(ij_cor+1,k,2))Then
         j=ind(i)
         ind(i)=ind(k)
         ind(k)=j
         vl_t=VL(:,k)
         VL(:,k)=VL(:,i)
         VL(:,i)=vl_t
         Do j=1,5
           tmp=E_bv(ij_cor+1,k,j)
           E_bv(ij_cor+1,k,j)=E_bv(ij_cor+1,i,j)
           E_bv(ij_cor+1,i,j)=tmp
         Enddo
       Endif
     Enddo
   Enddo

!  Whether or not one wants to sort over some additional parameter
   Call FindStr(Name_of_in,5,'#Sortir')
   read(5,*)i_sort
   If(i_sort/=0)Then
     read(5,*)v_cut,R_cut,Bv_cut,ge_cut
     v_cut=v_cut+Ryd*ij_cor*(ij_cor+1.)/(2.*w*R_cut*R_cut)
     i_vac=1
     Do i=1,np*ns
       If((E_bv(ij_cor+1,i,2)<v_cut).and.(E_bv(ij_cor+1,i,3)>Bv_cut)&
            .and.(E_bv(ij_cor+1,i,4)<ge_cut).and.E_bv(ij_cor+1,i,5)>-1.)then
         vl_t=VL(:,i_vac)
         VL(:,i_vac)=VL(:,i)
         VL(:,i)=vl_t
         Do l=1,5
           tmp=E_bv(ij_cor+1,i_vac,l)
           E_bv(ij_cor+1,i_vac,l)=E_bv(ij_cor+1,i,l)
           E_bv(ij_cor+1,i,l)=tmp
         Enddo
         i_vac=i_vac+1
       Endif
     Enddo
     nvj=i_vac-1
   Else
!    No any additional sorting
     nvj=np*ns
   Endif
   Close(5)

 End Subroutine Bv_sort

! -----------------------------------------------------------
!     Block of Auxilary procedures
! -----------------------------------------------------------

! -----------------------------------------------------------
!    Standart spline procedure
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
            end do
        CALL TRIDIA(ALPHA,BETA,GAMMA,B,CM(2:n),N-2)
      END
      SUBROUTINE TRIDIA(ALPHA,BETA,GAMMA,B,X,N)
        Integer i,n,j
        Real*8 ALPHA(1:n+2),BETA(1:n+2),GAMMA(1:n+2),B(1:n+2),X(1:n+1),RAP
        DO I=2,N
        RAP=BETA(I-1)/ALPHA(I-1)
        ALPHA(I)=ALPHA(I)-RAP*GAMMA(I-1)
        B(I)=B(I)-RAP*B(I-1)
        end do
        X(N)=B(N)/ALPHA(N)
        DO J=2,N
        I=N-J+1
        X(I)=(B(I)-GAMMA(I)*X(I+1))/ALPHA(I)
        END DO
      END
      Real*8 FUNCTION SPL(N,X,Y,M,T)
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
! -----------------------------------------------------------
! The procedure searches for the string "Name_of_in" in the input file
 Subroutine FindStr(Name_of_in,num_file,str_to_find)
   Integer num_file,ifile
   Character*7 blockname,str_to_find
   Character*10 Name_of_in
   Open(num_file,File=Name_of_in)
   do ifile = 1,1000
     read(num_file,'(A7)',end=5)blockname
     If (blockname.eq.str_to_find)Return
   Enddo
5  Close(num_file)
   Print *,'The string ',str_to_find,' is not found in the file ',Name_of_in
 End Subroutine FindStr
! -----------------------------------------------------------
 Character*3 function cex(nj)
   integer cnj,dnj,unj,nj
   cnj = nj/100
   dnj = (nj - cnj*100)/10
   unj = nj - cnj*100 - dnj*10
   cex = char(cnj+48)//char(dnj+48)//char(unj+48)
 end function cex
! -----------------------------------------------------------
