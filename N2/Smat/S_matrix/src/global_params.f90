Module global_params
    complex(kind=8), parameter :: ci = (0.d0, 1.d0)
    integer, parameter :: num_gamma_Q = 8,l_max = 4, n_max = 3, vib_num = 9, num_Q_ = 10, neut_vib_num = 3 
    integer, parameter :: energy_step_num = 10000, num_conv_points = 5000, dip_mult = 5
    real(8), parameter :: sigma = .0001
    real*8 :: ryd,au_eV,au_cm,au_to_sec
    real(8) :: pi,ev_au
    integer :: num_Q, num_v ,num_ev,num_scatE,neut_vib
    integer, allocatable, dimension(:,:) ::  qn_ev
    real(kind=8), allocatable, dimension(:) ::  geom,energies_ev,E_el_Rmat,energies_v
    Complex*16, allocatable, dimension(:,:,:,:) :: S_Q,S_Q_check
    real(8), parameter, dimension(num_Q_) ::  geom_ = (/1.9, 1.944, 1.989, 2.0333, 2.0778, 2.122, 2.1667, 2.211, 2.255, 2.3/) ! in atomic units
    integer, parameter :: dim_lm = (l_max+1)**2 * n_max
    real(8), allocatable, dimension(:,:) :: chan_ens
    integer, parameter  ::  num_tot = dim_lm*vib_num
    real(8), parameter  ::  ion_shift = -.00432 !eV
    
  
  End Module global_params