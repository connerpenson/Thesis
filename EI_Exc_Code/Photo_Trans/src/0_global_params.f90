module global_params

    integer, parameter  ::  e_max = 2, v_max = 3

    character(len=11)   ::  bsp_wf_dir = "./B-SPLINE/"
    character(len=7)    ::  dips_dir   = "./Dipo/"
    character(len=15)   ::  ens_dir    = "./Vib_Energies/" 

    real(8), parameter  ::  au_to_cm   = 219474.63
    real(8), parameter  ::  const      = 2.026e-6

end module global_params