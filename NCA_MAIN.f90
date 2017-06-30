module NCA_MAIN
  USE NCA_VARS_GLOBAL
  USE NCA_SETUP
  USE NCA_DIAG
  USE NCA_GREENS_FUNCTIONS
  !
  USE DMFT_FFTGF
  !
  implicit none

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine nca_init_solver(Hloc)
    complex(8),dimension(Nspin,Nspin,Norb,Norb),intent(in) :: Hloc
    write(*,"(A)")"INIT NCA SOLVER"
    !>debug
    print*,"entering nca_init_structure"
    !<debug
    call nca_init_structure()
    !>debug
    print*,"calling set_Hloc"
    !<debug
    call set_Hloc(Hloc)
    !>debug
    print*,"entering setup_pointers"
    !<debug
    call setup_pointers
    !>debug
    print*,"entering setup_eigenspace"
    !<debug
    call setup_eigenspace
    call nca_diagonalize_hlocal
    call nca_build_operators
    call nca_build_bare_propagator
  end subroutine nca_init_solver





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine nca_solver(Delta)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: Delta
    integer                                           :: ispin,iorb
    write(*,"(A)")"START NCA SOLVER"
    NcaDeltaAnd_iw = Delta
    do ispin=1,Nspin
       do iorb=1,Norb
          call fft_gf_iw2tau(NcaDeltaAnd_iw(ispin,ispin,iorb,iorb,:),&
               NcaDeltaAnd_tau(ispin,ispin,iorb,iorb,0:),beta)
       enddo
    enddo
    call nca_build_dressed_propagator
    call nca_build_impurity_gf
  end subroutine nca_solver

end module NCA_MAIN
