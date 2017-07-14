module NCA_MAIN
  USE NCA_VARS_GLOBAL
  USE NCA_AUX_FUNX
  USE NCA_SETUP
  USE NCA_DIAG
  USE NCA_GREENS_FUNCTIONS
  USE NCA_OBSERVABLES
  !
  USE SF_ARRAYS
  USE SF_IOTOOLS
  USE DMFT_FFTGF
  !
  implicit none

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine nca_init_solver(Hloc)
    complex(8),dimension(Nspin,Nspin,Norb,Norb),intent(in) :: Hloc
    write(LOGfile,"(A)")"INIT NCA SOLVER"
    call nca_init_structure()
    call set_Hloc(Hloc)
    call setup_pointers
    call nca_diagonalize_hlocal
    call nca_build_operators
    call nca_build_bare_propagator
  end subroutine nca_init_solver





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine nca_solve(Delta)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: Delta
    integer                                           :: ispin,iorb
    real(8) :: tau(0:Ltau)
    write(LOGfile,"(A)")"START NCA SOLVER"
    NcaDeltaAnd_iw = Delta
    do ispin=1,Nspin
       do iorb=1,Norb
          call fft_gf_iw2tau(NcaDeltaAnd_iw(ispin,ispin,iorb,iorb,:),&
               NcaDeltaAnd_tau(ispin,ispin,iorb,iorb,0:),beta)
       enddo
    enddo
    call nca_build_dressed_propagator !<-- get ncaR(tau)
    call nca_get_observables
    call nca_get_impurity_gf
  end subroutine nca_solve

end module NCA_MAIN
