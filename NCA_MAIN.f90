module NCA_MAIN
  USE FFTGF
  !
  USE NCA_INPUT_VARS
  USE NCA_VARS_GLOBAL
  USE NCA_AUX_FUNX
  USE NCA_DIAG
  USE NCA_GREENS_FUNCTIONS

  implicit none

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine nca_init_solver(Hunit)
    character(len=*),optional,intent(in)              :: Hunit
    character(len=100)                                :: Hunit_
    write(LOGfile,"(A)")"INIT NCA SOLVER"
    Hunit_='inputHLOC.in';if(present(Hunit))Hunit_=Hunit
    call nca_init_structure(Hunit_)
    call setup_pointers
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
    write(LOGfile,"(A)")"START NCA SOLVER"
    NcaDeltaAnd_iw = Delta
    do ispin=1,Nspin
       do iorb=1,Norb
          call fftgf_iw2tau(NcaDeltaAnd_iw(ispin,ispin,iorb,iorb,:),&
               NcaDeltaAnd_tau(ispin,ispin,iorb,iorb,0:),beta)
       enddo
    enddo
    call nca_build_dressed_propagator
    call nca_build_impurity_gf
  end subroutine nca_solver

end module NCA_MAIN
