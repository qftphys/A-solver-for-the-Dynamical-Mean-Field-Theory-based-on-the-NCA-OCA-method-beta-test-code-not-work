!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!|{ImpUP1,...,ImpUPN},BathUP>|{ImpDW1,...,ImpDWN},BathDW>
!########################################################################
module NCA_DIAG
  USE SF_CONSTANTS,only: one,xi,zero,pi
  USE SF_TIMER
  USE SF_IOTOOLS,  only: free_unit,reg,free_units,txtfy
  USE SF_ARRAYS,   only: arange,linspace
  USE SF_LINALG,   only: inv,eigh
  !
  USE NCA_INPUT_VARS
  USE NCA_VARS_GLOBAL
  USE NCA_AUX_FUNX
  USE NCA_HLOCAL
  implicit none
  private

  public :: nca_diagonalize_hlocal
  public :: nca_build_operators


contains





  !+-------------------------------------------------------------------+
  !                    FULL DIAGONALIZATION
  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+-------------------------------------------------------------------+
  subroutine nca_diagonalize_hlocal
    integer                  :: i,j,isector,dim,unit,dim_stride
    real(8),dimension(Nsect) :: e0 
    real(8)                  :: egs
    e0=0.d0
    write(*,"(A)")"Get Hamiltonian:"
    call start_timer
    dim_stride=0
    do isector=1,Nsect
       !call progress(isector,Nsect)
       dim=getdim(isector)
       call build_hlocal(isector,espace(isector)%H(:,:))
       call eigh(espace(isector)%H,espace(isector)%e,'V','U')
       e0(isector)=minval(espace(isector)%e)
    enddo
    call stop_timer
    !
    egs=minval(e0)
    forall(isector=1:Nsect)espace(isector)%e = espace(isector)%e - egs
    !Get the partition function Z and rescale energies
    zeta_function_diag=0d0
    do isector=1,Nsect
       dim=getdim(isector)
       do i=1,dim
          zeta_function_diag=zeta_function_diag+exp(-beta*espace(isector)%e(i))
       enddo
    enddo
    write(*,"(A)")"DIAG resume:"
    write(*,"(A,f20.12)")'egs  =',egs
    write(*,"(A,f20.12)")'Z    =',zeta_function_diag
    unit=free_unit()
    open(unit,file="Egs.nca")
    write(unit,*)egs
    close(unit)
    return
  end subroutine nca_diagonalize_hlocal




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine nca_build_operators
    integer :: ispin,iorb
    do ispin=1,2
       do iorb=1,Norb
          allocate(Coperator(ispin,iorb)%Op(NN,NN))
          allocate(CDGoperator(ispin,iorb)%Op(NN,NN))
          call build_c_operator(ispin,iorb,Coperator(ispin,iorb)%Op(:,:))
          CDGoperator(ispin,iorb)%Op(:,:)=transpose(Coperator(ispin,iorb)%Op(:,:))
       enddo
    enddo
  end subroutine nca_build_operators



end MODULE NCA_DIAG
