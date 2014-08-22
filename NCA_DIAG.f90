!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!|{ImpUP1,...,ImpUPN},BathUP>|{ImpDW1,...,ImpDWN},BathDW>
!########################################################################
module NCA_DIAG
  USE CONSTANTS,only: one,xi,zero,pi
  USE TIMER
  USE IOTOOLS,  only: free_unit,reg,free_units,txtfy
  USE ARRAYS,   only: arange,linspace
  USE MATRIX,   only: matrix_inverse,matrix_diagonalize
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
    real(8),dimension(NN,NN) :: Htotal
    e0=0.d0
    write(LOGfile,"(A)")"Get Hamiltonian:"
    call start_progress(LOGfile)
    dim_stride=0
    Htotal=0.d0
    do isector=1,Nsect
       !call progress(isector,Nsect)
       dim=getdim(isector)
       call build_hlocal(isector,espace(isector)%H(:,:))
       !embed H into Htotal
       Htotal(dim_stride+1:dim_stride+dim,dim_stride+1:dim_stride+dim) = espace(isector)%H(1:dim,1:dim)
       dim_stride=dim_stride+dim
       call matrix_diagonalize(espace(isector)%H,espace(isector)%e,'V','U')
       do i=1,dim
          write(LOGfile,"(100(F8.4,1x))")(espace(isector)%H(i,j),j=1,dim)
       enddo
       write(LOGfile,*)""
       e0(isector)=minval(espace(isector)%e)
    enddo
    call stop_progress
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
    write(LOGfile,"(A)")"DIAG resume:"
    write(LOGfile,"(A,f20.12)")'egs  =',egs
    write(LOGfile,"(A,f20.12)")'Z    =',zeta_function_diag
    unit=free_unit()
    open(unit,file="Egs.nca")
    write(unit,*)egs
    close(unit)
    write(LOGfile,"(A)")"Htotal:"
    do i=1,NN
       write(LOGfile,"(100(F8.4,1x))")(Htotal(i,j),j=1,NN)
    enddo
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
