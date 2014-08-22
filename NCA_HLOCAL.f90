MODULE NCA_HLOCAL
  USE CONSTANTS,only:zero
  USE NCA_INPUT_VARS
  USE NCA_VARS_GLOBAL
  USE NCA_AUX_FUNX
  implicit none
  private

  public :: build_hlocal
  public :: build_c_operator
  public :: build_cdg_operator

contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix DOUBLE PRECISION
  !+------------------------------------------------------------------+
  subroutine build_hlocal(isector,Hmat)
    real(8),dimension(:,:)             :: Hmat
    integer                            :: isector
    integer,dimension(:),allocatable   :: Hmap    !map of the Sector S to Hilbert space H
    integer,dimension(Ntot)            :: ib
    integer                            :: dim
    integer                            :: i,j,m,iorb,jorb,ispin,impi,ishift
    integer                            :: k1,k2,k3,k4
    real(8)                            :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)            :: nup,ndw
    real(8)                            :: htmp
    real(8),dimension(Nspin,Norb)      :: eloc
    logical                            :: Jcondition
    integer                            :: first_state,last_state
    !
    dim=getdim(isector)
    allocate(Hmap(dim))
    call build_sector(isector,Hmap)
    !
    if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "NCA_HLOCAL/build_hlocal: wrong dimensions in H"
    Hmat = 0.d0
    !
    ishift     = 0
    first_state= 1
    last_state = dim
    !
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=dreal(Hloc(ispin,ispin,iorb,iorb))
       enddo
    enddo
    !
    !-----------------------------------------------!
    !BUILD ED HAMILTONIAN AS A SPARSE MATRIX
    !this part is identical between d_ and c_ codes.
    include "nca_build_hlocal.f90"
    !-----------------------------------------------!
    !
    deallocate(Hmap)
    !
  end subroutine build_hlocal


end MODULE NCA_HLOCAL
