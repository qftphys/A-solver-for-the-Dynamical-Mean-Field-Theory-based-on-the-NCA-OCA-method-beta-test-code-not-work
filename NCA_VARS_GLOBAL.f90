MODULE NCA_VARS_GLOBAL
  USE NCA_INPUT_VARS
  implicit none


  !SIZE OF THE PROBLEM
  !Ns       = # of levels per spin
  !Nlevels  = 2*Ns = total #  of levels
  !Nfock = 2**Nlevels = 2**(2Ns) max size of the Hilbert space
  !Nsectors = # of sectors
  !=========================================================
  integer                                                 :: Ns
  integer                                                 :: Nlevels
  integer                                                 :: Nsectors
  integer                                                 :: Nfock


  !local part of the Hamiltonian
  !=========================================================
  complex(8),dimension(:,:,:,:),allocatable               :: impHloc           !local hamiltonian


  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  integer,allocatable,dimension(:,:)                      :: getsector
  integer,allocatable,dimension(:,:)                      :: getCsector
  integer,allocatable,dimension(:,:)                      :: getCDGsector
  integer,allocatable,dimension(:)                        :: getdim,getnup,getndw



  !TO BE REMOVED or MODIFIED
  type hilbert_space_operator
     real(8),dimension(:,:),allocatable :: Op
  end type hilbert_space_operator
  type(hilbert_space_operator),dimension(:,:),allocatable :: Coperator,CDGoperator
  real(8),dimension(:,:),allocatable                      :: EigBasis
  real(8),dimension(:),allocatable                        :: EigValues


  !SECTOR-TO-FOCK SPACE STRUCTURE
  !=========================================================
  type sector_map
     integer,dimension(:),allocatable :: map
  end type sector_map
  interface map_allocate
     module procedure                 :: map_allocate_scalar
     module procedure                 :: map_allocate_vector
  end interface map_allocate
  interface map_deallocate
     module procedure                 :: map_deallocate_scalar
     module procedure                 :: map_deallocate_vector
  end interface map_deallocate


  !Partition functions
  !=========================================================
  real(8)                                                 :: zeta_function



  !Anderson model Hybridization function Delta_nca(tau) : (Nspin,Nspin,Norb,Norb,0:Ltau)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:)             :: NcaDeltaAnd_iw
  real(8),allocatable,dimension(:,:,:,:,:)                :: NcaDeltaAnd_tau



  !Atomic limit bare and dressed propagato: ncaR0(tau),ncaR(tau): [Nfock][Nfock][0:Ltau]
  !=========================================================
  real(8),allocatable,dimension(:,:,:)                    :: ncaR0
  real(8),allocatable,dimension(:,:,:)                    :: ncaR


  !Impurity Green's function and Self-Energies: (Nspin,Nspin,Norb,Norb,:)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:)             :: impSmats
  complex(8),allocatable,dimension(:,:,:,:,:)             :: impStail
  complex(8),allocatable,dimension(:,:,:,:,:)             :: impGmats
  complex(8),allocatable,dimension(:,:,:,:,:)             :: impG0mats

  !Density and double occupancy
  !=========================================================
  real(8),dimension(:),allocatable                        ::  nca_dens
  real(8),dimension(:),allocatable                        ::  nca_dens_up
  real(8),dimension(:),allocatable                        ::  nca_dens_dw
  real(8),dimension(:),allocatable                        ::  nca_docc
  real(8),dimension(:),allocatable                        ::  nca_sz2
  real(8),dimension(:,:),allocatable                      ::  nca_zeta



contains


  subroutine map_allocate_scalar(H,N)
    type(sector_map)              :: H
    integer                       :: N
    allocate(H%map(N))
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:) :: H
    integer,dimension(size(H))    :: N
    integer                       :: i
    do i=1,size(H)
       allocate(H(i)%map(N(i)))
    enddo
  end subroutine map_allocate_vector


  subroutine map_deallocate_scalar(H)
    type(sector_map)              :: H
    deallocate(H%map)
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:) :: H
    integer                       :: i
    do i=1,size(H)
       deallocate(H(i)%map)
    enddo
  end subroutine map_deallocate_vector



END MODULE NCA_VARS_GLOBAL
