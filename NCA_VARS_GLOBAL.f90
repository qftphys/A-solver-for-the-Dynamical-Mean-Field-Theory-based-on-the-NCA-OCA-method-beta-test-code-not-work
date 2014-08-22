MODULE NCA_VARS_GLOBAL  
  implicit none

  !SIZE OF THE PROBLEM
  !Ns   = # of levels per spin
  !Ntot = 2*Ns = total #  of levels
  !NN   = 2**Ntot = 2**(2Ns) max size of the Hilbert space
  !Nbo  =# number of bath sites (all sites - impurity sites)
  !Nsect=# of sectors
  !=========================================================
  integer                                                 :: Ns,Ntot,NN
  integer                                                 :: Nsect


  !local part of the Hamiltonian
  !=========================================================
  complex(8),dimension(:,:,:,:),allocatable               :: Hloc           !local hamiltonian


  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  integer,allocatable,dimension(:,:)                      :: getsector
  integer,allocatable,dimension(:,:)                      :: getCsector
  integer,allocatable,dimension(:,:)                      :: getCDGsector
  integer,allocatable,dimension(:,:)                      :: impIndex
  integer,allocatable,dimension(:)                        :: getdim,getnup,getndw


  !derived types:
  !=========================================================  
  type complete_espace
     real(8),dimension(:),allocatable                     :: E
     real(8),dimension(:,:),allocatable                   :: H
  end type complete_espace
  type hilbert_space_operator
     real(8),dimension(:,:),allocatable                   :: Op
  end type hilbert_space_operator
  type(complete_espace),dimension(:),allocatable          :: Espace
  type(hilbert_space_operator),dimension(:,:),allocatable :: Coperator,CDGoperator


  !Partition functions
  !=========================================================
  real(8)                                                 :: zeta_function_diag
  real(8)                                                 :: zeta_function


  !Anderson model Hybridization function Delta_nca(tau) : (Nspin,Nspin,Norb,Norb,0:Ltau)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:)             :: NcaDeltaAnd_iw
  real(8),allocatable,dimension(:,:,:,:,:)                :: NcaDeltaAnd_tau


  !Impurity Green's function and Self-Energies: (Nspin,Nspin,Norb,Norb,:)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:)             :: impSmats
  complex(8),allocatable,dimension(:,:,:,:,:)             :: impGmats


END MODULE NCA_VARS_GLOBAL
