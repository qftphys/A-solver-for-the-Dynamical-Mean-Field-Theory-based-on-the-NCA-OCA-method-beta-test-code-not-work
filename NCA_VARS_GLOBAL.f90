MODULE NCA_VARS_GLOBAL
  USE SF_PARSE_INPUT
  USE SF_VERSION
  implicit none


  ! !GIT VERSION
  include "revision.inc"  !this file is generated at compilation time in the Makefile

  !SIZE OF THE PROBLEM
  !Ns       = # of levels per spin
  !Nlevels  = 2*Ns = total #  of levels
  !Nhilbert = 2**Nlevels = 2**(2Ns) max size of the Hilbert space
  !Nsectors = # of sectors
  !=========================================================
  integer                                                 :: Ns
  integer                                                 :: Nlevels
  integer                                                 :: Nsectors
  integer                                                 :: Nhilbert



  !local part of the Hamiltonian
  !=========================================================
  complex(8),dimension(:,:,:,:),allocatable               :: impHloc           !local hamiltonian


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
     real(8),dimension(:),pointer   :: E
     real(8),dimension(:,:),pointer :: H
  end type complete_espace
  type(complete_espace),dimension(:),allocatable          :: Espace


  !TO BE REMOVED or MODIFIED
  type hilbert_space_operator
     real(8),dimension(:,:),allocatable :: Op
     real(8),dimension(:),allocatable   :: val
     real(8),dimension(:),allocatable   :: col
  end type hilbert_space_operator
  type(hilbert_space_operator),dimension(:,:),allocatable :: Coperator,CDGoperator



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


  !Density and double occupancy
  !=========================================================
  real(8),dimension(:),allocatable                        ::  nca_dens
  real(8),dimension(:),allocatable                        ::  nca_dens_up,nca_dens_dw
  real(8),dimension(:),allocatable                        ::  nca_docc



  !input variables
  !=========================================================
  integer                                                 :: Norb                !Norb =# of impurity orbitals
  integer                                                 :: Nspin               !Nspin=# spin degeneracy (max 2)
  integer                                                 :: nloop               !max dmft loop variables
  real(8)                                                 :: Ust                 !intra-orbitals interactions
  real(8)                                                 :: Jh                  !J_Hund: Hunds' coupling constant 
  real(8)                                                 :: Jx                  !J_X: coupling constant for the spin-eXchange interaction term
  real(8)                                                 :: Jp                  !J_P: coupling constant for the Pair-hopping interaction term 
  real(8),dimension(3)                                    :: Uloc                !local interactions
  real(8)                                                 :: xmu                 !chemical potential
  real(8)                                                 :: beta                !inverse temperature
  integer                                                 :: Nsuccess            !
  logical                                                 :: Jhflag              !spin-exchange and pair-hopping flag.
  logical                                                 :: HFmode              !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)                                                 :: cutoff              !cutoff for spectral summation
  real(8)                                                 :: gs_threshold        !Energy threshold for ground state degeneracy loop up
  real(8)                                                 :: dmft_error          !dmft convergence threshold
  real(8)                                                 :: nca_error        !Energy threshold for ground state degeneracy loop up


  !Some parameters for function dimension:
  !=========================================================
  integer                                                 :: Lmats
  integer                                                 :: Ltau


  !LOG AND Hamiltonian UNITS
  !=========================================================
  integer,save                                            :: LOGfile



contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine nca_read_input(INPUTunit)
    character(len=*)                                      :: INPUTunit
    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(beta,"BETA",INPUTunit,default=500.d0,comment="Inverse temperature, at T=0 is used as a IR cut-off.")
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of impurity orbitals.")
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy (max 2)")
    call parse_input_variable(uloc,"ULOC",INPUTunit,default=[5.d0,0.d0,0.d0],comment="Values of the local interaction per orbital (max 3)")
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0,comment="Value of the inter-orbital interaction term")
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0,comment="Hunds coupling")
    call parse_input_variable(Jx,"JX",INPUTunit,default=0.d0,comment="S-E coupling")
    call parse_input_variable(Jp,"JP",INPUTunit,default=0.d0,comment="P-H coupling")
    call parse_input_variable(jhflag,"JHFLAG",INPUTunit,default=.false.,comment="Flag to include full SU(2) invariant term: spin-flip, pair-hopping.")
    call parse_input_variable(hfmode,"HFMODE",INPUTunit,default=.true.,comment="Flag to set the Hartree form of the interaction (n-1/2). see xmu.")
    !
    call parse_input_variable(nloop,"NLOOP",INPUTunit,default=100,comment="Max number of DMFT iterations.")
    call parse_input_variable(dmft_error,"DMFT_ERROR",INPUTunit,default=0.00001d0,comment="Error threshold for DMFT convergence")
    call parse_input_variable(nca_error,"NCA_ERROR",INPUTunit,default=1d-12,comment="Error threshold for the NCA cycle")
    call parse_input_variable(nsuccess,"NSUCCESS",INPUTunit,default=1,comment="Number of successive iterations below threshold for convergence")
    !
    call parse_input_variable(Lmats,"LMATS",INPUTunit,default=2000,comment="Number of Matsubara frequencies.")
    call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1000,comment="Number of imaginary time points.")
    call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9,comment="Spectrum cut-off, used to determine the number states to be retained.")
    call parse_input_variable(gs_threshold,"GS_THRESHOLD",INPUTunit,default=1.d-9,comment="Energy threshold for ground state degeneracy loop up")
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6,comment="LOG unit.")
    Ltau=max(int(beta),Ltau)
    !
    call save_input_file(INPUTunit)
    call scifor_version()
    call code_version(revision)
    !
  end subroutine nca_read_input






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
