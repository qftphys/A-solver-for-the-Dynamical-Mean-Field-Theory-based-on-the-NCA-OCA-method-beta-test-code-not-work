program nca_hm_square
  USE DMFT_NCA
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                                     :: iloop,Lk,Nso
  logical                                     :: converged
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:) :: Delta
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)     :: Hk
  complex(8),allocatable,dimension(:,:)       :: modelHloc
  complex(8),allocatable,dimension(:,:,:,:)   :: Hloc
  real(8),allocatable,dimension(:)            :: Wtk

  !variables for the model:
  integer                                     :: Nk,Nkpath
  real(8)                                     :: ts,wmixing
  character(len=32)                           :: finput
  character(len=32)                           :: hkfile
  logical :: bool
  !


  !Parse additional variables && read Input && read H(k)^2x2
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS","inputED.conf",default=0.25d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call nca_read_input(trim(finput))

  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  !
  if(Norb/=1)stop "Wrong setup from input file: Norb=1"
  Nso=Nspin*Norb



  !Allocate Weiss Field:
  allocate(Delta(Nspin,Nspin,Norb,Norb,Lmats));Delta=zero
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Hloc(Nspin,Nspin,Norb,Norb));Hloc=zero


  !Build the Hamiltonian on a grid or on a path
  call build_hk(trim(hkfile))
  Hloc = so2nn_reshape(modelHloc,Nspin,Norb)

  call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats,iprint=1)
  call dmft_delta(Gmats,Smats,Delta,Hloc,iprint=1)

  !Setup solver
  call nca_init_solver(Hloc)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     call nca_solve(Delta)

     Smats = impSmats

     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats,iprint=1)

     call dmft_delta(Gmats,Smats,Delta,Hloc,iprint=1)

     converged = check_convergence(Delta(1,1,1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo


  ! call dmft_kinetic_energy(Hk,Wtk,Smats)




contains




  !--------------------------------------------------------------------!
  !Lattice Hamitonian:
  !--------------------------------------------------------------------!
  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N,i
    complex(8),dimension(N,N) :: hk
    real(8)                   :: kx,ky
    !
    if(N/=Nso)stop "hk_model: N != Nso"
    !
    kx=kpoint(1)
    ky=kpoint(2)
    !
    hk(:,:) = zero
    !
    hk = -2d0*ts*(cos(kx)+cos(ky))*eye(N)
    !
  end function hk_model






  !---------------------------------------------------------------------
  !PURPOSE: get model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional             :: file
    integer                               :: i,j,ik
    integer                               :: ix,iy
    real(8)                               :: kx,ky  
    integer                               :: iorb,jorb
    integer                               :: isporb,jsporb
    integer                               :: ispin,jspin
    integer                               :: unit
    real(8),dimension(2)                  :: kvec
    real(8)                               :: blen,area_hex,area_rect,points_in,points_tot
    real(8),allocatable,dimension(:)      :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable    :: kpath

    Lk= Nk*Nk

    write(LOGfile,*)"Build H(k)    :",Lk
    write(LOGfile,*)"# of SO-bands :",Nso

    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk));Hk=zero
    allocate(wtk(Lk));Wtk=0d0

    call TB_set_bk([pi2,0d0],[0d0,pi2])

    call TB_build_model(Hk, hk_model, Nso, [Nk,Nk])

    Wtk = 1d0/Lk

    if(present(file))&
         call TB_write_hk(Hk, trim(file), &
         No = Nso,Nd = Norb,Np = 0,Nineq = 1,&
         Nkvec=[Nk,Nk])
    !
    allocate(modelHloc(Nso,Nso))
    modelHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs((modelHloc))<1.d-4)modelHloc=zero

    !path: G X M G
    allocate(kpath(4,2))
    kpath(1,:)=[0d0,0d0]
    kpath(2,:)=[ pi,0d0]
    kpath(3,:)=[ pi, pi]
    kpath(4,:)=[0d0,0d0]
    call TB_solve_model(hk_model,Nso,kpath,Nkpath,&
         colors_name=[red1],&
         points_name=[character(len=10) :: "G","X","M", "G"],&
         file="Eigenbands.nint")
  end subroutine build_hk








end program nca_hm_square



