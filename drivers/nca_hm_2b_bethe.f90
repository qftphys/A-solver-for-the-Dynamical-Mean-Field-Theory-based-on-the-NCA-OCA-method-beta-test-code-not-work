program nca_2b_bethe
  USE DMFT_NCA
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                       :: iloop,Le,Nso,Nlso,Nlat
  logical                                       :: converged
  integer                                       :: ispin,ilat!,i,j


  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:) :: Delta,Delta_prev
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats

  !hamiltonian input:

  complex(8),allocatable,dimension(:,:,:,:)   :: Hloc
  real(8),dimension(:,:),allocatable          :: Dbands
  real(8),dimension(:,:),allocatable          :: Ebands
  real(8),dimension(:),allocatable            :: H0,de

  !variables for the model:
  real(8),dimension(2)                          :: Wband
  real(8)                                       :: wmixing,Mh
  character(len=32)                             :: finput
  character(len=32)                             :: hkfile
  !


  !Parse additional variables && read Input && read H(k)^2x2
  call parse_cmd_variable(finput,"FINPUT",default='inputNCA.conf')
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Wband,"WBAND",finput,default=[1d0,1d0])
  call parse_input_variable(Mh,"MH",finput,default=0d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call nca_read_input(trim(finput))

  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  !
  if(Nspin/=1.OR.Norb/=2)stop "Wrong setup from input file: Nspin=1; Norb=2"
  Nso=Nspin*Norb
  Nlso=Nso



  !Allocate Weiss Field:
  allocate(Delta(Nspin,Nspin,Norb,Norb,Lmats));Delta=zero
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Hloc(Nspin,Nspin,Norb,Norb));Hloc=zero


  !Build the Hamiltonian on a grid or on a path
  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  allocate(de(Nso))
  Ebands(1,:) = linspace(-Wband(1),Wband(1),Le,mesh=de(1))
  Ebands(2,:) = linspace(-Wband(2),Wband(2),Le,mesh=de(2))
  !
  Dbands(1,:) = dens_bethe(Ebands(1,:),Wband(1))*de(1)
  Dbands(2,:) = dens_bethe(Ebands(2,:),Wband(2))*de(2)
  !
  allocate(H0(Nso))
  H0=[-Mh/2,Mh/2]
  call TB_write_Hloc(one*diag(H0))
  Hloc(1,1,:,:)=diag(H0)

  ! compute the local gf:
  call dmft_gloc_matsubara(Ebands,Dbands,H0,Gmats,Smats,iprint=1)
  call dmft_delta(Gmats,Smats,Delta,Hloc,iprint=1)

  !Setup solver
  call nca_init_solver(Hloc)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call nca_solve(Delta)!,iprint=1)

     Smats = impSmats

     ! compute the local gf:
     call dmft_gloc_matsubara(Ebands,Dbands,H0,Gmats,Smats,iprint=1)


     Delta_prev = Delta
     call dmft_delta(Gmats,Smats,Delta,Hloc,iprint=1)
     Delta = wmixing*Delta + wmixing*Delta_prev

     converged = check_convergence(Delta(1,1,1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo

  call dmft_kinetic_energy(Ebands,Dbands,H0,Smats(1,1,:,:,:))



end program nca_2b_bethe



