program nca_hm_bethe
  USE DMFT_NCA
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                                     :: dmft_loop,Nb,Le
  logical                                     :: converged
  real(8)                                     :: wband

  !The local hybridization function:
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  complex(8),allocatable,dimension(:,:,:,:,:) :: Delta,Smats,Sreal,Gmats,Greal
  character(len=16)                           :: finput
  real(8)                                     :: wmixing,Eout(2),de,dens
  real(8),allocatable                         :: Gtau(:)
  real(8),dimension(:,:,:),allocatable        :: He
  real(8),dimension(:),allocatable            :: Wte
  logical                                     :: sc_bethe

  call parse_cmd_variable(finput,"FINPUT",default='inputNCA.conf')
  call parse_input_variable(wband,"wband",finput,default=1d0)
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(sc_bethe,"SC_BETHE",finput,default=.false.)
  !
  call nca_read_input(trim(finput))

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")


  allocate(He(Nspin*Norb,Nspin*Norb,Le),Wte(Le))
  He(1,1,:)         = linspace(-Wband,Wband,Le,mesh=de)
  He(Nspin,Nspin,:) = linspace(-Wband,Wband,Le,mesh=de)
  Wte               = dens_bethe(He(1,1,:),wband)*de


  !Allocate Fields:
  allocate(Delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc=zero

  Smats=zero
  call dmft_gloc_matsubara(one*He,Wte,Gmats,Smats,iprint=1)
  Delta = wband**2*0.25d0*Gmats

  call nca_init_solver(Hloc)

  converged=.false. ; dmft_loop=0 
  do while(.not.converged.AND.dmft_loop<nloop)
     dmft_loop=dmft_loop+1
     call start_loop(dmft_loop,nloop,"DMFT-loop")

     call nca_solve(Delta)
     Smats = impSmats !??

     ! compute the local gf:
     call dmft_gloc_matsubara(one*He,Wte,Gmats,Smats,iprint=1)

     if(sc_bethe)then
        Delta = Wband**2*0.25d0*impGmats
     else
        call dmft_delta(Gmats,Smats,Delta,Hloc,iprint=1)
     endif


     converged = check_convergence(Delta(1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)

     call end_loop
  enddo

  call dmft_kinetic_energy(one*He,Wte,Smats)

  ! open(100,file="Delta_iw.nca")
  ! do i=1,Lmats
  !    write(100,*)pi/beta*(2*i-1),dimag(Delta(1,1,1,1,i)),dreal(Delta(1,1,1,1,i))
  ! enddo
  ! close(100)
end program nca_hm_bethe
