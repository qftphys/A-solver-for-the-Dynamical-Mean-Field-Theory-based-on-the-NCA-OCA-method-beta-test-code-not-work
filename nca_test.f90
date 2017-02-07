program nca_test
  USE DMFT_NCA
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter                           :: Nbath=7
  real(8),dimension(Nbath)                    :: epsk,vk
  complex(8),dimension(:,:,:,:,:),allocatable :: Delta
  integer                                     :: i,iorb
  !
  logical :: converged
  integer :: dmft_loop

  print*, "NCA_TEST:"
  call nca_read_input("inputNCA.conf")
  allocate(Delta(Nspin,Nspin,Norb,Norb,Lmats))
  Delta=zero
  open(100,file="Delta_iw.input")
  do i=1,Lmats
     Delta(1,1,1,1,i) = 0.5d0/(xi*pi/beta*(2*i-1) - Uloc(1)/2d0)+0.5d0/(xi*pi/beta*(2*i-1) + Uloc(1)/2d0)
     write(100,*)pi/beta*(2*i-1),dimag(Delta(1,1,1,1,i)),dreal(Delta(1,1,1,1,i))
  enddo

  call nca_init_solver()

  converged=.false. ; dmft_loop=0 
  do while(.not.converged)
     dmft_loop=dmft_loop+1
     call start_loop(dmft_loop,nloop,"DMFT-loop")
     call nca_solver(Delta)
     Delta = 0.25d0*impGmats
     converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     call end_loop
  enddo
  open(100,file="Delta_iw.nca")
  do i=1,Lmats
     write(100,*)pi/beta*(2*i-1),dimag(Delta(1,1,1,1,i)),dreal(Delta(1,1,1,1,i))
  enddo
  close(100)
end program nca_test
