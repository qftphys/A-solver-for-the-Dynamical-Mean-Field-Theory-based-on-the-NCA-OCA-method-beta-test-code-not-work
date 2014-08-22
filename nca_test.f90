program nca_test
  USE DMFT_NCA
  USE SCIFOR
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
  ! epsk=[-2.571482655102 ,-1.368374574376 ,-0.388770595499 ,0.003495419243 ,1.368693417935 ,0.382102611201 ,2.580412764301]
  ! vk=[0.337849045837 ,-0.101748602044 ,0.014286205613 ,-0.004336228118 ,-0.104680410833 ,0.011945938170 ,0.337108241410]
  ! forall(iorb=1:Norb,i=1:Lmats)Delta(1,1,iorb,iorb,i)=sum(vk**2/(xi*pi/beta*dble(2*i-1)-epsk))
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

end program nca_test
