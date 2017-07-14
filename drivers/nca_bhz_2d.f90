!                    MODEL Hamiltonian is:
!
! |     h^{2x2}(k)             &          hso^{2x2}(k)        |
! |     [hso^{2x2}]*(k)        &         [h^{2x2}]*(-k)       |
!
! h^{2x2}(k):=
!
! | m-(Cos{kx}+Cos{ky})         & \lambda*(Sin{kx}-i*Sin{ky}) |
! | \lambda*(Sin{kx}+i*Sin{ky}) & -m+(Cos{kx}+Cos{ky})        |
!
! hso^{2x2}(k):=
! | xi*rh*(sin(kx)-xi*sin(ky))  &         \delta              |
! |         -\delta             &             0               |
program nca_bhz
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE DMFT_NCA
  !
  implicit none
  integer                :: iloop,Lk,Nso
  logical                :: converged
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:)
  complex(8),allocatable :: Gmats(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:)
  real(8),allocatable    :: Wtk(:)
  real(8),allocatable    :: kxgrid(:),kygrid(:)
  integer,allocatable    :: ik2ix(:),ik2iy(:)
  !variables for the model:
  integer                :: Nk,Nkpath
  real(8)                :: mh,lambda,wmixing
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  logical                :: spinsym
  !

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  !
  call nca_read_input(trim(finput))
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")

  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))

  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))


  Smats=zero
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats,iprint=1)
  call dmft_delta(Gmats,Smats,Delta,Hloc=j2so(bhzHloc),iprint=0)


  call nca_init_solver(Hloc=j2so(bhzHloc))

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call nca_solve(Delta)
     Smats = impSmats


     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats,iprint=1)

     call dmft_delta(Gmats,Smats,Delta,Hloc=j2so(bhzHloc),iprint=1)

     converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo

  call dmft_kinetic_energy(Hk,Wtk,Smats)


contains



  !---------------------------------------------------------------------
  !PURPOSE: GET BHZ HAMILTONIAN (from the NonInteracting code)
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy
    real(8)                             :: kx,ky    
    integer                             :: iorb,jorb
    integer                             :: isporb,jsporb
    integer                             :: ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    complex(8),dimension(Nso,Nso,Lmats) :: Gmats
    real(8)                             :: wm(Lmats)

    call build_hk_GXMG()

    write(LOGfile,*)"Build H(k) for BHZ:"
    Lk=Nk**2
    write(*,*)"# of k-points     :",Lk
    write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
    allocate(wtk(Lk))

    call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
    call TB_build_model(Hk,hk_bhz,Nso,[Nk,Nk])
    wtk = 1d0/Lk
    if(present(file))then
       call TB_write_hk(Hk,trim(file),Nso,&
            Nd=Norb,&
            Np=1,   &
            Nineq=1,&
            Nkvec=[Nk,Nk])
    endif
    allocate(bhzHloc(Nso,Nso))
    bhzHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0d0
    call TB_write_Hloc(bhzHloc)

  end subroutine build_hk




  !---------------------------------------------------------------------
  !PURPOSE: GET THE BHZ HAMILTONIAN ALONG THE Gamma-X-M-Gamma path
  !---------------------------------------------------------------------
  subroutine build_hk_GXMG(kpath_)
    integer                            :: i,j
    integer                            :: Npts
    real(8),dimension(:,:),optional        :: kpath_
    real(8),dimension(:,:),allocatable :: kpath
    character(len=64)                      :: file
    !This routine build the H(k) along the GXMG path in BZ,
    !Hk(k) is constructed along this path.
    if(present(kpath_))then
       write(LOGfile,*)"Build H(k) BHZ along a given path:"
       Npts = size(kpath_,1)
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,size(kpath_,2)))
       kpath=kpath_
       file="Eig_path.nint"
    else
       write(LOGfile,*)"Build H(k) BHZ along the path GXMG:"
       Npts = 4
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_Gamma
       kpath(2,:)=kpoint_M1
       kpath(3,:)=kpoint_X1
       kpath(4,:)=kpoint_Gamma
       file="Eigenbands.nint"
    endif
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
    allocate(wtk(Lk))
    call TB_build_model(Hk,hk_bhz,Nso,kpath,Nkpath)
    wtk = 1d0/Lk
    call TB_Solve_model(hk_bhz,Nso,kpath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
         file=reg(file))
  end subroutine build_hk_GXMG


  !--------------------------------------------------------------------!
  !BHZ HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_bhz(kvec,N) result(hk)
    real(8),dimension(:)      :: kvec
    integer                   :: N
    complex(8),dimension(N,N) :: hk
    real(8)                   :: kx,ky
    if(N/=Nso)stop "hk_bhz error: N != Nspin*Norb == 4"
    kx=kvec(1)
    ky=kvec(2)
    Hk          = zero
    Hk(1:2,1:2) = hk_bhz2x2(kx,ky)
    Hk(3:4,3:4) = conjg(hk_bhz2x2(-kx,-ky))
  end function hk_bhz

  function hk_bhz2x2(kx,ky) result(hk)
    real(8)                   :: kx,ky,epsik
    complex(8),dimension(2,2) :: hk
    epsik   = cos(kx)+cos(ky)
    hk(1,1) = mh - epsik
    hk(2,2) =-mh + epsik
    hk(1,2) = lambda*(sin(kx)-xi*sin(ky))
    hk(2,1) = lambda*(sin(kx)+xi*sin(ky))
  end function hk_bhz2x2










  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  subroutine read_sigma(sigma)
    complex(8)        :: sigma(:,:,:,:,:)
    integer           :: iorb,ispin,i,L,unit
    real(8)           :: reS(Nspin),imS(Nspin),ww
    character(len=20) :: suffix
    if(size(sigma,1)/=Nspin)stop "read_sigma: error in dim 1. Nspin"
    if(size(sigma,3)/=Norb)stop "read_sigma: error in dim 3. Norb"
    L=size(sigma,5);print*,L
    if(L/=Lmats)stop "read_sigma: error in dim 5. Lmats"
    do iorb=1,Norb
       unit=free_unit()
       suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_iw.ed"
       write(*,*)"read from file=","impSigma"//reg(suffix)
       open(unit,file="impSigma"//reg(suffix),status='old')
       do i=1,L
          read(unit,"(F26.15,6(F26.15))")ww,(imS(ispin),reS(ispin),ispin=1,Nspin)
          forall(ispin=1:Nspin)sigma(ispin,ispin,iorb,iorb,i)=dcmplx(reS(ispin),imS(ispin))
       enddo
       close(unit)
    enddo
  end subroutine read_sigma





  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg,Nso) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    integer                                     :: Nso,i,j,iorb,jorb,ispin,jspin
    complex(8),dimension(Nso,Nso)               :: g

    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nso,Nso)               :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so



end program nca_bhz



