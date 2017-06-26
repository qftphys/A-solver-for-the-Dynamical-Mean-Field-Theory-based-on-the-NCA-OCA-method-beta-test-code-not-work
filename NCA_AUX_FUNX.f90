MODULE NCA_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg
  USE NCA_INPUT_VARS
  USE NCA_VARS_GLOBAL
  implicit none
  !private

  interface print_state_vector
     module procedure print_state_vector_ivec,print_state_vector_int
  end interface print_state_vector

  public :: print_Hloc
  !
  public :: nca_init_structure
  !
  public :: setup_pointers
  public :: build_sector
  public :: bdecomp
  public :: bjoin
  public :: print_state_vector
  !
  public :: c
  public :: cdg
  public :: build_c_operator
  public :: build_cdg_operator
  public :: binary_search
  !
  public :: setup_eigenspace
  public :: reset_eigenspace

contains





  !+------------------------------------------------------------------+
  !PURPOSE  : Init calculation
  !+------------------------------------------------------------------+
  subroutine nca_init_structure(Hunit)
    character(len=*)                         :: Hunit
    logical                                  :: control
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: reHloc         !local hamiltonian, real part 
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: imHloc         !local hamiltonian, imag part
    integer                                  :: NP,nup,ndw,iorb,jorb,ispin,jspin
    !
    !Norb=# of impurity orbitals
    !Ns=total number of sites
    Ns = Norb
    Ntot  = 2*Ns
    NN    = 2**Ntot
    !
    nup=Ns/2
    ndw=Ns-nup
    Nsect = (Ns+1)*(Ns+1)
    NP=get_sector_dimension(nup,ndw)
    !
    write(*,*)"Summary:"
    write(*,*)"--------------------------------------------"
    write(*,*)'Number of sites              = ',Norb
    write(*,*)'Max. Number of electrons     = ',Ntot
    write(*,*)'Hilbert space dimension      = ',NN
    write(*,*)'Number of sectors            = ',Nsect
    write(*,*)'Largest sector dimension     = ',NP
    write(*,*)"--------------------------------------------"

    allocate(Hloc(Nspin,Nspin,Norb,Norb))
    reHloc = 0.d0
    imHloc = 0.d0

    inquire(file=Hunit,exist=control)
    if(control)then
       write(*,*)"Reading Hloc from file: "//Hunit
       open(50,file=Hunit,status='old')
       do ispin=1,Nspin
          do iorb=1,Norb
             read(50,*)((reHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       do ispin=1,Nspin
          do iorb=1,Norb
             read(50,*)((imHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       close(50)
    else
       write(*,*)"Hloc file not found."
       write(*,*)"Hloc should be defined elsewhere..."
    endif
    Hloc = dcmplx(reHloc,imHloc)
    write(*,"(A)")"H_local:"
    call print_Hloc(Hloc)

    !Allocate pointers
    allocate(impIndex(2,Norb))
    allocate(getdim(Nsect),getnup(Nsect),getndw(Nsect))
    allocate(getsector(0:Ns,0:Ns))

    !Some check:
    if(Nspin>2)stop "Nspin > 2 ERROR. ask developer or develop your own on separate branch..."
    if(Norb>3)stop "Norb > 3 ERROR. ask developer or develop your own on separate branch..." 

    !allocate nca_delta
    allocate(NcaDeltaAnd_tau(Nspin,Nspin,Norb,Norb,0:Ltau))
    allocate(NcaDeltaAnd_iw(Nspin,Nspin,Norb,Norb,Lmats))

    !allocate functions
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGmats(Nspin,Nspin,Norb,Norb,Lmats))
 
    !allocate memory for operators
    allocate(Coperator(2,Norb))
    allocate(CDGoperator(2,Norb))

  end subroutine nca_init_structure









  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_Hloc(hloc,unit)
    integer,optional                            :: unit
    integer                                     :: iorb,jorb,ispin,jspin
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          write(*,"(20(A1,F7.3,A1,F7.3,A1,2x))")&
               (&
               (&
               '(',dreal(Hloc(ispin,jspin,iorb,jorb)),',',dimag(Hloc(ispin,jspin,iorb,jorb)),')',&
               jorb =1,Norb),&
               jspin=1,Nspin)
       enddo
    enddo
    if(present(unit))then
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(90F12.6)")((dreal(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(90F12.6)")((dimag(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
    endif
  end subroutine print_Hloc






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_pointers
    integer                          :: ispin,dim,isector
    integer                          :: nup,ndw,iorb
    write(*,"(A)")"Setting up pointers:"
    call start_timer
    isector=0
    do nup=0,Ns
       do ndw=0,Ns
          isector=isector+1
          getsector(nup,ndw)=isector
          getnup(isector)=nup
          getndw(isector)=ndw
          dim = get_sector_dimension(nup,ndw)
          getdim(isector)=dim
       enddo
    enddo
    call stop_timer

    do ispin=1,2
       do iorb=1,Norb
          impIndex(ispin,iorb)=iorb+(ispin-1)*Ns
       enddo
    enddo
  end subroutine setup_pointers





  !+------------------------------------------------------------------+
  !PURPOSE  : constructs the sectors by storing the map to the 
  !states i\in Hilbert_space from the states count in H_sector.
  !+------------------------------------------------------------------+
  !|ImpUP,BathUP>|ImpDW,BathDW >
  subroutine build_sector(isector,map)
    integer              :: i,isector,iup,idw,dim
    integer              :: nup,ndw
    integer              :: ivec(Ntot)
    integer,dimension(:) :: map
    dim=0
    nup = getnup(isector)
    ndw = getndw(isector)
    do i=1,NN
       call bdecomp(i,ivec)
       iup = sum(ivec(1:Ns))
       idw = sum(ivec(Ns+1:2*Ns))
       if(iup==nup.AND.idw==ndw)then
          dim           = dim+1 !count the states in the sector (n_up,n_dw)
          map(dim)      = i       !build the map to full space states
       endif
    enddo
  end subroutine build_sector







  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Ntot)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  subroutine bdecomp(i,ivec)
    integer :: ivec(Ntot)         
    integer :: l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i-1,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end subroutine bdecomp



  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Ntot) with the binary sequence 
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  subroutine bjoin(ivec,i)
    integer,dimension(ntot) :: ivec
    integer                 :: i,j
    i=1
    do j=1,Ntot
       i=i+ivec(j)*2**(j-1)
    enddo
  end subroutine bjoin




  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(m,i,j,sgn)
    integer :: ib(Ntot)
    integer :: i,j,m,km
    integer :: isg
    real(8) :: sgn
    call bdecomp(i,ib)
    if (ib(m)==0)then
       j=0
    else
       if(m==1)then
          j=i-1
       else
          km=sum(ib(1:m-1))
          isg=(-1)**km
          j=(i-2**(m-1))*isg
       endif
    endif
    sgn=dfloat(j)/dfloat(abs(j))
    j=abs(j)
  end subroutine c

  subroutine build_c_operator(ispin,iorb,Cmat)
    integer                  :: ispin,iorb
    real(8),dimension(NN,NN) :: Cmat
    integer                  :: imp,i,j
    integer                  :: ib(Ntot)
    real(8)                  :: c_
    !build <j|C|i>
    imp = impIndex(ispin,iorb)
    Cmat=0d0
    do i=1,NN
       call bdecomp(i,ib)
       if(ib(imp)==0)cycle
       call c(imp,i,j,c_)
       Cmat(i,j)=c_
    enddo
    return
  end subroutine build_c_operator



  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm+|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine cdg(m,i,j,sgn)
    integer :: ib(Ntot)
    integer :: i,j,m,km
    integer :: isg
    real(8) :: sgn
    call bdecomp(i,ib)
    if (ib(m)==1)then
       j=0
    else
       if(m==1)then
          j=i+1
       else
          km=sum(ib(1:m-1))
          isg=(-1)**km
          j=(i+2**(m-1))*isg
       endif
    endif
    sgn=dfloat(j)/dfloat(abs(j))
    j=abs(j)
  end subroutine cdg

  subroutine build_cdg_operator(ispin,iorb,Cmat)
    integer                  :: ispin,iorb
    real(8),dimension(NN,NN) :: Cmat
    integer                  :: imp,i,j
    integer                  :: ib(Ntot)
    real(8)                  :: c_
    !build <j|C^+|i>
    imp = impIndex(ispin,iorb)
    Cmat=0d0
    do i=1,NN
       call bdecomp(i,ib)
       if(ib(imp)==1)cycle
       call cdg(imp,i,j,c_)
       Cmat(i,j)=c_
    enddo
    return
  end subroutine build_cdg_operator




  !+------------------------------------------------------------------+
  !PURPOSE  : print a state vector |{up}>|{dw}>
  !+------------------------------------------------------------------+
  subroutine print_state_vector_ivec(ib,unit)
    integer,optional :: unit
    integer :: i,j,unit_
    integer :: ib(ntot)
    unit_=6;if(present(unit))unit_=unit
    call bjoin(ib,i)
    write(unit_,"(I3,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ib(j),j=1,Ns)
    write(unit_,"(A1,A1)",advance="no")">","|"
    write(unit_,"(10I1)",advance="no")(ib(ns+j),j=1,Ns)
    write(unit_,"(A1)",advance="yes")">"
  end subroutine print_state_vector_ivec
  !
  subroutine print_state_vector_int(i,unit)
    integer,optional :: unit
    integer :: i,j,unit_
    integer :: ib(ntot)
    unit_=6;if(present(unit))unit_=unit
    call bdecomp(i,ib)
    write(unit_,"(I3,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ib(j),j=1,Ns)
    write(unit_,"(A2)",advance="no")">|"
    write(unit_,"(10I1)",advance="no")(ib(ns+j),j=1,Ns)
    write(unit_,"(A1)",advance="yes")">"
  end subroutine print_state_vector_int



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial
  !+------------------------------------------------------------------+
  function get_sector_dimension(nup,ndw) result(dim)
    integer :: nup,ndw,dim,dimup,dimdw
    dimup=binomial(Ns,nup)
    dimdw=binomial(Ns,ndw)
    dim=dimup*dimdw
  end function get_sector_dimension



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial of an integer N!=1.2.3...(N-1).N
  !+------------------------------------------------------------------+
  recursive function factorial(n) result(f)
    integer            :: f
    integer,intent(in) :: n
    if(n<=0)then
       f=1
    else
       f=n*factorial(n-1)
    end if
  end function factorial



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor
  !+------------------------------------------------------------------+
  function binomial(n1,n2) result(nchoos)
    real(8) :: xh
    integer :: n1,n2,i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*real(n1+1-i,8)/real(i,8)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial




  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a value in an array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_eigenspace
    integer :: isector,dim
    if(allocated(espace)) deallocate(espace)
    allocate(espace(1:Nsect))
    do isector=1,Nsect
       dim=getdim(isector)
       allocate(espace(isector)%e(dim),espace(isector)%H(dim,dim))
    enddo
  end subroutine setup_eigenspace



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine reset_eigenspace
    integer :: isector
    forall(isector=1:Nsect)
       espace(isector)%E=0.d0
       espace(isector)%H=0.d0
    end forall
  end subroutine reset_eigenspace


END MODULE NCA_AUX_FUNX
