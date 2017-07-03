MODULE NCA_SETUP
  USE NCA_VARS_GLOBAL
  !
  USE SF_CONSTANTS, only:zero
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,str
  USE SF_MISC, only: assert_shape
  !
  implicit none
  !private

  interface set_Hloc
     module procedure :: set_Hloc_nn
     module procedure :: set_Hloc_so
  end interface set_Hloc


  public :: nca_init_structure
  !
  public :: setup_pointers
  public :: build_sector
  public :: bdecomp
  public :: binary_search
  !
  public :: set_Hloc
  public :: print_Hloc
  !
  public :: c
  public :: cdg
  public :: build_c_operator
  public :: build_cdg_operator
  public :: nca_build_operators
  !
  public :: setup_eigenspace
  public :: reset_eigenspace

contains





  !+------------------------------------------------------------------+
  !PURPOSE  : Init calculation
  !+------------------------------------------------------------------+
  subroutine nca_init_structure()
    integer :: NP,nup,ndw,iorb,jorb,ispin,jspin
    !
    !Norb=# of impurity orbitals
    !Ns=total number of sites
    Ns       = Norb
    !
    Nlevels  = 2*Ns
    !
    Nsectors = (Ns+1)*(Ns+1) !nup=0:Ns;ndw=0:Ns
    !    
    Nhilbert = 2**Nlevels
    !
    nup=Ns/2
    ndw=Ns-nup
    NP=get_sector_dimension(nup,ndw)
    !
    write(LOGfile,"(A)")"Summary:"
    write(LOGfile,"(A)")"--------------------------------------------"
    write(LOGfile,"(A,I15)")'# of levels/spin      = ',Ns
    write(LOGfile,"(A,I15)")'Total system size     = ',Nlevels
    write(LOGfile,"(A,I15)")'Fock space dimension  = ',Nhilbert
    write(LOGfile,"(A,I15)")'Number of sectors     = ',Nsectors
    write(LOGfile,"(A,2I15)")'Largest Sector(s)    = ',NP
    write(LOGfile,"(A)")"--------------------------------------------"


    allocate(impHloc(Nspin,Nspin,Norb,Norb))
    impHloc=zero


    !Allocate pointers
    allocate(impIndex(2,Norb))
    allocate(getdim(Nsectors),getnup(Nsectors),getndw(Nsectors))
    allocate(getsector(0:Ns,0:Ns))


    !Some check:
    if(Nspin>=2)stop "Nspin > 2 ERROR. ask developer or develop your own on separate branch..."
    if(Norb>3)stop "Norb > 3 ERROR. ask developer or develop your own on separate branch..." 


    !allocate nca_delta
    allocate(NcaDeltaAnd_tau(Nspin,Nspin,Norb,Norb,0:Ltau))
    allocate(NcaDeltaAnd_iw(Nspin,Nspin,Norb,Norb,Lmats))

    !allocate atomic propagators
    allocate(ncaR0(Nhilbert,Nhilbert,0:Ltau))
    allocate(ncaR(Nhilbert,Nhilbert,0:Ltau))

    !allocate observables:
    allocate(nca_dens(Norb))
    allocate(nca_dens_up(Norb))
    allocate(nca_dens_dw(Norb))
    allocate(nca_docc(Norb))


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
  subroutine set_Hloc_nn(Hloc)
    complex(8),dimension(:,:,:,:) :: Hloc ![Nspin][Nspin][Norb][Norb]
    call assert_shape(Hloc,[Nspin,Nspin,Norb,Norb],"set_Hloc_nn","Hloc")
    !
    impHloc = Hloc
    !
    write(LOGfile,"(A)")"Set Hloc: done"
    call print_Hloc(impHloc)
  end subroutine set_Hloc_nn

  subroutine set_Hloc_so(hloc)
    complex(8),dimension(:,:) :: Hloc ![Nspin*Norb][Nspin*Norb]
    integer :: iorb,jorb
    integer :: ispin,jspin
    integer :: is,js
    call assert_shape(Hloc,[Nspin*Norb,Nspin*Norb],"set_Hloc_so","Hloc")
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb 
                js = jorb + (jspin-1)*Norb 
                impHloc(ispin,jspin,iorb,jorb) = Hloc(is,js)
             enddo
          enddo
       enddo
    enddo
    !
    write(LOGfile,"(A)")"Set Hloc: done"
    call print_Hloc(impHloc)
  end subroutine set_Hloc_so



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_Hloc(hloc,file)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: hloc
    character(len=*),optional                   :: file
    integer                                     :: ilat,jlat
    integer                                     :: iorb,jorb
    integer                                     :: ispin,jspin
    integer                                     :: unit
    character(len=32)                           :: fmt
    !
    unit=LOGfile;
    !
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc to file :"//reg(file)
    endif
    write(fmt,"(A,I0,A)")"(",Nspin*Norb,"A)"
    do ispin=1,Nspin
       do iorb=1,Norb
          write(unit,fmt)((str(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
       enddo
    enddo
    write(unit,*)""
    if(present(file))close(unit)
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
  !|ImpUP,BathUP>|ImpDW,BathDW >
  !+------------------------------------------------------------------+
  subroutine build_sector(isector,Hup)
    integer                                      :: isector
    type(sector_map)                             :: Hup
    integer                                      :: nup,ndw,sz,nt,twoJz
    integer                                      :: nup_,ndw_,sz_,nt_
    integer                                      :: twoSz_,twoLz_
    integer                                      :: i,ibath,iorb
    integer                                      :: iup,idw
    integer                                      :: dim
    integer                                      :: ivec(Ns),jvec(Ns)
    nup = getNup(isector)
    ndw = getNdw(isector)
    dim = getDim(isector)
    call map_allocate(Hup,dim)
    dim=0
    do idw=0,2**Ns-1
       jvec  = bdecomp(idw,Ns)
       ndw_  = sum(jvec)
       if(ndw_ /= ndw)cycle
       do iup=0,2**Ns-1
          ivec  = bdecomp(iup,Ns)
          nup_  = sum(ivec)
          if(nup_ /= nup)cycle
          dim      = dim+1
          Hup%map(dim) = iup + idw*2**Ns
       enddo
    enddo
  end subroutine build_sector





  subroutine delete_sector(isector,Hup)!,Hdw)
    integer                   :: isector
    type(sector_map)          :: Hup
    ! type(sector_map),optional :: Hdw
    call map_deallocate(Hup)
  end subroutine delete_sector






  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp






  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates 
  !   |out>=C_pos|in>  OR  |out>=C^+_pos|in> ; 
  !   the sign of |out> has the phase convention, pos labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c

  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
  end subroutine cdg




  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=C/Cm+|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine build_c_operator(ispin,iorb,Cmat)
    integer                              :: ispin,iorb
    real(8),dimension(Nhilbert,Nhilbert) :: Cmat
    integer                              :: imp,i,j
    integer                              :: ib(Nlevels)
    real(8)                              :: c_
    !build <j|C|i>
    imp = impIndex(ispin,iorb)
    Cmat=0d0
    do i=0,Nhilbert-1
       ib =  bdecomp(i,2*Ns)
       if(ib(imp)==0)cycle
       call c(imp,i,j,c_)
       Cmat(i+1,j+1)=c_
    enddo
    return
  end subroutine build_c_operator


  subroutine build_cdg_operator(ispin,iorb,Cmat)
    integer                              :: ispin,iorb
    real(8),dimension(Nhilbert,Nhilbert) :: Cmat
    integer                              :: imp,i,j
    integer                              :: ib(Nlevels)
    real(8)                              :: c_
    !build <j|C|i>
    imp = impIndex(ispin,iorb)
    Cmat=0d0
    do i=0,Nhilbert-1
       ib =  bdecomp(i,2*Ns)
       if(ib(imp)==1)cycle
       call cdg(imp,i,j,c_)
       Cmat(i+1,j+1)=c_
    enddo
    return
  end subroutine build_cdg_operator


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine nca_build_operators
    integer :: ispin,iorb,NN
    NN = Nhilbert
    do ispin=1,2
       do iorb=1,Norb
          allocate(Coperator(ispin,iorb)%Op(NN,NN))
          allocate(CDGoperator(ispin,iorb)%Op(NN,NN))
          call build_c_operator(ispin,iorb,Coperator(ispin,iorb)%Op(:,:))
          CDGoperator(ispin,iorb)%Op(:,:)=transpose(Coperator(ispin,iorb)%Op(:,:))
       enddo
    enddo
  end subroutine nca_build_operators




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
  !PURPOSE  : calculate the binomial factor n1 over n2
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
       xh = xh*dble(n1+1-i)/dble(i)
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
    allocate(espace(1:Nsectors))
    do isector=1,Nsectors
       dim=getdim(isector)
       allocate(espace(isector)%e(dim),espace(isector)%H(dim,dim))
    enddo
  end subroutine setup_eigenspace



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine reset_eigenspace
    integer :: isector
    forall(isector=1:Nsectors)
       espace(isector)%E=0.d0
       espace(isector)%H=0.d0
    end forall
  end subroutine reset_eigenspace


END MODULE NCA_SETUP
