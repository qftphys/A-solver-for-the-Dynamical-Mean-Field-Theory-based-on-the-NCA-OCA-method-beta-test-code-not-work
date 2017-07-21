MODULE NCA_SETUP
  USE NCA_VARS_GLOBAL
  !
  implicit none
  private


  public :: nca_init_structure
  !
  public :: setup_pointers
  public :: build_sector
  public :: delete_sector
  public :: bdecomp
  public :: binary_search
  !
  public :: c
  public :: cdg
  !
  public :: nca_build_operators
  !
  public :: build_dens_operator
  public :: build_docc_operator
  public :: build_sz2_operator


  complex(8),parameter :: zero=dcmplx(0d0,0d0)


contains





  !+------------------------------------------------------------------+
  !PURPOSE  : Init calculation
  !+------------------------------------------------------------------+
  subroutine nca_init_structure()
    integer :: NP,nup,ndw,iorb,jorb,ispin,jspin
    !
    Ns       = Norb
    !
    Nlevels  = 2*Ns
    !
    Nsectors = (Ns+1)*(Ns+1) !nup=0:Ns;ndw=0:Ns
    !    
    Nfock    = 2**Nlevels
    !
    nup=Ns/2
    ndw=Ns-nup
    NP=get_sector_dimension(nup,ndw)
    !
    write(LOGfile,"(A)")"Summary:"
    write(LOGfile,"(A,I15)") '# of levels/spin      = ',Ns
    write(LOGfile,"(A,I15)") 'Total system size     = ',Nlevels
    write(LOGfile,"(A,I15)") 'Fock space dimension  = ',Nfock
    write(LOGfile,"(A,I15)") 'Number of sectors     = ',Nsectors
    write(LOGfile,"(A,2I15)")'Largest Sector(s)     = ',NP
    write(LOGfile,"(A)")"--------------------------------------------"


    allocate(impHloc(Nspin,Nspin,Norb,Norb))
    impHloc=zero


    !Allocate pointers
    allocate(getdim(Nsectors),getnup(Nsectors),getndw(Nsectors))
    allocate(getsector(0:Ns,0:Ns))

    allocate(getCsector(2,Nsectors));getCsector=0
    allocate(getCDGsector(2,Nsectors));getCDGsector=0

    !Some check:
    write(LOGfile,*) "WARNING: this code only works for multi-orbital diagonal problems, no mixed-GF are evaluated"
    call sleep(1)
    if(Nspin>2)stop "Nspin > 2 ERROR. ask developer or develop your own on separate branch..."
    if(Norb>2)stop "Norb > 2 ERROR. ask developer or develop your own on separate branch..." 


    !allocate nca_delta
    allocate(NcaDeltaAnd_tau(2,2,Norb,Norb,0:Ltau));NcaDeltaAnd_tau=0d0
    allocate(NcaDeltaAnd_iw(2,2,Norb,Norb,Lmats));NcaDeltaAnd_iw=zero

    !allocate atomic propagators
    allocate(ncaR0(Nfock,Nfock,0:Ltau));ncaR0=0d0
    allocate(ncaR(Nfock,Nfock,0:Ltau));ncaR=0d0

    allocate(EigBasis(Nfock,Nfock));EigBasis=0d0
    allocate(EigValues(Nfock));EigValues=0d0

    !allocate observables:
    allocate(nca_dens_spin(2,Norb));nca_dens_spin=0d0
    allocate(nca_dens(Norb));nca_dens=0d0
    allocate(nca_docc(Norb));nca_docc=0d0
    allocate(nca_sz2(Norb));nca_sz2=0d0
    allocate(nca_zeta(Nspin,Norb));nca_zeta=0d0


    !allocate functions
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impStail(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impG0mats(Nspin,Nspin,Norb,Norb,Lmats))

    !allocate memory for operators
    allocate(Coperator(2,Norb))
    allocate(CDGoperator(2,Norb))

  end subroutine nca_init_structure









  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_pointers
    integer                          :: ispin,dim,isector,jsector
    integer                          :: nup,ndw,iorb,iup,jup,idw,jdw
    write(*,"(A)")"Setting up pointers:"
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

    getCsector=0
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup-1;jdw=ndw;if(jup < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw-1;if(jdw < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(2,isector)=jsector
    enddo
    !
    getCDGsector=0
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup+1;jdw=ndw;if(jup > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers




  !+------------------------------------------------------------------+
  !PURPOSE  : constructs the sectors by storing the map to the 
  !states i\in Hilbert_space from the states count in H_sector.
  !|ImpUP,BathUP>|ImpDW,BathDW >
  !+------------------------------------------------------------------+
  subroutine build_sector(isector,H)
    integer                                      :: isector
    type(sector_map)                             :: H
    integer                                      :: nup,ndw
    integer                                      :: nup_,ndw_
    integer                                      :: istate
    integer                                      :: dim
    integer                                      :: ivec(2*Ns)
    nup = getNup(isector)
    ndw = getNdw(isector)
    dim = getDim(isector)
    call map_allocate(H,dim)
    dim=0
    do istate=1,Nfock
       ivec=bdecomp(istate,2*Ns)
       nup_=sum(ivec(1:Ns))
       ndw_=sum(ivec(Ns+1:2*Ns))
       if(nup_/=nup .OR. ndw_/=ndw)cycle
       dim      = dim+1
       H%map(dim) = istate
    enddo
  end subroutine build_sector

  subroutine delete_sector(isector,H)
    integer                   :: isector
    type(sector_map)          :: H
    call map_deallocate(H)
  end subroutine delete_sector






  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition, corresponding to the number i-1.
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: i,Ntot
    integer :: ivec(Ntot)
    integer :: bit,istate
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state istate=[0,...,2^Ntot-1]
    !remark: to respect the binary structure 0=|0...00...0> to Nfock=|1...11...1>
    !one needs to btest the i-1 state to shift the 1...Nfock to 0...2^Ntot-1
    istate=i-1
    do bit=1,Ntot
       busy=btest(istate,bit-1)
       ivec(bit)=0
       if(busy)ivec(bit)=1
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
    integer               :: l,instate
    instate = in-1
    if(.not.btest(instate,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(instate,l-1))fsgn=-fsgn
    enddo
    out = ibclr(instate,pos-1) + 1
  end subroutine c

  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l,instate
    instate = in - 1
    if(btest(instate,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(instate,l-1))fsgn=-fsgn
    enddo
    out = ibset(instate,pos-1) + 1
  end subroutine cdg








  !+-------------------------------------------------------------------+
  !PURPOSE: c_{iorb,ispin} = \sum_ab <a|c_{iorb,ispin}|b> |a><b|
  !+-------------------------------------------------------------------+
  subroutine build_c_operator(ispin,iorb,Cmat)
    integer                        :: ispin,iorb
    real(8),dimension(Nfock,Nfock) :: Cmat
    real(8)                        :: c_,Dab
    integer                        :: imp,i,j,ia,ib,iia,iib
    integer                        :: ivec(Nlevels)
    !
    imp = iorb+(ispin-1)*Norb
    Cmat=0d0
    do ib=1,Nfock               !loop over |b> states
       do ia=1,Nfock            !loop over <a| states
          Dab=0d0
          do iib=1,Nfock        !now loop over the components of |b>=\sum_iib u^b_iib |iib>
             ivec = bdecomp(iib,2*Ns)
             if(ivec(imp) == 0) cycle
             call c(imp,iib,iia,c_)
             Dab = Dab + c_*EigBasis(ia,iia)*EigBasis(iib,ib)
          enddo
          Cmat(ia,ib) = Dab
       enddo
    enddo
  end subroutine build_c_operator


  subroutine build_cdg_operator(ispin,iorb,Cmat)
    integer                        :: ispin,iorb
    real(8),dimension(Nfock,Nfock) :: Cmat
    real(8)                        :: c_,Dab
    integer                        :: imp,i,j,ia,ib,iia,iib
    integer                        :: ivec(Nlevels)
    !
    imp = iorb+(ispin-1)*Norb
    !
    Cmat=0d0
    !
    do ib=1,Nfock               !loop over |b> states
       do ia=1,Nfock            !loop over <a| states
          Dab=0d0
          do iib=1,Nfock        !now loop over the components of |b>=\sum_iib u^b_iib |iib>
             ivec = bdecomp(iib,2*Ns)
             if(ivec(imp) == 1) cycle
             call cdg(imp,iib,iia,c_)
             Dab = Dab + c_*EigBasis(ia,iia)*EigBasis(iib,ib)
          enddo
          Cmat(ia,ib) = Dab
       enddo
    enddo
  end subroutine build_cdg_operator


  subroutine build_dens_operator(ispin,iorb,Cmat)
    integer                        :: ispin,iorb
    real(8),dimension(Nfock,Nfock) :: Cmat
    real(8)                        :: dens,Dab
    integer                        :: imp,i,j,ia,ib,iia,iib
    integer                        :: ivec(Nlevels)
    !
    imp = iorb+(ispin-1)*Norb
    !
    Cmat=0d0
    !
    do ia=1,Nfock
       Dab=0d0
       do iia=1,Nfock
          ivec = bdecomp(iia,2*Ns)
          dens = dble(ivec(imp))
          Dab  = Dab + dens*EigBasis(ia,iia)*EigBasis(iia,ia)
       enddo
       Cmat(ia,ia) =  Dab
    enddo
  end subroutine build_dens_operator



  subroutine build_docc_operator(iorb,Cmat)
    integer                        :: iorb
    real(8),dimension(Nfock,Nfock) :: Cmat
    real(8)                        :: dens_up,dens_dw,docc,Dab
    integer                        :: imp_up,imp_dw,ia,ib,iia,iib
    integer                        :: ivec(Nlevels)
    !
    imp_up = iorb
    imp_dw = iorb + Ns
    !
    Cmat=0d0
    !
    do ia=1,Nfock
       Dab=0d0
       do iia=1,Nfock
          ivec = bdecomp(iia,2*Ns)
          dens_up = dble(ivec(imp_up))
          dens_dw = dble(ivec(imp_dw))
          docc    = dens_up*dens_dw
          Dab  = Dab + docc*EigBasis(ia,iia)*EigBasis(iia,ia)
       enddo
       Cmat(ia,ia) =  Dab
    enddo
  end subroutine build_docc_operator



  subroutine build_sz2_operator(iorb,Cmat)
    integer                        :: iorb
    real(8),dimension(Nfock,Nfock) :: Cmat
    real(8)                        :: dens_up,dens_dw,sz,sz2,Dab
    integer                        :: imp_up,imp_dw,ia,ib,iia,iib
    integer                        :: ivec(Nlevels)
    !
    imp_up = iorb
    imp_dw = iorb + Ns
    !
    Cmat=0d0
    !
    do ia=1,Nfock
       Dab=0d0
       do iia=1,Nfock
          ivec = bdecomp(iia,2*Ns)
          dens_up = dble(ivec(imp_up))
          dens_dw = dble(ivec(imp_dw))
          sz      = 0.5d0*(dens_up-dens_dw)
          sz2     = sz**2
          Dab  = Dab + sz2*EigBasis(ia,iia)*EigBasis(iia,ia)
       enddo
       Cmat(ia,ia) =  Dab
    enddo
  end subroutine build_sz2_operator







  !+------------------------------------------------------------------+
  !Purpose  : 
  !+------------------------------------------------------------------+
  subroutine nca_build_operators
    integer                              :: ispin,iorb,i,j,isector
    integer                              :: ivec(2*Ns),stride,dim
    do ispin=1,2
       do iorb=1,Norb
          allocate(Coperator(ispin,iorb)%Op(Nfock,Nfock))
          allocate(CDGoperator(ispin,iorb)%Op(Nfock,Nfock))
          Coperator(ispin,iorb)%Op(:,:)=0d0
          CDGoperator(ispin,iorb)%Op=0d0
          !
          call build_c_operator(ispin,iorb,Coperator(ispin,iorb)%Op)
          call build_cdg_operator(ispin,iorb,CDGoperator(ispin,iorb)%Op)
          ! CDGoperator(ispin,iorb)%Op = transpose(Coperator(ispin,iorb)%Op)
          !
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



END MODULE NCA_SETUP
