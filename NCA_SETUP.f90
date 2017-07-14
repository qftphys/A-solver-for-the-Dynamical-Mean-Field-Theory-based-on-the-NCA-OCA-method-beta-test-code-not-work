MODULE NCA_SETUP
  USE NCA_VARS_GLOBAL
  !
  USE SF_CONSTANTS, only:zero
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,str
  USE SF_MISC, only: assert_shape
  !
  implicit none
  private

  interface print_state_vector
     module procedure print_state_vector_ivec
  end interface print_state_vector



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
  public :: nca_build_operators

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
    Nfock = 2**Nlevels
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
    allocate(NcaDeltaAnd_tau(Nspin,Nspin,Norb,Norb,0:Ltau));NcaDeltaAnd_tau=0d0
    allocate(NcaDeltaAnd_iw(Nspin,Nspin,Norb,Norb,Lmats));NcaDeltaAnd_iw=zero

    !allocate atomic propagators
    allocate(ncaR0(Nfock,Nfock,0:Ltau));ncaR0=0d0
    allocate(ncaR(Nfock,Nfock,0:Ltau));ncaR=0d0

    allocate(EigBasis(Nfock,Nfock));EigBasis=0d0
    allocate(EigValues(Nfock));EigValues=0d0

    !allocate observables:
    allocate(nca_dens(Norb));nca_dens=0d0
    allocate(nca_dens_up(Norb));nca_dens_up=0d0
    allocate(nca_dens_dw(Norb));nca_dens_dw=0d0
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
  subroutine build_sector(isector,Hup)
    integer                                      :: isector
    type(sector_map)                             :: Hup
    integer                                      :: nup,ndw
    integer                                      :: nup_,ndw_
    integer                                      :: i
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





  subroutine delete_sector(isector,Hup)
    integer                   :: isector
    type(sector_map)          :: Hup
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
    do l=1,Ntot
       busy=btest(i,l-1)
       ivec(l)=0
       if(busy)ivec(l)=1
    enddo
  end function bdecomp



  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Nlevels) with the binary sequence 
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bjoin(ib,Ntot) result(i)
    integer                 :: Ntot
    integer,dimension(Ntot) :: ib
    integer                 :: i,j
    i=0
    do j=0,Ntot-1
       i=i+ib(j+1)*2**j
    enddo
  end function bjoin





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
  !PURPOSE  : build <j|C|i> with j,i \in FockBasis
  !+-------------------------------------------------------------------+
  subroutine build_c_operator_fock(ispin,iorb,Cmat)
    integer                              :: ispin,iorb
    real(8),dimension(Nfock,Nfock) :: Cmat
    integer                              :: imp,i,j,iup,idw
    integer                              :: ivec(Nlevels)
    real(8)                              :: c_,c_op,c_sgn
    imp = iorb+(ispin-1)*Ns
    Cmat=0d0
    do idw=0,2**Ns-1
       do iup=0,2**Ns-1
          i = iup + idw*2**Ns
          ivec = bdecomp(i,2*Ns)
          if(ivec(imp)==0)cycle
          call c(imp,i,j,c_)
          Cmat(i+1,j+1) = c_
       enddo
    enddo
  end subroutine build_c_operator_fock

  subroutine build_cdg_operator_fock(ispin,iorb,Cmat)
    integer                              :: ispin,iorb
    real(8),dimension(Nfock,Nfock) :: Cmat
    integer                              :: imp,i,j,iup,idw
    integer                              :: ivec(Nlevels)
    real(8)                              :: c_
    imp = iorb+(ispin-1)*Ns
    Cmat=0d0
    do idw=0,2**Ns-1
       do iup=0,2**Ns-1
          i = iup + idw*2**Ns
          ivec = bdecomp(i,2*Ns)
          if(ivec(imp)==1)cycle
          call cdg(imp,i,j,c_)
          Cmat(i+1,j+1) = c_
       enddo
    enddo
  end subroutine build_cdg_operator_fock





  !+-------------------------------------------------------------------+
  !PURPOSE:
  !+-------------------------------------------------------------------+
  subroutine build_c_operator_eig(ispin,iorb,Cmat)
    integer                        :: ispin,iorb
    real(8),dimension(Nfock,Nfock) :: Cmat
    real(8),dimension(Nfock)       :: avec,bvec
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
          do iib=0,Nfock-1        !now loop over the components of |b>=\sum_iib u^b_iib |iib>
             ivec = bdecomp(iib,2*Ns)
             if(ivec(imp) == 0) cycle
             call c(imp,iib,iia,c_)
             Dab = Dab + EigBasis(ia,iia+1)*c_*EigBasis(iib+1,ib)
          enddo
          Cmat(ia,ib) = Cmat(ia,ib) + Dab
       enddo
    enddo
  end subroutine build_c_operator_eig


  subroutine build_cdg_operator_eig(ispin,iorb,Cmat)
    integer                        :: ispin,iorb
    real(8),dimension(Nfock,Nfock) :: Cmat
    real(8),dimension(Nfock)       :: avec,bvec
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
          do iib=0,Nfock-1        !now loop over the components of |b>=\sum_iib u^b_iib |iib>
             ivec = bdecomp(iib,2*Ns)
             if(ivec(imp) == 1) cycle
             call cdg(imp,iib,iia,c_)
             Dab = Dab + EigBasis(ia,iia+1)*c_*EigBasis(iib+1,ib)
          enddo
          Cmat(ia,ib) = Cmat(ia,ib) + Dab
       enddo
    enddo
  end subroutine build_cdg_operator_eig




  !+------------------------------------------------------------------+
  !Purpose  : 
  !+------------------------------------------------------------------+
  subroutine nca_build_operators
    integer                              :: ispin,iorb,i,j,isector
    integer                              :: ivec(2*Ns),stride,dim

    do ispin=1,2
       do iorb=1,Norb
          allocate(Coperator(ispin,iorb)%Op(Nfock,Nfock))
          Coperator(ispin,iorb)%Op(:,:)=0d0
          call build_c_operator_eig(ispin,iorb,Coperator(ispin,iorb)%Op)
          !
          allocate(CDGoperator(ispin,iorb)%Op(Nfock,Nfock))
          CDGoperator(ispin,iorb)%Op=0d0
          call build_cdg_operator_eig(ispin,iorb,CDGoperator(ispin,iorb)%Op)
       enddo
    enddo
    !
    ! ispin=1
    ! do iorb=1,Norb
    !    write(*,*)""
    !    write(*,*)"C_l"//str(iorb)//"_s"//str(ispin)
    !    do i=1,Nfock
    !       ivec = bdecomp(i-1,2*Ns)
    !       write(800,"(I9,A1,100I1)",advance="no")i-1," | ",(ivec(j),j=1,2*Ns)
    !       write(800,"(1000F9.2)")(Coperator(ispin,iorb)%Op(i,j),j=1,Nfock)
    !    enddo
    !    write(*,*)""
    !    write(*,*)"Cdg_l"//str(iorb)//"_s"//str(ispin)
    !    do i=1,Nfock
    !       ivec = bdecomp(i-1,2*Ns)
    !       write(801,"(I9,A1,100I1)",advance="no")i-1," | ",(ivec(j),j=1,2*Ns)
    !       write(801,"(1000F9.2)")(CDGoperator(ispin,iorb)%Op(i,j),j=1,Nfock)
    !    enddo
    !    write(*,*)""
    ! enddo
    ! !
    ! ispin=2
    ! do iorb=1,Norb
    !    write(*,*)""
    !    write(*,*)"C_l"//str(iorb)//"_s"//str(ispin)
    !    do i=1,Nfock
    !       ivec = bdecomp(i-1,2*Ns)
    !       write(802,"(I9,A1,100I1)",advance="no")i-1," | ",(ivec(j),j=1,2*Ns)
    !       write(802,"(1000F9.2)")(Coperator(ispin,iorb)%Op(i,j),j=1,Nfock)
    !    enddo
    !    write(*,*)""
    !    write(*,*)"Cdg_l"//str(iorb)//"_s"//str(ispin)
    !    do i=1,Nfock
    !       ivec = bdecomp(i-1,2*Ns)
    !       write(803,"(I9,A1,100I1)",advance="no")i-1," | ",(ivec(j),j=1,2*Ns)
    !       write(803,"(1000F9.2)")(CDGoperator(ispin,iorb)%Op(i,j),j=1,Nfock)
    !    enddo
    !    write(*,*)""
    ! enddo


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
  !PURPOSE  : print a state vector |{up}>|{dw}>
  !+------------------------------------------------------------------+
  subroutine print_state_vector_ivec(ivec,unit)
    integer,intent(in) :: ivec(:)
    integer,optional   :: unit
    integer            :: unit_
    integer            :: i,j,Ntot
    character(len=2)   :: fbt
    character(len=16)  :: fmt
    unit_=6;if(present(unit))unit_=unit
    Ntot = size(ivec)
    i= bjoin(ivec,Ntot)
    write(unit_,"(I9,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="yes")">"
  end subroutine print_state_vector_ivec
  !



END MODULE NCA_SETUP
