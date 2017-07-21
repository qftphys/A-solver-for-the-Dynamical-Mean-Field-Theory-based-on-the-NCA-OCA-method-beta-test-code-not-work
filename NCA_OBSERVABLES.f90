MODULE NCA_OBSERVABLES
  USE NCA_VARS_GLOBAL
  USE NCA_SETUP
  !
  USE SF_IOTOOLS, only:str,free_unit
  !
  implicit none
  private 

  public :: nca_get_observables

  logical,save                       :: iolegend=.true.
  real(8),dimension(:),allocatable   :: nca_magz
  
contains


  subroutine nca_get_observables()
    integer                              :: ispin,iorb
    real(8),dimension(Nfock,Nfock) :: Matrix,Mat_up,Mat_dw
    !
    !
    nca_dens_spin = 0d0
    nca_dens      = 0d0
    nca_docc      = 0d0
    nca_sz2       = 0d0
    !
    !
    !Get N_up
    do iorb=1,Norb
       call build_dens_operator(1,iorb,Matrix)
       call nca_build_observable(Matrix,nca_dens_spin(1,iorb))
    enddo
    !
    !Get N_dw
    do iorb=1,Norb
       call build_dens_operator(2,iorb,Matrix)
       call nca_build_observable(Matrix,nca_dens_spin(2,iorb))
    enddo
    !
    !Get N=N_up + N_dw
    nca_dens = sum(nca_dens_spin,dim=1)
    !
    !Get double occupancy
    do iorb=1,Norb
       call build_docc_operator(iorb,Matrix)
       call nca_build_observable(Matrix,nca_docc(iorb))
    enddo

    do iorb=1,Norb
       call build_sz2_operator(iorb,Matrix)
       call nca_build_observable(Matrix,nca_sz2(iorb))
    enddo
    !
    !
    if(.not.allocated(nca_magz))allocate(nca_magz(Norb))
    nca_magz = nca_dens_spin(1,:) - nca_dens_spin(2,:)
    !
    write(LOGfile,"(A,10f18.12)")       "dens UP=",(nca_dens_spin(1,iorb),iorb=1,Norb)
    write(LOGfile,"(A,10f18.12)")       "dens DW=",(nca_dens_spin(2,iorb),iorb=1,Norb)
    write(LOGfile,"(A,10f18.12,f18.12)")"dens   =",(nca_dens(iorb),iorb=1,Norb),sum(nca_dens)
    write(LOGfile,"(A,10f18.12)")       "docc   =",(nca_docc(iorb),iorb=1,Norb)
    if(Nspin>1)write(LOGfile,"(A,10f18.12,A)")"mag    =",(nca_magz(iorb),iorb=1,Norb)
    !
    if(iolegend)call write_legend
    call write_observables()
    !
  end subroutine nca_get_observables


  subroutine nca_build_observable(Matrix,observable)
    real(8),dimension(Nfock,Nfock) :: Matrix
    real(8)                              :: observable
    integer :: is
    !
    observable=0d0
    !
    Matrix = matmul(ncaR(:,:,Ltau),Matrix)
    do is=1,Nfock
       observable= observable + Matrix(is,is)/zeta_function
    enddo
    !
  end subroutine nca_build_observable



  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_info.nca")
    write(unit,"(A1,90(A10,6X))")"#",&
         (str(iorb)//"dens_"//str(iorb),iorb=1,Norb),&
         (str(Norb+iorb)//"docc_"//str(iorb),iorb=1,Norb),&
         (str(2*Norb+iorb)//"nup_"//str(iorb),iorb=1,Norb),&
         (str(3*Norb+iorb)//"ndw_"//str(iorb),iorb=1,Norb),&
         (str(4*Norb+iorb)//"mag_"//str(iorb),iorb=1,Norb),&
         (str(5*Norb+iorb)//"sz2_"//str(iorb),iorb=1,Norb),&
         ((str(6*Norb+iorb)//"z_"//str(iorb)//"s"//str(ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="parameters_info.nca")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (str(2+iorb)//"U_"//str(iorb),iorb=1,Norb),&
         str(2+Norb+1)//"U'",str(2+Norb+2)//"Jh"
    close(unit)
    !
    iolegend=.false.
  end subroutine write_legend


  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    integer :: unit
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_all.nca",position='append')
    write(unit,"(90(F15.9,1X))")&
         (nca_dens(iorb),iorb=1,Norb),&
         (nca_docc(iorb),iorb=1,Norb),&
         (nca_dens_spin(1,iorb),iorb=1,Norb),&
         (nca_dens_spin(2,iorb),iorb=1,Norb),&
         (nca_magz(iorb),iorb=1,Norb),&
         (nca_sz2(iorb),iorb=1,Norb),&
         ((nca_zeta(ispin,iorb),iorb=1,Norb),ispin=1,Nspin)
    close(unit)    
    !
    unit = free_unit()
    open(unit,file="parameters_last.nca")
    write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
    close(unit)
    !
    unit = free_unit()
    open(unit,file="observables_last.nca")
    write(unit,"(90(F15.9,1X))")&
         (nca_dens(iorb),iorb=1,Norb),&
         (nca_docc(iorb),iorb=1,Norb),&
         (nca_dens_spin(1,iorb),iorb=1,Norb),&
         (nca_dens_spin(2,iorb),iorb=1,Norb),&
         (nca_magz(iorb),iorb=1,Norb),&
         (nca_sz2(iorb),iorb=1,Norb),&
         ((nca_zeta(ispin,iorb),iorb=1,Norb),ispin=1,Nspin)
    close(unit)         
  end subroutine write_observables


END MODULE NCA_OBSERVABLES










! subroutine nca_get_observables()
!   integer                              :: ispin,iorb,is,itau
!   real(8),dimension(Nfock,Nfock) :: Matrix
!   real(8),dimension(Nfock)          :: Rbeta
!   real(8),dimension(Norb)              :: dens
!   !
!   do ispin=1,Nspin
!      nca_dens = 0d0
!      do iorb=1,Norb
!         Matrix = matmul(CDGoperator(ispin,iorb)%Op,Coperator(ispin,iorb)%Op)
!         do is=1,Nfock
!            nca_dens(iorb)= nca_dens(iorb) + Rbeta(is)*Matrix(is,is)/zeta_function
!         enddo
!      enddo
!   enddo
!   print*,"Z     =",zeta_function
!   print*,"Nimp  =",2.d0*nca_dens
!   !
!   do ispin=1,Nspin
!      dens = 0d0
!      do iorb=1,Norb
!         Matrix = matmul(CDGoperator(ispin,iorb)%Op,Coperator(ispin,iorb)%Op)
!         Matrix = matmul(ncaR(:,:,Ltau),Matrix)
!         do is=1,Nfock
!            dens(iorb)= dens(iorb) + Matrix(is,is)/zeta_function
!         enddo
!      enddo
!   enddo
!   print*,"Ndens =",2d0*dens
!   !
! end subroutine nca_get_observables
