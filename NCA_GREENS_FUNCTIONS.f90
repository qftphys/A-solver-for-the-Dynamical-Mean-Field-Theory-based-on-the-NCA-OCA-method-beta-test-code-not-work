MODULE NCA_GREENS_FUNCTIONS
  USE NCA_VARS_GLOBAL
  USE NCA_SETUP
  !
  USE SF_CONSTANTS,only: one,xi,zero,pi
  USE SF_IOTOOLS,  only: free_unit,reg,str,splot
  USE SF_ARRAYS,   only: arange,linspace
  USE DMFT_FFTGF
  USE DMFT_FFTAUX
  !
  implicit none
  private 

  public :: nca_build_bare_propagator
  public :: nca_build_dressed_propagator
  public :: nca_get_impurity_gf


  
  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable     :: wm,tau

  !Strong coupling propagators R_{0,nca}(tau) & R_nca(tau): (2*Ntot,2*Ntot,0:Ltau)
  !=========================================================
  real(8),allocatable,dimension(:,:,:,:,:) :: impGtau




contains



  !>TODO:
  !1. Build up the NCA_R0 exploiting its block diagonal structure, i.e.
  ! evaliuate <j|exp(-tau*H)|i> for i,j in a given sector, and add up the
  ! results in a large Nhilbert**2 matrix nca_R0.
  !2.nca_R0, nca_Self, nca_R are the only large, i.e. Nhilbert**2, matrices to keep
  !3. Collapse the C,CDG operators into two Nhilbert vectors containing actual
  ! values and column indices of the very sparse operators. Overload the
  ! matmul operation to deal with these structure (type) to perform
  ! the mat-mat product in nca_get_sigma_nca.
  !
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine nca_build_bare_propagator
    integer                              :: is,itau,isector,stride,dim
    real(8),dimension(Nhilbert)          :: Eval
    real(8),dimension(Nhilbert,Nhilbert) :: Umat
    stride= 0
    Umat  = 0d0
    Eval  = 0d0
    call allocate_grids
    ncaR0=zero
    ncaR =zero
    do isector = 1 , Nsectors
       dim = getdim(isector)
       Umat(stride+1:stride+dim,stride+1:stride+dim) = espace(isector)%H(1:dim,1:dim)
       Eval(stride+1:stride+dim) = espace(isector)%E(1:dim)
       stride = stride+dim
    enddo
    forall(is=1:Nhilbert,itau=0:Ltau)ncaR0(is,is,itau) = exp(-tau(itau)*Eval(is))
    do itau=0,Ltau
       ncaR0(:,:,itau) = matmul(matmul(Umat,ncaR0(:,:,itau)),transpose(Umat))
    enddo
    ncaR = ncaR0
  end subroutine nca_build_bare_propagator





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !This routine performs the NCA self-consistency loop:
  !given R0(tau) and a guess for R(tau)
  !evaluates: 
  !Snca(tau) = \sum_ab [C_a * R(tau) * C^+_b * \Delta_{ab}(-tau) -&
  !                     C^+_a * R(tau) * C_b * \Delta_{ab}(tau)]
  !+------------------------------------------------------------------+
  subroutine nca_get_sigma_nca(ncaSigma)
    real(8),dimension(Nhilbert,Nhilbert,0:Ltau) :: ncaSigma
    real(8),dimension(Nhilbert,Nhilbert)        :: Matrix1,Matrix2
    integer                                     :: itau,ispin,iorb
    ncaSigma=0d0
    Matrix1=0d0
    Matrix2=0d0
    do itau=0,Ltau
       do ispin=1,Nspin
          do iorb=1,Norb
             Matrix1 = matmul(Coperator(ispin,iorb)%Op,ncaR(:,:,itau))
             Matrix1 = matmul(Matrix1,CDGoperator(ispin,iorb)%Op)
             Matrix1 =-Matrix1*NcaDeltaAnd_tau(ispin,ispin,iorb,iorb,Ltau-itau)
             !
             Matrix2 = matmul(CDGoperator(ispin,iorb)%Op,ncaR(:,:,itau))
             Matrix2 = matmul(Matrix2,Coperator(ispin,iorb)%Op)
             Matrix2 = Matrix2*NcaDeltaAnd_tau(ispin,ispin,iorb,iorb,itau)
             !
             ncaSigma(:,:,itau) = Matrix1 - Matrix2
          enddo
       enddo
    enddo
  end subroutine nca_get_sigma_nca





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !Given R(tau),R0(tau) and S(tau) solve the Dyson equation:
  !R(tau) = R0(tau) + int_0^tau d2 int_0_tau2 d1 R(tau-tau2)*S(tau2-tau1)*R0(tau1)
  !+------------------------------------------------------------------+
  subroutine nca_build_dressed_propagator
    real(8),dimension(Nhilbert,Nhilbert,0:Ltau) :: ncaSigma,SxR0,ncaRold
    real(8),dimension(Nhilbert,Nhilbert)        :: Matrix1
    real(8)                                     :: dtau,error
    integer                                     :: itau,is
    logical                                     :: converged
    converged=.false.
    dtau=tau(1)-tau(0)
    do while (.not.converged)
       ncaRold = ncaR
       !Get the Self-energy:
       call nca_get_sigma_nca(ncaSigma)
       !Get the polarization SxR0(tau) = int_0^tau d1 ncaSigma(tau-tau1)R0(tau1)
       do itau=0,Ltau
          SxR0(:,:,itau) = tau_convolution(ncaSigma,ncaR0,itau,Ltau)*dtau
       enddo
       !R(tau) = R0(tau) + int_0^tau d1 R(tau-tau1)SxR0(tau1)
       do itau=0,Ltau
          Matrix1 = tau_convolution(ncaR,SxR0,itau,Ltau)*dtau
          ncaR(:,:,itau) = ncaR0(:,:,itau) + Matrix1
       enddo
       error = sum(abs(ncaR(1,1,:)-ncaRold(1,1,:)))/sum(abs(ncaR(1,1,:)))
       converged = (error<=nca_error)
       print*,error,converged
    enddo
    !
    !
    zeta_function = 0d0
    do is=1,Nhilbert
       zeta_function=zeta_function+ncaR(is,is,Ltau)
    end do
    !
  end subroutine nca_build_dressed_propagator




  !+------------------------------------------------------------------+
  !PURPOSE  :  G(tau) = -Tr(R(beta-tau)*C*R(tau)*C+)/Z
  !+------------------------------------------------------------------+
  subroutine nca_get_impurity_gf
    integer                              :: ispin,iorb,is,itau
    real(8)                              :: trace
    real(8),dimension(Nhilbert,Nhilbert) :: Matrix
    real(8),dimension(0:3)               :: Ctail
    !
    zeta_function=0d0
    do is=1,Nhilbert
       zeta_function = zeta_function + ncaR(is,is,Ltau)
    end do
    !
    if(allocated(impGtau))deallocate(impGtau)
    allocate(impGtau(Nspin,Nspin,Norb,Norb,0:Ltau))
    !
    impGtau=0d0
    do itau=0,Ltau
       do ispin=1,Nspin
          do iorb=1,Norb
             Matrix =-matmul(ncaR(:,:,Ltau-itau),Coperator(ispin,iorb)%Op)
             Matrix = matmul(Matrix,ncaR(:,:,itau))
             Matrix = matmul(Matrix,CDGoperator(ispin,iorb)%Op)
             trace  = 0d0
             do is=1,Nhilbert
                trace=trace+Matrix(is,is)
             enddo
             impGtau(ispin,ispin,iorb,iorb,itau)=trace/zeta_function
          enddo
       enddo
    enddo
    !
    do ispin=1,Nspin
       do iorb=1,Norb
          Ctail =  tail_coeff_glat(Uloc(iorb),nca_dens(iorb)/2d0,xmu,0d0)
          call fft_gf_tau2iw(impGmats(ispin,ispin,iorb,iorb,:),impGtau(ispin,ispin,iorb,iorb,0:Ltau),beta,C=Ctail)
       enddo
    enddo
    !
    ! !
    call print_imp_gf
    !
  end subroutine nca_get_impurity_gf





  !+------------------------------------------------------------------+
  !PURPOSE  : Print normal Green's functions
  !+------------------------------------------------------------------+
  subroutine print_imp_gf
    integer                                           :: i,ispin,iorb
    complex(8)                                        :: fg0
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: impG0mats,Stail
    character(len=20)                                 :: suffix
    real(8)                                           :: n0,C0,C1
    !
    impSmats = zero
    n0 = nca_dens(1)
    C0=Uloc(1)*(n0-0.5d0)
    C1=Uloc(1)**2*n0*(1.d0-n0)
    !
    !Diagonal in both spin and orbital: this is ensured by the special *per impurity" bath structure
    !no intra-orbital hoopings
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Lmats
             fg0 = xi*wm(i) + xmu - impHloc(ispin,ispin,iorb,iorb) - ncaDeltaAnd_iw(ispin,ispin,iorb,iorb,i)
             impG0mats(ispin,ispin,iorb,iorb,i) = one/fg0
             impSmats(ispin,ispin,iorb,iorb,i)= fg0 - one/impGmats(ispin,ispin,iorb,iorb,i)
             Stail(ispin,ispin,iorb,iorb,i) = C0 + C1/(xi*wm(i))
          enddo
       enddo
    enddo
    !

    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//str(iorb)//str(iorb)//"_s"//str(ispin)
          call splot("impSigma"//reg(suffix)//"_iw.nca",wm,impSmats(ispin,ispin,iorb,iorb,:))
          call splot("impStail"//reg(suffix)//"_iw.nca",wm,Stail(ispin,ispin,iorb,iorb,:))
          call splot("impG"//reg(suffix)//"_iw.nca"    ,wm,impGmats(ispin,ispin,iorb,iorb,:))
          call splot("impG0"//reg(suffix)//"_iw.nca"   ,wm,impG0mats(ispin,ispin,iorb,iorb,:))
          call splot("impG"//reg(suffix)//"_tau.nca",tau,impGtau(ispin,ispin,iorb,iorb,:))
       enddo
    enddo
  end subroutine print_imp_gf








  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate and setup arrays with frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    allocate(wm(Lmats))
    allocate(tau(0:Ltau))
    wm     = pi/beta*dble(2*arange(1,Lmats)-1)
    tau(0:)= linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function tau_convolution(A,B,itau,L) result(C)
    integer                                  :: L
    real(8),dimension(Nhilbert,Nhilbert,0:L) :: A,B
    real(8),dimension(Nhilbert,Nhilbert)     :: C
    integer                                  :: itau
    real(8),dimension(Nhilbert,Nhilbert,0:L) :: AxB
    integer                                  :: itau1
    AxB=0d0
    do itau1=0,itau
       AxB(:,:,itau1)=matmul(A(:,:,itau-itau1),B(:,:,itau1))
    end do
    C = tau_trapz(AxB(:,:,0:),0,itau)
  end function tau_convolution


  !----------------------------------------------------------------------------
  !  This function calculates the sum
  !    \sum_{k=i}^{j} w_{k}^{i,j} f_{k}.
  !    w_{k}^{i,j} = 1/2  for k=i,j
  !                = 1    for i<k<j
  !                      OR
  !                = 1    for i<k<=j (half-edge)
  !----------------------------------------------------------------------------
  function tau_trapz(f,ia,ib) result(sum)
    real(8),dimension(Nhilbert,Nhilbert,0:Ltau),intent(in) :: f
    integer,intent(in)                                     :: ia, ib
    integer                                                :: k
    real(8),dimension(Nhilbert,Nhilbert)                   :: sum
    sum=0d0
    if(ia==ib)then
       return
    else
       sum=sum+0.5d0*f(:,:,ia)
       do k=ia+1,ib-1
          sum=sum+f(:,:,k)
       end do
       sum=sum+0.5d0*f(:,:,ib)
    end if
  end function tau_trapz




end MODULE NCA_GREENS_FUNCTIONS


! 
