MODULE NCA_GREENS_FUNCTIONS
  USE NCA_VARS_GLOBAL
  USE NCA_SETUP
  !
  USE SF_TIMER,    only: start_timer,stop_timer
  USE SF_LINALG,   only: diag,outerprod
  USE SF_CONSTANTS,only: one,xi,zero,pi
  USE SF_IOTOOLS,  only: free_unit,reg,str,splot
  USE SF_ARRAYS,   only: arange,linspace
  USE SF_SPECIAL,  only: fermi
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
  real(8),dimension(:),allocatable         :: wm,tau
  real(8),allocatable,dimension(:,:,:,:,:) :: impGtau

contains



  !+------------------------------------------------------------------+
  !PURPOSE  :
  ! G(w_n; E) = [iw_n − E]^−1 --> [f(E)-1]*exp(-tau*E)
  !+------------------------------------------------------------------+
  subroutine nca_build_bare_propagator
    integer                     :: iw,itau,istate
    !
    write(LOGfile,"(A)")"Evaluating bare atomic propagator R0(tau)"
    !
    call allocate_grids
    ncaR0=zero
    ncaR =zero
    !
    do istate=1,Nfock
       do itau=0,Ltau
          !this is: exp(-tau*Hloc) = sum_ab <b|exp(-tau*E_a)|a> |b><a| = sum_a exp(-tau*E_a) |a><a|
          ncaR0(istate,istate,itau) = exp(-tau(itau)*EigValues(istate))
          !this is the FFT of one/(iw - E_a):
          !ncaR0(istate,istate,itau) = exp(-tau(itau)*EigValues(istate))*(1d0-fermi(Eigvalues(istate),beta))
       enddo
       call splot("R0_tau.nca",tau,ncaR0(istate,istate,0:),append=.true.)
    enddo
    !
    ncaR = ncaR0
    !    
  end subroutine nca_build_bare_propagator





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !This routine performs the NCA self-consistency loop:
  !given R0(tau) and a guess for R(tau)
  !evaluates: 
  !Snca(tau) = \sum_ab [C_a * R(tau) * C^+_b . \Delta_{ab}(-tau) -&
  !                     C^+_a * R(tau) * C_b . \Delta_{ab}(tau)]
  !+------------------------------------------------------------------+
  subroutine nca_get_sigma_nca(ncaSigma)
    real(8),dimension(Nfock,Nfock,0:Ltau) :: ncaSigma
    real(8),dimension(Nfock,Nfock)        :: Matrix1,Matrix2
    integer                               :: itau,ispin,iorb
    integer :: i,j,k1
    ncaSigma=0d0
    !
    do itau=0,Ltau
       !
       do ispin=1,Nspin
          do iorb=1,Norb
             Matrix1=0d0
             Matrix2=0d0
             Matrix1 = matmul(matmul(Coperator(ispin,iorb)%Op,ncaR(:,:,itau)),CDGoperator(ispin,iorb)%Op)
             !
             Matrix2 = matmul(matmul(CDGoperator(ispin,iorb)%Op,ncaR(:,:,itau)),Coperator(ispin,iorb)%Op)
             !
             ncaSigma(:,:,itau) = ncaSigma(:,:,itau) - Matrix1*NcaDeltaAnd_tau(ispin,ispin,iorb,iorb,Ltau-itau) 
             ncaSigma(:,:,itau) = ncaSigma(:,:,itau) - Matrix2*NcaDeltaAnd_tau(ispin,ispin,iorb,iorb,itau)
             !
          enddo
       enddo
       !
    enddo

  end subroutine nca_get_sigma_nca





!+------------------------------------------------------------------+
!PURPOSE  : 
!Given R(tau),R0(tau) and S(tau) solve the Dyson equation:
!R(tau) = R0(tau) + int_0^tau d2 int_0_tau2 d1 R(tau-tau2)*S(tau2-tau1)*R0(tau1)
!+------------------------------------------------------------------+
subroutine nca_build_dressed_propagator
real(8),dimension(Nfock,Nfock,0:Ltau) :: ncaSigma,SxR0,ncaRold
real(8),dimension(Nfock,Nfock)        :: Matrix1
real(8)                                     :: dtau,error
integer                                     :: itau,is
logical                                     :: converged
converged=.false.
dtau=tau(1)-tau(0)
call start_timer()
do while (.not.converged)
 ncaRold = ncaR
 !
 !Get the Self-energy:
 call nca_get_sigma_nca(ncaSigma)
 !
 !Get the polarization SxR0(tau) = int_0^tau d1 ncaSigma(tau-tau1)R0(tau1)
 do itau=0,Ltau
    SxR0(:,:,itau) = tau_convolution(ncaSigma,ncaR0,itau)*dtau
 enddo
 !
 !R(tau) = R0(tau) + int_0^tau d1 R(tau-tau1)SxR0(tau1)
 do itau=0,Ltau
    Matrix1 = tau_convolution(ncaR,SxR0,itau)*dtau
    ncaR(:,:,itau) = ncaR0(:,:,itau) + Matrix1
 enddo
 !
 error = sum(abs(ncaR-ncaRold))/sum(abs(ncaR))
 !
 converged = (error<=nca_error)
 !
 write(LOGfile,*)error,converged
enddo
call stop_timer()
!
!
zeta_function = 0d0
do is=1,Nfock
 zeta_function=zeta_function+ncaR(is,is,Ltau)
end do
!
end subroutine nca_build_dressed_propagator




!+------------------------------------------------------------------+
!PURPOSE  :  G(tau) = -Tr(R(beta-tau)*C*R(tau)*C+)/Z
!+------------------------------------------------------------------+
subroutine nca_get_impurity_gf
integer                        :: ispin,iorb,is,itau,isector,dim,stride
real(8)                        :: trace
real(8),dimension(Nfock,Nfock) :: Matrix
real(8),dimension(0:3)         :: Ctail
!
if(allocated(impGtau))deallocate(impGtau)
allocate(impGtau(Nspin,Nspin,Norb,Norb,0:Ltau))
!
impGtau=0d0
do ispin=1,Nspin
 do iorb=1,Norb
    do itau=0,Ltau
       Matrix = 0d0
       Matrix =-matmul(ncaR(:,:,Ltau-itau),Coperator(ispin,iorb)%Op)
       Matrix = matmul(Matrix,ncaR(:,:,itau))
       Matrix = matmul(Matrix,CDGoperator(ispin,iorb)%Op)
       trace  = 0d0
       do is=1,Nfock
          trace=trace+Matrix(is,is)
       enddo
       impGtau(ispin,ispin,iorb,iorb,itau)=trace/zeta_function
    enddo
 enddo
enddo
!
do ispin=1,Nspin
 do iorb=1,Norb
    Ctail =  tail_coeff_glat(Uloc(iorb),nca_dens_spin(ispin,iorb),xmu,0d0)
    call fft_gf_tau2iw(impGmats(ispin,ispin,iorb,iorb,:),impGtau(ispin,ispin,iorb,iorb,0:Ltau),beta,C=Ctail)
 enddo
enddo
!
!
call Get_Self_Energy()
call print_imp_gf()
!
end subroutine nca_get_impurity_gf





!+------------------------------------------------------------------+
!PURPOSE  : Print normal Green's functions
!+------------------------------------------------------------------+
subroutine Get_Self_Energy
integer           :: i,ispin,iorb
complex(8)        :: fg0
real(8)           :: n0,C0,C1
real(8)           :: wm1,wm2
!
impSmats = zero
!
!Diagonal in both spin and orbital: this is ensured by the special *per impurity" bath structure
!no intra-orbital hoopings
do ispin=1,Nspin
 do iorb=1,Norb
    n0 = nca_dens(iorb)
    C0=Uloc(iorb)*(n0-0.5d0)
    C1=Uloc(iorb)**2*n0*(1.d0-n0)
    do i=1,Lmats
       fg0 = xi*wm(i) + xmu - impHloc(ispin,ispin,iorb,iorb) - ncaDeltaAnd_iw(ispin,ispin,iorb,iorb,i)
       impG0mats(ispin,ispin,iorb,iorb,i) = one/fg0
       impSmats(ispin,ispin,iorb,iorb,i)= fg0 - one/impGmats(ispin,ispin,iorb,iorb,i)
       impStail(ispin,ispin,iorb,iorb,i) = C0 + C1/(xi*wm(i))
    enddo
 enddo
enddo
!
wm1 = pi/beta ; wm2=3d0*pi/beta    
do ispin=1,Nspin
 do iorb=1,Norb
    nca_zeta(ispin,iorb)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ispin,ispin,iorb,iorb,1))/wm1 ))
 enddo
enddo
end subroutine Get_Self_Energy





!+------------------------------------------------------------------+
!PURPOSE  : Print normal Green's functions
!+------------------------------------------------------------------+
subroutine print_imp_gf
integer                                           :: i,ispin,iorb
character(len=20)                                 :: suffix
!
do ispin=1,Nspin
 do iorb=1,Norb
    suffix="_l"//str(iorb)//str(iorb)//"_s"//str(ispin)
    call splot("impSigma"//reg(suffix)//"_iw.nca",wm,impSmats(ispin,ispin,iorb,iorb,:))
    call splot("impStail"//reg(suffix)//"_iw.nca",wm,impStail(ispin,ispin,iorb,iorb,:))
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
function tau_convolution(A,B,itau) result(C)
real(8),dimension(Nfock,Nfock,0:Ltau) :: A,B
real(8),dimension(Nfock,Nfock)        :: C
real(8),dimension(Nfock,Nfock,0:Ltau) :: AxB
integer                                     :: itau
integer                                     :: itau1
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
real(8),dimension(Nfock,Nfock,0:Ltau),intent(in) :: f
integer,intent(in)                                     :: ia, ib
integer                                                :: k
real(8),dimension(Nfock,Nfock)                   :: sum
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
