module NCA_DIAG
  USE NCA_VARS_GLOBAL
  USE NCA_SETUP
  !
  USE SF_CONSTANTS,only: one,xi,zero,pi
  USE SF_TIMER
  USE SF_IOTOOLS,  only: free_unit,reg,str
  USE SF_ARRAYS,   only: arange,linspace
  USE SF_LINALG,   only: eigh
  USE SF_MISC,     only: assert_shape
  !
  implicit none
  private

  public :: nca_diagonalize_hlocal


contains





  !+-------------------------------------------------------------------+
  !                    FULL DIAGONALIZATION
  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+-------------------------------------------------------------------+
  subroutine nca_diagonalize_hlocal
    integer                     :: i,j,isector,dim,unit,nup,ndw
    real(8),dimension(Nsectors) :: e0 
    real(8)                     :: egs
    e0=0.d0
    write(LOGfile,"(A)")"Diagonalizing Hamiltonian:"
    call start_timer
    do isector=1,Nsectors
       nup  = getnup(isector)
       ndw  = getndw(isector)
       write(LOGfile,"(A,I4,A6,I2,A6,I2,A6,I15)")&
            "Solving sector:",isector,", nup:",nup,", ndw:",ndw,", dim=",getdim(isector)
       dim=getdim(isector)
       call build_hlocal(isector,espace(isector)%H(:,:))
       call eigh(espace(isector)%H,espace(isector)%e,'V','U')
       e0(isector)=minval(espace(isector)%e)
    enddo
    call stop_timer
    !
    egs=minval(e0)
    forall(isector=1:Nsectors)espace(isector)%e = espace(isector)%e - egs
    !Get the partition function Z and rescale energies
    zeta_function_diag=0d0
    do isector=1,Nsectors
       dim=getdim(isector)
       do i=1,dim
          zeta_function_diag=zeta_function_diag+exp(-beta*espace(isector)%e(i))
       enddo
    enddo
    write(LOGfile,"(A)")"DIAG summary:"
    write(LOGfile,"(A,f20.12)")'egs  =',egs
    write(LOGfile,"(A,f20.12)")'Z    =',zeta_function_diag
    unit=free_unit()
    open(unit,file="Egs.nca")
    write(unit,*)egs
    close(unit)
    return
  end subroutine nca_diagonalize_hlocal











  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix DOUBLE PRECISION
  !+------------------------------------------------------------------+
  subroutine build_hlocal(isector,Hmat)
    real(8),dimension(:,:)           :: Hmat
    integer                          :: isector
    integer,dimension(:),allocatable :: Hmap    !map of the Sector S to Hilbert space H
    integer,dimension(Nlevels)       :: ib
    integer                          :: dim
    integer                          :: i,j,m,iorb,jorb,ispin,impi,ishift
    integer                          :: k1,k2,k3,k4
    real(8)                          :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)          :: nup,ndw
    real(8)                          :: htmp
    real(8),dimension(Nspin,Norb)    :: eloc
    logical                          :: Jcondition
    type(sector_map)                 :: H
    !
    dim=getdim(isector)
    !
    call build_sector(isector,H)
    !
    call assert_shape(Hmat,[Dim,Dim],"Build_Hlocal","Hmat")
    !
    Hmat = 0d0
    !    
    !-----------------------------------------------!
    states: do i=1,Dim
       m = H%map(i)
       impi = i
       ib = bdecomp(m,2*Ns)
       !
       do iorb=1,Norb
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       htmp=0.d0
       !
       !LOCAL HAMILTONIAN PART:
       !local energies
       htmp = htmp - xmu*(sum(nup)+sum(ndw))
       !
       do iorb=1,Norb
          htmp = htmp + dreal(impHloc(1,1,iorb,iorb))*nup(iorb)
          htmp = htmp + dreal(impHloc(Nspin,Nspin,iorb,iorb))*ndw(iorb)
       enddo
       !
       Hmat(impi,i)=Hmat(impi,i)+htmp
       !
       !
       !Density-density interaction: same orbital, opposite spins
       !Euloc=\sum=i U_i*(n_u*n_d)_i
       htmp = zero
       do iorb=1,Norb
          htmp = htmp + Uloc(iorb)*nup(iorb)*ndw(iorb)
       enddo
       !
       if(Norb>1)then
          !density-density interaction: different orbitals, opposite spins:
          ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
          ! =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                htmp = htmp + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))
             enddo
          enddo
          !density-density interaction: different orbitals, parallel spins
          ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
          ! = \sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                htmp = htmp + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))
             enddo
          enddo
       endif
       !if using the Hartree-shifted chemical potential: mu=0 for half-filling
       !sum up the contributions of hartree terms:
       if(hfmode)then
          do iorb=1,Norb
             htmp = htmp - 0.5d0*Uloc(iorb)*(nup(iorb)+ndw(iorb)) + 0.25d0*uloc(iorb)
          enddo
          if(Norb>1)then
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   htmp=htmp-0.5d0*Ust*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))+0.25d0*Ust
                   htmp=htmp-0.5d0*(Ust-Jh)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))+0.25d0*(Ust-Jh)
                enddo
             enddo
          endif
       endif
       !
       Hmat(impi,i)=Hmat(impi,i)+htmp
       !
       !
       !
       !SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
       !S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
       !S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
       if(Norb>1.AND.Jhflag)then
          do iorb=1,Norb
             do jorb=1,Norb
                Jcondition=(&
                     (iorb/=jorb).AND.&
                     (ib(jorb)==1).AND.&
                     (ib(iorb+Ns)==1).AND.&
                     (ib(jorb+Ns)==0).AND.&
                     (ib(iorb)==0))
                if(Jcondition)then
                   call c(jorb,m,k1,sg1)
                   call c(iorb+Ns,k1,k2,sg2)
                   call cdg(jorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
                   j=binary_search(Hmap,k4)
                   htmp = Jx*sg1*sg2*sg3*sg4
                   !
                   if(j/=0)Hmat(impi,j)=Hmat(impi,j)+htmp
                   !
                endif
             enddo
          enddo
       endif
       !
       !PAIR-HOPPING (P-H) TERMS
       !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
       !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
       if(Norb>1.AND.Jhflag)then
          do iorb=1,Norb
             do jorb=1,Norb
                Jcondition=(&
                     (iorb/=jorb).AND.&
                     (ib(jorb)==1).AND.&
                     (ib(jorb+Ns)==1).AND.&
                     (ib(iorb+Ns)==0).AND.&
                     (ib(iorb)==0))
                if(Jcondition)then
                   call c(jorb,m,k1,sg1)
                   call c(jorb+Ns,k1,k2,sg2)
                   call cdg(iorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
                   j=binary_search(Hmap,k4)
                   htmp = Jp*sg1*sg2*sg3*sg4
                   !
                   if(j/=0)Hmat(impi,j)=Hmat(impi,j)+htmp
                   !
                endif
             enddo
          enddo
       endif
       !
       !
       !
       !IMPURITY LOCAL HYBRIDIZATION
       !Off-diagonal elements, i.e. non-local part
       !same spin:
       do iorb=1,Norb
          do jorb=1,Norb
             !this loop considers only the orbital off-diagonal terms
             !because iorb=jorb can not have simultaneously
             !occupation 0 and 1, as required by this if Jcondition:
             !UP
             Jcondition = &
                  (impHloc(1,1,iorb,jorb)/=zero) .AND. &
                  (ib(jorb)==1)                  .AND. &
                  (ib(iorb)==0)
             if (Jcondition) then
                call c(jorb,m,k1,sg1)
                call cdg(iorb,k1,k2,sg2)
                j = binary_search(H%map,k2)
                htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
                !
                Hmat(i,j) = Hmat(i,j) + htmp
                !
             endif
             !DW
             Jcondition = &
                  (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                  (ib(jorb+Ns)==1)                       .AND. &
                  (ib(iorb+Ns)==0)
             if (Jcondition) then
                call c(jorb+Ns,m,k1,sg1)
                call cdg(iorb+Ns,k1,k2,sg2)
                j = binary_search(H%map,k2)
                htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
                !
                Hmat(i,j) = Hmat(i,j) + htmp
                !
             endif
          enddo
       enddo
       !
    enddo states
    !-----------------------------------------------!
    !
    call delete_sector(isector,H)
    !
  end subroutine build_hlocal






end MODULE NCA_DIAG
