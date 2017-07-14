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
    integer                            :: i,j,isector,dim,unit,stride
    real(8),dimension(Nsectors)        :: e0 
    real(8)                            :: egs,zeta_function_diag
    real(8),dimension(:,:),allocatable :: Hmat
    real(8),dimension(:),allocatable   :: Eval
    !
    e0=0d0
    !
    write(LOGfile,"(A)")"Diagonalizing Hamiltonian:"
    call start_timer
    stride=0
    EigBasis=0d0
    EigValues=0d0
    do isector=1,Nsectors
       dim=getdim(isector)
       !
       allocate(Hmat(dim,dim));Hmat=0d0
       allocate(Eval(dim));Eval=0d0
       !
       call build_hlocal(isector,Hmat)
       call eigh(Hmat,Eval,'V','U')
       !
       e0(isector)=minval(Eval)
       !
       write(LOGfile,"(A,F15.8)")&
            "Solving sector: "//str(isector)//&
            ", nup:"//str(getnup(isector))//", ndw:"//str(getndw(isector))//&
            ", dim="//str(getdim(isector))//&
            ", Egs=",e0(isector)
       !
       EigBasis(stride+1:stride+dim,stride+1:stride+dim) = Hmat(:,:)
       EigValues(stride+1:stride+dim)                    = Eval
       !
       stride=stride+dim
       !
       deallocate(Hmat)
       deallocate(Eval)
    enddo
    call stop_timer
    !

    !Get the partition function Z and rescale energies
    Egs=minval(e0)
    EigValues = EigValues - Egs
    zeta_function_diag=sum(exp(-beta*EigValues(:)))
    !
    write(LOGfile,"(A)")"DIAG summary:"
    write(LOGfile,"(A,f20.12)")'Egs  =',Egs
    write(LOGfile,"(A,f20.12)")'Z    =',zeta_function_diag
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
       Hmat(i,i)=Hmat(i,i)+htmp
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
       Hmat(i,i)=Hmat(i,i)+htmp
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
                   if(j/=0)Hmat(i,j)=Hmat(i,j)+htmp
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
                   if(j/=0)Hmat(i,j)=Hmat(i,j)+htmp
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
