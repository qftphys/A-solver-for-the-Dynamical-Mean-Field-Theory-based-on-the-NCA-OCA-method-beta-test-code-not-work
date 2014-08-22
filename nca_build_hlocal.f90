  states: do i=first_state,last_state
     m=Hmap(i)
     impi = i-ishift
     call bdecomp(m,ib)
     htmp=0.d0
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo
     !
     !LOCAL HAMILTONIAN PART:
     !local energies
     htmp = htmp - xmu*(sum(nup)+sum(ndw))
     htmp = htmp + dot_product(eloc(1,:),nup) + dot_product(eloc(Nspin,:),ndw)
     !
     !Density-density interaction: same orbital, opposite spins
     !Euloc=\sum=i U_i*(n_u*n_d)_i
     htmp = htmp + dot_product(uloc,nup*ndw)!=\sum=i U_i*(n_u*n_d)_i
     if(Norb>1)then
        !density-density interaction: different orbitals, opposite spins
        !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
        !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))
           enddo
        enddo
        !density-density interaction: different orbitals, parallel spins
        !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
        !    "="\sum_ij (Ust-Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
        !    "="\sum_ij (Uloc-3*Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))
           enddo
        enddo
     endif
     !if using the Hartree-shifted chemical potential: mu=0 for half-filling
     !sum up the contributions of hartree terms:
     if(hfmode)then
        htmp=htmp - 0.5d0*dot_product(uloc,nup+ndw) + 0.25d0*sum(uloc)
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
                 htmp = Jh*sg1*sg2*sg3*sg4
                 !
                 Hmat(impi,j)=Hmat(impi,j)+htmp
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
                 htmp = Jh*sg1*sg2*sg3*sg4
                 !
                 Hmat(impi,j)=Hmat(impi,j)+htmp
                 !
              endif
           enddo
        enddo
     endif
     !
     !
     !
     !IMPURITY LOCAL HYBRIDIZATION
     do iorb=1,Norb
        do jorb=1,Norb
           !SPIN UP
           if((ib(iorb)==0).AND.(ib(jorb)==1))then
              call c(jorb,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j=binary_search(Hmap,k2)
              htmp = Hloc(1,1,iorb,jorb)*sg1*sg2
              !
              Hmat(impi,j)=Hmat(impi,j)+htmp
              !
           endif
           !SPIN DW
           if((ib(iorb+Ns)==0).AND.(ib(jorb+Ns)==1))then
              call c(jorb+Ns,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j=binary_search(Hmap,k2)
              htmp = Hloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
              !
              Hmat(impi,j)=Hmat(impi,j)+htmp
              !
           endif
        enddo
     enddo
     !
  enddo states

