MODULE NCA_AUX_FUNX
  USE NCA_VARS_GLOBAL
  !
  USE SF_CONSTANTS, only:zero
  USE SF_IOTOOLS, only:free_unit,reg,str
  USE SF_MISC, only: assert_shape
  !
  implicit none
  private

  interface set_Hloc
     module procedure :: set_Hloc_nn
     module procedure :: set_Hloc_so
  end interface set_Hloc

  interface lso2nnn_reshape
     module procedure d_nlso2nnn
     module procedure c_nlso2nnn
  end interface lso2nnn_reshape

  interface so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
  end interface so2nn_reshape

  interface nnn2lso_reshape
     module procedure d_nnn2nlso
     module procedure c_nnn2nlso
  end interface nnn2lso_reshape

  interface nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
  end interface nn2so_reshape


  public :: set_Hloc
  public :: print_Hloc
  !
  public :: lso2nnn_reshape
  public :: so2nn_reshape
  public :: nnn2lso_reshape
  public :: nn2so_reshape
  !
  ! public :: search_chemical_potential 

contains

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
          write(unit,fmt)((str(Hloc(ispin,jspin,iorb,jorb),d=3)//" ",jorb=1,Norb),jspin=1,Nspin)
       enddo
    enddo
    write(unit,*)""
    if(present(file))close(unit)
  end subroutine print_Hloc








  !##################################################################
  !                   RESHAPE ROUTINES
  !##################################################################
  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlso][Nlso] shape
  ! from/to the [Nlat][Nspin][Nspin][Norb][Norb] shape.
  ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]  !
  ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn
  function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                            :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn

  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                  :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlat][Nspin][Nspin][Norb][Norb] shape
  ! from/to the [Nlso][Nlso] shape.
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
  !+-----------------------------------------------------------------------------+!
  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso

  function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                            :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso

  function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                            :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                            :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso







  ! !##################################################################
  ! !##################################################################
  ! ! ROUTINES TO SEARCH CHEMICAL POTENTIAL UP TO SOME ACCURACY
  ! ! can be used to fix any other *var so that  *ntmp == nread
  ! !##################################################################
  ! !##################################################################

  ! !+------------------------------------------------------------------+
  ! !PURPOSE  : 
  ! !+------------------------------------------------------------------+
  ! subroutine search_chemical_potential(var,ntmp,converged)
  !   real(8),intent(inout) :: var
  !   real(8),intent(in)    :: ntmp
  !   logical,intent(inout) :: converged
  !   logical               :: bool
  !   real(8)               :: ndiff
  !   integer,save          :: count=0,totcount=0,i
  !   integer,save          :: nindex=0
  !   integer               :: nindex_old(3)
  !   real(8)               :: ndelta_old,nratio
  !   integer,save          :: nth_magnitude=-2,nth_magnitude_old=-2
  !   real(8),save          :: nth=1.d-2
  !   logical,save          :: ireduce=.true.
  !   integer               :: unit
  !   !
  !   ndiff=ntmp-nread
  !   nratio = 0.5d0;!nratio = 1.d0/(6.d0/11.d0*pi)
  !   !
  !   !check actual value of the density *ntmp* with respect to goal value *nread*
  !   count=count+1
  !   totcount=totcount+1
  !   if(count>2)then
  !      do i=1,2
  !         nindex_old(i+1)=nindex_old(i)
  !      enddo
  !   endif
  !   nindex_old(1)=nindex
  !   !
  !   if(ndiff >= nth)then
  !      nindex=-1
  !   elseif(ndiff <= -nth)then
  !      nindex=1
  !   else
  !      nindex=0
  !   endif
  !   !
  !   ndelta_old=ndelta
  !   bool=nindex/=0.AND.( (nindex+nindex_old(1)==0).OR.(nindex+sum(nindex_old(:))==0) )
  !   !if(nindex_old(1)+nindex==0.AND.nindex/=0)then !avoid loop forth and back
  !   if(bool)then
  !      ndelta=ndelta_old*nratio !decreasing the step
  !   else
  !      ndelta=ndelta_old
  !   endif
  !   !
  !   if(ndelta_old<1.d-9)then
  !      ndelta_old=0.d0
  !      nindex=0
  !   endif
  !   !update chemical potential
  !   var=var+dble(nindex)*ndelta
  !   !xmu=xmu+dble(nindex)*ndelta
  !   !
  !   !Print information
  !   write(LOGfile,"(A,f16.9,A,f15.9)")"n    = ",ntmp," /",nread
  !   if(nindex>0)then
  !      write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," ==>"
  !   elseif(nindex<0)then
  !      write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," <=="
  !   else
  !      write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," == "
  !   endif
  !   write(LOGfile,"(A,f15.9)")"var  = ",var
  !   write(LOGfile,"(A,ES16.9,A,ES16.9)")"dn   = ",ndiff,"/",nth
  !   unit=free_unit()
  !   open(unit,file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
  !   write(unit,*)var,ntmp,ndiff
  !   close(unit)
  !   !
  !   !check convergence within actual threshold
  !   !if reduce is activetd
  !   !if density is in the actual threshold
  !   !if DMFT is converged
  !   !if threshold is larger than nerror (i.e. this is not last loop)
  !   bool=ireduce.AND.(abs(ndiff)<nth).AND.converged.AND.(nth>nerr)
  !   if(bool)then
  !      nth_magnitude_old=nth_magnitude        !save old threshold magnitude
  !      nth_magnitude=nth_magnitude_old-1      !decrease threshold magnitude || floor(log10(abs(ntmp-nread)))
  !      nth=max(nerr,10.d0**(nth_magnitude))   !set the new threshold 
  !      count=0                                !reset the counter
  !      converged=.false.                      !reset convergence
  !      ndelta=ndelta_old*nratio                  !reduce the delta step
  !      !
  !   endif
  !   !
  !   !if density is not converged set convergence to .false.
  !   if(abs(ntmp-nread)>nth)converged=.false.
  !   !
  !   !check convergence for this threshold
  !   !!---if smallest threshold-- NO MORE
  !   !if reduce is active (you reduced the treshold at least once)
  !   !if # iterations > max number
  !   !if not yet converged
  !   !set threshold back to the previous larger one.
  !   !bool=(nth==nerr).AND.ireduce.AND.(count>niter).AND.(.not.converged)
  !   bool=ireduce.AND.(count>niter).AND.(.not.converged)
  !   if(bool)then
  !      ireduce=.false.
  !      nth=10.d0**(nth_magnitude_old)
  !   endif
  !   !
  !   write(LOGfile,"(A,I5)")"count= ",count
  !   write(LOGfile,"(A,L2)")"Converged=",converged
  !   print*,""
  !   !
  ! end subroutine search_chemical_potential


END MODULE NCA_AUX_FUNX
