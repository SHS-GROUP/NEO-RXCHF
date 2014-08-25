!=======================================================================
      program xcneo_mpi

!=======================================================================
      implicit none
      include "mpif.h"

#if MPI32
      integer*4 ierr,req,prov
#else
      integer ierr,req,prov
#endif
      parameter (req=MPI_THREAD_FUNNELED)

      call MPI_INIT_THREAD(req,prov,ierr)

      call xcneo_mpi_driver

      call MPI_FINALIZE(ierr)

      end
!=======================================================================
      subroutine xcneo_mpi_driver

!=======================================================================
      implicit none
      include 'mpif.h'

!-------Basis Set Info-------(
      integer,allocatable :: ELCAM(:,:)  ! Angular mom for electrons
      integer,allocatable :: NUCAM(:,:)   ! Angular mom for quantum nuclei
      double precision,allocatable :: ELCEX(:) ! Exponents: elec basis
      double precision,allocatable :: NUCEX(:)  ! Exponents: nuc basis
      double precision,allocatable :: ELCBFC(:,:) ! Basis centers: elec basis
      double precision,allocatable :: NUCBFC(:,:)  ! basis centers: nuc basis
      integer,allocatable :: AMPEB2C(:) ! Map primitive index to contracted
      double precision,allocatable :: AGEBFCC(:) ! Map prim index to contract coef
      double precision,allocatable :: AGNBFCC(:)  ! Nuclear contract coef
      integer,allocatable :: KPESTR(:)  ! Map contracted index to primitive start
      integer,allocatable :: KPEEND(:)  ! Map contracted index to primitive end
      double precision,allocatable :: zan(:) ! Classical nuclear charges
      double precision,allocatable :: cat(:,:) ! XYZ Coordinates of classical atoms
      double precision,allocatable :: bgem(:)
      double precision,allocatable :: ggem(:)
      integer nat
      integer ngtg
      integer nelec
      integer NAE               ! Number of regular electrons
      integer NBE               ! Number of special electrons
      integer NAE_alp
      integer NAE_beta
      integer NAalpE,NAbetE
      double precision pmass    ! Mass of nonelectron quantum particle 
!-------Basis Set Info-------)
      double precision,allocatable :: GM2ICR(:)
      double precision,allocatable :: GM2sICR(:)
      double precision,allocatable :: GM2exICR(:)
      double precision,allocatable :: GM3IC1(:)
      double precision,allocatable :: GM4ICR(:)
      double precision,allocatable :: GM2_1ICR(:)
      double precision,allocatable :: GM2_2ICR(:)
      double precision,allocatable :: GM3_1IC1(:)
      double precision,allocatable :: GM3_2IC1(:)
      integer SZG2ICR
      integer SZG3IC1
      integer SZG4IC
      integer i,j,idum,istat
      integer NUCST
      integer NG2CHK
      integer NG3CHK
      integer NG4CHK
      integer EXCHLEV
      logical LNEOHF
      logical LXCUHF
      logical LXCROHF
      logical LG2DSCF
      logical LG3DSCF
      logical LG2IC1 ! GAM2 Integrals In-Core algo-1: contracted
      logical LG3IC1 ! GAM3 Integrals In-Core algo-1
      logical LG3IC2 ! GAM3 Integrals In-Core algo-2
      logical LG4DSCF
      logical LG4IC  ! GAM4 Integrals In-Core 
      logical read_CE
      logical read_CP
      logical read_GAM2
      logical read_GAM3
      logical read_GAM4
      logical LGAM4
      logical LDBG
      logical LSOSCF
      logical LDIIS
      logical LSORXCHF
      integer OCBSETYP
      logical LRXCHF
      logical LRXCUHF

      double precision a2bohr,bohr2a
      parameter(bohr2a=0.529177249d+00)
      parameter(a2bohr=1.0d+00/0.529177249d+00)

      integer npebf,nebf,npbf

      integer ngee
      integer ng1,ng2,ng3,ng4
      integer ng1prm,ng2prm,ng3prm,ng4prm

      integer junk

#if MPI32
      integer*4 nproc4,rank4,namelen,ierr
      integer*4 ngtg_,nat_,npebf_,npbf_,nebf_
#else
      integer namelen,ierr
      integer ngtg_,nat_,npebf_,npbf_,nebf_
#endif
      integer nproc,rank
      character(len=25) :: procname

      double precision wtime,wtime1,wtime2

      integer nblocks
      integer blockrank

      logical foundrec,inperr
      character*70 exe
      character*30 fname
      character*72 line
      character*6  CATOMS
      character*7  CEBASIS
      character*7  CPBASIS
      character*1  istring

      namelist /sysinfo/ nat,nebf,npebf,npbf,ngtg,pmass,nelec,nucst
      namelist /control/ lneohf,lxcuhf,lxcrohf,lrxchf,lrxcuhf
      namelist /guessmo/ read_ce,read_cp
      namelist /intctrl/ read_gam2,read_gam3,read_gam4,
     x                   ng2chk,ng3chk,ng4chk,nblocks,blockrank
      namelist /scfctrl/ ldbg,lsoscf,ldiis,lsorxchf,ocbsetyp
      namelist /geminal/ bgem,ggem
      namelist /xcuhf  / nae,nbe
      namelist /rxchf  / nae,nbe,exchlev
      namelist /rxcuhf / nae_alp,nae_beta

#if MPI32
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc4,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank4,ierr)
      rank=int(rank4,kind=8)
      nproc=int(nproc4,kind=8)
#else
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
#endif
      call MPI_GET_PROCESSOR_NAME(procname,namelen,ierr)

      if(rank.eq.0) then
       write(*,2000)

       call getarg(0,exe)
       call getarg(1,fname)
       write(*,*) "Running executable: "
       write(*,*) "   "//trim(exe)
       write(*,*) "Input file: "
       write(*,*) "   "//trim(fname)
       write(*,*)
#if MPI32
       write(*,*) "Assuming 32-bit integer MPI library"
#else
       write(*,*) "Assuming 64-bit integer MPI library"
#endif
       write(*,*)
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      write(*,5000) "Process with rank",rank,
     x              "running on host",trim(procname)

! Initialize
      LNEOHF=.false.
      LXCUHF=.false.
      LXCROHF=.false.
      LRXCHF=.true.
      LRXCUHF=.false.
      read_CE=.false.
      read_CP=.false.
      read_GAM2=.false.
      read_GAM3=.false.
      read_GAM4=.false.
      NG2CHK=1
      NG3CHK=1
      NG4CHK=1
      LGAM4=.true.
      LG4DSCF=.false.
      LG4IC=.true.
      LG3DSCF=.false.
      LG3IC1=.true.
      LG3IC2=.false.
      LG2DSCF=.false.
      LG2IC1=.true.
      LDBG=.false.
      LSOSCF=.false.
      LDIIS=.false.
      LSORXCHF=.false.
      OCBSETYP=2
      CATOMS='&ATOMS'
      CEBASIS='&EBASIS'
      CPBASIS='&PBASIS'
      nblocks=1
      blockrank=0
      inperr=.false.

! Read from file with master process
      if(rank.eq.0) then

       open(unit=9,file=trim(fname)//'.inp')

!!!!!!!!!!!!!!! Read system info !!!!!!!!!!!!!! 
       read(9,nml=sysinfo)
       write(*,*)
       write(*,nml=sysinfo)
       write(*,*)

!!!!!!!!!!!!!!! Read control info !!!!!!!!!!!!!! 
       rewind(9)
       read(9,nml=control)
       write(*,*)
       write(*,nml=control)
       write(*,*)

       if(LNEOHF.and.LRXCHF) then
        LRXCHF=.false.
        write(*,*) "Overriding LRXCHF since LNEOHF=.TRUE."
       end if
       if(LNEOHF.and.LRXCUHF) then
        LRXCUHF=.false.
        LXCUHF=.false.
        write(*,*) "Overriding LRXCUHF since LNEOHF=.TRUE."
       end if

       if(LRXCHF.and.LRXCUHF) then
        LRXCHF=.false.
        write(*,*) "Overriding LRXCHF since LRXCUHF=.TRUE."
       end if

       if ((LRXCHF).and.(LXCUHF)) then
        write(*,*) "Cannot have LRXCHF and LXCUHF"
        write(*,*) "Exiting..."
        inperr=.true.
       end if

       if ((LRXCUHF).and.(LXCUHF)) then
        write(*,*) "Cannot have LRXCUHF and LXCUHF"
        write(*,*) "Exiting..."
        inperr=.true.
       end if

!!!!!!!!!!!!!!! Read guess info !!!!!!!!!!!!!! 
       rewind(9)
       read(9,nml=guessmo)
       write(*,*)
       write(*,nml=guessmo)
       write(*,*)

!!!!!!!!!!!!!!! Read integral info !!!!!!!!!!!!!! 
       rewind(9)
       read(9,nml=intctrl)
       write(*,*)
       write(*,nml=intctrl)
       write(*,*)
 
!!!!!!!!!!!!!!! Read SCF info !!!!!!!!!!!!!! 
       rewind(9)
       read(9,nml=scfctrl)
       write(*,*)
       write(*,nml=scfctrl)
       write(*,*)
       if(.not.LSORXCHF) then
        write(*,*) "Must use SORXCHF to use UHF for regulars"
        write(*,*) "Exiting...."
        call exit
       end if       
       if(LSOSCF.and.LDIIS) then
        write(*,*) "Cannot have both SOSCF and DIIS"
        write(*,*) "Proceeding with DIIS..."
        LSOSCF=.false.
       end if
       if(LSOSCF.and.LSORXCHF) then
        write(*,*) "Cannot have both SOSCF and SORXCHF"
        write(*,*) "Proceeding with SORXCHF..."
        LSOSCF=.false.
       end if
C       if(LDIIS.and.LSORXCHF) then
C        write(*,*) "Cannot have both DIIS and SORXCHF"
C        write(*,*) "Proceeding with SORXCHF..."
C        LDIIS=.false.
C       end if

!!!!!!!!!!!!!!! Read in geminal info !!!!!!!!!!!!!! 
       if(allocated(bgem)) deallocate(bgem)
       if(allocated(ggem)) deallocate(ggem)

       if(.not.(LNEOHF)) then
        allocate(bgem(ngtg))
        allocate(ggem(ngtg))
        rewind(9)
        read(9,nml=geminal)
       else
        ngtg=1
        allocate(bgem(ngtg))
        allocate(ggem(ngtg))
        bgem=1.0d+00
        ggem=1.0d-15
       end if

!!!!!!!!!!!!!!! Read in XCHF/NEO-HF openshell info !!!!!!!!!!!!!!!
       if(LXCUHF) then
        rewind(9)
        read(9,nml=xcuhf)
        write(*,*)
        write(*,nml=xcuhf)
        write(*,*)
        if((nae+nbe).ne.nelec) then
         write(*,*) "Num alpha electrons:",NAE
         write(*,*) "Num beta  electrons:",NBE
         write(*,*) "Num total electrons:",nelec
         write(*,*) "Not equal. Exiting..."
         inperr=.true.
        end if
       end if

!!!!!!!!!!!!!!! Read in RXCHF info !!!!!!!!!!!!!!!
       if((LRXCHF).or.(LRXCUHF)) then
        rewind(9)
        read(9,nml=rxchf)
        write(*,*)
        write(*,nml=rxchf)
        write(*,*)
        if((nae+nbe).ne.nelec) then
         write(*,*) "Num regular electrons:",NAE
         write(*,*) "Num special electrons:",NBE
         write(*,*) "Num total   electrons:",nelec
         write(*,*) "Not equal. Exiting..."
         inperr=.true.
        end if
        if(.not.(LRXCUHF)) then
         nae_alp=nae/2
         nae_beta=nae/2
        end if
       end if

       if (((LRXCHF).or.(LRXCUHF)).and.(nelec.ge.4).and.(exchlev.eq.2))
     x       then
        write(*,*) "Currently only nelec<4 is supported for RXCHF-fe"
        write(*,*) "Exiting..."
        inperr=.true.
       end if

       if (((LRXCHF).or.(LRXCUHF)).and.
     x          ((exchlev.lt.0).or.(exchlev.gt.2))) then
        write(*,*) "EXCHLEV must be between 0 and 2:"
        write(*,*) "        = 0 for RXCHF-ne"
        write(*,*) "        = 1 for RXCHF-ae"
        write(*,*) "        = 2 for RXCHF-fe"
        write(*,*) "Exiting..."
        inperr=.true.
       end if

       if ((EXCHLEV.gt.0).and.(LRXCHF).and.(NBE.eq.1)) then
        write(*,*) "EXCHLEV > 0 should only be used with LRXCUHF"
        write(*,*) "Exiting..."
        inperr=.true.
       end if

       if((LRXCUHF).and.(nae.eq.1)) then
        write(*,*) "For NAE=1 regular electron, use RXCHF-ne"
        write(*,*) "Exiting..."
        inperr=.true.
       end if

       if ((LRXCHF.or.LRXCUHF).and.(NBE.gt.1).and.(EXCHLEV.gt.1)) then
        write(*,*) "Only RXCHF-ne and RXCHF-ae are currently"
        write(*,*) "supported for more than one special electron."
        write(*,*) "Exiting..."
        inperr=.true.
       end if

       if ((LRXCHF.or.LRXCUHF).and.(NBE.gt.2).and.(EXCHLEV.eq.1)) then
        write(*,*) "WARNING: A balanced version of RXCHF-ae has only"
        write(*,*) "         been derived for up to two special"
        write(*,*) "         electrons: use with caution for NBE>2!"
       end if

       if((LRXCHF.or.LRXCUHF).and.((OCBSETYP.lt.0).or.(OCBSETYP.gt.3)))
     x  then
        write(*,*) "OCBSETYP must be one of:"
        write(*,*) " = 0 : OCBSE off"
        write(*,*) " = 1 : OCBSE"
        write(*,*) " = 2 : OCBSE2"
        write(*,*) " = 3 : OCBSE3"
        write(*,*) "Exiting..."
        inperr=.true.
       end if

!!!!!!!!!!!!!!! Read in RXCUHF info !!!!!!!!!!!!!!!
       if(LRXCUHF) then
        rewind(9)
        read(9,nml=rxcuhf)
        write(*,*)
        write(*,nml=rxcuhf)
        write(*,*)
        if((nae_beta).gt.nae_alp) then
         write(*,*) "Num regular alpha electrons:",NAE_alp
         write(*,*) "Num regular beta  electrons:",NAE_beta
         write(*,*) "Num of alpha electrons must be equal or greater"
         write(*,*) "than number of beta electrons. Exiting..."
         inperr=.true.
        end if
        if((nae_alp+nae_beta).ne.nae) then
         write(*,*) "Num regular alpha electrons:",NAE_alp
         write(*,*) "Num regular beta  electrons:",NAE_beta
         write(*,*) "Num regular       electrons:",NAE
         write(*,*) "Not equal. Exiting..."
         inperr=.true.
        end if
       end if

!!!!!!!!!!!!!!! Read in atomic coordinates !!!!!!!!!!!!!!!
       if(allocated(zan)) deallocate(zan)
       if(allocated(cat)) deallocate(cat)
       allocate(zan(nat))
       allocate(cat(3,nat))

       rewind(9)
       foundrec=.false.
       do while (.not.(foundrec))
         read(9,*) line
         if(trim(line).eq.CATOMS) then
          foundrec=.true.
         end if
       end do

       do i=1,nat
          read(9,*)zan(i),cat(1,i),cat(2,i),cat(3,i)
          cat(1,i)=a2bohr*cat(1,i)
          cat(2,i)=a2bohr*cat(2,i)
          cat(3,i)=a2bohr*cat(3,i)
       end do

!!!!!!!!!!!!!!! Read in electronic basis sets !!!!!!!!!!!!!!!
       if(allocated(AMPEB2C)) deallocate(AMPEB2C)
       if(allocated(ELCEX)) deallocate(ELCEX)
       if(allocated(AGEBFCC)) deallocate(AGEBFCC)
       if(allocated(ELCAM)) deallocate(ELCAM)
       if(allocated(ELCBFC)) deallocate(ELCBFC)
       allocate(AMPEB2C(npebf))
       allocate(ELCEX(npebf))
       allocate(AGEBFCC(npebf))
       allocate(ELCAM(npebf,3))
       allocate(ELCBFC(npebf,3))

       rewind(9)
       foundrec=.false.
       do while (.not.(foundrec))
         read(9,*) line
         if(trim(line).eq.CEBASIS) then
          foundrec=.true.
         end if
       end do

       do i=1,npebf
          read(9,*)idum,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),ELCAM(i,3),
     x  ELCEX(i),AGEBFCC(i),ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
          do j=1,3
             ELCBFC(i,j)=a2bohr*ELCBFC(i,j)
          end do
       end do

!!!!!!!!!!!!!!! Read in nuclear basis sets !!!!!!!!!!!!!!!
       if(allocated(NUCEX)) deallocate(NUCEX)
       if(allocated(AGNBFCC)) deallocate(AGNBFCC)
       if(allocated(NUCAM)) deallocate(NUCAM)
       if(allocated(NUCBFC)) deallocate(NUCBFC)
       allocate(NUCEX(npbf))
       allocate(AGNBFCC(npbf))
       allocate(NUCAM(npbf,3))
       allocate(NUCBFC(npbf,3))

       rewind(9)
       foundrec=.false.
       do while (.not.(foundrec))
         read(9,*) line
         if(trim(line).eq.CPBASIS) then
          foundrec=.true.
         end if
       end do

       do i=1,npbf
          read(9,*)idum,idum,NUCAM(i,1),NUCAM(i,2),NUCAM(i,3),
     x  NUCEX(i),AGNBFCC(i),NUCBFC(i,1),NUCBFC(i,2),NUCBFC(i,3)
          do j=1,3
             NUCBFC(i,j)=a2bohr*NUCBFC(i,j)
          end do
       end do

       close(9)

       if(allocated(KPESTR)) deallocate(KPESTR)
       allocate( KPESTR(nebf),stat=istat )
       if(allocated(KPEEND)) deallocate(KPEEND)
       allocate( KPEEND(nebf),stat=istat )
       call make_KPE(nebf,npebf,AMPEB2C,KPESTR,KPEEND)

      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call MPI_BCAST(inperr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      if(inperr) return

! Broadcast static variables including array dimensions to slave processes
      call MPI_BCAST(ngtg,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nebf,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(npebf,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(npbf,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nat,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)

! INTEGER8 and LOGICAL must be declared specificially since MPI_INTEGER
! and MPI_LOGICAL do not seem to get specified as long with compiler flags 
      call MPI_BCAST(pmass,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nelec,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NAE,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NBE,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NAE_alp,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NAE_beta,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NUCST,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LNEOHF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LXCUHF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LXCROHF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LRXCHF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LRXCUHF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(read_CE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(read_CP,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(read_GAM2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(read_GAM3,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(read_GAM4,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NG4CHK,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LGAM4,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG4DSCF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG4IC,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NG3CHK,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG3DSCF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG3IC1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG3IC2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NG2CHK,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG2DSCF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG2IC1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LDBG,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LSOSCF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LDIIS,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LSORXCHF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(OCBSETYP,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(EXCHLEV,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nblocks,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(blockrank,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Allocate arrays on slave processes
      if (rank.gt.0) then
       if(allocated(bgem)) deallocate(bgem)
       allocate( bgem(ngtg),stat=istat )
       if(allocated(ggem)) deallocate(ggem)
       allocate( ggem(ngtg),stat=istat )
       if(allocated(zan)) deallocate(zan)
       allocate( zan(nat),stat=istat )
       if(allocated(cat)) deallocate(cat)
       allocate( cat(3,nat),stat=istat )
       if(allocated(AMPEB2C)) deallocate(AMPEB2C)
       allocate( AMPEB2C(npebf),stat=istat )
       if(allocated(ELCEX)) deallocate(ELCEX)
       allocate( ELCEX(npebf),stat=istat )
       if(allocated(AGEBFCC)) deallocate(AGEBFCC)
       allocate( AGEBFCC(npebf),stat=istat )
       if(allocated(ELCAM)) deallocate(ELCAM)
       allocate( ELCAM(npebf,3),stat=istat )
       if(allocated(ELCBFC)) deallocate(ELCBFC)
       allocate( ELCBFC(npebf,3),stat=istat )
       if(allocated(NUCEX)) deallocate(NUCEX)
       allocate( NUCEX(npbf),stat=istat )
       if(allocated(AGNBFCC)) deallocate(AGNBFCC)
       allocate( AGNBFCC(npbf),stat=istat )
       if(allocated(NUCAM)) deallocate(NUCAM)
       allocate( NUCAM(npbf,3),stat=istat )
       if(allocated(NUCBFC)) deallocate(NUCBFC)
       allocate( NUCBFC(npbf,3),stat=istat )
       if(allocated(KPESTR)) deallocate(KPESTR)
       allocate( KPESTR(nebf),stat=istat )
       if(allocated(KPEEND)) deallocate(KPEEND)
       allocate( KPEEND(nebf),stat=istat )
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Broadcast array data
#if MPI32
      ngtg_=int(ngtg,kind=4)
      nat_=int(nat,kind=4)
      npebf_=int(npebf,kind=4)
      npbf_=int(npbf,kind=4)
      nebf_=int(nebf,kind=4)
#else
      ngtg_=ngtg
      nat_=nat
      npebf_=npebf
      npbf_=npbf
      nebf_=nebf
#endif
      call MPI_BCAST(bgem,ngtg_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ggem,ngtg_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(zan,nat_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cat,3*nat_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(AMPEB2C,npebf_,MPI_INTEGER8,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ELCEX,npebf_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(AGEBFCC,npebf_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ELCAM,3*npebf_,MPI_INTEGER8,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ELCBFC,3*npebf_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NUCEX,npbf_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(AGNBFCC,npbf_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NUCAM,3*npbf_,MPI_INTEGER8,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NUCBFC,3*npbf_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(KPESTR,nebf_,MPI_INTEGER8,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(KPEEND,nebf_,MPI_INTEGER8,
     x               0,MPI_COMM_WORLD,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(LNEOHF.or.LXCUHF) then
       if (rank.eq.0) 
     x    write(*,*) "MPI code only supports RXCHF-ne/ae with NBE>1"
       return
      end if

      if(.not.((LRXCHF.or.LRXCUHF).and.(NBE.ge.2).and.(EXCHLEV.lt.2)))
     x  then
       if (rank.eq.0) 
     x    write(*,*) "MPI code only supports RXCHF-ne/ae with NBE>1"
       return
      end if

      ngee=nebf*nebf*nebf*nebf

      ng1=nebf*nebf*npbf*npbf
      ng2=nebf*nebf*nebf*nebf*npbf*npbf
      ng3=nebf*nebf*nebf*nebf*nebf*nebf*npbf*npbf
      ng4=nebf*nebf*nebf*nebf*nebf*nebf*nebf*nebf*npbf*npbf

      ng1prm=npebf*npebf*npbf*npbf
      ng2prm=npebf*npebf*npebf*npebf*npbf*npbf
      ng3prm=npebf*npebf*npebf*npebf*npebf*npebf*npbf*npbf

      if (rank.eq.0) then

       write(*,*)
       write(*,*) "========================================"
       write(*,*)
       if (LRXCHF) then
        write(*,*) " Running closed shell RXCHF calculation"
       else
        write(*,*) " Running open shell RXCHF calculation"
       end if
       if (EXCHLEV.eq.2) then
        write(*,*) " performed at the RXCHF-fe level"
       else if (EXCHLEV.eq.1) then
        write(*,*) " performed at the RXCHF-ae level"
       else
        write(*,*) " performed at the RXCHF-ne level"
       end if
       write(*,*)
       write(*,*) "========================================"
       write(*,*)

       write(*,*)
       write(*,*) "GEMINAL PARAMETERS"
       write(*,*) "  k          b_k                 gamma_k"
       do i=1,ngtg
         write(*,6000) i,bgem(i),ggem(i)
       end do

       write(*,*) 
       write(*,*) 'ATOMIC COORDINATES (BOHR)'
       write(*,*) 
       write(*,*) 'CHARGE -X COORDINATE- -Y COORDINATE- -Z COORDINATE-'
       do i=1,nat
         write(*,7000) zan(i),cat(1,i),cat(2,i),cat(3,i)
       end do
       write(*,*) 
      
       write(*,*)
       write(*,*)' CHECK CONTRACTED ELECTRONIC BASIS FUNCTIONS '
       write(*,*)'CONT INDEX    KPESTR     KPEEND'
       do i=1,nebf
          write(*,8000) i,KPESTR(i),KPEEND(i)
       end do
 
       WRITE(*,*)
       WRITE(*,*)'ELECTRONIC BASIS FUNCTIONS:'
       WRITE(*,*)
       WRITE(*,*)'PRIM  CONT    ANG      EXPONENT CONTRACT  -X- -Y- -Z-'
       WRITE(*,*)'INDEX INDEX   MOM                  COEF'
       DO i=1,npebf
         WRITE(*,9000) i,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),ELCAM(i,3),
     x  ELCEX(i),AGEBFCC(i),ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
       END DO

      end if

! Normalize the contraction coefficients for elec basis functions
      call ELCNORM3(npebf,nebf,
     x              AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

      if (rank.eq.0) then
       WRITE(*,*)
       WRITE(*,*)'ELECTRONIC BASIS FUNCTIONS:'
       WRITE(*,*)'CONTRACT COEFF HAVE BEEN NORMALIZED'
       WRITE(*,*)
       WRITE(*,*)'PRIM  CONT    ANG      EXPONENT CONTRACT  -X- -Y- -Z-'
       WRITE(*,*)'INDEX INDEX   MOM                  COEF'
       DO i=1,npebf
         WRITE(*,9000) i,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),ELCAM(i,3),
     x  ELCEX(i),AGEBFCC(i),ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
       END DO

       WRITE(*,*)
       WRITE(*,*)'NUCLEAR BASIS FUNCTIONS:'
       WRITE(*,*)
       WRITE(*,*)'PRIM  CONT    ANG      EXPONENT CONTRACT  -X- -Y- -Z-'
       WRITE(*,*)'INDEX INDEX   MOM                  COEF'
       DO i=1,npbf
         WRITE(*,9000) i,i,NUCAM(i,1),NUCAM(i,2),NUCAM(i,3),
     x  NUCEX(i),AGNBFCC(i),NUCBFC(i,1),NUCBFC(i,2),NUCBFC(i,3)
       END DO
      end if

! Normalize the contraction coefficients for nuc basis functions
      call NUCNORM(npbf,AGNBFCC,NUCEX,NUCAM,NUCBFC)

      if (rank.eq.0) then
       WRITE(*,*)
       WRITE(*,*)'NUCLEAR BASIS FUNCTIONS:'
       WRITE(*,*)'CONTRACT COEFF HAVE BEEN NORMALIZED'
       WRITE(*,*)
       WRITE(*,*)'PRIM  CONT    ANG      EXPONENT CONTRACT  -X- -Y- -Z-'
       WRITE(*,*)'INDEX INDEX   MOM                  COEF'
       DO i=1,npbf
         WRITE(*,9000) i,i,NUCAM(i,1),NUCAM(i,2),NUCAM(i,3),
     x  NUCEX(i),AGNBFCC(i),NUCBFC(i,1),NUCBFC(i,2),NUCBFC(i,3)
       END DO
      end if

      ng1prm=(npebf**2)*npbf**2
      ng2prm=(npebf**4)*npbf**2
      ng3prm=(npebf**6)*npbf**2

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Calculate inexpensive integrals and store on disk using only master process
!   - requires use of global scratch directories
      if (rank.eq.0) then

       call class_nuc_rep(nat,zan,cat)

       call elec_ovlap(npebf,nebf,nebf*nebf,
     x                 AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

       call check_elec_ovlap(nebf)

       call nuc_ovlap(npbf,npbf*npbf,AGNBFCC,NUCEX,NUCAM,NUCBFC)

       call check_nuc_ovlap(npbf)

       write(*,*)
       write(*,*)'**************************************'
       write(*,*)'    Computing GAM_epcore Integrals    '

       call calc_GAM_epcore(nebf,npebf,npbf,nebf*nebf,npbf*npbf,
     x                      nat,pmass,zan,cat,
     x                      AMPEB2C,AGEBFCC,AGNBFCC,
     x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

       write(*,*)'    Computing GAM_ep     Integrals    '

       call calc_GAM_ep(nebf,npebf,npbf,ng1,
     x                  AMPEB2C,AGEBFCC,AGNBFCC,
     x                  ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

       write(*,*)'    Computing GAM_ee     Integrals    '

       call calc_GAM_ee(nebf,npebf,ngee,
     x                  AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

       write(*,*)'**************************************'
       write(*,*)

      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call RXCHF_driver_MPI(nproc,rank,
     x                      nblocks,blockrank,
     x                      nelec,nae,nbe,nae_alp,nae_beta,nucst,
     x                      nebf,npebf,npbf,nat,ngtg,
     x                      ng1,ng2,ng3,ng4,ngee,
     x                      ng1prm,ng2prm,ng3prm,
     x                      pmass,cat,zan,bgem,ggem,
     x                      KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                      NG2CHK,NG3CHK,NG4CHK,
     x                      read_CE,read_CP,
     x                      read_GAM2,read_GAM3,read_GAM4,
     x                      LG2IC1,LG3IC1,LG4IC,
     x                      LG2DSCF,LG3DSCF,LG4DSCF,
     x                      LRXCUHF,
     x                      LSOSCF,LDIIS,LSORXCHF,OCBSETYP,
     x                      LDBG,EXCHLEV)

      write(*,*) "Process is finished:",rank

 5000 format(1X,A,1X,I4,1X,A,1X,A)
 6000 format(1X,I3,G20.8,G20.8)
 7000 format(1X,F5.1,3(2X,F13.8))
 8000 format(1X,I3,I6,I5)
 9000 format(1X,I3,I6,I5,I3,I3,F12.6,F10.6,F10.6,F10.6,F10.6)

 2000 format(/,
     x' ============================================================',/,
     x'        ______     _ _   _ _____ ___            ______       ',/,
     x'       / / / /    |   \ | | ____/ _ \           \ \ \ \      ',/,
     x'      / / / /     | |  \| |  _|| | | |_____      \ \ \ \     ',/,
     x'      \ \ \ \     | | |\  | |__| |_|  _____|     / / / /     ',/,
     x'       \_\_\_\    |  _| \_|_____\___/           /_/_/_/      ',/,
     x'                  |_|                                        ',/,
     x'                                                             ',/,
     x' ============================================================'/)

      return
      end

