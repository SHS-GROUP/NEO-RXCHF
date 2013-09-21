!=======================================================================
      program xcneo_mpi

!=======================================================================
      implicit none
      include "mpif.h"

      integer ierr,req,prov
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
      double precision,allocatable :: bcoef1(:)
      double precision,allocatable :: gamma1(:)
      integer nat
      integer ngtg1
      integer nelec
      integer NAE               ! Number of regular electrons
      integer NBE               ! Number of special electrons
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
      logical LCMF
      logical LSOSCF
      logical LRXCHF
      logical LRXCUHF
      logical LOCBSE

      double precision a2bohr,bohr2a
      parameter(bohr2a=0.529177249d+00)
      parameter(a2bohr=1.0d+00/0.529177249d+00)

      integer npebf,nebf,npbf
      integer nebf2,npbf2,NPR,NEBFLT
      integer NPRA,NPRB

      integer ngee
      integer ng1,ng2,ng3,ng4
      integer ng1prm,ng2prm,ng3prm,ng4prm

      integer junk

      integer*4 nproc4,rank4,namelen,ierr
      integer nproc,rank
      character(len=25) :: procname

      double precision wtime,wtime1,wtime2

      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc4,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank4,ierr)
      call MPI_GET_PROCESSOR_NAME(procname,namelen,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      rank=int(rank4,kind=8)
      nproc=int(nproc4,kind=8)
      write(*,7000) "Process with rank",rank,
     x              "running on host",trim(procname)

! Read from file with master process
      if(rank.eq.0) then
         open(unit=9,file='basis_definition.inp')

         read(9,*) ngtg1
         if(allocated(bcoef1)) deallocate(bcoef1)
         allocate( bcoef1(ngtg1),stat=istat )
         if(allocated(gamma1)) deallocate(gamma1)
         allocate( gamma1(ngtg1),stat=istat )
         do i=1,ngtg1
            read(9,*)bcoef1(i),gamma1(i)
         end do

         read(9,*)nat
         if(allocated(zan)) deallocate(zan)
         allocate( zan(nat),stat=istat )
         if(allocated(cat)) deallocate(cat)
         allocate( cat(3,nat),stat=istat )
         do i=1,nat
            read(9,*)zan(i),cat(1,i),cat(2,i),cat(3,i)
            cat(1,i)=a2bohr*cat(1,i)
            cat(2,i)=a2bohr*cat(2,i)
            cat(3,i)=a2bohr*cat(3,i)
         end do

         read(9,*)nebf
         read(9,*)npebf
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
         do i=1,npebf
            read(9,*)idum,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),ELCAM(i,3),
     x    ELCEX(i),AGEBFCC(i),ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
            do j=1,3
               ELCBFC(i,j)=a2bohr*ELCBFC(i,j)
            end do
         end do

         read(9,*)npbf
         if(allocated(NUCEX)) deallocate(NUCEX)
         allocate( NUCEX(npbf),stat=istat )
         if(allocated(AGNBFCC)) deallocate(AGNBFCC)
         allocate( AGNBFCC(npbf),stat=istat )
         if(allocated(NUCAM)) deallocate(NUCAM)
         allocate( NUCAM(npbf,3),stat=istat )
         if(allocated(NUCBFC)) deallocate(NUCBFC)
         allocate( NUCBFC(npbf,3),stat=istat )
         do i=1,npbf
            read(9,*)idum,idum,NUCAM(i,1),NUCAM(i,2),NUCAM(i,3),
     x    NUCEX(i),AGNBFCC(i),NUCBFC(i,1),NUCBFC(i,2),NUCBFC(i,3)
            do j=1,3
               NUCBFC(i,j)=a2bohr*NUCBFC(i,j)
            end do
         end do
         if(allocated(KPESTR)) deallocate(KPESTR)
         allocate( KPESTR(nebf),stat=istat )
         if(allocated(KPEEND)) deallocate(KPEEND)
         allocate( KPEEND(nebf),stat=istat )
         call make_KPE(nebf,npebf,AMPEB2C,KPESTR,KPEEND)

         read(9,*)pmass
         read(9,*)nelec
         read(9,*)NAE
         read(9,*)NBE
         read(9,*) NUCST
         read(9,*) LNEOHF
         read(9,*) LXCUHF
         read(9,*) LXCROHF
         read(9,*) LRXCHF
         read(9,*) LRXCUHF
         read(9,*) read_CE
         read(9,*) read_CP
         read(9,*) read_GAM2
         read(9,*) read_GAM3
         read(9,*) read_GAM4
         read(9,*) NG4CHK
         read(9,*) LGAM4
         read(9,*) LG4DSCF
         read(9,*) LG4IC
         read(9,*) NG3CHK
         read(9,*) LG3DSCF
         read(9,*) LG3IC1
         read(9,*) LG3IC2
         read(9,*) NG2CHK
         read(9,*) LG2DSCF
         read(9,*) LG2IC1
         read(9,*) LCMF
         read(9,*) LSOSCF
         read(9,*) LOCBSE
         read(9,*) EXCHLEV ! 0=RXCHF-ne; 1=RXCHF-ae; 2=RXCHF-fe

         close(9)

      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Broadcast static variables including array dimensions to slave processes
      call MPI_BCAST(ngtg1,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nebf,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(npebf,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(npbf,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nat,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)

! INTEGER8 and LOGICAL8 must be declared specificially since MPI_INTEGER
! and MPI_LOGICAL do not seem to get specified as long with compiler flags 
      call MPI_BCAST(pmass,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nelec,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NAE,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NBE,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NUCST,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LNEOHF,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LXCUHF,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LXCROHF,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LRXCHF,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LRXCUHF,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(read_CE,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(read_CP,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(read_GAM2,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(read_GAM3,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(read_GAM4,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NG4CHK,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LGAM4,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG4DSCF,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG4IC,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NG3CHK,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG3DSCF,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG3IC1,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG3IC2,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NG2CHK,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG2DSCF,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LG2IC1,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LCMF,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LSOSCF,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LOCBSE,1,MPI_LOGICAL8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(EXCHLEV,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)


      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Allocate arrays on slave processes
      if (rank.gt.0) then
       if(allocated(bcoef1)) deallocate(bcoef1)
       allocate( bcoef1(ngtg1),stat=istat )
       if(allocated(gamma1)) deallocate(gamma1)
       allocate( gamma1(ngtg1),stat=istat )
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
      call MPI_BCAST(bcoef1,ngtg1,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(gamma1,ngtg1,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(zan,nat,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cat,3*nat,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(AMPEB2C,npebf,MPI_INTEGER8,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ELCEX,npebf,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(AGEBFCC,npebf,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ELCAM,3*npebf,MPI_INTEGER8,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ELCBFC,3*npebf,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NUCEX,npbf,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(AGNBFCC,npbf,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NUCAM,3*npbf,MPI_INTEGER8,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NUCBFC,3*npbf,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(KPESTR,nebf,MPI_INTEGER8,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(KPEEND,nebf,MPI_INTEGER8,
     x               0,MPI_COMM_WORLD,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(LNEOHF.or.LXCUHF.or.LRXCUHF) then
       if (rank.eq.0) 
     x    write(*,*) "MPI code only supports RXCHF-ne/ae with NBE>1"
       return
      end if

      if(.not.(LRXCHF.and.(NBE.ge.2).and.(EXCHLEV.lt.2))) then
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
       write(*,*)'nat   =',nat
       write(*,*)'npebf =',npebf
       write(*,*)'nebf  =',nebf
       write(*,*)'npbf  =',npbf
       write(*,*)'ngtg1 =',ngtg1
       write(*,*)'ngee  =',ngee
       write(*,*)'ng1   =',ng1
       write(*,*)'ng2   =',ng2
       write(*,*)'ng3   =',ng3
       write(*,*)'ng4   =',ng4
       write(*,*)'ng1prm=',ng1prm
       write(*,*)'ng2prm=',ng2prm
       write(*,*)'ng3prm=',ng3prm
       write(*,*)
       write(*,*)'PMASS   =',PMASS
       write(*,*)'nelec   =',nelec
       write(*,*)'NAE     =',NAE,'= total number of regular electrons'
       write(*,*)'NBE     =',NBE,'= number of special electrons'
       write(*,*)'NUCST   =',NUCST
       write(*,*)'LNEOHF  =',LNEOHF
       write(*,*)'LXCUHF  =',LXCUHF
       write(*,*)'LXCROHF =',LXCROHF
       write(*,*)'LRXCHF  =',LRXCHF
       write(*,*)'LRXCUHF  =',LRXCUHF
       write(*,*)'read_CE =',read_CE
       write(*,*)'read_CP =',read_CP
       write(*,*)'READ_GAM2=',read_GAM2
       write(*,*)'READ_GAM3=',read_GAM3
       write(*,*)'READ_GAM4=',read_GAM4
       write(*,*)'NG4CHK  =',NG4CHK
       write(*,*)'LGAM4   =',LGAM4
       write(*,*)'LG4DSCF =',LG4DSCF
       write(*,*)'LG4IC   =',LG4IC 
       write(*,*)'NG3CHK  =',NG3CHK
       write(*,*)'LG3DSCF =',LG3DSCF
       write(*,*)'LG3IC1  =',LG3IC1
       write(*,*)'LG3IC2  =',LG3IC2
       write(*,*)'NG2CHK  =',NG2CHK
       write(*,*)'LG2DSCF =',LG2DSCF
       write(*,*)'LG2IC1  =',LG2IC1
       write(*,*)'LCMF    =',LCMF
       write(*,*)'LSOSCF  =',LSOSCF
       write(*,*)'LOCBSE  =',LOCBSE
       write(*,*)'EXCHLEV =',EXCHLEV,' (=2: fe; =1: ae; =0: ne)'
       write(*,*) "Geminal parameters: k, b_k, gamm_k"
          do i=1,ngtg1
             write(*,*) i,bcoef1(i),gamma1(i)
          end do
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
     x                      nelec,nae,nbe,nucst,
     x                      nebf,npebf,npbf,nat,ngtg1,
     x                      ng1,ng2,ng3,ng4,ngee,
     x                      ng1prm,ng2prm,ng3prm,
     x                      pmass,cat,zan,bcoef1,gamma1,
     x                      KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                      NG2CHK,NG3CHK,NG4CHK,
     x                      read_CE,read_CP,
     x                      read_GAM2,read_GAM3,read_GAM4,
     x                      LG2IC1,LG3IC1,LG4IC,
     x                      LG2DSCF,LG3DSCF,LG4DSCF,
     x                      LSOSCF,LOCBSE,LCMF,EXCHLEV)

      write(*,*) "Process is finished:",rank

 7000 format(1X,A,1X,I4,1X,A,1X,A)
 8000 format(1X,I3,I6,I5)
 9000 format(1X,I3,I6,I5,I3,I3,F12.6,F10.6,F10.6,F10.6,F10.6)

      return
      end

