C=======================================================================
      subroutine RXCHF_OMG4_MPI(nproc,rank,
     x                          Nchunks,nebf,npbf,
     x                          ng4,
     x                          DAE,DBE,DP,
     x                          pinds,einds,GM4,
     x                          FAE,FBE,FP,E_OMG4)

C Add interaction OMG4 contributions to all Fock matrices
C=======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

C Input variables
      integer           nproc,rank
      integer           Nchunks
      integer           ng4        ! Number of nontrivial integrals on MPI proc
      integer           nebf,npbf
      integer(kind=4)   pinds(ng4)
      integer(kind=8)   einds(ng4)
      double precision  GM4(ng4)  ! GAM4 integrals
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FAE(nebf,nebf)
      double precision  FBE(nebf,nebf)
      double precision  FP(npbf,npbf)
      double precision  E_OMG4

C Local variables
      integer           i
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4
      double precision  val
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  TFAE(nebf,nebf)
      double precision  TFBE(nebf,nebf)
      double precision  TFP(npbf,npbf)
      double precision  TE_OMG4
      double precision  three

#if MPI32
      integer*4 ierr,nebf2,npbf2
#else
      integer ierr,nebf2,npbf2
#endif

      three=3.0d+00

C Initialize
      XFAE      = 0.0d+00
      XFBE      = 0.0d+00
      XFP       = 0.0d+00
      E_OMG4    = 0.0d+00
      TFAE      = 0.0d+00
      TFBE      = 0.0d+00
      TFP       = 0.0d+00
      TE_OMG4   = 0.0d+00

#if MPI32
      nebf2=int(nebf*nebf,kind=4)
      npbf2=int(npbf*npbf,kind=4)
#else
      nebf2=nebf*nebf
      npbf2=npbf*npbf
#endif

C Form intermediate Fock matrices

!$omp parallel 
!$ompx shared(three)
!$ompx shared(ng4)
!$ompx shared(nebf,npbf)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(pinds)
!$ompx shared(einds)
!$ompx shared(GM4)
!$ompx private(i) 
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(iec4,jec4)
!$ompx private(val)
!$ompx reduction(+:TFAE)
!$ompx reduction(+:TFBE)
!$ompx reduction(+:TFP)
!$ompx reduction(+:TE_OMG4)

!$omp do SCHEDULE(RUNTIME)
      do i=1,ng4

         call bitunpack(pinds(i),einds(i),ip,jp,
     x                  iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4)

! GM4(ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4) contributes to
! F^e(ie1,je1),F^q(ie2,je2),F^p(ip,jp)
         val=GM4(i)

         TFAE(iec1,jec1)=TFAE(iec1,jec1)+DP(ip,jp)*DBE(iec2,jec2)*
     x                   DBE(iec3,jec3)*DBE(iec4,jec4)*val
         TFP(ip,jp)=TFP(ip,jp)+DAE(iec1,jec1)*DBE(iec2,jec2)*
     x              DBE(iec3,jec3)*DBE(iec4,jec4)*val
         TFBE(iec2,jec2)=TFBE(iec2,jec2)+three*DAE(iec1,jec1)*
     x                   DBE(iec3,jec3)*DBE(iec4,jec4)*DP(ip,jp)*val
         TE_OMG4=TE_OMG4+DP(ip,jp)*DAE(iec1,jec1)*
     x             DBE(iec2,jec2)*DBE(iec3,jec3)*DBE(iec4,jec4)*val

      end do
!$omp end do
!$omp end parallel      

C Update Fock matrices
C Use additional auxiliary array since MPI_IN_PLACE fails with MPICH2
      call MPI_ALLREDUCE(TFAE(1,1),XFAE(1,1),nebf2,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TFBE(1,1),XFBE(1,1),nebf2,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TFP(1,1),XFP(1,1),npbf2,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TE_OMG4,E_OMG4,1,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)

      call add2fock(nebf,XFAE,FAE)
      call add2fock(nebf,XFBE,FBE)
      call add2fock(npbf,XFP,FP)

      return
      end
C=======================================================================
      subroutine XCHF_OMG4_MPI(nproc,rank,
     x                         Nchunks,nebf,npbf,
     x                         ng4,
     x                         DE,DP,
     x                         pinds,einds,GM4,
     x                         FE,FP,E_OMG4)

C Add XCHF OMG4 contributions to all Fock matrices
C=======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

C Input variables
      integer           nproc,rank
      integer           Nchunks
      integer           ng4        ! Number of nontrivial integrals on MPI proc
      integer           nebf,npbf
      integer(kind=4)   pinds(ng4)
      integer(kind=8)   einds(ng4)
      double precision  GM4(ng4)  ! GAM4 integrals
      double precision  DE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FE(nebf,nebf)
      double precision  FP(npbf,npbf)
      double precision  E_OMG4

C Local variables
      integer           i
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4
      double precision  val
      double precision  XFE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  TFE(nebf,nebf)
      double precision  TFP(npbf,npbf)
      double precision  TE_OMG4
      double precision  four

#if MPI32
      integer*4 ierr,nebf2,npbf2
#else
      integer ierr,nebf2,npbf2
#endif

      four=4.0d+00

C Initialize
      XFE       = 0.0d+00
      XFP       = 0.0d+00
      E_OMG4    = 0.0d+00
      TFE       = 0.0d+00
      TFP       = 0.0d+00
      TE_OMG4   = 0.0d+00

#if MPI32
      nebf2=int(nebf*nebf,kind=4)
      npbf2=int(npbf*npbf,kind=4)
#else
      nebf2=nebf*nebf
      npbf2=npbf*npbf
#endif

C Form intermediate Fock matrices

!$omp parallel 
!$ompx shared(four)
!$ompx shared(ng4)
!$ompx shared(nebf,npbf)
!$ompx shared(DE,DP)
!$ompx shared(pinds)
!$ompx shared(einds)
!$ompx shared(GM4)
!$ompx private(i) 
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(iec4,jec4)
!$ompx private(val)
!$ompx reduction(+:TFE)
!$ompx reduction(+:TFP)
!$ompx reduction(+:TE_OMG4)

!$omp do SCHEDULE(RUNTIME)
      do i=1,ng4

         call bitunpack(pinds(i),einds(i),ip,jp,
     x                  iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4)

! GM4(ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4) contributes to
! F^e(ie1,je1),F^p(ip,jp)
         val=GM4(i)

         TFE(iec1,jec1)=TFE(iec1,jec1)+DP(ip,jp)*DE(iec2,jec2)*
     x                   DE(iec3,jec3)*DE(iec4,jec4)*four*val
         TFP(ip,jp)=TFP(ip,jp)+DE(iec1,jec1)*DE(iec2,jec2)*
     x              DE(iec3,jec3)*DE(iec4,jec4)*val
         TE_OMG4=TE_OMG4+DP(ip,jp)*DE(iec1,jec1)*
     x             DE(iec2,jec2)*DE(iec3,jec3)*DE(iec4,jec4)*val

      end do
!$omp end do
!$omp end parallel      

C Update Fock matrices
C Use additional auxiliary array since MPI_IN_PLACE fails with MPICH2
      call MPI_ALLREDUCE(TFE(1,1),XFE(1,1),nebf2,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TFP(1,1),XFP(1,1),npbf2,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TE_OMG4,E_OMG4,1,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)

      call add2fock(nebf,XFE,FE)
      call add2fock(npbf,XFP,FP)

      return
      end

