C=======================================================================
      subroutine RXCHF_OMG2_MPI(nproc,rank,
     x                          Nchunks,nebf,npbf,
     x                          ng2,ng2loc,
     x                          DAE,DBE,DP,
     x                          GM2_1,
     x                          FAE,FBE,FP,E_OMG2)

C Add interaction OMG2 contributions to all Fock matrices
C=======================================================================
      implicit none
      include 'mpif.h'

C Input variables
      integer           nproc,rank
      integer           Nchunks
      integer           ng2        ! Total number of integrals
      integer           ng2loc     ! Number of integrals on MPI proc
      integer           nebf,npbf
      double precision  GM2_1(ng2loc)  ! GAM2   integrals
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FAE(nebf,nebf)
      double precision  FBE(nebf,nebf)
      double precision  FP(npbf,npbf)
      double precision  E_OMG2

C Local variables
      integer           istat,ichunk,istart,iend,ng2_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  TFAE(nebf,nebf)
      double precision  TFBE(nebf,nebf)
      double precision  TFP(npbf,npbf)
      double precision  TE_OMG2

      integer,allocatable :: loop_map(:,:)

      integer   mpistart,mpiend,arrstart
#if MPI32
      integer*4 ierr,nebf2,npbf2
#else
      integer ierr,nebf2,npbf2
#endif

C Initialize
      XFAE      = 0.0d+00
      XFBE      = 0.0d+00
      XFP       = 0.0d+00
      E_OMG2    = 0.0d+00
      TFAE      = 0.0d+00
      TFBE      = 0.0d+00
      TFP       = 0.0d+00
      TE_OMG2   = 0.0d+00

#if MPI32
      nebf2=int(nebf*nebf,kind=4)
      npbf2=int(npbf*npbf,kind=4)
#else
      nebf2=nebf*nebf
      npbf2=npbf*npbf
#endif

! Each process has ng2/nproc integrals according to rank
! Last process also has ng2%nproc remaining integrals
      call get_mpi_range(ng2,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng2

C Chop up calculation of OMG2 terms
      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng2/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng2_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng2_seg,6),stat=istat )

C Nested loop compression for this chunk:
         Loopi=0
         imas=0
         do ip=1,npbf
         do jp=1,npbf
            do iec1=1,nebf
            do jec1=1,nebf
               do iec2=1,nebf
               do jec2=1,nebf

                  imas=imas+1 ! imas is master_index
                  if(imas.ge.istart.and.imas.le.iend) then
                     Loopi=Loopi+1
                     loop_map(Loopi,1)=jec2
                     loop_map(Loopi,2)=iec2
                     loop_map(Loopi,3)=jec1
                     loop_map(Loopi,4)=iec1
                     loop_map(Loopi,5)=jp
                     loop_map(Loopi,6)=ip

! Save coordinates of first value for array index shifts
                     if ((ichunk.eq.1).and.(Loopi.eq.1)) then
                        call index_GAM_2PK(nebf,npbf,
     x                           ip,jp,iec1,jec1,iec2,jec2,arrstart)
                     end if

                  end if

               end do
               end do
            end do
            end do
         end do
         end do

         call RXCHF_OMG2_thread_MPI(istart,iend,ng2_seg,ng2loc,
     x                              nebf,npbf,
     x                              loop_map,arrstart,
     x                              DAE,DBE,DP,
     x                              GM2_1,
     x                              TFAE,TFBE,TFP,TE_OMG2)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

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
      call MPI_ALLREDUCE(TE_OMG2,E_OMG2,1,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)

      call add2fock(nebf,XFAE,FAE)
      call add2fock(nebf,XFBE,FBE)
      call add2fock(npbf,XFP,FP)

      return
      end
C=======================================================================
      subroutine RXCHF_OMG2ex_MPI(nproc,rank,
     x                            Nchunks,nebf,npbf,
     x                            ng2,ng2loc,
     x                            DAE,DBE,DP,
     x                            GM2_1,GM2_2,
     x                            FAE,FBE,FP,E_OMG2)

C Add interaction OMG2 (incl exchange) contributions to all Fock matrices
C=======================================================================
      implicit none
      include 'mpif.h'

C Input variables
      integer           nproc,rank
      integer           Nchunks
      integer           ng2        ! Total number of integrals
      integer           ng2loc     ! Number of integrals on MPI proc
      integer           nebf,npbf
      double precision  GM2_1(ng2loc)  ! GAM2   integrals
      double precision  GM2_2(ng2loc)  ! GAM2ex integrals
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FAE(nebf,nebf)
      double precision  FBE(nebf,nebf)
      double precision  FP(npbf,npbf)
      double precision  E_OMG2

C Local variables
      integer           istat,ichunk,istart,iend,ng2_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  TFAE(nebf,nebf)
      double precision  TFBE(nebf,nebf)
      double precision  TFP(npbf,npbf)
      double precision  TE_OMG2

      integer,allocatable :: loop_map(:,:)

      integer   mpistart,mpiend,arrstart
#if MPI32
      integer*4 ierr,nebf2,npbf2
#else
      integer ierr,nebf2,npbf2
#endif

C Initialize
      XFAE      = 0.0d+00
      XFBE      = 0.0d+00
      XFP       = 0.0d+00
      E_OMG2    = 0.0d+00
      TFAE      = 0.0d+00
      TFBE      = 0.0d+00
      TFP       = 0.0d+00
      TE_OMG2   = 0.0d+00

#if MPI32
      nebf2=int(nebf*nebf,kind=4)
      npbf2=int(npbf*npbf,kind=4)
#else
      nebf2=nebf*nebf
      npbf2=npbf*npbf
#endif

! Each process has ng2/nproc integrals according to rank
! Last process also has ng2%nproc remaining integrals
      call get_mpi_range(ng2,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng2

C Chop up calculation of OMG2 terms
      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng2/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng2_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng2_seg,6),stat=istat )

C Nested loop compression for this chunk:
         Loopi=0
         imas=0
         do ip=1,npbf
         do jp=1,npbf
            do iec1=1,nebf
            do jec1=1,nebf
               do iec2=1,nebf
               do jec2=1,nebf

                  imas=imas+1 ! imas is master_index
                  if(imas.ge.istart.and.imas.le.iend) then
                     Loopi=Loopi+1
                     loop_map(Loopi,1)=jec2
                     loop_map(Loopi,2)=iec2
                     loop_map(Loopi,3)=jec1
                     loop_map(Loopi,4)=iec1
                     loop_map(Loopi,5)=jp
                     loop_map(Loopi,6)=ip

! Save coordinates of first value for array index shifts
                     if ((ichunk.eq.1).and.(Loopi.eq.1)) then
                        call index_GAM_2PK(nebf,npbf,
     x                           ip,jp,iec1,jec1,iec2,jec2,arrstart)
                     end if

                  end if

               end do
               end do
            end do
            end do
         end do
         end do

         call RXCHF_OMG2ex_thread_MPI(istart,iend,ng2_seg,ng2loc,
     x                                nebf,npbf,
     x                                loop_map,arrstart,
     x                                DAE,DBE,DP,
     x                                GM2_1,GM2_2,
     x                                TFAE,TFBE,TFP,TE_OMG2)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

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
      call MPI_ALLREDUCE(TE_OMG2,E_OMG2,1,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)

      call add2fock(nebf,XFAE,FAE)
      call add2fock(nebf,XFBE,FBE)
      call add2fock(npbf,XFP,FP)

      return
      end
C=======================================================================
      subroutine XCHF_OMG2_MPI(nproc,rank,
     x                         Nchunks,nebf,npbf,
     x                         ng2,ng2loc,
     x                         DE,DP,
     x                         GM2_1,GM2s,
     x                         FE,SE,FP,SP,E_OMG2,S_OMG2)

C Add XCHF OMG2 contributions to all Fock matrices
C=======================================================================
      implicit none
      include 'mpif.h'

C Input variables
      integer           nproc,rank
      integer           Nchunks
      integer           ng2        ! Total number of integrals
      integer           ng2loc     ! Number of integrals on MPI proc
      integer           nebf,npbf
      double precision  GM2_1(ng2loc)  ! GAM2   integrals
      double precision  GM2s(ng2loc)   ! GAM2s  integrals
      double precision  DE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FE(nebf,nebf)
      double precision  SE(nebf,nebf)
      double precision  FP(npbf,npbf)
      double precision  SP(npbf,npbf)
      double precision  E_OMG2
      double precision  S_OMG2

C Local variables
      integer           istat,ichunk,istart,iend,ng2_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2
      double precision  XFE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  TFE(nebf,nebf)
      double precision  TFP(npbf,npbf)
      double precision  TSE(nebf,nebf)
      double precision  TSP(npbf,npbf)
      double precision  TE_OMG2
      double precision  TS_OMG2

      integer,allocatable :: loop_map(:,:)

      integer   mpistart,mpiend,arrstart
#if MPI32
      integer*4 ierr,nebf2,npbf2
#else
      integer ierr,nebf2,npbf2
#endif

C Initialize
      XFE      = 0.0d+00
      SE       = 0.0d+00
      XFP      = 0.0d+00
      SP       = 0.0d+00
      E_OMG2   = 0.0d+00
      S_OMG2   = 0.0d+00
      TFE      = 0.0d+00
      TSE      = 0.0d+00
      TFP      = 0.0d+00
      TSP      = 0.0d+00
      TE_OMG2  = 0.0d+00
      TS_OMG2  = 0.0d+00

#if MPI32
      nebf2=int(nebf*nebf,kind=4)
      npbf2=int(npbf*npbf,kind=4)
#else
      nebf2=nebf*nebf
      npbf2=npbf*npbf
#endif

! Each process has ng2/nproc integrals according to rank
! Last process also has ng2%nproc remaining integrals
      call get_mpi_range(ng2,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng2

C Chop up calculation of OMG2 terms
      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng2/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng2_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng2_seg,6),stat=istat )

C Nested loop compression for this chunk:
         Loopi=0
         imas=0
         do ip=1,npbf
         do jp=1,npbf
            do iec1=1,nebf
            do jec1=1,nebf
               do iec2=1,nebf
               do jec2=1,nebf

                  imas=imas+1 ! imas is master_index
                  if(imas.ge.istart.and.imas.le.iend) then
                     Loopi=Loopi+1
                     loop_map(Loopi,1)=jec2
                     loop_map(Loopi,2)=iec2
                     loop_map(Loopi,3)=jec1
                     loop_map(Loopi,4)=iec1
                     loop_map(Loopi,5)=jp
                     loop_map(Loopi,6)=ip

! Save coordinates of first value for array index shifts
                     if ((ichunk.eq.1).and.(Loopi.eq.1)) then
                        call index_GAM_2PK(nebf,npbf,
     x                           ip,jp,iec1,jec1,iec2,jec2,arrstart)
                     end if

                  end if

               end do
               end do
            end do
            end do
         end do
         end do

         call XCHF_OMG2_thread_MPI(istart,iend,ng2_seg,ng2loc,
     x                             nebf,npbf,
     x                             loop_map,arrstart,
     x                             DE,DP,
     x                             GM2_1,GM2s,
     x                             TFE,TSE,TFP,TSP,TE_OMG2,TS_OMG2)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrices
C Use additional auxiliary array since MPI_IN_PLACE fails with MPICH2
      call MPI_ALLREDUCE(TFE(1,1),XFE(1,1),nebf2,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TSE(1,1),SE(1,1),nebf2,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TFP(1,1),XFP(1,1),npbf2,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TSP(1,1),SP(1,1),npbf2,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TE_OMG2,E_OMG2,1,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TS_OMG2,S_OMG2,1,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)

      call add2fock(nebf,XFE,FE)
      call add2fock(npbf,XFP,FP)

      return
      end
C======================================================================
      subroutine RXCHF_OMG2_thread_MPI(istart,iend,ng2_seg,ng2,
     x                                 nebf,npbf,
     x                                 loop_map,arrstart,
     x                                 DAE,DBE,DP,
     x                                 GM2_1,
     x                                 XFAE,XFBE,XFP,E_OMG2)
C======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng2_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng2
      integer           arrstart
      integer           loop_map(ng2_seg,6)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM2_1(ng2)  ! GAM2   integrals

C Output variables
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  E_OMG2

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1
      integer           iec2,jec2
      integer           imap,ia
      double precision  val
      double precision  wtime

C GM2ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron indicies

!$omp parallel 
!$ompx shared(loop_map,arrstart)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM2_1)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(ia)
!$ompx private(val)
!$ompx reduction(+:XFAE)
!$ompx reduction(+:XFBE)
!$ompx reduction(+:XFP)
!$ompx reduction(+:E_OMG2)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

       imap=iLp-istart+1
       jec2=loop_map(imap,1)
       iec2=loop_map(imap,2)
       jec1=loop_map(imap,3)
       iec1=loop_map(imap,4)
       jp =loop_map(imap,5)
       ip =loop_map(imap,6)

       call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec1,iec2,jec2,ia)

! GM2_1(ip,jp,ie1,je1,ie2,je2) contributes to F^e(ie1,je1),F^q(ie2,je2),F^p(ip,jp)
       val=GM2_1(ia-arrstart+1)

       XFAE(iec1,jec1)=XFAE(iec1,jec1)+DP(ip,jp)*DBE(iec2,jec2)*val
       XFP(ip,jp)=XFP(ip,jp)+DAE(iec1,jec1)*DBE(iec2,jec2)*val
       XFBE(iec2,jec2)=XFBE(iec2,jec2)+DP(ip,jp)*DAE(iec1,jec1)*val
       E_OMG2=E_OMG2+DP(ip,jp)*DAE(iec1,jec1)*DBE(iec2,jec2)*val

      end do
!$omp end do
!$omp end parallel      

      return
      end
C======================================================================
      subroutine RXCHF_OMG2ex_thread_MPI(istart,iend,ng2_seg,ng2,
     x                                   nebf,npbf,
     x                                   loop_map,arrstart,
     x                                   DAE,DBE,DP,
     x                                   GM2_1,GM2_2,
     x                                   XFAE,XFBE,XFP,E_OMG2)
C======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng2_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng2
      integer           arrstart
      integer           loop_map(ng2_seg,6)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM2_1(ng2)  ! GAM2   integrals
      double precision  GM2_2(ng2)  ! GAM2ex integrals

C Output variables
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  E_OMG2

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1
      integer           iec2,jec2
      integer           imap,ia
      double precision  val
      double precision  wtime
      double precision  half
      half=0.50d+00

C GM2ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron indicies

!$omp parallel 
!$ompx shared(half)
!$ompx shared(loop_map,arrstart)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM2_1,GM2_2)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(ia)
!$ompx private(val)
!$ompx reduction(+:XFAE)
!$ompx reduction(+:XFBE)
!$ompx reduction(+:XFP)
!$ompx reduction(+:E_OMG2)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

        imap=iLp-istart+1
        jec2=loop_map(imap,1)
        iec2=loop_map(imap,2)
        jec1=loop_map(imap,3)
        iec1=loop_map(imap,4)
        jp =loop_map(imap,5)
        ip =loop_map(imap,6)

        call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec1,iec2,jec2,ia)

! GM2_1(ip,jp,ie1,je1,ie2,je2) contributes to F^e(ie1,je1),F^q(ie2,je2),F^p(ip,jp)
        val=GM2_1(ia-arrstart+1)

        XFAE(iec1,jec1)=XFAE(iec1,jec1)+DP(ip,jp)*DBE(iec2,jec2)*val
        XFP(ip,jp)=XFP(ip,jp)+DAE(iec1,jec1)*DBE(iec2,jec2)*val
        XFBE(iec2,jec2)=XFBE(iec2,jec2)+DP(ip,jp)*DAE(iec1,jec1)*val
        E_OMG2=E_OMG2+DP(ip,jp)*DAE(iec1,jec1)*DBE(iec2,jec2)*val

! GM2_2(ip,jp,ie1,je1,ie2,je2) contributes to F^e(ie1,je2),F^q(ie2,je1),F^p(ip,jp)
        val=-half*GM2_2(ia-arrstart+1)

        XFAE(iec1,jec2)=XFAE(iec1,jec2)+DP(ip,jp)*DBE(iec2,jec1)*val
        XFP(ip,jp)=XFP(ip,jp)+DAE(iec1,jec2)*DBE(iec2,jec1)*val
        XFBE(iec2,jec1)=XFBE(iec2,jec1)+DP(ip,jp)*DAE(iec1,jec2)*val
        E_OMG2=E_OMG2+DP(ip,jp)*DAE(iec1,jec2)*DBE(iec2,jec1)*val
         
      end do
!$omp end do
!$omp end parallel      

      return
      end
C======================================================================
      subroutine XCHF_OMG2_thread_MPI(istart,iend,ng2_seg,ng2,
     x                                nebf,npbf,
     x                                loop_map,arrstart,
     x                                DE,DP,
     x                                GM2_1,GM2s,
     x                                XFE,SE,XFP,SP,E_OMG2,S_OMG2)
C======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng2_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng2
      integer           arrstart
      integer           loop_map(ng2_seg,6)
      double precision  DE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM2_1(ng2)  ! GAM2   integrals
      double precision  GM2s(ng2)   ! GAM2s  integrals

C Output variables
      double precision  XFE(nebf,nebf)
      double precision  SE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  SP(npbf,npbf)
      double precision  E_OMG2
      double precision  S_OMG2

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1
      integer           iec2,jec2
      integer           imap,ia
      double precision  val1,val2
      double precision  wtime
      double precision  two,half
      half=0.50d+00
      two=2.0d+00

C GM2ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron indicies

!$omp parallel 
!$ompx shared(two,half)
!$ompx shared(loop_map,arrstart)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx shared(DE,DP)
!$ompx shared(GM2_1,GM2s)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(ia)
!$ompx private(val1,val2)
!$ompx reduction(+:XFE)
!$ompx reduction(+:SE)
!$ompx reduction(+:XFP)
!$ompx reduction(+:SP)
!$ompx reduction(+:E_OMG2)
!$ompx reduction(+:S_OMG2)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

       imap=iLp-istart+1
       jec2=loop_map(imap,1)
       iec2=loop_map(imap,2)
       jec1=loop_map(imap,3)
       iec1=loop_map(imap,4)
       jp =loop_map(imap,5)
       ip =loop_map(imap,6)

       call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec1,iec2,jec2,ia)

! GM2_1(ip,jp,ie1,je1,ie2,je2) contributes
! twice to F^p(ip,jp) and once to F^e(ie1,je1),F^e(ie1,je2)
       val1=GM2_1(ia-arrstart+1)
       val2=-half*GM2_1(ia-arrstart+1)

       XFE(iec1,jec1)=XFE(iec1,jec1)+DP(ip,jp)*DE(iec2,jec2)*two*val1
       XFP(ip,jp)=XFP(ip,jp)+DE(iec1,jec1)*DE(iec2,jec2)*val1
       E_OMG2=E_OMG2+DP(ip,jp)*DE(iec1,jec1)*DE(iec2,jec2)*val1

       XFE(iec1,jec2)=XFE(iec1,jec2)+DP(ip,jp)*DE(iec2,jec1)*two*val2
       XFP(ip,jp)=XFP(ip,jp)+DE(iec1,jec2)*DE(iec2,jec1)*val2
       E_OMG2=E_OMG2+DP(ip,jp)*DE(iec1,jec2)*DE(iec2,jec1)*val2

! GM2s(ip,jp,ie1,je1,ie2,je2) contributes
! twice to S^p(ip,jp) and once to S^e(ie1,je1),S^e(ie1,je2)
       val1=GM2s(ia-arrstart+1)
       val2=-half*GM2s(ia-arrstart+1)

       SE(iec1,jec1)=SE(iec1,jec1)+DP(ip,jp)*DE(iec2,jec2)*val1
       SP(ip,jp)=SP(ip,jp)+DE(iec1,jec1)*DE(iec2,jec2)*val1
       S_OMG2=S_OMG2+DP(ip,jp)*DE(iec1,jec1)*DE(iec2,jec2)*val1

       SE(iec1,jec2)=SE(iec1,jec2)+DP(ip,jp)*DE(iec2,jec1)*val2
       SP(ip,jp)=SP(ip,jp)+DE(iec1,jec2)*DE(iec2,jec1)*val2
       S_OMG2=S_OMG2+DP(ip,jp)*DE(iec1,jec2)*DE(iec2,jec1)*val2

      end do
!$omp end do
!$omp end parallel      

      return
      end

