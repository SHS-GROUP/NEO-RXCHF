C=======================================================================
      subroutine RXCHF_OMG3_MPI(nproc,rank,
     x                          Nchunks,nebf,npbf,
     x                          ng3,ng3loc,
     x                          DAE,DBE,DP,
     x                          GM3_1,
     x                          FAE,FBE,FP,E_OMG3)

C Add interaction OMG3 contributions to all Fock matrices
C=======================================================================
      implicit none
      include 'mpif.h'

C Input variables
      integer           nproc,rank
      integer           Nchunks
      integer           ng3        ! Total number of integrals
      integer           ng3loc     ! Number of integrals on MPI proc
      integer           nebf,npbf
      double precision  GM3_1(ng3loc)  ! GAM3    integrals
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FAE(nebf,nebf)
      double precision  FBE(nebf,nebf)
      double precision  FP(npbf,npbf)
      double precision  E_OMG3

C Local variables
      integer           istat,ichunk,istart,iend,ng3_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  TFAE(nebf,nebf)
      double precision  TFBE(nebf,nebf)
      double precision  TFP(npbf,npbf)
      double precision  TE_OMG3

      integer,allocatable :: loop_map(:,:)

      integer   mpistart,mpiend,arrstart
      integer ierr

C Initialize
      XFAE      = 0.0d+00
      XFBE      = 0.0d+00
      XFP       = 0.0d+00
      E_OMG3    = 0.0d+00
      TFAE      = 0.0d+00
      TFBE      = 0.0d+00
      TFP       = 0.0d+00
      TE_OMG3   = 0.0d+00

! Each process has ng3/nproc integrals according to rank
! Last process also has ng3%nproc remaining integrals
      call get_mpi_range(ng3,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng3

C Chop up calculation of OMG3 terms
      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng3/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng3_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )

C Nested loop compression for this chunk:
         Loopi=0
         imas=0
         do ip=1,npbf
         do jp=1,npbf
            do iec1=1,nebf
            do jec1=1,nebf
               do iec2=1,nebf
               do jec2=1,nebf
                  do iec3=1,nebf
                  do jec3=1,nebf

                         imas=imas+1 ! imas is master_index
                         if(imas.ge.istart.and.imas.le.iend) then
                            Loopi=Loopi+1
                            loop_map(Loopi,1)=jec3
                            loop_map(Loopi,2)=iec3
                            loop_map(Loopi,3)=jec2
                            loop_map(Loopi,4)=iec2
                            loop_map(Loopi,5)=jec1
                            loop_map(Loopi,6)=iec1
                            loop_map(Loopi,7)=jp
                            loop_map(Loopi,8)=ip

! Save coordinates of first value for array index shifts
                            if ((ichunk.eq.1).and.(Loopi.eq.1)) then
                               call index_GAM_3PK(nebf,npbf,
     x                     ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,arrstart)
                            end if

                         end if

                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

         call RXCHF_OMG3_thread_MPI(istart,iend,ng3_seg,ng3loc,
     x                              nebf,npbf,
     x                              loop_map,arrstart,
     x                              DAE,DBE,DP,
     x                              GM3_1,
     x                              TFAE,TFBE,TFP,TE_OMG3)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
C Use additional auxiliary array since MPI_IN_PLACE fails with MPICH2
      call MPI_ALLREDUCE(TFAE(1,1),XFAE(1,1),nebf*nebf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TFBE(1,1),XFBE(1,1),nebf*nebf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TFP(1,1),XFP(1,1),npbf*npbf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TE_OMG3,E_OMG3,1,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)

      call add2fock(nebf,XFAE,FAE)
      call add2fock(nebf,XFBE,FBE)
      call add2fock(npbf,XFP,FP)

      return
      end
C=======================================================================
      subroutine RXCHF_OMG3ex_MPI(nproc,rank,
     x                            Nchunks,nebf,npbf,
     x                            ng3,ng3loc,
     x                            DAE,DBE,DP,
     x                            GM3_1,GM3_2,GM3_3,
     x                            FAE,FBE,FP,E_OMG3)

C Add interaction OMG3 (incl exchange) contributions to all Fock matrices
C=======================================================================
      implicit none
      include 'mpif.h'

C Input variables
      integer           nproc,rank
      integer           Nchunks
      integer           ng3        ! Total number of integrals
      integer           ng3loc     ! Number of integrals on MPI proc
      integer           nebf,npbf
      double precision  GM3_1(ng3loc)  ! GAM3    integrals
      double precision  GM3_2(ng3loc)  ! GAM3ex1 integrals
      double precision  GM3_3(ng3loc)  ! GAM3ex2 integrals
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FAE(nebf,nebf)
      double precision  FBE(nebf,nebf)
      double precision  FP(npbf,npbf)
      double precision  E_OMG3

C Local variables
      integer           istat,ichunk,istart,iend,ng3_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  TFAE(nebf,nebf)
      double precision  TFBE(nebf,nebf)
      double precision  TFP(npbf,npbf)
      double precision  TE_OMG3

      integer,allocatable :: loop_map(:,:)

      integer   mpistart,mpiend,arrstart
      integer ierr

C Initialize
      XFAE      = 0.0d+00
      XFBE      = 0.0d+00
      XFP       = 0.0d+00
      E_OMG3    = 0.0d+00
      TFAE      = 0.0d+00
      TFBE      = 0.0d+00
      TFP       = 0.0d+00
      TE_OMG3   = 0.0d+00

! Each process has ng3/nproc integrals according to rank
! Last process also has ng3%nproc remaining integrals
      call get_mpi_range(ng3,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng3

C Chop up calculation of OMG3 terms
      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng3/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng3_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )

C Nested loop compression for this chunk:
         Loopi=0
         imas=0
         do ip=1,npbf
         do jp=1,npbf
            do iec1=1,nebf
            do jec1=1,nebf
               do iec2=1,nebf
               do jec2=1,nebf
                  do iec3=1,nebf
                  do jec3=1,nebf

                         imas=imas+1 ! imas is master_index
                         if(imas.ge.istart.and.imas.le.iend) then
                            Loopi=Loopi+1
                            loop_map(Loopi,1)=jec3
                            loop_map(Loopi,2)=iec3
                            loop_map(Loopi,3)=jec2
                            loop_map(Loopi,4)=iec2
                            loop_map(Loopi,5)=jec1
                            loop_map(Loopi,6)=iec1
                            loop_map(Loopi,7)=jp
                            loop_map(Loopi,8)=ip

! Save coordinates of first value for array index shifts
                            if ((ichunk.eq.1).and.(Loopi.eq.1)) then
                               call index_GAM_3PK(nebf,npbf,
     x                     ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,arrstart)
                            end if

                         end if

                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

         call RXCHF_OMG3ex_thread_MPI(istart,iend,ng3_seg,ng3loc,
     x                                nebf,npbf,
     x                                loop_map,arrstart,
     x                                DAE,DBE,DP,
     x                                GM3_1,GM3_2,GM3_3,
     x                                TFAE,TFBE,TFP,TE_OMG3)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
C Use additional auxiliary array since MPI_IN_PLACE fails with MPICH2
      call MPI_ALLREDUCE(TFAE(1,1),XFAE(1,1),nebf*nebf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TFBE(1,1),XFBE(1,1),nebf*nebf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TFP(1,1),XFP(1,1),npbf*npbf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TE_OMG3,E_OMG3,1,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)

      call add2fock(nebf,XFAE,FAE)
      call add2fock(nebf,XFBE,FBE)
      call add2fock(npbf,XFP,FP)

      return
      end
C=======================================================================
      subroutine XCHF_OMG3_MPI(nproc,rank,
     x                         Nchunks,nebf,npbf,
     x                         ng3,ng3loc,
     x                         DE,DP,
     x                         GM3_1,
     x                         FE,FP,E_OMG3)

C Add XCHF OMG3 contributions to all Fock matrices
C=======================================================================
      implicit none
      include 'mpif.h'

C Input variables
      integer           nproc,rank
      integer           Nchunks
      integer           ng3        ! Total number of integrals
      integer           ng3loc     ! Number of integrals on MPI proc
      integer           nebf,npbf
      double precision  GM3_1(ng3loc)  ! GAM3    integrals
      double precision  DE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FE(nebf,nebf)
      double precision  FP(npbf,npbf)
      double precision  E_OMG3

C Local variables
      integer           istat,ichunk,istart,iend,ng3_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      double precision  XFE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  TFE(nebf,nebf)
      double precision  TFP(npbf,npbf)
      double precision  TE_OMG3

      integer,allocatable :: loop_map(:,:)

      integer   mpistart,mpiend,arrstart
      integer ierr

C Initialize
      XFE       = 0.0d+00
      XFP       = 0.0d+00
      E_OMG3    = 0.0d+00
      TFE       = 0.0d+00
      TFP       = 0.0d+00
      TE_OMG3   = 0.0d+00

! Each process has ng3/nproc integrals according to rank
! Last process also has ng3%nproc remaining integrals
      call get_mpi_range(ng3,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng3

C Chop up calculation of OMG3 terms
      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng3/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng3_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )

C Nested loop compression for this chunk:
         Loopi=0
         imas=0
         do ip=1,npbf
         do jp=1,npbf
            do iec1=1,nebf
            do jec1=1,nebf
               do iec2=1,nebf
               do jec2=1,nebf
                  do iec3=1,nebf
                  do jec3=1,nebf

                         imas=imas+1 ! imas is master_index
                         if(imas.ge.istart.and.imas.le.iend) then
                            Loopi=Loopi+1
                            loop_map(Loopi,1)=jec3
                            loop_map(Loopi,2)=iec3
                            loop_map(Loopi,3)=jec2
                            loop_map(Loopi,4)=iec2
                            loop_map(Loopi,5)=jec1
                            loop_map(Loopi,6)=iec1
                            loop_map(Loopi,7)=jp
                            loop_map(Loopi,8)=ip

! Save coordinates of first value for array index shifts
                            if ((ichunk.eq.1).and.(Loopi.eq.1)) then
                               call index_GAM_3PK(nebf,npbf,
     x                     ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,arrstart)
                            end if

                         end if

                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

         call XCHF_OMG3_thread_MPI(istart,iend,ng3_seg,ng3loc,
     x                             nebf,npbf,
     x                             loop_map,arrstart,
     x                             DE,DP,
     x                             GM3_1,
     x                             TFE,TFP,TE_OMG3)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
C Use additional auxiliary array since MPI_IN_PLACE fails with MPICH2
      call MPI_ALLREDUCE(TFE(1,1),XFE(1,1),nebf*nebf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TFP(1,1),XFP(1,1),npbf*npbf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(TE_OMG3,E_OMG3,1,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)

      call add2fock(nebf,XFE,FE)
      call add2fock(npbf,XFP,FP)

      return
      end
C=======================================================================
      subroutine RXCHF_OMG3_thread_MPI(istart,iend,ng3_seg,ng3,
     x                                 nebf,npbf,
     x                                 loop_map,arrstart,
     x                                 DAE,DBE,DP,
     x                                 GM3_1,
     x                                 XFAE,XFBE,XFP,E_OMG3)
C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng3_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng3
      integer           arrstart
      integer           loop_map(ng3_seg,8)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM3_1(ng3)  ! GAM3    integrals

C Output variables
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  E_OMG3

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1  !
      integer           iec2,jec2  ! Contracted elec basis function indices
      integer           iec3,jec3  !
      integer           imap,ia
      double precision  val1,val2
      double precision  wtime
      double precision  half,two
      half=0.50d+00
      two=2.0d+00

C GM3ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron 1 indicies
C    iec3,jec3 : special electron 2 indicies

!$omp parallel 
!$ompx shared(half,two)
!$ompx shared(loop_map,arrstart)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM3_1)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia) 
!$ompx private(val1,val2)
!$ompx reduction(+:XFAE)
!$ompx reduction(+:XFBE)
!$ompx reduction(+:XFP)
!$ompx reduction(+:E_OMG3)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

        imap=iLp-istart+1
        jec3=loop_map(imap,1)
        iec3=loop_map(imap,2)
        jec2=loop_map(imap,3)
        iec2=loop_map(imap,4)
        jec1=loop_map(imap,5)
        iec1=loop_map(imap,6)
        jp =loop_map(imap,7)
        ip =loop_map(imap,8)

        call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)

! GM3_1(ip,jp,ie1,je1,ie2,je2,ie3,je3) contributes
! twice to F^e(ie1,je1),F^p(ip,jp) and once to F^q(ie2,je2),F^q(ie2,je3)
        val1=GM3_1(ia-arrstart+1)
        val2=-half*GM3_1(ia-arrstart+1)

        XFAE(iec1,jec1)=XFAE(iec1,jec1)+val1*
     x                  DBE(iec2,jec2)*DBE(iec3,jec3)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val1*
     x             DAE(iec1,jec1)*DBE(iec2,jec2)*DBE(iec3,jec3)
        XFBE(iec2,jec2)=XFBE(iec2,jec2)+two*val1*
     x                  DAE(iec1,jec1)*DBE(iec3,jec3)*DP(ip,jp)
        E_OMG3=E_OMG3+DP(ip,jp)*DAE(iec1,jec1)*
     x            DBE(iec2,jec2)*DBE(iec3,jec3)*val1

        XFAE(iec1,jec1)=XFAE(iec1,jec1)+val2*
     x                  DBE(iec2,jec3)*DBE(iec3,jec2)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val2*
     x             DAE(iec1,jec1)*DBE(iec2,jec3)*DBE(iec3,jec2)
        XFBE(iec2,jec3)=XFBE(iec2,jec3)+two*val2*
     x                  DAE(iec1,jec1)*DBE(iec3,jec2)*DP(ip,jp)
        E_OMG3=E_OMG3+DP(ip,jp)*DAE(iec1,jec1)*
     x            DBE(iec2,jec3)*DBE(iec3,jec2)*val2

      end do
!$omp end do
!$omp end parallel      

      return
      end
C=======================================================================
      subroutine RXCHF_OMG3ex_thread_MPI(istart,iend,ng3_seg,ng3,
     x                                   nebf,npbf,
     x                                   loop_map,arrstart,
     x                                   DAE,DBE,DP,
     x                                   GM3_1,GM3_2,GM3_3,
     x                                   XFAE,XFBE,XFP,E_OMG3)
C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng3_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng3
      integer           arrstart
      integer           loop_map(ng3_seg,8)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM3_1(ng3)  ! GAM3    integrals
      double precision  GM3_2(ng3)  ! GAM3ex1 integrals
      double precision  GM3_3(ng3)  ! GAM3ex2 integrals

C Output variables
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  E_OMG3

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1  !
      integer           iec2,jec2  ! Contracted elec basis function indices
      integer           iec3,jec3  !
      integer           imap,ia
      double precision  val1,val2
      double precision  wtime
      double precision  half,two
      half=0.50d+00
      two=2.0d+00

C GM3ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron 1 indicies
C    iec3,jec3 : special electron 2 indicies

!$omp parallel 
!$ompx shared(half,two)
!$ompx shared(loop_map,arrstart)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM3_1,GM3_2,GM3_3)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia) 
!$ompx private(val1,val2)
!$ompx reduction(+:XFAE)
!$ompx reduction(+:XFBE)
!$ompx reduction(+:XFP)
!$ompx reduction(+:E_OMG3)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

        imap=iLp-istart+1
        jec3=loop_map(imap,1)
        iec3=loop_map(imap,2)
        jec2=loop_map(imap,3)
        iec2=loop_map(imap,4)
        jec1=loop_map(imap,5)
        iec1=loop_map(imap,6)
        jp =loop_map(imap,7)
        ip =loop_map(imap,8)

        call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)

! GM3_1(ip,jp,ie1,je1,ie2,je2,ie3,je3) contributes
! twice to F^e(ie1,je1),F^p(ip,jp) and once to F^q(ie2,je2),F^q(ie2,je3)
        val1=GM3_1(ia-arrstart+1)
        val2=-half*GM3_1(ia-arrstart+1)

        XFAE(iec1,jec1)=XFAE(iec1,jec1)+val1*
     x                  DBE(iec2,jec2)*DBE(iec3,jec3)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val1*
     x             DAE(iec1,jec1)*DBE(iec2,jec2)*DBE(iec3,jec3)
        XFBE(iec2,jec2)=XFBE(iec2,jec2)+two*val1*
     x                  DAE(iec1,jec1)*DBE(iec3,jec3)*DP(ip,jp)
        E_OMG3=E_OMG3+DP(ip,jp)*DAE(iec1,jec1)*
     x            DBE(iec2,jec2)*DBE(iec3,jec3)*val1

        XFAE(iec1,jec1)=XFAE(iec1,jec1)+val2*
     x                  DBE(iec2,jec3)*DBE(iec3,jec2)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val2*
     x             DAE(iec1,jec1)*DBE(iec2,jec3)*DBE(iec3,jec2)
        XFBE(iec2,jec3)=XFBE(iec2,jec3)+two*val2*
     x                  DAE(iec1,jec1)*DBE(iec3,jec2)*DP(ip,jp)
        E_OMG3=E_OMG3+DP(ip,jp)*DAE(iec1,jec1)*
     x            DBE(iec2,jec3)*DBE(iec3,jec2)*val2

! GM3_2(ip,jp,ie1,je1,ie2,je2,ie3,je3) contributes
! twice to F^e(ie1,je2),F^p(ip,jp) and
! once to F^q(ie2,je1),F^q(ie2,je3),F^q(ie3,je3),F^q(ie3,je1)
        val1=-half*GM3_2(ia-arrstart+1)
        val2=half*half*GM3_2(ia-arrstart+1)

        XFAE(iec1,jec2)=XFAE(iec1,jec2)+val1*
     x                  DBE(iec2,jec1)*DBE(iec3,jec3)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val1*
     x             DAE(iec1,jec2)*DBE(iec2,jec1)*DBE(iec3,jec3)
        XFBE(iec2,jec1)=XFBE(iec2,jec1)+val1*
     x                  DAE(iec1,jec2)*DBE(iec3,jec3)*DP(ip,jp)
        E_OMG3=E_OMG3+DP(ip,jp)*DAE(iec1,jec2)*
     x            DBE(iec2,jec1)*DBE(iec3,jec3)*val1

        XFAE(iec1,jec2)=XFAE(iec1,jec2)+val2*
     x                  DBE(iec2,jec3)*DBE(iec3,jec1)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val2*
     x             DAE(iec1,jec2)*DBE(iec2,jec3)*DBE(iec3,jec1)
        XFBE(iec2,jec3)=XFBE(iec2,jec3)+val2*
     x                  DAE(iec1,jec2)*DBE(iec3,jec1)*DP(ip,jp)
        E_OMG3=E_OMG3+DP(ip,jp)*DAE(iec1,jec2)*
     x            DBE(iec2,jec3)*DBE(iec3,jec1)*val2

        XFBE(iec3,jec3)=XFBE(iec3,jec3)+val1*
     x                  DAE(iec1,jec2)*DBE(iec2,jec1)*DP(ip,jp)

        XFBE(iec3,jec1)=XFBE(iec3,jec1)+val2*
     x                  DAE(iec1,jec2)*DBE(iec2,jec3)*DP(ip,jp)

! GM3_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) contributes
! twice to F^e(ie1,je3),F^p(ip,jp) and
! once to F^q(ie2,je2),F^q(ie2,je1),F^q(ie3,je1),F^q(ie3,je2)
        val1=-half*GM3_3(ia-arrstart+1)
        val2=half*half*GM3_3(ia-arrstart+1)

        XFAE(iec1,jec3)=XFAE(iec1,jec3)+val1*
     x                  DBE(iec2,jec2)*DBE(iec3,jec1)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val1*
     x             DAE(iec1,jec3)*DBE(iec2,jec2)*DBE(iec3,jec1)
        XFBE(iec2,jec2)=XFBE(iec2,jec2)+val1*
     x                  DAE(iec1,jec3)*DBE(iec3,jec1)*DP(ip,jp)
        E_OMG3=E_OMG3+DP(ip,jp)*DAE(iec1,jec3)*
     x            DBE(iec2,jec2)*DBE(iec3,jec1)*val1

        XFAE(iec1,jec3)=XFAE(iec1,jec3)+val2*
     x                  DBE(iec2,jec1)*DBE(iec3,jec2)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val2*
     x             DAE(iec1,jec3)*DBE(iec2,jec1)*DBE(iec3,jec2)
        XFBE(iec2,jec1)=XFBE(iec2,jec1)+val2*
     x                  DAE(iec1,jec3)*DBE(iec3,jec2)*DP(ip,jp)
        E_OMG3=E_OMG3+DP(ip,jp)*DAE(iec1,jec3)*
     x            DBE(iec2,jec1)*DBE(iec3,jec2)*val2

        XFBE(iec3,jec1)=XFBE(iec3,jec1)+val1*
     x                  DAE(iec1,jec3)*DBE(iec2,jec2)*DP(ip,jp)

        XFBE(iec3,jec2)=XFBE(iec3,jec2)+val2*
     x                  DAE(iec1,jec3)*DBE(iec2,jec1)*DP(ip,jp)

      end do
!$omp end do
!$omp end parallel      

      return
      end
C=======================================================================
      subroutine XCHF_OMG3_thread_MPI(istart,iend,ng3_seg,ng3,
     x                                nebf,npbf,
     x                                loop_map,arrstart,
     x                                DE,DP,
     x                                GM3_1,
     x                                XFE,XFP,E_OMG3)
C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng3_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng3
      integer           arrstart
      integer           loop_map(ng3_seg,8)
      double precision  DE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM3_1(ng3)  ! GAM3    integrals

C Output variables
      double precision  XFE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  E_OMG3

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1  !
      integer           iec2,jec2  ! Contracted elec basis function indices
      integer           iec3,jec3  !
      integer           imap,ia
      double precision  val1,val2,val3
      double precision  wtime
      double precision  half,two,three
      half=0.50d+00
      two=2.0d+00
      three=3.0d+00

C GM3ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron 1 indicies
C    iec3,jec3 : special electron 2 indicies

!$omp parallel 
!$ompx shared(half,two,three)
!$ompx shared(loop_map,arrstart)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DE,DP)
!$ompx shared(GM3_1)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia) 
!$ompx private(val1,val2,val3)
!$ompx reduction(+:XFE)
!$ompx reduction(+:XFP)
!$ompx reduction(+:E_OMG3)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

        imap=iLp-istart+1
        jec3=loop_map(imap,1)
        iec3=loop_map(imap,2)
        jec2=loop_map(imap,3)
        iec2=loop_map(imap,4)
        jec1=loop_map(imap,5)
        iec1=loop_map(imap,6)
        jp =loop_map(imap,7)
        ip =loop_map(imap,8)

        call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)

! GM3_1(ip,jp,ie1,je1,ie2,je2,ie3,je3) contributes
! six times to F^p(ip,jp) and twice to F^e(ie1,je1),F^e(ie1,je2),F^e(ie1,je3)
        val1=GM3_1(ia-arrstart+1)
        val2=-half*GM3_1(ia-arrstart+1)
        val3=half*half*GM3_1(ia-arrstart+1)

        XFE(iec1,jec1)=XFE(iec1,jec1)+three*val1*
     x                  DE(iec2,jec2)*DE(iec3,jec3)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val1*
     x             DE(iec1,jec1)*DE(iec2,jec2)*DE(iec3,jec3)
        E_OMG3=E_OMG3+DP(ip,jp)*DE(iec1,jec1)*
     x            DE(iec2,jec2)*DE(iec3,jec3)*val1

        XFE(iec1,jec2)=XFE(iec1,jec2)+three*val2*
     x                  DE(iec2,jec1)*DE(iec3,jec3)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val2*
     x             DE(iec1,jec2)*DE(iec2,jec1)*DE(iec3,jec3)
        E_OMG3=E_OMG3+DP(ip,jp)*DE(iec1,jec2)*
     x            DE(iec2,jec1)*DE(iec3,jec3)*val2

        XFE(iec1,jec2)=XFE(iec1,jec2)+three*val3*
     x                  DE(iec2,jec3)*DE(iec3,jec1)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val3*
     x             DE(iec1,jec2)*DE(iec2,jec3)*DE(iec3,jec1)
        E_OMG3=E_OMG3+DP(ip,jp)*DE(iec1,jec2)*
     x            DE(iec2,jec3)*DE(iec3,jec1)*val3

        XFE(iec1,jec1)=XFE(iec1,jec1)+three*val2*
     x                  DE(iec2,jec3)*DE(iec3,jec2)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val2*
     x             DE(iec1,jec1)*DE(iec2,jec3)*DE(iec3,jec2)
        E_OMG3=E_OMG3+DP(ip,jp)*DE(iec1,jec1)*
     x            DE(iec2,jec3)*DE(iec3,jec2)*val2

        XFE(iec1,jec3)=XFE(iec1,jec3)+three*val3*
     x                  DE(iec2,jec1)*DE(iec3,jec2)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val3*
     x             DE(iec1,jec3)*DE(iec2,jec1)*DE(iec3,jec2)
        E_OMG3=E_OMG3+DP(ip,jp)*DE(iec1,jec3)*
     x            DE(iec2,jec1)*DE(iec3,jec2)*val3

        XFE(iec1,jec3)=XFE(iec1,jec3)+three*val2*
     x                  DE(iec2,jec2)*DE(iec3,jec1)*DP(ip,jp)
        XFP(ip,jp)=XFP(ip,jp)+val2*
     x             DE(iec1,jec3)*DE(iec2,jec2)*DE(iec3,jec1)
        E_OMG3=E_OMG3+DP(ip,jp)*DE(iec1,jec3)*
     x            DE(iec2,jec2)*DE(iec3,jec1)*val2

      end do
!$omp end do
!$omp end parallel      

      return
      end

