C=======================================================================
      subroutine RXCHF_OMG4_MPI(nproc,rank,
     x                          Nchunks,nebf,npbf,
     x                          ng4,ng4loc,
     x                          DAE,DBE,DP,
     x                          GM4,
     x                          FAE,FBE,FP,E_OMG4)

C Add interaction OMG4 contributions to all Fock matrices
C=======================================================================
      implicit none
      include 'mpif.h'

C Input variables
      integer           nproc,rank
      integer           Nchunks
      integer           ng4        ! Total number of integrals
      integer           ng4loc     ! Number of integrals on MPI proc
      integer           nebf,npbf
      double precision  GM4(ng4loc)  ! GAM4 integrals
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FAE(nebf,nebf)
      double precision  FBE(nebf,nebf)
      double precision  FP(npbf,npbf)
      double precision  E_OMG4

C Local variables
      integer           istat,ichunk,istart,iend,ng4_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)

      integer,allocatable :: loop_map(:,:)

      integer   mpistart,mpiend,arrstart
      integer*4 ierr

C Initialize
      XFAE      = 0.0d+00
      XFBE      = 0.0d+00
      XFP       = 0.0d+00
      E_OMG4    = 0.0d+00

! Each process has ng4/nproc integrals according to rank
! Last process also has ng4%nproc remaining integrals
      call get_mpi_range(ng4,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng4

C Chop up calculation of OMG4 terms
      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng4/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )

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
                     do iec4=1,nebf
                     do jec4=1,nebf

                            imas=imas+1 ! imas is master_index
                            if(imas.ge.istart.and.imas.le.iend) then
                               Loopi=Loopi+1
                               loop_map(Loopi,1)=jec4
                               loop_map(Loopi,2)=iec4
                               loop_map(Loopi,3)=jec3
                               loop_map(Loopi,4)=iec3
                               loop_map(Loopi,5)=jec2
                               loop_map(Loopi,6)=iec2
                               loop_map(Loopi,7)=jec1
                               loop_map(Loopi,8)=iec1
                               loop_map(Loopi,9)=jp
                               loop_map(Loopi,10)=ip

! Save coordinates of first value for array index shifts
                     if ((ichunk.eq.1).and.(Loopi.eq.1)) then
                        call index_GAM_4PK2(nebf,npbf,ip,jp,iec1,jec1,
     x                           iec2,jec2,iec3,jec3,iec4,jec4,arrstart)
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
         end do
         end do

         call RXCHF_OMG4_thread_MPI(istart,iend,ng4_seg,ng4loc,
     x                              nebf,npbf,
     x                              loop_map,arrstart,
     x                              DAE,DBE,DP,
     x                              GM4,
     x                              XFAE,XFBE,XFP,E_OMG4)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrices
      call MPI_ALLREDUCE(MPI_IN_PLACE,XFAE,nebf*nebf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,XFBE,nebf*nebf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,XFP,npbf*npbf,
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
     x                         ng4,ng4loc,
     x                         DE,DP,
     x                         GM4,
     x                         FE,FP,E_OMG4)

C Add XCHF OMG4 contributions to all Fock matrices
C=======================================================================
      implicit none
      include 'mpif.h'

C Input variables
      integer           nproc,rank
      integer           Nchunks
      integer           ng4        ! Total number of integrals
      integer           ng4loc     ! Number of integrals on MPI proc
      integer           nebf,npbf
      double precision  GM4(ng4loc)  ! GAM4 integrals
      double precision  DE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FE(nebf,nebf)
      double precision  FP(npbf,npbf)
      double precision  E_OMG4

C Local variables
      integer           istat,ichunk,istart,iend,ng4_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4
      double precision  XFE(nebf,nebf)
      double precision  XFP(npbf,npbf)

      integer,allocatable :: loop_map(:,:)

      integer   mpistart,mpiend,arrstart
      integer*4 ierr

C Initialize
      XFE       = 0.0d+00
      XFP       = 0.0d+00
      E_OMG4    = 0.0d+00

! Each process has ng4/nproc integrals according to rank
! Last process also has ng4%nproc remaining integrals
      call get_mpi_range(ng4,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng4

C Chop up calculation of OMG4 terms
      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng4/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )

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
                     do iec4=1,nebf
                     do jec4=1,nebf

                            imas=imas+1 ! imas is master_index
                            if(imas.ge.istart.and.imas.le.iend) then
                               Loopi=Loopi+1
                               loop_map(Loopi,1)=jec4
                               loop_map(Loopi,2)=iec4
                               loop_map(Loopi,3)=jec3
                               loop_map(Loopi,4)=iec3
                               loop_map(Loopi,5)=jec2
                               loop_map(Loopi,6)=iec2
                               loop_map(Loopi,7)=jec1
                               loop_map(Loopi,8)=iec1
                               loop_map(Loopi,9)=jp
                               loop_map(Loopi,10)=ip

! Save coordinates of first value for array index shifts
                     if ((ichunk.eq.1).and.(Loopi.eq.1)) then
                        call index_GAM_4PK2(nebf,npbf,ip,jp,iec1,jec1,
     x                           iec2,jec2,iec3,jec3,iec4,jec4,arrstart)
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
         end do
         end do

         call XCHF_OMG4_thread_MPI(istart,iend,ng4_seg,ng4loc,
     x                             nebf,npbf,
     x                             loop_map,arrstart,
     x                             DE,DP,
     x                             GM4,
     x                             XFE,XFP,E_OMG4)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrices
      call MPI_ALLREDUCE(MPI_IN_PLACE,XFE,nebf*nebf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,XFP,npbf*npbf,
     x                   MPI_DOUBLE_PRECISION,MPI_SUM,
     x                   MPI_COMM_WORLD,ierr)

      call add2fock(nebf,XFE,FE)
      call add2fock(npbf,XFP,FP)

      return
      end
C=======================================================================
      subroutine RXCHF_OMG4_thread_MPI(istart,iend,ng4_seg,ng4,
     x                                 nebf,npbf,
     x                                 loop_map,arrstart,
     x                                 DAE,DBE,DP,
     x                                 GM4,
     x                                 XFAE,XFBE,XFP,E_OMG4)
C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng4_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng4
      integer           arrstart
      integer           loop_map(ng4_seg,10)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM4(ng4)  ! GAM4   integrals

C Output variables
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  E_OMG4

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1  !
      integer           iec2,jec2  ! Contracted elec basis function indices
      integer           iec3,jec3  !
      integer           iec4,jec4  !
      integer           imap,ia
      double precision  val
      double precision  wtime
      double precision  three
      three=3.0d+00

C GM4ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron 1 indicies
C    iec3,jec3 : special electron 2 indicies
C    iec4,jec4 : special electron 3 indicies

!$omp parallel 
!$ompx shared(three)
!$ompx shared(loop_map,arrstart)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng4_seg)
!$ompx shared(ng4)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM4)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(iec4,jec4)
!$ompx private(ia)
!$ompx private(val)
!$ompx reduction(+:XFAE)
!$ompx reduction(+:XFBE)
!$ompx reduction(+:XFP)
!$ompx reduction(+:E_OMG4)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec4=loop_map(imap,1)
         iec4=loop_map(imap,2)
         jec3=loop_map(imap,3)
         iec3=loop_map(imap,4)
         jec2=loop_map(imap,5)
         iec2=loop_map(imap,6)
         jec1=loop_map(imap,7)
         iec1=loop_map(imap,8)
         jp =loop_map(imap,9)
         ip =loop_map(imap,10)

         call index_GAM_4PK2(nebf,npbf,ip,jp,iec1,jec1,
     x                       iec2,jec2,iec3,jec3,iec4,jec4,ia)

! GM4(ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4) contributes to
! F^e(ie1,je1),F^q(ie2,je2),F^p(ip,jp)
         val=GM4(ia-arrstart+1)

         XFAE(iec1,jec1)=XFAE(iec1,jec1)+DP(ip,jp)*DBE(iec2,jec2)*
     x                   DBE(iec3,jec3)*DBE(iec4,jec4)*val
         XFP(ip,jp)=XFP(ip,jp)+DAE(iec1,jec1)*DBE(iec2,jec2)*
     x              DBE(iec3,jec3)*DBE(iec4,jec4)*val
         XFBE(iec2,jec2)=XFBE(iec2,jec2)+three*DAE(iec1,jec1)*
     x                   DBE(iec3,jec3)*DBE(iec4,jec4)*DP(ip,jp)*val
         E_OMG4=E_OMG4+DP(ip,jp)*DAE(iec1,jec1)*
     x             DBE(iec2,jec2)*DBE(iec3,jec3)*DBE(iec4,jec4)*val

      end do
!$omp end do
!$omp end parallel      

      return
      end
C=======================================================================
      subroutine XCHF_OMG4_thread_MPI(istart,iend,ng4_seg,ng4,
     x                                nebf,npbf,
     x                                loop_map,arrstart,
     x                                DE,DP,
     x                                GM4,
     x                                XFE,XFP,E_OMG4)
C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng4_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng4
      integer           arrstart
      integer           loop_map(ng4_seg,10)
      double precision  DE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM4(ng4)  ! GAM4   integrals

C Output variables
      double precision  XFE(nebf,nebf)
      double precision  XFP(npbf,npbf)
      double precision  E_OMG4

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1  !
      integer           iec2,jec2  ! Contracted elec basis function indices
      integer           iec3,jec3  !
      integer           iec4,jec4  !
      integer           imap,ia
      double precision  val
      double precision  wtime
      double precision  four
      four=4.0d+00

C GM4ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron 1 indicies
C    iec3,jec3 : special electron 2 indicies
C    iec4,jec4 : special electron 3 indicies

!$omp parallel 
!$ompx shared(four)
!$ompx shared(loop_map,arrstart)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng4_seg)
!$ompx shared(ng4)
!$ompx shared(DE,DP)
!$ompx shared(GM4)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(iec4,jec4)
!$ompx private(ia)
!$ompx private(val)
!$ompx reduction(+:XFE)
!$ompx reduction(+:XFP)
!$ompx reduction(+:E_OMG4)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec4=loop_map(imap,1)
         iec4=loop_map(imap,2)
         jec3=loop_map(imap,3)
         iec3=loop_map(imap,4)
         jec2=loop_map(imap,5)
         iec2=loop_map(imap,6)
         jec1=loop_map(imap,7)
         iec1=loop_map(imap,8)
         jp =loop_map(imap,9)
         ip =loop_map(imap,10)

         call index_GAM_4PK2(nebf,npbf,ip,jp,iec1,jec1,
     x                       iec2,jec2,iec3,jec3,iec4,jec4,ia)

! GM4(ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4) contributes to
! F^e(ie1,je1),F^p(ip,jp)
         val=GM4(ia-arrstart+1)

         XFE(iec1,jec1)=XFE(iec1,jec1)+DP(ip,jp)*DE(iec2,jec2)*
     x                   DE(iec3,jec3)*DE(iec4,jec4)*four*val
         XFP(ip,jp)=XFP(ip,jp)+DE(iec1,jec1)*DE(iec2,jec2)*
     x              DE(iec3,jec3)*DE(iec4,jec4)*val
         E_OMG4=E_OMG4+DP(ip,jp)*DE(iec1,jec1)*
     x             DE(iec2,jec2)*DE(iec3,jec3)*DE(iec4,jec4)*val

      end do
!$omp end do
!$omp end parallel      

      return
      end

