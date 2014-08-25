!======================================================================
      subroutine RXCHF_GAM4_MPI(nproc,rank,
     x                          Nchunks,ne,np,ngee,ng2,ng4,
     x                          ng2loc,ng4loc,
     x                          GM2s)

! Calculates INT_GAM4 integrals (split onto MPI procs)
!  - requires each proc to store all XCHF_GAM2s integrals in memory
!======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ngee,ne,np
      integer ng4,ng2       ! Total number of integrals
      integer ng2loc,ng4loc ! Number of integrals on this proc
      double precision GM2s(ng2loc) ! GAM2s integrals on this MPI proc
      
! Local Variables
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,i
      integer ipstart,jpstart
      integer ie1start,je1start,ie2start,je2start
      integer ie3start,je3start,ie4start,je4start
      double precision GAM_ee(ngee)
      double precision GM4(ng4loc)  ! INT OMG4 integrals (symm)

      integer nproc,rank
      integer mpistart,mpiend,arrstart

      integer nints,int_count

      integer(kind=4),  allocatable :: protinds(:)
      integer(kind=8),  allocatable :: elecinds(:)
      double precision, allocatable :: nonzeroints(:)

      double precision, allocatable :: TGM2s(:),TGM4(:)
#if MPI32
      integer*4 ierr
      integer*4 ng2loc_,ng4loc_
      integer*4 ng2locarr(nproc),ng4locarr(nproc),displarr(nproc)
#else
      integer ierr
      integer ng2loc_,ng4loc_
      integer ng2locarr(nproc),ng4locarr(nproc),displarr(nproc)
#endif

      double precision zero
      parameter(zero=0.0d+00)

      character*4  istring              ! File I/O variables

      double precision wtime,wtime1,wtime2

! Have each process calculate ng4/nproc integrals according to rank
! Have last process calculate ng4%nproc remaining integrals
      call get_mpi_range(ng4,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng4

      if (rank.eq.0) then
       write(*,1000) ng4,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      GM4=0.0d+00

      if (rank.eq.0) call read_GAM_ee(ne,ngee,GAM_ee) 
      call MPI_BCAST(GAM_ee,ngee,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)

! Form global GAM2s on each process
      if(allocated(TGM2s)) deallocate(TGM2s)
      allocate(TGM2s(ng2))
      TGM2s=zero

! Get number of elements calculated by each proc
#if MPI32
      ng2loc_=int(ng2loc,kind=4)
      call MPI_ALLGATHER(ng2loc_,1,MPI_INTEGER,
     x                   ng2locarr(1),1,MPI_INTEGER,
     x                   MPI_COMM_WORLD,ierr)
#else
      ng2loc_=ng2loc
      call MPI_ALLGATHER(ng2loc_,1,MPI_INTEGER8,
     x                   ng2locarr(1),1,MPI_INTEGER8,
     x                   MPI_COMM_WORLD,ierr)
#endif

! Get displacements for array storage
      displarr(1)=0
      do i=2,nproc
        displarr(i)=displarr(i-1)+ng2locarr(i-1)
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Form global GM2s on each proc
      call MPI_ALLGATHERV(GM2s(1),ng2loc_,MPI_DOUBLE_PRECISION,
     x                    TGM2s(1),ng2locarr,displarr,
     x                    MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

C Variables for file I/O
      write(istring,'(I4.4)') rank

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2000) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      wtime = MPI_WTIME()

      nints=0

      do ip=1,np
      do jp=1,np
        do ie1=1,ne
        do je1=1,ne
          do ie2=1,ne
          do je2=1,ne
            do ie3=1,ne
            do je3=1,ne
              do ie4=1,ne
              do je4=1,ne

                nints=nints+1
                if(nints.eq.mpistart) then
                 ipstart=ip
                 jpstart=jp
                 ie1start=ie1
                 je1start=je1
                 ie2start=ie2
                 je2start=je2
                 ie3start=ie3
                 je3start=je3
                 ie4start=ie4
                 je4start=je4
                 GO TO 100
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

  100 CONTINUE

! Loop through to calc both types of ints and get # of each that are nonzero
      call RXCHF_GAM4_thread_MPI(ne,np,ngee,ng2,ng4loc,
     x                           ipstart,jpstart,
     x                           ie1start,je1start,
     x                           ie2start,je2start,
     x                           ie3start,je3start,
     x                           ie4start,je4start,
     x                           GAM_ee,TGM2s,GM4,
     x                           int_count)

      wtime1 = MPI_WTIME() - wtime
      write(*,2011) rank,wtime1

! Loop again to store nontrivial INT_GAM4 and indices
      if(allocated(protinds)) deallocate(protinds)
      allocate(protinds(int_count))
      if(allocated(elecinds)) deallocate(elecinds)
      allocate(elecinds(int_count))
      if(allocated(nonzeroints)) deallocate(nonzeroints)
      allocate(nonzeroints(int_count))
      protinds=0
      elecinds=0
      nonzeroints=0.0d+00

      call GAM4_nonzero_MPI(ne,np,ng4loc,int_count,
     x                      ipstart,jpstart,
     x                      ie1start,je1start,
     x                      ie2start,je2start,
     x                      ie3start,je3start,
     x                      ie4start,je4start,
     x                      GM4,protinds,elecinds,nonzeroints)

      call writebitint(nproc,rank,int_count,17,
     x                 "INT_GAM4-"//istring//".ufm",
     x                 protinds,elecinds,nonzeroints)

      if(allocated(nonzeroints)) deallocate(nonzeroints)
      if(allocated(elecinds)) deallocate(elecinds)
      if(allocated(protinds)) deallocate(protinds)

      wtime2 = MPI_WTIME() - wtime

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2010)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      write(*,2011) rank,wtime2
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng4loc:",mpistart,mpiend,ng4loc
!      do i=1,ng4loc
!       write(*,9001) GM4(i)
!      end do

      if(allocated(TGM2s)) deallocate(TGM2s)

! Construct global array on master process for testing
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if(allocated(TGM4)) deallocate(TGM4)
!      if (rank.eq.0) then
!       allocate(TGM4(ng4))
!      else
!       allocate(TGM4(1))
!      end if
!      TGM4=zero
!
!      ng4loc4=int(ng4loc,kind=4)
!
!! Get number of elements calculated by each proc
!      call MPI_GATHER(ng4loc4,1,MPI_INTEGER,
!     x                ng4locarr(1),1,MPI_INTEGER,
!     x                0,MPI_COMM_WORLD,ierr)
!
!! Get displacements for array storage
!      if (rank.eq.0) then
!        displarr(1)=0
!        do i=2,nproc
!          displarr(i)=displarr(i-1)+ng4locarr(i-1)
!        end do
!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM4 on root
!      call MPI_GATHERV(GM4(1),ng4loc,MPI_DOUBLE_PRECISION,
!     x                 TGM4(1),ng4locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!      if (rank.eq.0) then
!!       write(*,*) "concatenated ng4"
!!       do i=1,ng4
!!        write(*,9001) TGM4(i)
!!       end do
!!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if (rank.eq.0) then
!       open(unit=20,file="INT_GAM4.ufm",form="unformatted")
!       write(20) TGM4
!       close(20)
!       write(*,*) "INT_GAM4 written to disk"
!      end if
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!      if(allocated(TGM4)) deallocate(TGM4)

 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 5-PARTICLE INTEGRALS        |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 5-PARTICLE INTEGRALS: ',1X,I15/
     x        8X,'  NUMBER OF BLOCKS (USER DEFINED): ',1X,I15/
     x        8X,'                          ',/,
     x        8X,'  COMPUTATIONAL RESOURCES:',/,
     x        8X,'  ------------------------',/)
                       
 1500 FORMAT( 8X,'      MPI PROCESSES:',1X,I4/
     x        8X,'        OMP THREADS:',1X,I4/)

 2000 FORMAT(/8X,'  ALL INTEGRAL CALCULATION TIMINGS:',/,
     x        8X,' -------------------------------------')

 2001 FORMAT( 8X,' PROCESS ',1X,I4,1X,' BLOCK ',1X,I4,1X,F10.2)

 2010 FORMAT(/8X,'  INTEGRAL CALC+SYMM TIMINGS:',/,
     x        8X,'  ---------------------------')

 2011 FORMAT( 8X,'   PROCESS ',1X,I4,1X,F10.2)

 9001 FORMAT(1X,1(F20.10))

      return
      end
!======================================================================
      subroutine XCHF_GAM4_MPI(nproc,rank,
     x                         Nchunks,ne,np,ngee,ng2,ng4,
     x                         ng2loc,ng4loc,
     x                         GM2s)

! Calculates INT_GAM4 and XCHF_GAM4 integrals (split onto MPI procs)
!  - requires each proc to store all XCHF_GAM2s integrals in memory
!======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ngee,ne,np
      integer ng4,ng2       ! Total number of integrals
      integer ng2loc,ng4loc ! Number of integrals on this proc
      double precision GM2s(ng2loc) ! GAM2s integrals on this MPI proc
      
! Local Variables
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,i
      integer ipstart,jpstart
      integer ie1start,je1start,ie2start,je2start
      integer ie3start,je3start,ie4start,je4start
      double precision GAM_ee(ngee)
      double precision GM4_1(ng4loc)  ! INT  OMG4 integrals (symm)
      double precision GM4_2(ng4loc)  ! XCHF OMG4 integrals (symm)

      integer nproc,rank
      integer mpistart,mpiend,arrstart

      integer nints,int_count,xchf_count

      integer(kind=4),  allocatable :: protinds(:)
      integer(kind=8),  allocatable :: elecinds(:)
      double precision, allocatable :: nonzeroints(:)

      double precision, allocatable :: TGM2s(:),TGM4(:)
#if MPI32
      integer*4 ierr
      integer*4 ng2loc_,ng4loc_
      integer*4 ng2locarr(nproc),ng4locarr(nproc),displarr(nproc)
#else
      integer ierr
      integer ng2loc_,ng4loc_
      integer ng2locarr(nproc),ng4locarr(nproc),displarr(nproc)
#endif

      double precision zero
      parameter(zero=0.0d+00)

      character*4  istring              ! File I/O variables

      double precision wtime,wtime1,wtime2

! Have each process calculate ng4/nproc integrals according to rank
! Have last process calculate ng4%nproc remaining integrals
      call get_mpi_range(ng4,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng4

      if (rank.eq.0) then
       write(*,1000) ng4,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      GM4_1=0.0d+00
      GM4_2=0.0d+00

      if (rank.eq.0) call read_GAM_ee(ne,ngee,GAM_ee) 
      call MPI_BCAST(GAM_ee,ngee,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)

! Form global GAM2s on each process
      if(allocated(TGM2s)) deallocate(TGM2s)
      allocate(TGM2s(ng2))
      TGM2s=zero

! Get number of elements calculated by each proc
#if MPI32
      ng2loc_=int(ng2loc,kind=4)
      call MPI_ALLGATHER(ng2loc_,1,MPI_INTEGER,
     x                   ng2locarr(1),1,MPI_INTEGER,
     x                   MPI_COMM_WORLD,ierr)
#else
      ng2loc_=ng2loc
      call MPI_ALLGATHER(ng2loc_,1,MPI_INTEGER8,
     x                   ng2locarr(1),1,MPI_INTEGER8,
     x                   MPI_COMM_WORLD,ierr)
#endif

! Get displacements for array storage
      displarr(1)=0
      do i=2,nproc
        displarr(i)=displarr(i-1)+ng2locarr(i-1)
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Form global GM2s on each proc
      call MPI_ALLGATHERV(GM2s(1),ng2loc_,MPI_DOUBLE_PRECISION,
     x                    TGM2s(1),ng2locarr,displarr,
     x                    MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

C Variables for file I/O
      write(istring,'(I4.4)') rank

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2000) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      wtime = MPI_WTIME()

      nints=0

      do ip=1,np
      do jp=1,np
        do ie1=1,ne
        do je1=1,ne
          do ie2=1,ne
          do je2=1,ne
            do ie3=1,ne
            do je3=1,ne
              do ie4=1,ne
              do je4=1,ne

                nints=nints+1
                if(nints.eq.mpistart) then
                 ipstart=ip
                 jpstart=jp
                 ie1start=ie1
                 je1start=je1
                 ie2start=ie2
                 je2start=je2
                 ie3start=ie3
                 je3start=je3
                 ie4start=ie4
                 je4start=je4
                 GO TO 200
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

  200 CONTINUE

! Loop through to calc both types of ints and get # of each that are nonzero
      call XCHF_GAM4_thread_MPI(ne,np,ngee,ng2,ng4loc,
     x                          ipstart,jpstart,
     x                          ie1start,je1start,
     x                          ie2start,je2start,
     x                          ie3start,je3start,
     x                          ie4start,je4start,
     x                          GAM_ee,TGM2s,GM4_1,GM4_2,
     x                          int_count,xchf_count)

      wtime1 = MPI_WTIME() - wtime
      write(*,2011) rank,wtime1

! Loop again to store nontrivial INT_GAM4 and indices
      if(allocated(protinds)) deallocate(protinds)
      allocate(protinds(int_count))
      if(allocated(elecinds)) deallocate(elecinds)
      allocate(elecinds(int_count))
      if(allocated(nonzeroints)) deallocate(nonzeroints)
      allocate(nonzeroints(int_count))
      protinds=0
      elecinds=0
      nonzeroints=0.0d+00

      call GAM4_nonzero_MPI(ne,np,ng4loc,int_count,
     x                      ipstart,jpstart,
     x                      ie1start,je1start,
     x                      ie2start,je2start,
     x                      ie3start,je3start,
     x                      ie4start,je4start,
     x                      GM4_1,protinds,elecinds,nonzeroints)

      call writebitint(nproc,rank,int_count,17,
     x                 "INT_GAM4-"//istring//".ufm",
     x                 protinds,elecinds,nonzeroints)

      if(allocated(nonzeroints)) deallocate(nonzeroints)
      if(allocated(elecinds)) deallocate(elecinds)
      if(allocated(protinds)) deallocate(protinds)

! Loop again to store nontrivial XCHF_GAM4 and indices
      if(allocated(protinds)) deallocate(protinds)
      allocate(protinds(xchf_count))
      if(allocated(elecinds)) deallocate(elecinds)
      allocate(elecinds(xchf_count))
      if(allocated(nonzeroints)) deallocate(nonzeroints)
      allocate(nonzeroints(xchf_count))
      protinds=0
      elecinds=0
      nonzeroints=0.0d+00

      call GAM4_nonzero_MPI(ne,np,ng4loc,xchf_count,
     x                      ipstart,jpstart,
     x                      ie1start,je1start,
     x                      ie2start,je2start,
     x                      ie3start,je3start,
     x                      ie4start,je4start,
     x                      GM4_2,protinds,elecinds,nonzeroints)

      call writebitint(nproc,rank,xchf_count,18,
     x                 "XCHF_GAM4-"//istring//".ufm",
     x                 protinds,elecinds,nonzeroints)

      if(allocated(nonzeroints)) deallocate(nonzeroints)
      if(allocated(elecinds)) deallocate(elecinds)
      if(allocated(protinds)) deallocate(protinds)

      wtime2 = MPI_WTIME() - wtime

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2010)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      write(*,2011) rank,wtime2
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng4loc:",mpistart,mpiend,ng4loc
!      do i=1,ng4loc
!       write(*,9001) GM4(i)
!      end do

      if(allocated(TGM2s)) deallocate(TGM2s)

! Construct global array on master process for testing
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if(allocated(TGM4)) deallocate(TGM4)
!      if (rank.eq.0) then
!       allocate(TGM4(ng4))
!      else
!       allocate(TGM4(1))
!      end if
!      TGM4=zero
!
!      ng4loc4=int(ng4loc,kind=4)
!
!! Get number of elements calculated by each proc
!      call MPI_GATHER(ng4loc4,1,MPI_INTEGER,
!     x                ng4locarr(1),1,MPI_INTEGER,
!     x                0,MPI_COMM_WORLD,ierr)
!
!! Get displacements for array storage
!      if (rank.eq.0) then
!        displarr(1)=0
!        do i=2,nproc
!          displarr(i)=displarr(i-1)+ng4locarr(i-1)
!        end do
!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM4 on root
!      call MPI_GATHERV(GM4(1),ng4loc,MPI_DOUBLE_PRECISION,
!     x                 TGM4(1),ng4locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!      if (rank.eq.0) then
!!       write(*,*) "concatenated ng4"
!!       do i=1,ng4
!!        write(*,9001) TGM4(i)
!!       end do
!!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if (rank.eq.0) then
!       open(unit=20,file="XCHF_GAM4.ufm",form="unformatted")
!       write(20) TGM4
!       close(20)
!       write(*,*) "XCHF_GAM4 written to disk"
!      end if
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!      if(allocated(TGM4)) deallocate(TGM4)

 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 5-PARTICLE INTEGRALS        |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 5-PARTICLE INTEGRALS: ',1X,I15/
     x        8X,'  NUMBER OF BLOCKS (USER DEFINED): ',1X,I15/
     x        8X,'                          ',/,
     x        8X,'  COMPUTATIONAL RESOURCES:',/,
     x        8X,'  ------------------------',/)
                       
 1500 FORMAT( 8X,'      MPI PROCESSES:',1X,I4/
     x        8X,'        OMP THREADS:',1X,I4/)

 2000 FORMAT(/8X,'  All INTEGRAL CALCULATION TIMINGS:',/,
     x        8X,' -------------------------------------')

 2001 FORMAT( 8X,' PROCESS ',1X,I4,1X,' BLOCK ',1X,I4,1X,F10.2)

 2010 FORMAT(/8X,'  INTEGRAL CALC+SYMM TIMINGS:',/,
     x        8X,'  ---------------------------')

 2011 FORMAT( 8X,'   PROCESS ',1X,I4,1X,F10.2)

 9001 FORMAT(1X,1(F20.10))

      return
      end
!======================================================================
      subroutine RXCHF_GAM4_thread_MPI(ne,np,ngee,ng2,ng4,
     x                                ipstart,jpstart,
     x                                ie1start,je1start,
     x                                ie2start,je2start,
     x                                ie3start,je3start,
     x                                ie4start,je4start,
     x                                GAM_ee,GAM_2s,GAM_4,
     x                                int_count)
!
!======================================================================
      implicit none
      include 'omp_lib.h'

! Input variables
      integer ne    ! Number of contracted electronic basis functions
      integer np    ! Number of nuclear basis functions
      integer ngee  ! Number of contracted 2-electron integrals
      integer ng2   ! Number of contracted 3-particle integrals
      integer ng4   ! Number of integrals calc by MPI process
      integer ipstart  !
      integer jpstart  !
      integer ie1start !
      integer ie2start !
      integer ie3start ! Initial indices for this MPI process
      integer ie4start !
      integer je1start !
      integer je2start !
      integer je3start !
      integer je4start !
      double precision GAM_ee(ngee)   ! Array storage of 2e-integrals
      double precision GAM_2s(ng2)    ! Array storage of 3-part overlaps

! Output variables
      double precision GAM_4(ng4)      ! Array storage of 5-particle ints
      integer int_count               ! # of nontrivial INT_GAM4 integrals

! Local variables
      integer i
      integer ip
      integer jp
      integer ie1
      integer ie2
      integer ie3
      integer ie4
      integer je1
      integer je2
      integer je3
      integer je4

      double precision ans
      double precision tol

      tol=1.0d-15

      int_count=0

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(ne,np)
!$ompx shared(ngee)
!$ompx shared(ng2)
!$ompx shared(ng4)
!$ompx shared(ipstart,jpstart) 
!$ompx shared(ie1start,je1start) 
!$ompx shared(ie2start,je2start) 
!$ompx shared(ie3start,je3start) 
!$ompx shared(ie4start,je4start) 
!$ompx shared(tol)
!$ompx shared(GAM_2s)
!$ompx shared(GAM_ee)
!$ompx shared(GAM_4)
!$ompx private(i)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(ie3,je3) 
!$ompx private(ie4,je4) 
!$ompx private(ans)
!$ompx reduction(+:int_count)

!$omp do SCHEDULE(RUNTIME)
      do i=1,ng4

         call getgam4inds(i,np,ne,ipstart,jpstart,
     x                    ie1start,je1start,ie2start,je2start,
     x                    ie3start,je3start,ie4start,je4start,
     x                    ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4)

         call RXCHFmult_symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

         GAM_4(i)=ans
         if(abs(ans).gt.tol) int_count=int_count+1

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

      return
      end 
!======================================================================
      subroutine XCHF_GAM4_thread_MPI(ne,np,ngee,ng2,ng4,
     x                                ipstart,jpstart,
     x                                ie1start,je1start,
     x                                ie2start,je2start,
     x                                ie3start,je3start,
     x                                ie4start,je4start,
     x                                GAM_ee,GAM_2s,GAM4_1,GAM4_2,
     x                                int_count,xchf_count)
!
!======================================================================
      implicit none
      include 'omp_lib.h'

! Input variables
      integer ne       ! Number of contracted electronic basis functions
      integer np       ! Number of nuclear basis functions
      integer ngee     ! Number of contracted 2-electron integrals
      integer ng2      ! Number of contracted 3-particle integrals
      integer ng4      ! Number of integrals calc by MPI process
      integer ipstart  !
      integer jpstart  !
      integer ie1start !
      integer ie2start !
      integer ie3start ! Initial indices for this MPI process
      integer ie4start !
      integer je1start !
      integer je2start !
      integer je3start !
      integer je4start !
      double precision GAM_ee(ngee)   ! Array storage of 2e-integrals
      double precision GAM_2s(ng2)    ! Array storage of 3-part overlaps

! Output variables
      double precision GAM4_1(ng4)    ! Array storage of 5-particle ints (INT)
      double precision GAM4_2(ng4)    ! Array storage of 5-particle ints (XCHF)
      integer int_count               ! # of nontrivial INT_GAM4 integrals
      integer xchf_count              ! # of nontrivial XCHF_GAM4 integrals

! Local variables
      integer i
      integer ip
      integer jp
      integer ie1
      integer ie2
      integer ie3
      integer ie4
      integer je1
      integer je2
      integer je3
      integer je4

      double precision ans
      double precision four
      double precision tol

      four=4.0d+00
      tol=1.0d-15

      int_count=0
      xchf_count=0

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(ne,np)
!$ompx shared(ngee)
!$ompx shared(ng2)
!$ompx shared(ng4)
!$ompx shared(ipstart,jpstart) 
!$ompx shared(ie1start,je1start) 
!$ompx shared(ie2start,je2start) 
!$ompx shared(ie3start,je3start) 
!$ompx shared(ie4start,je4start) 
!$ompx shared(four)
!$ompx shared(tol)
!$ompx shared(GAM_2s)
!$ompx shared(GAM_ee)
!$ompx shared(GAM4_1)
!$ompx shared(GAM4_2)
!$ompx private(i) 
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(ie3,je3) 
!$ompx private(ie4,je4) 
!$ompx private(ans)
!$ompx reduction(+:int_count)
!$ompx reduction(+:xchf_count)

!$omp do SCHEDULE(RUNTIME)
      do i=1,ng4

         call getgam4inds(i,np,ne,ipstart,jpstart,
     x                    ie1start,je1start,ie2start,je2start,
     x                    ie3start,je3start,ie4start,je4start,
     x                    ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4)

         call RXCHFmult_symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

         GAM4_1(i)=ans
         if(abs(ans).gt.tol) int_count=int_count+1

         call symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

         ans=ans/(four*four)
         GAM4_2(i)=ans
         if(abs(ans).gt.tol) xchf_count=xchf_count+1

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

      return
      end 
!======================================================================
      subroutine GAM4_nonzero_MPI(ne,np,ng4,ng4_nz,
     x                            ipstart,jpstart,
     x                            ie1start,je1start,
     x                            ie2start,je2start,
     x                            ie3start,je3start,
     x                            ie4start,je4start,
     x                            GAM4,pinds,einds,ints_nz)
!
!======================================================================
      implicit none

! Input variables
      integer ne       ! Number of contracted electronic basis functions
      integer np       ! Number of nuclear basis functions
      integer ng4      ! Tot number of integrals calc by MPI process
      integer ng4_nz   ! Number of nonzero integrals calc by MPI process
      integer ipstart  !
      integer jpstart  !
      integer ie1start !
      integer ie2start !
      integer ie3start ! Initial indices for this MPI process
      integer ie4start !
      integer je1start !
      integer je2start !
      integer je3start !
      integer je4start !
      double precision GAM4(ng4)   ! All GAM4 integrals calc by this MPI proc

! Output variables
      double precision ints_nz(ng4_nz) ! Nonzero integrals
      double precision pinds(ng4_nz)   ! Packed proton indices
      double precision einds(ng4_nz)   ! Packed electron indices

! Local variables
      integer i
      integer ip
      integer jp
      integer ie1
      integer ie2
      integer ie3
      integer ie4
      integer je1
      integer je2
      integer je3
      integer je4

      integer(kind=4)  packp
      integer(kind=8)  packe
      integer          currind

      double precision tol

      tol=1.0d-15

      currind=0

      do i=1,ng4

         if(abs(GAM4(i)).gt.tol) then

          currind=currind+1
          call getgam4inds(i,np,ne,ipstart,jpstart,
     x                     ie1start,je1start,ie2start,je2start,
     x                     ie3start,je3start,ie4start,je4start,
     x                     ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4)

          call bitpack(ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,
     x                 packp,packe)

          pinds(currind)=packp
          einds(currind)=packe
          ints_nz(currind)=GAM4(i)

         end if

      end do

      return
      end 


      subroutine getgam4inds(currind,npbf,nebf,
     x                       ips,jps,i1s,j1s,i2s,j2s,i3s,j3s,i4s,j4s,
     x                       ip,jp,i1,j1,i2,j2,i3,j3,i4,j4)
      implicit none

! Input variables
      integer currind
      integer npbf,nebf
      integer ips,jps
      integer i1s,j1s
      integer i2s,j2s
      integer i3s,j3s
      integer i4s,j4s

! Output variables
      integer ip,jp
      integer i1,j1
      integer i2,j2
      integer i3,j3
      integer i4,j4

! Local variables
      integer iwork

      iwork=currind
      ip=ips
      jp=jps
      i1=i1s
      j1=j1s
      i2=i2s
      j2=j2s
      i3=i3s
      j3=j3s
      i4=i4s
      j4=j4s

      do while (iwork.gt.0)

        if (iwork.le.(nebf-j4)) then
         j4=j4+iwork
         iwork=0
        else
         j4=1
         iwork=iwork-(nebf-j4)

         if (i4.lt.nebf) then
          i4=i4+1
         else
          i4=1

          if (j3.lt.nebf) then
           j3=j3+1
          else
           j3=1

           if (i3.lt.nebf) then
            i3=i3+1
           else
            i3=1

            if (j2.lt.nebf) then
             j2=j2+1
            else
             j2=1

             if (i2.lt.nebf) then
              i2=i2+1
             else
              i2=1

              if (j1.lt.nebf) then
               j1=j1+1
              else
               j1=1

               if (i1.lt.nebf) then
                i1=i1+1
               else
                i1=1

                if (jp.lt.npbf) then
                 jp=jp+1
                else
                 jp=1

                 if (ip.lt.npbf) then
                  ip=ip+1
                 else
                  write(*,*) "ERROR IN UPDATEINDS: index overflow"
                 end if

                end if

               end if

              end if

             end if

            end if

           end if

          end if

         end if

        end if

      end do

      return
      end

