!======================================================================
      subroutine RXCHF_GAM4_MPI(nproc,rank,
     x                          Nchunks,ne,np,ngee,ng2,ng4,
     x                          ng2loc,ng4loc,
     x                          GM2s,GM4)

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
      double precision GM4(ng4loc)  ! INT OMG4 integrals (symm)
      double precision GM2s(ng2loc) ! GAM2s integrals on this MPI proc
      
! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,i
      integer,allocatable :: loop_map(:,:)
      double precision,allocatable :: GAM_ee(:)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

      double precision, allocatable :: TGM2s(:),TGM4(:)
      integer*4 ng2loc4,ng4loc4
      integer*4 ng2locarr(nproc),ng4locarr(nproc),displarr(nproc)

      double precision zero
      parameter(zero=0.0d+00)

      double precision wtime,wtime2

! Have each process calculate ng4/nproc integrals according to rank
! Have last process calculate ng4%nproc remaining integrals
      call get_mpi_range(ng4,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng4
      write(*,*) "rank,mpistart,mpiend,ng4loc:",
     x           rank,mpistart,mpiend,ng4loc

      if (rank.eq.0) then
       write(*,1000) ng4,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      GM4=0.0d+00

!----READ-GAM_ee--INTO-MEMORY-------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      allocate( GAM_ee(ngee),stat=istat )
      if (rank.eq.0) call read_GAM_ee(ne,ngee,GAM_ee) 
      call MPI_BCAST(GAM_ee,ngee,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
!----READ-GAM_ee-INTO-MEMORY-------------------------------)

! Form global GAM2s on each process
      if(allocated(TGM2s)) deallocate(TGM2s)
      allocate(TGM2s(ng2))
      TGM2s=zero

      ng2loc4=int(ng2loc,kind=4)

! Get number of elements calculated by each proc
      call MPI_ALLGATHER(ng2loc4,1,MPI_INTEGER,
     x                   ng2locarr(1),1,MPI_INTEGER,
     x                   MPI_COMM_WORLD,ierr)

! Get displacements for array storage
      displarr(1)=0
      do i=2,nproc
        displarr(i)=displarr(i-1)+ng2locarr(i-1)
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Form global GM2s on each proc
      call MPI_ALLGATHERV(GM2s(1),ng2loc,MPI_DOUBLE_PRECISION,
     x                    TGM2s(1),ng2locarr,displarr,
     x                    MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(
      wtime = MPI_WTIME()

      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng4/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )

! Nested loop compression for this chunk:
         Loopi=0
         imas=0
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

                         imas=imas+1 ! imas is master_index
                         if(imas.ge.istart.and.imas.le.iend) then
                            Loopi=Loopi+1
                            loop_map(Loopi,1)=je4
                            loop_map(Loopi,2)=ie4
                            loop_map(Loopi,3)=je3
                            loop_map(Loopi,4)=ie3
                            loop_map(Loopi,5)=je2
                            loop_map(Loopi,6)=ie2
                            loop_map(Loopi,7)=je1
                            loop_map(Loopi,8)=ie1
                            loop_map(Loopi,9)=jp
                            loop_map(Loopi,10)=ip

! Save coordinates of first value for array transfer later
                            if ((ichunk.eq.1).and.(Loopi.eq.1)) then
                               call index_GAM_4PK2(ne,np,
     x                                             ip,jp,
     x                                             ie1,je1,
     x                                             ie2,je2,
     x                                             ie3,je3,
     x                                             ie4,je4,arrstart)
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

         call RXCHF_GAM4_thread_MPI(ne,np,ngee,ng2,ng4loc,ng4_seg,
     x                              istart,iend,loop_map,arrstart,
     x                              GAM_ee,TGM2s,GM4)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)
!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

      wtime2 = MPI_WTIME() - wtime
      write(*,2000)rank,wtime2

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng4loc:",mpistart,mpiend,ng4loc
!      do i=1,ng4loc
!       write(*,9001) GM4(i)
!      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(allocated(TGM2s)) deallocate(TGM2s)
      if(allocated(GAM_ee)) deallocate(GAM_ee)

! Construct global array on master process for testing
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(allocated(TGM4)) deallocate(TGM4)
      if (rank.eq.0) then
       allocate(TGM4(ng4))
      else
       allocate(TGM4(1))
      end if
      TGM4=zero

      ng4loc4=int(ng4loc,kind=4)

! Get number of elements calculated by each proc
      call MPI_GATHER(ng4loc4,1,MPI_INTEGER,
     x                ng4locarr(1),1,MPI_INTEGER,
     x                0,MPI_COMM_WORLD,ierr)

! Get displacements for array storage
      if (rank.eq.0) then
        displarr(1)=0
        do i=2,nproc
          displarr(i)=displarr(i-1)+ng4locarr(i-1)
        end do
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Form global GM4 on root
      call MPI_GATHERV(GM4(1),ng4loc,MPI_DOUBLE_PRECISION,
     x                 TGM4(1),ng4locarr,displarr,
     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if (rank.eq.0) then
!       write(*,*) "concatenated ng4"
!       do i=1,ng4
!        write(*,9001) TGM4(i)
!       end do
!      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (rank.eq.0) then
       open(unit=20,file="INT_GAM4.ufm",form="unformatted")
       write(20) TGM4
       close(20)
       write(*,*) "INT_GAM4 written to disk"
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(allocated(TGM4)) deallocate(TGM4)

 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 5-PARTICLE INTEGRALS        |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 5-PARTICLE INTEGRALS: ',1X,I12/
     x        8X,'  NUMBER OF BLOCKS (USER DEFINED): ',1X,I12/
     x        8X,'                          ',/,
     x        8X,'  COMPUTATIONAL RESOURCES:',/,
     x        8X,'  ------------------------',/)
                       
 1500 FORMAT( 8X,'      MPI PROCESSES:',1X,I3/
     x        8X,'        OMP THREADS:',1X,I3/)

 2000 FORMAT(8X,'    TIME FOR PROCESS ',1X,I4,1X,F10.2)

 9001 FORMAT(1X,1(F20.10))

      return
      end
!======================================================================
      subroutine XCHF_GAM4_MPI(nproc,rank,
     x                         Nchunks,ne,np,ngee,ng2,ng4,
     x                         ng2loc,ng4loc,
     x                         GM2s,GM4)

! Calculates XCHF_GAM4 integrals (split onto MPI procs)
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
      double precision GM4(ng4loc)  ! XCHF OMG4 integrals (symm)
      double precision GM2s(ng2loc) ! GAM2s integrals on this MPI proc
      
! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,i
      integer,allocatable :: loop_map(:,:)
      double precision,allocatable :: GAM_ee(:)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

      double precision, allocatable :: TGM2s(:),TGM4(:)
      integer*4 ng2loc4,ng4loc4
      integer*4 ng2locarr(nproc),ng4locarr(nproc),displarr(nproc)

      double precision zero
      parameter(zero=0.0d+00)

      double precision wtime,wtime2

! Have each process calculate ng4/nproc integrals according to rank
! Have last process calculate ng4%nproc remaining integrals
      call get_mpi_range(ng4,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng4
      write(*,*) "rank,mpistart,mpiend,ng4loc:",
     x           rank,mpistart,mpiend,ng4loc

      if (rank.eq.0) then
       write(*,1000) ng4,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      GM4=0.0d+00

!----READ-GAM_ee--INTO-MEMORY-------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      allocate( GAM_ee(ngee),stat=istat )
      if (rank.eq.0) call read_GAM_ee(ne,ngee,GAM_ee) 
      call MPI_BCAST(GAM_ee,ngee,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
!----READ-GAM_ee-INTO-MEMORY-------------------------------)

! Form global GAM2s on each process
      if(allocated(TGM2s)) deallocate(TGM2s)
      allocate(TGM2s(ng2))
      TGM2s=zero

      ng2loc4=int(ng2loc,kind=4)

! Get number of elements calculated by each proc
      call MPI_ALLGATHER(ng2loc4,1,MPI_INTEGER,
     x                   ng2locarr(1),1,MPI_INTEGER,
     x                   MPI_COMM_WORLD,ierr)

! Get displacements for array storage
      displarr(1)=0
      do i=2,nproc
        displarr(i)=displarr(i-1)+ng2locarr(i-1)
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Form global GM2s on each proc
      call MPI_ALLGATHERV(GM2s(1),ng2loc,MPI_DOUBLE_PRECISION,
     x                    TGM2s(1),ng2locarr,displarr,
     x                    MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(
      wtime = MPI_WTIME()

      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng4/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )

! Nested loop compression for this chunk:
         Loopi=0
         imas=0
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

                         imas=imas+1 ! imas is master_index
                         if(imas.ge.istart.and.imas.le.iend) then
                            Loopi=Loopi+1
                            loop_map(Loopi,1)=je4
                            loop_map(Loopi,2)=ie4
                            loop_map(Loopi,3)=je3
                            loop_map(Loopi,4)=ie3
                            loop_map(Loopi,5)=je2
                            loop_map(Loopi,6)=ie2
                            loop_map(Loopi,7)=je1
                            loop_map(Loopi,8)=ie1
                            loop_map(Loopi,9)=jp
                            loop_map(Loopi,10)=ip

! Save coordinates of first value for array transfer later
                            if ((ichunk.eq.1).and.(Loopi.eq.1)) then
                               call index_GAM_4PK2(ne,np,
     x                                             ip,jp,
     x                                             ie1,je1,
     x                                             ie2,je2,
     x                                             ie3,je3,
     x                                             ie4,je4,arrstart)
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

         call XCHF_GAM4_thread_MPI(ne,np,ngee,ng2,ng4loc,ng4_seg,
     x                             istart,iend,loop_map,arrstart,
     x                             GAM_ee,TGM2s,GM4)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)
!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

      wtime2 = MPI_WTIME() - wtime
      write(*,2000)rank,wtime2

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng4loc:",mpistart,mpiend,ng4loc
!      do i=1,ng4loc
!       write(*,9001) GM4(i)
!      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(allocated(TGM2s)) deallocate(TGM2s)
      if(allocated(GAM_ee)) deallocate(GAM_ee)

! Construct global array on master process for testing
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(allocated(TGM4)) deallocate(TGM4)
      if (rank.eq.0) then
       allocate(TGM4(ng4))
      else
       allocate(TGM4(1))
      end if
      TGM4=zero

      ng4loc4=int(ng4loc,kind=4)

! Get number of elements calculated by each proc
      call MPI_GATHER(ng4loc4,1,MPI_INTEGER,
     x                ng4locarr(1),1,MPI_INTEGER,
     x                0,MPI_COMM_WORLD,ierr)

! Get displacements for array storage
      if (rank.eq.0) then
        displarr(1)=0
        do i=2,nproc
          displarr(i)=displarr(i-1)+ng4locarr(i-1)
        end do
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Form global GM4 on root
      call MPI_GATHERV(GM4(1),ng4loc,MPI_DOUBLE_PRECISION,
     x                 TGM4(1),ng4locarr,displarr,
     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if (rank.eq.0) then
!       write(*,*) "concatenated ng4"
!       do i=1,ng4
!        write(*,9001) TGM4(i)
!       end do
!      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (rank.eq.0) then
       open(unit=20,file="XCHF_GAM4.ufm",form="unformatted")
       write(20) TGM4
       close(20)
       write(*,*) "XCHF_GAM4 written to disk"
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(allocated(TGM4)) deallocate(TGM4)

 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 5-PARTICLE INTEGRALS        |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 5-PARTICLE INTEGRALS: ',1X,I12/
     x        8X,'  NUMBER OF BLOCKS (USER DEFINED): ',1X,I12/
     x        8X,'                          ',/,
     x        8X,'  COMPUTATIONAL RESOURCES:',/,
     x        8X,'  ------------------------',/)
                       
 1500 FORMAT( 8X,'      MPI PROCESSES:',1X,I3/
     x        8X,'        OMP THREADS:',1X,I3/)

 2000 FORMAT(8X,'    TIME FOR PROCESS ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,'    TIME FOR RESIDUAL ',1X,F10.2)

 9001 FORMAT(1X,1(F20.10))

      return
      end
!======================================================================
      subroutine RXCHF_GAM4_thread_MPI(ne,np,ngee,ng2,ng4,ng4_seg,
     x                                 istart,iend,loop_map,arrstart,
     x                                 GAM_ee,GAM_2s,GAM_4)
!
!======================================================================
      implicit none
      include 'omp_lib.h'

      integer istart,iend,imap
      integer ne    ! Number of contracted electronic basis functions
      integer np    ! Number of nuclear basis functions
      integer ngee  ! Number of contracted 2-electron integrals
      integer ng2   ! Number of contracted 3-particle integrals
      integer ng4_seg ! dimension of chunk of contracted 5-particle integrals
      integer ng4   ! Number of integrals calc by MPI process
      integer arrstart ! Index of first integral

      integer loop_map(ng4_seg,10)
      double precision GAM_ee(ngee)   ! Array storage of 2e-integrals
      double precision GAM_2s(ng2)    ! Array storage of 3-part overlaps
      double precision GAM_4(ng4)     ! Array storage of 5-particle ints

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
!---OPENMP-RELATED-VARIABLES-----)

! Local variables
      integer ia
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
!CWS_int_stat
!     integer iii
!     integer jjj
!     integer icount1
!     integer icount2
!     double precision threshold1
!     double precision threshold2
!CWS_int_stat

      double precision ans

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map,arrstart)
!$ompx shared(istart,iend)
!$ompx shared(ne,np)
!$ompx shared(ngee)
!$ompx shared(ng2)
!$ompx shared(ng4_seg)
!$ompx shared(GAM_2s)
!$ompx shared(GAM_ee)
!$ompx shared(GAM_4)
!$ompx shared(IFIL)
!$ompx private(iLp) 
!$ompx private(imap) 
!$ompx private(ia)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(ie3,je3) 
!$ompx private(ie4,je4) 
!$ompx private(ans)
!$ompx private(id)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         je4=loop_map(imap,1) 
         ie4=loop_map(imap,2) 
         je3=loop_map(imap,3) 
         ie3=loop_map(imap,4) 
         je2=loop_map(imap,5) 
         ie2=loop_map(imap,6) 
         je1=loop_map(imap,7) 
         ie1=loop_map(imap,8) 
         jp =loop_map(imap,9) 
         ip =loop_map(imap,10)

         call index_GAM_4PK2(ne,np,
     x                       ip,jp,
     x                       ie1,je1,
     x                       ie2,je2,
     x                       ie3,je3,
     x                       ie4,je4,ia)


         call RXCHFmult_symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

         GAM_4(ia-arrstart+1)=ans

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

      return
      end 
!======================================================================
      subroutine XCHF_GAM4_thread_MPI(ne,np,ngee,ng2,ng4,ng4_seg,
     x                                istart,iend,loop_map,arrstart,
     x                                GAM_ee,GAM_2s,GAM_4)
!
!======================================================================
      implicit none
      include 'omp_lib.h'

      integer istart,iend,imap
      integer ne    ! Number of contracted electronic basis functions
      integer np    ! Number of nuclear basis functions
      integer ngee  ! Number of contracted 2-electron integrals
      integer ng2   ! Number of contracted 3-particle integrals
      integer ng4_seg ! dimension of chunk of contracted 5-particle integrals
      integer ng4   ! Number of integrals calc by MPI process
      integer arrstart ! Index of first integral

      integer loop_map(ng4_seg,10)
      double precision GAM_ee(ngee)   ! Array storage of 2e-integrals
      double precision GAM_2s(ng2)    ! Array storage of 3-part overlaps
      double precision GAM_4(ng4)     ! Array storage of 5-particle ints

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
!---OPENMP-RELATED-VARIABLES-----)

! Local variables
      integer ia
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
!CWS_int_stat
!     integer iii
!     integer jjj
!     integer icount1
!     integer icount2
!     double precision threshold1
!     double precision threshold2
!CWS_int_stat

      double precision ans
      double precision four

      four=4.0d+00

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map,arrstart)
!$ompx shared(istart,iend)
!$ompx shared(ne,np)
!$ompx shared(ngee)
!$ompx shared(ng2)
!$ompx shared(ng4_seg)
!$ompx shared(four)
!$ompx shared(GAM_2s)
!$ompx shared(GAM_ee)
!$ompx shared(GAM_4)
!$ompx shared(IFIL)
!$ompx private(iLp) 
!$ompx private(imap) 
!$ompx private(ia)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(ie3,je3) 
!$ompx private(ie4,je4) 
!$ompx private(ans)
!$ompx private(id)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         je4=loop_map(imap,1) 
         ie4=loop_map(imap,2) 
         je3=loop_map(imap,3) 
         ie3=loop_map(imap,4) 
         je2=loop_map(imap,5) 
         ie2=loop_map(imap,6) 
         je1=loop_map(imap,7) 
         ie1=loop_map(imap,8) 
         jp =loop_map(imap,9) 
         ip =loop_map(imap,10)

         call index_GAM_4PK2(ne,np,
     x                       ip,jp,
     x                       ie1,je1,
     x                       ie2,je2,
     x                       ie3,je3,
     x                       ie4,je4,ia)


         call symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

         ans=ans/(four*four)
         GAM_4(ia-arrstart+1)=ans

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

      return
      end 

