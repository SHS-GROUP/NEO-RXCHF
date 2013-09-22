!======================================================================
      subroutine RXCHF_GAM4_MPI(nproc,rank,
     x                          Nchunks,ne,np,ngee,ng2,ng4,
     x                          GAM_2s,GAM_4)

!======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ng4,ng2,ngee,ne,np
      double precision GAM_4(ng4)
      double precision GAM_2s(ng2)
      
! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4
      integer,allocatable :: loop_map(:,:)
      double precision,allocatable :: GAM_ee(:)
      double precision wtime,wtime2

      double precision GM4_MPI(ng4/nproc)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

! Have each process calculate ng4/nproc integrals according to rank then allgather
! Have each process calculate ng4%nproc remaining integrals
!  - saves waiting for master and eliminates one broadcast

      call get_mpi_range(ng4,nproc,rank,mpistart,mpiend)

      if (rank.eq.0) then
       write(*,1000) ng4,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      wtime = omp_get_wtime()

!----READ-GAM_ee--INTO-MEMORY-------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      allocate( GAM_ee(ngee),stat=istat )
      if (rank.eq.0) call read_GAM_ee(ne,ngee,GAM_ee) 
      call MPI_BCAST(GAM_ee,ngee,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
!----READ-GAM_ee-AND-GAM_2s-INTO-MEMORY-------------------------------)

      GM4_MPI=0.0d+00

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(
      wtime = MPI_WTIME()

      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng3/nproc integrals
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

         call RXCHFmult_thread_gam4_IC(ne,np,ngee,ng2,ng4,ng4_seg,
     x                       istart,iend,loop_map,GAM_ee,GAM_2s,GAM_4)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)
!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

      wtime2 = MPI_WTIME() - wtime
      write(*,2000)rank,wtime2

! Here, a contiguous block of ng4/nproc integrals are stored in GAM4 arr
! Pass this to a ng4/nproc-dimensional array to prepare for allgather
      call copy_arr(ng4/nproc,GAM_4(arrstart),GM4_MPI)

! Pass and broadcast these major chunks of arrays to all processes
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHER(GM4_MPI(1),ng4/nproc,MPI_DOUBLE_PRECISION,
     x                   GAM_4(1),ng4/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM4 allgather"

! At this stage, ng4%nproc last elements are missing from arrays
      if (mod(ng4,nproc).ne.0) then

         wtime = MPI_WTIME()

! Have threads chop calculation of ng4%nproc integrals
         call loop_size(nproc*(ng4/nproc)+1,ng4,1,0,
     x                  istart,iend)
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
                            if (Loopi.eq.1) then
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

         call RXCHFmult_thread_gam4_IC(ne,np,ngee,ng2,ng4,ng4_seg,
     x                       istart,iend,loop_map,GAM_ee,GAM_2s,GAM_4)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
         if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

         wtime2 = MPI_WTIME() - wtime
         if (rank.eq.0) write(*,3000) wtime2

!         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!         write(*,*) "ng4+resid:",ng4
!         do istart=1,ng4
!         write(*,9001) GAM_4(istart)
!         end do

      end if ! resid exists

      if(allocated(GAM_ee)) deallocate(GAM_ee)

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
      subroutine XCHF_GAM4_MPI(nproc,rank,
     s                         Nchunks,ne,np,ngee,ng2,ng4,
     x                         GAM_2s,GAM_4)

!======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ng4,ng2,ngee,ne,np
      double precision GAM_4(ng4)
      double precision GAM_2s(ng2)
      
! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4
      integer,allocatable :: loop_map(:,:)
      double precision,allocatable :: GAM_ee(:)
      double precision wtime,wtime2

      double precision GM4_MPI(ng4/nproc)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

! Have each process calculate ng4/nproc integrals according to rank then allgather
! Have each process calculate ng4%nproc remaining integrals
!  - saves waiting for master and eliminates one broadcast

      call get_mpi_range(ng4,nproc,rank,mpistart,mpiend)

      if (rank.eq.0) then
       write(*,1000) ng4,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      wtime = omp_get_wtime()

!----READ-GAM_ee--INTO-MEMORY-------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      allocate( GAM_ee(ngee),stat=istat )
      if (rank.eq.0) call read_GAM_ee(ne,ngee,GAM_ee) 
      call MPI_BCAST(GAM_ee,ngee,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
!----READ-GAM_ee-AND-GAM_2s-INTO-MEMORY-------------------------------)

      GM4_MPI=0.0d+00

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(
      wtime = MPI_WTIME()

      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng3/nproc integrals
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

         call thread_gam4_IC(ne,np,ngee,ng2,ng4,ng4_seg,istart,iend,
     x                       loop_map,GAM_ee,GAM_2s,GAM_4)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)
!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

      wtime2 = MPI_WTIME() - wtime
      write(*,2000)rank,wtime2

! Here, a contiguous block of ng4/nproc integrals are stored in GAM4 arr
! Pass this to a ng4/nproc-dimensional array to prepare for allgather
      call copy_arr(ng4/nproc,GAM_4(arrstart),GM4_MPI)

! Pass and broadcast these major chunks of arrays to all processes
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHER(GM4_MPI(1),ng4/nproc,MPI_DOUBLE_PRECISION,
     x                   GAM_4(1),ng4/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM4 allgather"

! At this stage, ng4%nproc last elements are missing from arrays
      if (mod(ng4,nproc).ne.0) then

         wtime = MPI_WTIME()

! Have threads chop calculation of ng4%nproc integrals
         call loop_size(nproc*(ng4/nproc)+1,ng4,1,0,
     x                  istart,iend)
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
                            if (Loopi.eq.1) then
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

         call thread_gam4_IC(ne,np,ngee,ng2,ng4,ng4_seg,istart,iend,
     x                       loop_map,GAM_ee,GAM_2s,GAM_4)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
         if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

         wtime2 = MPI_WTIME() - wtime
         if (rank.eq.0) write(*,3000) wtime2

!         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!         write(*,*) "ng4+resid:",ng4
!         do istart=1,ng4
!         write(*,9001) GAM_4(istart)
!         end do

      end if ! resid exists

      if(allocated(GAM_ee)) deallocate(GAM_ee)


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

