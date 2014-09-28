!======================================================================
      subroutine RXCHF_GAM4_INT_MPI(nproc,rank,
     x                              Nchunks,nebf,nebfBE,npbf,
     x                              ngee,ng2,ng4,
     x                              ng2loc,ng4loc,
     x                              GM2s)

! Calculates INT_GAM4 integrals (split onto MPI procs)
!  - requires each proc to store all XCHF_GAM2s integrals in memory
!======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ngee,nebf,nebfBE,npbf
      integer ng4,ng2       ! Total number of integrals
      integer ng2loc,ng4loc ! Number of integrals on this proc
      double precision GM2s(ng2loc) ! GAM2s integrals on this MPI proc
      
! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,i
      integer ips,jps,ie1s,je1s,ie2s,je2s,ie3s,je3s,ie4s,je4s

      integer(kind=4) packp
      integer(kind=8) packe

      integer(kind=4),allocatable :: prot_map(:)
      integer(kind=8),allocatable :: elec_map(:)

      double precision, allocatable :: GM4_chk(:) ! INT OMG4 integrals (symm)

      double precision GAM_ee(ngee)

      integer nproc,rank
      integer mpistart,mpiend,arrstart

      integer int_count,nzcount

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

      double precision wtime,wtime1,wtime2

! Have each process calculate ng4/nproc integrals according to rank
! Have last process calculate ng4%nproc remaining integrals
      call get_mpi_range(ng4,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng4

      if (rank.eq.0) then
       write(*,1000) ng4,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      if (rank.eq.0) call read_GAM_ee(nebf,ngee,GAM_ee) 
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

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2000) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      wtime = MPI_WTIME()

      int_count=0

      do ichunk=1,Nchunks

        wtime1 = MPI_WTIME()

! Have threads chop calculation of mpiend-mpistart+1=ng4/nproc integrals
        call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
        ng4_seg=1+iend-istart

        if(allocated(prot_map)) deallocate(prot_map)
        allocate(prot_map(ng4_seg))
        if(allocated(elec_map)) deallocate(elec_map)
        allocate(elec_map(ng4_seg))
        if(allocated(GM4_chk)) deallocate(GM4_chk)
        allocate(GM4_chk(ng4_seg))
        GM4_chk=0.0d+00

! Calculate starting indices to not waste loop iterating
        call RXCHF_index_GAM4_unpk(nebf,nebfBE,npbf,istart,ips,jps,
     x                             ie1s,je1s,ie2s,je2s,
     x                             ie3s,je3s,ie4s,je4s)

! Initialize for loop
        ip=ips
        jp=jps
        ie1=ie1s
        je1=je1s
        ie2=ie2s
        je2=je2s
        ie3=ie3s
        je3=je3s
        ie4=ie4s
        je4=je4s
        Loopi=0
        imas=istart

        do while (imas.le.iend)

          if((imas.ge.istart).and.(imas.le.iend)) then

           Loopi=Loopi+1
           call bitpack(ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,
     x                  packp,packe)
           prot_map(Loopi)=packp
           elec_map(Loopi)=packe

          end if

          call RXCHF_update_GAM4_inds(nebf,nebfBE,npbf,
     x                                ip,jp,ie1,je1,ie2,je2,
     x                                ie3,je3,ie4,je4,istat)
          if((istat.ne.0).and.(imas.ne.iend)) then
           write(*,*) "ERROR with updating GAM4 indices"
           return
          end if

          imas=imas+1

        end do

        call RXCHF_GAM4_thread_MPI(nebf,nebfBE,npbf,ngee,ng2,ng4_seg,
     x                             prot_map,elec_map,
     x                             GAM_ee,TGM2s,
     x                             GM4_chk,nzcount)
        int_count=int_count+nzcount

        call RXCHF_GAM4_nz_MPI(nproc,rank,
     x                         ichunk,ng4_seg,nzcount,prot_map,elec_map,
     x                         GM4_chk)

        wtime2 = MPI_WTIME() - wtime1
        write(*,2001) rank,ichunk,wtime2

      end do !end loop over chunks

      if(allocated(GM4_chk)) deallocate(GM4_chk)
      if(allocated(elec_map)) deallocate(elec_map)
      if(allocated(prot_map)) deallocate(prot_map)

C Concatenate all nonzero integrals and associated indices from various chunk files and write out
      call RXCHF_GAM4_concat_MPI(nproc,rank,
     x                           Nchunks,int_count)

      wtime2 = MPI_WTIME() - wtime
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2010)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      write(*,2011) rank,wtime2,ng4loc,
     x              int_count,dble(int_count)/dble(ng4loc)*1.0d+02
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

 2000 FORMAT(/8X,'  INTEGRAL BLOCK CALCULATION TIMINGS:',/,
     x        8X,' -------------------------------------')

 2001 FORMAT( 8X,' PROCESS ',1X,I4,1X,' BLOCK ',1X,I4,1X,F10.2)

 2010 FORMAT(/8X,'   SUMMARY (TIMING, # CALC, # NONZERO, % SURV)',/,
     x        8X,' -----------------------------------------------')

 2011 FORMAT( 5X,' PROCESS ',1X,I4,1X,F10.2,2(1X,I10),1X,F6.2)

 9001 FORMAT(1X,1(F20.10))

      return
      end
!======================================================================
      subroutine RXCHF_GAM4_XCHF_MPI(nproc,rank,
     x                               Nchunks,nebf,nebfBE,npbf,ngee,
     x                               ng2,ng4,
     x                               ng2loc,ng4loc,
     x                               GM2s)

! Calculates XCHF_GAM4 integrals (split onto MPI procs)
!  - requires each proc to store all XCHF_GAM2s integrals in memory
!======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ngee,nebf,nebfBE,npbf
      integer ng4,ng2       ! Total number of integrals
      integer ng2loc,ng4loc ! Number of integrals on this proc
      double precision GM2s(ng2loc) ! GAM2s integrals on this MPI proc
      
! Local Variables
      integer ngeeBE
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,i
      integer ips,jps,ie1s,je1s,ie2s,je2s,ie3s,je3s,ie4s,je4s

      integer(kind=4) packp
      integer(kind=8) packe

      integer(kind=4),allocatable :: prot_map(:)
      integer(kind=8),allocatable :: elec_map(:)

      double precision, allocatable :: GM4_chk(:) ! INT  OMG4 integrals (symm)

      double precision GAM_ee(ngee)
      double precision, allocatable :: GAM_eeBE(:)

      integer nproc,rank
      integer mpistart,mpiend,arrstart

      integer int_count,nzcount

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

      double precision wtime,wtime1,wtime2

      ngeeBE=nebfBE*nebfBE*nebfBE*nebfBE
      if(allocated(GAM_eeBE)) deallocate(GAM_eeBE)
      allocate(GAM_eeBE(ngeeBE))
      GAM_eeBE=zero

! Have each process calculate ng4/nproc integrals according to rank
! Have last process calculate ng4%nproc remaining integrals
      call get_mpi_range(ng4,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng4

      if (rank.eq.0) then
       write(*,1000) ng4,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      if (rank.eq.0) then
!       write(*,*) "nebf,nebfBE,ngee,ngeeBE:",
!     x             nebf,nebfBE,ngee,ngeeBE
!       write(*,*) "reading gam_ee..."
       call read_GAM_ee(nebf,ngee,GAM_ee) 
!       write(*,*) "packing gam_ee..."
       call RXCHF_pack_gamee(nebf,nebfBE,ngee,ngeeBE,GAM_ee,GAM_eeBE)
!       write(*,*) "done"
      end if
      call MPI_BCAST(GAM_eeBE,ngeeBE,MPI_DOUBLE_PRECISION,
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

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2000) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      wtime = MPI_WTIME()

      int_count=0

      do ichunk=1,Nchunks

        wtime1 = MPI_WTIME()

! Have threads chop calculation of mpiend-mpistart+1=ng4/nproc integrals
        call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
        ng4_seg=1+iend-istart

        if(allocated(prot_map)) deallocate(prot_map)
        allocate(prot_map(ng4_seg))
        if(allocated(elec_map)) deallocate(elec_map)
        allocate(elec_map(ng4_seg))
        if(allocated(GM4_chk)) deallocate(GM4_chk)
        allocate(GM4_chk(ng4_seg))
        GM4_chk=0.0d+00

! Calculate starting indices to not waste loop iterating
        call index_GAM4_unpk(nebfBE,npbf,istart,ips,jps,ie1s,je1s,
     x                       ie2s,je2s,ie3s,je3s,ie4s,je4s)

! Initialize for loop
        ip=ips
        jp=jps
        ie1=ie1s
        je1=je1s
        ie2=ie2s
        je2=je2s
        ie3=ie3s
        je3=je3s
        ie4=ie4s
        je4=je4s
        Loopi=0
        imas=istart

        do while (imas.le.iend)

          if((imas.ge.istart).and.(imas.le.iend)) then

           Loopi=Loopi+1
           call bitpack(ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,
     x                  packp,packe)
           prot_map(Loopi)=packp
           elec_map(Loopi)=packe

          end if

          call update_GAM4_inds(nebfBE,npbf,
     x                          ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,
     x                          istat)
          if((istat.ne.0).and.(imas.ne.iend)) then
           write(*,*) "ERROR with updating GAM4 indices"
           return
          end if

          imas=imas+1

        end do

        call XCHF_GAM4_thread_MPI(nebfBE,npbf,ngeeBE,ng2,ng4_seg,
     x                            prot_map,elec_map,
     x                            GAM_eeBE,TGM2s,
     x                            GM4_chk,nzcount)
        int_count=int_count+nzcount

        call RXCHF_GAM4_nz_MPI(nproc,rank,
     x                         ichunk,ng4_seg,nzcount,prot_map,elec_map,
     x                         GM4_chk)

        wtime2 = MPI_WTIME() - wtime1
        write(*,2001) rank,ichunk,wtime2

      end do !end loop over chunks

      if(allocated(GM4_chk)) deallocate(GM4_chk)
      if(allocated(elec_map)) deallocate(elec_map)
      if(allocated(prot_map)) deallocate(prot_map)

C Concatenate all nonzero integrals and associated indices from various chunk files and write out
      call XCHF_GAM4_concat_MPI(nproc,rank,
     x                          Nchunks,int_count)

      wtime2 = MPI_WTIME() - wtime
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2010)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      write(*,2011) rank,wtime2,ng4loc,
     x              int_count,dble(int_count)/dble(ng4loc)*1.0d+02
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng4loc:",mpistart,mpiend,ng4loc
!      do i=1,ng4loc
!       write(*,9001) GM4(i)
!      end do

      if(allocated(TGM2s)) deallocate(TGM2s)
      if(allocated(GAM_eeBE)) deallocate(GAM_eeBE)

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

 2000 FORMAT(/8X,'  INTEGRAL BLOCK CALCULATION TIMINGS:',/,
     x        8X,' -------------------------------------')

 2001 FORMAT( 8X,' PROCESS ',1X,I4,1X,' BLOCK ',1X,I4,1X,F10.2)

 2010 FORMAT(/8X,'   SUMMARY (TIMING, # CALC, # NONZERO, % SURV)',/,
     x        8X,' -----------------------------------------------')

 2011 FORMAT( 5X,' PROCESS ',1X,I4,1X,F10.2,2(1X,I10),1X,F6.2)

 9001 FORMAT(1X,1(F20.10))

      return
      end
!======================================================================
      subroutine RXCHF_GAM4_thread_MPI(ne,neBE,np,ngee,ng2,ng4,
     x                                 prot_map,elec_map,
     x                                 GAM_ee,GAM_2s,
     x                                 GM4_chk,nzcount)
!
!======================================================================
      implicit none
      include 'omp_lib.h'

      integer ne   ! Number of contracted electronic basis functions
      integer neBE ! Number of contracted special electronic basis functions
      integer np   ! Number of nuclear basis functions
      integer ngee ! Number of contracted 2-electron integrals
      integer ng2  ! Number of contracted 3-particle integrals
      integer ng4  ! dimension of chunk of contracted 5-particle integrals

      integer(kind=4)  prot_map(ng4) ! Array of packed proton indices
      integer(kind=8)  elec_map(ng4) ! Array of packed electron indices
      double precision GAM_ee(ngee)  ! Array storage of 2e-integrals
      double precision GAM_2s(ng2)   ! Array storage of 3-part overlaps

      double precision GM4_chk(ng4)  ! GAM4 integrals
      integer          nzcount       ! # of nontrivial INT_GAM4 integrals

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

      nzcount=0

!$omp parallel 
!$ompx shared(prot_map,elec_map)
!$ompx shared(GM4_chk)
!$ompx shared(ne,neBE,np)
!$ompx shared(ngee)
!$ompx shared(ng2)
!$ompx shared(ng4)
!$ompx shared(GAM_2s)
!$ompx shared(GAM_ee)
!$ompx shared(tol)
!$ompx private(i) 
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(ie3,je3) 
!$ompx private(ie4,je4) 
!$ompx private(ans)
!$ompx reduction(+:nzcount)

!$omp do SCHEDULE(RUNTIME)
      do i=1,ng4

         call bitunpack(prot_map(i),elec_map(i),ip,jp,
     x                  ie1,je1,ie2,je2,ie3,je3,ie4,je4)

         call RXCHFmult_symm_gam4(ne,neBE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

         if(abs(ans).gt.tol) then
          nzcount=nzcount+1
          GM4_chk(i)=ans
         end if

      end do
!$omp end do
!$omp end parallel      

      return
      end 
!======================================================================
      subroutine XCHF_GAM4_thread_MPI(ne,np,ngee,ng2,ng4,
     x                                prot_map,elec_map,
     x                                GAM_ee,GAM_2s,
     x                                GM4_chk,nzcount)
!
!======================================================================
      implicit none
      include 'omp_lib.h'

      integer ne    ! Number of contracted electronic basis functions
      integer np    ! Number of nuclear basis functions
      integer ngee  ! Number of contracted 2-electron integrals
      integer ng2   ! Number of contracted 3-particle integrals
      integer ng4  ! dimension of chunk of contracted 5-particle integrals

      integer(kind=4)  prot_map(ng4)  ! Array of packed proton indices
      integer(kind=8)  elec_map(ng4)  ! Array of packed electron indices
      double precision GAM_ee(ngee)   ! Array storage of 2e-integrals
      double precision GAM_2s(ng2)    ! Array storage of 3-part overlaps

      double precision GM4_chk(ng4)   ! XCHF_GAM4 integrals
      integer          nzcount        ! # of nontrivial XCHF_GAM4 integrals

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

      tol=1.0d-15
      four=4.0d+00

      nzcount=0

!$omp parallel 
!$ompx shared(prot_map,elec_map)
!$ompx shared(GM4_chk)
!$ompx shared(ne,np)
!$ompx shared(ngee)
!$ompx shared(ng2)
!$ompx shared(ng4)
!$ompx shared(four)
!$ompx shared(GAM_2s)
!$ompx shared(GAM_ee)
!$ompx shared(tol)
!$ompx private(i) 
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(ie3,je3) 
!$ompx private(ie4,je4) 
!$ompx private(ans)
!$ompx reduction(+:nzcount)

!$omp do SCHEDULE(RUNTIME)
      do i=1,ng4

         call bitunpack(prot_map(i),elec_map(i),ip,jp,
     x                  ie1,je1,ie2,je2,ie3,je3,ie4,je4)

         call symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

         ans=ans/(four*four)
         if(abs(ans).gt.tol) then
          nzcount=nzcount+1
          GM4_chk(i)=ans
         end if

      end do
!$omp end do
!$omp end parallel      

      return
      end 
!======================================================================
      subroutine RXCHF_GAM4_nz_MPI(nproc,rank,
     x                             ichunk,ng4,ng4_nz,protinds,elecinds,
     x                             ints_all)
!
! Packs all nonzero integrals from ints_all into local ints_nz eliminating
! the zero entries and determines prot/elec indices in pinds,einds
! The three local arrays are written to disk
!
! Input:
!   ints_all(ng4)   : Array of whole chunk size of integrals
!
! Output:
!   GAM4_[rank]_[chunk].ufm : unformatted file with writes
!                             1) # nonzero integrals for this chunk
!                             2) pinds array as one write
!                             3) einds array as one write
!                             4) ints_nz array as one write
!======================================================================
      implicit none

! Input variables
      integer nproc,rank
      integer ichunk                  ! Index of chunk being processed
      integer ng4                     ! Dimension of chunk of contracted 5-particle integrals
      integer ng4_nz                  ! Number of nonzero integrals calc by MPI process
      integer(kind=4)  protinds(ng4)  ! All proton indices (packed)
      integer(kind=8)  elecinds(ng4)  ! All electron indices (packed)
      double precision ints_all(ng4)  ! All integrals

! Local variables
      integer i
      integer currind

      integer(kind=4)  pinds(ng4_nz)   ! Nonzero proton indices (packed)
      integer(kind=8)  einds(ng4_nz)   ! Nonzero electron indices (packed)
      double precision ints_nz(ng4_nz) ! Nonzero integrals

      character*4 istring,jstring      ! File I/O variables

      double precision tol

      tol=1.0d-15

C Variables for file I/O
      write(istring,'(I4.4)') rank
      write(jstring,'(I4.4)') ichunk

C Pack nonzero integrals and corresponding indices
      currind=0
      do i=1,ng4
        if(abs(ints_all(i)).gt.tol) then
         currind=currind+1
         pinds(currind)=protinds(i)
         einds(currind)=elecinds(i)
         ints_nz(currind)=ints_all(i)
        end if
      end do

      if(currind.ne.ng4_nz) then
       write(*,*) "ERROR in RXCHF_GAM4_nz_MPI:",rank
       write(*,*) "Found diff # of nz integrals on second pass:",
     x            currind,ng4_nz
      end if

C Write to file
      call writebitint(nproc,rank,ng4_nz,18,
     x                 "GAM4_"//istring//"_"//jstring//".ufm",
     x                 pinds,einds,ints_nz)

      return
      end 
!======================================================================
      subroutine RXCHF_GAM4_concat_MPI(nproc,rank,
     x                                 nchunk,nz_tot)
!
!======================================================================
      implicit none

! Input variables
      integer nproc,rank
      integer nchunk     ! Number of chunks written
      integer nz_tot     ! Number of total nonzero integrals

! Local variables
      integer i
      integer currind
      integer nz_chk     ! Number of nonzero integrals for given chunk

      integer(kind=4)  pinds(nz_tot)   ! All nonzero proton indices (packed)
      integer(kind=8)  einds(nz_tot)   ! All nonzero electron indices (packed)
      double precision ints_nz(nz_tot) ! All nonzero integrals

      integer     unitno               !
      character*4 istring,jstring      ! File I/O variables

C Variables for file I/O
      write(istring,'(I4.4)') rank
      unitno=200+rank

      currind=1

      do i=1,Nchunk

        write(jstring,'(I4.4)') i

C Read integrals from this chunk
        open(unit=unitno,file="GAM4_"//istring//"_"//jstring//".ufm",
     x       form="unformatted")
        read(unitno) nz_chk
        call RXCHF_readbit_mpi(nproc,rank,nz_chk,unitno,
     x                         pinds(currind),einds(currind),
     x                         ints_nz(currind))
        close(unitno,status="delete")

C Update index for next read
        currind=currind+nz_chk

      end do

C Write to file
      call writebitint(nproc,rank,nz_tot,17,
     x                 "INT_GAM4-"//istring//".ufm",
     x                 pinds,einds,ints_nz)

      return
      end 
!======================================================================
      subroutine XCHF_GAM4_concat_MPI(nproc,rank,
     x                                nchunk,nz_tot)
!
!======================================================================
      implicit none

! Input variables
      integer nproc,rank
      integer nchunk     ! Number of chunks written
      integer nz_tot     ! Number of total nonzero integrals

! Local variables
      integer i
      integer currind
      integer nz_chk     ! Number of nonzero integrals for given chunk

      integer(kind=4)  pinds(nz_tot)   ! All nonzero proton indices (packed)
      integer(kind=8)  einds(nz_tot)   ! All nonzero electron indices (packed)
      double precision ints_nz(nz_tot) ! All nonzero integrals

      integer     unitno               !
      character*4 istring,jstring      ! File I/O variables

C Variables for file I/O
      write(istring,'(I4.4)') rank
      unitno=200+rank

      currind=1

      do i=1,Nchunk

        write(jstring,'(I4.4)') i

C Read integrals from this chunk
        open(unit=unitno,file="GAM4_"//istring//"_"//jstring//".ufm",
     x       form="unformatted")
        read(unitno) nz_chk
        call RXCHF_readbit_mpi(nproc,rank,nz_chk,unitno,
     x                         pinds(currind),einds(currind),
     x                         ints_nz(currind))
        close(unitno,status="delete")

C Update index for next read
        currind=currind+nz_chk

      end do

C Write to file
      call writebitint(nproc,rank,nz_tot,18,
     x                 "XCHF_GAM4-"//istring//".ufm",
     x                 pinds,einds,ints_nz)

      return
      end 
!======================================================================
      subroutine RXCHF_index_GAM4_unpk(ne,neBE,np,ipk,ip,jp,ie1,je1,
     x                                 ie2,je2,ie3,je3,ie4,je4)
!
! Calculates all proton and electron indices corresponding to some
! flattened index [ipk] (essentially inverse of RXCHF_index_GAM_4PK)
!======================================================================
      implicit none

! Input variables
      integer ne        ! Number of total electron basis functions
      integer neBE      ! Number of special electron basis functions
      integer np        ! Number of nuclear  basis functions
      integer ipk       ! Starting flattened index

! Output variables
      integer ip,jp     ! Starting proton indices
      integer ie1,je1   ! Starting electron indices
      integer ie2,je2   !
      integer ie3,je3   !
      integer ie4,je4   !

! Local variables
      integer i

      i=ipk-1

      je4=1+mod(i,neBE)
      i=i/neBE

      ie4=1+mod(i,neBE)
      i=i/neBE

      je3=1+mod(i,neBE)
      i=i/neBE

      ie3=1+mod(i,neBE)
      i=i/neBE

      je2=1+mod(i,neBE)
      i=i/neBE

      ie2=1+mod(i,neBE)
      i=i/neBE

      je1=1+mod(i,ne)
      i=i/ne

      ie1=1+mod(i,ne)
      i=i/ne

      jp=1+mod(i,np)
      i=i/np

      ip=1+i

      return
      end
!======================================================================
      subroutine index_GAM4_unpk(ne,np,ipk,ip,jp,ie1,je1,
     x                           ie2,je2,ie3,je3,ie4,je4)
!
! Calculates all proton and electron indices corresponding to some
! flattened index [ipk] (essentially inverse of index_GAM_4PK)
!======================================================================
      implicit none

! Input variables
      integer ne        ! Number of electron basis functions
      integer np        ! Number of nuclear  basis functions
      integer ipk       ! Starting flattened index

! Output variables
      integer ip,jp     ! Starting proton indices
      integer ie1,je1   ! Starting electron indices
      integer ie2,je2   !
      integer ie3,je3   !
      integer ie4,je4   !

! Local variables
      integer i

      i=ipk-1

      je4=1+mod(i,ne)
      i=i/ne

      ie4=1+mod(i,ne)
      i=i/ne

      je3=1+mod(i,ne)
      i=i/ne

      ie3=1+mod(i,ne)
      i=i/ne

      je2=1+mod(i,ne)
      i=i/ne

      ie2=1+mod(i,ne)
      i=i/ne

      je1=1+mod(i,ne)
      i=i/ne

      ie1=1+mod(i,ne)
      i=i/ne

      jp=1+mod(i,np)
      i=i/np

      ip=1+i

      return
      end
!======================================================================
      subroutine RXCHF_update_GAM4_inds(ne,neBE,np,
     x                                  ip,jp,ie1,je1,ie2,je2,
     x                                  ie3,je3,ie4,je4,ierr)
!
! Increments GAM4 indices for restricted basis set packing
!======================================================================
      implicit none

! Input variables
      integer ne        ! Number of electron basis functions
      integer neBE      ! Number of special electron basis functions
      integer np        ! Number of nuclear  basis functions

! Input/Output variables
      integer ip,jp     ! Starting proton indices
      integer ie1,je1   ! Starting electron indices
      integer ie2,je2   !
      integer ie3,je3   !
      integer ie4,je4   !

! Output variables
      integer ierr      ! =0 if clean exit; =1 if overflow

      ierr=0

      if (je4.lt.neBE) then
       je4=je4+1
      else
       je4=1

       if (ie4.lt.neBE) then
        ie4=ie4+1
       else
        ie4=1

        if (je3.lt.neBE) then
         je3=je3+1
        else
         je3=1

         if (ie3.lt.neBE) then
          ie3=ie3+1
         else
          ie3=1

          if (je2.lt.neBE) then
           je2=je2+1
          else
           je2=1
          
           if (ie2.lt.neBE) then
            ie2=ie2+1
           else
            ie2=1

            if (je1.lt.ne) then
             je1=je1+1
            else
             je1=1
           
             if (ie1.lt.ne) then
              ie1=ie1+1
             else
              ie1=1

              if (jp.lt.np) then
               jp=jp+1
              else
               jp=1
           
               if (ip.lt.np) then
                ip=ip+1
               else
                ip=1
                ierr=1
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

      return
      end
!======================================================================
      subroutine update_GAM4_inds(ne,np,
     x                            ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,
     x                            ierr)
!
! Increments GAM4 indices
!======================================================================
      implicit none

! Input variables
      integer ne        ! Number of electron basis functions
      integer np        ! Number of nuclear  basis functions

! Input/Output variables
      integer ip,jp     ! Starting proton indices
      integer ie1,je1   ! Starting electron indices
      integer ie2,je2   !
      integer ie3,je3   !
      integer ie4,je4   !

! Output variables
      integer ierr      ! =0 if clean exit; =1 if overflow

      ierr=0

      if (je4.lt.ne) then
       je4=je4+1
      else
       je4=1

       if (ie4.lt.ne) then
        ie4=ie4+1
       else
        ie4=1

        if (je3.lt.ne) then
         je3=je3+1
        else
         je3=1

         if (ie3.lt.ne) then
          ie3=ie3+1
         else
          ie3=1

          if (je2.lt.ne) then
           je2=je2+1
          else
           je2=1
          
           if (ie2.lt.ne) then
            ie2=ie2+1
           else
            ie2=1

            if (je1.lt.ne) then
             je1=je1+1
            else
             je1=1
           
             if (ie1.lt.ne) then
              ie1=ie1+1
             else
              ie1=1

              if (jp.lt.np) then
               jp=jp+1
              else
               jp=1
           
               if (ip.lt.np) then
                ip=ip+1
               else
                ip=1
                ierr=1
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

      return
      end
!======================================================================
      subroutine RXCHF_pack_gamee(nebf,nebfBE,ngee,ngeeBE,GAM,GAMBE)
!
! Packs two electron integrals (nebf^4) into those for special electrons
! (nebfBE^4) from GAM -> GAMBE
!======================================================================
      implicit none

! Input variables
      integer nebf                   ! Number of all electron basis functions
      integer nebfBE                 ! Number of restricted electron basis functions
      integer ngee                   ! Number of 2-e integrals for all electron basis
      integer ngeeBE                 ! Number of 2-e integrals for restricted electron basis
      double precision GAM(ngee)     ! 2-e integrals for all electron basis

! Output variables
      double precision GAMBE(ngeeBE) ! 2-e integrals for restricted electron basis

! Local variables
      integer ie1,je1
      integer ie2,je2
      integer iold,inew

      do ie1=1,nebfBE
      do je1=1,nebfBE
        do ie2=1,nebfBE
        do je2=1,nebfBE
          call pack_4D(nebf,nebf,nebf,je2,ie2,je1,ie1,iold)
          call pack_4D(nebfBE,nebfBE,nebfBE,je2,ie2,je1,ie1,inew)
          GAMBE(inew)=GAM(iold)
        end do
        end do
      end do
      end do

      return
      end

