!=======================================================================
      subroutine RXCHF_GAM3_MPI(nproc,rank,
     x                          Nchunks,nebf,npebf,npbf,
     x                          ng3,ng3loc,ng3prm,nat,ngtg1,
     x                          pmass,cat,zan,bcoef1,gamma1,
     x                          KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                          ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                          GM3_1,GM3_2)

!=======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ng3             ! Total number of integrals
      integer ng3loc          ! Number of integrals for MPI proc to calc
      integer nebf,npebf,npbf,ng3prm
      integer nat,ngtg1
!-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
      integer KPESTR(nebf)  ! Map contracted index to primitive start
      integer KPEEND(nebf)  ! Map contracted index to primitive end
!-------Basis Set Info-------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)

! Variables Returned
      double precision GM3_1(ng3loc)  ! XCHF OMG3  integrals (full symm)
      double precision GM3_2(ng3loc)  ! INT  OMG3  integrals (partial symm)

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,i
      integer,allocatable :: loop_map(:,:)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

      integer*4 tag_2
      integer*4 sendrank,recvrank
      integer*4, allocatable :: reqs(:),statarr(:,:)
      integer ia_first,ia_last,nmsgs,reqscount
      logical locblock

      double precision, allocatable :: XGM3_2(:)

      integer ia
      integer ia_123
      integer ia_132
      integer ia_213
      integer ia_231
      integer ia_312
      integer ia_321

      double precision x123
      double precision x132
      double precision x213
      double precision x321
      double precision xxxx

      double precision zero,half,six
      parameter(zero=0.0d+00,half=0.5d+00,six=6.0d+00)

      double precision wtime
      double precision wtime2

! Testing Variables
!      double precision, allocatable :: TGM3_1(:) ! Testing arrays
!      double precision, allocatable :: TGM3_2(:)
!
!      integer*4 ng3loc4
!      integer*4 ng3locarr(nproc),displarr(nproc)


! Have each process calculate ng3/nproc integrals according to rank
! Have last process calculate ng3%nproc remaining integrals
      call get_mpi_range(ng3,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng3

      if (rank.eq.0) then
       write(*,1000) ng3,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      GM3_1=0.0d+00
      GM3_2=0.0d+00

!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      wtime = MPI_WTIME()

      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng3/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng3_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )

! Nested loop compression for this chunk:
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

! Save coordinates of first value for array transfer later
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

         call RXCHF_GAM3_thread_MPI(istart,iend,ng3_seg,ng3loc,
     x                         nebf,npebf,npbf,nat,ngtg1,
     x                         pmass,cat,zan,bcoef1,gamma1,
     x                         loop_map,arrstart,GM3_1,GM3_2,
     x                         KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                         ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------)
!-----CLEAN-UP-MEMORY-------------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

      wtime2 = MPI_WTIME() - wtime

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(rank.eq.0) write(*,2000) 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      write(*,2001) rank,wtime2

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng3loc:",mpistart,mpiend,ng3loc
!      do i=1,ng3loc
!       write(*,9001) GM3_1(i),GM3_2(i)
!      end do

!--------------------SYMMETRIZE----------------------------------------(
C Symmetrized integrals in GM3_1ICR (XCHF integrals)
C  - completed using separate routine to reduce memory requirements
C Symmetrized integrals in GM3_2ICR (interaction integrals)
C  - completed simultaneously

      wtime = MPI_WTIME() 

! Allocate storage for temporary arrays to store ia_ji integrals
      if(allocated(XGM3_2)) deallocate(XGM3_2)
      allocate(XGM3_2(ng3loc)) ! stores ia_132 integrals
      XGM3_2=0.0d+00

! Get appropriate ia_ji values for XGM3_2
      if (rank.eq.0) write(*,*) "Symmetrizing GM3_2"

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf

! Get indices of first/last integrals to be passed for fixed ip,jp,iec1,jec1
            call index_GAM_3PK(nebf,npbf,
     x                         ip,jp,iec1,jec1,
     x                         1,1,1,1,ia_first)
            call index_GAM_3PK(nebf,npbf,
     x                         ip,jp,iec1,jec1,
     x                         nebf,nebf,nebf,nebf,ia_last)
! Determine how many messages (send+recv) will be passed by this MPI proc
! Only pass through MPI in case of block being split over >1 MPI proc
            if((mpistart.le.ia_first).and.(mpiend.ge.ia_last)) then
             nmsgs=2*(ia_last-ia_first+1)
             locblock=.true.
            else
             locblock=.false.
             if((mpistart.ge.ia_first).and.(mpistart.le.ia_last)) then
              if (mpiend.le.ia_last) then
               nmsgs=2*(mpiend-mpistart+1)
              else
               nmsgs=2*(ia_last-mpistart+1)
              end if
             else if((mpiend.ge.ia_first).and.(mpiend.le.ia_last)) then
              if (mpistart.le.ia_first) then
               nmsgs=2*(mpiend-ia_first+1)
              else
               nmsgs=2*(mpiend-mpistart+1)
              end if
             else
              nmsgs=1
              locblock=.true.
             end if
            end if
! Allocate and initialize arrays for MPI requests
            if(allocated(reqs)) deallocate(reqs)
            allocate(reqs(nmsgs))
            reqs=MPI_REQUEST_NULL
            reqscount=0

C            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C            if(rank.eq.0) write(*,1001) "index:",rank,ip,jp,iec1,jec1
C            write(*,*) "nmsgs:",rank,nmsgs
C            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec3,jec3,iec2,jec2,ia)
                             ia_132=ia

! ia_123 index is on this MPI process
       if ((ia_123.ge.mpistart).and.(ia_123.le.mpiend)) then

        if (locblock) then
         XGM3_2(ia_123-arrstart+1)=GM3_2(ia_132-arrstart+1)
        else
! MPI receive ia_ji integrals and store in X* arrs at ia_ji index
         call pack_4D(nebf,nebf,nebf,
     x                jec3,iec3,jec2,iec2,ia)
         tag_2=int(ia,kind=4)
         reqscount=reqscount+1

! Get ia_132 for GM3_2
         call get_mpi_proc(ng3,nproc,ia_132,recvrank)
         call MPI_IRECV(XGM3_2(ia_123-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,recvrank,tag_2,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

! ia_132 index is on this MPI process
       if ((ia_132.ge.mpistart).and.(ia_132.le.mpiend)) then

        if (.not.(locblock)) then
! Get rank of MPI process with ia_123 integrals
         call get_mpi_proc(ng3,nproc,ia_123,sendrank)

! MPI send ia_132 integrals for GM3_2
         call pack_4D(nebf,nebf,nebf,
     x                jec3,iec3,jec2,iec2,ia)
         tag_2=int(ia,kind=4)
         reqscount=reqscount+1
         call MPI_ISEND(GM3_2(ia_132-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,sendrank,tag_2,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)

        end if

       end if

               end do
               end do
            end do
            end do

            if(allocated(statarr)) deallocate(statarr)
            allocate(statarr(MPI_STATUS_SIZE,nmsgs))
            call MPI_WAITALL(nmsgs,reqs,statarr,ierr)
            if (ierr.ne.0) write(*,*) "Trouble with reqs waitall"

         end do
         end do
      end do
      end do

      if(allocated(statarr)) deallocate(reqs)
      if(allocated(reqs)) deallocate(reqs)

! Symmetrize integrals locally
      do i=1,ng3loc

        x123=GM3_2(i)
        x132=XGM3_2(i)
        GM3_2(i)=(x123+x132)*half

      end do

! Done with GM3_2 integrals
      if(allocated(XGM3_2)) deallocate(XGM3_2)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Symmetrize GM3_1 integrals
!      if (rank.eq.0) write(*,*) "Symmetrizing GM3_1"

      call RXCHF_GAM3ex_MPI_symm(nproc,rank,
     x                           ng3,ng3loc,
     x                           mpistart,mpiend,arrstart,
     x                           nebf,npbf,
     x                           GM3_1)

!      write(*,*) "start,end,ng3loc:",mpistart,mpiend,ng3loc
!      do i=1,ng3loc
!       write(*,9001) GM3_1(i),GM3_2(i)
!      end do

      wtime2 = MPI_WTIME() - wtime

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(rank.eq.0) write(*,3000) 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      write(*,3001) rank,wtime2
!--------------------SYMMETRIZE----------------------------------------)

! Construct global arrays on master process for testing
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if(allocated(TGM3_1)) deallocate(TGM3_1)
!      if(allocated(TGM3_2)) deallocate(TGM3_2)
!      if (rank.eq.0) then
!       allocate(TGM3_1(ng3))
!       allocate(TGM3_2(ng3))
!      else
!       allocate(TGM3_1(1))
!       allocate(TGM3_2(1))
!      end if
!      TGM3_1=zero
!      TGM3_2=zero
!
!      ng3loc4=int(ng3loc,kind=4)
!
!! Get number of elements calculated by each proc
!      call MPI_GATHER(ng3loc4,1,MPI_INTEGER,
!     x                ng3locarr(1),1,MPI_INTEGER,
!     x                0,MPI_COMM_WORLD,ierr)
!
!! Get displacements for array storage
!      if (rank.eq.0) then
!        displarr(1)=0
!        do i=2,nproc
!          displarr(i)=displarr(i-1)+ng3locarr(i-1)
!        end do
!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM3_1 on root
!      call MPI_GATHERV(GM3_1(1),ng3loc,MPI_DOUBLE_PRECISION,
!     x                 TGM3_1(1),ng3locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM3_2 on root
!      call MPI_GATHERV(GM3_2(1),ng3loc,MPI_DOUBLE_PRECISION,
!     x                 TGM3_2(1),ng3locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!!      if (rank.eq.0) then
!!       write(*,*) "concatenated ng3"
!!       do i=1,ng3
!!        write(*,9001) TGM3_1(i),TGM3_2(i)
!!       end do
!!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if (rank.eq.0) then
!       open(unit=20,file="XCHF_GAM3.ufm",form="unformatted")
!       write(20) TGM3_1
!       close(20)
!       write(*,*) "XCHF_GAM3 written to disk"
!       open(unit=21,file="INT_GAM3.ufm",form="unformatted")
!       write(21) TGM3_2
!       close(21)
!       write(*,*) "INT_GAM3 written to disk"
!      end if
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!      if(allocated(TGM3_2)) deallocate(TGM3_2)
!      if(allocated(TGM3_1)) deallocate(TGM3_1)

 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 4-PARTICLE INTEGRALS        |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 4-PARTICLE INTEGRALS: ',1X,I12/
     x        8X,'  NUMBER OF BLOCKS (USER DEFINED): ',1X,I12/
     x        8X,'                          ',/,
     x        8X,'  COMPUTATIONAL RESOURCES:',/,
     x        8X,'  ------------------------',/)

 1500 FORMAT( 8X,'      MPI PROCESSES:',1X,I3/
     x        8X,'        OMP THREADS:',1X,I3/)

 2000 FORMAT(/8X,'  INTEGRAL CALCULATION TIMINGS:',/,
     x        8X,'  -----------------------------')

 2001 FORMAT( 8X,'    PROCESS ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,' INTEGRAL SYMMETRIZATION TIMINGS:',/,
     x        8X,' --------------------------------')

 3001 FORMAT( 8X,'    PROCESS ',1X,I4,1X,F10.2)

 9001 FORMAT(1X,2(F20.10))

      return
      end
!=======================================================================
      subroutine RXCHF_GAM3ex_MPI(nproc,rank,
     x                            Nchunks,nebf,npebf,npbf,
     x                            ng3,ng3loc,ng3prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            GM3_1,GM3_2,
     x                            GM3_3,GM3_4)

!=======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ng3             ! Total number of integrals
      integer ng3loc          ! Number of integrals for MPI proc to calc
      integer nebf,npebf,npbf,ng3prm
      integer nat,ngtg1
!-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
      integer KPESTR(nebf)  ! Map contracted index to primitive start
      integer KPEEND(nebf)  ! Map contracted index to primitive end
!-------Basis Set Info-------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)

! Variables Returned
      double precision GM3_1(ng3loc)  ! XCHF OMG3    integrals (full symm)
      double precision GM3_2(ng3loc)  ! INT  OMG3    integrals (partial symm)
      double precision GM3_3(ng3loc)  ! INT  OMG3ex1 integrals (partial symm)
      double precision GM3_4(ng3loc)  ! INT  OMG3ex2 integrals (partial symm)

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,i
      integer,allocatable :: loop_map(:,:)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

      integer tag_2,tag_3,tag_4
      integer sendrank,recvrank
      integer*4, allocatable :: reqs(:)

      double precision, allocatable :: XGM3_2(:)
      double precision, allocatable :: XGM3_3(:)
      double precision, allocatable :: XGM3_4(:)

      integer ia
      integer ia_123
      integer ia_132
      integer ia_213
      integer ia_231
      integer ia_312
      integer ia_321

      double precision x123
      double precision x132
      double precision x213
      double precision x321
      double precision xall

      double precision zero,half,six
      parameter(zero=0.0d+00,half=0.5d+00,six=6.0d+00)

      double precision wtime
      double precision wtime2

! Testing Variables
!      double precision, allocatable :: TGM3_1(:) ! Testing arrays
!      double precision, allocatable :: TGM3_2(:)
!      double precision, allocatable :: TGM3_3(:)
!      double precision, allocatable :: TGM3_4(:)
!
!      integer*4 ng3loc4
!      integer*4 ng3locarr(nproc),displarr(nproc)


! Have each process calculate ng3/nproc integrals according to rank
! Have last process calculate ng3%nproc remaining integrals
      call get_mpi_range(ng3,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng3

      if (rank.eq.0) then
       write(*,1000) ng3,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      GM3_1=0.0d+00
      GM3_2=0.0d+00
      GM3_3=0.0d+00
      GM3_4=0.0d+00

!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      wtime = MPI_WTIME()

      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng3/nproc integrals
         call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
         ng3_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )

! Nested loop compression for this chunk:
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

         call RXCHF_GAM3ex_thread_MPI(istart,iend,ng3_seg,ng3loc,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           loop_map,arrstart,
     x                           GM3_1,GM3_2,GM3_3,GM3_4,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------)
!-----CLEAN-UP-MEMORY-------------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

      wtime2 = MPI_WTIME() - wtime

      if(rank.eq.0) write(*,2000) 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      write(*,2001) rank,wtime2

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng3loc:",mpistart,mpiend,ng3loc
!      do i=1,ng3loc
!       write(*,9001) GM3_1(i),GM3_2(i),
!     x               GM3_3(i),GM3_4(i)
!      end do

!--------------------SYMMETRIZE----------------------------------------(
C Symmetrized integrals in GM3_1ICR (XCHF integrals)
C  - completed using separate routine to reduce memory requirements
C Symmetrized integrals in GM3_2ICR (interaction integrals)
C Symmetrized integrals in GM3_3ICR (exchange integrals)
C Symmetrized integrals in GM3_4ICR (exchange integrals)
C  - completed sequentially in this routine

      wtime = MPI_WTIME() 

! Get appropriate ia_ji values for XGM3_2
      if (rank.eq.0) write(*,*) "Symmetrizing GM3_2"

! Allocate storage for temporary arrays to store ia_ji integrals
      if(allocated(XGM3_2)) deallocate(XGM3_2)
      allocate(XGM3_2(ng3loc)) ! stores ia_132 integrals
      XGM3_2=0.0d+00

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf

! Get indices of first/last integrals to be passed for fixed ip,jp,iec1,jec1
            call index_GAM_3PK(nebf,npbf,
     x                         ip,jp,iec1,jec1,
     x                         1,1,1,1,ia_first)
            call index_GAM_3PK(nebf,npbf,
     x                         ip,jp,iec1,jec1,
     x                         nebf,nebf,nebf,nebf,ia_last)
! Determine how many messages (send+recv) will be passed by this MPI proc
! Only pass through MPI in case of block being split over >1 MPI proc
            if((mpistart.le.ia_first).and.(mpiend.ge.ia_last)) then
             nmsgs=2*(ia_last-ia_first+1)
             locblock=.true.
            else
             locblock=.false.
             if((mpistart.ge.ia_first).and.(mpistart.le.ia_last)) then
              if (mpiend.le.ia_last) then
               nmsgs=2*(mpiend-mpistart+1)
              else
               nmsgs=2*(ia_last-mpistart+1)
              end if
             else if((mpiend.ge.ia_first).and.(mpiend.le.ia_last)) then
              if (mpistart.le.ia_first) then
               nmsgs=2*(mpiend-ia_first+1)
              else
               nmsgs=2*(mpiend-mpistart+1)
              end if
             else
              nmsgs=1
              locblock=.true.
             end if
            end if
! Allocate and initialize arrays for MPI requests
            if(allocated(reqs)) deallocate(reqs)
            allocate(reqs(nmsgs))
            reqs=MPI_REQUEST_NULL
            reqscount=0

C            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C            if(rank.eq.0) write(*,1001) "index:",rank,ip,jp,iec1,jec1
C            write(*,*) "nmsgs:",rank,nmsgs
C            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec3,jec3,iec2,jec2,ia)
                             ia_132=ia

! ia_123 index is on this MPI process
       if ((ia_123.ge.mpistart).and.(ia_123.le.mpiend)) then

        if (locblock) then
         XGM3_2(ia_123-arrstart+1)=GM3_2(ia_132-arrstart+1)
        else
! MPI receive ia_ji integrals and store in X* arrs at ia_ji index
         call pack_4D(nebf,nebf,nebf,
     x                jec3,iec3,jec2,iec2,ia)
         tag_2=int(ia,kind=4)
         reqscount=reqscount+1

! Get ia_132 for GM3_2
         call get_mpi_proc(ng3,nproc,ia_132,recvrank)
         call MPI_IRECV(XGM3_2(ia_123-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,recvrank,tag_2,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

! ia_132 index is on this MPI process
       if ((ia_132.ge.mpistart).and.(ia_132.le.mpiend)) then

        if (.not.(locblock)) then
! Get rank of MPI process with ia_123 integrals
         call get_mpi_proc(ng3,nproc,ia_123,sendrank)

! MPI send ia_132 integrals for GM3_2
         call pack_4D(nebf,nebf,nebf,
     x                jec3,iec3,jec2,iec2,ia)
         tag_2=int(ia,kind=4)
         reqscount=reqscount+1
         call MPI_ISEND(GM3_2(ia_132-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,sendrank,tag_2,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

               end do
               end do
            end do
            end do

            if(allocated(statarr)) deallocate(statarr)
            allocate(statarr(MPI_STATUS_SIZE,nmsgs))
            call MPI_WAITALL(nmsgs,reqs,statarr,ierr)
            if (ierr.ne.0) write(*,*) "Trouble with reqs waitall"

         end do
         end do
      end do
      end do

      if(allocated(statarr)) deallocate(reqs)
      if(allocated(reqs)) deallocate(reqs)

! Symmetrize integrals locally
      do i=1,ng3loc

        x123=GM3_2(i)
        x132=XGM3_2(i)
        GM3_2(i)=(x123+x132)*half

      end do

! Done with GM3_2 integrals
      if(allocated(XGM3_2)) deallocate(XGM3_2)

! Get appropriate ia_ji values for XGM3_3
      if (rank.eq.0) write(*,*) "Symmetrizing GM3_3"

! Allocate storage for temporary arrays to store ia_ji integrals
      if(allocated(XGM3_3)) deallocate(XGM3_3)
      allocate(XGM3_3(ng3loc)) ! stores ia_213 integrals
      XGM3_3=0.0d+00

      do ip=1,npbf
      do jp=1,npbf
         do iec3=1,nebf
         do jec3=1,nebf

! Get indices of first/last integrals to be passed for fixed ip,jp,iec3,jec3
            call index_GAM_3PK(nebf,npbf,
     x                         ip,jp,1,1,
     x                         1,1,iec3,jec3,ia_first)
            call index_GAM_3PK(nebf,npbf,
     x                         ip,jp,nebf,nebf,
     x                         nebf,nebf,iec3,jec3,ia_last)
! Determine how many messages (send+recv) will be passed by this MPI proc
! Only pass through MPI in case of block being split over >1 MPI proc
            if((mpistart.le.ia_first).and.(mpiend.ge.ia_last)) then
             nmsgs=2*(ia_last-ia_first+1)
             locblock=.true.
            else
             locblock=.false.
             if((mpistart.ge.ia_first).and.(mpistart.le.ia_last)) then
              if (mpiend.le.ia_last) then
               nmsgs=2*(mpiend-mpistart+1)
              else
               nmsgs=2*(ia_last-mpistart+1)
              end if
             else if((mpiend.ge.ia_first).and.(mpiend.le.ia_last)) then
              if (mpistart.le.ia_first) then
               nmsgs=2*(mpiend-ia_first+1)
              else
               nmsgs=2*(mpiend-mpistart+1)
              end if
             else
              nmsgs=1
              locblock=.true.
             end if
            end if
! Allocate and initialize arrays for MPI requests
            if(allocated(reqs)) deallocate(reqs)
            allocate(reqs(nmsgs))
            reqs=MPI_REQUEST_NULL
            reqscount=0

C            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C            if(rank.eq.0) write(*,1001) "index:",rank,ip,jp,iec1,jec1
C            write(*,*) "nmsgs:",rank,nmsgs
C            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            do iec1=1,nebf
            do jec1=1,nebf
               do iec2=1,nebf
               do jec2=1,nebf

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,iec3,jec3,ia)
                             ia_213=ia

! ia_123 index is on this MPI process
       if ((ia_123.ge.mpistart).and.(ia_123.le.mpiend)) then

        if (locblock) then
         XGM3_3(ia_123-arrstart+1)=GM3_3(ia_213-arrstart+1)
        else
! MPI receive ia_ji integrals and store in X* arrs at ia_ji index
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_3=int(ia,kind=4)
         reqscount=reqscount+1

! Get ia_213 for GM3_3
        call get_mpi_proc(ng3,nproc,ia_213,recvrank)
        call MPI_IRECV(XGM3_3(ia_123-arrstart+1),1,
     x                 MPI_DOUBLE_PRECISION,recvrank,tag_3,
     x                 MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

! ia_213 index is on this MPI process
       if ((ia_213.ge.mpistart).and.(ia_213.le.mpiend)) then

        if (.not.(locblock)) then
! Get rank of MPI process with ia_123 integrals
         call get_mpi_proc(ng3,nproc,ia_123,sendrank)

! MPI send ia_213 integrals for GM3_3
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_3=int(ia,kind=4)
         reqscount=reqscount+1
         call MPI_ISEND(GM3_3(ia_213-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,sendrank,tag_3,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

               end do
               end do
            end do
            end do

            if(allocated(statarr)) deallocate(statarr)
            allocate(statarr(MPI_STATUS_SIZE,nmsgs))
            call MPI_WAITALL(nmsgs,reqs,statarr,ierr)
            if (ierr.ne.0) write(*,*) "Trouble with reqs waitall"

         end do
         end do
      end do
      end do

      if(allocated(statarr)) deallocate(reqs)
      if(allocated(reqs)) deallocate(reqs)

! Symmetrize integrals locally
      do i=1,ng3loc

        x123=GM3_3(i)
        x213=XGM3_3(i)
        GM3_3(i)=(x123+x213)*half

      end do

! Done with GM3_3 integrals
      if(allocated(XGM3_3)) deallocate(XGM3_3)

! Get appropriate ia_ji values for XGM3_4
      if (rank.eq.0) write(*,*) "Symmetrizing GM3_4"

! Allocate storage for temporary arrays to store ia_ji integrals
      if(allocated(XGM3_4)) deallocate(XGM3_4)
      allocate(XGM3_4(ng3loc)) ! stores ia_321 integrals
      XGM3_4=0.0d+00

      do ip=1,npbf
      do jp=1,npbf
         do iec2=1,nebf
         do jec2=1,nebf

! Get indices of first/last integrals to be passed for fixed ip,jp,iec2,jec2
            call index_GAM_3PK(nebf,npbf,
     x                         ip,jp,1,1,
     x                         iec2,iec2,1,1,ia_first)
            call index_GAM_3PK(nebf,npbf,
     x                         ip,jp,nebf,nebf,
     x                         iec2,iec2,nebf,nebf,ia_last)
! Determine how many messages (send+recv) will be passed by this MPI proc
! Only pass through MPI in case of block being split over >1 MPI proc
            if((mpistart.le.ia_first).and.(mpiend.ge.ia_last)) then
             nmsgs=2*(ia_last-ia_first+1)
             locblock=.true.
            else
             locblock=.false.
             if((mpistart.ge.ia_first).and.(mpistart.le.ia_last)) then
              if (mpiend.le.ia_last) then
               nmsgs=2*(mpiend-mpistart+1)
              else
               nmsgs=2*(ia_last-mpistart+1)
              end if
             else if((mpiend.ge.ia_first).and.(mpiend.le.ia_last)) then
              if (mpistart.le.ia_first) then
               nmsgs=2*(mpiend-ia_first+1)
              else
               nmsgs=2*(mpiend-mpistart+1)
              end if
             else
              nmsgs=1
              locblock=.true.
             end if
            end if
! Allocate and initialize arrays for MPI requests
            if(allocated(reqs)) deallocate(reqs)
            allocate(reqs(nmsgs))
            reqs=MPI_REQUEST_NULL
            reqscount=0

C            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C            if(rank.eq.0) write(*,1001) "index:",rank,ip,jp,iec2,jec2
C            write(*,*) "nmsgs:",rank,nmsgs
C            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            do iec1=1,nebf
            do jec1=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec3,jec3,iec2,jec2,iec1,jec1,ia)
                             ia_321=ia

! ia_123 index is on this MPI process
       if ((ia_123.ge.mpistart).and.(ia_123.le.mpiend)) then

        if (locblock) then
         XGM3_4(ia_123-arrstart+1)=GM3_4(ia_321-arrstart+1)
        else
! MPI receive ia_ji integrals and store in X* arrs at ia_ji index
         call pack_4D(nebf,nebf,nebf,
     x                jec3,iec3,jec1,iec1,ia)
         tag_4=int(ia,kind=4)
         reqscount=reqscount+1

! Get ia_321 for GM3_4
         call get_mpi_proc(ng3,nproc,ia_321,recvrank)
         call MPI_IRECV(XGM3_4(ia_123-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,recvrank,tag_4,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

! ia_321 index is on this MPI process
       if ((ia_321.ge.mpistart).and.(ia_321.le.mpiend)) then

        if (.not.(locblock)) then
! Get rank of MPI process with ia_123 integrals
         call get_mpi_proc(ng3,nproc,ia_123,sendrank)

! MPI send ia_321 integrals for GM3_4
         call pack_4D(nebf,nebf,nebf,
     x                jec3,iec3,jec1,iec1,ia)
         tag_4=int(ia,kind=4)
         reqscount=reqscount+1
         call MPI_ISEND(GM3_4(ia_321-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,sendrank,tag_4,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

               end do
               end do
            end do
            end do

            if(allocated(statarr)) deallocate(statarr)
            allocate(statarr(MPI_STATUS_SIZE,nmsgs))
            call MPI_WAITALL(nmsgs,reqs,statarr,ierr)
            if (ierr.ne.0) write(*,*) "Trouble with reqs waitall"

         end do
         end do
      end do
      end do

      if(allocated(statarr)) deallocate(reqs)
      if(allocated(reqs)) deallocate(reqs)

! Symmetrize integrals locally
      do i=1,ng3loc

        x123=GM3_4(i)
        x321=XGM3_4(i)
        GM3_4(i)=(x123+x321)*half

      end do

! Done with GM3_4 integrals
      if(allocated(XGM3_4)) deallocate(XGM3_4)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Symmetrize GM3_1 integrals
      if (rank.eq.0) write(*,*) "Symmetrizing GM3_1"

      call RXCHF_GAM3ex_MPI_symm(nproc,rank,
     x                           ng3,ng3loc,
     x                           mpistart,mpiend,arrstart,
     x                           nebf,npbf,
     x                           GM3_1)

!      write(*,*) "start,end,ng3loc:",mpistart,mpiend,ng3loc
!      do i=1,ng3loc
!       write(*,9001) GM3_1(i),GM3_2(i),
!     x               GM3_3(i),GM3_4(i)
!      end do

      wtime2 = MPI_WTIME() - wtime

      if(rank.eq.0) write(*,3000) 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      write(*,3001) rank,wtime2
!--------------------SYMMETRIZE----------------------------------------)

! Construct global arrays on master process for testing
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if(allocated(TGM3_1)) deallocate(TGM3_1)
!      if(allocated(TGM3_2)) deallocate(TGM3_2)
!      if(allocated(TGM3_3)) deallocate(TGM3_3)
!      if(allocated(TGM3_4)) deallocate(TGM3_4)
!      if (rank.eq.0) then
!       allocate(TGM3_1(ng3))
!       allocate(TGM3_2(ng3))
!       allocate(TGM3_3(ng3))
!       allocate(TGM3_4(ng3))
!      else
!       allocate(TGM3_1(1))
!       allocate(TGM3_2(1))
!       allocate(TGM3_3(1))
!       allocate(TGM3_4(1))
!      end if
!      TGM3_1=zero
!      TGM3_2=zero
!      TGM3_3=zero
!      TGM3_4=zero
!
!      ng3loc4=int(ng3loc,kind=4)
!
!! Get number of elements calculated by each proc
!      call MPI_GATHER(ng3loc4,1,MPI_INTEGER,
!     x                ng3locarr(1),1,MPI_INTEGER,
!     x                0,MPI_COMM_WORLD,ierr)
!
!! Get displacements for array storage
!      if (rank.eq.0) then
!        displarr(1)=0
!        do i=2,nproc
!          displarr(i)=displarr(i-1)+ng3locarr(i-1)
!        end do
!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM3_1 on root
!      call MPI_GATHERV(GM3_1(1),ng3loc,MPI_DOUBLE_PRECISION,
!     x                 TGM3_1(1),ng3locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM3_2 on root
!      call MPI_GATHERV(GM3_2(1),ng3loc,MPI_DOUBLE_PRECISION,
!     x                 TGM3_2(1),ng3locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM3_3 on root
!      call MPI_GATHERV(GM3_3(1),ng3loc,MPI_DOUBLE_PRECISION,
!     x                 TGM3_3(1),ng3locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM3_4 on root
!      call MPI_GATHERV(GM3_4(1),ng3loc,MPI_DOUBLE_PRECISION,
!     x                 TGM3_4(1),ng3locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!!      if (rank.eq.0) then
!!       write(*,*) "concatenated ng3"
!!       do i=1,ng3
!!        write(*,9001) TGM3_1(i),TGM3_2(i),
!!     x                TGM3_3(i),TGM3_4(i)
!!       end do
!!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if (rank.eq.0) then
!       open(unit=20,file="XCHF_GAM3.ufm",form="unformatted")
!       write(20) TGM3_1
!       close(20)
!       write(*,*) "XCHF_GAM3 written to disk"
!       open(unit=21,file="INT_GAM3.ufm",form="unformatted")
!       write(21) TGM3_2
!       close(21)
!       write(*,*) "INT_GAM3 written to disk"
!       open(unit=22,file="INT_GAM3ex1.ufm",form="unformatted")
!       write(22) TGM3_3
!       close(22)
!       write(*,*) "INT_GAM3ex1 written to disk"
!       open(unit=23,file="INT_GAM3ex2.ufm",form="unformatted")
!       write(23) TGM3_4
!       close(23)
!       write(*,*) "INT_GAM3ex2 written to disk"
!      end if
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!      if(allocated(TGM3_4)) deallocate(TGM3_4)
!      if(allocated(TGM3_3)) deallocate(TGM3_3)
!      if(allocated(TGM3_2)) deallocate(TGM3_2)
!      if(allocated(TGM3_1)) deallocate(TGM3_1)

 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 4-PARTICLE INTEGRALS        |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 4-PARTICLE INTEGRALS: ',1X,I12/
     x        8X,'  NUMBER OF BLOCKS (USER DEFINED): ',1X,I12/
     x        8X,'                          ',/,
     x        8X,'  COMPUTATIONAL RESOURCES:',/,
     x        8X,'  ------------------------',/)

 1500 FORMAT( 8X,'      MPI PROCESSES:',1X,I3/
     x        8X,'        OMP THREADS:',1X,I3/)

 2000 FORMAT(/8X,'  INTEGRAL CALCULATION TIMINGS:',/,
     x        8X,'  -----------------------------')

 2001 FORMAT( 8X,'    PROCESS ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,' INTEGRAL SYMMETRIZATION TIMINGS:',/,
     x        8X,' --------------------------------')

 3001 FORMAT( 8X,'    PROCESS ',1X,I4,1X,F10.2)

 9001 FORMAT(1X,4(F20.10))

      return
      end
C=======================================================================
      subroutine RXCHF_GAM3_thread_MPI(istart,iend,ng3_seg,ng3,
     x                          nebf,npebf,npbf,nat,ngtg1,
     x                          pmass,cat,zan,bcoef1,gamma1,
     x                          loop_map,arrstart,GM3_1,GM3_2,
     x                          KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                          ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input Variables
      integer istart,iend,ng3_seg
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas
      integer ng3      ! Number of integrals calc by MPI process
      integer arrstart ! Index of first integral

C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
      integer KPESTR(nebf)  ! Map contracted index to primitive start
      integer KPEEND(nebf)  ! Map contracted index to primitive end
C-------Basis Set Info-------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)
      integer loop_map(ng3_seg,8)

! Variables Returned
      double precision GM3_1(ng3),GM3_2(ng3)

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer iec3,jec3  !
      integer imap,ia
      double precision OMG3_1,OMG3_2

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(istart,iend)
!$ompx shared(loop_map,arrstart)
!$ompx shared(GM3_1,GM3_2)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(OMG3_1,OMG3_2)
!$ompx private(ia) 

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

         call RXCHFmult_contract_omega3_conv(ip,jp,iec1,jec1,iec2,jec2,
     x                            iec3,jec3,nebf,npebf,npbf,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            OMG3_1,OMG3_2)

         call index_GAM_3PK(nebf,npbf,
     x          ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)

         GM3_1(ia-arrstart+1)=OMG3_1
         GM3_2(ia-arrstart+1)=OMG3_2

      end do
!$omp end do
!$omp end parallel      
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

      return
      end
C=======================================================================
      subroutine RXCHF_GAM3ex_thread_MPI(istart,iend,ng3_seg,ng3,
     x                            nebf,npebf,npbf,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            loop_map,arrstart,
     x                            GM3_1,GM3_2,GM3_3,GM3_4,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input Variables
      integer istart,iend,ng3_seg
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas
      integer ng3      ! Number of integrals calc by MPI process
      integer arrstart ! Index of first integral

C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
      integer KPESTR(nebf)  ! Map contracted index to primitive start
      integer KPEEND(nebf)  ! Map contracted index to primitive end
C-------Basis Set Info-------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)
      integer loop_map(ng3_seg,8)

! Variables Returned
      double precision GM3_1(ng3),GM3_2(ng3)
      double precision GM3_3(ng3),GM3_4(ng3)

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer iec3,jec3  !
      integer imap,ia
      double precision OMG3_1,OMG3_2,OMG3_3,OMG3_4

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(istart,iend)
!$ompx shared(loop_map,arrstart)
!$ompx shared(GM3_1,GM3_2,GM3_3,GM3_4)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(OMG3_1,OMG3_2,OMG3_3,OMG3_4)
!$ompx private(ia) 

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

        call RXCHFmult_contract_omega3_convex(ip,jp,iec1,jec1,iec2,jec2,
     x                            iec3,jec3,nebf,npebf,npbf,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            OMG3_1,OMG3_2,OMG3_3,OMG3_4)

         call index_GAM_3PK(nebf,npbf,
     x          ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)

         GM3_1(ia-arrstart+1)=OMG3_1
         GM3_2(ia-arrstart+1)=OMG3_2
         GM3_3(ia-arrstart+1)=OMG3_3
         GM3_4(ia-arrstart+1)=OMG3_4

      end do
!$omp end do
!$omp end parallel      
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

      return
      end
C=======================================================================
      subroutine RXCHF_GAM3ex_MPI_symm(nproc,rank,
     x                                 ng3,ng3loc,
     x                                 mpistart,mpiend,arrstart,
     x                                 nebf,npbf,
     x                                 GM3_1)

C=======================================================================
      implicit none
      include 'mpif.h'

C Input Variables
      integer nproc     ! Number of MPI processes
      integer rank      ! Rank of this MPI proc
      integer ng3       ! Total number of integrals
      integer ng3loc    ! Number of integrals for MPI proc to calc
      integer mpistart  ! Start index covered by this MPI proc
      integer mpiend    ! End index covered by this MPI proc
      integer arrstart  ! Index shift for this MPI proc
      integer nebf      ! Number contracted electronic basis functions
      integer npbf      ! Number nuclear basis functions
      integer*4 ierr

! Output Variables
      double precision GM3_1(ng3loc)  ! On output, symmetrized integrals

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer iec3,jec3  !
      integer sendrank,recvrank

      double precision XGM3_1(ng3loc)    ! Stores reduction of integrals
      double precision XGM3_1aux(ng3loc) ! Stores individual set of integrals

      integer tag_1
      integer*4 reqs_1(2*ng3loc)

      integer i,ia
      integer ia_123
      integer ia_132
      integer ia_213
      integer ia_231
      integer ia_312
      integer ia_321

      double precision x123
      double precision xall

      double precision zero,six
      parameter(zero=0.0d+00,six=6.0d+00)

!--------------------SYMMETRIZE----------------------------------------(
C Symmetrized integrals in GM3_1ICR (XCHF integrals)
C  - completed using six passes to reduce memory requirements

      XGM3_1=zero

! On first pass, get ia_132 set for XGM3_1
!      if (rank.eq.0) write(*,*) "   ... starting first pass"
      XGM3_1aux=zero
      reqs_1=MPI_REQUEST_NULL

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec3,jec3,iec2,jec2,ia)
                             ia_132=ia

       if ((ia_123.ge.mpistart).and.(ia_123.le.mpiend)) then
        tag_1=ia_123
        call get_mpi_proc(ng3,nproc,ia_132,recvrank)
        call MPI_IRECV(XGM3_1aux(ia_123-arrstart+1),1,
     x                 MPI_DOUBLE_PRECISION,
     x                 recvrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_123-arrstart+1),ierr)
       end if

       if ((ia_132.ge.mpistart).and.(ia_132.le.mpiend)) then
        call get_mpi_proc(ng3,nproc,ia_123,sendrank)
        tag_1=ia_123
        call MPI_ISEND(GM3_1(ia_132-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 sendrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_132-arrstart+1+ng3loc),ierr)
       end if

               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      call MPI_WAITALL(2*ng3loc,reqs_1,MPI_STATUSES_IGNORE,ierr)
      if (ierr.ne.0) write(*,*) "Trouble with ia_132 waitall"

! Add ia_132 data to XGM3_1
      call add_to_arr(ng3loc,XGM3_1aux,XGM3_1)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! On second pass, get ia_213 set for XGM3_1
!      if (rank.eq.0) write(*,*) "   ... starting second pass"
      XGM3_1aux=zero
      reqs_1=MPI_REQUEST_NULL

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,iec3,jec3,ia)
                             ia_213=ia

       if ((ia_123.ge.mpistart).and.(ia_123.le.mpiend)) then
        tag_1=ia_123
        call get_mpi_proc(ng3,nproc,ia_213,recvrank)
        call MPI_IRECV(XGM3_1aux(ia_123-arrstart+1),1,
     x                 MPI_DOUBLE_PRECISION,
     x                 recvrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_123-arrstart+1),ierr)
       end if

       if ((ia_213.ge.mpistart).and.(ia_213.le.mpiend)) then
        call get_mpi_proc(ng3,nproc,ia_123,sendrank)
        tag_1=ia_123
        call MPI_ISEND(GM3_1(ia_213-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 sendrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_213-arrstart+1+ng3loc),ierr)
       end if

               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      call MPI_WAITALL(2*ng3loc,reqs_1,MPI_STATUSES_IGNORE,ierr)
      if (ierr.ne.0) write(*,*) "Trouble with ia_213 waitall"

! Add ia_213 data to XGM3_1
      call add_to_arr(ng3loc,XGM3_1aux,XGM3_1)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! On third pass, get ia_231 set for XGM3_1
!      if (rank.eq.0) write(*,*) "   ... starting third pass"
      XGM3_1aux=zero
      reqs_1=MPI_REQUEST_NULL

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec3,jec3,iec1,jec1,ia)
                             ia_231=ia

       if ((ia_123.ge.mpistart).and.(ia_123.le.mpiend)) then
        tag_1=ia_123
        call get_mpi_proc(ng3,nproc,ia_231,recvrank)
        call MPI_IRECV(XGM3_1aux(ia_123-arrstart+1),1,
     x                 MPI_DOUBLE_PRECISION,
     x                 recvrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_123-arrstart+1),ierr)
       end if

       if ((ia_231.ge.mpistart).and.(ia_231.le.mpiend)) then
        call get_mpi_proc(ng3,nproc,ia_123,sendrank)
        tag_1=ia_123
        call MPI_ISEND(GM3_1(ia_231-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 sendrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_231-arrstart+1+ng3loc),ierr)
       end if

               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      call MPI_WAITALL(2*ng3loc,reqs_1,MPI_STATUSES_IGNORE,ierr)
      if (ierr.ne.0) write(*,*) "Trouble with ia_231 waitall"

! Add ia_231 data to XGM3_1
      call add_to_arr(ng3loc,XGM3_1aux,XGM3_1)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! On fourth pass, get ia_312 set for XGM3_1
!      if (rank.eq.0) write(*,*) "   ... starting fourth pass"
      XGM3_1aux=zero
      reqs_1=MPI_REQUEST_NULL

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec3,jec3,iec1,jec1,iec2,jec2,ia)
                             ia_312=ia

       if ((ia_123.ge.mpistart).and.(ia_123.le.mpiend)) then
        tag_1=ia_123
        call get_mpi_proc(ng3,nproc,ia_312,recvrank)
        call MPI_IRECV(XGM3_1aux(ia_123-arrstart+1),1,
     x                 MPI_DOUBLE_PRECISION,
     x                 recvrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_123-arrstart+1),ierr)
       end if

       if ((ia_312.ge.mpistart).and.(ia_312.le.mpiend)) then
        call get_mpi_proc(ng3,nproc,ia_123,sendrank)
        tag_1=ia_123
        call MPI_ISEND(GM3_1(ia_312-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 sendrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_312-arrstart+1+ng3loc),ierr)
       end if

               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      call MPI_WAITALL(2*ng3loc,reqs_1,MPI_STATUSES_IGNORE,ierr)
      if (ierr.ne.0) write(*,*) "Trouble with ia_312 waitall"

! Add ia_312 data to XGM3_1
      call add_to_arr(ng3loc,XGM3_1aux,XGM3_1)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! On fifth pass, get ia_321 set for XGM3_1
!      if (rank.eq.0) write(*,*) "   ... starting fifth pass"
      XGM3_1aux=zero
      reqs_1=MPI_REQUEST_NULL

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec3,jec3,iec2,jec2,iec1,jec1,ia)
                             ia_321=ia

       if ((ia_123.ge.mpistart).and.(ia_123.le.mpiend)) then
        tag_1=ia_123
        call get_mpi_proc(ng3,nproc,ia_321,recvrank)
        call MPI_IRECV(XGM3_1aux(ia_123-arrstart+1),1,
     x                 MPI_DOUBLE_PRECISION,
     x                 recvrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_123-arrstart+1),ierr)
       end if

       if ((ia_321.ge.mpistart).and.(ia_321.le.mpiend)) then
        call get_mpi_proc(ng3,nproc,ia_123,sendrank)
        tag_1=ia_123
        call MPI_ISEND(GM3_1(ia_321-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 sendrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_321-arrstart+1+ng3loc),ierr)
       end if

               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      call MPI_WAITALL(2*ng3loc,reqs_1,MPI_STATUSES_IGNORE,ierr)
      if (ierr.ne.0) write(*,*) "Trouble with ia_321 waitall"

! Add ia_321 data to XGM3_1
      call add_to_arr(ng3loc,XGM3_1aux,XGM3_1)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Symmetrize GM3_1 locally
      do i=1,ng3loc
        x123=GM3_1(i)
        xall=XGM3_1(i)
        GM3_1(i)=(x123+xall)/six
      end do

      return
      end


