!=======================================================================
      subroutine RXCHF_GAM2_MPI(nproc,rank,
     x                          Nchunks,nebf,npebf,npbf,
     x                          ng2,ng2loc,ng2prm,nat,ngtg1,
     x                          pmass,cat,zan,bcoef1,gamma1,
     x                          KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                          ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                          GM2_1,GM2_2,GM2s)

!=======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ng2             ! Total number of integrals
      integer ng2loc          ! Number of integrals for MPI proc to calc
      integer nebf,npebf,npbf,ng2prm
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
      double precision GM2_1(ng2loc)  ! XCHF OMG2   integrals (symm)
      double precision GM2_2(ng2loc)  ! INT  OMG2   integrals (unsymm)
      double precision GM2s(ng2loc)   ! XCHF OMG2s  integrals (symm)

! Local Variables
      integer istat,ichunk,istart,iend,ng2_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,i,j
      integer,allocatable :: loop_map(:,:)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

      integer*4 tag_1,tag_s
      integer*4 sendrank,recvrank
C      integer*4,allocatable :: req(:)
C      integer*4,pointer :: req(:)
      integer*4, allocatable :: reqs(:)
C      integer*4, allocatable, target :: reqs(:)
      integer ia_first,ia_last,nmsgs,reqscount
      logical flag
      logical locblock

      double precision, allocatable :: XGM2_1(:)
      double precision, allocatable :: XGM2s(:)

      integer ia
      integer ia_12
      integer ia_21

      double precision x12,x21

      double precision zero,half,two
      parameter(zero=0.0d+00,half=0.5d+00,two=2.0d+00)

      double precision wtime
      double precision wtime2

! Testing Variables
!      double precision, allocatable :: TGM2_1(:) ! Testing arrays
!      double precision, allocatable :: TGM2_2(:)
!      double precision, allocatable :: TGM2s(:)
!
!      integer*4 ng2loc4
!      integer*4 ng2locarr(nproc),displarr(nproc)


! Have each process calculate ng2/nproc integrals according to rank
! Have last process calculate ng2%nproc remaining integrals
      call get_mpi_range(ng2,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng2

      if (rank.eq.0) then
       write(*,1000) ng2,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      GM2_1=0.0d+00
      GM2_2=0.0d+00
      GM2s=0.0d+00

!-----CHOP-UP-THE-CALCULATION-OF-GAM_2--------------------------------(
      wtime = MPI_WTIME()

      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng2/nproc integrals
            call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
            ng2_seg=1+iend-istart

            if(allocated(loop_map)) deallocate(loop_map)
            allocate( loop_map(ng2_seg,6),stat=istat )

! Nested loop compression for this chunk:
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

! Save coordinates of first value for array transfer later
                         if ((ichunk.eq.1).and.(Loopi.eq.1)) then
                            call index_GAM_2PK(nebf,npbf,
     x                               ip,jp,iec1,jec1,iec2,jec2,arrstart)
                         end if

                      end if

                  end do
                  end do
               end do
               end do
            end do
            end do

         call RXCHF_GAM2_thread_MPI(istart,iend,ng2_seg,ng2loc,
     x                         nebf,npebf,npbf,nat,ngtg1,
     x                         pmass,cat,zan,bcoef1,gamma1,
     x                         loop_map,arrstart,
     x                         GM2_1,GM2_2,GM2s,
     x                         KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                         ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_2--------------------------------)
!-----CLEAN-UP-MEMORY-------------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

      wtime2 = MPI_WTIME() - wtime

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(rank.eq.0) write(*,2000) 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      write(*,2001) rank,wtime2

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng2loc:",mpistart,mpiend,ng2loc
!      do i=1,ng2loc
!       write(*,9001) GM2_1(i),GM2_2(i),
!     x               GM2s(i)
!      end do

!--------------------SYMMETRIZE----------------------------------------(
C Symmetrized integrals in GM2_1ICR (XCHF integrals)
C Unsymmetrized integrals in GM2_2ICR (interaction integrals)
C Symmetrized integrals in GM2sICR (XCHF integrals)

      wtime = MPI_WTIME() 

! Symmetrize GM2_1 integrals

! Allocate storage for temporary arrays to store ia_21 integrals
      if(allocated(XGM2_1)) deallocate(XGM2_1)
      allocate(XGM2_1(ng2loc))
      XGM2_1=0.0d+00

! Get appropriate ia_ji values for XGM2_1
      if (rank.eq.0) write(*,*) "Symmetrizing GM2_1"

      do ip=1,npbf
      do jp=1,npbf

! Get indices of first/last integrals to be passed for fixed ip,jp
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,1,1,1,1,ia_first)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,nebf,nebf,nebf,nebf,ia_last)
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
C         if(locblock) then
C          write(*,*) "locblock",rank,ia_first,ia_last,
C     x               mpistart,mpiend,nmsgs
C         else
C          write(*,*) "not locblock,first,last:",rank,ia_first,ia_last,
C     x               mpistart,mpiend,nmsgs
C         end if

C         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C         if(rank.eq.0) write(*,1001) "index:",rank,ip,jp
C         write(*,*) "nmsgs:",rank,nmsgs
C         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf

         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia_12)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,ia_21)

! ia_12 index is on this MPI process
       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then

        if (locblock) then
         XGM2_1(ia_12-arrstart+1)=GM2_1(ia_21-arrstart+1)
        else
! MPI receive ia_21 integrals and store in X* arrs at ia_12 index
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_1=int(ia,kind=4)
         reqscount=reqscount+1

! Get ia_21 for GM2_1
         call get_mpi_proc(ng2,nproc,ia_21,recvrank)
         call MPI_IRECV(XGM2_1(ia_12-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,recvrank,tag_1,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
C         write(*,*) "add req",rank,reqscount,reqs(reqscount),"recv"
        end if

       end if

! ia_21 index is on this MPI process
       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then

        if (.not.(locblock)) then
! Get rank of MPI process with ia_12 integrals
         call get_mpi_proc(ng2,nproc,ia_12,sendrank)

! MPI send ia_21 integrals for GM2_1
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_1=int(ia,kind=4)
         reqscount=reqscount+1
         call MPI_ISEND(GM2_1(ia_21-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,sendrank,tag_1,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
C         write(*,*) "add req",rank,reqscount,reqs(reqscount),"send"

        end if

       end if

            end do
            end do
         end do
         end do

         if (.not.locblock) then
C          call MPI_TESTALL(nmsgs,reqs,flag,MPI_STATUSES_IGNORE,ierr)
C          write(*,*) "testall flag:",rank,flag
C          call MPI_WAITALL(nmsgs,reqs,MPI_STATUSES_IGNORE,ierr)
C          if (ierr.ne.0) write(*,*) "Trouble with GM2_1 waitall"

C          write(*,*) "nmsgs,reqscount:",rank,nmsgs,reqscount
C          allocate(req(nmsgs),stat=allocstat)
C          write(*,*) "stat:",rank,allocstat
C          do i=1,nmsgs
C            do j=1,nmsgs
C              req(j)=reqs(j)
C            end do
C            write(*,*) "loop:",rank,i
CC            if(allocated(req)) deallocate(req)
CC            req(:)=reqs(:)
C            write(*,*) "i,reqs 1:",rank,i,reqs(1),reqs(2),reqs(3)
C            write(*,*) "i,req 1:",rank,i,req(1),req(2),req(3)
C            if (req(i).ne.MPI_REQUEST_NULL) then
CC             call MPI_WAIT(req,stat,ierr)
CC           write(*,*) "allocated 1:",rank,allocated(reqs),allocated(req)
C             call MPI_TEST(req(i),flag,stat1,ierr)
CC           write(*,*) "allocated 2:",rank,allocated(reqs),allocated(req)
C             write(*,*) "i,flag:",rank,i,flag
C             write(*,*) "i,req 2:",rank,i,req(1),req(2),req(3)
C             write(*,*) "i,reqs 2:",rank,i,reqs(1),reqs(2),reqs(3)
C            end if
C          end do

C          write(*,*) "after all:",rank!,reqs
          call RXCHF_MPI_WAITALL2(rank,nmsgs,reqs)
         end if

      end do
      end do

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      if(allocated(XGM2_1)) deallocate(XGM2_1)
C      RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C      write(*,*) "after loop:",rank
C      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C      if(allocated(reqs)) deallocate(reqs)
C      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C      write(*,*) "after loop2:",rank
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Symmetrize integrals locally
      do i=1,ng2loc

        x12=GM2_1(i)
        x21=XGM2_1(i)
        GM2_1(i)=(x12+x21)/two

      end do

! Done with GM2_1 integrals
      if(allocated(XGM2_1)) deallocate(XGM2_1)

! Symmetrize GM2s integrals

! Allocate storage for temporary arrays to store ia_21 integrals
      if(allocated(XGM2s)) deallocate(XGM2s)
      allocate(XGM2s(ng2loc))
      XGM2s=0.0d+00

! Get appropriate ia_ji values for XGM2s
      if (rank.eq.0) write(*,*) "Symmetrizing GM2s"

      do ip=1,npbf
      do jp=1,npbf

! Get indices of first/last integrals to be passed for fixed ip,jp
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,1,1,1,1,ia_first)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,nebf,nebf,nebf,nebf,ia_last)
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

C         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C         if(rank.eq.0) write(*,1001) "index:",rank,ip,jp
C         write(*,*) "nmsgs:",rank,nmsgs
C         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf

         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia_12)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,ia_21)

! ia_12 index is on this MPI process
       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then

        if (locblock) then
         XGM2s(ia_12-arrstart+1)=GM2s(ia_21-arrstart+1)
        else
! MPI receive ia_21 integrals and store in X* arrs at ia_12 index
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_s=int(ia,kind=4)
         reqscount=reqscount+1

! Get ia_21 for GM2s
         call get_mpi_proc(ng2,nproc,ia_21,recvrank)
         call MPI_IRECV(XGM2s(ia_12-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,recvrank,tag_s,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

! ia_21 index is on this MPI process
       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then

        if (.not.(locblock)) then
! Get rank of MPI process with ia_12 integrals
         call get_mpi_proc(ng2,nproc,ia_12,sendrank)

! MPI send ia_21 integrals for GM2s
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_s=int(ia,kind=4)
         reqscount=reqscount+1
         call MPI_ISEND(GM2s(ia_21-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,sendrank,tag_s,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)

        end if

       end if

            end do
            end do
         end do
         end do

         if (.not.locblock) then
CC          call MPI_WAITALL(nmsgs,reqs,MPI_STATUSES_IGNORE,ierr)
CC          if (ierr.ne.0) write(*,*) "Trouble with GM2s waitall"
C          flag=.false.
C          write(*,*) "nmsgs,reqscount:",rank,nmsgs,reqscount
C          do i=1,nmsgs
C            write(*,*) "loop:",rank,i
CC            if(allocated(req)) deallocate(req)
C            allocate(req(nmsgs),stat=allocstat)
C            write(*,*) "stat:",rank,i,allocstat
CC            req=reqs(i)
C            req(:)=reqs(:)
C            write(*,*) "i,reqs 1:",rank,i,reqs(1),reqs(2),reqs(3)
C            write(*,*) "i,req 1:",rank,i,req(1),req(2),req(3)
C            if (req(i).ne.MPI_REQUEST_NULL) then
CC            call MPI_WAIT(req,stat,ierr)
CC           write(*,*) "allocated 1:",rank,allocated(reqs),allocated(req)
CC             call MPI_TEST(req(i),flag,stat,ierr)
CC           write(*,*) "allocated 2:",rank,allocated(reqs),allocated(req)
C             write(*,*) "i,flag:",rank,i,flag
CC             write(*,*) "i,reqs 2:",rank,i,reqs(1),reqs(2),reqs(3)
CC             write(*,*) "i,req 2:",rank,i,req(1),req(2),req(3)
C            end if
C          end do
C          write(*,*) "after all:",rank!,reqs
          call RXCHF_MPI_WAITALL2(rank,nmsgs,reqs)
         end if

      end do
      end do

      if(allocated(reqs)) deallocate(reqs)

! Symmetrize integrals locally
      do i=1,ng2loc

        x12=GM2s(i)
        x21=XGM2s(i)
        GM2s(i)=(x12+x21)/two

      end do

! Done with GM2s integrals
      if(allocated(XGM2s)) deallocate(XGM2s)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng2loc:",mpistart,mpiend,ng2loc
!      do i=1,ng2loc
!       write(*,9001) GM2_1(i),GM2_2(i),
!     x               GM2s(i)
!      end do

      wtime2 = MPI_WTIME() - wtime

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(rank.eq.0) write(*,3000) 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      write(*,3001) rank,wtime2
!--------------------SYMMETRIZE----------------------------------------)

! Construct global array on master process for testing
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if(allocated(TGM2_1)) deallocate(TGM2_1)
!      if(allocated(TGM2_2)) deallocate(TGM2_2)
!      if(allocated(TGM2s)) deallocate(TGM2s)
!      if (rank.eq.0) then
!       allocate(TGM2_1(ng2))
!       allocate(TGM2_2(ng2))
!       allocate(TGM2s(ng2))
!      else
!       allocate(TGM2_1(1))
!       allocate(TGM2_2(1))
!       allocate(TGM2s(1))
!      end if
!      TGM2_1=zero
!      TGM2_2=zero
!      TGM2s=zero
!
!      ng2loc4=int(ng2loc,kind=4)
!
!! Get number of elements calculated by each proc
!      call MPI_GATHER(ng2loc4,1,MPI_INTEGER,
!     x                ng2locarr(1),1,MPI_INTEGER,
!     x                0,MPI_COMM_WORLD,ierr)
!
!! Get displacements for array storage
!      if (rank.eq.0) then
!        displarr(1)=0
!        do i=2,nproc
!          displarr(i)=displarr(i-1)+ng2locarr(i-1)
!        end do
!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM2_1 on root
!      call MPI_GATHERV(GM2_1(1),ng2loc,MPI_DOUBLE_PRECISION,
!     x                 TGM2_1(1),ng2locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM2_2 on root
!      call MPI_GATHERV(GM2_2(1),ng2loc,MPI_DOUBLE_PRECISION,
!     x                 TGM2_2(1),ng2locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM2s on root
!      call MPI_GATHERV(GM2s(1),ng2loc,MPI_DOUBLE_PRECISION,
!     x                 TGM2s(1),ng2locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!!      if (rank.eq.0) then
!!       write(*,*) "concatenated ng2"
!!       do i=1,ng2
!!        write(*,9001) TGM2_1(i),TGM2_2(i),
!!     x                TGM2s(i)
!!       end do
!!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if (rank.eq.0) then
!       open(unit=20,file="XCHF_GAM2.ufm",form="unformatted")
!       write(20) TGM2_1
!       close(20)
!       write(*,*) "XCHF_GAM2 written to disk"
!       open(unit=21,file="INT_GAM2.ufm",form="unformatted")
!       write(21) TGM2_2
!       close(21)
!       write(*,*) "INT_GAM2 written to disk"
!       open(unit=22,file="XCHF_GAM2s.ufm",form="unformatted")
!       write(22) TGM2s
!       close(22)
!       write(*,*) "XCHF_GAM2s written to disk"
!      end if
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!      if(allocated(TGM2s)) deallocate(TGM2s)
!      if(allocated(TGM2_2)) deallocate(TGM2_2)
!      if(allocated(TGM2_1)) deallocate(TGM2_1)

 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 3-PARTICLE INTEGRALS        |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 3-PARTICLE INTEGRALS: ',1X,I12/
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

 9001 FORMAT(1X,3(F20.10))

      return
      end
!=======================================================================
      subroutine RXCHF_GAM2ex_MPI(nproc,rank,
     x                            Nchunks,nebf,npebf,npbf,
     x                            ng2,ng2loc,ng2prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            GM2_1,GM2_2,GM2_3,GM2s)

!=======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ng2             ! Total number of integrals
      integer ng2loc          ! Number of integrals for MPI proc to calc
      integer nebf,npebf,npbf,ng2prm
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
      double precision GM2_1(ng2loc)  ! XCHF OMG2   integrals (symm)
      double precision GM2_2(ng2loc)  ! INT  OMG2   integrals (unsymm)
      double precision GM2_3(ng2loc)  ! INT  OMG2ex integrals (symm)
      double precision GM2s(ng2loc)   ! XCHF OMG2s  integrals (symm)

! Local Variables
      integer istat,ichunk,istart,iend,ng2_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,i
      integer,allocatable :: loop_map(:,:)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

      integer*4 tag_1,tag_3,tag_s
      integer*4 sendrank,recvrank
      integer*4, allocatable :: reqs(:)
      integer ia_first,ia_last,nmsgs,reqscount
      logical locblock

      double precision, allocatable :: XGM2_1(:)
      double precision, allocatable :: XGM2_3(:)
      double precision, allocatable :: XGM2s(:)

      integer ia
      integer ia_12
      integer ia_21

      double precision x12,x21

      double precision zero,half,two
      parameter(zero=0.0d+00,half=0.5d+00,two=2.0d+00)

      double precision wtime
      double precision wtime2

! Testing Variables
!      double precision, allocatable :: TGM2_1(:) ! Testing arrays
!      double precision, allocatable :: TGM2_2(:)
!      double precision, allocatable :: TGM2_3(:)
!      double precision, allocatable :: TGM2s(:)
!
!      integer*4 ng2loc4
!      integer*4 ng2locarr(nproc),displarr(nproc)


! Have each process calculate ng2/nproc integrals according to rank
! Have last process calculate ng2%nproc remaining integrals
      call get_mpi_range(ng2,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ng2

      if (rank.eq.0) then
       write(*,1000) ng2,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      GM2_1=0.0d+00
      GM2_2=0.0d+00
      GM2_3=0.0d+00
      GM2s=0.0d+00

!-----CHOP-UP-THE-CALCULATION-OF-GAM_2--------------------------------(
      wtime = MPI_WTIME()

      do ichunk=1,Nchunks

! Have threads chop calculation of mpiend-mpistart+1=ng2/nproc integrals
            call loop_size(mpistart,mpiend,Nchunks,ichunk-1,istart,iend)
            ng2_seg=1+iend-istart

            if(allocated(loop_map)) deallocate(loop_map)
            allocate( loop_map(ng2_seg,6),stat=istat )

! Nested loop compression for this chunk:
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
     x                               ip,jp,iec1,jec1,iec2,jec2,arrstart)
                         end if

                      end if

                  end do
                  end do
               end do
               end do
            end do
            end do

         call RXCHF_GAM2ex_thread_MPI(istart,iend,ng2_seg,ng2loc,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           loop_map,arrstart,
     x                           GM2_1,GM2_2,GM2_3,GM2s,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_2--------------------------------)
!-----CLEAN-UP-MEMORY-------------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

      wtime2 = MPI_WTIME() - wtime

      if(rank.eq.0) write(*,2000) 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      write(*,2001) rank,wtime2

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng2loc:",mpistart,mpiend,ng2loc
!      do i=1,ng2loc
!       write(*,9001) GM2_1(i),GM2_2(i),
!     x               GM2_3(i),GM2s(i)
!      end do

!--------------------SYMMETRIZE----------------------------------------(
C Symmetrized integrals in GM2_1ICR (XCHF integrals)
C Unsymmetrized integrals in GM2_2ICR (interaction integrals)
C Symmetrized integrals in GM2_3ICR (exchange integrals)
C Symmetrized integrals in GM2sICR (XCHF integrals)

      wtime = MPI_WTIME() 

! Symmetrize GM2_1 integrals

! Allocate storage for temporary arrays to store ia_21 integrals
      if(allocated(XGM2_1)) deallocate(XGM2_1)
      allocate(XGM2_1(ng2loc))
      XGM2_1=0.0d+00

! Get appropriate ia_ji values for XGM2_1
      if (rank.eq.0) write(*,*) "Symmetrizing GM2_1"

      do ip=1,npbf
      do jp=1,npbf

! Get indices of first/last integrals to be passed for fixed ip,jp
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,1,1,1,1,ia_first)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,nebf,nebf,nebf,nebf,ia_last)
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

C         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C         if(rank.eq.0) write(*,1001) "index:",rank,ip,jp
C         write(*,*) "nmsgs:",rank,nmsgs
C         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf

         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia_12)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,ia_21)

! ia_12 index is on this MPI process
       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then

        if (locblock) then
         XGM2_1(ia_12-arrstart+1)=GM2_1(ia_21-arrstart+1)
        else
! MPI receive ia_21 integrals and store in X* arrs at ia_12 index
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_1=int(ia,kind=4)
         reqscount=reqscount+1

! Get ia_21 for GM2_1
         call get_mpi_proc(ng2,nproc,ia_21,recvrank)
         call MPI_IRECV(XGM2_1(ia_12-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,recvrank,tag_1,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

! ia_21 index is on this MPI process
       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then

        if (.not.(locblock)) then
! Get rank of MPI process with ia_12 integrals
         call get_mpi_proc(ng2,nproc,ia_12,sendrank)

! MPI send ia_21 integrals for GM2_1
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_1=int(ia,kind=4)
         reqscount=reqscount+1
         call MPI_ISEND(GM2_1(ia_21-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,sendrank,tag_1,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)

        end if

       end if

            end do
            end do
         end do
         end do

         if (.not.locblock) then
C          call MPI_WAITALL(nmsgs,reqs,MPI_STATUSES_IGNORE,ierr)
C          if (ierr.ne.0) write(*,*) "Trouble with GM2_1 waitall"
C          do i=1,nmsgs
C            req=reqs(i)
C            call MPI_WAIT(req,stat,ierr)
C          end do
          call RXCHF_MPI_WAITALL2(rank,nmsgs,reqs)
         end if

      end do
      end do

      if(allocated(reqs)) deallocate(reqs)

! Symmetrize integrals locally
      do i=1,ng2loc

        x12=GM2_1(i)
        x21=XGM2_1(i)
        GM2_1(i)=(x12+x21)/two

      end do

! Done with GM2_1 integrals
      if(allocated(XGM2_1)) deallocate(XGM2_1)

! Symmetrize GM2_3 integrals

! Allocate storage for temporary arrays to store ia_21 integrals
      if(allocated(XGM2_3)) deallocate(XGM2_3)
      allocate(XGM2_3(ng2loc))
      XGM2_3=0.0d+00

! Get appropriate ia_ji values for XGM2_3
      if (rank.eq.0) write(*,*) "Symmetrizing GM2_3"

      do ip=1,npbf
      do jp=1,npbf

! Get indices of first/last integrals to be passed for fixed ip,jp
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,1,1,1,1,ia_first)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,nebf,nebf,nebf,nebf,ia_last)
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

C         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C         if(rank.eq.0) write(*,1001) "index:",rank,ip,jp
C         write(*,*) "nmsgs:",rank,nmsgs
C         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf

         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia_12)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,ia_21)

! ia_12 index is on this MPI process
       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then

        if (locblock) then
         XGM2_3(ia_12-arrstart+1)=GM2_3(ia_21-arrstart+1)
        else
! MPI receive ia_21 integrals and store in X* arrs at ia_12 index
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_3=int(ia,kind=4)
         reqscount=reqscount+1

! Get ia_21 for GM2_3
         call get_mpi_proc(ng2,nproc,ia_21,recvrank)
         call MPI_IRECV(XGM2_3(ia_12-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,recvrank,tag_3,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

! ia_21 index is on this MPI process
       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then

        if (.not.(locblock)) then
! Get rank of MPI process with ia_12 integrals
         call get_mpi_proc(ng2,nproc,ia_12,sendrank)

! MPI send ia_21 integrals for GM2_3
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_3=int(ia,kind=4)
         reqscount=reqscount+1
         call MPI_ISEND(GM2_3(ia_21-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,sendrank,tag_3,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)

        end if

       end if

            end do
            end do
         end do
         end do

         if (.not.locblock) then
C          call MPI_WAITALL(nmsgs,reqs,MPI_STATUSES_IGNORE,ierr)
C          if (ierr.ne.0) write(*,*) "Trouble with GM2_3 waitall"
C          do i=1,nmsgs
C            req=reqs(i)
C            call MPI_WAIT(req,stat,ierr)
C          end do
          call RXCHF_MPI_WAITALL2(rank,nmsgs,reqs)
         end if

      end do
      end do

      if(allocated(reqs)) deallocate(reqs)

! Symmetrize integrals locally
      do i=1,ng2loc

        x12=GM2_3(i)
        x21=XGM2_3(i)
        GM2_3(i)=(x12+x21)/two

      end do

! Done with GM2_3 integrals
      if(allocated(XGM2_3)) deallocate(XGM2_3)

! Symmetrize GM2s integrals

! Allocate storage for temporary arrays to store ia_21 integrals
      if(allocated(XGM2s)) deallocate(XGM2s)
      allocate(XGM2s(ng2loc))
      XGM2s=0.0d+00

! Get appropriate ia_ji values for XGM2s
      if (rank.eq.0) write(*,*) "Symmetrizing GM2s"

      do ip=1,npbf
      do jp=1,npbf

! Get indices of first/last integrals to be passed for fixed ip,jp
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,1,1,1,1,ia_first)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,nebf,nebf,nebf,nebf,ia_last)
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

C         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C         if(rank.eq.0) write(*,1001) "index:",rank,ip,jp
C         write(*,*) "nmsgs:",rank,nmsgs
C         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf

         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia_12)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,ia_21)

! ia_12 index is on this MPI process
       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then

        if (locblock) then
         XGM2s(ia_12-arrstart+1)=GM2s(ia_21-arrstart+1)
        else
! MPI receive ia_21 integrals and store in X* arrs at ia_12 index
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_s=int(ia,kind=4)
         reqscount=reqscount+1

! Get ia_21 for GM2s
         call get_mpi_proc(ng2,nproc,ia_21,recvrank)
         call MPI_IRECV(XGM2s(ia_12-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,recvrank,tag_s,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

! ia_21 index is on this MPI process
       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then

        if (.not.(locblock)) then
! Get rank of MPI process with ia_12 integrals
         call get_mpi_proc(ng2,nproc,ia_12,sendrank)

! MPI send ia_21 integrals for GM2s
         call pack_4D(nebf,nebf,nebf,
     x                jec2,iec2,jec1,iec1,ia)
         tag_s=int(ia,kind=4)
         reqscount=reqscount+1
         call MPI_ISEND(GM2s(ia_21-arrstart+1),1,
     x                  MPI_DOUBLE_PRECISION,sendrank,tag_s,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)

        end if

       end if

            end do
            end do
         end do
         end do

         if (.not.locblock) then
C          call MPI_WAITALL(nmsgs,reqs,MPI_STATUSES_IGNORE,ierr)
C          if (ierr.ne.0) write(*,*) "Trouble with GM2s waitall"
C          do i=1,nmsgs
C            req=reqs(i)
C            call MPI_WAIT(req,stat,ierr)
C          end do
          call RXCHF_MPI_WAITALL2(rank,nmsgs,reqs)
         end if

      end do
      end do

      if(allocated(reqs)) deallocate(reqs)

! Symmetrize integrals locally
      do i=1,ng2loc

        x12=GM2s(i)
        x21=XGM2s(i)
        GM2s(i)=(x12+x21)/two

      end do

! Done with GM2s integrals
      if(allocated(XGM2s)) deallocate(XGM2s)

      wtime2 = MPI_WTIME() - wtime

      if(rank.eq.0) write(*,3000) 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      write(*,3001) rank,wtime2
!--------------------SYMMETRIZE----------------------------------------)

! Construct global array on master process for testing
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if(allocated(TGM2_1)) deallocate(TGM2_1)
!      if(allocated(TGM2_2)) deallocate(TGM2_2)
!      if(allocated(TGM2_3)) deallocate(TGM2_3)
!      if(allocated(TGM2s)) deallocate(TGM2s)
!      if (rank.eq.0) then
!       allocate(TGM2_1(ng2))
!       allocate(TGM2_2(ng2))
!       allocate(TGM2_3(ng2))
!       allocate(TGM2s(ng2))
!      else
!       allocate(TGM2_1(1))
!       allocate(TGM2_2(1))
!       allocate(TGM2_3(1))
!       allocate(TGM2s(1))
!      end if
!      TGM2_1=zero
!      TGM2_2=zero
!      TGM2_3=zero
!      TGM2s=zero
!
!      ng2loc4=int(ng2loc,kind=4)
!
!! Get number of elements calculated by each proc
!      call MPI_GATHER(ng2loc4,1,MPI_INTEGER,
!     x                ng2locarr(1),1,MPI_INTEGER,
!     x                0,MPI_COMM_WORLD,ierr)
!
!! Get displacements for array storage
!      if (rank.eq.0) then
!        displarr(1)=0
!        do i=2,nproc
!          displarr(i)=displarr(i-1)+ng2locarr(i-1)
!        end do
!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM2_1 on root
!      call MPI_GATHERV(GM2_1(1),ng2loc,MPI_DOUBLE_PRECISION,
!     x                 TGM2_1(1),ng2locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM2_2 on root
!      call MPI_GATHERV(GM2_2(1),ng2loc,MPI_DOUBLE_PRECISION,
!     x                 TGM2_2(1),ng2locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM2_3 on root
!      call MPI_GATHERV(GM2_3(1),ng2loc,MPI_DOUBLE_PRECISION,
!     x                 TGM2_3(1),ng2locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!! Form global GM2s on root
!      call MPI_GATHERV(GM2s(1),ng2loc,MPI_DOUBLE_PRECISION,
!     x                 TGM2s(1),ng2locarr,displarr,
!     x                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!!      if (rank.eq.0) then
!!       write(*,*) "concatenated ng2"
!!       do i=1,ng2
!!        write(*,9001) TGM2_1(i),TGM2_2(i),
!!     x                TGM2_3(i),TGM2s(i)
!!       end do
!!      end if
!
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      if (rank.eq.0) then
!       open(unit=20,file="XCHF_GAM2.ufm",form="unformatted")
!       write(20) TGM2_1
!       close(20)
!       write(*,*) "XCHF_GAM2 written to disk"
!       open(unit=21,file="INT_GAM2.ufm",form="unformatted")
!       write(21) TGM2_2
!       close(21)
!       write(*,*) "INT_GAM2 written to disk"
!       open(unit=22,file="INT_GAM2ex.ufm",form="unformatted")
!       write(22) TGM2_3
!       close(22)
!       write(*,*) "INT_GAM2ex written to disk"
!       open(unit=23,file="XCHF_GAM2s.ufm",form="unformatted")
!       write(23) TGM2s
!       close(23)
!       write(*,*) "XCHF_GAM2s written to disk"
!      end if
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!      if(allocated(TGM2s)) deallocate(TGM2s)
!      if(allocated(TGM2_3)) deallocate(TGM2_3)
!      if(allocated(TGM2_2)) deallocate(TGM2_2)
!      if(allocated(TGM2_1)) deallocate(TGM2_1)

 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 3-PARTICLE INTEGRALS        |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 3-PARTICLE INTEGRALS: ',1X,I12/
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
!=======================================================================
      subroutine RXCHF_GAM2_thread_MPI(istart,iend,ng2_seg,ng2,
     x                          nebf,npebf,npbf,nat,ngtg1,
     x                          pmass,cat,zan,bcoef1,gamma1,
     x                          loop_map,arrstart,
     x                          GM2_1,GM2_2,GM2s,
     x                          KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                          ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

!=======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng2_seg
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas
      integer ng2      ! Number of integrals calc by MPI process
      integer arrstart ! Index of first integral

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
      integer loop_map(ng2_seg,6)

! Variables Returned
      double precision GM2_1(ng2),GM2_2(ng2)
      double precision GM2s(ng2)

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap,ia
      double precision OMG2_1,OMG2_2,OMG2s

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(istart,iend)
!$ompx shared(loop_map,arrstart)
!$ompx shared(GM2_1,GM2_2,GM2s)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(OMG2_1,OMG2_2,OMG2s)
!$ompx private(ia)
!$ompx private(id)

!     id= omp_get_thread_num()
!     write(*,*)' Hello from process ',id
!     if(id.eq.0) then
!        write(*,*)'Threads in use', omp_get_num_threads()
!     end if

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec2=loop_map(imap,1)
         iec2=loop_map(imap,2)
         jec1=loop_map(imap,3)
         iec1=loop_map(imap,4)
         jp =loop_map(imap,5)
         ip =loop_map(imap,6)

         call RXCHFmult_contract_omega2_conv(ip,jp,iec1,jec1,iec2,jec2,
     x                            nebf,npebf,npbf,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            OMG2_1,OMG2_2,OMG2s)


         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia)

         GM2_1(ia-arrstart+1)=OMG2_1
         GM2_2(ia-arrstart+1)=OMG2_2
         GM2s(ia-arrstart+1)=OMG2s

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

      return
      end
!=======================================================================
      subroutine RXCHF_GAM2ex_thread_MPI(istart,iend,ng2_seg,ng2,
     x                            nebf,npebf,npbf,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            loop_map,arrstart,
     x                            GM2_1,GM2_2,GM2_3,GM2s,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

!=======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng2_seg
      integer npebf    ! Number primitive electronic basis functions
      integer nebf     ! Number contracted electronic basis functions
      integer npbf     ! Number nuclear basis functions
      integer nat      ! Number of atoms
      integer ngtg1    ! Number BGammas
      integer ng2      ! Number of integrals calc by MPI process
      integer arrstart ! Index of first integral

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
      integer loop_map(ng2_seg,6)

! Variables Returned
      double precision GM2_1(ng2),GM2_2(ng2)
      double precision GM2_3(ng2),GM2s(ng2)

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap,ia
      double precision OMG2_1,OMG2_2,OMG2_3,OMG2s

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(istart,iend)
!$ompx shared(loop_map,arrstart)
!$ompx shared(GM2_1,GM2_2,GM2_3,GM2s)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(OMG2_1,OMG2_2,OMG2_3,OMG2s)
!$ompx private(ia)
!$ompx private(id)

!     id= omp_get_thread_num()
!     write(*,*)' Hello from process ',id
!     if(id.eq.0) then
!        write(*,*)'Threads in use', omp_get_num_threads()
!     end if

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec2=loop_map(imap,1)
         iec2=loop_map(imap,2)
         jec1=loop_map(imap,3)
         iec1=loop_map(imap,4)
         jp =loop_map(imap,5)
         ip =loop_map(imap,6)

        call RXCHFmult_contract_omega2_convex(ip,jp,iec1,jec1,iec2,jec2,
     x                            nebf,npebf,npbf,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            OMG2_1,OMG2_2,OMG2_3,OMG2s)


         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia)

         GM2_1(ia-arrstart+1)=OMG2_1
         GM2_2(ia-arrstart+1)=OMG2_2
         GM2_3(ia-arrstart+1)=OMG2_3
         GM2s(ia-arrstart+1)=OMG2s

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

      return
      end

