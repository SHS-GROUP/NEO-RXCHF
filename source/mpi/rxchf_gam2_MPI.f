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
      integer ip,jp,iec1,jec1,iec2,jec2,i
      integer,allocatable :: loop_map(:,:)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

      integer tag_1,tag_s
      integer sendrank,recvrank
      integer*4, allocatable :: reqs_1(:)
      integer*4, allocatable :: reqs_s(:)

      double precision, allocatable :: XGM2_1(:)
      double precision, allocatable :: XGM2s(:)

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

! Allocate storage for temporary arrays to store ia_21 integrals
      if(allocated(XGM2_1)) deallocate(XGM2_1)
      allocate(XGM2_1(ng2loc))
      if(allocated(XGM2s)) deallocate(XGM2s)
      allocate(XGM2s(ng2loc))

! Allocate and initialize arrays for MPI requests
      if(allocated(reqs_1)) deallocate(reqs_1)
      allocate(reqs_1(2*ng2loc))
      if(allocated(reqs_s)) deallocate(reqs_s)
      allocate(reqs_s(2*ng2loc))
      reqs_1=MPI_REQUEST_NULL
      reqs_s=MPI_REQUEST_NULL

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf

!  GAM_2 Symmetrization:
!  Determine packing indices for XGAM_2 integral matrices

!              As Packed-->       XGAM_2(je2,ie2,je1,ie1,jp,ip)
!                           XGAM_2(ip,jp,ie1,je1,ie2,je2,) 
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia_12)
!              As Packed-->       XGAM_2(je1,ie1,je2,ie2,jp,ip)
!                           XGAM_2(ip,jp,ie2,je2,ie1,je1) 
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,ia_21)

C RXCHFmult(
C    index 1: regular electron
C    index 2: special electron
C    index 3: proton
C )

       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then
! ia_12 index is on this MPI process

! Get rank of MPI process with ia_21 integrals
        call get_mpi_proc(ng2,nproc,ia_21,recvrank)
! MPI receive ia_21 integrals and store in X* arrs at ia_12 index
        tag_1=ia_12     ! Check if it still works for longints
        tag_s=ia_12+ng2 !
        call MPI_IRECV(XGM2_1(ia_12-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 recvrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_12-arrstart+1),ierr)
        call MPI_IRECV(XGM2s(ia_12-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 recvrank,tag_s,MPI_COMM_WORLD,
     x                 reqs_s(ia_12-arrstart+1),ierr)

       end if

       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then
! ia_21 index is on this MPI process

! Get rank of MPI process with ia_12 integrals
        call get_mpi_proc(ng2,nproc,ia_12,sendrank)
! MPI send ia_21 integrals
        tag_1=ia_12     ! Check if it still works for longints
        tag_s=ia_12+ng2 !
        call MPI_ISEND(GM2_1(ia_21-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 sendrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_21-arrstart+1+ng2loc),ierr)
        call MPI_ISEND(GM2s(ia_21-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 sendrank,tag_s,MPI_COMM_WORLD,
     x                 reqs_s(ia_21-arrstart+1+ng2loc),ierr)

       end if

            end do
            end do
         end do
         end do
      end do
      end do

! Wait for all messages to be sent and received by all processes
      call MPI_WAITALL(2*ng2loc,reqs_1,MPI_STATUSES_IGNORE,ierr)
      if (ierr.ne.0) write(*,*) "Trouble with reqs_1 waitall"
      call MPI_WAITALL(2*ng2loc,reqs_s,MPI_STATUSES_IGNORE,ierr)
      if (ierr.ne.0) write(*,*) "Trouble with reqs_s waitall"

      if(allocated(reqs_s)) deallocate(reqs_s)
      if(allocated(reqs_1)) deallocate(reqs_1)

! Symmetrize integrals locally
      do i=1,ng2loc

        x12=GM2_1(i)
        x21=XGM2_1(i)
        GM2_1(i)=(x12+x21)/two

        x12=GM2s(i)
        x21=XGM2s(i)
        GM2s(i)=(x12+x21)/two

      end do

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng2loc:",mpistart,mpiend,ng2loc
!      do i=1,ng2loc
!       write(*,9001) GM2_1(i),GM2_2(i),
!     x               GM2s(i)
!      end do

      if(allocated(XGM2s)) deallocate(XGM2s)
      if(allocated(XGM2_1)) deallocate(XGM2_1)

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

      integer tag_1,tag_3,tag_s
      integer sendrank,recvrank
      integer*4, allocatable :: reqs_1(:)
      integer*4, allocatable :: reqs_3(:)
      integer*4, allocatable :: reqs_s(:)

      double precision, allocatable :: XGM2_1(:)
      double precision, allocatable :: XGM2_3(:)
      double precision, allocatable :: XGM2s(:)

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

! Allocate storage for temporary arrays to store ia_21 integrals
      if(allocated(XGM2_1)) deallocate(XGM2_1)
      allocate(XGM2_1(ng2loc))
      if(allocated(XGM2_3)) deallocate(XGM2_3)
      allocate(XGM2_3(ng2loc))
      if(allocated(XGM2s)) deallocate(XGM2s)
      allocate(XGM2s(ng2loc))

! Allocate and initialize arrays for MPI requests
      if(allocated(reqs_1)) deallocate(reqs_1)
      allocate(reqs_1(2*ng2loc))
      if(allocated(reqs_3)) deallocate(reqs_3)
      allocate(reqs_3(2*ng2loc))
      if(allocated(reqs_s)) deallocate(reqs_s)
      allocate(reqs_s(2*ng2loc))
      reqs_1=MPI_REQUEST_NULL
      reqs_3=MPI_REQUEST_NULL
      reqs_s=MPI_REQUEST_NULL

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf

!  GAM_2 Symmetrization:
!  Determine packing indices for XGAM_2 integral matrices

!              As Packed-->       XGAM_2(je2,ie2,je1,ie1,jp,ip)
!                           XGAM_2(ip,jp,ie1,je1,ie2,je2,) 
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia_12)
!              As Packed-->       XGAM_2(je1,ie1,je2,ie2,jp,ip)
!                           XGAM_2(ip,jp,ie2,je2,ie1,je1) 
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,ia_21)

C RXCHFmult(
C    index 1: regular electron
C    index 2: special electron
C    index 3: proton
C )

       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then
! ia_12 index is on this MPI process

! Get rank of MPI process with ia_21 integrals
        call get_mpi_proc(ng2,nproc,ia_21,recvrank)
! MPI receive ia_21 integrals and store in X* arrs at ia_12 index
        tag_1=ia_12       !
        tag_3=ia_12+ng2   ! Check if it still works for longints
        tag_s=ia_12+ng2*2 !
        call MPI_IRECV(XGM2_1(ia_12-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 recvrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_12-arrstart+1),ierr)
        call MPI_IRECV(XGM2_3(ia_12-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 recvrank,tag_3,MPI_COMM_WORLD,
     x                 reqs_3(ia_12-arrstart+1),ierr)
        call MPI_IRECV(XGM2s(ia_12-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 recvrank,tag_s,MPI_COMM_WORLD,
     x                 reqs_s(ia_12-arrstart+1),ierr)

       end if

       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then
! ia_21 index is on this MPI process

! Get rank of MPI process with ia_12 integrals
        call get_mpi_proc(ng2,nproc,ia_12,sendrank)
! MPI send ia_21 integrals
        tag_1=ia_12       !
        tag_3=ia_12+ng2   ! Check if it still works for longints
        tag_s=ia_12+ng2*2 !
        call MPI_ISEND(GM2_1(ia_21-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 sendrank,tag_1,MPI_COMM_WORLD,
     x                 reqs_1(ia_21-arrstart+1+ng2loc),ierr)
        call MPI_ISEND(GM2_3(ia_21-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 sendrank,tag_3,MPI_COMM_WORLD,
     x                 reqs_3(ia_21-arrstart+1+ng2loc),ierr)
        call MPI_ISEND(GM2s(ia_21-arrstart+1),1,MPI_DOUBLE_PRECISION,
     x                 sendrank,tag_s,MPI_COMM_WORLD,
     x                 reqs_s(ia_21-arrstart+1+ng2loc),ierr)

       end if

            end do
            end do
         end do
         end do
      end do
      end do

! Wait for all messages to be sent and received by all processes
      call MPI_WAITALL(2*ng2loc,reqs_1,MPI_STATUSES_IGNORE,ierr)
      if (ierr.ne.0) write(*,*) "Trouble with reqs_1 waitall"
      call MPI_WAITALL(2*ng2loc,reqs_3,MPI_STATUSES_IGNORE,ierr)
      if (ierr.ne.0) write(*,*) "Trouble with reqs_3 waitall"
      call MPI_WAITALL(2*ng2loc,reqs_s,MPI_STATUSES_IGNORE,ierr)
      if (ierr.ne.0) write(*,*) "Trouble with reqs_s waitall"

      if(allocated(reqs_s)) deallocate(reqs_s)
      if(allocated(reqs_3)) deallocate(reqs_3)
      if(allocated(reqs_1)) deallocate(reqs_1)

! Symmetrize integrals locally
      do i=1,ng2loc

        x12=GM2_1(i)
        x21=XGM2_1(i)
        GM2_1(i)=(x12+x21)/two

        x12=GM2_3(i)
        x21=XGM2_3(i)
        GM2_3(i)=(x12+x21)/two

        x12=GM2s(i)
        x21=XGM2s(i)
        GM2s(i)=(x12+x21)/two

      end do

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng2loc:",mpistart,mpiend,ng2loc
!      do i=1,ng2loc
!       write(*,9001) GM2_1(i),GM2_2(i),
!     x               GM2_3(i),GM2s(i)
!      end do

      if(allocated(XGM2s)) deallocate(XGM2s)
      if(allocated(XGM2_3)) deallocate(XGM2_3)
      if(allocated(XGM2_1)) deallocate(XGM2_1)

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

