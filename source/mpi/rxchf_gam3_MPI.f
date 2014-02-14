!=======================================================================
      subroutine RXCHF_GAM3_MPI(nproc,rank,
     x                          nblocks,blockrank,
     x                          blockstart,blockend,
     x                          Nchunks,nebf,npebf,npbf,
     x                          ng3,ng3loc,ng3prm,nat,ngtg1,
     x                          pmass,cat,zan,bcoef1,gamma1,
     x                          KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                          ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                          LSYMM,GM3_1,GM3_2)

!=======================================================================
      implicit none
      include 'mpif.h'

! Input Variables
      integer nblocks
      integer blockrank
      integer blockstart,blockend
      integer Nchunks
      integer ng3             ! Total number of integrals
      integer ng3loc          ! Number of integrals for MPI proc to calc
      integer nebf,npebf,npbf,ng3prm
      integer nat,ngtg1
      logical LSYMM           ! Flag to symmetrize GM3_1
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
      integer ierr
      integer mpistart,mpiend,arrstart

      character*4  istring              ! File I/O variables

      double precision zero,half,six
      parameter(zero=0.0d+00,half=0.5d+00,six=6.0d+00)

      double precision wtime
      double precision wtime1
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

      mpistart=mpistart+blockstart-1
      mpiend=mpiend+blockstart-1
C      write(*,*) "nblocks,blockrank,blockstart,blockend,
C     x            mpistart,mpiend:",rank,nblocks,blockrank,
C     x            blockstart,blockend,mpistart,mpiend

      if (rank.eq.0) then
       write(*,1000) ng3,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      GM3_1=0.0d+00
      GM3_2=0.0d+00

C Variables for file I/O
      write(istring,'(I4.4)') rank

!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2000) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      wtime = MPI_WTIME()

      do ichunk=1,Nchunks

         wtime1 = MPI_WTIME()

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

         call RXCHF_writeint_mpi(nproc,rank,ng3loc,18,
     x                           "XCHF_GAM3-"//istring//".ufm",
     x                           GM3_1)
         call RXCHF_writeint_mpi(nproc,rank,ng3loc,17,
     x                           "INT_GAM3-"//istring//".ufm",
     x                           GM3_2)

         wtime2 = MPI_WTIME() - wtime1
         write(*,2001) rank,ichunk,wtime2

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------)

      if(allocated(loop_map)) deallocate(loop_map)

      wtime2 = MPI_WTIME() - wtime

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2010) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      write(*,2011) rank,wtime2
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng3loc:",mpistart,mpiend,ng3loc
!      do i=1,ng3loc
!       write(*,9001) GM3_1(i),GM3_2(i)
!      end do

! Symmetrize only if all integrals have been calculated
      if (nblocks.eq.1) then

!--------------------SYMMETRIZE----------------------------------------(
C Symmetrized integrals in GM3_1ICR (XCHF integrals) (if LSYMM)
C Symmetrized integrals in GM3_2ICR (interaction integrals)

      wtime = MPI_WTIME() 

      if (rank.eq.0) write(*,*) "Symmetrizing GM3_2"
      call gam3_symm1(nproc,rank,
     x                nebf,npbf,ng3,ng3loc,GM3_2)

      if (LSYMM) then
C       if (rank.eq.0) write(*,*) "Symmetrizing GM3_1"
C       call RXCHF_GAM3_MPI_symm(nproc,rank,
C     x                          ng3,ng3loc,
C     x                          mpistart,mpiend,arrstart,
C     x                          nebf,npbf,
C     x                          GM3_1)
       if (rank.eq.0) write(*,*) "GM3_1 symm not supported yet"
      end if
      if (rank.eq.0) write(*,*) "Not symmetrizing GM3_1"

      wtime2 = MPI_WTIME() - wtime

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(rank.eq.0) write(*,3000) 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      write(*,3001) rank,wtime2
!--------------------SYMMETRIZE----------------------------------------)

      else
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       if(rank.eq.0) then
        write(*,*) "  NOT SYMMETRIZING GAM3!!!  "
       end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      end if

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

 2000 FORMAT(/8X,'  INTEGRAL BLOCK CALCULATION TIMINGS:',/,
     x        8X,' -------------------------------------')

 2001 FORMAT( 8X,' PROCESS ',1X,I4,1X,' BLOCK ',1X,I4,1X,F10.2)

 2010 FORMAT(/8X,'  TOTAL INTEGRAL CALCULATION TIMINGS:',/,
     x        8X,'  -----------------------------------')

 2011 FORMAT( 8X,'       PROCESS ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,' INTEGRAL SYMMETRIZATION TIMINGS:',/,
     x        8X,' --------------------------------')

 3001 FORMAT( 8X,'     PROCESS ',1X,I4,1X,F10.2)

 9001 FORMAT(1X,2(F20.10))

      return
      end
!=======================================================================
      subroutine RXCHF_GAM3ex_MPI(nproc,rank,
     x                            nblocks,blockrank,
     x                            blockstart,blockend,
     x                            Nchunks,nebf,npebf,npbf,
     x                            ng3,ng3loc,ng3prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            LSYMM,GM3_1,GM3_2,
     x                            GM3_3,GM3_4)

!=======================================================================
      implicit none
      include 'mpif.h'

! Input Variables
      integer nblocks
      integer blockrank
      integer blockstart,blockend
      integer Nchunks
      integer ng3             ! Total number of integrals
      integer ng3loc          ! Number of integrals for MPI proc to calc
      integer nebf,npebf,npbf,ng3prm
      integer nat,ngtg1
      logical LSYMM           ! Flag to symmetrize GM3_1
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
      integer ierr
      integer mpistart,mpiend,arrstart

      character*4  istring              ! File I/O variables

      double precision zero,half,six
      parameter(zero=0.0d+00,half=0.5d+00,six=6.0d+00)

      double precision wtime
      double precision wtime1
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

      mpistart=mpistart+blockstart-1
      mpiend=mpiend+blockstart-1
C      write(*,*) "nblocks,blockrank,blockstart,blockend,
C     x            mpistart,mpiend:",rank,nblocks,blockrank,
C     x            blockstart,blockend,mpistart,mpiend

      if (rank.eq.0) then
       write(*,1000) ng3,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      GM3_1=0.0d+00
      GM3_2=0.0d+00
      GM3_3=0.0d+00
      GM3_4=0.0d+00

C Variables for file I/O
      write(istring,'(I4.4)') rank

!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2000) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      wtime = MPI_WTIME()

      do ichunk=1,Nchunks

         wtime1 = MPI_WTIME()

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

         call RXCHF_writeint_mpi(nproc,rank,ng3loc,18,
     x                           "XCHF_GAM3-"//istring//".ufm",
     x                           GM3_1)
         call RXCHF_writeint_mpi(nproc,rank,ng3loc,17,
     x                           "INT_GAM3-"//istring//".ufm",
     x                           GM3_2)
         call RXCHF_writeint_mpi(nproc,rank,ng3loc,20,
     x                           "INT_GAM3ex1-"//istring//".ufm",
     x                           GM3_3)
         call RXCHF_writeint_mpi(nproc,rank,ng3loc,20,
     x                           "INT_GAM3ex2-"//istring//".ufm",
     x                           GM3_4)

         wtime2 = MPI_WTIME() - wtime1
         write(*,2001) rank,ichunk,wtime2

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------)

      if(allocated(loop_map)) deallocate(loop_map)

      wtime2 = MPI_WTIME() - wtime

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(rank.eq.0) write(*,2010) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      write(*,2011) rank,wtime2
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "start,end,ng3loc:",mpistart,mpiend,ng3loc
!      do i=1,ng3loc
!       write(*,9001) GM3_1(i),GM3_2(i),
!     x               GM3_3(i),GM3_4(i)
!      end do

! Symmetrize only if all integrals have been calculated
      if (nblocks.eq.1) then

!--------------------SYMMETRIZE----------------------------------------(
C Symmetrized integrals in GM3_1ICR (XCHF integrals) (if LSYMM)
C Symmetrized integrals in GM3_2ICR (interaction integrals)
C Symmetrized integrals in GM3_3ICR (exchange integrals)
C Symmetrized integrals in GM3_4ICR (exchange integrals)

      wtime = MPI_WTIME()

      if (rank.eq.0) write(*,*) "Symmetrizing GM3_2"
      call gam3_symm1(nproc,rank,
     x                nebf,npbf,ng3,ng3loc,GM3_2)

      if (rank.eq.0) write(*,*) "Symmetrizing GM3_3"
      call gam3_symm2(nproc,rank,
     x                nebf,npbf,ng3,ng3loc,GM3_3)

      if (rank.eq.0) write(*,*) "Symmetrizing GM3_4"
      call gam3_symm3(nproc,rank,
     x                nebf,npbf,ng3,ng3loc,GM3_4)

      if (LSYMM) then
C       if (rank.eq.0) write(*,*) "Symmetrizing GM3_1"
C       call RXCHF_GAM3_MPI_symm(nproc,rank,
C     x                          ng3,ng3loc,
C     x                          mpistart,mpiend,arrstart,
C     x                          nebf,npbf,
C     x                          GM3_1)
       if (rank.eq.0) write(*,*) "GM3_1 symm not supported yet"
      end if
      if (rank.eq.0) write(*,*) "Not symmetrizing GM3_1"

      wtime2 = MPI_WTIME() - wtime

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(rank.eq.0) write(*,3000) 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      write(*,3001) rank,wtime2

!--------------------SYMMETRIZE----------------------------------------)

      else
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       if(rank.eq.0) then
        write(*,*) "  NOT SYMMETRIZING GAM3!!!  "
       end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      end if

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
!      if (rank.eq.0) then
!       write(*,*) "concatenated ng3"
!       do i=1,ng3
!        write(*,9001) TGM3_1(i),TGM3_2(i),
!     x                TGM3_3(i),TGM3_4(i)
!       end do
!      end if
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

 2000 FORMAT(/8X,'  INTEGRAL BLOCK CALCULATION TIMINGS:',/,
     x        8X,' -------------------------------------')

 2001 FORMAT( 8X,' PROCESS ',1X,I4,1X,' BLOCK ',1X,I4,1X,F10.2)

 2010 FORMAT(/8X,'  TOTAL INTEGRAL CALCULATION TIMINGS:',/,
     x        8X,'  -----------------------------------')

 2011 FORMAT( 8X,'       PROCESS ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,' INTEGRAL SYMMETRIZATION TIMINGS:',/,
     x        8X,' --------------------------------')

 3001 FORMAT( 8X,'     PROCESS ',1X,I4,1X,F10.2)

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
      subroutine gam3_symm1(nproc,rank,
     x                      nebf,npbf,ntotints,nlocints,GAM3)

C=======================================================================
      implicit none
      include "mpif.h"

! Input variables
      integer nproc,rank
      integer nebf,npbf
      integer ntotints,nlocints

! Input/output variables
      double precision GAM3(nlocints)

! Local variables
      integer i
      integer mpistart,mpiend
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3
      integer ia_first,ia_last
      integer ia,ia_12,ia_21
      integer tag,reqscount,ierr
      integer sendrank,recvrank
      integer reqs(2*nlocints)
      double precision x12,x21
      double precision xGAM3(nlocints)

      xGAM3=0.0d+00

      call get_mpi_range(ntotints,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ntotints

! Get appropriate ia_ji values for xGAM3

      do ip=1,npbf
      do jp=1,npbf
         do ie1=1,nebf
         do je1=1,nebf

            reqscount=0

            do ie2=1,nebf
            do je2=1,nebf
               do ie3=1,nebf
               do je3=1,nebf

             call index_GAM_3PK(nebf,npbf,
     x                          ip,jp,ie1,je1,ie2,je2,ie3,je3,ia_12)
             call index_GAM_3PK(nebf,npbf,
     x                          ip,jp,ie1,je1,ie3,je3,ie2,je2,ia_21)

       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then

         call get_mpi_proc(ntotints,nproc,ia_21,recvrank)

         if (recvrank.eq.rank) then
          xGAM3(ia_12-mpistart+1)=GAM3(ia_21-mpistart+1)
         else
          call pack_4D(nebf,nebf,nebf,
     x                 je3,ie3,je2,ie2,ia)
          tag=ia
          reqscount=reqscount+1
          call MPI_IRECV(xGAM3(ia_12-mpistart+1),1,
     x                   MPI_DOUBLE_PRECISION,recvrank,tag,
     x                   MPI_COMM_WORLD,reqs(reqscount),ierr)
         end if

       end if

       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then

         call get_mpi_proc(ntotints,nproc,ia_12,sendrank)

         if (sendrank.ne.rank) then
          call pack_4D(nebf,nebf,nebf,
     x                 je3,ie3,je2,ie2,ia)
          tag=ia
          reqscount=reqscount+1
          call MPI_ISEND(GAM3(ia_21-mpistart+1),1,
     x                   MPI_DOUBLE_PRECISION,sendrank,tag,
     x                   MPI_COMM_WORLD,reqs(reqscount),ierr)
         end if

       end if

               end do
               end do
            end do
            end do

            if(reqscount.eq.0) then
             reqscount=1
             reqs(1)=MPI_REQUEST_NULL
            end if
            call MPI_WAITALL(reqscount,reqs(1:reqscount),
     x                       MPI_STATUSES_IGNORE,ierr)
            if (ierr.ne.0) write(*,*) "Trouble with GAM3 waitall"
C            if(reqscount.ne.1) then
C             write(*,*) "Success with waitall:",
C     x                  rank,ip,jp,ie1,je1,reqscount
C            end if

         end do
         end do
      end do
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Symmetrize integrals locally
      do i=1,nlocints
        x12=GAM3(i)
        x21=xGAM3(i)
        GAM3(i)=(x12+x21)/2.0d+00
      end do

      return
      end
C=======================================================================
      subroutine gam3_symm2(nproc,rank,
     x                      nebf,npbf,ntotints,nlocints,GAM3)

C=======================================================================
      implicit none
      include "mpif.h"

! Input variables
      integer nproc,rank
      integer nebf,npbf
      integer ntotints,nlocints

! Input/output variables
      double precision GAM3(nlocints)

! Local variables
      integer i
      integer mpistart,mpiend
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3
      integer ia_first,ia_last
      integer ia,ia_12,ia_21
      integer tag,reqscount,ierr
      integer sendrank,recvrank
      integer reqs(2*nlocints)
      double precision x12,x21
      double precision xGAM3(nlocints)

      xGAM3=0.0d+00

      call get_mpi_range(ntotints,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ntotints

! Get appropriate ia_ji values for xGAM3

      do ip=1,npbf
      do jp=1,npbf
         do ie3=1,nebf
         do je3=1,nebf

            reqscount=0

            do ie1=1,nebf
            do je1=1,nebf
               do ie2=1,nebf
               do je2=1,nebf

             call index_GAM_3PK(nebf,npbf,
     x                          ip,jp,ie1,je1,ie2,je2,ie3,je3,ia_12)
             call index_GAM_3PK(nebf,npbf,
     x                          ip,jp,ie2,je2,ie1,je1,ie3,je3,ia_21)

       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then

         call get_mpi_proc(ntotints,nproc,ia_21,recvrank)

         if (recvrank.eq.rank) then
          xGAM3(ia_12-mpistart+1)=GAM3(ia_21-mpistart+1)
         else
          call pack_4D(nebf,nebf,nebf,
     x                 je2,ie2,je1,ie1,ia)
          tag=ia
          reqscount=reqscount+1
          call MPI_IRECV(xGAM3(ia_12-mpistart+1),1,
     x                   MPI_DOUBLE_PRECISION,recvrank,tag,
     x                   MPI_COMM_WORLD,reqs(reqscount),ierr)
         end if

       end if

       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then

         call get_mpi_proc(ntotints,nproc,ia_12,sendrank)

         if (sendrank.ne.rank) then
          call pack_4D(nebf,nebf,nebf,
     x                 je2,ie2,je1,ie1,ia)
          tag=ia
          reqscount=reqscount+1
          call MPI_ISEND(GAM3(ia_21-mpistart+1),1,
     x                   MPI_DOUBLE_PRECISION,sendrank,tag,
     x                   MPI_COMM_WORLD,reqs(reqscount),ierr)
         end if

       end if

               end do
               end do
            end do
            end do

            if(reqscount.eq.0) then
             reqscount=1
             reqs(1)=MPI_REQUEST_NULL
            end if
            call MPI_WAITALL(reqscount,reqs(1:reqscount),
     x                       MPI_STATUSES_IGNORE,ierr)
            if (ierr.ne.0) write(*,*) "Trouble with GAM3 waitall"
C            if(reqscount.ne.1) then
C             write(*,*) "Success with waitall:",
C     x                  rank,ip,jp,ie3,je3,reqscount
C            end if

         end do
         end do
      end do
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Symmetrize integrals locally
      do i=1,nlocints
        x12=GAM3(i)
        x21=xGAM3(i)
        GAM3(i)=(x12+x21)/2.0d+00
      end do

      return
      end
C=======================================================================
      subroutine gam3_symm3(nproc,rank,
     x                      nebf,npbf,ntotints,nlocints,GAM3)

C=======================================================================
      implicit none
      include "mpif.h"

! Input variables
      integer nproc,rank
      integer nebf,npbf
      integer ntotints,nlocints

! Input/output variables
      double precision GAM3(nlocints)

! Local variables
      integer i
      integer mpistart,mpiend
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3
      integer ia_first,ia_last
      integer ia,ia_12,ia_21
      integer tag,reqscount,ierr
      integer sendrank,recvrank
      integer reqs(2*nlocints)
      double precision x12,x21
      double precision xGAM3(nlocints)

      xGAM3=0.0d+00

      call get_mpi_range(ntotints,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ntotints

! Get appropriate ia_ji values for xGAM3

      do ip=1,npbf
      do jp=1,npbf
         do ie2=1,nebf
         do je2=1,nebf

            reqscount=0

            do ie1=1,nebf
            do je1=1,nebf
               do ie3=1,nebf
               do je3=1,nebf

             call index_GAM_3PK(nebf,npbf,
     x                          ip,jp,ie1,je1,ie2,je2,ie3,je3,ia_12)
             call index_GAM_3PK(nebf,npbf,
     x                          ip,jp,ie3,je3,ie2,je2,ie1,je1,ia_21)

       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then

         call get_mpi_proc(ntotints,nproc,ia_21,recvrank)

         if (recvrank.eq.rank) then
          xGAM3(ia_12-mpistart+1)=GAM3(ia_21-mpistart+1)
         else
          call pack_4D(nebf,nebf,nebf,
     x                 je3,ie3,je1,ie1,ia)
          tag=ia
          reqscount=reqscount+1
          call MPI_IRECV(xGAM3(ia_12-mpistart+1),1,
     x                   MPI_DOUBLE_PRECISION,recvrank,tag,
     x                   MPI_COMM_WORLD,reqs(reqscount),ierr)
         end if

       end if

       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then

         call get_mpi_proc(ntotints,nproc,ia_12,sendrank)

         if (sendrank.ne.rank) then
          call pack_4D(nebf,nebf,nebf,
     x                 je3,ie3,je1,ie1,ia)
          tag=ia
          reqscount=reqscount+1
          call MPI_ISEND(GAM3(ia_21-mpistart+1),1,
     x                   MPI_DOUBLE_PRECISION,sendrank,tag,
     x                   MPI_COMM_WORLD,reqs(reqscount),ierr)
         end if

       end if

               end do
               end do
            end do
            end do

            if(reqscount.eq.0) then
             reqscount=1
             reqs(1)=MPI_REQUEST_NULL
            end if
            call MPI_WAITALL(reqscount,reqs(1:reqscount),
     x                       MPI_STATUSES_IGNORE,ierr)
            if (ierr.ne.0) write(*,*) "Trouble with GAM3 waitall"
C            if(reqscount.ne.1) then
C             write(*,*) "Success with waitall:",
C     x                  rank,ip,jp,ie2,je2,reqscount
C            end if

         end do
         end do
      end do
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Symmetrize integrals locally
      do i=1,nlocints
        x12=GAM3(i)
        x21=xGAM3(i)
        GAM3(i)=(x12+x21)/2.0d+00
      end do

      return
      end

