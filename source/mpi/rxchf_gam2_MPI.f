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
      integer ierr
      integer mpistart,mpiend,arrstart

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

      if(allocated(loop_map)) deallocate(loop_map)

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

      if (rank.eq.0) write(*,*) "Symmetrizing GM2_1"
      call gam2_symm(nproc,rank,nebf,npbf,ng2,ng2loc,GM2_1)

      if (rank.eq.0) write(*,*) "Symmetrizing GM2s"
      call gam2_symm(nproc,rank,nebf,npbf,ng2,ng2loc,GM2s)

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
      integer ierr
      integer mpistart,mpiend,arrstart

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

      if(allocated(loop_map)) deallocate(loop_map)

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

      if (rank.eq.0) write(*,*) "Symmetrizing GM2_1"
      call gam2_symm(nproc,rank,nebf,npbf,ng2,ng2loc,GM2_1)

      if (rank.eq.0) write(*,*) "Symmetrizing GM2_3"
      call gam2_symm(nproc,rank,nebf,npbf,ng2,ng2loc,GM2_3)

      if (rank.eq.0) write(*,*) "Symmetrizing GM2s"
      call gam2_symm(nproc,rank,nebf,npbf,ng2,ng2loc,GM2s)

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
!=======================================================================
      subroutine gam2_symm(nproc,rank,
     x                     nebf,npbf,ntotints,nlocints,GAM2)

!=======================================================================
      implicit none
      include "mpif.h"

! Input variables
      integer nproc,rank
      integer nebf,npbf
      integer ntotints,nlocints

! Input/Output variables
      double precision GAM2(nlocints)

! Local variables
      integer i
      integer mpistart,mpiend
      integer ip,jp,ie1,je1,ie2,je2
      integer ia_first,ia_last
      integer ia,ia_12,ia_21
      integer tag,reqscount,ierr
      integer sendrank,recvrank
      integer reqs(2*nlocints)
      double precision x12,x21
      double precision xGAM2(nlocints)

      xGAM2=0.0d+00

      call get_mpi_range(ntotints,nproc,rank,mpistart,mpiend)
      if(rank.eq.(nproc-1)) mpiend=ntotints

! Get appropriate ia_ji values for xGAM2

      do ip=1,npbf
      do jp=1,npbf

         reqscount=0

         do ie1=1,nebf
         do je1=1,nebf
            do ie2=1,nebf
            do je2=1,nebf

         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,ie1,je1,ie2,je2,ia_12)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,ie2,je2,ie1,je1,ia_21)

       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then

        call get_mpi_proc(ntotints,nproc,ia_21,recvrank)

        if (recvrank.eq.rank) then
         xGAM2(ia_12-mpistart+1)=GAM2(ia_21-mpistart+1)
        else
         call pack_4D(nebf,nebf,nebf,
     x                je2,ie2,je1,ie1,ia)
         tag=ia
         reqscount=reqscount+1
         call MPI_IRECV(xGAM2(ia_12-mpistart+1),1,
     x                  MPI_DOUBLE_PRECISION,recvrank,tag,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
        end if

       end if

       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then

        call get_mpi_proc(ntotints,nproc,ia_12,sendrank)
        if (sendrank.ne.rank) then
         call pack_4D(nebf,nebf,nebf,
     x                je2,ie2,je1,ie1,ia)
         tag=ia
         reqscount=reqscount+1
         call MPI_ISEND(GAM2(ia_21-mpistart+1),1,
     x                  MPI_DOUBLE_PRECISION,sendrank,tag,
     x                  MPI_COMM_WORLD,reqs(reqscount),ierr)
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
     x                    MPI_STATUSES_IGNORE,ierr)
         if (ierr.ne.0) write(*,*) "Trouble with GAM2 waitall"
C         if(reqscount.ne.1) then
C          write(*,*) "Success with waitall:",rank,ip,jp,reqscount
C         end if

      end do
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Symmetrize integrals locally
      do i=1,nlocints
        x12=GAM2(i)
        x21=xGAM2(i)
        GAM2(i)=(x12+x21)/2.0d+00
      end do

      return
      end

