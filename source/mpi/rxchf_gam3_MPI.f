!=======================================================================
      subroutine RXCHF_GAM3_MPI(nproc,rank,
     x                          Nchunks,nebf,npebf,npbf,
     x                          ng3,ng3prm,nat,ngtg1,
     x                          pmass,cat,zan,bcoef1,gamma1,
     x                          KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                          ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                          GM3_1ICR,GM3_2ICR)

!=======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ng3,nebf,npebf,npbf,ng3prm
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
      double precision GM3_1ICR(ng3)  ! XCHF OMG3  integrals (full symm)
      double precision GM3_2ICR(ng3)  ! INT  OMG3  integrals (partial symm)

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)

      double precision,allocatable :: XG3_1ICR(:), XG3_2ICR(:)

      double precision GM3_1MPI(ng3/nproc)  ! XCHF OMG3    integrals (full symm)
      double precision GM3_2MPI(ng3/nproc)  ! INT  OMG3    integrals (partial symm)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

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
      double precision x231
      double precision x312
      double precision x321
      double precision xxxx

      double precision zero,half,six
      parameter(zero=0.0d+00,half=0.5d+00,six=6.0d+00)

      double precision wtime
      double precision wtime2

! Have each process calculate ng3/nproc integrals according to rank then allgather
! Have each process calculate ng3%nproc remaining integrals
!  - saves waiting for master and eliminates one broadcast

      call get_mpi_range(ng3,nproc,rank,mpistart,mpiend)

      if (rank.eq.0) then
       write(*,1000) ng3,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      if(allocated(XG3_1ICR)) deallocate(XG3_1ICR)
      allocate( XG3_1ICR(ng3),stat=istat )
      if(allocated(XG3_2ICR)) deallocate(XG3_2ICR)
      allocate( XG3_2ICR(ng3),stat=istat )

      XG3_1ICR=0.0d+00
      XG3_2ICR=0.0d+00

      GM3_1MPI=0.0d+00
      GM3_2MPI=0.0d+00

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

         call RXCHFmult_thread_gam3_IC1(istart,iend,ng3_seg,ng3,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        loop_map,XG3_1ICR,XG3_2ICR,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------)
!-----CLEAN-UP-MEMORY-------------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

      wtime2 = MPI_WTIME() - wtime
      write(*,2000)rank,wtime2

! Here, a contiguous block of ng3/nproc integrals are stored in X* arrs
! Pass this to a ng3/nproc-dimensional array to prepare for allgather
      call copy_arr(ng3/nproc,XG3_1ICR(arrstart),GM3_1MPI)
      call copy_arr(ng3/nproc,XG3_2ICR(arrstart),GM3_2MPI)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "ng3:",ng3
!      do istart=1,ng3
!       write(*,9001) XG3_1ICR(istart),XG3_2ICR(istart)
!      end do
!      write(*,*) "ng3/nproc:",ng3/nproc
!      do istart=1,ng3/nproc
!       write(*,9001) GM3_1MPI(istart),GM3_2MPI(istart)
!      end do

! Pass and broadcast these major chunks of arrays to all processes
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHER(GM3_1MPI(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   XG3_1ICR(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM3_1 allgather"
      call MPI_ALLGATHER(GM3_2MPI(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   XG3_2ICR(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM3_2 allgather"

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "ng3 second:",ng3
!      do istart=1,ng3
!       write(*,9001) XG3_1ICR(istart),XG3_2ICR(istart)
!      end do

! At this stage, ng3%nproc last elements are missing from arrays
      if (mod(ng3,nproc).ne.0) then

         wtime = MPI_WTIME()

! Have threads chop calculation of ng3%nproc integrals
         call loop_size(nproc*(ng3/nproc)+1,ng3,1,0,
     x                  istart,iend)
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
                        if (Loopi.eq.1) then
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

         call RXCHFmult_thread_gam3_IC1(istart,iend,ng3_seg,ng3,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        loop_map,XG3_1ICR,XG3_2ICR,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

!-----CLEAN-UP-MEMORY-------------------------------------------------(
         if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

         wtime2 = MPI_WTIME() - wtime
         if (rank.eq.0) write(*,3000) wtime2

!         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!         write(*,*) "ng3+resid:",ng3
!         do istart=1,ng3
!         write(*,9001) XG3_1ICR(istart),XG3_2ICR(istart)
!         end do

      end if ! resid exists

!--------------------SYMMETRIZE----------------------------------------(
      wtime = MPI_WTIME() 

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

C  GAM_3 Symmetrization:
C  Determine packing indices for XGAM_3 integral matrices

C              As Packed-->       XGAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
!                            read(905,REC=ia_123) x123
                             ! x123=GAM_3(ia_123)

C              As Packed-->       XGAM_3(je2,ie2,je3,ie3,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec3,jec3,iec2,jec2,ia)
                             ia_132=ia
!                            read(905,REC=ia_132) x132
                             ! x132=GAM_3(ia_132)

C              As Packed-->       XGAM_3(je3,ie3,je1,ie1,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,iec3,jec3,ia)
                             ia_213=ia
!                            read(905,REC=ia_213) x213
                             ! x213=GAM_3(ia_213)

C              As Packed-->       XGAM_3(je1,ie1,je3,ie3,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec3,jec3,iec1,jec1,ia)
                             ia_231=ia
!                            read(905,REC=ia_231) x231
                             ! x231=GAM_3(ia_231)

C              As Packed-->       XGAM_3(je2,ie2,je1,ie1,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec3,jec3,iec1,jec1,iec2,jec2,ia)
                             ia_312=ia
!                            read(905,REC=ia_312) x312
                             ! x312=GAM_3(ia_312)

C              As Packed-->       XGAM_3(je1,ie1,je2,ie2,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec3,jec3,iec2,jec2,iec1,jec1,ia)
                             ia_321=ia
!                            read(905,REC=ia_321) x321
                             ! x321=GAM_3(ia_321)

C RXCHFmult(
C    index 1: regular electron
C    index 2: special electron 1
C    index 3: special electron 2
C    index 4: proton
C )

C Fully symmetrized integrals in GM3_1ICR
                       x123=XG3_1ICR(ia_123)
                       x132=XG3_1ICR(ia_132)
                       x213=XG3_1ICR(ia_213)
                       x231=XG3_1ICR(ia_231)
                       x312=XG3_1ICR(ia_312)
                       x321=XG3_1ICR(ia_321)
                       xxxx=(x123+x132+x213+x231+x312+x321)/six

                       GM3_1ICR(ia_123)=xxxx 
                       GM3_1ICR(ia_132)=xxxx 
                       GM3_1ICR(ia_213)=xxxx 
                       GM3_1ICR(ia_231)=xxxx 
                       GM3_1ICR(ia_312)=xxxx 
                       GM3_1ICR(ia_321)=xxxx 

C Partially symmetrized integrals in GM3_2ICR
                       x123=XG3_2ICR(ia_123)
                       x132=XG3_2ICR(ia_132)
                       xxxx=(x123+x132)*half

                       GM3_2ICR(ia_123)=xxxx 
                       GM3_2ICR(ia_132)=xxxx 


               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      wtime2 = MPI_WTIME() - wtime

      if(allocated(XG3_2ICR)) deallocate(XG3_2ICR)
      if(allocated(XG3_1ICR)) deallocate(XG3_1ICR)

      if (rank.eq.0) write(*,4000) wtime2
!--------------------SYMMETRIZE----------------------------------------)


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

 2000 FORMAT(8X,'    TIME FOR PROCESS ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,'    TIME FOR RESIDUAL ',1X,F10.2)

 4000 FORMAT(8X,'      TIME TO SYMMETRIZE INTEGRALS:',1X,F12.4/)

 9001 FORMAT(1X,2(F20.10))

      return
      end
!=======================================================================
      subroutine RXCHF_GAM3ex_MPI(nproc,rank,
     x                            Nchunks,nebf,npebf,npbf,
     x                            ng3,ng3prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            GM3_1ICR,GM3_2ICR,
     x                            GM3_3ICR,GM3_4ICR)

!=======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ng3,nebf,npebf,npbf,ng3prm
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
      double precision GM3_1ICR(ng3)  ! XCHF OMG3    integrals (full symm)
      double precision GM3_2ICR(ng3)  ! INT  OMG3    integrals (partial symm)
      double precision GM3_3ICR(ng3)  ! INT  OMG3ex1 integrals (partial symm)
      double precision GM3_4ICR(ng3)  ! INT  OMG3ex2 integrals (partial symm)

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)

      double precision,allocatable :: XG3_1ICR(:),XG3_2ICR(:)
      double precision,allocatable :: XG3_3ICR(:),XG3_4ICR(:)

      double precision GM3_1MPI(ng3/nproc)  ! XCHF OMG3    integrals (full symm)
      double precision GM3_2MPI(ng3/nproc)  ! INT  OMG3    integrals (partial symm)
      double precision GM3_3MPI(ng3/nproc)  ! INT  OMG3ex1 integrals (partial symm)
      double precision GM3_4MPI(ng3/nproc)  ! INT  OMG3ex2 integrals (partial symm)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

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
      double precision x231
      double precision x312
      double precision x321
      double precision xxxx

      double precision zero,half,six
      parameter(zero=0.0d+00,half=0.5d+00,six=6.0d+00)

      double precision wtime
      double precision wtime2

! Have each process calculate ng3/nproc integrals according to rank then allgather
! Have each process calculate ng3%nproc remaining integrals
!  - saves waiting for master and eliminates one broadcast

      call get_mpi_range(ng3,nproc,rank,mpistart,mpiend)

      if (rank.eq.0) then
       write(*,1000) ng3,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      if(allocated(XG3_1ICR)) deallocate(XG3_1ICR)
      allocate( XG3_1ICR(ng3),stat=istat )
      if(allocated(XG3_2ICR)) deallocate(XG3_2ICR)
      allocate( XG3_2ICR(ng3),stat=istat )
      if(allocated(XG3_3ICR)) deallocate(XG3_3ICR)
      allocate( XG3_3ICR(ng3),stat=istat )
      if(allocated(XG3_4ICR)) deallocate(XG3_4ICR)
      allocate( XG3_4ICR(ng3),stat=istat )

      XG3_1ICR=0.0d+00
      XG3_2ICR=0.0d+00
      XG3_3ICR=0.0d+00
      XG3_4ICR=0.0d+00

      GM3_1MPI=0.0d+00
      GM3_2MPI=0.0d+00
      GM3_3MPI=0.0d+00
      GM3_4MPI=0.0d+00

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

         call RXCHFmult_thread_gam3_IC1ex(istart,iend,ng3_seg,ng3,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        loop_map,XG3_1ICR,XG3_2ICR,
     x                        XG3_3ICR,XG3_4ICR,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------)
!-----CLEAN-UP-MEMORY-------------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

      wtime2 = MPI_WTIME() - wtime
      write(*,2000)rank,wtime2

! Here, a contiguous block of ng3/nproc integrals are stored in X* arrs
! Pass this to a ng3/nproc-dimensional array to prepare for allgather
      call copy_arr(ng3/nproc,XG3_1ICR(arrstart),GM3_1MPI)
      call copy_arr(ng3/nproc,XG3_2ICR(arrstart),GM3_2MPI)
      call copy_arr(ng3/nproc,XG3_3ICR(arrstart),GM3_3MPI)
      call copy_arr(ng3/nproc,XG3_4ICR(arrstart),GM3_4MPI)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "ng3:",ng3
!      do istart=1,ng3
!       write(*,9001) XG3_1ICR(istart),XG3_2ICR(istart),
!     x               XG3_3ICR(istart),XG3_4ICR(istart)
!      end do
!      write(*,*) "ng3/nproc:",ng3/nproc
!      do istart=1,ng3/nproc
!       write(*,9001) GM3_1MPI(istart),GM3_2MPI(istart),
!     x               GM3_3MPI(istart),GM3_4MPI(istart)
!      end do

! Pass and broadcast these major chunks of arrays to all processes
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHER(GM3_1MPI(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   XG3_1ICR(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM3_1 allgather"
      call MPI_ALLGATHER(GM3_2MPI(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   XG3_2ICR(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM3_2 allgather"
      call MPI_ALLGATHER(GM3_3MPI(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   XG3_3ICR(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM3_3 allgather"
      call MPI_ALLGATHER(GM3_4MPI(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   XG3_4ICR(1),ng3/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM3_4 allgather"

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "ng3 second:",ng3
!      do istart=1,ng3
!       write(*,9001) XG3_1ICR(istart),XG3_2ICR(istart),
!     x               XG3_3ICR(istart),XG3_4ICR(istart)
!      end do

! At this stage, ng3%nproc last elements are missing from arrays
      if (mod(ng3,nproc).ne.0) then

         wtime = MPI_WTIME()

! Have threads chop calculation of ng3%nproc integrals
         call loop_size(nproc*(ng3/nproc)+1,ng3,1,0,
     x                  istart,iend)
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
                        if (Loopi.eq.1) then
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

         call RXCHFmult_thread_gam3_IC1ex(istart,iend,ng3_seg,ng3,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        loop_map,XG3_1ICR,XG3_2ICR,
     x                        XG3_3ICR,XG3_4ICR,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

!-----CLEAN-UP-MEMORY-------------------------------------------------(
         if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

         wtime2 = MPI_WTIME() - wtime
         if (rank.eq.0) write(*,3000) wtime2

!         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!         write(*,*) "ng3+resid:",ng3
!         do istart=1,ng3
!         write(*,9001) XG3_1ICR(istart),XG3_2ICR(istart),
!     x                 XG3_3ICR(istart),XG3_4ICR(istart)
!         end do

      end if ! resid exists

!--------------------SYMMETRIZE----------------------------------------(
      wtime = MPI_WTIME() 

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

C  GAM_3 Symmetrization:
C  Determine packing indices for XGAM_3 integral matrices

C              As Packed-->       XGAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
!                            read(905,REC=ia_123) x123
                             ! x123=GAM_3(ia_123)

C              As Packed-->       XGAM_3(je2,ie2,je3,ie3,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec3,jec3,iec2,jec2,ia)
                             ia_132=ia
!                            read(905,REC=ia_132) x132
                             ! x132=GAM_3(ia_132)

C              As Packed-->       XGAM_3(je3,ie3,je1,ie1,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,iec3,jec3,ia)
                             ia_213=ia
!                            read(905,REC=ia_213) x213
                             ! x213=GAM_3(ia_213)

C              As Packed-->       XGAM_3(je1,ie1,je3,ie3,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec3,jec3,iec1,jec1,ia)
                             ia_231=ia
!                            read(905,REC=ia_231) x231
                             ! x231=GAM_3(ia_231)

C              As Packed-->       XGAM_3(je2,ie2,je1,ie1,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec3,jec3,iec1,jec1,iec2,jec2,ia)
                             ia_312=ia
!                            read(905,REC=ia_312) x312
                             ! x312=GAM_3(ia_312)

C              As Packed-->       XGAM_3(je1,ie1,je2,ie2,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec3,jec3,iec2,jec2,iec1,jec1,ia)
                             ia_321=ia
!                            read(905,REC=ia_321) x321
                             ! x321=GAM_3(ia_321)

C RXCHFmult(
C    index 1: regular electron
C    index 2: special electron 1
C    index 3: special electron 2
C    index 4: proton
C )

C Fully symmetrized integrals in GM3_1ICR
                       x123=XG3_1ICR(ia_123)
                       x132=XG3_1ICR(ia_132)
                       x213=XG3_1ICR(ia_213)
                       x231=XG3_1ICR(ia_231)
                       x312=XG3_1ICR(ia_312)
                       x321=XG3_1ICR(ia_321)
                       xxxx=(x123+x132+x213+x231+x312+x321)/six

                       GM3_1ICR(ia_123)=xxxx 
                       GM3_1ICR(ia_132)=xxxx 
                       GM3_1ICR(ia_213)=xxxx 
                       GM3_1ICR(ia_231)=xxxx 
                       GM3_1ICR(ia_312)=xxxx 
                       GM3_1ICR(ia_321)=xxxx 

C Partially symmetrized integrals in GM3_2ICR
                       x123=XG3_2ICR(ia_123)
                       x132=XG3_2ICR(ia_132)
                       xxxx=(x123+x132)*half

                       GM3_2ICR(ia_123)=xxxx 
                       GM3_2ICR(ia_132)=xxxx 

C Partially symmetrized integrals in GM3_3ICR
                       x123=XG3_3ICR(ia_123)
                       x213=XG3_3ICR(ia_213)
                       xxxx=(x123+x213)*half

                       GM3_3ICR(ia_123)=xxxx 
                       GM3_3ICR(ia_213)=xxxx 

C Partially symmetrized integrals in GM3_4ICR
                       x123=XG3_4ICR(ia_123)
                       x321=XG3_4ICR(ia_321)
                       xxxx=(x123+x321)*half

                       GM3_4ICR(ia_123)=xxxx 
                       GM3_4ICR(ia_321)=xxxx 

               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      wtime2 = MPI_WTIME() - wtime

      if(allocated(XG3_4ICR)) deallocate(XG3_4ICR)
      if(allocated(XG3_3ICR)) deallocate(XG3_3ICR)
      if(allocated(XG3_2ICR)) deallocate(XG3_2ICR)
      if(allocated(XG3_1ICR)) deallocate(XG3_1ICR)

      if (rank.eq.0) write(*,4000) wtime2
!--------------------SYMMETRIZE----------------------------------------)


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

 2000 FORMAT(8X,'    TIME FOR PROCESS ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,'    TIME FOR RESIDUAL ',1X,F10.2)

 4000 FORMAT(8X,'      TIME TO SYMMETRIZE INTEGRALS:',1X,F12.4/)

 9001 FORMAT(1X,4(F20.10))

      return
      end

