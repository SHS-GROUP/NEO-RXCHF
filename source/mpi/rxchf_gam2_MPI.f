!=======================================================================
      subroutine RXCHF_GAM2_MPI(nproc,rank,
     x                          Nchunks,nebf,npebf,npbf,
     x                          ng2,ng2prm,nat,ngtg1,
     x                          pmass,cat,zan,bcoef1,gamma1,
     x                          KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                          ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                          GM2_1ICR,GM2_2ICR,GM2sICR)

!=======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ng2,nebf,npebf,npbf,ng2prm
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
      double precision GM2_1ICR(ng2)  ! XCHF OMG2  integrals (symm)
      double precision GM2_2ICR(ng2)  ! INT  OMG2  integrals (unsymm)
      double precision GM2sICR(ng2)   ! XCHF OMG2s integrals (symm)

! Local Variables
      integer istat,ichunk,istart,iend,ng2_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2
      integer,allocatable :: loop_map(:,:)

      double precision,allocatable :: XG2_1ICR(:),XG2_2ICR(:),XG2sICR(:)

      double precision GM2_1MPI(ng2/nproc)  ! XCHF OMG2   integrals (symm)
      double precision GM2_2MPI(ng2/nproc)  ! INT  OMG2   integrals (unsymm)
      double precision GM2sMPI(ng2/nproc)   ! XCHF OMG2s  integrals (symm)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

      integer ia_12
      integer ia_21

      double precision x12,y12
      double precision x21,y21
      double precision xx
      double precision yy

      double precision zero,half,two
      parameter(zero=0.0d+00,half=0.5d+00,two=2.0d+00)

      double precision wtime
      double precision wtime2

! Have each process calculate ng2/nproc integrals according to rank then allgather
! Have each process calculate ng2%nproc remaining integrals
!  - saves waiting for master and eliminates one broadcast

      call get_mpi_range(ng2,nproc,rank,mpistart,mpiend)

      if (rank.eq.0) then
       write(*,1000) ng2,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      if(allocated(XG2_1ICR)) deallocate(XG2_1ICR)
      allocate( XG2_1ICR(ng2),stat=istat )
      if(allocated(XG2_2ICR)) deallocate(XG2_2ICR)
      allocate( XG2_2ICR(ng2),stat=istat )
      if(allocated(XG2sICR)) deallocate(XG2sICR)
      allocate( XG2sICR(ng2),stat=istat )

      XG2_1ICR=0.0d+00
      XG2_2ICR=0.0d+00
      XG2sICR=0.0d+00

      GM2_1MPI=0.0d+00
      GM2_2MPI=0.0d+00
      GM2sMPI=0.0d+00

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

         call RXCHFmult_thread_gam2_IC1(istart,iend,ng2_seg,ng2,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           loop_map,XG2_1ICR,XG2_2ICR,XG2sICR,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_2--------------------------------)
!-----CLEAN-UP-MEMORY-------------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

      wtime2 = MPI_WTIME() - wtime
      write(*,2000)rank,wtime2

! Here, a contiguous block of ng2/nproc integrals are stored in X* arrs
! Pass this to a ng2/nproc-dimensional array to prepare for allgather
      call copy_arr(ng2/nproc,XG2_1ICR(arrstart),GM2_1MPI)
      call copy_arr(ng2/nproc,XG2_2ICR(arrstart),GM2_2MPI)
      call copy_arr(ng2/nproc,XG2sICR(arrstart),GM2sMPI)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "ng2:",ng2
!      do istart=1,ng2
!       write(*,9001) XG2_1ICR(istart),XG2_2ICR(istart),
!     x               XG2sICR(istart)
!      end do
!      write(*,*) "ng2/nproc:",ng2/nproc
!      do istart=1,ng2/nproc
!       write(*,9001) GM2_1MPI(istart),GM2_2MPI(istart),
!     x               GM2sMPI(istart)
!      end do

! Pass and broadcast these major chunks of arrays to all processes
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHER(GM2_1MPI(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   XG2_1ICR(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM2_1 allgather"
      call MPI_ALLGATHER(GM2_2MPI(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   XG2_2ICR(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM2_2 allgather"
      call MPI_ALLGATHER(GM2sMPI(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   XG2sICR(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM2s allgather"

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "ng2 second:",ng2
!      do istart=1,ng2
!       write(*,9001) XG2_1ICR(istart),XG2_2ICR(istart),
!     x               XG2sICR(istart)
!      end do

! At this stage, ng2%nproc last elements are missing from arrays
      if (mod(ng2,nproc).ne.0) then

         wtime = MPI_WTIME()

! Have threads chop calculation of ng2%nproc integrals
         call loop_size(nproc*(ng2/nproc)+1,ng2,1,0,
     x                  istart,iend)
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
                      if (Loopi.eq.1) then
                         call index_GAM_2PK(nebf,npbf,
     x                            ip,jp,iec1,jec1,iec2,jec2,arrstart)
                      end if

                   end if

               end do
               end do
            end do
            end do
         end do
         end do

         call RXCHFmult_thread_gam2_IC1(istart,iend,ng2_seg,ng2,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        loop_map,XG2_1ICR,XG2_2ICR,XG2sICR,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

!-----CLEAN-UP-MEMORY-------------------------------------------------(
         if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

         wtime2 = MPI_WTIME() - wtime
         if (rank.eq.0) write(*,3000) wtime2

!         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!         write(*,*) "ng2+resid:",ng2
!         do istart=1,ng2
!         write(*,9001) XG2_1ICR(istart),XG2_2ICR(istart),
!     x                 XG2sICR(istart)
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

C Symmetrized integrals in GM2_1ICR (XCHF integrals)
                       x12=XG2_1ICR(ia_12)
                       x21=XG2_1ICR(ia_21)
                       GM2_1ICR(ia_12)=(x12+x21)/two

C Unsymmetrized integrals in GM2_2ICR (interaction integrals)
                       x12=XG2_2ICR(ia_12)
                       GM2_2ICR(ia_12)=x12

C Symmetrized integrals in GM2sICR (XCHF integrals)
                       x12=XG2sICR(ia_12)
                       x21=XG2sICR(ia_21)
                       GM2sICR(ia_12)=(x12+x21)/two

               end do
               end do
            end do
            end do
         end do
         end do

         wtime2 = MPI_WTIME() - wtime

         if(allocated(XG2sICR)) deallocate(XG2sICR)
         if(allocated(XG2_2ICR)) deallocate(XG2_2ICR)
         if(allocated(XG2_1ICR)) deallocate(XG2_1ICR)

         if (rank.eq.0) write(*,4000) wtime2
!--------------------SYMMETRIZE----------------------------------------)


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

 2000 FORMAT(8X,'    TIME FOR PROCESS ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,'    TIME FOR RESIDUAL ',1X,F10.2)

 4000 FORMAT(8X,'      TIME TO SYMMETRIZE INTEGRALS:',1X,F12.4/)

 9001 FORMAT(1X,3(F20.10))

      return
      end
!=======================================================================
      subroutine RXCHF_GAM2ex_MPI(nproc,rank,
     x                            Nchunks,nebf,npebf,npbf,
     x                            ng2,ng2prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            GM2_1ICR,GM2_2ICR,GM2_3ICR,GM2sICR)

!=======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer Nchunks
      integer ng2,nebf,npebf,npbf,ng2prm
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
      double precision GM2_1ICR(ng2)  ! XCHF OMG2   integrals (symm)
      double precision GM2_2ICR(ng2)  ! INT  OMG2   integrals (unsymm)
      double precision GM2_3ICR(ng2)  ! INT  OMG2ex integrals (symm)
      double precision GM2sICR(ng2)   ! XCHF OMG2s  integrals (symm)

! Local Variables
      integer istat,ichunk,istart,iend,ng2_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2
      integer,allocatable :: loop_map(:,:)

      double precision,allocatable :: XG2_1ICR(:),XG2_2ICR(:)
      double precision,allocatable :: XG2_3ICR(:),XG2sICR(:)

      double precision GM2_1MPI(ng2/nproc)  ! XCHF OMG2   integrals (symm)
      double precision GM2_2MPI(ng2/nproc)  ! INT  OMG2   integrals (unsymm)
      double precision GM2_3MPI(ng2/nproc)  ! INT  OMG2ex integrals (symm)
      double precision GM2sMPI(ng2/nproc)   ! XCHF OMG2s  integrals (symm)

      integer nproc,rank
      integer*4 ierr
      integer mpistart,mpiend,arrstart

      integer ia_12
      integer ia_21

      double precision x12,y12
      double precision x21,y21
      double precision xx
      double precision yy

      double precision zero,half,two
      parameter(zero=0.0d+00,half=0.5d+00,two=2.0d+00)

      double precision wtime
      double precision wtime2

! Have each process calculate ng2/nproc integrals according to rank then allgather
! Have each process calculate ng2%nproc remaining integrals
!  - saves waiting for master and eliminates one broadcast

      call get_mpi_range(ng2,nproc,rank,mpistart,mpiend)

      if (rank.eq.0) then
       write(*,1000) ng2,nchunks
       write(*,1500) nproc,omp_get_max_threads()
      end if

      if(allocated(XG2_1ICR)) deallocate(XG2_1ICR)
      allocate( XG2_1ICR(ng2),stat=istat )
      if(allocated(XG2_2ICR)) deallocate(XG2_2ICR)
      allocate( XG2_2ICR(ng2),stat=istat )
      if(allocated(XG2_3ICR)) deallocate(XG2_3ICR)
      allocate( XG2_3ICR(ng2),stat=istat )
      if(allocated(XG2sICR)) deallocate(XG2sICR)
      allocate( XG2sICR(ng2),stat=istat )

      XG2_1ICR=0.0d+00
      XG2_2ICR=0.0d+00
      XG2_3ICR=0.0d+00
      XG2sICR=0.0d+00

      GM2_1MPI=0.0d+00
      GM2_2MPI=0.0d+00
      GM2_3MPI=0.0d+00
      GM2sMPI=0.0d+00

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

         call RXCHFmult_thread_gam2_IC1ex(istart,iend,ng2_seg,ng2,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           loop_map,XG2_1ICR,XG2_2ICR,
     x                           XG2_3ICR,XG2sICR,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_2--------------------------------)
!-----CLEAN-UP-MEMORY-------------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

      wtime2 = MPI_WTIME() - wtime
      write(*,2000)rank,wtime2

! Here, a contiguous block of ng2/nproc integrals are stored in X* arrs
! Pass this to a ng2/nproc-dimensional array to prepare for allgather
      call copy_arr(ng2/nproc,XG2_1ICR(arrstart),GM2_1MPI)
      call copy_arr(ng2/nproc,XG2_2ICR(arrstart),GM2_2MPI)
      call copy_arr(ng2/nproc,XG2_3ICR(arrstart),GM2_3MPI)
      call copy_arr(ng2/nproc,XG2sICR(arrstart),GM2sMPI)

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "ng2:",ng2
!      do istart=1,ng2
!       write(*,9001) XG2_1ICR(istart),XG2_2ICR(istart),
!     x               XG2_3ICR(istart),XG2sICR(istart)
!      end do
!      write(*,*) "ng2/nproc:",ng2/nproc
!      do istart=1,ng2/nproc
!       write(*,9001) GM2_1MPI(istart),GM2_2MPI(istart),
!     x               GM2_3MPI(istart),GM2sMPI(istart)
!      end do

! Pass and broadcast these major chunks of arrays to all processes
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHER(GM2_1MPI(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   XG2_1ICR(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM2_1 allgather"
      call MPI_ALLGATHER(GM2_2MPI(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   XG2_2ICR(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM2_2 allgather"
      call MPI_ALLGATHER(GM2_3MPI(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   XG2_3ICR(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM2_3 allgather"
      call MPI_ALLGATHER(GM2sMPI(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   XG2sICR(1),ng2/nproc,MPI_DOUBLE_PRECISION,
     x                   MPI_COMM_WORLD,ierr)
      if (ierr.ne.0) write (*,*) "Trouble with GM2s allgather"

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      write(*,*) "ng2 second:",ng2
!      do istart=1,ng2
!       write(*,9001) XG2_1ICR(istart),XG2_2ICR(istart),
!     x               XG2_3ICR(istart),XG2sICR(istart)
!      end do

! At this stage, ng2%nproc last elements are missing from arrays
      if (mod(ng2,nproc).ne.0) then

         wtime = MPI_WTIME()

! Have threads chop calculation of ng2%nproc integrals
         call loop_size(nproc*(ng2/nproc)+1,ng2,1,0,
     x                  istart,iend)
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
                      if (Loopi.eq.1) then
                         call index_GAM_2PK(nebf,npbf,
     x                            ip,jp,iec1,jec1,iec2,jec2,arrstart)
                      end if

                   end if

               end do
               end do
            end do
            end do
         end do
         end do

         call RXCHFmult_thread_gam2_IC1ex(istart,iend,ng2_seg,ng2,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        loop_map,XG2_1ICR,XG2_2ICR,
     x                        XG2_3ICR,XG2sICR,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

!-----CLEAN-UP-MEMORY-------------------------------------------------(
         if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

         wtime2 = MPI_WTIME() - wtime
         if (rank.eq.0) write(*,3000) wtime2

!         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!         write(*,*) "ng2+resid:",ng2
!         do istart=1,ng2
!         write(*,9001) XG2_1ICR(istart),XG2_2ICR(istart),
!     x                 XG2_3ICR(istart),XG2sICR(istart)
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

C Symmetrized integrals in GM2_1ICR (XCHF integrals)
                       x12=XG2_1ICR(ia_12)
                       x21=XG2_1ICR(ia_21)
                       GM2_1ICR(ia_12)=(x12+x21)/two

C Unsymmetrized integrals in GM2_2ICR (interaction integrals)
                       x12=XG2_2ICR(ia_12)
                       GM2_2ICR(ia_12)=x12

C Symmetrized integrals in GM2_3ICR (exchange integrals)
                       x12=XG2_3ICR(ia_12)
                       x21=XG2_3ICR(ia_21)
                       GM2_3ICR(ia_12)=(x12+x21)/two

C Symmetrized integrals in GM2sICR (XCHF integrals)
                       x12=XG2sICR(ia_12)
                       x21=XG2sICR(ia_21)
                       GM2sICR(ia_12)=(x12+x21)/two

               end do
               end do
            end do
            end do
         end do
         end do

         wtime2 = MPI_WTIME() - wtime

         if(allocated(XG2sICR)) deallocate(XG2sICR)
         if(allocated(XG2_3ICR)) deallocate(XG2_3ICR)
         if(allocated(XG2_2ICR)) deallocate(XG2_2ICR)
         if(allocated(XG2_1ICR)) deallocate(XG2_1ICR)

         if (rank.eq.0) write(*,4000) wtime2
!--------------------SYMMETRIZE----------------------------------------)


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

 2000 FORMAT(8X,'    TIME FOR PROCESS ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,'    TIME FOR RESIDUAL ',1X,F10.2)

 4000 FORMAT(8X,'      TIME TO SYMMETRIZE INTEGRALS:',1X,F12.4/)

 9001 FORMAT(1X,4(F20.10))

      return
      end

