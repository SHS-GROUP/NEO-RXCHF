!=======================================================================
      subroutine RXCHF_GAM3_INT_MPI(nproc,rank,
     x                              nblocks,blockrank,
     x                              blockstart,blockend,
     x                              Nchunks,
     x                              nebf,npebf,nebfBE,npebfBE,npbf,
     x                              ng3,ng3loc,ng3prm,nat,ngtg1,
     x                              pmass,cat,zan,bcoef1,gamma1,
     x                              KPESTR,KPEEND,
     x                              AMPEB2C,AGEBFCC,
     x                              ELCEX,ELCAM,ELCBFC,
     x                              KPESTR_be,KPEEND_be,
     x                              AMPEB2C_be,AGEBFCC_be,
     x                              ELCEX_be,ELCAM_be,ELCBFC_be,
     x                              AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                              INT_GAM3)

!=======================================================================
      implicit none
      include 'mpif.h'
      include 'omp_lib.h'

! Input Variables
      integer nblocks
      integer blockrank
      integer blockstart,blockend
      integer Nchunks
      integer ng3             ! Total number of integrals
      integer ng3loc          ! Number of integrals for MPI proc to calc
      integer nebf,npebf,npbf,ng3prm
      integer npebfBE! Number primitive special electronic basis functions
      integer nebfBE ! Number contracted special electronic basis functions
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
! Special electron basis
      integer ELCAM_be(npebfBE,3)           ! 
      double precision ELCEX_be(npebfBE)    ! 
      double precision ELCBFC_be(npebfBE,3) ! 
      integer AMPEB2C_be(npebfBE)           ! Analogs for special electron basis
      double precision AGEBFCC_be(npebfBE)  ! 
      integer KPESTR_be(nebfBE)             ! 
      integer KPEEND_be(nebfBE)             ! 
!-------Basis Set Info-------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)

! Variables Returned
      double precision INT_GAM3(ng3loc)  ! INT  OMG3  integrals (partial symm)

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,i
      integer,allocatable :: loop_map(:,:)

      integer nproc,rank
#if MPI32
      integer*4 ierr
#else
      integer ierr
#endif
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

      INT_GAM3=0.0d+00

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
               do iec2=1,nebfBE
               do jec2=1,nebfBE
                  do iec3=1,nebfBE
                  do jec3=1,nebfBE

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
                           call RXCHFmult_GAM_3PK(nebf,nebfBE,npbf,
     x                                            ip,jp,
     x                                            iec1,jec1,
     x                                            iec2,jec2,
     x                                            iec3,jec3,arrstart)
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

         call RXCHF_GAM3_thread_INT_MPI(istart,iend,ng3_seg,ng3loc,
     x                         nebf,npebf,nebfBE,npebfBE,npbf,nat,ngtg1,
     x                         pmass,cat,zan,bcoef1,gamma1,
     x                         loop_map,arrstart,INT_GAM3,
     x                         KPESTR,KPEEND,
     x                         AMPEB2C,AGEBFCC,
     x                         ELCEX,ELCAM,ELCBFC,
     x                         KPESTR_be,KPEEND_be,
     x                         AMPEB2C_be,AGEBFCC_be,
     x                         ELCEX_be,ELCAM_be,ELCBFC_be,
     x                         AGNBFCC,NUCEX,NUCAM,NUCBFC)

         call RXCHF_writeint_mpi(nproc,rank,ng3loc,17,
     x                           "INT_GAM3-"//istring//".ufm",
     x                           INT_GAM3)

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

      wtime = MPI_WTIME() 

      if (rank.eq.0) write(*,*) "Symmetrizing INT_GAM3"
      call gam3_symm_132_restr(nproc,rank,
     x                         nebf,nebfBE,npbf,ng3,ng3loc,
     x                         mpistart,mpiend,INT_GAM3)

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
     x        8X,'   NUMBER OF 4-PARTICLE INTEGRALS: ',1X,I15/
     x        8X,'  NUMBER OF BLOCKS (USER DEFINED): ',1X,I15/
     x        8X,'                          ',/,
     x        8X,'  COMPUTATIONAL RESOURCES:',/,
     x        8X,'  ------------------------',/)

 1500 FORMAT( 8X,'      MPI PROCESSES:',1X,I4/
     x        8X,'        OMP THREADS:',1X,I4/)

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
C=======================================================================
      subroutine RXCHF_GAM3_thread_INT_MPI(istart,iend,ng3_seg,ng3,
     x                         nebf,npebf,nebfBE,npebfBE,npbf,nat,ngtg1,
     x                         pmass,cat,zan,bcoef1,gamma1,
     x                         loop_map,arrstart,INT_GAM3,
     x                         KPESTR,KPEEND,
     x                         AMPEB2C,AGEBFCC,
     x                         ELCEX,ELCAM,ELCBFC,
     x                         KPESTR_be,KPEEND_be,
     x                         AMPEB2C_be,AGEBFCC_be,
     x                         ELCEX_be,ELCAM_be,ELCBFC_be,
     x                         AGNBFCC,NUCEX,NUCAM,NUCBFC)

C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input Variables
      integer istart,iend,ng3_seg
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npebfBE! Number primitive special electronic basis functions
      integer nebfBE ! Number contracted special electronic basis functions
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
! Special electron basis
      integer ELCAM_be(npebfBE,3)           ! 
      double precision ELCEX_be(npebfBE)    ! 
      double precision ELCBFC_be(npebfBE,3) ! 
      integer AMPEB2C_be(npebfBE)           ! Analogs for special electron basis
      double precision AGEBFCC_be(npebfBE)  ! 
      integer KPESTR_be(nebfBE)             ! 
      integer KPEEND_be(nebfBE)             ! 
C-------Basis Set Info-------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)
      integer loop_map(ng3_seg,8)

! Variables Returned
      double precision INT_GAM3(ng3)

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer iec3,jec3  !
      integer imap,ia
      double precision OMG3

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
!$ompx shared(INT_GAM3)
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
!$ompx private(OMG3)
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

         call RXCHFmult_contract_omega3_INT(ip,jp,iec1,jec1,iec2,jec2,
     x                            iec3,jec3,nebf,npebf,nebfBE,npebfBE,
     x                            npbf,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,
     x                            AMPEB2C,AGEBFCC,
     x                            ELCEX,ELCAM,ELCBFC,
     x                            KPESTR_be,KPEEND_be,
     x                            AMPEB2C_be,AGEBFCC_be,
     x                            ELCEX_be,ELCAM_be,ELCBFC_be,
     x                            AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                            OMG3)

         call RXCHFmult_GAM_3PK(nebf,nebfBE,npbf,
     x                          ip,jp,
     x                          iec1,jec1,
     x                          iec2,jec2,
     x                          iec3,jec3,ia)

         INT_GAM3(ia-arrstart+1)=OMG3

      end do
!$omp end do
!$omp end parallel      
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

      return
      end
C=======================================================================
      subroutine gam3_symm_132_restr(nproc,rank,
     x                               nebf,nebfBE,npbf,ntotints,nlocints,
     x                               mpistart,mpiend,GAM3)

C Symmetrizes INT_GAM3 integrals averaging
C   < ip ie1 ie2 ie3 | O3 | ip je1 je2 je3 > = Omega_3(1,2,3)
C   < ip ie1 ie3 ie2 | O3 | ip je1 je3 je2 > = Omega_3(1,3,2)
C=======================================================================
      implicit none
      include "mpif.h"

! Input variables
      integer nproc,rank
      integer nebf,nebfBE,npbf
      integer ntotints,nlocints
      integer mpistart,mpiend

! Input/output variables
      double precision GAM3(nlocints)

! Local variables
      integer i
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3
      integer ia_first,ia_last
      integer ia,ia_12,ia_21
#if MPI32
      integer*4 tag,reqscount,ierr
      integer*4 sendrank,recvrank,rank_
      integer*4 reqs(2*nlocints)
#else
      integer tag,reqscount,ierr
      integer sendrank,recvrank,rank_
      integer reqs(2*nlocints)
#endif
      double precision x12,x21
      double precision xGAM3(nlocints)

#if MPI32
      rank_=int(rank,kind=4)
#else
      rank_=rank
#endif

      xGAM3=0.0d+00

! Get appropriate ia_ji values for xGAM3

      do ip=1,npbf
      do jp=1,npbf
         do ie1=1,nebf
         do je1=1,nebf

            reqscount=0

            do ie2=1,nebfBE
            do je2=1,nebfBE
               do ie3=1,nebfBE
               do je3=1,nebfBE

                 call RXCHFmult_GAM_3PK(nebf,nebfBE,npbf,
     x                                  ip,jp,
     x                                  ie1,je1,
     x                                  ie2,je2,
     x                                  ie3,je3,ia_12)
                 call RXCHFmult_GAM_3PK(nebf,nebfBE,npbf,
     x                                  ip,jp,
     x                                  ie1,je1,
     x                                  ie3,je3,
     x                                  ie2,je2,ia_21)

       if ((ia_12.ge.mpistart).and.(ia_12.le.mpiend)) then

         call get_mpi_proc(ntotints,nproc,ia_21,recvrank)

         if (recvrank.eq.rank_) then
          xGAM3(ia_12-mpistart+1)=GAM3(ia_21-mpistart+1)
         else
          call pack_4D(nebfBE,nebfBE,nebfBE,
     x                 je3,ie3,je2,ie2,ia)
#if MPI32
          tag=int(ia,kind=4)
#else
          tag=ia
#endif
          reqscount=reqscount+1
          call MPI_IRECV(xGAM3(ia_12-mpistart+1),1,
     x                   MPI_DOUBLE_PRECISION,recvrank,tag,
     x                   MPI_COMM_WORLD,reqs(reqscount),ierr)
         end if

       end if

       if ((ia_21.ge.mpistart).and.(ia_21.le.mpiend)) then

         call get_mpi_proc(ntotints,nproc,ia_12,sendrank)

         if (sendrank.ne.rank_) then
          call pack_4D(nebfBE,nebfBE,nebfBE,
     x                 je3,ie3,je2,ie2,ia)
#if MPI32
          tag=int(ia,kind=4)
#else
          tag=ia
#endif
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

