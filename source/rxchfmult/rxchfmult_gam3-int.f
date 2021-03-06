!=======================================================================
      subroutine RXCHFmult_GAM3_INT(Nchunks,
     x                              nebf,npebf,nebfBE,npebfBE,npbf,
     x                              ng3,nat,ngtg1,
     x                              pmass,cat,zan,bcoef1,gamma1,
     x                              KPESTR,KPEEND,
     x                              AMPEB2C,AGEBFCC,
     x                              ELCEX,ELCAM,ELCBFC,
     x                              KPESTR_be,KPEEND_be,
     x                              AMPEB2C_be,AGEBFCC_be,
     x                              ELCEX_be,ELCAM_be,ELCBFC_be,
     x                              AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                              INT_GAM3)

C Calculates four-particle interaction integrals
!=======================================================================
      implicit none
      include 'omp_lib.h'
! Input Variables
      integer Nchunks
      integer ng3,nebf,npebf,nebfBE,npebfBE,npbf
      integer nat,ngtg1
!-------Basis Set Info-------(
      integer ELCAM(npebf,3)                ! Angular mom for electrons
      integer NUCAM(npbf,3)                 ! Angular mom for quantum nuclei
      double precision ELCEX(npebf)         ! Exponents: elec basis
      double precision NUCEX(npbf)          ! Exponents: nuc basis
      double precision ELCBFC(npebf,3)      ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)       ! basis centers: nuc basis
      integer AMPEB2C(npebf)                ! Map primitive index to contracted
      double precision AGEBFCC(npebf)       ! Map prim index to contract coef
      double precision AGNBFCC(npbf)        ! Nuclear contract coef
      integer KPESTR(nebf)                  ! Map contracted index to primitive start
      integer KPEEND(nebf)                  ! Map contracted index to primitive end
! Special electron basis
      integer ELCAM_be(npebfBE,3)           ! 
      double precision ELCEX_be(npebfBE)    ! 
      double precision ELCBFC_be(npebfBE,3) ! 
      integer AMPEB2C_be(npebfBE)           ! Analogs for special electron basis
      double precision AGEBFCC_be(npebfBE)  ! 
      integer KPESTR_be(nebfBE)             ! 
      integer KPEEND_be(nebfBE)             ! 
!-------Basis Set Info-------)
      double precision pmass                ! Mass of nonelectron quantum particle 
      double precision zan(nat)             ! Classical nuclear charges
      double precision cat(3,nat)           ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)

! Variables Returned
      double precision INT_GAM3(ng3)  ! INT  OMG3  integrals (partial symm)

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
      double precision,allocatable :: XG3ICR(:)

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


         write(*,1000) ng3,nchunks,omp_get_num_procs(),
     xomp_get_max_threads(),1
         wtime = omp_get_wtime()

         if(allocated(XG3ICR)) deallocate(XG3ICR)
         allocate( XG3ICR(ng3),stat=istat )

!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      do ichunk=1,Nchunks

         wtime2 = omp_get_wtime()

         call loop_size(1,ng3,Nchunks,ichunk-1,istart,iend)
         ng3_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )

C RXCHFmult(
C    index 1: regular electron   (iec1,jec1)
C    index 2: special electron 1 (iec2,jec2)
C    index 3: special electron 2 (iec3,jec3)
C    index 4: proton             (  ip,jp  )
C )

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
                     end if

                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

         call RXCHFmult_thread_gam3_INT(istart,iend,ng3_seg,ng3,
     x                        nebf,npebf,nebfBE,npebfBE,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        loop_map,XG3ICR,
     x                        KPESTR,KPEEND,
     x                        AMPEB2C,AGEBFCC,
     x                        ELCEX,ELCAM,ELCBFC,
     x                        KPESTR_be,KPEEND_be,
     x                        AMPEB2C_be,AGEBFCC_be,
     x                        ELCEX_be,ELCAM_be,ELCBFC_be,
     x                        AGNBFCC,NUCEX,NUCAM,NUCBFC)


         wtime2 = omp_get_wtime() - wtime2
         write(*,2000)ichunk,wtime2

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------)

!-----CLEAN-UP-MEMORY-------------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

!--------------------SYMMETRIZE----------------------------------------(

      wtime = omp_get_wtime() - wtime
      write(*,3000)wtime
      wtime = omp_get_wtime() 

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebfBE
            do jec2=1,nebfBE
               do iec3=1,nebfBE
               do jec3=1,nebfBE

C  GAM_3 Symmetrization:
C  Determine packing indices for XGAM_3 integral matrices

C              As Packed-->       XGAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) 
                 call RXCHFmult_GAM_3PK(nebf,nebfBE,npbf,
     x                                  ip,jp,
     x                                  iec1,jec1,
     x                                  iec2,jec2,
     x                                  iec3,jec3,ia)

                 ia_123=ia

C              As Packed-->       XGAM_3(je2,ie2,je3,ie3,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2) 
                 call RXCHFmult_GAM_3PK(nebf,nebfBE,npbf,
     x                                  ip,jp,
     x                                  iec1,jec1,
     x                                  iec3,jec3,
     x                                  iec2,jec2,ia)
                 ia_132=ia

C Partially symmetrized integrals in INT_GAM3
                 x123=XG3ICR(ia_123)
                 x132=XG3ICR(ia_132)
                 xxxx=(x123+x132)*half

                 INT_GAM3(ia_123)=xxxx 
                 INT_GAM3(ia_132)=xxxx 

               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      wtime = omp_get_wtime() - wtime

      if(allocated(XG3ICR)) deallocate(XG3ICR)

      write(*,4000) wtime
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
     x        8X,'  ------------------------',/,
     x        8X,'     CORES PER NODE:',1X,I3/
     x        8X,'          AVAILABLE:',1X,I3/
     x        8X,'    NUMBER OF NODES:',1X,I3/)

 2000 FORMAT(8X,'    TIME TO EVALUATE BLOCK ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,'  TIMING SUMMARY FOR 4-PARTICLE INTEGRALS:',/,
     x        8X,'  ----------------------------------------',/,
     x        8X,'    TIME TO EVALUATE ALL INTEGRALS:',1X,F12.4)

 4000 FORMAT(8X,'      TIME TO SYMMETRIZE INTEGRALS:',1X,F12.4/)


      return
      end
C=======================================================================
      subroutine RXCHFmult_thread_gam3_INT(istart,iend,ng3_seg,ng3,
     x                         nebf,npebf,nebfBE,npebfBE,npbf,nat,ngtg1,
     x                         pmass,cat,zan,bcoef1,gamma1,
     x                         loop_map,XG3ICR,
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
      integer npebfBE! Number primitive electronic basis functions
      integer nebfBE ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas
      integer ng3

      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer iec3,jec3  !
!-------Basis Set Info-------(
      integer ELCAM(npebf,3)                ! Angular mom for electrons
      integer NUCAM(npbf,3)                 ! Angular mom for quantum nuclei
      double precision ELCEX(npebf)         ! Exponents: elec basis
      double precision NUCEX(npbf)          ! Exponents: nuc basis
      double precision ELCBFC(npebf,3)      ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)       ! basis centers: nuc basis
      integer AMPEB2C(npebf)                ! Map primitive index to contracted
      double precision AGEBFCC(npebf)       ! Map prim index to contract coef
      double precision AGNBFCC(npbf)        ! Nuclear contract coef
      integer KPESTR(nebf)                  ! Map contracted index to primitive start
      integer KPEEND(nebf)                  ! Map contracted index to primitive end
! Special electron basis
      integer ELCAM_be(npebfBE,3)           ! 
      double precision ELCEX_be(npebfBE)    ! 
      double precision ELCBFC_be(npebfBE,3) ! 
      integer AMPEB2C_be(npebfBE)           ! Analogs for special electron basis
      double precision AGEBFCC_be(npebfBE)  ! 
      integer KPESTR_be(nebfBE)             ! 
      integer KPEEND_be(nebfBE)             ! 
!-------Basis Set Info-------)
      double precision pmass                ! Mass of nonelectron quantum particle 
      double precision zan(nat)             ! Classical nuclear charges
      double precision cat(3,nat)           ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)
      integer loop_map(ng3_seg,8)
      double precision XG3ICR(ng3)

! Variables Returned

! Local Variables
      integer imap,ia
      double precision OMG3

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!!!!!!!!!!!!$ompx shared(GAM_3)
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(XG3ICR)
!$ompx shared(KPESTR,KPEEND)
!$ompx shared(AMPEB2C,AGEBFCC)
!$ompx shared(ELCEX,ELCAM,ELCBFC)
!$ompx shared(KPESTR_be,KPEEND_be)
!$ompx shared(AMPEB2C_be,AGEBFCC_be)
!$ompx shared(ELCEX_be,ELCAM_be,ELCBFC_be)
!$ompx shared(AGNBFCC,NUCEX,NUCAM,NUCBFC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,nebfBE,npebfBE,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(istart,iend)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia) 
!$ompx private(OMG3)
!$ompx private(id)

!     id= omp_get_thread_num()
!     write(*,*)' Hello from process ',id
!     if(id.eq.0) then
!        write(*,*)'Threads in use', omp_get_num_threads()
!     end if

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

         call RXCHFmult_contract_omega3_INT(ip,jp,
     x                            iec1,jec1,iec2,jec2,iec3,jec3,
     x                            nebf,npebf,nebfBE,npebfBE,
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

!        As Packed--> XGAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
!        XGAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) 
         call RXCHFmult_GAM_3PK(nebf,nebfBE,npbf,
     x                          ip,jp,
     x                          iec1,jec1,
     x                          iec2,jec2,
     x                          iec3,jec3,ia)
         XG3ICR(ia)=OMG3

      end do
!$omp end do
!$omp end parallel      
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

C---OPENMP-TIMING------------------------------------------------------(
!     wtime = omp_get_wtime() - wtime
!     write(*,*)'TIME TO CALCULATE GAM_3 INTEGRALS: ',wtime
!     write(*,*)'ISTART=',istart
!     write(*,*)'IEND  =',iend
C---OPENMP-TIMING------------------------------------------------------)


      return
      end
C======================================================================
      subroutine RXCHFmult_contract_omega3_INT(ip,jp,
     x                            iec1,jec1,iec2,jec2,iec3,jec3,
     x                            nebf,npebf,nebfBE,npebfBE,
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

C======================================================================
      implicit none
c     include 'omp_lib.h'
c     include 'mpif.h'

C Input Variables
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npebfBE! Number primitive electronic basis functions
      integer nebfBE ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas

      integer ip,jp
      integer iec1,jec1                     !
      integer iec2,jec2                     ! Contracted elec basis function indices
      integer iec3,jec3                     !
!-------Basis Set Info-------(
      integer ELCAM(npebf,3)                ! Angular mom for electrons
      integer NUCAM(npbf,3)                 ! Angular mom for quantum nuclei
      double precision ELCEX(npebf)         ! Exponents: elec basis
      double precision NUCEX(npbf)          ! Exponents: nuc basis
      double precision ELCBFC(npebf,3)      ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)       ! basis centers: nuc basis
      integer AMPEB2C(npebf)                ! Map primitive index to contracted
      double precision AGEBFCC(npebf)       ! Map prim index to contract coef
      double precision AGNBFCC(npbf)        ! Nuclear contract coef
      integer KPESTR(nebf)                  ! Map contracted index to primitive start
      integer KPEEND(nebf)                  ! Map contracted index to primitive end
! Special electron basis
      integer ELCAM_be(npebfBE,3)           ! 
      double precision ELCEX_be(npebfBE)    ! 
      double precision ELCBFC_be(npebfBE,3) ! 
      integer AMPEB2C_be(npebfBE)           ! Analogs for special electron basis
      double precision AGEBFCC_be(npebfBE)  ! 
      integer KPESTR_be(nebfBE)             ! 
      integer KPEEND_be(nebfBE)             ! 
!-------Basis Set Info-------)
      double precision pmass                ! Mass of nonelectron quantum particle 
      double precision zan(nat)             ! Classical nuclear charges
      double precision cat(3,nat)           ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)

C Variables Returned
      double precision OMG3

C Local Variables
      integer ie1,je1
      integer ie2,je2
      integer ie3,je3
      integer ie1_start
      integer ie2_start
      integer ie3_start
      integer je1_start
      integer je2_start
      integer je3_start
      integer ie1_end
      integer ie2_end
      integer ie3_end
      integer je1_end
      integer je2_end
      integer je3_end

      double precision Cof_ie1,Cof_je1
      double precision Cof_ie2,Cof_je2
      double precision Cof_ie3,Cof_je3
      double precision Cof_ip,Cof_jp
C--------------------------------(
C Basis set-related local variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1,Amat1(3) 
      double precision A2,Amat2(3) 
      double precision A3,Amat3(3) 
      double precision A4,Amat4(3) 
      double precision B1,Bmat1(3) 
      double precision B2,Bmat2(3) 
      double precision B3,Bmat3(3) 
      double precision B4,Bmat4(3) 
C--------------------------------)
      double precision ans


      ie1_start=KPESTR(iec1)
      ie2_start=KPESTR_be(iec2)
      ie3_start=KPESTR_be(iec3)

      je1_start=KPESTR(jec1)
      je2_start=KPESTR_be(jec2)
      je3_start=KPESTR_be(jec3)

      ie1_end=KPEEND(iec1)
      ie2_end=KPEEND_be(iec2)
      ie3_end=KPEEND_be(iec3)

      je1_end=KPEEND(jec1)
      je2_end=KPEEND_be(jec2)
      je3_end=KPEEND_be(jec3)

      OMG3=0.0d+00

      do ie1=ie1_start,ie1_end
      do je1=je1_start,je1_end

        do ie2=ie2_start,ie2_end
        do je2=je2_start,je2_end

          do ie3=ie3_start,ie3_end
          do je3=je3_start,je3_end
               
C             Get Basis set info:
              A1=ELCEX(ie1)
              I1=ELCAM(ie1,1)
              J1=ELCAM(ie1,2)
              K1=ELCAM(ie1,3)
              Amat1(1)=ELCBFC(ie1,1)
              Amat1(2)=ELCBFC(ie1,2)
              Amat1(3)=ELCBFC(ie1,3)

              A2=ELCEX_be(ie2)
              I2=ELCAM_be(ie2,1)
              J2=ELCAM_be(ie2,2)
              K2=ELCAM_be(ie2,3)
              Amat2(1)=ELCBFC_be(ie2,1)
              Amat2(2)=ELCBFC_be(ie2,2)
              Amat2(3)=ELCBFC_be(ie2,3)

              A3=ELCEX_be(ie3)
              I3=ELCAM_be(ie3,1)
              J3=ELCAM_be(ie3,2)
              K3=ELCAM_be(ie3,3)
              Amat3(1)=ELCBFC_be(ie3,1)
              Amat3(2)=ELCBFC_be(ie3,2)
              Amat3(3)=ELCBFC_be(ie3,3)

              A4=NUCEX(ip)
              I4=NUCAM(ip,1)
              J4=NUCAM(ip,2)
              K4=NUCAM(ip,3)
              Amat4(1)=NUCBFC(ip,1)
              Amat4(2)=NUCBFC(ip,2)
              Amat4(3)=NUCBFC(ip,3)

              B1=ELCEX(je1)
              L1=ELCAM(je1,1)
              M1=ELCAM(je1,2)
              N1=ELCAM(je1,3)
              Bmat1(1)=ELCBFC(je1,1)
              Bmat1(2)=ELCBFC(je1,2)
              Bmat1(3)=ELCBFC(je1,3)

              B2=ELCEX_be(je2)
              L2=ELCAM_be(je2,1)
              M2=ELCAM_be(je2,2)
              N2=ELCAM_be(je2,3)
              Bmat2(1)=ELCBFC_be(je2,1)
              Bmat2(2)=ELCBFC_be(je2,2)
              Bmat2(3)=ELCBFC_be(je2,3)

              B3=ELCEX_be(je3)
              L3=ELCAM_be(je3,1)
              M3=ELCAM_be(je3,2)
              N3=ELCAM_be(je3,3)
              Bmat3(1)=ELCBFC_be(je3,1)
              Bmat3(2)=ELCBFC_be(je3,2)
              Bmat3(3)=ELCBFC_be(je3,3)

              B4=NUCEX(jp)
              L4=NUCAM(jp,1)
              M4=NUCAM(jp,2)
              N4=NUCAM(jp,3)
              Bmat4(1)=NUCBFC(jp,1)
              Bmat4(2)=NUCBFC(jp,2)
              Bmat4(3)=NUCBFC(jp,3)

C  Get primitive Electron Basis Function Contraction Coefficients 
              Cof_ie1=AGEBFCC(ie1)
              Cof_ie2=AGEBFCC_be(ie2)
              Cof_ie3=AGEBFCC_be(ie3)
              Cof_je1=AGEBFCC(je1)
              Cof_je2=AGEBFCC_be(je2)
              Cof_je3=AGEBFCC_be(je3)
C  Get Nuclear Basis Function Contraction Coefficients
              Cof_ip=AGNBFCC(ip)
              Cof_jp=AGNBFCC(jp)

C---------------------OMG_123------------------------------------------(
              call RXCHFmult_INT_GAM3_MD(I1,J1,K1,A1,Amat1,
     x                                   I2,J2,K2,A2,Amat2,
     x                                   I3,J3,K3,A3,Amat3,
     x                                   I4,J4,K4,A4,Amat4,
     x                                   L1,M1,N1,B1,Bmat1,
     x                                   L2,M2,N2,B2,Bmat2,
     x                                   L3,M3,N3,B3,Bmat3,
     x                                   L4,M4,N4,B4,Bmat4,
     x                                   nat,ngtg1,
     x                                   pmass,cat,zan,
     x                                   bcoef1,gamma1,
     x                                   ans)

!                       call underflow(ans)

                        OMG3=OMG3+ans
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2
     x                      *Cof_ie3*Cof_je3

C---------------------OMG_123------------------------------------------)

          end do
          end do
        end do
        end do
      end do
      end do


      return
      end
C======================================================================
      subroutine RXCHFmult_INT_GAM3_MD(I1,J1,K1,A1,Amat1,
     x                                 I2,J2,K2,A2,Amat2,
     x                                 I3,J3,K3,A3,Amat3,
     x                                 I4,J4,K4,A4,Amat4,
     x                                 L1,M1,N1,B1,Bmat1,
     x                                 L2,M2,N2,B2,Bmat2,
     x                                 L3,M3,N3,B3,Bmat3,
     x                                 L4,M4,N4,B4,Bmat4,
     x                                 nat,ngtg1,
     x                                 pmass,cat,zan,
     x                                 bcoef1,gamma1,
     x                                 ans)

C Adapted ../gam_3_OMP.f to account for interaction terms
C======================================================================
      implicit none

C Input Variables
      integer nat
      integer ngtg1
      double precision pmass
      double precision zan(nat)
      double precision cat(3,nat)
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)

      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1,Amat1(3) 
      double precision A2,Amat2(3) 
      double precision A3,Amat3(3) 
      double precision A4,Amat4(3) 
      double precision B1,Bmat1(3) 
      double precision B2,Bmat2(3) 
      double precision B3,Bmat3(3) 
      double precision B4,Bmat4(3) 

C Variables Returned
      double precision ans   ! INT OMG3 contribution

C Local Variables
      integer iii    ! Index for looping over natoms
      integer ik,il  ! Indices for geminal loops
      integer iat

      double precision gamA14
      double precision gamA24
      double precision gamA34
      double precision gamB14
      double precision gamB24
      double precision gamB34
      double precision gamA
      double precision gamB

      double precision cmat(3)
      double precision znuc

      double precision gVEE
      double precision xgVEE
      double precision xx,yy,zz
      double precision zero,half,one,two,four
      parameter(zero=0.0d+00,one=1.0d+00,two=2.0d+00,four=4.0d+00)
      parameter(half=0.5d+00)

!      double precision gHEg
      double precision gVEEg1
      double precision gVEEg2
      double precision gVEPg1  ! INT_GAM3 only
!      double precision xgHEg
      double precision xgVEEg1 
      double precision xgVEEg2 
      double precision xgVEPg1 ! INT_GAM3 only

      double precision xmass,coulomb_sign
      double precision xke,Vc
      double precision val_vec 
    
      integer NQUAD_coul
      integer NQUAD_ovlap


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     BASIS FUNCTIONS: ASSIGN CENTERS, EXPONENTS, ANGULAR MOM.
C

C     *****ie1 ::  electron 1 bra *****
C     *****ie2 ::  electron 2 bra *****      
C     *****ie3 ::  electron 3 bra *****      
C     *****ip  ::  proton bra     *****

C     *****je1 ::  electron 1 ket *****
C     *****je2 ::  electron 2 ket *****
C     *****je3 ::  electron 3 ket *****
C     *****jp  ::  proton ket     *****

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C      gVEE=zero
C
C      DO IK=1,NGTG1
C
C         gamA=gamma1(ik)
C         gamB=zero
C
CC  xgVEE
CC  --- g(1,p)VEE(2,3)---
C         call G2_MD_xggs(I1,J1,K1,A1,Amat1,
C     1                   I4,J4,K4,A4,Amat4,
C     2                   L1,M1,N1,B1,Bmat1,
C     3                   L4,M4,N4,B4,Bmat4,
C     4                   gamA,gamB,xx)
C
Cc        call pgiovlap(I1,J1,K1,A1,Amat1,
Cc    x                 I4,J4,K4,A4,Amat4,
Cc    x                 L1,M1,N1,B1,Bmat1,
Cc    x                 L4,M4,N4,B4,Bmat4,
Cc    x                 gamA,gamB,xx)
C
C
C         call gfvee(I2,J2,K2,A2,Amat2,
C     x              I3,J3,K3,A3,Amat3,
C     x              L2,M2,N2,B2,Bmat2,
C     x              L3,M3,N3,B3,Bmat3,
C     x              yy)
C
C         call underflow(xx)
C         call underflow(yy)
C
C         xgVEE=xx*yy
C         gVEE=gVEE+(bcoef1(ik)*xgVEE)
C
CC End 1 gamma loop
C      end do
C Begin 2 gamma loop

!      gHEg=zero
      gVEEg1=zero
      gVEEg2=zero
      gVEPg1=zero ! INT_GAM3 only

      DO IK=1,NGTG1
         DO IL=1,NGTG1

C>>>>>>>>>>>>>>>>>>>>  xgHEg <<<<<<<<<<<<<<<<<<<<
c           ndim=3
c           natom=nat
c           xmass=one
c           coulomb_sign=-one
c           call o3_Hcore_val(NDIM,NATOM,
c    x                        xmass,zan,cmat,
c    x                        coulomb_sign,
c    x                        I3,J3,K3,A3,Amat3,
c    x                        L3,M3,N3,B3,Bmat3,
c    x                        xx)
CCCCCC-rather than call o3_Hcore_val evaluate 
C everything right here...

!            xmass=one
!            coulomb_sign=-one

!            call gfke(I3,J3,K3,A3,Amat3,
!     x                L3,M3,N3,B3,Bmat3,
!     x                xmass,xke)
!            Vc = ZERO
!            do iat=1,nat
!                 Cmat(1)=cat(1,iat)
!                 Cmat(2)=cat(2,iat)
!                 Cmat(3)=cat(3,iat)
!                 call gfvec(I3,J3,K3,A3,Amat3,
!     x                      L3,M3,N3,B3,Bmat3,
!     x                      Cmat,val_vec)
!                 Vc = Vc + (zan(iat)*val_vec)
!             end do
c            hcore = xke + (Vc*coulomb_sign)
!             xx = xke + (Vc*coulomb_sign)
CCCCCC

c           gamA14 = gamma1(ik)
c           gamA24 = ZERO
c           gamA34 = ZERO
c           gamB14 = ZERO
c           gamB24 = gamma1(il)
c           gamB34 = ZERO

!            gamA = gamma1(ik)
!            gamB = gamma1(il)

c           call G4ovlap_typ1(I1,J1,K1,A1,Amat1,
c           call G3ovlap(I1,J1,K1,A1,Amat1,
c    x                   I2,J2,K2,A2,Amat2,
c    x                   I4,J4,K4,A4,Amat4,
c    x                   L1,M1,N1,B1,Bmat1,
c    x                   L2,M2,N2,B2,Bmat2,
c    x                   L4,M4,N4,B4,Bmat4,
c    x                   ZERO,gamA,ZERO,
c    x                   ZERO,ZERO,gamB,yy)
cc   4                gamA12,gamA13,gamA23,
cc   4                gamB12,gamB13,gamB23,sval)
!            call G3_MD_xggs(I1,J1,K1,A1,Amat1,
!     x                      I2,J2,K2,A2,Amat2,
!     x                      I4,J4,K4,A4,Amat4,
!     x                      L1,M1,N1,B1,Bmat1,
!     x                      L2,M2,N2,B2,Bmat2,
!     x                      L4,M4,N4,B4,Bmat4,
!     x                      ZERO,gamA,ZERO,
!     x                      ZERO,ZERO,gamB,yy)

C RXCHFmult(  reorder indicies for INT_GAM3 - should not affect XCHF_GAM3
C    index 1: regular electron
C    index 2: special electron 1
C    index 3: special electron 2
C    index 4: proton
C )
C           --g(e2,p1) V^{ep}(e1,p1) g(e3,p1)--
            coulomb_sign  = -ONE
            gamA14 = gamma1(ik)
            gamA24 = ZERO
            gamA34 = ZERO
            gamB14 = ZERO
            gamB24 = gamma1(il)
            gamB34 = ZERO
!           call rys_G4vee_r34(NQUAD_coul,NQUAD_ovlap,
c           call interface_G4vee_r34(NQUAD_coul,NQUAD_ovlap,
c    x                               I1,J1,K1,A1,Amat1,
c    x                               I2,J2,K2,A2,Amat2,
c    x                               I3,J3,K3,A3,Amat3,
c    x                               I4,J4,K4,A4,Amat4,
c    x                               L1,M1,N1,B1,Bmat1,
c    x                               L2,M2,N2,B2,Bmat2,
c    x                               L3,M3,N3,B3,Bmat3,
c    x                               L4,M4,N4,B4,Bmat4,
c    x                               gamA14,gamA24,gamA34,
c    x                               gamB14,gamB24,gamB34,zz)
            call G4_MD_xgVepg(I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        I1,J1,K1,A1,Amat1,
     *                        I4,J4,K4,A4,Amat4,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gamA14,gamB24,
     *                        zz)
c           call G4Vep_AUX_g14g24V34(I1,J1,K1,A1,Amat1,
c    *                               I2,J2,K2,A2,Amat2,
c    *                               I3,J3,K3,A3,Amat3,
c    *                               I4,J4,K4,A4,Amat4,
c    *                               L1,M1,N1,B1,Bmat1,
c    *                               L2,M2,N2,B2,Bmat2,
c    *                               L3,M3,N3,B3,Bmat3,
c    *                               L4,M4,N4,B4,Bmat4,
c    *                               gamA14,zero,
c    *                               zero,gamB24,
c    *                               zz)

!             call underflow(xx)
!             call underflow(yy)
             call underflow(zz)

!             xgHEg = ((xx*yy) + (zz*coulomb_sign) )
           xgVEPg1 = zz*coulomb_sign


C>>>>>>>>>>>>>>>>>>>>  xgVEEg1 <<<<<<<<<<<<<<<<<<<<
C RXCHFmult( reorder indicies for INT_GAM3 - should not affect XCHF_GAM3
C    index 1: regular electron
C    index 2: special electron 1
C    index 3: special electron 2
C    index 4: proton
C )
C  --- g(2,p)VEE(1,3)g(2,p)---
            gamA=gamma1(ik)
            gamB=gamma1(il)
            call G2_MD_xggs(I2,J2,K2,A2,Amat2,
     1                      I4,J4,K4,A4,Amat4,
     2                      L2,M2,N2,B2,Bmat2,
     3                      L4,M4,N4,B4,Bmat4,
     4                      gamA,gamB,xx)

c           call pgiovlap(I1,J1,K1,A1,Amat1,
c    x                    I4,J4,K4,A4,Amat4,
c    x                    L1,M1,N1,B1,Bmat1,
c    x                    L4,M4,N4,B4,Bmat4,
c    x                    gamA,gamB,xx)


            call gfvee(I1,J1,K1,A1,Amat1,
     x                 I3,J3,K3,A3,Amat3,
     x                 L1,M1,N1,B1,Bmat1,
     x                 L3,M3,N3,B3,Bmat3,
     x                 yy)


            call underflow(xx)
            call underflow(yy)

            xgVEEg1=xx*yy

C>>>>>>>>>>>>>>>>>>>>  xgVEEg2 <<<<<<<<<<<<<<<<<<<<
C ans = <GA(1)GA(2)GA(3)GB(4)|g(1,4)g(2,4)g(3,4)/(r1-r2)|GB(1)GB(2)GB(3)GB(4)>
C ans = <ie1 ie2 ie3 ip|g(1,p)g(2,p)g(3,p)/(r1-r2)|je1 je2 je3 jp>
            gamA14=gamma1(ik)
            gamA24=zero
            gamA34=zero
            gamB14=zero
            gamB24=zero
            gamB34=gamma1(il)
CCWS- 11-08-2010(:  There is a problem with the rys_G4Vee_r12
C                   AC routine, so call a CS routine to 
C                   evaluate the xgVEEg2 integral:
c           NQUAD_coul=5
c           NQUAD_ovlap=5
c           call rys_G4Vee_r12(NQUAD_coul,NQUAD_ovlap,
c    x                              I1,J1,K1,A1,Amat1,
c    x                              I2,J2,K2,A2,Amat2,
c    x                              I3,J3,K3,A3,Amat3,
c    x                              I4,J4,K4,A4,Amat4,
c    x                              L1,M1,N1,B1,Bmat1,
c    x                              L2,M2,N2,B2,Bmat2,
c    x                              L3,M3,N3,B3,Bmat3,
c    x                              L4,M4,N4,B4,Bmat4,
c    x                              gamA14,gamA24,gamA34,
c    x                              gamB14,gamB24,gamB34,
c    x                              xgVEEg2)
C RXCHFmult(  reorder indicies for INT_GAM2 - should not affect XCHF_GAM2
C    index 1: regular electron
C    index 2: special electron 1
C    index 3: special electron 2
C    index 4: proton
C  --- g(2,p)VEE(1,2)g(3,p)---
            call G4_MD_xgVeeg(I2,J2,K2,A2,Amat2,
     *                        I1,J1,K1,A1,Amat1,
     *                        I3,J3,K3,A3,Amat3,
     *                        I4,J4,K4,A4,Amat4,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gamA14,gamB34,
     *                        xgVEEg2)
CCWS- 11-08-2010 )
c           call G4Vee_AUX_g14g34V12(I1,J1,K1,A1,Amat1,
c    *                               I2,J2,K2,A2,Amat2,
c    *                               I3,J3,K3,A3,Amat3,
c    *                               I4,J4,K4,A4,Amat4,
c    *                               L1,M1,N1,B1,Bmat1,
c    *                               L2,M2,N2,B2,Bmat2,
c    *                               L3,M3,N3,B3,Bmat3,
c    *                               L4,M4,N4,B4,Bmat4,
c    *                               gamA14,gamB34,
c    *                               xgVEEg2)


C  Sum terms against bcoeff
!            gHEg  =gHEg  +(bcoef1(ik)*bcoef1(il)*xgHEg)
            gVEEg1=gVEEg1+(bcoef1(ik)*bcoef1(il)*xgVEEg1)
            gVEEg2=gVEEg2+(bcoef1(ik)*bcoef1(il)*xgVEEg2)
            gVEPg1=gVEPg1+(bcoef1(ik)*bcoef1(il)*xgVEPg1) ! INT_GAM3 only

C End 2 gamma loop
         end do
      end do


C Total integral build
!      ans1 = gHEg + gVEEg1*half + four*gVEEg2*half
      ans = gVEPg1 + gVEEg1 + two*gVEEg2
CCWS-debug
c                    write(*,*)'========='
c                    write(*,*)'ie1=',ie1
c                    write(*,*)'je1=',je1
c                    write(*,*)'ie2=',ie2
c                    write(*,*)'je2=',je2
c                    write(*,*)'ie3=',ie3
c                    write(*,*)'je3=',je3
c                    write(*,*)'gHEg   =',gHEg
c                    write(*,*)'gVEEg1 =',gVEEg1
c                    write(*,*)'gVEEg2 =',gVEEg2
c                    write(*,*)'ans    =',ans
c                    write(*,*)'========='
CCWS-debug



      return
      end

