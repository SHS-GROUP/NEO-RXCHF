!=======================================================================
      subroutine GAM3_IC2(Nchunks,nebf,npebf,npbf,
     x                    ng3,ng3prm,nat,ngtg1,
     x                    pmass,cat,zan,bcoef1,gamma1,
     x                    KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                    ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                    GM3ICR)

!=======================================================================
      implicit none
      include 'omp_lib.h'
! Input Variables
      integer Nchunks
      integer ng3,nebf,npebf,npbf,ng3prm
      integer nat,ngtg1
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

! Variables Returned
      double precision GM3ICR(ng3)

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer ie1,je1,ie2,je2,ie3,je3
      integer,allocatable :: loop_map(:,:)
      double precision,allocatable :: XG3ICR(:)
!     double precision,allocatable :: GAM_3(:)
!     double precision,allocatable :: TMPARY(:)

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


         write(*,1000) ng3prm,nchunks,omp_get_num_procs(),
     xomp_get_max_threads(),1
         wtime = omp_get_wtime()

! Memory to hold un-symmetrized GAM3 integrals
         if(allocated(XG3ICR)) deallocate(XG3ICR)
         allocate( XG3ICR(ng3),stat=istat )

!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      do ichunk=1,Nchunks

         wtime2 = omp_get_wtime()

         call loop_size(1,ng3prm,Nchunks,ichunk-1,istart,iend)
         ng3_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )

! Nested loop compression for this chunk:
         Loopi=0
         imas=0
         do ip=1,npbf
         do jp=1,npbf
            do ie1=1,npebf
            do je1=1,npebf
               do ie2=1,npebf
               do je2=1,npebf
                  do ie3=1,npebf
                  do je3=1,npebf

                     imas=imas+1 ! imas is master_index
                     if(imas.ge.istart.and.imas.le.iend) then
                        Loopi=Loopi+1
                        loop_map(Loopi,1)=je3
                        loop_map(Loopi,2)=ie3
                        loop_map(Loopi,3)=je2
                        loop_map(Loopi,4)=ie2
                        loop_map(Loopi,5)=je1
                        loop_map(Loopi,6)=ie1
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

         call thread_gam3_IC2(istart,iend,ng3_seg,ng3,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        loop_map,XG3ICR,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)


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
            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

C  GAM_2 Symmetrization:
C  Determine packing indices for XGAM_3 integral matrices

C              As Packed-->       XGAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
!                            read(905,REC=ia_123) x123
                             ! x123=GAM_3(ia_123)
                               x123=XG3ICR(ia_123)

C              As Packed-->       XGAM_3(je2,ie2,je3,ie3,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec3,jec3,iec2,jec2,ia)
                             ia_132=ia
!                            read(905,REC=ia_132) x132
                             ! x132=GAM_3(ia_132)
                               x132=XG3ICR(ia_132)

C              As Packed-->       XGAM_3(je3,ie3,je1,ie1,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,iec3,jec3,ia)
                             ia_213=ia
!                            read(905,REC=ia_213) x213
                             ! x213=GAM_3(ia_213)
                               x213=XG3ICR(ia_213)

C              As Packed-->       XGAM_3(je1,ie1,je3,ie3,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec3,jec3,iec1,jec1,ia)
                             ia_231=ia
!                            read(905,REC=ia_231) x231
                             ! x231=GAM_3(ia_231)
                               x231=XG3ICR(ia_231)

C              As Packed-->       XGAM_3(je2,ie2,je1,ie1,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec3,jec3,iec1,jec1,iec2,jec2,ia)
                             ia_312=ia
!                            read(905,REC=ia_312) x312
                             ! x312=GAM_3(ia_312)
                               x312=XG3ICR(ia_312)

C              As Packed-->       XGAM_3(je1,ie1,je2,ie2,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec3,jec3,iec2,jec2,iec1,jec1,ia)
                             ia_321=ia
!                            read(905,REC=ia_321) x321
                             ! x321=GAM_3(ia_321)
                               x321=XG3ICR(ia_321)

                       xxxx=(x123+x132+x213+x231+x312+x321)/six


c                      write(*,*)'xxxx=',xxxx
cc                     GAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3)=xxxx 
c                      GAM_3PK(ia_123)=xxxx 
c                      call put_GAM3(ia_123,ng3,xxxx)
!                      write(805,REC=ia_123) xxxx
                       GM3ICR(ia_123)=xxxx 

cc                     GAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2)=xxxx
c                      GAM_3PK(ia_132)=xxxx 
c                      call put_GAM3(ia_132,ng3,xxxx)
!                      write(805,REC=ia_132) xxxx
                       GM3ICR(ia_132)=xxxx 

cc                     GAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3)=xxxx
c                      GAM_3PK(ia_213)=xxxx 
c                      call put_GAM3(ia_213,ng3,xxxx)
!                      write(805,REC=ia_213) xxxx
                       GM3ICR(ia_213)=xxxx 

cc                     GAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1)=xxxx
c                      GAM_3PK(ia_231)=xxxx 
c                      call put_GAM3(ia_231,ng3,xxxx)
!                      write(805,REC=ia_231) xxxx
                       GM3ICR(ia_231)=xxxx 

cc                     GAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2)=xxxx
c                      GAM_3PK(ia_312)=xxxx 
c                      call put_GAM3(ia_312,ng3,xxxx)
!                      write(805,REC=ia_312) xxxx
                       GM3ICR(ia_312)=xxxx 

cc                     GAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1)=xxxx
c                      GAM_3PK(ia_321)=xxxx 
c                      call put_GAM3(ia_321,ng3,xxxx)
!                      write(805,REC=ia_321) xxxx
                       GM3ICR(ia_321)=xxxx 


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
     x        6X,'|            --IN-CORE APPROACH2--            |',/,
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
      subroutine thread_gam3_IC2(istart,iend,ng3_seg,ng3,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           loop_map,XG3ICR,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

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
      integer ng3

      integer ip,jp
      integer ie1,je1  !
      integer ie2,je2  ! Primitive elec basis function indices
      integer ie3,je3  !
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
!     double precision GAM_3(ng3_seg)

! Variables Returned
      double precision XG3ICR(ng3)

! Local Variables
      integer imap,ia
      double precision ans
!--------------------------------(
! Basis set-related local variables
      integer iec1,jec1,iec2,jec2,iec3,jec3
      double precision Cof_ie1,Cof_je1
      double precision Cof_ie2,Cof_je2
      double precision Cof_ie3,Cof_je3
      double precision Cof_ip,Cof_jp

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
!--------------------------------)

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(istart,iend)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1)
!$ompx private(ie2,je2)
!$ompx private(ie3,je3)
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia) 
!$ompx private(A1,I1,J1,K1,Amat1)
!$ompx private(A2,I2,J2,K2,Amat2)
!$ompx private(A3,I3,J3,K3,Amat3)
!$ompx private(A4,I4,J4,K4,Amat4)
!$ompx private(B1,L1,M1,N1,Bmat1)
!$ompx private(B2,L2,M2,N2,Bmat2)
!$ompx private(B3,L3,M3,N3,Bmat3)
!$ompx private(B4,L4,M4,N4,Bmat4)
!$ompx private(Cof_ie1,Cof_je1)
!$ompx private(Cof_ie2,Cof_je2)
!$ompx private(Cof_ie3,Cof_je3)
!$ompx private(Cof_ip,Cof_jp)
!$ompx private(ans)
!$ompx private(id)
!$ompx reduction(+:XG3ICR)
!!!!!!!!!!!!!x$ompx shared(XG3ICR)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         je3=loop_map(imap,1)
         ie3=loop_map(imap,2)
         je2=loop_map(imap,3)
         ie2=loop_map(imap,4)
         je1=loop_map(imap,5)
         ie1=loop_map(imap,6)
         jp =loop_map(imap,7)
         ip =loop_map(imap,8)

C Get Basis set info:
         A1=ELCEX(ie1)
         I1=ELCAM(ie1,1)
         J1=ELCAM(ie1,2)
         K1=ELCAM(ie1,3)
         Amat1(1)=ELCBFC(ie1,1)
         Amat1(2)=ELCBFC(ie1,2)
         Amat1(3)=ELCBFC(ie1,3)

         A2=ELCEX(ie2)
         I2=ELCAM(ie2,1)
         J2=ELCAM(ie2,2)
         K2=ELCAM(ie2,3)
         Amat2(1)=ELCBFC(ie2,1)
         Amat2(2)=ELCBFC(ie2,2)
         Amat2(3)=ELCBFC(ie2,3)

         A3=ELCEX(ie3)
         I3=ELCAM(ie3,1)
         J3=ELCAM(ie3,2)
         K3=ELCAM(ie3,3)
         Amat3(1)=ELCBFC(ie3,1)
         Amat3(2)=ELCBFC(ie3,2)
         Amat3(3)=ELCBFC(ie3,3)

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

         B2=ELCEX(je2)
         L2=ELCAM(je2,1)
         M2=ELCAM(je2,2)
         N2=ELCAM(je2,3)
         Bmat2(1)=ELCBFC(je2,1)
         Bmat2(2)=ELCBFC(je2,2)
         Bmat2(3)=ELCBFC(je2,3)

         B3=ELCEX(je3)
         L3=ELCAM(je3,1)
         M3=ELCAM(je3,2)
         N3=ELCAM(je3,3)
         Bmat3(1)=ELCBFC(je3,1)
         Bmat3(2)=ELCBFC(je3,2)
         Bmat3(3)=ELCBFC(je3,3)

         B4=NUCEX(jp)
         L4=NUCAM(jp,1)
         M4=NUCAM(jp,2)
         N4=NUCAM(jp,3)
         Bmat4(1)=NUCBFC(jp,1)
         Bmat4(2)=NUCBFC(jp,2)
         Bmat4(3)=NUCBFC(jp,3)

         call xcalc_GAM3_MD(I1,J1,K1,A1,Amat1,
     x                      I2,J2,K2,A2,Amat2,
     x                      I3,J3,K3,A3,Amat3,
     x                      I4,J4,K4,A4,Amat4,
     x                      L1,M1,N1,B1,Bmat1,
     x                      L2,M2,N2,B2,Bmat2,
     x                      L3,M3,N3,B3,Bmat3,
     x                      L4,M4,N4,B4,Bmat4,
     x                      nat,ngtg1,
     x                      pmass,cat,zan,
     x                      bcoef1,gamma1,
     x                      ans)

!                       call underflow(ans)

!  Map from primitive BF indices to contracted indices
                        iec1=AMPEB2C(ie1)
                        iec2=AMPEB2C(ie2)
                        iec3=AMPEB2C(ie3)
                        jec1=AMPEB2C(je1)
                        jec2=AMPEB2C(je2)
                        jec3=AMPEB2C(je3)
!  Get primitive Electron Basis Function Contraction Coefficients 
                        Cof_ie1=AGEBFCC(ie1)
                        Cof_ie2=AGEBFCC(ie2)
                        Cof_ie3=AGEBFCC(ie3)
                        Cof_je1=AGEBFCC(je1)
                        Cof_je2=AGEBFCC(je2)
                        Cof_je3=AGEBFCC(je3)
!  Get Nuclear Basis Function Contraction Coefficients
                        Cof_ip=AGNBFCC(ip)
                        Cof_jp=AGNBFCC(jp)
!  Map the 8-index contracted integral to 1-D:
                        call index_GAM_3PK(nebf,npbf,
     x                                     ip,jp,
     x                                     iec1,jec1,
     x                                     iec2,jec2,
     x                                     iec3,jec3,ia)

                        XG3ICR(ia)=XG3ICR(ia)+ans
     x                            *Cof_ip*Cof_jp
     x                            *Cof_ie1*Cof_je1
     x                            *Cof_ie2*Cof_je2
     x                            *Cof_ie3*Cof_je3


!        call index_GAM_3PK(nebf,npbf,
!    x          ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)

      end do
!$omp end do
!$omp end parallel      
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
