!=======================================================================
      subroutine FAE_OMG3(Nchunks,nebf,npbf,ng3,NAE,NBE,
     x                    DAE,DBE,DP,GM3ICR,FAE,E_AE_OMG3)

! Calculate Alpha and Beta Fock matrices
! and contribution to total energy for 4-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng3,nebf,npbf
      integer NAE,NBE
      double precision GM3ICR(ng3)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FAE(nebf,nebf)
!     double precision FBE(nebf,nebf)
!     double precision FP(npbf,npbf)
      double precision E_AE_OMG3

! Local Variables
      logical LAB,LBB,LAA
      integer istat,ichunk,istart,iend,ng3_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
      double precision XFAE(nebf,nebf)
!     double precision XFBE(nebf,nebf)
!     double precision XFP(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_AE_OMG3=0.0d+00
      XFAE=0.0d+00
!     XFBE=0.0d+00
!     XFP =0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

      LAB=.FALSE.
      LBB=.FALSE.
      LAA=.FALSE.

! Determine if there will be mixed Alpha-Beta contributions to the 
! Alpha Fock matrices:
!     LAB = (NBE.ge.1)
      LAB = ( NBE.ge.1 .and. NAE.ge.2 )
! Determine if there will be a Beta-Beta contribution to the 
! Alpha Fock matrices
      LBB = (NBE.ge.2)
! Determine if there will be a Alpha-Alpha contribution to the 
! Alpha Fock matrices
!     LAA = (NAE.ge.2)
      LAA = (NAE.ge.3)

!     LAB=.TRUE.
!     LBB=.TRUE.
!     LAA=.TRUE.

!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      do ichunk=1,Nchunks

         call loop_size(1,ng3,Nchunks,ichunk-1,istart,iend)

! Segment of ng3:
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
                         end if

                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

         call thread_FAEOMG3(LAB,LBB,LAA,NAE,istart,iend,ng3_seg,ng3,
     x                       nebf,npbf,loop_map,DAE,DBE,DP,GM3ICR,
     x                       XFAE,E_AE_OMG3)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
      call add2fock(nebf,XFAE,FAE)
!     call add2fock(nebf,XFBE,FBE)
!     call add2fock(npbf,XFP,FP)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!======================================================================
      subroutine thread_FAEOMG3(LAB,LBB,LAA,NAE,istart,iend,ng3_seg,ng3,
     x                       nebf,npbf,loop_map,DAE,DBE,DP,GM3ICR,
     x                       XFAE,E_AE_OMG3)
!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      logical LAB,LBB,LAA
      integer istart,iend,ng3_seg
      integer NAE
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng3

      integer loop_map(ng3_seg,8)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM3ICR(ng3)

! Variables Returned
      double precision XFAE(nebf,nebf)
!     double precision XFBE(nebf,nebf)
!     double precision XFP(npbf,npbf)
      double precision E_AE_OMG3

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer iec3,jec3  !
      integer imap,ia
      integer ia_123
      integer ia_213
      integer ia_312
      integer ia_132
      integer ia_231
      integer ia_321

      double precision OMG6,OMG23,OMG12,OMG13 
      double precision XAA,XBB,XAB,XBA 

      double precision EFAC

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

!     EFAC=1.0d+00
      EFAC=3.0d+00
!     if(NAE.eq.1) EFAC=1.0d+00

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(LAB,LBB,LAA)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(EFAC)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM3ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia_123,ia_213,ia_312,ia_132,ia_231,ia_321)
!$ompx private(OMG6,OMG23,OMG12,OMG13)
!$ompx private(XAA,XBB,XAB,XBA)
!$ompx reduction(+:E_AE_OMG3)
!$ompx reduction(+:XFAE)

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

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia_123)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec2,iec2,jec1,iec3,jec3,ia_213)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec3,iec2,jec1,iec3,jec2,ia_312)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec3,iec3,jec2,ia_132)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec2,iec2,jec3,iec3,jec1,ia_231)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec3,iec2,jec2,iec3,jec1,ia_321)

         OMG6=0.0d+00
         OMG23=0.0d+00
         OMG12=0.0d+00
         OMG13=0.0d+00

         if(LAA) OMG6=
     x       ( GM3ICR(ia_123)     
     x        -GM3ICR(ia_213)
     x        +GM3ICR(ia_312)
     x        -GM3ICR(ia_132)
     x        +GM3ICR(ia_231)
     x        -GM3ICR(ia_321) )

!        P23 OMG3_132
!        P12 OMG3_213
!        P13 OMG3_321
         if(LBB) OMG23=( GM3ICR(ia_123) - GM3ICR(ia_132) )
         if(LAB) OMG12=( GM3ICR(ia_123) - GM3ICR(ia_213) )
         if(LAB) OMG13=( GM3ICR(ia_123) - GM3ICR(ia_321) )

!----------------Form-Alpha-Fock-Matrix-------------------(

         XAA=DAE(iec2,jec2)*DAE(iec3,jec3)*OMG6
         XBB=DBE(iec2,jec2)*DBE(iec3,jec3)*OMG23
         XAB=DAE(iec2,jec2)*DBE(iec3,jec3)*OMG12
         XBA=DBE(iec2,jec2)*DAE(iec3,jec3)*OMG13

!        XFAE(iec1,jec1)=XFAE(iec1,jec1)+
!    x         3.0d+00*DP(ip,jp)*(XAA+XBB+XAB+XBA)

         XFAE(iec1,jec1)=XFAE(iec1,jec1)+
     x         EFAC*DP(ip,jp)*(XAA+XBB+XAB+XBA)

!----------------Form-Alpha-Fock-Matrix-------------------)

! Energy
         E_AE_OMG3 = E_AE_OMG3 +
     x         DP(ip,jp)*DAE(iec1,jec1)*(XAA+XBB+XAB+XBA)
         

         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!=======================================================================
      subroutine UFP_OMG3(Nchunks,nebf,npbf,ng3,NAE,NBE,
     x                    DAE,DBE,DP,GM3ICR,FP,E_P_OMG3)

! Calculate Alpha and Beta contributions
! to QM particle Fock matrix for 4-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng3,nebf,npbf
      integer NAE,NBE
      double precision GM3ICR(ng3)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
!     double precision FAE(nebf,nebf)
!     double precision FBE(nebf,nebf)
      double precision FP(npbf,npbf)
      double precision E_P_OMG3

! Local Variables
      logical LAAA,LBBB,LAAB,LBBA
      integer istat,ichunk,istart,iend,ng3_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
!     double precision XFAE(nebf,nebf)
!     double precision XFBE(nebf,nebf)
      double precision XFP(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_P_OMG3=0.0d+00
      XFP =0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

      LAAA=.FALSE.
      LBBB=.FALSE.
      LAAB=.FALSE.
      LBBA=.FALSE.

! ? Alpha-Alpha-Alpha contribution to the Alpha Fock matrices
      LAAA = (NAE.ge.3)
! ? Beta-Beta-Beta contribution to the Alpha Fock matrices
      LBBB = (NBE.ge.3)
! ? Alpha-Alpha-Beta contribution
      LAAB = ( (NAE.ge.2) .and. (NBE.ge.1) )
! ? Beta-Beta-Alpha contribution
      LBBA = ( (NBE.ge.2) .and. (NAE.ge.1) )

!     LAAA=.true.
!     LBBB=.true.
!     LAAB=.true.
!     LBBA=.true.

!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      do ichunk=1,Nchunks

         call loop_size(1,ng3,Nchunks,ichunk-1,istart,iend)

! Segment of ng3:
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
                         end if

                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

         call thread_UFPOMG3(LAAA,LBBB,LAAB,LBBA,
     x                       istart,iend,ng3_seg,ng3,nebf,npbf,
     x                       loop_map,DAE,DBE,DP,GM3ICR,XFP,E_P_OMG3)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
!     call add2fock(nebf,XFAE,FAE)
!     call add2fock(nebf,XFBE,FBE)
      call add2fock(npbf,XFP,FP)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!=======================================================================
      subroutine thread_UFPOMG3(LAAA,LBBB,LAAB,LBBA,
     x                          istart,iend,ng3_seg,ng3,nebf,npbf,
     x                          loop_map,DAE,DBE,DP,GM3ICR,XFP,E_P_OMG3)
!=======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      logical LAAA,LBBB,LAAB,LBBA
      integer istart,iend,ng3_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng3

      integer loop_map(ng3_seg,8)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM3ICR(ng3)

! Variables Returned
!     double precision XFAE(nebf,nebf)
!     double precision XFBE(nebf,nebf)
      double precision XFP(npbf,npbf)
      double precision E_P_OMG3

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer iec3,jec3  !
      integer imap,ia
      integer ia_123
      integer ia_213
      integer ia_312
      integer ia_132
      integer ia_231
      integer ia_321

      double precision OMG6,OMG23,OMG12,OMG13 
!     double precision XAA,XBB,XAB,XBA 
!     double precision XAAA_BBB,XAAB_BBA,XABA_BAB,XABB_BAA 
      double precision CAAA,CBBB,CAAB,CBBA,CABA,CBAB,CABB,CBAA

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(LAAA,LBBB,LAAB,LBBA)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM3ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia_123,ia_213,ia_312,ia_132,ia_231,ia_321)
!$ompx private(OMG6,OMG23,OMG12,OMG13)
!$ompx private(CAAA,CBBB,CAAB,CBBA,CABA,CBAB,CABB,CBAA)
!$ompx reduction(+:XFP)
!$ompx reduction(+:E_P_OMG3)

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

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia_123)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec2,iec2,jec1,iec3,jec3,ia_213)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec3,iec2,jec1,iec3,jec2,ia_312)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec3,iec3,jec2,ia_132)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec2,iec2,jec3,iec3,jec1,ia_231)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec3,iec2,jec2,iec3,jec1,ia_321)

         OMG6=
     x       ( GM3ICR(ia_123)     
     x        -GM3ICR(ia_213)
     x        +GM3ICR(ia_312)
     x        -GM3ICR(ia_132)
     x        +GM3ICR(ia_231)
     x        -GM3ICR(ia_321) )

!        P23 OMG3_132
!        P12 OMG3_213
!        P13 OMG3_321
         OMG23=( GM3ICR(ia_123) - GM3ICR(ia_132) )
         OMG12=( GM3ICR(ia_123) - GM3ICR(ia_213) )
         OMG13=( GM3ICR(ia_123) - GM3ICR(ia_321) )

!----------------Form-QM-Partice-Fock-Matrix-------------------------(
CAAA,CBBB,CAAB,CBBA,CABA,CBAB,CABB,CBAA
      CAAA=0.0d+00
      CBBB=0.0d+00 
      CAAB=0.0d+00 
      CBBA=0.0d+00 
      CABA=0.0d+00 
      CBAB=0.0d+00 
      CABB=0.0d+00 
      CBAA=0.0d+00 
 
      if(LAAA) CAAA=DAE(iec1,jec1)*DAE(iec2,jec2)*DAE(iec3,jec3)

      if(LBBB) CBBB=DBE(iec1,jec1)*DBE(iec2,jec2)*DBE(iec3,jec3)
       
      if(LAAB) then 
         CAAB=DAE(iec1,jec1)*DAE(iec2,jec2)*DBE(iec3,jec3)
         CABA=DAE(iec1,jec1)*DBE(iec2,jec2)*DAE(iec3,jec3)
         CBAA=DBE(iec1,jec1)*DAE(iec2,jec2)*DAE(iec3,jec3)
      end if

      if(LBBA) then
         CBBA=DBE(iec1,jec1)*DBE(iec2,jec2)*DAE(iec3,jec3)
         CBAB=DBE(iec1,jec1)*DAE(iec2,jec2)*DBE(iec3,jec3)
         CABB=DAE(iec1,jec1)*DBE(iec2,jec2)*DBE(iec3,jec3)
      end if

      XFP(ip,jp)=XFP(ip,jp)
     x          +(CAAA+CBBB)*OMG6
     x          +(CAAB+CBBA)*OMG12
     x          +(CABA+CBAB)*OMG13
     x          +(CABB+CBAA)*OMG23

      E_P_OMG3 = E_P_OMG3 + DP(ip,jp)*(
     x          +(CAAA+CBBB)*OMG6
     x          +(CAAB+CBBA)*OMG12
     x          +(CABA+CBAB)*OMG13
     x          +(CABB+CBAA)*OMG23 )


!      XAAA_BBB=(DAE(iec1,jec1)*DAE(iec2,jec2)*DAE(iec3,jec3))
!    x         +(DBE(iec1,jec1)*DBE(iec2,jec2)*DBE(iec3,jec3))*OMG6

!      XAAB_BBA=(DAE(iec1,jec1)*DAE(iec2,jec2)*DBE(iec3,jec3))
!    x         +(DBE(iec1,jec1)*DBE(iec2,jec2)*DAE(iec3,jec3))*OMG12

!      XABA_BAB=(DAE(iec1,jec1)*DBE(iec2,jec2)*DAE(iec3,jec3))
!    x         +(DBE(iec1,jec1)*DAE(iec2,jec2)*DBE(iec3,jec3))*OMG13

!      XABB_BAA=(DAE(iec1,jec1)*DBE(iec2,jec2)*DBE(iec3,jec3))
!    x         +(DBE(iec1,jec1)*DAE(iec2,jec2)*DAE(iec3,jec3))*OMG23

!      XFP(ip,jp)=XFP(ip,jp)+(XAAA_BBB+XAAB_BBA+XABA_BAB+XABB_BAA)
!----------------Form-QM-Partice-Fock-Matrix-------------------------)


         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end

