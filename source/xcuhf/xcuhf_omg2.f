!=======================================================================
      subroutine FAE_OMG2(Nchunks,NAE,NBE,nebf,npbf,ng2,
     x                    DAE,DBE,DETOT,DP,GM2ICR,GM2SICR,
     x                    FAE,SAE2,E_AE_OMG2,S_AE_OMG2)

! Calculate Alpha or Beta Fock matrices
! and contribution to total energy for 3-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng2,nebf,npbf
      integer NAE
      integer NBE
      double precision GM2ICR(ng2)
      double precision GM2SICR(ng2)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DETOT(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FAE(nebf,nebf)
      double precision SAE2(nebf,nebf)
!     double precision FBE(nebf,nebf)
!     double precision SBE(nebf,nebf)
!     double precision FP(npbf,npbf)
!     double precision SP(npbf,npbf)
      double precision E_AE_OMG2
      double precision S_AE_OMG2

! Local Variables
      logical LAA,LAB
      integer istat,ichunk,istart,iend,ng2_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2
      integer,allocatable :: loop_map(:,:)
      double precision XFAE(nebf,nebf)
!     double precision XSAE(nebf,nebf)
!     double precision XFBE(nebf,nebf)
!     double precision XSBE(nebf,nebf)
!     double precision XFP(npbf,npbf)
!     double precision XSP(npbf,npbf)

!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_AE_OMG2=0.0d+00
      S_AE_OMG2=0.0d+00
      XFAE=0.0d+00
      SAE2=0.0d+00
!     XFBE=0.0d+00
!     XSBE=0.0d+00
!     XFP =0.0d+00
!     XSP =0.0d+00

! HACK to form total electronic density:
!     DETOT=0.0d+00
!     call add2fock(nebf,DAE,DETOT)
!     call add2fock(nebf,DBE,DETOT)
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

      LAA=.FALSE.
      LAB=.FALSE.
      if( NAE.ge.2 ) LAA=.TRUE.
      if( NAE.ge.1 .and. NBE.ge.1 ) LAB=.TRUE.
!     LAA=.TRUE.
!     LAB=.TRUE.

!      write(*,*)'OMG2: NAE=  LAA=',NAE,LAA
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      do ichunk=1,Nchunks

         call loop_size(1,ng2,Nchunks,ichunk-1,istart,iend)

! Segment of ng3:
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
                  end if

               end do
               end do
            end do
            end do
         end do
         end do

!        call thread_FAEOMG2(istart,iend,NAE,ng2_seg,ng2,nebf,npbf,
!    x                      loop_map,DAE,DETOT,DP,GM2ICR,GM2SICR,
!    x                      XFAE,SAE2,E_AE_OMG2,S_AE_OMG2)

         call thread_FAEOMG2(LAA,LAB,istart,iend,NAE,
     x                      ng2_seg,ng2,nebf,npbf,
     x                      loop_map,DAE,DBE,DETOT,DP,GM2ICR,GM2SICR,
     x                      XFAE,SAE2,E_AE_OMG2,S_AE_OMG2)
      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
      call add2fock(nebf,XFAE,FAE)
!     call add2fock(nebf,XSAE,SAE)
!     call add2fock(nebf,XFBE,FBE)
!     call add2fock(nebf,XSBE,SBE)
!     call add2fock(npbf,XFP,FP)
!     call add2fock(npbf,XSP,SP)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!======================================================================
      subroutine thread_FAEOMG2(LAA,LAB,istart,iend,NAE,
     x                      ng2_seg,ng2,nebf,npbf,
     x                      loop_map,DAE,DBE,DETOT,DP,GM2ICR,GM2SICR,
     x                      XFAE,XSAE,E_AE_OMG2,S_AE_OMG2)
!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      logical LAA,LAB
      integer istart,iend,ng2_seg
      integer NAE
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng2

      integer loop_map(ng2_seg,6)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DETOT(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM2ICR(ng2)
      double precision GM2SICR(ng2)

! Variables Returned
      double precision XFAE(nebf,nebf)
      double precision XSAE(nebf,nebf)
!     double precision XFBE(nebf,nebf)
!     double precision XSBE(nebf,nebf)
!     double precision XFP(npbf,npbf)
!     double precision XSP(npbf,npbf)
      double precision E_AE_OMG2
      double precision S_AE_OMG2

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap,ia
      integer ia_12
      integer ia_21

      double precision OMG12,OMG21 
      double precision SOMG12,SOMG21 
      double precision COEF1,COEF2
      double precision EFAC

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

!     EFAC=1.0d+00
      EFAC=2.0d+00
!     if(NAE.eq.1) EFAC=1.0d+00
!     LAA=.true.

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(LAA,LAB)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx shared(EFAC)
!$ompx shared(DAE,DBE,DETOT,DP)
!$ompx shared(GM2ICR)
!$ompx shared(GM2SICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(ia_12,ia_21)
!$ompx private(OMG12,OMG21)
!$ompx private(COEF1,COEF2)
!$ompx private(SOMG12,SOMG21)
!$ompx reduction(+:XFAE)
!$ompx reduction(+:XSAE)
!$ompx reduction(+:E_AE_OMG2)
!$ompx reduction(+:S_AE_OMG2)

!xomp do SCHEDULE(RUNTIME)
!     do iLp=1,12
!        if(LAA) write(*,*)'LAA is true'
!        if(.not.LAA) write(*,*)'LAA is false'
!     end do
!xomp end do 

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec2=loop_map(imap,1)
         iec2=loop_map(imap,2)
         jec1=loop_map(imap,3)
         iec1=loop_map(imap,4)
         jp =loop_map(imap,5)
         ip =loop_map(imap,6)

         call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec1,iec2,jec2,ia_12)
         call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec2,iec2,jec1,ia_21)

         OMG12=GM2ICR(ia_12)
         OMG21=GM2ICR(ia_21)

         SOMG12=GM2SICR(ia_12)
         SOMG21=GM2SICR(ia_21)

!----------------Form-Alpha-Fock-Matrix-------------------(
!        XFAE(iec1,jec1)=XFAE(iec1,jec1) + 2.0d+00 *
!    x                  ( DP(ip,jp)*DETOT(iec2,jec2)*OMG12
!    x                  - DP(ip,jp)*DAE(iec2,jec2)*OMG21 )
         
!        XSAE(iec1,jec1)=XSAE(iec1,jec1) + 2.0d+00 *
!    x                  ( DP(ip,jp)*DETOT(iec2,jec2)*SOMG12
!    x                  - DP(ip,jp)*DAE(iec2,jec2)*SOMG21 )


         COEF1=0.0d+00
         COEF2=0.0d+00
         if(LAA) COEF1=DAE(iec2,jec2)*DP(ip,jp)
         if(LAB) COEF2=DBE(iec2,jec2)*DP(ip,jp)
!        COEF1=DAE(iec2,jec2)*DP(ip,jp)
!        COEF2=DBE(iec2,jec2)*DP(ip,jp)
         
         XFAE(iec1,jec1)=XFAE(iec1,jec1) + EFAC *
     x                  ( COEF1*(OMG12-OMG21)
     x                  + COEF2*OMG12 )

         XSAE(iec1,jec1)=XSAE(iec1,jec1) + EFAC *
     x                  ( COEF1*(SOMG12-SOMG21)
     x                  + COEF2*SOMG12 )
!----------------Form-Alpha-Fock-Matrix-------------------)

! Energy and Overlap
!        E_AE_OMG2=E_AE_OMG2+
!    x                  + DP(ip,jp)*DAE(iec1,jec1)
!    x     *( DETOT(iec2,jec2)*OMG12 - DAE(iec2,jec2)*OMG21 )
         
!        S_AE_OMG2=S_AE_OMG2+
!    x                  + DP(ip,jp)*DAE(iec1,jec1)
!    x     *( DETOT(iec2,jec2)*SOMG12 - DAE(iec2,jec2)*SOMG21 )

         E_AE_OMG2=E_AE_OMG2 + DAE(iec1,jec1) *
     x                  ( COEF1*(OMG12-OMG21)
     x                  + COEF2*OMG12 )

         S_AE_OMG2=S_AE_OMG2 + DAE(iec1,jec1) *
     x                  ( COEF1*(SOMG12-SOMG21)
     x                  + COEF2*SOMG12 )
 
         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!=======================================================================
      subroutine UFP_OMG2(Nchunks,nebf,npbf,ng2,NAE,NBE,
     x                    DAE,DBE,DETOT,DP,GM2ICR,GM2SICR,FP,SP2,
     x                    E_P_OMG2,S_P_OMG2)

! Calculate QM particle Fock matrices
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng2,nebf,npbf
      integer NAE,NBE
      double precision GM2ICR(ng2)
      double precision GM2SICR(ng2)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision DETOT(nebf,nebf)

! Variables Returned
!     double precision FAE(nebf,nebf)
!     double precision SAE(nebf,nebf)
!     double precision FBE(nebf,nebf)
!     double precision SBE(nebf,nebf)
      double precision FP(npbf,npbf)
      double precision SP2(npbf,npbf)
      double precision E_P_OMG2
      double precision S_P_OMG2
!     double precision E_gam2

! Local Variables
      logical LAB
      logical LAA
      logical LBB
      integer istat,ichunk,istart,iend,ng2_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2
      integer,allocatable :: loop_map(:,:)
!     double precision XFAE(nebf,nebf)
!     double precision XSAE(nebf,nebf)
!     double precision XFBE(nebf,nebf)
!     double precision XSBE(nebf,nebf)
      double precision XFP(npbf,npbf)
!     double precision XSP(npbf,npbf)

!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
!     E_gam2=0.0d+00
!     XFAE=0.0d+00
!     XSAE=0.0d+00
!     XFBE=0.0d+00
!     XSBE=0.0d+00
      XFP =0.0d+00
      SP2 =0.0d+00
      E_P_OMG2 = 0.0d+00
      S_P_OMG2 = 0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

      LAA=.FALSE.
      LBB=.FALSE.
      LAB=.FALSE.

! Determine control for mixed alpha-beta contributions to the 
! QM particle Fock matrices.  If NAE and NBE are both one or more
! then the mixed contributions exist.
      LAB = ( (NAE.ge.1) .and. (NBE.ge.1) )
!!!!
      LAA = ( NAE.ge.2 )
      LBB = ( NBE.ge.2 )
!!!!
!     if( (NAE.gt.0) .and. (NBE.gt.0) ) LAB=.TRUE.
!     LAA=.TRUE.
!     LBB=.TRUE.
!     LAB=.TRUE.

!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      do ichunk=1,Nchunks

         call loop_size(1,ng2,Nchunks,ichunk-1,istart,iend)

! Segment of ng3:
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
                  end if

               end do
               end do
            end do
            end do
         end do
         end do

         call thread_UFPOMG2(LAA,LBB,LAB,istart,iend,
     x                       ng2_seg,ng2,nebf,npbf,
     x                       loop_map,DAE,DBE,DETOT,DP,GM2ICR,GM2SICR,
     x                       XFP,SP2,E_P_OMG2,S_P_OMG2)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
!     call add2fock(nebf,XFAE,FAE)
!     call add2fock(nebf,XSAE,SAE)
!     call add2fock(nebf,XFBE,FBE)
!     call add2fock(nebf,XSBE,SBE)
      call add2fock(npbf,XFP,FP)
!     call add2fock(npbf,XSP,SP)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!======================================================================
!     subroutine thread_UFPOMG2(LAB,istart,iend,ng2_seg,ng2,nebf,npbf,
!    x                        loop_map,DAE,DBE,DETOT,DP,GM2ICR,GM2SICR,
!    x                         XFP,XSP,E_P_OMG2,S_P_OMG2)
      subroutine thread_UFPOMG2(LAA,LBB,LAB,istart,iend,
     x                          ng2_seg,ng2,nebf,npbf,
     x                         loop_map,DAE,DBE,DETOT,DP,GM2ICR,GM2SICR,
     x                          XFP,XSP,E_P_OMG2,S_P_OMG2)

!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      logical LAA
      logical LBB
      logical LAB
      integer istart,iend,ng2_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng2

      integer loop_map(ng2_seg,6)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DETOT(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM2ICR(ng2)
      double precision GM2SICR(ng2)

! Variables Returned
!     double precision XFAE(nebf,nebf)
!     double precision XSAE(nebf,nebf)
!     double precision XFBE(nebf,nebf)
!     double precision XSBE(nebf,nebf)
      double precision XFP(npbf,npbf)
      double precision XSP(npbf,npbf)
      double precision E_P_OMG2
      double precision S_P_OMG2
!     double precision E_gam2

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap,ia
      integer ia_12
      integer ia_21

      double precision OMG12,OMG21 
      double precision SOMG12,SOMG21 
      double precision COEF1,COEF2,COEF3

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(LAA,LBB,LAB)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx shared(DAE,DBE,DETOT,DP)
!$ompx shared(GM2ICR)
!$ompx shared(GM2SICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(ia_12,ia_21)
!$ompx private(OMG12,OMG21)
!$ompx private(SOMG12,SOMG21)
!$ompx private(COEF1,COEF2,COEF3)
!$ompx reduction(+:XFP)
!$ompx reduction(+:XSP)
!$ompx reduction(+:E_P_OMG2)
!$ompx reduction(+:S_P_OMG2)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec2=loop_map(imap,1)
         iec2=loop_map(imap,2)
         jec1=loop_map(imap,3)
         iec1=loop_map(imap,4)
         jp =loop_map(imap,5)
         ip =loop_map(imap,6)

         call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec1,iec2,jec2,ia_12)
         call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec2,iec2,jec1,ia_21)

         OMG12=GM2ICR(ia_12)
         OMG21=GM2ICR(ia_21)

         SOMG12=GM2SICR(ia_12)
         SOMG21=GM2SICR(ia_21)

!-------------------Form-QM-Particle-Fock-Matrix------------------------(

       COEF1=0.0d+00
       COEF2=0.0d+00
       COEF3=0.0d+00

       if(LAB) then
          COEF1=DAE(iec1,jec1)*DBE(iec2,jec2)
     x         +DBE(iec1,jec1)*DAE(iec2,jec2)
       end if
       if(LAA) COEF2=DAE(iec1,jec1)*DAE(iec2,jec2)
       if(LBB) COEF3=DBE(iec1,jec1)*DBE(iec2,jec2)

       XFP(ip,jp)=XFP(ip,jp)
     x           +COEF1*OMG12
     x           +(COEF2+COEF3)*(OMG12-OMG21)

       XSP(ip,jp)=XSP(ip,jp)
     x           +COEF1*SOMG12
     x           +(COEF2+COEF3)*(SOMG12-SOMG21)

! Energy and Overlap
       E_P_OMG2=E_P_OMG2+DP(ip,jp)*(
     x            COEF1*OMG12
     x           +(COEF2+COEF3)*(OMG12-OMG21) )

       S_P_OMG2=S_P_OMG2+DP(ip,jp)*(
     x            COEF1*SOMG12
     x           +(COEF2+COEF3)*(SOMG12-SOMG21) )

!-------------------Form-QM-Particle-Fock-Matrix------------------------)

! Energy and Overlap
!      E_P_OMG2 = E_P_OMG2+DP(ip,jp)*(COEF1*OMG12-(COEF2+COEF3)*OMG21)
!      S_P_OMG2 = S_P_OMG2+DP(ip,jp)*(COEF1*SOMG12-(COEF2+COEF3)*SOMG21)

!!!!!!!!!!!!-previous
!-------------------Form-QM-Particle-Fock-Matrix------------------------(
!
!      COEF1=DETOT(iec1,jec1)*DETOT(iec2,jec2)
!
!      COEF2=0.0d+00
!
!      if(LAB) then
!         COEF2=DAE(iec1,jec1)*DBE(iec2,jec2)
!    x         +DBE(iec1,jec1)*DAE(iec2,jec2)
!      end if
!
!      XFP(ip,jp)=XFP(ip,jp)+COEF1*OMG12-COEF2*OMG21
!      XSP(ip,jp)=XSP(ip,jp)+COEF1*SOMG12-COEF2*SOMG21
!
!-------------------Form-QM-Particle-Fock-Matrix------------------------)
!
! Energy and Overlap
!      E_P_OMG2 = E_P_OMG2 + DP(ip,jp)*(COEF1*OMG12-COEF2*OMG21)
!      S_P_OMG2 = S_P_OMG2 + DP(ip,jp)*(COEF1*SOMG12-COEF2*SOMG21)
!!!!!!!!!!!!-previous

         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end

