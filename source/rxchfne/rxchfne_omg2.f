!=======================================================================
      subroutine RXCHFne_FAE_OMG2(Nchunks,nebf,npbf,ng2,
     x                          DAE,DBE,DP,GM2ICR,FAE,E_AE_OMG2)

! Calculate regular electron or special electron Fock matrices
! and contribution to total energy for 3-particle terms (same for both)
! AO contracted integrals are stored in-core.
! DAE represents density matrix for electrons coinciding with FAE
! DBE represents density matrix for electrons not coinciding with FAE
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng2,nebf,npbf
      double precision GM2ICR(ng2) ! Unsymmetrized OMG2 integrals
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FAE(nebf,nebf)
!     double precision FBE(nebf,nebf)
!     double precision SBE(nebf,nebf)
!     double precision FP(npbf,npbf)
!     double precision SP(npbf,npbf)
      double precision E_AE_OMG2

! Local Variables
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
      XFAE=0.0d+00
!     XFBE=0.0d+00
!     XSBE=0.0d+00
!     XFP =0.0d+00
!     XSP =0.0d+00

! HACK to form total electronic density:
!     DETOT=0.0d+00
!     call add2fock(nebf,DAE,DETOT)
!     call add2fock(nebf,DBE,DETOT)
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

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

         call RXCHFne_thread_FAEOMG2(istart,iend,ng2_seg,ng2,nebf,npbf,
     x                      loop_map,DAE,DBE,DP,GM2ICR,
     x                      XFAE,E_AE_OMG2)

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
      subroutine RXCHFne_thread_FAEOMG2(istart,iend,ng2_seg,ng2,nebf,
     x                      npbf,loop_map,DAE,DBE,DP,GM2ICR,
     x                      XFAE,E_AE_OMG2)
!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng2_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng2

      integer loop_map(ng2_seg,6)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM2ICR(ng2) ! Unsymmetrized OMG2 integrals

! Variables Returned
      double precision XFAE(nebf,nebf)
      double precision E_AE_OMG2

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap,ia
      integer ia_12
      integer ia_21

      double precision val
!     double precision COEF1,COEF2

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM2ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(ia_12,ia_21)
!$ompx private(val)
!$ompx reduction(+:XFAE)
!$ompx reduction(+:E_AE_OMG2)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec2=loop_map(imap,1)
         iec2=loop_map(imap,2)
         jec1=loop_map(imap,3)
         iec1=loop_map(imap,4)
         jp =loop_map(imap,5)
         ip =loop_map(imap,6)

C ARS( particle 1: special e ; particle 2: regular e ; index 3: prot )
         call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec1,iec2,jec2,ia_12)

         val=GM2ICR(ia_12)

!----------------Form-Alpha-Fock-Matrix-------------------(

         XFAE(iec2,jec2)=XFAE(iec2,jec2)+DP(ip,jp)*DBE(iec1,jec1)*val
         
!----------------Form-Alpha-Fock-Matrix-------------------)

! Energy
         E_AE_OMG2=E_AE_OMG2+DP(ip,jp)*DAE(iec2,jec2)*DBE(iec1,jec1)*val
         
         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!=======================================================================
      subroutine RXCHFne_FP_OMG2(Nchunks,nebf,npbf,ng2,NAE,NBE,
     x                    DAE,DBE,DP,GM2ICR,FP,E_P_OMG2)

! Calculate QM particle Fock matrices
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng2,nebf,npbf
      integer NAE,NBE
      double precision GM2ICR(ng2) ! Unsymmetrized OMG2 integrals
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FP(npbf,npbf)
      double precision E_P_OMG2

! Local Variables
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
      E_P_OMG2 = 0.0d+00

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

         call RXCHFne_thread_FPOMG2(istart,iend,ng2_seg,ng2,nebf,npbf,
     x                      loop_map,DAE,DBE,DP,GM2ICR,XFP,E_P_OMG2)

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
      subroutine RXCHFne_thread_FPOMG2(istart,iend,ng2_seg,ng2,nebf,
     x                    npbf,loop_map,DAE,DBE,DP,GM2ICR,XFP,E_P_OMG2)

!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng2_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng2

      integer loop_map(ng2_seg,6)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM2ICR(ng2) ! Unsymmetrized OMG2 integrals

! Variables Returned
      double precision XFP(npbf,npbf)
      double precision E_P_OMG2

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap,ia
      integer ia_12
      integer ia_21
      double precision val

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM2ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(ia_12,ia_21)
!$ompx private(val)
!$ompx reduction(+:XFP)
!$ompx reduction(+:E_P_OMG2)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec2=loop_map(imap,1)
         iec2=loop_map(imap,2)
         jec1=loop_map(imap,3)
         iec1=loop_map(imap,4)
         jp =loop_map(imap,5)
         ip =loop_map(imap,6)

C ARS( particle 1: special e ; particle 2: regular e ; index 3: prot )
         call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec1,iec2,jec2,ia_12)

         val=GM2ICR(ia_12)

!-------------------Form-QM-Particle-Fock-Matrix------------------------(

       XFP(ip,jp)=XFP(ip,jp)+DAE(iec2,jec2)*DBE(iec1,jec1)*val

!-------------------Form-QM-Particle-Fock-Matrix------------------------)

! Energy
       E_P_OMG2=E_P_OMG2+DP(ip,jp)*DAE(iec2,jec2)*DBE(iec1,jec1)*val

         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!=======================================================================
      subroutine RXCHFne_FBE_OMG2(Nchunks,nebf,npbf,ng2,NAE,NBE,
     x                    DAE,DBE,DP,GM2ICR,FBE,E_BE_OMG2)

! Calculate special electronic Fock matrices
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng2,nebf,npbf
      integer NAE,NBE
      double precision GM2ICR(ng2) ! Unsymmetrized OMG2 integrals
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FBE(nebf,nebf)
      double precision E_BE_OMG2

! Local Variables
      integer istat,ichunk,istart,iend,ng2_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2
      integer,allocatable :: loop_map(:,:)
!     double precision XFAE(nebf,nebf)
!     double precision XSAE(nebf,nebf)
      double precision XFBE(nebf,nebf)
!     double precision XSBE(nebf,nebf)
!     double precision XFP(npbf,npbf)
!     double precision XSP(npbf,npbf)

!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
!     E_gam2=0.0d+00
!     XFAE=0.0d+00
!     XSAE=0.0d+00
!     XFBE=0.0d+00
!     XSBE=0.0d+00
      XFBE =0.0d+00
      E_BE_OMG2 = 0.0d+00

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

         call RXCHFne_thread_FBEOMG2(istart,iend,ng2_seg,ng2,nebf,npbf,
     x                      loop_map,DAE,DBE,DP,GM2ICR,XFBE,E_BE_OMG2)

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
      call add2fock(nebf,XFBE,FBE)
!     call add2fock(npbf,XSP,SP)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!======================================================================
      subroutine RXCHFne_thread_FBEOMG2(istart,iend,ng2_seg,ng2,nebf,
     x                   npbf,loop_map,DAE,DBE,DP,GM2ICR,XFBE,E_BE_OMG2)

!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng2_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng2

      integer loop_map(ng2_seg,6)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM2ICR(ng2) ! Unsymmetrized OMG2 integrals

! Variables Returned
      double precision XFBE(nebf,nebf)
      double precision E_BE_OMG2

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap,ia
      integer ia_12
      integer ia_21
      double precision val

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM2ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(ia_12,ia_21)
!$ompx private(val)
!$ompx reduction(+:XFBE)
!$ompx reduction(+:E_BE_OMG2)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec2=loop_map(imap,1)
         iec2=loop_map(imap,2)
         jec1=loop_map(imap,3)
         iec1=loop_map(imap,4)
         jp =loop_map(imap,5)
         ip =loop_map(imap,6)

C ARS( particle 1: special e ; particle 2: regular e ; index 3: prot )
         call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec1,iec2,jec2,ia_12)

         val=GM2ICR(ia_12)

!-------------------Form-Special-Electron-Fock-Matrix---------(

       XFBE(iec1,jec1)=XFBE(iec1,jec1)+DAE(iec2,jec2)*DP(ip,jp)*val

!-------------------Form-Special-Electron-Fock-Matrix---------)

! Energy
       E_BE_OMG2=E_BE_OMG2+DP(ip,jp)*DAE(iec2,jec2)*DBE(iec1,jec1)*val

         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end

