!=======================================================================
      subroutine RXCHF_FAE_OMG3(Nchunks,nebf,npbf,ng3,
     x                    DAE,DBE,DP,GM3_1ICR,GM3_2ICR,GM3_3ICR,
     x                    FAE,E_AE_OMG3)

! Calculate Alpha and Beta contributions
! to QM particle Fock matrix for 4-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng3,nebf,npbf
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3),
     x                 GM3_3ICR(ng3)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FAE(nebf,nebf)
!     double precision FBE(nebf,nebf)
!     double precision FP(npbf,npbf)
      double precision E_AE_OMG3

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
      double precision XFAE(nebf,nebf)
!     double precision XFBE(nebf,nebf)
!     double precision XFP(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_AE_OMG3=0.0d+00
      XFAE =0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

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

         call RXCHF_thread_FAEOMG3(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                   loop_map,DAE,DBE,DP,GM3_1ICR,GM3_2ICR,GM3_3ICR,
     x                   XFAE,E_AE_OMG3)

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
!=======================================================================
      subroutine RXCHF_thread_FAEOMG3(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                   loop_map,DAE,DBE,DP,GM3_1ICR,GM3_2ICR,GM3_3ICR,
     x                   XFAE,E_AE_OMG3)
!=======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng3_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng3

      integer loop_map(ng3_seg,8)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3),
     x                 GM3_3ICR(ng3)

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

      double precision val1,val2,val3
      double precision half

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

      half=0.50d+00

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(half)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM3_1ICR,GM3_2ICR,GM3_3ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia_123,ia_213,ia_312,ia_132,ia_231,ia_321)
!$ompx private(val1,val2,val3)
!$ompx reduction(+:XFAE)
!$ompx reduction(+:E_AE_OMG3)

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

         val1=GM3_1ICR(ia_123)-half*GM3_1ICR(ia_132)
         val2=GM3_2ICR(ia_213)-half*GM3_2ICR(ia_312)
         val3=GM3_3ICR(ia_321)-half*GM3_3ICR(ia_231)

      XFAE(iec2,jec2)=XFAE(iec2,jec2)+
     x                DAE(iec3,jec3)*DBE(iec1,jec1)*DP(ip,jp)
     x               *(val1-val2+val3)

      E_AE_OMG3 = E_AE_OMG3 + half*DP(ip,jp)*DBE(iec1,jec1)
     x               *DAE(iec2,jec2)*DAE(iec3,jec3)
     x               *(val1-val2+val3)

         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end

!=======================================================================
      subroutine RXCHF_FP_OMG3(Nchunks,nebf,npbf,ng3,
     x                    DAE,DBE,DP,GM3_1ICR,GM3_2ICR,GM3_3ICR,
     x                    FP,E_P_OMG3)

! Calculate Alpha and Beta contributions
! to QM particle Fock matrix for 4-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng3,nebf,npbf
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3),
     x                 GM3_3ICR(ng3)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
!     double precision FAE(nebf,nebf)
!     double precision FBE(nebf,nebf)
      double precision FP(npbf,npbf)
      double precision E_P_OMG3

! Local Variables
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

         call RXCHF_thread_FPOMG3(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                   loop_map,DAE,DBE,DP,GM3_1ICR,GM3_2ICR,GM3_3ICR,
     x                   XFP,E_P_OMG3)

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
      subroutine RXCHF_thread_FPOMG3(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                   loop_map,DAE,DBE,DP,GM3_1ICR,GM3_2ICR,GM3_3ICR,
     x                   XFP,E_P_OMG3)
!=======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng3_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng3

      integer loop_map(ng3_seg,8)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3),
     x                 GM3_3ICR(ng3)

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

      double precision val1,val2,val3
      double precision half

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

      half=0.50d+00

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(half)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM3_1ICR,GM3_2ICR,GM3_3ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia_123,ia_213,ia_312,ia_132,ia_231,ia_321)
!$ompx private(val1,val2,val3)
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

         val1=GM3_1ICR(ia_123)-half*GM3_1ICR(ia_132)
         val2=GM3_2ICR(ia_213)-half*GM3_2ICR(ia_312)
         val3=GM3_3ICR(ia_321)-half*GM3_3ICR(ia_231)

      XFP(ip,jp)=XFP(ip,jp)+half
     x               *DAE(iec2,jec2)*DAE(iec3,jec3)*DBE(iec1,jec1)
     x               *(val1-val2+val3)

      E_P_OMG3 = E_P_OMG3 + half*DP(ip,jp)*DBE(iec1,jec1)
     x               *DAE(iec2,jec2)*DAE(iec3,jec3)
     x               *(val1-val2+val3)

         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!=======================================================================
      subroutine RXCHF_FBE_OMG3(Nchunks,nebf,npbf,ng3,
     x                    DAE,DBE,DP,GM3_1ICR,GM3_2ICR,GM3_3ICR,
     x                    FBE,E_BE_OMG3)

! Calculate Alpha and Beta contributions
! to QM particle Fock matrix for 4-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng3,nebf,npbf
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3),
     x                 GM3_3ICR(ng3)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
!     double precision FAE(nebf,nebf)
      double precision FBE(nebf,nebf)
!     double precision FP(npbf,npbf)
      double precision E_BE_OMG3

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
!     double precision XFAE(nebf,nebf)
      double precision XFBE(nebf,nebf)
!     double precision XFP(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_BE_OMG3=0.0d+00
      XFBE =0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

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

         call RXCHF_thread_FBEOMG3(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                   loop_map,DAE,DBE,DP,GM3_1ICR,GM3_2ICR,GM3_3ICR,
     x                   XFBE,E_BE_OMG3)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
!     call add2fock(nebf,XFAE,FAE)
      call add2fock(nebf,XFBE,FBE)
!     call add2fock(npbf,XFP,FP)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!=======================================================================
      subroutine RXCHF_thread_FBEOMG3(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                   loop_map,DAE,DBE,DP,GM3_1ICR,GM3_2ICR,GM3_3ICR,
     x                   XFBE,E_BE_OMG3)
!=======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng3_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng3

      integer loop_map(ng3_seg,8)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3),
     x                 GM3_3ICR(ng3)

! Variables Returned
!     double precision XFAE(nebf,nebf)
      double precision XFBE(nebf,nebf)
!     double precision XFP(npbf,npbf)
      double precision E_BE_OMG3

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

      double precision val1,val2,val3
      double precision half

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

      half=0.50d+00

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(half)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM3_1ICR,GM3_2ICR,GM3_3ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia_123,ia_213,ia_312,ia_132,ia_231,ia_321)
!$ompx private(val1,val2,val3)
!$ompx reduction(+:XFBE)
!$ompx reduction(+:E_BE_OMG3)

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

         val1=GM3_1ICR(ia_123)-half*GM3_1ICR(ia_132)
         val2=GM3_2ICR(ia_213)-half*GM3_2ICR(ia_312)
         val3=GM3_3ICR(ia_321)-half*GM3_3ICR(ia_231)

      XFBE(iec1,jec1)=XFBE(iec1,jec1)+half
     x               *DAE(iec2,jec2)*DAE(iec3,jec3)*DP(ip,jp)
     x               *(val1-val2+val3)

      E_BE_OMG3 = E_BE_OMG3 + half*DP(ip,jp)*DBE(iec1,jec1)
     x               *DAE(iec2,jec2)*DAE(iec3,jec3)
     x               *(val1-val2+val3)

         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end

