!=======================================================================
      subroutine RXCHF_FAE_OMG4(Nchunks,nebf,npbf,ng4,
     x                    DAE,DBE,DP,GM4ICR,FAE,E_AE_OMG4)

! Calculate QM Particle Fock matrices
! and contribution to total energy for 5-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng4,nebf,npbf
      double precision GM4ICR(ng4)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FAE(npbf,npbf)
      double precision E_AE_OMG4

! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4
      integer,allocatable :: loop_map(:,:)
      double precision XFAE(nebf,nebf)
!     double precision XFP(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_AE_OMG4=0.0d+00
      XFAE =0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(
      do ichunk=1,Nchunks

         call loop_size(1,ng4,Nchunks,ichunk-1,istart,iend)

! Segment of ng4:
         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )

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
                    do iec4=1,nebf
                    do jec4=1,nebf
                         imas=imas+1 ! imas is master_index
                         if(imas.ge.istart.and.imas.le.iend) then
                            Loopi=Loopi+1
                            loop_map(Loopi,1)=jec4
                            loop_map(Loopi,2)=iec4
                            loop_map(Loopi,3)=jec3
                            loop_map(Loopi,4)=iec3
                            loop_map(Loopi,5)=jec2
                            loop_map(Loopi,6)=iec2
                            loop_map(Loopi,7)=jec1
                            loop_map(Loopi,8)=iec1
                            loop_map(Loopi,9)=jp
                            loop_map(Loopi,10)=ip
                         end if

                    end do
                    end do
                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

         call RXCHF_thread_FAEOMG4(istart,iend,ng4_seg,ng4,nebf,npbf,
     x                       loop_map,DAE,DBE,DP,GM4ICR,XFAE,E_AE_OMG4)

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
      subroutine RXCHF_thread_FAEOMG4(istart,iend,ng4_seg,ng4,nebf,npbf,
     x                       loop_map,DAE,DBE,DP,GM4ICR,XFAE,E_AE_OMG4)

!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng4_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng4

      integer loop_map(ng4_seg,10)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM4ICR(ng4)

! Variables Returned
      double precision XFAE(npbf,npbf)
      double precision E_AE_OMG4

! Local Variables
      integer ip,jp
      integer ie1,je1  !
      integer ie2,je2  ! Contracted elec basis function indices
      integer ie3,je3  !
      integer ie4,je4  !
      integer imap,ia_3124,ia_3142,ia_2134,ia_2143
      double precision val1,val2,val3,val4
      double precision half

!---OPENMP-RELATED-VARIABLES-----(
!     integer IFIL
!     integer id
      integer loopi,iLP
!     double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

      half=0.50d+00

!---OPENMP-TIMING------------------------------------------------------(
!     wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(half)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng4_seg)
!$ompx shared(ng4)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM4ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1)
!$ompx private(ie2,je2)
!$ompx private(ie3,je3)
!$ompx private(ie4,je4)
!$ompx private(ia_3124,ia_3142,ia_2134,ia_2143)
!$ompx private(val1,val2,val3,val4)
!$ompx reduction(+:XFAE)
!$ompx reduction(+:E_AE_OMG4)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         je4=loop_map(imap,1)
         ie4=loop_map(imap,2)
         je3=loop_map(imap,3)
         ie3=loop_map(imap,4)
         je2=loop_map(imap,5)
         ie2=loop_map(imap,6)
         je1=loop_map(imap,7)
         ie1=loop_map(imap,8)
         jp =loop_map(imap,9)
         ip =loop_map(imap,10)


         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je1,ie3,je2,ie4,je4,ia_3124)
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je1,ie3,je4,ie4,je2,ia_3142)
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je1,ie3,je3,ie4,je4,ia_2134)
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je1,ie3,je4,ie4,je3,ia_2143)

         val1=GM4ICR(ia_3124)
         val2=GM4ICR(ia_3142)
         val3=GM4ICR(ia_2134)
         val4=GM4ICR(ia_2143)

         XFAE(ie2,je2)=XFAE(ie2,je2)+half*DP(ip,jp)*DBE(ie1,je1)
     x                  *DAE(ie3,je3)*DAE(ie4,je4)
     x                  *(val1-half*val2-val3+half*val4)

         E_AE_OMG4 = E_AE_OMG4 + half*half*DP(ip,jp)*DBE(ie1,je1)
     x                  *DAE(ie2,je2)*DAE(ie3,je3)*DAE(ie4,je4)
     x                  *(val1-half*val2-val3+half*val4)


      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!=======================================================================
      subroutine RXCHF_FP_OMG4(Nchunks,nebf,npbf,ng4,
     x                    DAE,DBE,DP,GM4ICR,FP,E_P_OMG4)

! Calculate QM Particle Fock matrices
! and contribution to total energy for 5-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng4,nebf,npbf
      double precision GM4ICR(ng4)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FP(npbf,npbf)
      double precision E_P_OMG4

! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4
      integer,allocatable :: loop_map(:,:)
!     double precision XFAE(nebf,nebf)
      double precision XFP(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_P_OMG4=0.0d+00
      XFP =0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(
      do ichunk=1,Nchunks

         call loop_size(1,ng4,Nchunks,ichunk-1,istart,iend)

! Segment of ng4:
         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )

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
                    do iec4=1,nebf
                    do jec4=1,nebf
                         imas=imas+1 ! imas is master_index
                         if(imas.ge.istart.and.imas.le.iend) then
                            Loopi=Loopi+1
                            loop_map(Loopi,1)=jec4
                            loop_map(Loopi,2)=iec4
                            loop_map(Loopi,3)=jec3
                            loop_map(Loopi,4)=iec3
                            loop_map(Loopi,5)=jec2
                            loop_map(Loopi,6)=iec2
                            loop_map(Loopi,7)=jec1
                            loop_map(Loopi,8)=iec1
                            loop_map(Loopi,9)=jp
                            loop_map(Loopi,10)=ip
                         end if

                    end do
                    end do
                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

         call RXCHF_thread_FPOMG4(istart,iend,ng4_seg,ng4,nebf,npbf,
     x                       loop_map,DAE,DBE,DP,GM4ICR,XFP,E_P_OMG4)

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
!======================================================================
      subroutine RXCHF_thread_FPOMG4(istart,iend,ng4_seg,ng4,nebf,npbf,
     x                       loop_map,DAE,DBE,DP,GM4ICR,XFP,E_P_OMG4)

!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng4_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng4

      integer loop_map(ng4_seg,10)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM4ICR(ng4)

! Variables Returned
      double precision XFP(npbf,npbf)
      double precision E_P_OMG4

! Local Variables
      integer ip,jp
      integer ie1,je1  !
      integer ie2,je2  ! Contracted elec basis function indices
      integer ie3,je3  !
      integer ie4,je4  !
      integer imap,ia_3124,ia_3142,ia_2134,ia_2143
      double precision val1,val2,val3,val4
      double precision half

!---OPENMP-RELATED-VARIABLES-----(
!     integer IFIL
!     integer id
      integer loopi,iLP
!     double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

      half=0.50d+00

!---OPENMP-TIMING------------------------------------------------------(
!     wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(half)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng4_seg)
!$ompx shared(ng4)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM4ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1)
!$ompx private(ie2,je2)
!$ompx private(ie3,je3)
!$ompx private(ie4,je4)
!$ompx private(ia_3124,ia_3142,ia_2134,ia_2143)
!$ompx private(val1,val2,val3,val4)
!$ompx reduction(+:XFP)
!$ompx reduction(+:E_P_OMG4)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         je4=loop_map(imap,1)
         ie4=loop_map(imap,2)
         je3=loop_map(imap,3)
         ie3=loop_map(imap,4)
         je2=loop_map(imap,5)
         ie2=loop_map(imap,6)
         je1=loop_map(imap,7)
         ie1=loop_map(imap,8)
         jp =loop_map(imap,9)
         ip =loop_map(imap,10)


         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je1,ie3,je2,ie4,je4,ia_3124)
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je1,ie3,je4,ie4,je2,ia_3142)
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je1,ie3,je3,ie4,je4,ia_2134)
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je1,ie3,je4,ie4,je3,ia_2143)

         val1=GM4ICR(ia_3124)
         val2=GM4ICR(ia_3142)
         val3=GM4ICR(ia_2134)
         val4=GM4ICR(ia_2143)

         XFP(ip,jp)=XFP(ip,jp)+half*half*DBE(ie1,je1)
     x                  *DAE(ie2,je2)*DAE(ie3,je3)*DAE(ie4,je4)
     x                  *(val1-half*val2-val3+half*val4)

         E_P_OMG4 = E_P_OMG4 + half*half*DBE(ie1,je1)*DP(ip,jp)
     x                  *DAE(ie2,je2)*DAE(ie3,je3)*DAE(ie4,je4)
     x                  *(val1-half*val2-val3+half*val4)


      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!=======================================================================
      subroutine RXCHF_FBE_OMG4(Nchunks,nebf,npbf,ng4,
     x                    DAE,DBE,DP,GM4ICR,FBE,E_BE_OMG4)

! Calculate QM Particle Fock matrices
! and contribution to total energy for 5-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng4,nebf,npbf
      double precision GM4ICR(ng4)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FBE(npbf,npbf)
      double precision E_BE_OMG4

! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4
      integer,allocatable :: loop_map(:,:)
      double precision XFBE(nebf,nebf)
!     double precision XFP(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_BE_OMG4=0.0d+00
      XFBE =0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(
      do ichunk=1,Nchunks

         call loop_size(1,ng4,Nchunks,ichunk-1,istart,iend)

! Segment of ng4:
         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )

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
                    do iec4=1,nebf
                    do jec4=1,nebf
                         imas=imas+1 ! imas is master_index
                         if(imas.ge.istart.and.imas.le.iend) then
                            Loopi=Loopi+1
                            loop_map(Loopi,1)=jec4
                            loop_map(Loopi,2)=iec4
                            loop_map(Loopi,3)=jec3
                            loop_map(Loopi,4)=iec3
                            loop_map(Loopi,5)=jec2
                            loop_map(Loopi,6)=iec2
                            loop_map(Loopi,7)=jec1
                            loop_map(Loopi,8)=iec1
                            loop_map(Loopi,9)=jp
                            loop_map(Loopi,10)=ip
                         end if

                    end do
                    end do
                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

         call RXCHF_thread_FBEOMG4(istart,iend,ng4_seg,ng4,nebf,npbf,
     x                       loop_map,DAE,DBE,DP,GM4ICR,XFBE,E_BE_OMG4)

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
!======================================================================
      subroutine RXCHF_thread_FBEOMG4(istart,iend,ng4_seg,ng4,nebf,npbf,
     x                       loop_map,DAE,DBE,DP,GM4ICR,XFBE,E_BE_OMG4)

!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng4_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng4

      integer loop_map(ng4_seg,10)

      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM4ICR(ng4)

! Variables Returned
      double precision XFBE(npbf,npbf)
      double precision E_BE_OMG4

! Local Variables
      integer ip,jp
      integer ie1,je1  !
      integer ie2,je2  ! Contracted elec basis function indices
      integer ie3,je3  !
      integer ie4,je4  !
      integer imap,ia_3124,ia_3142,ia_2134,ia_2143
      double precision val1,val2,val3,val4
      double precision half

!---OPENMP-RELATED-VARIABLES-----(
!     integer IFIL
!     integer id
      integer loopi,iLP
!     double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

      half=0.50d+00

!---OPENMP-TIMING------------------------------------------------------(
!     wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(half)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng4_seg)
!$ompx shared(ng4)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM4ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1)
!$ompx private(ie2,je2)
!$ompx private(ie3,je3)
!$ompx private(ie4,je4)
!$ompx private(ia_3124,ia_3142,ia_2134,ia_2143)
!$ompx private(val1,val2,val3,val4)
!$ompx reduction(+:XFBE)
!$ompx reduction(+:E_BE_OMG4)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         je4=loop_map(imap,1)
         ie4=loop_map(imap,2)
         je3=loop_map(imap,3)
         ie3=loop_map(imap,4)
         je2=loop_map(imap,5)
         ie2=loop_map(imap,6)
         je1=loop_map(imap,7)
         ie1=loop_map(imap,8)
         jp =loop_map(imap,9)
         ip =loop_map(imap,10)


         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je1,ie3,je2,ie4,je4,ia_3124)
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je1,ie3,je4,ie4,je2,ia_3142)
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je1,ie3,je3,ie4,je4,ia_2134)
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je1,ie3,je4,ie4,je3,ia_2143)

         val1=GM4ICR(ia_3124)
         val2=GM4ICR(ia_3142)
         val3=GM4ICR(ia_2134)
         val4=GM4ICR(ia_2143)

         XFBE(ie1,je1)=XFBE(ie1,je1)+half*half*DP(ip,jp)
     x                  *DAE(ie2,je2)*DAE(ie3,je3)*DAE(ie4,je4)
     x                  *(val1-half*val2-val3+half*val4)

         E_BE_OMG4 = E_BE_OMG4 + half*half*DBE(ie1,je1)*DP(ip,jp)
     x                  *DAE(ie2,je2)*DAE(ie3,je3)*DAE(ie4,je4)
     x                  *(val1-half*val2-val3+half*val4)


      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
