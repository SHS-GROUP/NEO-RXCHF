!=======================================================================
      subroutine RXCUHF_FE_OMG3(Nchunks,nebf,npbf,ng3,
     x                          DAalpE,DAbetE,DAtotE,DBE,DP,
     x                          GM3_1ICR,GM3_2ICR,
     x                          FAalpE,FAbetE,E_AE_OMG3)

! Calculate Alpha and Beta contributions
! to QM particle Fock matrix for 4-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng3,nebf,npbf
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3)
      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
      double precision DAtotE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FAalpE(nebf,nebf)
      double precision FAbetE(nebf,nebf)
      double precision E_AE_OMG3

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
      double precision XFAalpE(nebf,nebf)
      double precision XFAbetE(nebf,nebf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      XFAalpE =0.0d+00
      XFAbetE =0.0d+00
      E_AE_OMG3=0.0d+00
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

         call RXCUHF_thread_FEOMG3(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                             loop_map,DAalpE,DAbetE,DAtotE,DBE,DP,
     x                             GM3_1ICR,GM3_2ICR,
     x                             XFAalpE,XFAbetE,E_AE_OMG3)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
      call add2fock(nebf,XFAalpE,FAalpE)
      call add2fock(nebf,XFAbetE,FAbetE)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!=======================================================================
      subroutine RXCUHF_thread_FEOMG3(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                             loop_map,DAalpE,DAbetE,DAtotE,DBE,DP,
     x                                GM3_1ICR,GM3_2ICR,
     x                                XFAalpE,XFAbetE,E_AE_OMG3)
!=======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng3_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng3

      integer loop_map(ng3_seg,8)

      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
      double precision DAtotE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3)

! Variables Returned
      double precision XFAalpE(nebf,nebf)
      double precision XFAbetE(nebf,nebf)
      double precision E_AE_OMG3

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer iec3,jec3  !
      integer imap,ia
      integer ia_123_123
      integer ia_123_213
      integer ia_123_312
      integer ia_123_132
      integer ia_123_231
      integer ia_123_321
      integer ib_123_213
      integer ib_123_312
      integer ia_132_312
      integer ia_132_213
      integer ib_132_312
      integer ib_132_213

      double precision val1a,val1b,val2a,val2b,val2c,val2d
      double precision half,two

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

      half=0.50d+00
      two=2.0d+00

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(half,two)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DAalpE,DAbetE,DAtotE,DBE,DP)
!$ompx shared(GM3_1ICR,GM3_2ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia_123_123,ia_123_213,ia_123_312)
!$ompx private(ia_123_132,ia_123_231,ia_123_321)
!$ompx private(ib_123_213,ib_123_312)
!$ompx private(ia_132_312,ia_132_213)
!$ompx private(ib_132_312,ib_132_213)
!$ompx private(val1a,val1b,val2a,val2b,val2c,val2d)
!$ompx reduction(+:XFAalpE)
!$ompx reduction(+:XFAbetE)
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
     x                  ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia_123_123)
         call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec2,iec2,jec1,iec3,jec3,ia_123_213)
         call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec3,iec2,jec1,iec3,jec2,ia_123_312)
         call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec1,iec2,jec3,iec3,jec2,ia_123_132)
         call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec2,iec2,jec3,iec3,jec1,ia_123_231)
         call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec3,iec2,jec2,iec3,jec1,ia_123_321)

         call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,iec2,jec2,jec1,iec3,jec3,ib_123_213)
         call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec3,jec2,jec1,iec3,iec2,ib_123_312)
         call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec3,iec3,jec1,iec2,jec2,ia_132_312)
         call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec2,iec3,jec1,iec2,jec3,ia_132_213)
         call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,jec3,iec3,jec1,jec2,iec2,ib_132_312)
         call index_GAM_3PK(nebf,npbf,
     x                  ip,jp,iec1,iec2,iec3,jec1,jec2,jec3,ib_132_213)

         val1a=GM3_1ICR(ia_123_123)
         val1b=GM3_1ICR(ia_123_132)
         val2a=(GM3_2ICR(ia_123_213)+GM3_2ICR(ib_123_213))*half
         val2b=(GM3_2ICR(ia_123_312)+GM3_2ICR(ib_123_312))*half
         val2c=(GM3_2ICR(ia_132_312)+GM3_2ICR(ib_132_312))*half
         val2d=(GM3_2ICR(ia_132_213)+GM3_2ICR(ib_132_213))*half

      XFAalpE(iec2,jec2)=XFAalpE(iec2,jec2)+DBE(iec1,jec1)*DP(ip,jp)*
     x                     ( DAtotE(iec3,jec3)*val1a
     x                      -DAalpE(iec3,jec3)*val1b
     x                      -DAtotE(iec3,jec3)*val2a
     x                      +DAalpE(iec3,jec3)*val2b
     x                      -DAalpE(iec3,jec3)*val2c
     x                      +DAalpE(iec3,jec3)*val2d )

      XFAbetE(iec2,jec2)=XFAbetE(iec2,jec2)+DBE(iec1,jec1)*DP(ip,jp)*
     x                     ( DAtotE(iec3,jec3)*val1a
     x                      -DAbetE(iec3,jec3)*val1b
     x                      -DAalpE(iec3,jec3)*val2c )

      E_AE_OMG3=E_AE_OMG3+DP(ip,jp)*DBE(iec1,jec1)*
     x      ( half*DAtotE(iec2,jec2)*DAtotE(iec3,jec3)*val1a
     x       -half*(DAalpE(iec2,jec2)*DAalpE(iec3,jec3)+
     x              DAbetE(iec2,jec2)*DAbetE(iec3,jec3))*val1b
     x       -DAalpE(iec2,jec2)*DAtotE(iec3,jec3)*GM3_2ICR(ia_123_213)
     x       +DAalpE(iec2,jec2)*DAalpE(iec3,jec3)*GM3_2ICR(ia_123_312) )

         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end

!=======================================================================
      subroutine RXCUHF_FP_OMG3(Nchunks,nebf,npbf,ng3,
     x                          DAalpE,DAbetE,DAtotE,DBE,DP,
     x                          GM3_1ICR,GM3_2ICR,
     x                          FP,E_P_OMG3)

! Calculate Alpha and Beta contributions
! to QM particle Fock matrix for 4-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng3,nebf,npbf
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3)
      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
      double precision DAtotE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FP(npbf,npbf)
      double precision E_P_OMG3

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
      double precision XFP(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      XFP =0.0d+00
      E_P_OMG3=0.0d+00
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

         call RXCUHF_thread_FPOMG3(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                             loop_map,DAalpE,DAbetE,DAtotE,DBE,DP,
     x                             GM3_1ICR,GM3_2ICR,
     x                             XFP,E_P_OMG3)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
      call add2fock(npbf,XFP,FP)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!=======================================================================
      subroutine RXCUHF_thread_FPOMG3(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                             loop_map,DAalpE,DAbetE,DAtotE,DBE,DP,
     x                                GM3_1ICR,GM3_2ICR,
     x                                XFP,E_P_OMG3)
!=======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng3_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng3

      integer loop_map(ng3_seg,8)

      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
      double precision DAtotE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3)

! Variables Returned
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

      double precision val1a,val1b,val2a,val2b
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
!$ompx shared(DAalpE,DAbetE,DAtotE,DBE,DP)
!$ompx shared(GM3_1ICR,GM3_2ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia_123,ia_213,ia_312,ia_132,ia_231,ia_321)
!$ompx private(val1a,val1b,val2a,val2b)
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

         val1a=GM3_1ICR(ia_123)
         val1b=GM3_1ICR(ia_132)
         val2a=GM3_2ICR(ia_213)
         val2b=GM3_2ICR(ia_312)

      XFP(ip,jp)=XFP(ip,jp)+DBE(iec1,jec1)*
     x                ( half*DAtotE(iec2,jec2)*DAtotE(iec3,jec3)*val1a
     x                 -half*(DAalpE(iec2,jec2)*DAalpE(iec3,jec3)+
     x                        DAbetE(iec2,jec2)*DAbetE(iec3,jec3))*val1b
     x                 -DAalpE(iec2,jec2)*DAtotE(iec3,jec3)*val2a
     x                 +DAalpE(iec2,jec2)*DAalpE(iec3,jec3)*val2b )

      E_P_OMG3=E_P_OMG3+DP(ip,jp)*DBE(iec1,jec1)*
     x                ( half*DAtotE(iec2,jec2)*DAtotE(iec3,jec3)*val1a
     x                 -half*(DAalpE(iec2,jec2)*DAalpE(iec3,jec3)+
     x                        DAbetE(iec2,jec2)*DAbetE(iec3,jec3))*val1b
     x                 -DAalpE(iec2,jec2)*DAtotE(iec3,jec3)*val2a
     x                 +DAalpE(iec2,jec2)*DAalpE(iec3,jec3)*val2b )

         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!=======================================================================
      subroutine RXCUHF_FBE_OMG3(Nchunks,nebf,npbf,ng3,
     x                           DAalpE,DAbetE,DAtotE,DBE,DP,
     x                           GM3_1ICR,GM3_2ICR,
     x                           FBE,E_BE_OMG3)

! Calculate Alpha and Beta contributions
! to QM particle Fock matrix for 4-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng3,nebf,npbf
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3)
      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
      double precision DAtotE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FBE(nebf,nebf)
      double precision E_BE_OMG3

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
      double precision XFBE(nebf,nebf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      XFBE =0.0d+00
      E_BE_OMG3=0.0d+00
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

         call RXCUHF_thread_FBEOMG3(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                             loop_map,DAalpE,DAbetE,DAtotE,DBE,DP,
     x                              GM3_1ICR,GM3_2ICR,
     x                              XFBE,E_BE_OMG3)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
      call add2fock(nebf,XFBE,FBE)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!=======================================================================
      subroutine RXCUHF_thread_FBEOMG3(istart,iend,ng3_seg,ng3,
     x                                 nebf,npbf,loop_map,
     x                                 DAalpE,DAbetE,DAtotE,DBE,DP,
     x                                 GM3_1ICR,GM3_2ICR,
     x                                 XFBE,E_BE_OMG3)
!=======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng3_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng3

      integer loop_map(ng3_seg,8)

      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
      double precision DAtotE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM3_1ICR(ng3),GM3_2ICR(ng3)

! Variables Returned
      double precision XFBE(nebf,nebf)
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

      double precision val1a,val1b,val2a,val2b
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
!$ompx shared(DAalpE,DAbetE,DAtotE,DBE,DP)
!$ompx shared(GM3_1ICR,GM3_2ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia_123,ia_213,ia_312,ia_132,ia_231,ia_321)
!$ompx private(val1a,val1b,val2a,val2b)
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

         val1a=GM3_1ICR(ia_123)
         val1b=GM3_1ICR(ia_132)
         val2a=GM3_2ICR(ia_213)
         val2b=GM3_2ICR(ia_312)

      XFBE(iec1,jec1)=XFBE(iec1,jec1)+DP(ip,jp)*
     x                ( half*DAtotE(iec2,jec2)*DAtotE(iec3,jec3)*val1a
     x                 -half*(DAalpE(iec2,jec2)*DAalpE(iec3,jec3)+
     x                        DAbetE(iec2,jec2)*DAbetE(iec3,jec3))*val1b
     x                 -DAalpE(iec2,jec2)*DAtotE(iec3,jec3)*val2a
     x                 +DAalpE(iec2,jec2)*DAalpE(iec3,jec3)*val2b )

      E_BE_OMG3=E_BE_OMG3+DP(ip,jp)*DBE(iec1,jec1)*
     x                ( half*DAtotE(iec2,jec2)*DAtotE(iec3,jec3)*val1a
     x                 -half*(DAalpE(iec2,jec2)*DAalpE(iec3,jec3)+
     x                        DAbetE(iec2,jec2)*DAbetE(iec3,jec3))*val1b
     x                 -DAalpE(iec2,jec2)*DAtotE(iec3,jec3)*val2a
     x                 +DAalpE(iec2,jec2)*DAalpE(iec3,jec3)*val2b )

         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end

