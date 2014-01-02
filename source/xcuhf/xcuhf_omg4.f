!=======================================================================
      subroutine FAE_OMG4(Nchunks,nebf,npbf,ng4,NAE,NBE,
     x                    DAE,DBE,DP,GM4ICR,FAE,E_AE_OMG4)

! Calculate Alpha and Beta Fock matrices
! and contribution to total energy for 5-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng4,nebf,npbf
      integer NAE,NBE
      double precision GM4ICR(ng4)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FAE(nebf,nebf)
      double precision E_AE_OMG4

! Local Variables
      logical LAAA,LAAB,LBBA,LBBB
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4
      integer,allocatable :: loop_map(:,:)
      double precision XFAE(nebf,nebf)
!     double precision XFP(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_AE_OMG4=0.0d+00
      XFAE=0.0d+00
!     XFP =0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

      LAAA = .FALSE.
      LAAB = .FALSE.
      LBBA = .FALSE.
      LBBB = .FALSE.

! ? Alpha-Alpha-Alpha contribution to the Alpha Fock matrices
      LAAA = (NAE.ge.4)
! ? Beta-Beta-Beta contribution to the Alpha Fock matrices
      LBBB = (NBE.ge.3)
! ? Alpha-Alpha-Beta contribution
      LAAB = ( (NAE.ge.3) .and. (NBE.ge.1) )
! ? Beta-Beta-Alpha contribution
      LBBA = ( (NBE.ge.2) .and. (NAE.ge.2) )

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

         call thread_FAEOMG4(LAAA,LAAB,LBBA,LBBB,
     x                       istart,iend,ng4_seg,ng4,
     x                       nebf,npbf,loop_map,DAE,DBE,DP,GM4ICR,
     x                       XFAE,E_AE_OMG4)

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
      subroutine thread_FAEOMG4(LAAA,LAAB,LBBA,LBBB,
     x                       istart,iend,ng4_seg,ng4,
     x                       nebf,npbf,loop_map,DAE,DBE,DP,GM4ICR,
     x                       XFAE,E_AE_OMG4)

!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      logical LAAA,LAAB,LBBA,LBBB
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
      double precision XFAE(nebf,nebf)
      double precision E_AE_OMG4

! Local Variables
      integer ip,jp
      integer ie1,je1  !
      integer ie2,je2  ! Contracted elec basis function indices
      integer ie3,je3  !
      integer ie4,je4  !
      integer imap,ia
      double precision ans
      double precision four
!     parameter(four=4.0d+00)

!---OPENMP-RELATED-VARIABLES-----(
!     integer IFIL
!     integer id
      integer loopi,iLP
!     double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

      four=4.0d+00

!---OPENMP-TIMING------------------------------------------------------(
!     wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(LAAA,LAAB,LBBA,LBBB)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng4_seg)
!$ompx shared(ng4)
!$ompx shared(four)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM4ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1)
!$ompx private(ie2,je2)
!$ompx private(ie3,je3)
!$ompx private(ie4,je4)
!$ompx private(ia)
!$ompx private(ans)
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


         call calc_FAEOMG4(LAAA,LAAB,LBBA,LBBB,
     x                     nebf,npbf,ng4,ip,jp,
     x                     ie1,ie2,ie3,ie4,
     x                     je1,je2,je3,je4,
     x                     DAE,DBE,DP,GM4ICR,ans)

         XFAE(ie1,je1) = XFAE(ie1,je1) + four*ans

         E_AE_OMG4 = E_AE_OMG4 + DAE(ie1,je1)*ans

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!=======================================================================
      subroutine calc_FAEOMG4(LAAA,LAAB,LBBA,LBBB,
     x                        nebf,npbf,ng4,ip,jp,
     x                        ie1,ie2,ie3,ie4,
     x                        je1,je2,je3,je4,
     x                        DAE,DBE,DP,GM4ICR,ans)
!=======================================================================
      implicit none

! Input Variables
      logical LAAA,LAAB,LBBA,LBBB
      integer nebf,npbf,ng4,ip,jp
      integer ie1,ie2,ie3,ie4
      integer je1,je2,je3,je4
      double precision DAE(nebf,nebf),DBE(nebf,nebf),DP(npbf,npbf)
      double precision GM4ICR(ng4)

! Variables Returned
      double precision ans

! Local Variables
      integer ia
      double precision OMG
      double precision OMG_23      
      double precision OMG_12      
      double precision OMG_13      
      double precision OMG_24      
      double precision OMG_14      
      double precision OMG_34      
      double precision OMG_12_23   
      double precision OMG_13_23   
      double precision OMG_12_24   
      double precision OMG_12_14   
      double precision OMG_12_34   
      double precision OMG_13_34   
      double precision OMG_14_34   
      double precision OMG_13_24   
      double precision OMG_14_23   
      double precision OMG_23_34   
      double precision OMG_24_34   
      double precision OMG_12_23_34
      double precision OMG_12_23_24
      double precision OMG_12_13_34
      double precision OMG_13_23_24
      double precision OMG_14_24_34
      double precision OMG_14_23_34

      double precision XAAA,AAA_OMG
      double precision XAAB,AAB_OMG
      double precision XABA,ABA_OMG
      double precision XBAA,BAA_OMG
      double precision XABB,ABB_OMG
      double precision XBAB,BAB_OMG
      double precision XBBA,BBA_OMG
      double precision XBBB,BBB_OMG


! OMG
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je2,ie3,je3,ie4,je4,ia)
         OMG=GM4ICR(ia)
! OMG-23
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je3,ie3,je2,ie4,je4,ia)
         OMG_23=GM4ICR(ia)
! OMG-12
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je1,ie3,je3,ie4,je4,ia)
         OMG_12=GM4ICR(ia)
! OMG-12-23
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je3,ie3,je1,ie4,je4,ia)
         OMG_12_23=GM4ICR(ia)
! OMG-13-23
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je1,ie3,je2,ie4,je4,ia)
         OMG_13_23=GM4ICR(ia)
! OMG-13
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je2,ie3,je1,ie4,je4,ia)
         OMG_13=GM4ICR(ia)
! OMG-24
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je4,ie3,je3,ie4,je2,ia)
         OMG_24=GM4ICR(ia)
! OMG-12-24
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je4,ie3,je3,ie4,je1,ia)
         OMG_12_24=GM4ICR(ia)
! OMG-12-14
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je1,ie3,je3,ie4,je2,ia)
         OMG_12_14=GM4ICR(ia)
! OMG-14
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je2,ie3,je3,ie4,je1,ia)
         OMG_14=GM4ICR(ia)
! OMG-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je2,ie3,je4,ie4,je3,ia)
         OMG_34=GM4ICR(ia)
! OMG-12-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je1,ie3,je4,ie4,je3,ia)
         OMG_12_34=GM4ICR(ia)
! OMG-13-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je2,ie3,je4,ie4,je1,ia)
         OMG_13_34=GM4ICR(ia)
! OMG-14-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je2,ie3,je1,ie4,je3,ia)
         OMG_14_34=GM4ICR(ia)
! OMG-13-24
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je4,ie3,je1,ie4,je2,ia)
         OMG_13_24=GM4ICR(ia)
! OMG-14-23
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je3,ie3,je2,ie4,je1,ia)
         OMG_14_23=GM4ICR(ia)
! OMG-23-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je3,ie3,je4,ie4,je2,ia)
         OMG_23_34=GM4ICR(ia)
! OMG-24-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je4,ie3,je2,ie4,je3,ia)
         OMG_24_34=GM4ICR(ia)

!--------------- Triple permutations -------------------(
! OMG_12_23_34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je3,ie3,je4,ie4,je1,ia)
         OMG_12_23_34=GM4ICR(ia)
! OMG_12_23_24
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je4,ie3,je1,ie4,je3,ia)
         OMG_12_23_24=GM4ICR(ia)
! OMG_12_13_34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je1,ie3,je4,ie4,je2,ia)
         OMG_12_13_34=GM4ICR(ia)
! OMG_13_23_24
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je4,ie3,je2,ie4,je1,ia)
         OMG_13_23_24=GM4ICR(ia)
! OMG_14_24_34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je1,ie3,je2,ie4,je3,ia)
         OMG_14_24_34=GM4ICR(ia)
! OMG_14_23_34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je3,ie3,je1,ie4,je2,ia)
         OMG_14_23_34=GM4ICR(ia)
!--------------- Triple permutations -------------------)

! AAA
      XAAA=0.0d+00
      ! Get 24-permutation-term:  OMG_24P
      if(LAAA) AAA_OMG = OMG
     x                 - OMG_23         
     x                 - OMG_12         
     x                 - OMG_13         
     x                 - OMG_24         
     x                 - OMG_14         
     x                 - OMG_34         
     x                 + OMG_12_23      
     x                 + OMG_13_23      
     x                 + OMG_12_24      
     x                 + OMG_12_14      
     x                 + OMG_12_34      
     x                 + OMG_13_34      
     x                 + OMG_14_34      
     x                 + OMG_13_24      
     x                 + OMG_14_23      
     x                 + OMG_23_34      
     x                 + OMG_24_34      
     x                 - OMG_12_23_34   
     x                 - OMG_12_23_24   
     x                 - OMG_12_13_34   
     x                 - OMG_13_23_24   
     x                 - OMG_14_24_34   
     x                 - OMG_14_23_34

      if(LAAA) XAAA = DAE(ie2,je2)*DAE(ie3,je3)*DAE(ie4,je4)*AAA_OMG

! AAB
! OMG - OMG-23 - OMG-12 + OMG-12-23 + OMG-13-23 - OMG-13
      XAAB=0.0d+00
      if(LAAB) AAB_OMG = 
     x    (   OMG
     x      - OMG_23
     x      - OMG_12
     x      + OMG_12_23
     x      + OMG_13_23
     x      - OMG_13   )

      if(LAAB) XAAB = DAE(ie2,je2)*DAE(ie3,je3)*DBE(ie4,je4)*AAB_OMG

! ABA
! OMG - OMG-24 - OMG-12 + OMG-12-24 + OMG-12-14 - OMG-14
      XABA=0.0d+00
      if(LAAB) ABA_OMG = 
     x    (   OMG
     x      - OMG_24
     x      - OMG_12
     x      + OMG_12_24
     x      + OMG_12_14
     x      - OMG_14   )

      if(LAAB) XABA = DAE(ie2,je2)*DBE(ie3,je3)*DAE(ie4,je4)*ABA_OMG

! BAA
! OMG - OMG-34 - OMG-13 + OMG-13-34 + OMG-14-34 - OMG-14
      XBAA=0.0d+00
      if(LAAB) BAA_OMG = 
     x    (   OMG
     x      - OMG_34
     x      - OMG_13
     x      + OMG_13_34
     x      + OMG_14_34
     x      - OMG_14   )

      if(LAAB) XBAA = DBE(ie2,je2)*DAE(ie3,je3)*DAE(ie4,je4)*BAA_OMG

! ABB
! OMG - OMG-34 - OMG-12 + OMG-12-34
      XABB=0.0d+00
      if(LBBA) ABB_OMG = 
     x    (   OMG
     x      - OMG_34
     x      - OMG_12
     x      + OMG_12_34 )

      if(LBBA) XABB = DAE(ie2,je2)*DBE(ie3,je3)*DBE(ie4,je4)*ABB_OMG

! BAB
! OMG - OMG-24 - OMG-13 + OMG-13-24
      XBAB=0.0d+00
      if(LBBA) BAB_OMG = 
     x    (   OMG
     x      - OMG_24
     x      - OMG_13
     x      + OMG_13_24 )

      if(LBBA) XBAB = DBE(ie2,je2)*DAE(ie3,je3)*DBE(ie4,je4)*BAB_OMG

! BBA
! OMG - OMG-23 - OMG-14 + OMG-14-23
      XBBA=0.0d+00
      if(LBBA) BBA_OMG = 
     x    (   OMG
     x      - OMG_23
     x      - OMG_14
     x      + OMG_14_23 )

      if(LBBA) XBBA = DBE(ie2,je2)*DBE(ie3,je3)*DAE(ie4,je4)*BBA_OMG

! BBB
! OMG - OMG-34 - OMG-23 + OMG-23-34 + OMG-24-34 - OMG-24
      XBBB=0.0d+00
      if(LBBB) BBB_OMG = 
     x    (   OMG
     x      - OMG_34
     x      - OMG_23
     x      + OMG_23_34
     x      + OMG_24_34
     x      - OMG_24   )

      if(LBBB) XBBB = DBE(ie2,je2)*DBE(ie3,je3)*DBE(ie4,je4)*BBB_OMG


      ans = DP(ip,jp)*(XAAA+XAAB+XABA+XBAA+XABB+XBAB+XBBA+XBBB)


      return
      end
!=======================================================================
      subroutine UFP_OMG4(Nchunks,nebf,npbf,ng4,NAE,NBE,
     x                    DAE,DBE,DP,GM4ICR,FP,E_P_OMG4)

! Calculate QM Particle Fock matrices
! and contribution to total energy for 5-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng4,nebf,npbf
      integer NAE,NBE
      double precision GM4ICR(ng4)
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FP(npbf,npbf)
      double precision E_P_OMG4

! Local Variables
      logical L4A,L4B,L3A1B,L1A3B,L2A2B
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

      L4A   = .FALSE. 
      L4B   = .FALSE. 
      L3A1B = .FALSE. 
      L1A3B = .FALSE. 
      L2A2B = .FALSE. 

! ? Alpha-Alpha-Alpha-Alpha contribution to the Alpha Fock matrices
      L4A = (NAE.ge.4)
! ? Beta-Beta-Beta-Beta contribution to the Alpha Fock matrices
      L4B = (NBE.ge.4)
! ? Alpha-Alpha-Alpha-Beta contribution
      L3A1B = ( (NAE.ge.3) .and. (NBE.ge.1) )
! ? Beta-Beta-Beta-Alpha contributions
      L1A3B = ( (NAE.ge.1) .and. (NBE.ge.3) )
! ? Alpha-Alpha-Beta-Beta contributions
      L2A2B = ( (NAE.ge.2) .and. (NBE.ge.2) )

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

         call thread_UFPOMG4(L4A,L4B,L3A1B,L1A3B,L2A2B,
     x                       istart,iend,ng4_seg,ng4,
     x                       nebf,npbf,loop_map,DAE,DBE,DP,GM4ICR,
     x                       XFP,E_P_OMG4)

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
      subroutine thread_UFPOMG4(L4A,L4B,L3A1B,L1A3B,L2A2B,
     x                       istart,iend,ng4_seg,ng4,
     x                       nebf,npbf,loop_map,DAE,DBE,DP,GM4ICR,
     x                       XFP,E_P_OMG4)

!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      logical L4A,L4B,L3A1B,L1A3B,L2A2B
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
      integer imap,ia
      double precision ans
!     double precision four
!     parameter(four=4.0d+00)

!---OPENMP-RELATED-VARIABLES-----(
!     integer IFIL
!     integer id
      integer loopi,iLP
!     double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!---OPENMP-TIMING------------------------------------------------------(
!     wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(L4A,L4B,L3A1B,L1A3B,L2A2B)
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
!$ompx private(ia)
!$ompx private(ans)
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


         call calc_UFPOMG4(L4A,L4B,L3A1B,L1A3B,L2A2B,
     x                     nebf,npbf,ng4,ip,jp,
     x                     ie1,ie2,ie3,ie4,
     x                     je1,je2,je3,je4,
     x                     DAE,DBE,DP,GM4ICR,ans)

         XFP(ip,jp) = XFP(ip,jp) + ans

         E_P_OMG4 = E_P_OMG4 + DP(ip,jp)*ans


      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!=======================================================================
      subroutine calc_UFPOMG4(L4A,L4B,L3A1B,L1A3B,L2A2B,
     x                        nebf,npbf,ng4,ip,jp,
     x                        ie1,ie2,ie3,ie4,
     x                        je1,je2,je3,je4,
     x                        DAE,DBE,DP,GM4ICR,ans)
!=======================================================================
      implicit none

! Input Variables
      logical L4A,L4B,L3A1B,L1A3B,L2A2B
      integer nebf,npbf,ng4,ip,jp
      integer ie1,ie2,ie3,ie4
      integer je1,je2,je3,je4
      double precision DAE(nebf,nebf),DBE(nebf,nebf),DP(npbf,npbf)
      double precision GM4ICR(ng4)

! Variables Returned
      double precision ans

! Local Variables
      integer ia
      double precision OMG
      double precision OMG_23      
      double precision OMG_12      
      double precision OMG_13      
      double precision OMG_24      
      double precision OMG_14      
      double precision OMG_34      
      double precision OMG_12_23   
      double precision OMG_13_23   
      double precision OMG_12_24   
      double precision OMG_12_14   
      double precision OMG_12_34   
      double precision OMG_13_34   
      double precision OMG_14_34   
      double precision OMG_13_24   
      double precision OMG_14_23   
      double precision OMG_23_34   
      double precision OMG_24_34   
      double precision OMG_12_23_34
      double precision OMG_12_23_24
      double precision OMG_12_13_34
      double precision OMG_13_23_24
      double precision OMG_14_24_34
      double precision OMG_14_23_34

      double precision XAAAA,XAAAB,XAABA,XABAA,XAABB,XABAB,XABBA,XABBB 
      double precision XBBBB,XBBBA,XBBAB,XBABB,XBBAA,XBABA,XBAAB,XBAAA 

      double precision AAAA_OMG 
      double precision AAAB_OMG 
      double precision AABA_OMG 
      double precision ABAA_OMG 
      double precision AABB_OMG 
      double precision ABAB_OMG 
      double precision ABBA_OMG 
      double precision ABBB_OMG 
      double precision BBBB_OMG 
      double precision BBBA_OMG 
      double precision BBAB_OMG 
      double precision BABB_OMG 
      double precision BBAA_OMG 
      double precision BABA_OMG 
      double precision BAAB_OMG 
      double precision BAAA_OMG 


! OMG
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je2,ie3,je3,ie4,je4,ia)
         OMG=GM4ICR(ia)
! OMG-23
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je3,ie3,je2,ie4,je4,ia)
         OMG_23=GM4ICR(ia)
! OMG-12
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je1,ie3,je3,ie4,je4,ia)
         OMG_12=GM4ICR(ia)
! OMG-12-23
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je3,ie3,je1,ie4,je4,ia)
         OMG_12_23=GM4ICR(ia)
! OMG-13-23
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je1,ie3,je2,ie4,je4,ia)
         OMG_13_23=GM4ICR(ia)
! OMG-13
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je2,ie3,je1,ie4,je4,ia)
         OMG_13=GM4ICR(ia)
! OMG-24
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je4,ie3,je3,ie4,je2,ia)
         OMG_24=GM4ICR(ia)
! OMG-12-24
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je4,ie3,je3,ie4,je1,ia)
         OMG_12_24=GM4ICR(ia)
! OMG-12-14
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je1,ie3,je3,ie4,je2,ia)
         OMG_12_14=GM4ICR(ia)
! OMG-14
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je2,ie3,je3,ie4,je1,ia)
         OMG_14=GM4ICR(ia)
! OMG-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je2,ie3,je4,ie4,je3,ia)
         OMG_34=GM4ICR(ia)
! OMG-12-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je1,ie3,je4,ie4,je3,ia)
         OMG_12_34=GM4ICR(ia)
! OMG-13-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je2,ie3,je4,ie4,je1,ia)
         OMG_13_34=GM4ICR(ia)
! OMG-14-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je2,ie3,je1,ie4,je3,ia)
         OMG_14_34=GM4ICR(ia)
! OMG-13-24
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je4,ie3,je1,ie4,je2,ia)
         OMG_13_24=GM4ICR(ia)
! OMG-14-23
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je3,ie3,je2,ie4,je1,ia)
         OMG_14_23=GM4ICR(ia)
! OMG-23-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je3,ie3,je4,ie4,je2,ia)
         OMG_23_34=GM4ICR(ia)
! OMG-24-34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je1,ie2,je4,ie3,je2,ie4,je3,ia)
         OMG_24_34=GM4ICR(ia)

!--------------- Triple permutations -------------------(
! OMG_12_23_34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je3,ie3,je4,ie4,je1,ia)
         OMG_12_23_34=GM4ICR(ia)
! OMG_12_23_24
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je2,ie2,je4,ie3,je1,ie4,je3,ia)
         OMG_12_23_24=GM4ICR(ia)
! OMG_12_13_34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je1,ie3,je4,ie4,je2,ia)
         OMG_12_13_34=GM4ICR(ia)
! OMG_13_23_24
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je3,ie2,je4,ie3,je2,ie4,je1,ia)
         OMG_13_23_24=GM4ICR(ia)
! OMG_14_24_34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je1,ie3,je2,ie4,je3,ia)
         OMG_14_24_34=GM4ICR(ia)
! OMG_14_23_34
         call index_GAM_4PK2(nebf,npbf,ip,jp,
     x                       ie1,je4,ie2,je3,ie3,je1,ie4,je2,ia)
         OMG_14_23_34=GM4ICR(ia)
!--------------- Triple permutations -------------------)

! AAAA / BBBB
      XAAAA=0.0d+00
      XBBBB=0.0d+00
      ! Get 24-permutation-term:  OMG_24P
      if(L4A .OR. L4B) AAAA_OMG = OMG
     x                          - OMG_23         
     x                          - OMG_12         
     x                          - OMG_13         
     x                          - OMG_24         
     x                          - OMG_14         
     x                          - OMG_34         
     x                          + OMG_12_23      
     x                          + OMG_13_23      
     x                          + OMG_12_24      
     x                          + OMG_12_14      
     x                          + OMG_12_34      
     x                          + OMG_13_34      
     x                          + OMG_14_34      
     x                          + OMG_13_24      
     x                          + OMG_14_23      
     x                          + OMG_23_34      
     x                          + OMG_24_34      
     x                          - OMG_12_23_34   
     x                          - OMG_12_23_24   
     x                          - OMG_12_13_34   
     x                          - OMG_13_23_24   
     x                          - OMG_14_24_34   
     x                          - OMG_14_23_34

      if(L4A) XAAAA = DAE(ie1,je1)*DAE(ie2,je2)
     x              * DAE(ie3,je3)*DAE(ie4,je4)*AAAA_OMG

      if(L4B) XBBBB = DBE(ie1,je1)*DBE(ie2,je2)
     x              * DBE(ie3,je3)*DBE(ie4,je4)*AAAA_OMG

! AAAB / BBBA
! OMG - OMG-23 - OMG-12 + OMG-12-23 + OMG-13-23 - OMG-13
      XAAAB=0.0d+00
      XBBBA=0.0d+00
      if(L3A1B .OR. L1A3B) AAAB_OMG = 
     x    (   OMG
     x      - OMG_23
     x      - OMG_12
     x      + OMG_12_23
     x      + OMG_13_23
     x      - OMG_13   )

      if(L3A1B) XAAAB = DAE(ie1,je1)*DAE(ie2,je2)
     x                * DAE(ie3,je3)*DBE(ie4,je4)*AAAB_OMG

      if(L1A3B) XBBBA = DBE(ie1,je1)*DBE(ie2,je2)
     x                * DBE(ie3,je3)*DAE(ie4,je4)*AAAB_OMG

! AABA / BBAB
! OMG - OMG-24 - OMG-12 + OMG-12-24 + OMG-12-14 - OMG-14
      XAABA=0.0d+00
      XBBAB=0.0d+00
      if(L3A1B .OR. L1A3B) AABA_OMG = 
     x    (   OMG
     x      - OMG_24
     x      - OMG_12
     x      + OMG_12_24
     x      + OMG_12_14
     x      - OMG_14   )

      if(L3A1B) XAABA = DAE(ie1,je1)*DAE(ie2,je2)
     x                * DBE(ie3,je3)*DAE(ie4,je4)*AABA_OMG

      if(L1A3B) XBBAB = DBE(ie1,je1)*DBE(ie2,je2)
     x                * DAE(ie3,je3)*DBE(ie4,je4)*AABA_OMG

! ABAA / BABB
! OMG - OMG-34 - OMG-13 + OMG-13-34 + OMG-14-34 - OMG-14
      XABAA=0.0d+00
      XBABB=0.0d+00
      if(L3A1B .OR. L1A3B) ABAA_OMG = 
     x    (   OMG
     x      - OMG_34
     x      - OMG_13
     x      + OMG_13_34
     x      + OMG_14_34
     x      - OMG_14   )

      if(L3A1B) XABAA = DAE(ie1,je1)*DBE(ie2,je2)
     x                * DAE(ie3,je3)*DAE(ie4,je4)*ABAA_OMG

      if(L1A3B) XBABB = DBE(ie1,je1)*DAE(ie2,je2)
     x                * DBE(ie3,je3)*DBE(ie4,je4)*ABAA_OMG

! AABB / BBAA
! OMG - OMG-34 - OMG-12 + OMG-12-34
      XAABB=0.0d+00
      XBBAA=0.0d+00
      if(L2A2B) AABB_OMG = 
     x    (   OMG
     x      - OMG_34
     x      - OMG_12
     x      + OMG_12_34 )

      if(L2A2B) XAABB = DAE(ie1,je1)*DAE(ie2,je2)
     x                * DBE(ie3,je3)*DBE(ie4,je4)*AABB_OMG

      if(L2A2B) XBBAA = DBE(ie1,je1)*DBE(ie2,je2)
     x                * DAE(ie3,je3)*DAE(ie4,je4)*AABB_OMG

! ABAB / BABA
! OMG - OMG-24 - OMG-13 + OMG-13-24
      XABAB=0.0d+00
      XBABA=0.0d+00
      if(L2A2B) ABAB_OMG = 
     x    (   OMG
     x      - OMG_24
     x      - OMG_13
     x      + OMG_13_24 )

      if(L2A2B) XABAB = DAE(ie1,je1)*DBE(ie2,je2)
     x                * DAE(ie3,je3)*DBE(ie4,je4)*ABAB_OMG

      if(L2A2B) XBABA = DBE(ie1,je1)*DAE(ie2,je2)
     x                * DBE(ie3,je3)*DAE(ie4,je4)*ABAB_OMG

! ABBA / BAAB
! OMG - OMG-23 - OMG-14 + OMG-14-23
      XABBA=0.0d+00
      XBAAB=0.0d+00
      if(L2A2B) ABBA_OMG = 
     x    (   OMG
     x      - OMG_23
     x      - OMG_14
     x      + OMG_14_23 )

      if(L2A2B) XABBA = DAE(ie1,je1)*DBE(ie2,je2)
     x                * DBE(ie3,je3)*DAE(ie4,je4)*ABBA_OMG

      if(L2A2B) XBAAB = DBE(ie1,je1)*DAE(ie2,je2)
     x                * DAE(ie3,je3)*DBE(ie4,je4)*ABBA_OMG

! ABBB / BBBA
! OMG - OMG-34 - OMG-23 + OMG-23-34 + OMG-24-34 - OMG-24
      XABBB=0.0d+00
      XBAAA=0.0d+00
      if(L1A3B .OR. L3A1B) ABBB_OMG = 
     x    (   OMG
     x      - OMG_34
     x      - OMG_23
     x      + OMG_23_34
     x      + OMG_24_34
     x      - OMG_24   )

      if(L1A3B) XABBB = DAE(ie1,je1)*DBE(ie2,je2)
     x                * DBE(ie3,je3)*DBE(ie4,je4)*ABBB_OMG

      if(L3A1B) XBAAA = DBE(ie1,je1)*DAE(ie2,je2)
     x                * DAE(ie3,je3)*DAE(ie4,je4)*ABBB_OMG


!     ans = DP(ip,jp)*
!    x ( XAAAA + XAAAB + XAABA + XABAA + XAABB + XABAB + XABBA + XABBB 
!    x + XBBBB + XBBBA + XBBAB + XBABB + XBBAA + XBABA + XBAAB + XBAAA )

      ans = 
     x ( XAAAA + XAAAB + XAABA + XABAA + XAABB + XABAB + XABBA + XABBB 
     x + XBBBB + XBBBA + XBBAB + XBABB + XBBAA + XBABA + XBAAB + XBAAA )


      return
      end
!======================================================================
      subroutine OMG4_ICR(Nchunks,ne,np,ngee,ng2,ng4,
     x                    GAM_2s,GAM_4)

!======================================================================
      implicit none
      include 'omp_lib.h'
! Input Variables
      integer Nchunks
      integer ng4,ng2,ngee,ne,np
      double precision GAM_4(ng4)
      double precision GAM_2s(ng2)
      
! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4
      integer,allocatable :: loop_map(:,:)
      double precision,allocatable :: GAM_ee(:)
!     double precision,allocatable :: GAM_2s(:)
      double precision wtime,wtime2


         write(*,1000) ng4,nchunks,omp_get_num_procs(),
     xomp_get_max_threads(),1
C         wtime = omp_get_wtime()

!----READ-GAM_ee-AND-GAM_2s-INTO-MEMORY-------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      allocate( GAM_ee(ngee),stat=istat )
!     write(*,*) 'allocate GAM_ee: ',istat
      call read_GAM_ee(ne,ngee,GAM_ee) 
!----READ-GAM_ee-AND-GAM_2s-INTO-MEMORY-------------------------------)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(

      do ichunk=1,Nchunks

C         wtime2 = omp_get_wtime()

         call loop_size(1,ng4,Nchunks,ichunk-1,istart,iend)

         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )
!        write(*,*) 'allocate loop_map: ',istat

! Nested loop compression for this chunk:
         Loopi=0
         imas=0
         do ip=1,np
         do jp=1,np
            do ie1=1,ne
            do je1=1,ne
               do ie2=1,ne
               do je2=1,ne
                  do ie3=1,ne
                  do je3=1,ne
                    do ie4=1,ne
                    do je4=1,ne

                         imas=imas+1 ! imas is master_index
                         if(imas.ge.istart.and.imas.le.iend) then
                            Loopi=Loopi+1
                            loop_map(Loopi,1)=je4
                            loop_map(Loopi,2)=ie4
                            loop_map(Loopi,3)=je3
                            loop_map(Loopi,4)=ie3
                            loop_map(Loopi,5)=je2
                            loop_map(Loopi,6)=ie2
                            loop_map(Loopi,7)=je1
                            loop_map(Loopi,8)=ie1
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

         call thread_omg4_IC(ne,np,ngee,ng2,ng4,ng4_seg,istart,iend,
     x                       loop_map,GAM_ee,GAM_2s,GAM_4)

C         wtime2 = omp_get_wtime() - wtime2
C         write(*,2000)ichunk,wtime2

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
!     if(allocated(GAM_2s)) deallocate(GAM_2s)
      if(allocated(loop_map)) deallocate(loop_map)
!     if(allocated(GAM_4)) deallocate(GAM_4)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

C      wtime = omp_get_wtime() - wtime
C      write(*,3000)wtime


 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 5-PARTICLE INTEGRALS        |',/,
     x        6X,'|         FOR NEO-XCUHF CALCULATION           |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 5-PARTICLE INTEGRALS: ',1X,I12/
     x        8X,'  NUMBER OF BLOCKS (USER DEFINED): ',1X,I12/
     x        8X,'                          ',/,
     x        8X,'  COMPUTATIONAL RESOURCES:',/,
     x        8X,'  ------------------------',/,
     x        8X,'     CORES PER NODE:',1X,I3/
     x        8X,'          AVAILABLE:',1X,I3/
     x        8X,'    NUMBER OF NODES:',1X,I3/)
                       
 2000 FORMAT(8X,'    TIME TO EVALUATE BLOCK ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,'  TIMING SUMMARY FOR 5-PARTICLE INTEGRALS:',/,
     x        8X,'  ----------------------------------------',/,
     x        8X,'    TIME TO EVALUATE ALL INTEGRALS:',1X,F12.4)


      return
      end
!======================================================================
      subroutine thread_omg4_IC(ne,np,ngee,ng2,ng4,ng4_seg,istart,iend,
     x                          loop_map,GAM_ee,GAM_2s,GAM_4)
 
!======================================================================
      implicit none
      include 'omp_lib.h'

      integer istart,iend,imap
      integer ne    ! Number of contracted electronic basis functions
      integer np    ! Number of nuclear basis functions
      integer ngee  ! Number of contracted 2-electron integrals
      integer ng2   ! Number of contracted 3-particle integrals
      integer ng4   ! Number of contracted 5-particle integrals
      integer ng4_seg ! dimension of chunk of contracted 5-particle integrals

      integer loop_map(ng4_seg,10)
      double precision GAM_ee(ngee)   ! Array storage of 2e-integrals
      double precision GAM_2s(ng2)    ! Array storage of 3-part overlaps
      double precision GAM_4(ng4)     ! Array storage of 5-particle ints

! Local variables
      integer iLP
      integer ia
      integer ip
      integer jp
      integer ie1
      integer ie2
      integer ie3
      integer ie4
      integer je1
      integer je2
      integer je3
      integer je4

      double precision ans


!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(ne,np)
!$ompx shared(ngee)
!$ompx shared(ng2)
!$ompx shared(ng4_seg)
!$ompx shared(GAM_2s)
!$ompx shared(GAM_ee)
!$ompx shared(GAM_4)
!$ompx private(iLp) 
!$ompx private(imap) 
!$ompx private(ia)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(ie3,je3) 
!$ompx private(ie4,je4) 
!$ompx private(ans)

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

         call index_GAM_4PK2(ne,np,
     x                       ip,jp,
     x                       ie1,je1,
     x                       ie2,je2,
     x                       ie3,je3,
     x                       ie4,je4,ia)


!        call symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
!    x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

!     ii  i   jj  j   kk  k   ll  l   ip  jp
!     ie1,je1,ie2,je2,ie3,je3,ie4,je4,ip1,jp1,ans)

!     call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
!    x                         ip,jp,ii,jj,kk,ll,i,j,k,l,x1)

         call i10(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                         ip,jp,ie1,ie2,ie3,ie4,
     x                               je1,je2,je3,je4,ans)

         GAM_4(ia)=ans

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

     
      return
      end 
