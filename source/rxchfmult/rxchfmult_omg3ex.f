C=======================================================================
      subroutine RXCHFmult_FAE_OMG3ex(Nchunks,nebf,npbf,ng3,
     x                                DAE,DBE,DP,
     x                                GM3ICR,GM3ex1ICR,GM3ex2ICR,
     x                                FAE,E_AE_OMG3)

C Calculate regular electron Fock matrix contribution
C to total energy for 4-particle terms
C AO contracted integrals are stored in-core in GM3ICR
C DAE : regular electronic density matrix
C DBE : special electronic density matrix
C=======================================================================
      implicit none

C Input variables
      integer           Nchunks
      integer           ng3,nebf,npbf
      double precision  GM3ICR(ng3)
      double precision  GM3ex1ICR(ng3)
      double precision  GM3ex2ICR(ng3)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FAE(nebf,nebf)
      double precision  E_AE_OMG3

C Local variables
      integer           istat,ichunk,istart,iend,ng3_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      double precision  XFAE(nebf,nebf)

      integer,allocatable :: loop_map(:,:)

C Initialize
      XFAE      = 0.0d+00
      E_AE_OMG3 = 0.0d+00

C Chop up calculation of OMG3 terms
      do ichunk=1,Nchunks

         call loop_size(1,ng3,Nchunks,ichunk-1,istart,iend)

         ng3_seg=1+iend-istart ! segment of ng3

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )

C Nested loop compression for this chunk:
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

         call RXCHFmult_thread_FAE_OMG3ex(istart,iend,ng3_seg,
     x                                    ng3,nebf,npbf,
     x                                    loop_map,DAE,DBE,DP,
     x                                    GM3ICR,GM3ex1ICR,GM3ex2ICR,
     x                                    XFAE,E_AE_OMG3)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
      call add2fock(nebf,XFAE,FAE)

      return
      end
C=======================================================================
      subroutine RXCHFmult_thread_FAE_OMG3ex(istart,iend,ng3_seg,
     x                                       ng3,nebf,npbf,
     x                                       loop_map,DAE,DBE,DP,
     x                                       GM3ICR,GM3ex1ICR,GM3ex2ICR,
     x                                       XFAE,E_AE_OMG3)
C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng3_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng3
      integer           loop_map(ng3_seg,8)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM3ICR(ng3)
      double precision  GM3ex1ICR(ng3)
      double precision  GM3ex2ICR(ng3)

C Output variables
      double precision  XFAE(nebf,nebf)
      double precision  E_AE_OMG3

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1  !
      integer           iec2,jec2  ! Contracted elec basis function indices
      integer           iec3,jec3  !
      integer           imap,ia
      integer           ia_123
      integer           ia_132
      integer           ia_213
      integer           ia_312
      integer           ia_321
      integer           ia_231
      double precision  val,val2,val3
      double precision  wtime
      double precision  half

      half=0.50d+00

C GM3ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron 1 indicies
C    iec3,jec3 : special electron 2 indicies

!$omp parallel 
!$ompx shared(half)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM3ICR,GM3ex1ICR,GM3ex2ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia_123,ia_132)
!$ompx private(ia_213,ia_312)
!$ompx private(ia_321,ia_231)
!$ompx private(val,val2,val3)
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
     x                      ip,jp,iec1,jec1,iec2,jec3,iec3,jec2,ia_132)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec2,iec2,jec1,iec3,jec3,ia_213)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec3,iec2,jec1,iec3,jec2,ia_312)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec3,iec2,jec2,iec3,jec1,ia_321)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec2,iec2,jec3,iec3,jec1,ia_231)

         val=GM3ICR(ia_123)-half*GM3ICR(ia_132)
         val2=GM3ex1ICR(ia_213)-half*GM3ex1ICR(ia_312)
         val3=GM3ex2ICR(ia_321)-half*GM3ex2ICR(ia_231)

         XFAE(iec1,jec1)=XFAE(iec1,jec1)+(val-half*(val2+val3))*
     x                   DBE(iec2,jec2)*DBE(iec3,jec3)*DP(ip,jp)
         E_AE_OMG3=E_AE_OMG3+DP(ip,jp)*DAE(iec1,jec1)*
     x             DBE(iec2,jec2)*DBE(iec3,jec3)*(val-half*(val2+val3))

      end do
!$omp end do
!$omp end parallel      

      return
      end
C=======================================================================
      subroutine RXCHFmult_FP_OMG3ex(Nchunks,nebf,npbf,ng3,
     x                               DAE,DBE,DP,
     x                               GM3ICR,GM3ex1ICR,GM3ex2ICR,
     x                               FP,E_P_OMG3)

C Calculate quantum particle Fock matrix contribution
C to total energy for 4-particle terms
C AO contracted integrals are stored in-core in GM3ICR
C DAE : regular electronic density matrix
C DBE : special electronic density matrix
C=======================================================================
      implicit none

C Input variables
      integer           Nchunks
      integer           ng3,nebf,npbf
      double precision  GM3ICR(ng3)
      double precision  GM3ex1ICR(ng3)
      double precision  GM3ex2ICR(ng3)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FP(npbf,npbf)
      double precision  E_P_OMG3

C Local variables
      integer           istat,ichunk,istart,iend,ng3_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      double precision  XFP(npbf,npbf)

      integer,allocatable :: loop_map(:,:)

C Initialize
      XFP      = 0.0d+00
      E_P_OMG3 = 0.0d+00

C Chop up calculation of OMG3 terms
      do ichunk=1,Nchunks

         call loop_size(1,ng3,Nchunks,ichunk-1,istart,iend)

         ng3_seg=1+iend-istart ! segment of ng3

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )

C Nested loop compression for this chunk:
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

         call RXCHFmult_thread_FP_OMG3ex(istart,iend,ng3_seg,
     x                                   ng3,nebf,npbf,
     x                                   loop_map,DAE,DBE,DP,
     x                                   GM3ICR,GM3ex1ICR,GM3ex2ICR,
     x                                   XFP,E_P_OMG3)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
      call add2fock(npbf,XFP,FP)

      return
      end
C=======================================================================
      subroutine RXCHFmult_thread_FP_OMG3ex(istart,iend,ng3_seg,
     x                                    ng3,nebf,npbf,
     x                                    loop_map,DAE,DBE,DP,
     x                                    GM3ICR,GM3ex1ICR,GM3ex2ICR,
     x                                    XFP,E_P_OMG3)
C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng3_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng3
      integer           loop_map(ng3_seg,8)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM3ICR(ng3)
      double precision  GM3ex1ICR(ng3)
      double precision  GM3ex2ICR(ng3)

C Output variables
      double precision  XFP(npbf,npbf)
      double precision  E_P_OMG3

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1  !
      integer           iec2,jec2  ! Contracted elec basis function indices
      integer           iec3,jec3  !
      integer           imap,ia
      integer           ia_123
      integer           ia_132
      integer           ia_213
      integer           ia_312
      integer           ia_321
      integer           ia_231
      double precision  val,val2,val3
      double precision  wtime
      double precision  half

      half=0.50d+00

C GM3ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron 1 indicies
C    iec3,jec3 : special electron 2 indicies

!$omp parallel 
!$ompx shared(half)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM3ICR,GM3ex1ICR,GM3ex2ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia_123,ia_132)
!$ompx private(ia_213,ia_312)
!$ompx private(ia_321,ia_231)
!$ompx private(val,val2,val3)
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
     x                      ip,jp,iec1,jec1,iec2,jec3,iec3,jec2,ia_132)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec2,iec2,jec1,iec3,jec3,ia_213)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec3,iec2,jec1,iec3,jec2,ia_312)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec3,iec2,jec2,iec3,jec1,ia_321)
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec2,iec2,jec3,iec3,jec1,ia_231)

         val=GM3ICR(ia_123)-half*GM3ICR(ia_132)
         val2=GM3ex1ICR(ia_213)-half*GM3ex1ICR(ia_312)
         val3=GM3ex2ICR(ia_321)-half*GM3ex2ICR(ia_231)

         XFP(ip,jp)=XFP(ip,jp)+(val-half*(val2+val3))*
     x              DAE(iec1,jec1)*DBE(iec2,jec2)*DBE(iec3,jec3)
         E_P_OMG3=E_P_OMG3+DP(ip,jp)*DAE(iec1,jec1)*
     x            DBE(iec2,jec2)*DBE(iec3,jec3)*(val-half*(val2+val3))

      end do
!$omp end do
!$omp end parallel      

      return
      end
C=======================================================================
      subroutine RXCHFmult_FBE_OMG3ex(Nchunks,nebf,npbf,ng3,
     x                                DAE,DBE,DP,
     x                                GM3ICR,GM3ex1ICR,GM3ex2ICR,
     x                                FBE,E_BE_OMG3)

C Calculate special electron Fock matrix contribution
C to total energy for 4-particle terms
C AO contracted integrals are stored in-core in GM3ICR
C DAE : regular electronic density matrix
C DBE : special electronic density matrix
C=======================================================================
      implicit none

C Input variables
      integer           Nchunks
      integer           ng3,nebf,npbf
      double precision  GM3ICR(ng3)
      double precision  GM3ex1ICR(ng3)
      double precision  GM3ex2ICR(ng3)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FBE(nebf,nebf)
      double precision  E_BE_OMG3

C Local variables
      integer           istat,ichunk,istart,iend,ng3_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      double precision  XFBE(nebf,nebf)

      integer,allocatable :: loop_map(:,:)

C Initialize
      XFBE      = 0.0d+00
      E_BE_OMG3 = 0.0d+00

C Chop up calculation of OMG3 terms
      do ichunk=1,Nchunks

         call loop_size(1,ng3,Nchunks,ichunk-1,istart,iend)

         ng3_seg=1+iend-istart ! segment of ng3

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )

C Nested loop compression for this chunk:
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

         call RXCHFmult_thread_FBE_OMG3ex(istart,iend,ng3_seg,
     x                                    ng3,nebf,npbf,
     x                                    loop_map,DAE,DBE,DP,
     x                                    GM3ICR,GM3ex1ICR,GM3ex2ICR,
     x                                    XFBE,E_BE_OMG3)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
      call add2fock(nebf,XFBE,FBE)

      return
      end
C=======================================================================
      subroutine RXCHFmult_thread_FBE_OMG3ex(istart,iend,ng3_seg,
     x                                     ng3,nebf,npbf,
     x                                     loop_map,DAE,DBE,DP,
     x                                     GM3ICR,GM3ex1ICR,GM3ex2ICR,
     x                                     XFBE,E_BE_OMG3)
C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng3_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng3
      integer           loop_map(ng3_seg,8)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM3ICR(ng3)
      double precision  GM3ex1ICR(ng3)
      double precision  GM3ex2ICR(ng3)

C Output variables
      double precision  XFBE(nebf,nebf)
      double precision  E_BE_OMG3

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1  !
      integer           iec2,jec2  ! Contracted elec basis function indices
      integer           iec3,jec3  !
      integer           imap,ia
      integer           ia_123_123
      integer           ia_123_132
      integer           ia_123_213
      integer           ia_123_312
      integer           ia_123_321
      integer           ia_123_231
      integer           ia_132_213
      integer           ia_132_312
      integer           ia_132_321
      integer           ia_132_231
      double precision  val,val2,val3
      double precision  wtime
      double precision  half,two

      half=0.50d+00
      two=2.0d+00

C GM3ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron 1 indicies
C    iec3,jec3 : special electron 2 indicies

!$omp parallel 
!$ompx shared(half,two)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM3ICR,GM3ex1ICR,GM3ex2ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia_123_123,ia_123_132)
!$ompx private(ia_123_213,ia_123_312)
!$ompx private(ia_123_321,ia_123_231)
!$ompx private(ia_132_213,ia_132_312)
!$ompx private(ia_132_321,ia_132_231)
!$ompx private(val,val2,val3)
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
     x                   ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia_123_123)
         call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec1,jec1,iec2,jec3,iec3,jec2,ia_123_132)
         call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec1,jec2,iec2,jec1,iec3,jec3,ia_123_213)
         call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec1,jec3,iec2,jec1,iec3,jec2,ia_123_312)
         call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec1,jec3,iec2,jec2,iec3,jec1,ia_123_321)
         call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec1,jec2,iec2,jec3,iec3,jec1,ia_123_231)
         call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec1,jec2,iec3,jec1,iec2,jec3,ia_132_213)
         call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec1,jec3,iec3,jec1,iec2,jec2,ia_132_312)
         call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec1,jec3,iec3,jec2,iec2,jec1,ia_132_321)
         call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec1,jec2,iec3,jec3,iec2,jec1,ia_132_231)

         val=GM3ICR(ia_123_123)-half*GM3ICR(ia_123_132)
         val2=GM3ex1ICR(ia_123_213)-half*GM3ex1ICR(ia_123_312)
     x       +GM3ex2ICR(ia_123_321)-half*GM3ex2ICR(ia_123_231)
         val3=GM3ex1ICR(ia_132_312)-half*GM3ex1ICR(ia_132_213)
     x       +GM3ex2ICR(ia_132_231)-half*GM3ex2ICR(ia_132_321)

         XFBE(iec2,jec2)=XFBE(iec2,jec2)+(two*val-half*(val2+val3))*
     x                   DAE(iec1,jec1)*DBE(iec3,jec3)*DP(ip,jp)
         E_BE_OMG3=E_BE_OMG3+DP(ip,jp)*DAE(iec1,jec1)*
     x             DBE(iec2,jec2)*DBE(iec3,jec3)*(val-half*val2)

      end do
!$omp end do
!$omp end parallel      

      return
      end

