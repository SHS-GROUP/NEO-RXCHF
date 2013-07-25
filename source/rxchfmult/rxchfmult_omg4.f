C=======================================================================
      subroutine RXCHFmult_FAE_OMG4(Nchunks,nebf,npbf,ng4,
     x                              DAE,DBE,DP,GM4ICR,
     x                              FAE,E_AE_OMG4)

C Calculate regular electron Fock matrix contribution
C to total energy for 5-particle terms
C AO contracted integrals are stored in-core in GM4ICR
C DAE : regular electronic density matrix
C DBE : special electronic density matrix
C=======================================================================
      implicit none

C Input variables
      integer           Nchunks
      integer           ng4,nebf,npbf
      double precision  GM4ICR(ng4)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FAE(nebf,nebf)
      double precision  E_AE_OMG4

C Local variables
      integer           istat,ichunk,istart,iend,ng4_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4
      double precision  XFAE(nebf,nebf)

      integer,allocatable :: loop_map(:,:)

C Initialize
      XFAE      = 0.0d+00
      E_AE_OMG4 = 0.0d+00

C Chop up calculation of OMG4 terms
      do ichunk=1,Nchunks

         call loop_size(1,ng4,Nchunks,ichunk-1,istart,iend)

         ng4_seg=1+iend-istart ! segment of ng4

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )

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

         call RXCHFmult_thread_FAE_OMG4(istart,iend,ng4_seg,
     x                                  ng4,nebf,npbf,
     x                                  loop_map,DAE,DBE,DP,GM4ICR,
     x                                  XFAE,E_AE_OMG4)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
      call add2fock(nebf,XFAE,FAE)

      return
      end
C=======================================================================
      subroutine RXCHFmult_thread_FAE_OMG4(istart,iend,ng4_seg,
     x                                     ng4,nebf,npbf,
     x                                     loop_map,DAE,DBE,DP,GM4ICR,
     x                                     XFAE,E_AE_OMG4)
C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng4_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng4
      integer           loop_map(ng4_seg,10)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM4ICR(ng4)

C Output variables
      double precision  XFAE(nebf,nebf)
      double precision  E_AE_OMG4

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1  !
      integer           iec2,jec2  ! Contracted elec basis function indices
      integer           iec3,jec3  !
      integer           iec4,jec4  !
      integer           imap,ia
      double precision  val
      double precision  wtime

C GM4ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron 1 indicies
C    iec3,jec3 : special electron 2 indicies
C    iec4,jec4 : special electron 3 indicies

!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng4_seg)
!$ompx shared(ng4)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM4ICR)
!$ompx private(iLp) 
!$ompx private(imap,ia)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(iec4,jec4)
!$ompx private(val)
!$ompx reduction(+:XFAE)
!$ompx reduction(+:E_AE_OMG4)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec4=loop_map(imap,1)
         iec4=loop_map(imap,2)
         jec3=loop_map(imap,3)
         iec3=loop_map(imap,4)
         jec2=loop_map(imap,5)
         iec2=loop_map(imap,6)
         jec1=loop_map(imap,7)
         iec1=loop_map(imap,8)
         jp =loop_map(imap,9)
         ip =loop_map(imap,10)

         call index_GAM_4PK2(nebf,npbf,ip,jp,iec1,jec1,
     x                       iec2,jec2,iec3,jec3,iec4,jec4,ia)

         val = GM4ICR(ia)

         XFAE(iec1,jec1)=XFAE(iec1,jec1)+DP(ip,jp)*DBE(iec2,jec2)*
     x                   DBE(iec3,jec3)*DBE(iec4,jec4)*val
         E_AE_OMG4=E_AE_OMG4+DP(ip,jp)*DAE(iec1,jec1)*
     x             DBE(iec2,jec2)*DBE(iec3,jec3)*DBE(iec4,jec4)*val

      end do
!$omp end do
!$omp end parallel      

      return
      end
C=======================================================================
      subroutine RXCHFmult_FP_OMG4(Nchunks,nebf,npbf,ng4,
     x                             DAE,DBE,DP,GM4ICR,
     x                             FP,E_P_OMG4)

C Calculate quantum particle Fock matrix contribution
C to total energy for 5-particle terms
C AO contracted integrals are stored in-core in GM4ICR
C DAE : regular electronic density matrix
C DBE : special electronic density matrix
C=======================================================================
      implicit none

C Input variables
      integer           Nchunks
      integer           ng4,nebf,npbf
      double precision  GM4ICR(ng4)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FP(npbf,npbf)
      double precision  E_P_OMG4

C Local variables
      integer           istat,ichunk,istart,iend,ng4_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4
      double precision  XFP(npbf,npbf)

      integer,allocatable :: loop_map(:,:)

C Initialize
      XFP      = 0.0d+00
      E_P_OMG4 = 0.0d+00

C Chop up calculation of OMG4 terms
      do ichunk=1,Nchunks

         call loop_size(1,ng4,Nchunks,ichunk-1,istart,iend)

         ng4_seg=1+iend-istart ! segment of ng4

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )

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

         call RXCHFmult_thread_FP_OMG4(istart,iend,ng4_seg,
     x                                 ng4,nebf,npbf,
     x                                 loop_map,DAE,DBE,DP,GM4ICR,
     x                                 XFP,E_P_OMG4)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
      call add2fock(npbf,XFP,FP)

      return
      end
C=======================================================================
      subroutine RXCHFmult_thread_FP_OMG4(istart,iend,ng4_seg,
     x                                    ng4,nebf,npbf,
     x                                    loop_map,DAE,DBE,DP,GM4ICR,
     x                                    XFP,E_P_OMG4)
C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng4_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng4
      integer           loop_map(ng4_seg,10)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM4ICR(ng4)

C Output variables
      double precision  XFP(npbf,npbf)
      double precision  E_P_OMG4

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1  !
      integer           iec2,jec2  ! Contracted elec basis function indices
      integer           iec3,jec3  !
      integer           iec4,jec4  !
      integer           imap,ia
      double precision  val
      double precision  wtime

C GM4ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron 1 indicies
C    iec3,jec3 : special electron 2 indicies
C    iec4,jec4 : special electron 3 indicies

!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng4_seg)
!$ompx shared(ng4)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM4ICR)
!$ompx private(iLp) 
!$ompx private(imap,ia)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(iec4,jec4)
!$ompx private(val)
!$ompx reduction(+:XFP)
!$ompx reduction(+:E_P_OMG4)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec4=loop_map(imap,1)
         iec4=loop_map(imap,2)
         jec3=loop_map(imap,3)
         iec3=loop_map(imap,4)
         jec2=loop_map(imap,5)
         iec2=loop_map(imap,6)
         jec1=loop_map(imap,7)
         iec1=loop_map(imap,8)
         jp =loop_map(imap,9)
         ip =loop_map(imap,10)

         call index_GAM_4PK2(nebf,npbf,ip,jp,iec1,jec1,
     x                       iec2,jec2,iec3,jec3,iec4,jec4,ia)

         val = GM4ICR(ia)

         XFP(ip,jp)=XFP(ip,jp)+DAE(iec1,jec1)*DBE(iec2,jec2)*
     x              DBE(iec3,jec3)*DBE(iec4,jec4)*val
         E_P_OMG4=E_P_OMG4+DP(ip,jp)*DAE(iec1,jec1)*
     x            DBE(iec2,jec2)*DBE(iec3,jec3)*DBE(iec4,jec4)*val

      end do
!$omp end do
!$omp end parallel      

      return
      end
C=======================================================================
      subroutine RXCHFmult_FBE_OMG4(Nchunks,nebf,npbf,ng4,
     x                              DAE,DBE,DP,GM4ICR,
     x                              FBE,E_BE_OMG4)

C Calculate special electron Fock matrix contribution
C to total energy for 5-particle terms
C AO contracted integrals are stored in-core in GM4ICR
C DAE : regular electronic density matrix
C DBE : special electronic density matrix
C=======================================================================
      implicit none

C Input variables
      integer           Nchunks
      integer           ng4,nebf,npbf
      double precision  GM4ICR(ng4)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FBE(nebf,nebf)
      double precision  E_BE_OMG4

C Local variables
      integer           istat,ichunk,istart,iend,ng4_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,iec4,jec4
      double precision  XFBE(nebf,nebf)

      integer,allocatable :: loop_map(:,:)

C Initialize
      XFBE      = 0.0d+00
      E_BE_OMG4 = 0.0d+00

C Chop up calculation of OMG4 terms
      do ichunk=1,Nchunks

         call loop_size(1,ng4,Nchunks,ichunk-1,istart,iend)

         ng4_seg=1+iend-istart ! segment of ng4

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )

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

         call RXCHFmult_thread_FBE_OMG4(istart,iend,ng4_seg,
     x                                  ng4,nebf,npbf,
     x                                  loop_map,DAE,DBE,DP,GM4ICR,
     x                                  XFBE,E_BE_OMG4)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
      call add2fock(nebf,XFBE,FBE)

      return
      end
C=======================================================================
      subroutine RXCHFmult_thread_FBE_OMG4(istart,iend,ng4_seg,
     x                                     ng4,nebf,npbf,
     x                                     loop_map,DAE,DBE,DP,GM4ICR,
     x                                     XFBE,E_BE_OMG4)
C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng4_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng4
      integer           loop_map(ng4_seg,10)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM4ICR(ng4)

C Output variables
      double precision  XFBE(nebf,nebf)
      double precision  E_BE_OMG4

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1  !
      integer           iec2,jec2  ! Contracted elec basis function indices
      integer           iec3,jec3  !
      integer           iec4,jec4  !
      integer           imap,ia
      integer           ia_1234
      integer           ia_1243
      integer           ia_1324
      integer           ia_1342
      integer           ia_1423
      integer           ia_1432
      double precision  val
      double precision  wtime
      double precision  three

      three=3.0d+00

C GM4ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron 1 indicies
C    iec3,jec3 : special electron 2 indicies
C    iec4,jec4 : special electron 3 indicies

!$omp parallel 
!$ompx shared(three)
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng4_seg)
!$ompx shared(ng4)
!$ompx shared(DAE,DBE,DP)
!$ompx shared(GM4ICR)
!$ompx private(iLp) 
!$ompx private(imap,ia)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(iec4,jec4)
!$ompx private(val)
!$ompx reduction(+:XFBE)
!$ompx reduction(+:E_BE_OMG4)

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
         jec4=loop_map(imap,1)
         iec4=loop_map(imap,2)
         jec3=loop_map(imap,3)
         iec3=loop_map(imap,4)
         jec2=loop_map(imap,5)
         iec2=loop_map(imap,6)
         jec1=loop_map(imap,7)
         iec1=loop_map(imap,8)
         jp =loop_map(imap,9)
         ip =loop_map(imap,10)

         call index_GAM_4PK2(nebf,npbf,ip,jp,iec1,jec1,
     x                       iec2,jec2,iec3,jec3,iec4,jec4,ia)

         val = GM4ICR(ia)

         XFBE(iec2,jec2)=XFBE(iec2,jec2)+three*DAE(iec1,jec1)*
     x                   DBE(iec3,jec3)*DBE(iec4,jec4)*DP(ip,jp)*val
         E_BE_OMG4=E_BE_OMG4+DP(ip,jp)*DAE(iec1,jec1)*
     x             DBE(iec2,jec2)*DBE(iec3,jec3)*DBE(iec4,jec4)*val

      end do
!$omp end do
!$omp end parallel      

      return
      end

