C=======================================================================
      subroutine RXCHFmult_FAE_OMG2(Nchunks,nebf,npbf,ng2,
     x                              DAE,DBE,DP,GM2ICR,
     x                              FAE,E_AE_OMG2)

C Calculate regular electron Fock matrix contribution
C to total energy for 3-particle terms
C AO contracted integrals are stored in-core in GM2ICR
C DAE : regular electronic density matrix
C DBE : special electronic density matrix
C=======================================================================
      implicit none

C Input variables
      integer           Nchunks
      integer           ng2,nebf,npbf
      double precision  GM2ICR(ng2)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FAE(nebf,nebf)
      double precision  E_AE_OMG2

C Local variables
      integer           istat,ichunk,istart,iend,ng2_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2
      double precision  XFAE(nebf,nebf)

      integer,allocatable :: loop_map(:,:)

C Initialize
      XFAE      = 0.0d+00
      E_AE_OMG2 = 0.0d+00

C Chop up calculation of OMG2 terms
      do ichunk=1,Nchunks

         call loop_size(1,ng2,Nchunks,ichunk-1,istart,iend)

         ng2_seg=1+iend-istart ! segment of ng2

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng2_seg,6),stat=istat )

C Nested loop compression for this chunk:
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

         call RXCHFmult_thread_FAE_OMG2(istart,iend,ng2_seg,
     x                                  ng2,nebf,npbf,
     x                                  loop_map,DAE,DBE,DP,GM2ICR,
     x                                  XFAE,E_AE_OMG2)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
      call add2fock(nebf,XFAE,FAE)

      return
      end
C======================================================================
      subroutine RXCHFmult_thread_FAE_OMG2(istart,iend,ng2_seg,
     x                                     ng2,nebf,npbf,
     x                                     loop_map,DAE,DBE,DP,GM2ICR,
     x                                     XFAE,E_AE_OMG2)
C======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng2_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng2
      integer           loop_map(ng2_seg,6)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM2ICR(ng2)

C Output variables
      double precision  XFAE(nebf,nebf)
      double precision  E_AE_OMG2

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1
      integer           iec2,jec2
      integer           imap,ia
      integer           ia_12
      double precision  val
      double precision  wtime

C GM2ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron indicies

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
!$ompx private(ia_12)
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

         call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec1,iec2,jec2,ia_12)
         val=GM2ICR(ia_12)

         XFAE(iec1,jec1)=XFAE(iec1,jec1)+DP(ip,jp)*DBE(iec2,jec2)*val
         E_AE_OMG2=E_AE_OMG2+DP(ip,jp)*DAE(iec1,jec1)*DBE(iec2,jec2)*val
         
      end do
!$omp end do
!$omp end parallel      

      return
      end
C=======================================================================
      subroutine RXCHFmult_FP_OMG2(Nchunks,nebf,npbf,ng2,
     x                             DAE,DBE,DP,GM2ICR,
     x                             FP,E_P_OMG2)

C Calculate quantum particle Fock matrix contribution
C to total energy for 3-particle terms
C AO contracted integrals are stored in-core in GM2ICR
C DAE : regular electronic density matrix
C DBE : special electronic density matrix
C=======================================================================
      implicit none

C Input variables
      integer           Nchunks
      integer           ng2,nebf,npbf
      double precision  GM2ICR(ng2)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FP(npbf,npbf)
      double precision  E_P_OMG2

C Local variables
      integer           istat,ichunk,istart,iend,ng2_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2
      double precision  XFP(npbf,npbf)

      integer,allocatable :: loop_map(:,:)

C Initialize
      XFP      = 0.0d+00
      E_P_OMG2 = 0.0d+00

C Chop up calculation of OMG2 terms
      do ichunk=1,Nchunks

         call loop_size(1,ng2,Nchunks,ichunk-1,istart,iend)

         ng2_seg=1+iend-istart ! segment of ng2

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng2_seg,6),stat=istat )

C Nested loop compression for this chunk:
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

         call RXCHFmult_thread_FP_OMG2(istart,iend,ng2_seg,
     x                                 ng2,nebf,npbf,
     x                                 loop_map,DAE,DBE,DP,GM2ICR,
     x                                 XFP,E_P_OMG2)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
      call add2fock(npbf,XFP,FP)

      return
      end
C======================================================================
      subroutine RXCHFmult_thread_FP_OMG2(istart,iend,ng2_seg,
     x                                    ng2,nebf,npbf,
     x                                    loop_map,DAE,DBE,DP,GM2ICR,
     x                                    XFP,E_P_OMG2)
C======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng2_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng2
      integer           loop_map(ng2_seg,6)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM2ICR(ng2)

C Output variables
      double precision  XFP(npbf,npbf)
      double precision  E_P_OMG2

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1
      integer           iec2,jec2
      integer           imap,ia
      integer           ia_12
      double precision  val
      double precision  wtime

C GM2ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron indicies

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
!$ompx private(ia_12)
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

         call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec1,iec2,jec2,ia_12)
         val=GM2ICR(ia_12)

         XFP(ip,jp)=XFP(ip,jp)+DAE(iec1,jec1)*DBE(iec2,jec2)*val
         E_P_OMG2=E_P_OMG2+DP(ip,jp)*DAE(iec1,jec1)*DBE(iec2,jec2)*val

      end do
!$omp end do
!$omp end parallel      

      return
      end
C=======================================================================
      subroutine RXCHFmult_FBE_OMG2(Nchunks,nebf,npbf,ng2,
     x                              DAE,DBE,DP,GM2ICR,
     x                              FBE,E_BE_OMG2)

C Calculate special electron Fock matrix contribution
C to total energy for 3-particle terms
C AO contracted integrals are stored in-core in GM2ICR
C DAE : regular electronic density matrix
C DBE : special electronic density matrix
C=======================================================================
      implicit none

C Input variables
      integer           Nchunks
      integer           ng2,nebf,npbf
      double precision  GM2ICR(ng2)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)

C Output variables
      double precision  FBE(nebf,nebf)
      double precision  E_BE_OMG2

C Local variables
      integer           istat,ichunk,istart,iend,ng2_seg
      integer           Loopi,imas
      integer           ip,jp,iec1,jec1,iec2,jec2
      double precision  XFBE(nebf,nebf)

      integer,allocatable :: loop_map(:,:)

C Initialize
      XFBE      = 0.0d+00
      E_BE_OMG2 = 0.0d+00

C Chop up calculation of OMG2 terms
      do ichunk=1,Nchunks

         call loop_size(1,ng2,Nchunks,ichunk-1,istart,iend)

         ng2_seg=1+iend-istart ! segment of ng2

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng2_seg,6),stat=istat )

C Nested loop compression for this chunk:
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

         call RXCHFmult_thread_FBE_OMG2(istart,iend,ng2_seg,
     x                                  ng2,nebf,npbf,
     x                                  loop_map,DAE,DBE,DP,GM2ICR,
     x                                  XFBE,E_BE_OMG2)

      end do !end loop over chunks

      if(allocated(loop_map)) deallocate(loop_map)

C Update Fock matrix
      call add2fock(nebf,XFBE,FBE)

      return
      end
C======================================================================
      subroutine RXCHFmult_thread_FBE_OMG2(istart,iend,ng2_seg,
     x                                     ng2,nebf,npbf,
     x                                     loop_map,DAE,DBE,DP,GM2ICR,
     x                                     XFBE,E_BE_OMG2)
C======================================================================
      implicit none
      include 'omp_lib.h'

C Input variables
      integer           istart,iend,ng2_seg
      integer           nebf   ! Number contracted electronic basis functions
      integer           npbf   ! Number nuclear basis functions
      integer           ng2
      integer           loop_map(ng2_seg,6)
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM2ICR(ng2)

C Output variables
      double precision  XFBE(nebf,nebf)
      double precision  E_BE_OMG2

C Local variables
      integer           IFIL
      integer           id
      integer           loopi,iLP
      integer           ip,jp
      integer           iec1,jec1
      integer           iec2,jec2
      integer           imap,ia
      integer           ia_12
      double precision  val
      double precision  wtime

C GM2ICR stored as:
C        ip,jp : QM particle indicies
C    iec1,jec1 : regular electron indicies
C    iec2,jec2 : special electron indicies

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
!$ompx private(ia_12)
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

         call index_GAM_2PK(nebf,npbf,ip,jp,iec1,jec1,iec2,jec2,ia_12)
         val=GM2ICR(ia_12)

         XFBE(iec2,jec2)=XFBE(iec2,jec2)+DAE(iec1,jec1)*DP(ip,jp)*val
         E_BE_OMG2=E_BE_OMG2+DP(ip,jp)*DAE(iec1,jec1)*DBE(iec2,jec2)*val

      end do
!$omp end do
!$omp end parallel      

      return
      end

