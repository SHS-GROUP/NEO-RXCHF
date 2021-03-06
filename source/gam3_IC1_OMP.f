!=======================================================================
      subroutine GAM3_IC1(Nchunks,nebf,npebf,npbf,
     x                    ng3,ng3prm,nat,ngtg1,
     x                    pmass,cat,zan,bcoef1,gamma1,
     x                    KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                    ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                    GM3ICR)

!=======================================================================
      implicit none
      include 'omp_lib.h'
! Input Variables
      integer Nchunks
      integer ng3,nebf,npebf,npbf,ng3prm
      integer nat,ngtg1
!-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
      integer KPESTR(nebf)  ! Map contracted index to primitive start
      integer KPEEND(nebf)  ! Map contracted index to primitive end
!-------Basis Set Info-------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)

! Variables Returned
      double precision GM3ICR(ng3)

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
      double precision,allocatable :: XG3ICR(:)
!     double precision,allocatable :: GAM_3(:)
!     double precision,allocatable :: TMPARY(:)

      integer ia
      integer ia_123
      integer ia_132
      integer ia_213
      integer ia_231
      integer ia_312
      integer ia_321

      double precision x123
      double precision x132
      double precision x213
      double precision x231
      double precision x312
      double precision x321
      double precision xxxx

      double precision zero,half,six
      parameter(zero=0.0d+00,half=0.5d+00,six=6.0d+00)

      double precision wtime
      double precision wtime2


         write(*,1000) ng3,nchunks,omp_get_num_procs(),
     xomp_get_max_threads(),1
         wtime = omp_get_wtime()

! Memory to hold un-symmetrized GAM3 integrals
         if(allocated(XG3ICR)) deallocate(XG3ICR)
         allocate( XG3ICR(ng3),stat=istat )

!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      do ichunk=1,Nchunks

         wtime2 = omp_get_wtime()

         call loop_size(1,ng3,Nchunks,ichunk-1,istart,iend)
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

         call thread_gam3_IC1(istart,iend,ng3_seg,ng3,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        loop_map,XG3ICR,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)


         wtime2 = omp_get_wtime() - wtime2
         write(*,2000)ichunk,wtime2

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------)

!-----CLEAN-UP-MEMORY-------------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-MEMORY-------------------------------------------------)

!--------------------SYMMETRIZE----------------------------------------(

      wtime = omp_get_wtime() - wtime
      write(*,3000)wtime
      wtime = omp_get_wtime() 

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

C  GAM_2 Symmetrization:
C  Determine packing indices for XGAM_3 integral matrices

C              As Packed-->       XGAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
!                            read(905,REC=ia_123) x123
                             ! x123=GAM_3(ia_123)
                               x123=XG3ICR(ia_123)

C              As Packed-->       XGAM_3(je2,ie2,je3,ie3,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec3,jec3,iec2,jec2,ia)
                             ia_132=ia
!                            read(905,REC=ia_132) x132
                             ! x132=GAM_3(ia_132)
                               x132=XG3ICR(ia_132)

C              As Packed-->       XGAM_3(je3,ie3,je1,ie1,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,iec3,jec3,ia)
                             ia_213=ia
!                            read(905,REC=ia_213) x213
                             ! x213=GAM_3(ia_213)
                               x213=XG3ICR(ia_213)

C              As Packed-->       XGAM_3(je1,ie1,je3,ie3,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec3,jec3,iec1,jec1,ia)
                             ia_231=ia
!                            read(905,REC=ia_231) x231
                             ! x231=GAM_3(ia_231)
                               x231=XG3ICR(ia_231)

C              As Packed-->       XGAM_3(je2,ie2,je1,ie1,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec3,jec3,iec1,jec1,iec2,jec2,ia)
                             ia_312=ia
!                            read(905,REC=ia_312) x312
                             ! x312=GAM_3(ia_312)
                               x312=XG3ICR(ia_312)

C              As Packed-->       XGAM_3(je1,ie1,je2,ie2,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec3,jec3,iec2,jec2,iec1,jec1,ia)
                             ia_321=ia
!                            read(905,REC=ia_321) x321
                             ! x321=GAM_3(ia_321)
                               x321=XG3ICR(ia_321)

                       xxxx=(x123+x132+x213+x231+x312+x321)/six


c                      write(*,*)'xxxx=',xxxx
cc                     GAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3)=xxxx 
c                      GAM_3PK(ia_123)=xxxx 
c                      call put_GAM3(ia_123,ng3,xxxx)
!                      write(805,REC=ia_123) xxxx
                       GM3ICR(ia_123)=xxxx 

cc                     GAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2)=xxxx
c                      GAM_3PK(ia_132)=xxxx 
c                      call put_GAM3(ia_132,ng3,xxxx)
!                      write(805,REC=ia_132) xxxx
                       GM3ICR(ia_132)=xxxx 

cc                     GAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3)=xxxx
c                      GAM_3PK(ia_213)=xxxx 
c                      call put_GAM3(ia_213,ng3,xxxx)
!                      write(805,REC=ia_213) xxxx
                       GM3ICR(ia_213)=xxxx 

cc                     GAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1)=xxxx
c                      GAM_3PK(ia_231)=xxxx 
c                      call put_GAM3(ia_231,ng3,xxxx)
!                      write(805,REC=ia_231) xxxx
                       GM3ICR(ia_231)=xxxx 

cc                     GAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2)=xxxx
c                      GAM_3PK(ia_312)=xxxx 
c                      call put_GAM3(ia_312,ng3,xxxx)
!                      write(805,REC=ia_312) xxxx
                       GM3ICR(ia_312)=xxxx 

cc                     GAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1)=xxxx
c                      GAM_3PK(ia_321)=xxxx 
c                      call put_GAM3(ia_321,ng3,xxxx)
!                      write(805,REC=ia_321) xxxx
                       GM3ICR(ia_321)=xxxx 


               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      wtime = omp_get_wtime() - wtime

      if(allocated(XG3ICR)) deallocate(XG3ICR)

      write(*,4000) wtime
!--------------------SYMMETRIZE----------------------------------------)


 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 4-PARTICLE INTEGRALS        |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 4-PARTICLE INTEGRALS: ',1X,I12/
     x        8X,'  NUMBER OF BLOCKS (USER DEFINED): ',1X,I12/
     x        8X,'                          ',/,
     x        8X,'  COMPUTATIONAL RESOURCES:',/,
     x        8X,'  ------------------------',/,
     x        8X,'     CORES PER NODE:',1X,I3/
     x        8X,'          AVAILABLE:',1X,I3/
     x        8X,'    NUMBER OF NODES:',1X,I3/)

 2000 FORMAT(8X,'    TIME TO EVALUATE BLOCK ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,'  TIMING SUMMARY FOR 4-PARTICLE INTEGRALS:',/,
     x        8X,'  ----------------------------------------',/,
     x        8X,'    TIME TO EVALUATE ALL INTEGRALS:',1X,F12.4)

 4000 FORMAT(8X,'      TIME TO SYMMETRIZE INTEGRALS:',1X,F12.4/)


      return
      end
C=======================================================================
      subroutine thread_gam3_IC1(istart,iend,ng3_seg,ng3,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           loop_map,XG3ICR,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input Variables
      integer istart,iend,ng3_seg
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas
      integer ng3

      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer iec3,jec3  !
C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
      integer KPESTR(nebf)  ! Map contracted index to primitive start
      integer KPEEND(nebf)  ! Map contracted index to primitive end
C-------Basis Set Info-------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)
      integer loop_map(ng3_seg,8)
      double precision XG3ICR(ng3)
!     double precision GAM_3(ng3_seg)

! Variables Returned

! Local Variables
      integer imap,ia
      double precision OMG3

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!!!!!!!!!!!!$ompx shared(GAM_3)
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(XG3ICR)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(istart,iend)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(ia) 
!$ompx private(OMG3)
!$ompx private(id)

!     id= omp_get_thread_num()
!     write(*,*)' Hello from process ',id
!     if(id.eq.0) then
!        write(*,*)'Threads in use', omp_get_num_threads()
!     end if

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

         call contract_omega3_conv(ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,
     x                             nebf,npebf,npbf,nat,ngtg1,
     x                             pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            OMG3)

!        GAM_3(imap)=OMG3
!        As Packed--> XGAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
!        XGAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) 
         call index_GAM_3PK(nebf,npbf,
     x          ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
         XG3ICR(ia)=OMG3

      end do
!$omp end do
!$omp end parallel      
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

C---OPENMP-TIMING------------------------------------------------------(
!     wtime = omp_get_wtime() - wtime
!     write(*,*)'TIME TO CALCULATE GAM_3 INTEGRALS: ',wtime
!     write(*,*)'ISTART=',istart
!     write(*,*)'IEND  =',iend
C---OPENMP-TIMING------------------------------------------------------)


      return
      end
!=======================================================================
      subroutine EG3IC1(Nchunks,nebf,npbf,ng3,
     x                  DE,DP,GM3ICR,focke,fockp,E_gam3)

!=======================================================================
      implicit none
! Input Variables
      integer Nchunks
      integer ng3,nebf,npbf
      double precision GM3ICR(ng3)
      double precision DE(nebf,nebf),DP(npbf,npbf)
! Variables Returned
      double precision focke(nebf,nebf),fockp(npbf,npbf)
      double precision E_gam3

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
      double precision xfocke(nebf,nebf),xfockp(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_gam3=0.0d+00
      xfocke=0.0d+00
      xfockp=0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      do ichunk=1,Nchunks

         call loop_size(1,ng3,Nchunks,ichunk-1,istart,iend)
!        write(*,*)'after call loop size'
!        write(*,*)'NG3=',ng3
!        write(*,*)'Nchunks=',Nchunks
!        write(*,*)'ichunk=',ichunk
!        write(*,*)'istart=',istart
!        write(*,*)'iend=',iend

! Segment of ng3:
         ng3_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )
!        write(*,*) 'allocate loop_map: ',istat

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

!        write(*,*)'done with nested loop compression'

         call thread_EG3ICR(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                      loop_map,DE,DP,GM3ICR,
     x                      xfocke,xfockp,E_gam3)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!     if(allocated(GAM_3)) deallocate(GAM_3)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
      call add2fock(nebf,xfocke,focke)
      call add2fock(npbf,xfockp,fockp)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end

C======================================================================
      subroutine thread_EG3ICR(istart,iend,ng3_seg,ng3,nebf,npbf,
     x                         loop_map,DE,DP,GM3ICR,
     x                         xfocke,xfockp,E_gam3)

C======================================================================
      implicit none
      include 'omp_lib.h'

C Input Variables
      integer istart,iend,ng3_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng3

      integer loop_map(ng3_seg,8)

      double precision DE(nebf,nebf),DP(npbf,npbf)
      double precision GM3ICR(ng3)

! Variables Returned
      double precision xfocke(nebf,nebf),xfockp(npbf,npbf)
      double precision E_gam3

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
      double precision GAM3
      double precision OMG3
      double precision OMG3_123       
      double precision OMG3_213
      double precision OMG3_312
      double precision OMG3_132
      double precision OMG3_231
      double precision OMG3_321


!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!---OPENMP-TIMING------------------------------------------------------(
!     wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(nebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DE,DP)
!$ompx shared(GM3ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(OMG3)
!$ompx private(GAM3)
!$ompx private(id)
!$ompx private(OMG3_123,OMG3_213,OMG3_312,OMG3_132,OMG3_231,OMG3_321)
!$ompx private(ia_123,ia_213,ia_312,ia_132,ia_231,ia_321)
!$ompx private(ia)
!$ompx reduction(+:E_gam3)
!$ompx reduction(+:xfocke)
!$ompx reduction(+:xfockp)

      id= omp_get_thread_num()
!     write(*,*)' Hello from process ',id
!     if(id.eq.0) then
!        write(*,*)'Threads in use', omp_get_num_threads()
!     end if

!$omp do
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

!--------------------------FOR-IN-CORE-DSCF----------------------------(
!        call contract_omega3(
!                           ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,
!    x                        OMG3_123)

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
                             OMG3_123=GM3ICR(ia_123)

!        call contract_omega3(
!                           ip,jp,iec1,jec2,iec2,jec1,iec3,jec3,ia)
!    x                        OMG3_213)

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec2,iec2,jec1,iec3,jec3,ia)
                             ia_213=ia
                             OMG3_213=GM3ICR(ia_213)


!        call contract_omega3(
!                           ip,jp,iec1,jec3,iec2,jec1,iec3,jec2,ia)
!    x                        OMG3_312)

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec3,iec2,jec1,iec3,jec2,ia)
                             ia_312=ia
                             OMG3_312=GM3ICR(ia_312)


!        call contract_omega3(
!                           ip,jp,iec1,jec1,iec2,jec3,iec3,jec2,ia)
!    x                        OMG3_132)

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec3,iec3,jec2,ia)
                             ia_132=ia
                             OMG3_132=GM3ICR(ia_132)


!        call contract_omega3(
!                           ip,jp,iec1,jec2,iec2,jec3,iec3,jec1,ia)
!    x                        OMG3_231)

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec2,iec2,jec3,iec3,jec1,ia)
                             ia_231=ia
                             OMG3_231=GM3ICR(ia_231)


!        call contract_omega3(
!                           ip,jp,iec1,jec3,iec2,jec2,iec3,jec1,ia)
!    x                        OMG3_321)

         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec3,iec2,jec2,iec3,jec1,ia)
                             ia_321=ia
                             OMG3_321=GM3ICR(ia_321)



         GAM3=8.0d+00*OMG3_123
     x       -4.0d+00*OMG3_213
     x       +2.0d+00*OMG3_312
     x       -4.0d+00*OMG3_132
     x       +2.0d+00*OMG3_231
     x       -4.0d+00*OMG3_321

         GAM3=GAM3/8.0d+00

            xfocke(iec1,jec1)=xfocke(iec1,jec1)+
     *             DP(ip,jp)*DE(iec2,jec2)*DE(iec3,jec3)*3.0d+00*GAM3

            xfockp(ip,jp)=xfockp(ip,jp)+
     *             DE(iec1,jec1)*DE(iec2,jec2)*DE(iec3,jec3)*GAM3

            E_gam3=E_gam3+
     x            DE(iec1,jec1)*DE(iec2,jec2)*DE(iec3,jec3)*
     x            DP(ip,jp)*GAM3

!--------------------------FOR-IN-CORE-DSCF----------------------------)

         end do
!$omp end do
!$omp end parallel      
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
