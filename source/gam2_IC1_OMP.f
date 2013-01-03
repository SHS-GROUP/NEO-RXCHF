!=======================================================================
      subroutine GAM2_IC1(Nchunks,nebf,npebf,npbf,
     x                    ng2,ng2prm,nat,ngtg1,
     x                    pmass,cat,zan,bcoef1,gamma1,
     x                    KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                    ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                    GM2ICR,GM2sICR)

!=======================================================================
      implicit none
      include 'omp_lib.h'
! Input Variables
      integer Nchunks
      integer ng2,nebf,npebf,npbf,ng2prm
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
      double precision GM2ICR(ng2)
      double precision GM2sICR(ng2)

! Local Variables
      integer istat,ichunk,istart,iend,ng2_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2
      integer,allocatable :: loop_map(:,:)

      double precision,allocatable :: XG2ICR(:)
      double precision,allocatable :: XG2sICR(:)

!     integer ia
      integer ia_12
      integer ia_21

      double precision x12,y12
      double precision x21,y21
      double precision xx
      double precision yy

      double precision zero,half,two
      parameter(zero=0.0d+00,half=0.5d+00,two=2.0d+00)

      double precision wtime
      double precision wtime2

         write(*,1000) ng2,nchunks,omp_get_num_procs(),
     xomp_get_max_threads(),1
         wtime = omp_get_wtime()

! Memory to hold un-symmetrized GAM2 integrals
      if(allocated(XG2ICR)) deallocate(XG2ICR)
      allocate( XG2ICR(ng2),stat=istat )

      if(allocated(XG2sICR)) deallocate(XG2sICR)
      allocate( XG2sICR(ng2),stat=istat )

!-----CHOP-UP-THE-CALCULATION-OF-GAM_2--------------------------------(
      do ichunk=1,Nchunks

            wtime2 = omp_get_wtime()

            call loop_size(1,ng2,Nchunks,ichunk-1,istart,iend)
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

! Determine the size of the subloop to be computed on this node:
!        call loop_size(1,ng2_seg,nodes,myid,my1st,mylast)

!        call thread_gam2_conv(my1st,mylast,ng2_seg,ng2,
         call thread_gam2_IC1(istart,iend,ng2_seg,ng2,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        loop_map,XG2ICR,XG2sICR,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

            wtime2 = omp_get_wtime() - wtime2
            write(*,2000)ichunk,wtime2

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_2--------------------------------)

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

!  GAM_2 Symmetrization:
!  Determine packing indices for XGAM_2 integral matrices

!              As Packed-->       XGAM_2(je2,ie2,je1,ie1,jp,ip)
!                           XGAM_2(ip,jp,ie1,je1,ie2,je2,) 
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia_12)

                             x12=XG2ICR(ia_12)
                             y12=XG2sICR(ia_12)


!              As Packed-->       XGAM_2(je1,ie1,je2,ie2,jp,ip)
!                           XGAM_2(ip,jp,ie2,je2,ie1,je1) 
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,ia_21)

                             x21=XG2ICR(ia_21)
                             y21=XG2sICR(ia_21)


                       GM2ICR(ia_12) =(x12+x21)/two
                       GM2sICR(ia_12)=(y12+y21)/two

               end do
               end do
            end do
            end do
         end do
         end do

         wtime = omp_get_wtime() - wtime

         if(allocated(XG2ICR)) deallocate(XG2ICR)
         if(allocated(XG2sICR)) deallocate(XG2sICR)

         write(*,4000) wtime
!--------------------SYMMETRIZE----------------------------------------)


 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 3-PARTICLE INTEGRALS        |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 3-PARTICLE INTEGRALS: ',1X,I12/
     x        8X,'  NUMBER OF BLOCKS (USER DEFINED): ',1X,I12/
     x        8X,'                          ',/,
     x        8X,'  COMPUTATIONAL RESOURCES:',/,
     x        8X,'  ------------------------',/,
     x        8X,'     CORES PER NODE:',1X,I3/
     x        8X,'          AVAILABLE:',1X,I3/
     x        8X,'    NUMBER OF NODES:',1X,I3/)

 2000 FORMAT(8X,'    TIME TO EVALUATE BLOCK ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,'  TIMING SUMMARY FOR 3-PARTICLE INTEGRALS:',/,
     x        8X,'  ----------------------------------------',/,
     x        8X,'    TIME TO EVALUATE ALL INTEGRALS:',1X,F12.4)

 4000 FORMAT(8X,'      TIME TO SYMMETRIZE INTEGRALS:',1X,F12.4/)


      return
      end
!=======================================================================
      subroutine thread_gam2_IC1(istart,iend,ng2_seg,ng2,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           loop_map,XG2ICR,XG2sICR,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

!=======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng2_seg
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas
      integer ng2

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
      integer loop_map(ng2_seg,6)

! Variables Returned
      double precision XG2ICR(ng2)
      double precision XG2sICR(ng2)

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap,ia
      double precision OMG2
      double precision OMG2s

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)


!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(istart,iend)
!$ompx shared(loop_map)
!$ompx shared(XG2ICR)
!$ompx shared(XG2sICR)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(OMG2)
!$ompx private(OMG2s)
!$ompx private(ia)
!$ompx private(id)

!     id= omp_get_thread_num()
!     write(*,*)' Hello from process ',id
!     if(id.eq.0) then
!        write(*,*)'Threads in use', omp_get_num_threads()
!     end if

!$omp do SCHEDULE(RUNTIME)
      do iLP=istart,iend

         imap=iLp-istart+1
!        imap=iLp
         jec2=loop_map(imap,1)
         iec2=loop_map(imap,2)
         jec1=loop_map(imap,3)
         iec1=loop_map(imap,4)
         jp =loop_map(imap,5)
         ip =loop_map(imap,6)

         call contract_omega2_conv(ip,jp,iec1,jec1,iec2,jec2,
     x                             nebf,npebf,npbf,nat,ngtg1,
     x                             pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            OMG2,OMG2s)


!        XG2ICR(imap)=OMG2
!        XG2sICR(imap)=OMG2s

         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia)

         XG2ICR(ia)=OMG2
         XG2sICR(ia)=OMG2s

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!=======================================================================
      subroutine EG2IC1(Nchunks,nebf,npbf,ng2,
     x                  DE,DP,GM2ICR,GM2SICR,
     x                  FE,SE,FP,SP,E_gam2,S_gam2)

! Calculate Fock matrices
! and contribution to total energy for 3-particle terms.
! AO contracted integrals are stored in-core.
!=======================================================================
      implicit none

! Input Variables
      integer Nchunks
      integer ng2,nebf,npbf
      double precision GM2ICR(ng2)
      double precision GM2SICR(ng2)
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FE(nebf,nebf)
      double precision FP(npbf,npbf)
      double precision SE(nebf,nebf)
      double precision SP(npbf,npbf)
      double precision E_gam2
      double precision S_gam2

! Local Variables
      integer istat,ichunk,istart,iend,ng2_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2
      integer,allocatable :: loop_map(:,:)
      double precision XFE(nebf,nebf)
!     double precision XSE(nebf,nebf)
      double precision XFP(npbf,npbf)
!     double precision XSP(npbf,npbf)

!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_gam2=0.0d+00
      S_gam2=0.0d+00
      XFE=0.0d+00
      SE=0.0d+00
      XFP =0.0d+00
      SP =0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_2--------------------------------(
      do ichunk=1,Nchunks

         call loop_size(1,ng2,Nchunks,ichunk-1,istart,iend)

! Segment of ng2:
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

         call thread_EG2ICR(istart,iend,ng2_seg,ng2,nebf,npbf,
     x                      loop_map,DE,DP,GM2ICR,GM2SICR,
     x                      XFE,SE,XFP,SP,E_gam2,S_gam2)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
      call add2fock(nebf,XFE,FE)
!     call add2fock(nebf,XSE,SE)
      call add2fock(npbf,XFP,FP)
!     call add2fock(npbf,XSP,SP)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!======================================================================
      subroutine thread_EG2ICR(istart,iend,ng2_seg,ng2,nebf,npbf,
     x                         loop_map,DE,DP,GM2ICR,GM2SICR,
     x                         XFE,SE,XFP,SP,E_gam2,S_gam2)

!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng2_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng2

      integer loop_map(ng2_seg,6)

      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM2ICR(ng2)
      double precision GM2SICR(ng2)

! Variables Returned
      double precision XFE(nebf,nebf)
      double precision SE(nebf,nebf)
      double precision XFP(npbf,npbf)
      double precision SP(npbf,npbf)
      double precision E_gam2
      double precision S_gam2

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap
      integer ia_12
      integer ia_21

      double precision OMG12,OMG21
      double precision SOMG12,SOMG21
      double precision val_gam2,val_gam2s
      double precision half,two
      parameter(half=0.5d+00,two=2.0d+00)

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
!$ompx shared(DE,DP)
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
!$ompx private(val_gam2,val_gam2s)
!$ompx reduction(+:XFE)
!$ompx reduction(+:SE)
!$ompx reduction(+:XFP)
!$ompx reduction(+:SP)
!$ompx reduction(+:E_gam2)
!$ompx reduction(+:S_gam2)

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
         val_gam2=OMG12-half*OMG21

         SOMG12=GM2SICR(ia_12)
         SOMG21=GM2SICR(ia_21)
         val_gam2s=SOMG12-half*SOMG21

!----------------Form-Elec-Fock-Matrix-------------------(
         XFE(iec1,jec1)=XFE(iec1,jec1)
     x                  + DP(ip,jp)*DE(iec2,jec2)*two*val_gam2

! Multiplication by 2 is taken care of in XCfock_correction
! This is a sloppy way of doing this and should be changed
!        XSE(iec1,jec1)=XSE(iec1,jec1)
!    x                  + DP(ip,jp)*DE(iec2,iec2)*two*SOMG12
         SE(iec1,jec1)=SE(iec1,jec1)
     x                 + DP(ip,jp)*DE(iec2,jec2)*val_gam2s
!----------------Form-Elec-Fock-Matrix-------------------)

!-------------------Form-QM-Partice-Fock-Matrix-------------------------(
       XFP(ip,jp)=XFP(ip,jp)+DE(iec1,jec1)*DE(iec2,jec2)*val_gam2
       SP(ip,jp)=SP(ip,jp)+DE(iec1,jec1)*DE(iec2,jec2)*val_gam2s
!-------------------Form-QM-Partice-Fock-Matrix-------------------------)

         E_gam2=E_gam2+DE(iec1,jec1)*DE(iec2,jec2)*DP(ip,jp)*val_gam2
         S_gam2=S_gam2+DE(iec1,jec1)*DE(iec2,jec2)*DP(ip,jp)*val_gam2s

         end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end

