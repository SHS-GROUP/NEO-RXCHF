!=======================================================================
      subroutine RXCHFmult_GAM2_IC1(Nchunks,
     x                            nebf1,npebf1,nebf2,npebf2,npbf,
     x                            ng2,ng2s,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR1,KPEEND1,
     x                            AMPEB2C1,AGEBFCC1,
     x                            ELCEX1,ELCAM1,ELCBFC1,
     x                            KPESTR2,KPEEND2,
     x                            AMPEB2C2,AGEBFCC2,
     x                            ELCEX2,ELCAM2,ELCBFC2,
     x                            AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                            GM2ICR,GM2sICR)

!=======================================================================
      implicit none
      include 'omp_lib.h'
! Input Variables
      integer Nchunks
      integer ng2,ng2s,ng2prm,npbf
      integer nebf1,npebf1
      integer nebf2,npebf2
      integer nat,ngtg1
!-------Basis Set Info-------(
      integer ELCAM1(npebf,3)  ! Angular mom for electrons
      integer ELCAM2(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX1(npebf) ! Exponents: elec basis
      double precision ELCEX2(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC1(npebf,3) ! Basis centers: elec basis
      double precision ELCBFC2(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C1(npebf) ! Map primitive index to contracted
      double precision AGEBFCC1(npebf) ! Map prim index to contract coef
      integer AMPEB2C2(npebf) ! Map primitive index to contracted
      double precision AGEBFCC2(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
      integer KPESTR1(nebf)  ! Map contracted index to primitive start
      integer KPEEND1(nebf)  ! Map contracted index to primitive end
      integer KPESTR2(nebf)  ! Map contracted index to primitive start
      integer KPEEND2(nebf)  ! Map contracted index to primitive end
!-------Basis Set Info-------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)

! Variables Returned
      double precision GM2ICR(ng2)
      double precision GM2sICR(ng2s)

! Local Variables
      integer istat,ichunk,istart,iend,ng2_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2
      integer,allocatable :: loop_map(:,:)

      double precision,allocatable :: XG2ICR(:),XG2sICR(:)

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

      if(allocated(XG2ICR)) deallocate(XG2ICR)
      allocate( XG2ICR(ng2),stat=istat )
      if(allocated(XG2sICR)) deallocate(XG2sICR)
      allocate( XG2sICR(ng2s),stat=istat )

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
               do iec1=1,nebf1
               do jec1=1,nebf1
                  do iec2=1,nebf2
                  do jec2=1,nebf2

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

         call RXCHFmult_thread_gam2_IC1(istart,iend,ng2_seg,ng2,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           loop_map,XG2_1ICR,XG2_2ICR,XG2sICR,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

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
            do iec1=1,nebf1
            do jec1=1,nebf1
               do iec2=1,nebf2
               do jec2=1,nebf2

!  GAM_2 Symmetrization:
!  Determine packing indices for XGAM_2 integral matrices

!              As Packed-->       XGAM_2(je2,ie2,je1,ie1,jp,ip)
!                           XGAM_2(ip,jp,ie1,je1,ie2,je2,) 
C ARS( particle 1: special e ; particle 2: regular e ; index 3: prot )
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia_12)
!              As Packed-->       XGAM_2(je1,ie1,je2,ie2,jp,ip)
!                           XGAM_2(ip,jp,ie2,je2,ie1,je1) 
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,ia_21)

C RXCHFmult(
C    index 1: regular electron
C    index 2: special electron
C    index 3: proton
C )

C Symmetrized integrals in GM2_1ICR (XCHF integrals)
                       x12=XG2_1ICR(ia_12)
                       x21=XG2_1ICR(ia_21)
                       GM2_1ICR(ia_12)=(x12+x21)/two

C Unsymmetrized integrals in GM2_2ICR (interaction integrals)
                       x12=XG2_2ICR(ia_12)
                       GM2_2ICR(ia_12)=x12

C Symmetrized integrals in GM2sICR (XCHF integrals)
                       x12=XG2sICR(ia_12)
                       x21=XG2sICR(ia_21)
                       GM2sICR(ia_12)=(x12+x21)/two

               end do
               end do
            end do
            end do
         end do
         end do

         wtime = omp_get_wtime() - wtime

         if(allocated(XG2sICR)) deallocate(XG2sICR)
         if(allocated(XG2ICR)) deallocate(XG2ICR)

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
      subroutine RXCHFmult_thread_gam2_IC1(istart,iend,ng2_seg,ng2,
     x                            nebf,npebf,npbf,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            loop_map,XG2_1ICR,XG2_2ICR,XG2sICR,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

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
      double precision XG2_1ICR(ng2),XG2_2ICR(ng2),XG2sICR(ng2)

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap,ia
      double precision OMG2_1,OMG2_2,OMG2s

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
!$ompx shared(XG2_1ICR,XG2_2ICR,XG2sICR)
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
!$ompx private(OMG2_1,OMG2_2,OMG2s)
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

         call RXCHFmult_contract_omega2_conv(ip,jp,iec1,jec1,iec2,jec2,
     x                            nebf,npebf,npbf,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            OMG2_1,OMG2_2,OMG2s)


!         XG2ICR(imap)=OMG2

         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia)

         XG2_1ICR(ia)=OMG2_1
         XG2_2ICR(ia)=OMG2_2
         XG2sICR(ia)=OMG2s

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
!======================================================================
      subroutine RXCHFmult_contract_omega2_conv(ip,jp,
     x                           iec1,jec1,iec2,jec2,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                           OMG2_1,OMG2_2,OMG2s)

!======================================================================
      implicit none
c     include 'omp_lib.h'
c     include 'mpif.h'

C Input Variables
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas

      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
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

C Variables Returned
      double precision OMG2_1,OMG2_2,OMG2s

C Local Variables
      integer ie1,je1
      integer ie2,je2
      integer ie1_start
      integer ie2_start
      integer je1_start
      integer je2_start
      integer ie1_end
      integer ie2_end
      integer je1_end
      integer je2_end

      double precision Cof_ie1,Cof_je1
      double precision Cof_ie2,Cof_je2
      double precision Cof_ip,Cof_jp
C--------------------------------(
C Basis set-related local variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3

      double precision A1,Amat1(3) 
      double precision A2,Amat2(3) 
      double precision A3,Amat3(3) 
      double precision B1,Bmat1(3) 
      double precision B2,Bmat2(3) 
      double precision B3,Bmat3(3) 
C--------------------------------)
      double precision ans
      double precision ansE1,ansE2,ansS


      ie1_start=KPESTR(iec1)
      ie2_start=KPESTR(iec2)

      je1_start=KPESTR(jec1)
      je2_start=KPESTR(jec2)

      ie1_end=KPEEND(iec1)
      ie2_end=KPEEND(iec2)

      je1_end=KPEEND(jec1)
      je2_end=KPEEND(jec2)

      OMG2_1=0.0d+00
      OMG2_2=0.0d+00
      OMG2s=0.0d+00

      do ie1=ie1_start,ie1_end
       do je1=je1_start,je1_end
        do ie2=ie2_start,ie2_end
         do je2=je2_start,je2_end
               
C             Get Basis set info:
              A1=ELCEX(ie1)
              I1=ELCAM(ie1,1)
              J1=ELCAM(ie1,2)
              K1=ELCAM(ie1,3)
              Amat1(1)=ELCBFC(ie1,1)
              Amat1(2)=ELCBFC(ie1,2)
              Amat1(3)=ELCBFC(ie1,3)

              A2=ELCEX(ie2)
              I2=ELCAM(ie2,1)
              J2=ELCAM(ie2,2)
              K2=ELCAM(ie2,3)
              Amat2(1)=ELCBFC(ie2,1)
              Amat2(2)=ELCBFC(ie2,2)
              Amat2(3)=ELCBFC(ie2,3)

              A3=NUCEX(ip)
              I3=NUCAM(ip,1)
              J3=NUCAM(ip,2)
              K3=NUCAM(ip,3)
              Amat3(1)=NUCBFC(ip,1)
              Amat3(2)=NUCBFC(ip,2)
              Amat3(3)=NUCBFC(ip,3)

              B1=ELCEX(je1)
              L1=ELCAM(je1,1)
              M1=ELCAM(je1,2)
              N1=ELCAM(je1,3)
              Bmat1(1)=ELCBFC(je1,1)
              Bmat1(2)=ELCBFC(je1,2)
              Bmat1(3)=ELCBFC(je1,3)

              B2=ELCEX(je2)
              L2=ELCAM(je2,1)
              M2=ELCAM(je2,2)
              N2=ELCAM(je2,3)
              Bmat2(1)=ELCBFC(je2,1)
              Bmat2(2)=ELCBFC(je2,2)
              Bmat2(3)=ELCBFC(je2,3)

              B3=NUCEX(jp)
              L3=NUCAM(jp,1)
              M3=NUCAM(jp,2)
              N3=NUCAM(jp,3)
              Bmat3(1)=NUCBFC(jp,1)
              Bmat3(2)=NUCBFC(jp,2)
              Bmat3(3)=NUCBFC(jp,3)

C  Get primitive Electron Basis Function Contraction Coefficients 
              Cof_ie1=AGEBFCC(ie1)
              Cof_ie2=AGEBFCC(ie2)
              Cof_je1=AGEBFCC(je1)
              Cof_je2=AGEBFCC(je2)
C  Get Nuclear Basis Function Contraction Coefficients
              Cof_ip=AGNBFCC(ip)
              Cof_jp=AGNBFCC(jp)

C ARS( particle 1: special e ; particle 2: regular e ; index 3: prot )
C---------------------OMG_12-------------------------------------------(
              call RXCHFmult_xcalc_GAM2_MD(I1,J1,K1,A1,Amat1,
     x                                     I2,J2,K2,A2,Amat2,
     x                                     I3,J3,K3,A3,Amat3,
     x                                     L1,M1,N1,B1,Bmat1,
     x                                     L2,M2,N2,B2,Bmat2,
     x                                     L3,M3,N3,B3,Bmat3,
     x                                     nat,ngtg1,
     x                                     pmass,cat,zan,
     x                                     bcoef1,gamma1,
     x                                     ansE1,ansE2,ansS)

!                       call underflow(ans)

                        OMG2_1=OMG2_1+ansE1
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2

                        OMG2_2=OMG2_2+ansE2
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2

                        OMG2s=OMG2s+ansS
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2
C---------------------OMG_12-------------------------------------------)

         end do
        end do
       end do
      end do


      return
      end

