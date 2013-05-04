!=======================================================================
      subroutine RXCHFne_GAM2_IC1(Nchunks,nebf,npebf,npbf,
     x                    ng2,ng2prm,nat,ngtg1,
     x                    pmass,cat,zan,bcoef1,gamma1,
     x                    KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                    ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                    GM2ICR,GM2exICR)

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
      double precision GM2ICR(ng2)  ! Regular OMG2 terms
      double precision GM2exICR(ng2) ! Exchange OMG2 terms

! Local Variables
      integer istat,ichunk,istart,iend,ng2_seg
      integer my1st,mylast
      integer iLp,imap
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2
      integer,allocatable :: loop_map(:,:)

      double precision,allocatable :: XG2ICR(:),XG2exICR(:)

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
C ARS( gam2 testing
C      integer ARSiLP,ARSimap
C      character*3 chunkstr
C )

         write(*,1000) ng2,nchunks,omp_get_num_procs(),
     xomp_get_max_threads(),1
         wtime = omp_get_wtime()

! Memory to hold un-symmetrized GAM2 integrals
      if(allocated(XG2ICR)) deallocate(XG2ICR)
      allocate( XG2ICR(ng2),stat=istat )
      if(allocated(XG2exICR)) deallocate(XG2exICR)
      allocate( XG2exICR(ng2),stat=istat )

!-----CHOP-UP-THE-CALCULATION-OF-GAM_2--------------------------------(
      do ichunk=1,Nchunks

            wtime2 = omp_get_wtime()

            call loop_size(1,ng2,Nchunks,ichunk-1,istart,iend)
            ng2_seg=1+iend-istart

C ARS( gam2 testing
C        write(*,*) "ichunk,istart,iend:"
C        write(*,*) ichunk,istart,iend
C )
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
C ARS( gam2 testing
C        if (ichunk.eq.1) then
C         write(*,*) "start output"
C         open(unit=97,file="params.fmt")
C         write(97,*) loopi,istart,iend,ng2_seg,ng2
C         write(97,*) nebf,npebf,npbf,nat,ngtg1,
C     x             pmass,cat,zan,bcoef1,gamma1,
C     x             KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
C     x             ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC
C         write(97,*) loop_map
C         close(97)
C         write(*,*) "end output"
C        end if
C )
         call RXCHFne_thread_gam2_IC1(istart,iend,ng2_seg,ng2,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        loop_map,XG2ICR,XG2exICR,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

            wtime2 = omp_get_wtime() - wtime2
            write(*,2000)ichunk,wtime2
C ARS( gam2 testing
C        write(chunkstr,'(I3.3)') ichunk
C        open(unit=98,file="gam2-ichunk"//chunkstr//".ufm",
C     x       form="unformatted",access="direct",recl=8)
C        do ARSiLP=istart,iend
C          ARSimap=ARSiLp-istart+1
C          write(98,rec=ARSimap) XG2ICR(ARSimap)
C        end do
C        close(98)
C )

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
C ARS( particle 1: special e ; particle 2: regular e ; index 3: prot )
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia_12)
         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec2,jec2,iec1,jec1,ia_21)


C Store unsymmetrized integrals in GM2ICR
                       x12=XG2ICR(ia_12)
                       GM2ICR(ia_12)=x12


!              As Packed-->       XGAM_2(je1,ie1,je2,ie2,jp,ip)
!                           XGAM_2(ip,jp,ie2,je2,ie1,je1) 

C Store symmetrized integrals in GM2ICR
                       x12=XG2exICR(ia_12)
                       x21=XG2exICR(ia_21)
                       GM2exICR(ia_12)=(x12+x21)/two

               end do
               end do
            end do
            end do
         end do
         end do

         wtime = omp_get_wtime() - wtime

         if(allocated(XG2ICR)) deallocate(XG2ICR)
         if(allocated(XG2exICR)) deallocate(XG2exICR)

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
      subroutine RXCHFne_thread_gam2_IC1(istart,iend,ng2_seg,ng2,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           loop_map,XG2ICR,XG2exICR,
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
      double precision XG2ICR(ng2),XG2exICR(ng2)

! Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap,ia
      double precision OMG2,OMG2ex

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
!$ompx shared(XG2ICR,XG2exICR)
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
!$ompx private(OMG2,OMG2ex)
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

         call RXCHFne_contract_omega2_conv(ip,jp,iec1,jec1,iec2,jec2,
     x                             nebf,npebf,npbf,nat,ngtg1,
     x                             pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            OMG2,OMG2ex)


!         XG2ICR(imap)=OMG2

         call index_GAM_2PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,ia)

         XG2ICR(ia)=OMG2
         XG2exICR(ia)=OMG2ex

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end
