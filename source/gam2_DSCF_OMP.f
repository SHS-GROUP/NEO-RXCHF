!=======================================================================
      subroutine GAM2_DSCF(Nchunks,nebf,npebf,npbf,
     x                     ng2,ng2prm,nat,ngtg1,
     x                     pmass,cat,zan,bcoef1,gamma1,
     x                     KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                     ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                     DE,DP,focke,fockp,se2,sp2,E_gam2,S_gam2)

!=======================================================================
      implicit none
! Input Variables
      integer Nchunks
      integer ng2,nebf,npebf,npbf,ng2prm
      integer nat,ngtg1
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

      double precision DE(nebf,nebf),DP(npbf,npbf)
! Variables Returned
      double precision focke(nebf,nebf),fockp(npbf,npbf)
      double precision se2(nebf,nebf),sp2(npbf,npbf)
      double precision E_gam2
      double precision S_gam2

! Local Variables
      integer istat,ichunk,istart,iend,ng2_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2
      integer,allocatable :: loop_map(:,:)
!     double precision,allocatable :: GAM_3(:)
      double precision xfocke(nebf,nebf),xfockp(npbf,npbf)
      double precision xse2(nebf,nebf),xsp2(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_gam2=0.0d+00
      xfocke=0.0d+00
      xfockp=0.0d+00
      S_gam2=0.0d+00
      xse2=0.0d+00
      xsp2=0.0d+00
!     call zero_out(nebf,xfocke)
!     call zero_out(npbf,xfockp)
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_2--------------------------------(
      do ichunk=1,Nchunks

         call loop_size(1,ng2,Nchunks,ichunk-1,istart,iend)
!        write(*,*)'after call loop size'
!        write(*,*)'NG2=',ng2
!        write(*,*)'Nchunks=',Nchunks
!        write(*,*)'ichunk=',ichunk
!        write(*,*)'istart=',istart
!        write(*,*)'iend=',iend

! Segment of ng2:
         ng2_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng2_seg,6),stat=istat )
!        write(*,*) 'allocate loop_map: ',istat

!        if(allocated(GAM_2)) deallocate(GAM_2)
!        allocate( GAM_2(ng2_seg),stat=istat )
!        write(*,*) 'allocate GAM_2: ',istat

! Nested loop compression for this chunk:
! There should be a better way to do this
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

!        write(*,*)'done with nested loop compression'

         call thread_gam2(istart,iend,ng2_seg,ng2,
     x                    nebf,npebf,npbf,nat,ngtg1,
     x                    pmass,cat,zan,bcoef1,gamma1,loop_map,
     x                    KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                    ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                    DE,DP,xfocke,xfockp,E_gam2,
     x                    xse2,xsp2,S_gam2)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!     if(allocated(GAM_2)) deallocate(GAM_2)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
      call add2fock(nebf,xfocke,focke)
      call add2fock(npbf,xfockp,fockp)
      call add2fock(nebf,xse2,se2)
      call add2fock(npbf,xsp2,sp2)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end

C======================================================================
      subroutine thread_gam2(istart,iend,ng2_seg,ng2,
     x                       nebf,npebf,npbf,nat,ngtg1,
     x                       pmass,cat,zan,bcoef1,gamma1,loop_map,
     x                       KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                       DE,DP,xfocke,xfockp,E_gam2,
     x                       xse2,xsp2,S_gam2)

C======================================================================
      implicit none
      include 'omp_lib.h'

C Input Variables
      integer istart,iend,ng2_seg
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas
      integer ng2

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
      integer loop_map(ng2_seg,6)
      double precision DE(nebf,nebf),DP(npbf,npbf)
C Variables Returned
!     double precision GAM_2(ng2_seg)
      double precision xfocke(nebf,nebf),xfockp(npbf,npbf)
      double precision xse2(nebf,nebf),xsp2(npbf,npbf)
      double precision E_gam2
      double precision S_gam2

C Local Variables
      integer ip,jp
      integer iec1,jec1  !
      integer iec2,jec2  ! Contracted elec basis function indices
      integer imap,ia
      double precision GAM2
      double precision GAM2s
      double precision OMG2
      double precision OMG2s
      double precision OMG2_12       
      double precision OMG2_21
      double precision OMG2s_12       
      double precision OMG2s_21
      double precision half,two

!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

      half=5.0d-01
      two=2.0d+00

!---OPENMP-TIMING------------------------------------------------------(
!     wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf,ng2_seg)
!$ompx shared(ng2)
!$ompx shared(half,two)
!$ompx shared(DE,DP)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(OMG2)
!$ompx private(OMG2s)
!$ompx private(GAM2)
!$ompx private(GAM2s)
!$ompx private(id)
!$ompx private(OMG2_12,OMG2_21)
!$ompx private(OMG2s_12,OMG2s_21)
!$ompx reduction(+:E_gam2)
!$ompx reduction(+:S_gam2)
!$ompx reduction(+:xfocke)
!$ompx reduction(+:xfockp)
!$ompx reduction(+:xse2)
!$ompx reduction(+:xsp2)
!XXXX!$ompx shared(GAM_2)

      id= omp_get_thread_num()
!     write(*,*)' Hello from process ',id
!     if(id.eq.0) then
!        write(*,*)'Threads in use', omp_get_num_threads()
!     end if

!$omp do
      do iLP=istart,iend

         imap=iLp-istart+1
         jec2=loop_map(imap,1)
         iec2=loop_map(imap,2)
         jec1=loop_map(imap,3)
         iec1=loop_map(imap,4)
         jp =loop_map(imap,5)
         ip =loop_map(imap,6)

C--------------------------FOR-DSCF------------------------------------(
         call contract_omega2(ip,jp,iec1,jec1,iec2,jec2,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        OMG2_12,OMG2s_12)

         call contract_omega2(ip,jp,iec1,jec2,iec2,jec1,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        OMG2_21,OMG2s_21)



         GAM2=OMG2_12-half*OMG2_21

         xfocke(iec1,jec1)=xfocke(iec1,jec1)+
     *          DP(ip,jp)*DE(iec2,jec2)*two*GAM2

         xfockp(ip,jp)=xfockp(ip,jp)+
     *          DE(iec1,jec1)*DE(iec2,jec2)*GAM2

         E_gam2=E_gam2+DE(iec1,jec1)*DE(iec2,jec2)*DP(ip,jp)*GAM2


         GAM2s=OMG2s_12-half*OMG2s_21

         xse2(iec1,jec1)=xse2(iec1,jec1)+
     *                DP(ip,jp)*DE(iec2,jec2)*GAM2s

         xsp2(ip,jp)=xsp2(ip,jp)+
     *                DE(iec1,jec1)*DE(iec2,jec2)*GAM2s

         S_gam2=S_gam2+DE(iec1,jec1)*DE(iec2,jec2)*DP(ip,jp)*GAM2s

C--------------------------FOR-DSCF------------------------------------)

         end do
!$omp end do
!$omp end parallel      
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

C---OPENMP-TIMING------------------------------------------------------(
!     wtime = omp_get_wtime() - wtime
!     write(*,*)'TIME TO CALCULATE GAM_3 INTEGRALS: ',wtime
C---OPENMP-TIMING------------------------------------------------------)

CCWS-IO
!        write(*,*)
!        write(*,*)'FINISHED CALCULATING CONTRACTED GAM_2 INTEGRALS:'
!        write(*,*)'ISTART=',istart
!        write(*,*)'IEND  =',iend
!        write(*,*)'WRITING THEM TO DISK...'
!        write(*,*)


      return
      end

C======================================================================
      subroutine contract_omega2(ip,jp,iec1,jec1,iec2,jec2,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                           OMG2,OMG2s)

C======================================================================
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
      double precision OMG2
      double precision OMG2s

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
      double precision ansE,ansS
      double precision OMG2_12
      double precision OMG2_21
      double precision OMG2s_12
      double precision OMG2s_21
      double precision half
      parameter(half=5.0d-01)


      ie1_start=KPESTR(iec1)
      ie2_start=KPESTR(iec2)

      je1_start=KPESTR(jec1)
      je2_start=KPESTR(jec2)

      ie1_end=KPEEND(iec1)
      ie2_end=KPEEND(iec2)

      je1_end=KPEEND(jec1)
      je2_end=KPEEND(jec2)

      OMG2_12=0.0d+00
      OMG2_21=0.0d+00

      OMG2s_12=0.0d+00
      OMG2s_21=0.0d+00

      A3=NUCEX(ip)
      I3=NUCAM(ip,1)
      J3=NUCAM(ip,2)
      K3=NUCAM(ip,3)
      Amat3(1)=NUCBFC(ip,1)
      Amat3(2)=NUCBFC(ip,2)
      Amat3(3)=NUCBFC(ip,3)

      B3=NUCEX(jp)
      L3=NUCAM(jp,1)
      M3=NUCAM(jp,2)
      N3=NUCAM(jp,3)
      Bmat3(1)=NUCBFC(jp,1)
      Bmat3(2)=NUCBFC(jp,2)
      Bmat3(3)=NUCBFC(jp,3)

C  Get Nuclear Basis Function Contraction Coefficients
      Cof_ip=AGNBFCC(ip)
      Cof_jp=AGNBFCC(jp)

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

C  Get primitive Electron Basis Function Contraction Coefficients 
              Cof_ie1=AGEBFCC(ie1)
              Cof_ie2=AGEBFCC(ie2)
              Cof_je1=AGEBFCC(je1)
              Cof_je2=AGEBFCC(je2)

C---------------------OMG_12-------------------------------------------(
              call xcalc_GAM2_MD(I1,J1,K1,A1,Amat1,
     x                           I2,J2,K2,A2,Amat2,
     x                           I3,J3,K3,A3,Amat3,
     x                           L1,M1,N1,B1,Bmat1,
     x                           L2,M2,N2,B2,Bmat2,
     x                           L3,M3,N3,B3,Bmat3,
     x                           nat,ngtg1,
     x                           pmass,cat,zan,
     x                           bcoef1,gamma1,
     x                           ansE,ansS)

                        call underflow(ansE)
                        call underflow(ansS)

                        OMG2_12=OMG2_12+ansE
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2

                        OMG2s_12=OMG2s_12+ansS
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2
C---------------------OMG_12-------------------------------------------)
C---------------------OMG_21-------------------------------------------(
              call xcalc_GAM2_MD(I2,J2,K2,A2,Amat2,
     x                           I1,J1,K1,A1,Amat1,
     x                           I3,J3,K3,A3,Amat3,
     x                           L2,M2,N2,B2,Bmat2,
     x                           L1,M1,N1,B1,Bmat1,
     x                           L3,M3,N3,B3,Bmat3,
     x                           nat,ngtg1,
     x                           pmass,cat,zan,
     x                           bcoef1,gamma1,
     x                           ansE,ansS)

                        call underflow(ansE)
                        call underflow(ansS)

                        OMG2_21=OMG2_21+ansE
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2

                        OMG2s_21=OMG2s_21+ansS
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2
C---------------------OMG_21-------------------------------------------)

         end do
        end do
       end do
      end do

C  Symmetrize 
      OMG2=half*(OMG2_12+OMG2_21)
      OMG2s=half*(OMG2s_12+OMG2s_21)


      return
      end

