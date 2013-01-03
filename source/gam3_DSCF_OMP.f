!=======================================================================
      subroutine GAM3_DSCF(Nchunks,nebf,npebf,npbf,
     x                     ng3,ng3prm,nat,ngtg1,
     x                     pmass,cat,zan,bcoef1,gamma1,
     x                     KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                     ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                     DE,DP,focke,fockp,E_gam3)

!=======================================================================
      implicit none
! Input Variables
      integer Nchunks
      integer ng3,nebf,npebf,npbf,ng3prm
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
      double precision E_gam3

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
!     double precision,allocatable :: GAM_3(:)
      double precision xfocke(nebf,nebf),xfockp(npbf,npbf)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_gam3=0.0d+00
      xfocke=0.0d+00
      xfockp=0.0d+00
!     call zero_out(nebf,xfocke)
!     call zero_out(npbf,xfockp)
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

!        if(allocated(GAM_3)) deallocate(GAM_3)
!        allocate( GAM_3(ng3_seg),stat=istat )
!        write(*,*) 'allocate GAM_3: ',istat

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

         call thread_gam3(istart,iend,ng3_seg,ng3,
     x                    nebf,npebf,npbf,nat,ngtg1,
     x                    pmass,cat,zan,bcoef1,gamma1,loop_map,
     x                    KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                    ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                    DE,DP,xfocke,xfockp,E_gam3)

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
      subroutine thread_gam3(istart,iend,ng3_seg,ng3,
     x                       nebf,npebf,npbf,nat,ngtg1,
     x                       pmass,cat,zan,bcoef1,gamma1,loop_map,
     x                       KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                       DE,DP,xfocke,xfockp,E_gam3)

C======================================================================
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
      double precision DE(nebf,nebf),DP(npbf,npbf)
C Variables Returned
!     double precision GAM_3(ng3_seg)
      double precision xfocke(nebf,nebf),xfockp(npbf,npbf)
      double precision E_gam3

C Local Variables
      integer imap,ia
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

!        write(*,*)
!        write(*,*)'**************************************'
!        write(*,*)'    Computing GAM_3 Integrals    '
!        write(*,*)
!        write(*,*)'nebf      =',nebf
!        write(*,*)'npbf      =',npbf
!        write(*,*)'Total ng3 =',ng3
!        write(*,*)'ISTART    =',istart
!        write(*,*)'IEND      =',iend
!        write(*,*)' Available processors: ',omp_get_num_procs()
!        write(*,*)' Available threads     ',omp_get_max_threads()
!        write(*,*)' Threads in use        ',omp_get_num_threads()
!        write(*,*)'**************************************'
!        write(*,*)

C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf,ng3_seg)
!$ompx shared(ng3)
!$ompx shared(DE,DP)
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
!$ompx reduction(+:E_gam3)
!$ompx reduction(+:xfocke)
!$ompx reduction(+:xfockp)
!XXXX!$ompx shared(GAM_3)

      id= omp_get_thread_num()
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

!        call contract_omega3(ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,
!    x                        nebf,npebf,npbf,nat,ngtg1,
!    x                        pmass,cat,zan,bcoef1,gamma1,
!    x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
!    x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
!    x                        OMG3)

!        GAM_3(imap)=OMG3
C--------------------------FOR-DSCF------------------------------------(
         call contract_omega3(ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        OMG3_123)

         call contract_omega3(ip,jp,iec1,jec2,iec2,jec1,iec3,jec3,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        OMG3_213)

         call contract_omega3(ip,jp,iec1,jec3,iec2,jec1,iec3,jec2,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        OMG3_312)

         call contract_omega3(ip,jp,iec1,jec1,iec2,jec3,iec3,jec2,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        OMG3_132)

         call contract_omega3(ip,jp,iec1,jec2,iec2,jec3,iec3,jec1,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        OMG3_231)

         call contract_omega3(ip,jp,iec1,jec3,iec2,jec2,iec3,jec1,
     x                        nebf,npebf,npbf,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        OMG3_321)


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
!        write(*,*)'FINISHED CALCULATING CONTRACTED GAM_3 INTEGRALS:'
!        write(*,*)'ISTART=',istart
!        write(*,*)'IEND  =',iend
!        write(*,*)'WRITING THEM TO DISK...'
!        write(*,*)


      return
      end

C======================================================================
      subroutine contract_omega3(ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,
     x                           nebf,npebf,npbf,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                           OMG3)

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

C Variables Returned
      double precision OMG3

C Local Variables
      integer ie1,je1
      integer ie2,je2
      integer ie3,je3
      integer ie1_start
      integer ie2_start
      integer ie3_start
      integer je1_start
      integer je2_start
      integer je3_start
      integer ie1_end
      integer ie2_end
      integer ie3_end
      integer je1_end
      integer je2_end
      integer je3_end

      double precision Cof_ie1,Cof_je1
      double precision Cof_ie2,Cof_je2
      double precision Cof_ie3,Cof_je3
      double precision Cof_ip,Cof_jp
C--------------------------------(
C Basis set-related local variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1,Amat1(3) 
      double precision A2,Amat2(3) 
      double precision A3,Amat3(3) 
      double precision A4,Amat4(3) 
      double precision B1,Bmat1(3) 
      double precision B2,Bmat2(3) 
      double precision B3,Bmat3(3) 
      double precision B4,Bmat4(3) 
C--------------------------------)
      double precision ans
      double precision OMG3_123
      double precision OMG3_132 
      double precision OMG3_213
      double precision OMG3_231 
      double precision OMG3_312 
      double precision OMG3_321
      double precision six
      parameter(six=6.0d+00)


      ie1_start=KPESTR(iec1)
      ie2_start=KPESTR(iec2)
      ie3_start=KPESTR(iec3)

      je1_start=KPESTR(jec1)
      je2_start=KPESTR(jec2)
      je3_start=KPESTR(jec3)

      ie1_end=KPEEND(iec1)
      ie2_end=KPEEND(iec2)
      ie3_end=KPEEND(iec3)

      je1_end=KPEEND(jec1)
      je2_end=KPEEND(jec2)
      je3_end=KPEEND(jec3)

!     OMG3=0.0d+00
      OMG3_123=0.0d+00
      OMG3_132=0.0d+00
      OMG3_213=0.0d+00
      OMG3_231=0.0d+00
      OMG3_312=0.0d+00
      OMG3_321=0.0d+00

      do ie1=ie1_start,ie1_end
       do je1=je1_start,je1_end
        do ie2=ie2_start,ie2_end
         do je2=je2_start,je2_end
          do ie3=ie3_start,ie3_end
           do je3=je3_start,je3_end
               
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

              A3=ELCEX(ie3)
              I3=ELCAM(ie3,1)
              J3=ELCAM(ie3,2)
              K3=ELCAM(ie3,3)
              Amat3(1)=ELCBFC(ie3,1)
              Amat3(2)=ELCBFC(ie3,2)
              Amat3(3)=ELCBFC(ie3,3)

              A4=NUCEX(ip)
              I4=NUCAM(ip,1)
              J4=NUCAM(ip,2)
              K4=NUCAM(ip,3)
              Amat4(1)=NUCBFC(ip,1)
              Amat4(2)=NUCBFC(ip,2)
              Amat4(3)=NUCBFC(ip,3)

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

              B3=ELCEX(je3)
              L3=ELCAM(je3,1)
              M3=ELCAM(je3,2)
              N3=ELCAM(je3,3)
              Bmat3(1)=ELCBFC(je3,1)
              Bmat3(2)=ELCBFC(je3,2)
              Bmat3(3)=ELCBFC(je3,3)

              B4=NUCEX(jp)
              L4=NUCAM(jp,1)
              M4=NUCAM(jp,2)
              N4=NUCAM(jp,3)
              Bmat4(1)=NUCBFC(jp,1)
              Bmat4(2)=NUCBFC(jp,2)
              Bmat4(3)=NUCBFC(jp,3)

C  Get primitive Electron Basis Function Contraction Coefficients 
              Cof_ie1=AGEBFCC(ie1)
              Cof_ie2=AGEBFCC(ie2)
              Cof_ie3=AGEBFCC(ie3)
              Cof_je1=AGEBFCC(je1)
              Cof_je2=AGEBFCC(je2)
              Cof_je3=AGEBFCC(je3)
C  Get Nuclear Basis Function Contraction Coefficients
              Cof_ip=AGNBFCC(ip)
              Cof_jp=AGNBFCC(jp)

C---------------------OMG_123------------------------------------------(
              call xcalc_GAM3_MD(I1,J1,K1,A1,Amat1,
     x                           I2,J2,K2,A2,Amat2,
     x                           I3,J3,K3,A3,Amat3,
     x                           I4,J4,K4,A4,Amat4,
     x                           L1,M1,N1,B1,Bmat1,
     x                           L2,M2,N2,B2,Bmat2,
     x                           L3,M3,N3,B3,Bmat3,
     x                           L4,M4,N4,B4,Bmat4,
     x                           nat,ngtg1,
     x                           pmass,cat,zan,
     x                           bcoef1,gamma1,
     x                           ans)

c                       call underflow(ans)

                        OMG3_123=OMG3_123+ans
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2
     x                      *Cof_ie3*Cof_je3
C---------------------OMG_123------------------------------------------)
C---------------------OMG_132------------------------------------------(
              call xcalc_GAM3_MD(I1,J1,K1,A1,Amat1,
     x                           I3,J3,K3,A3,Amat3,
     x                           I2,J2,K2,A2,Amat2,
     x                           I4,J4,K4,A4,Amat4,
     x                           L1,M1,N1,B1,Bmat1,
     x                           L3,M3,N3,B3,Bmat3,
     x                           L2,M2,N2,B2,Bmat2,
     x                           L4,M4,N4,B4,Bmat4,
     x                           nat,ngtg1,
     x                           pmass,cat,zan,
     x                           bcoef1,gamma1,
     x                           ans)

c                       call underflow(ans)

                        OMG3_132=OMG3_132+ans
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2
     x                      *Cof_ie3*Cof_je3
C---------------------OMG_132------------------------------------------)
C---------------------OMG_213------------------------------------------(
              call xcalc_GAM3_MD(I2,J2,K2,A2,Amat2,
     x                           I1,J1,K1,A1,Amat1,
     x                           I3,J3,K3,A3,Amat3,
     x                           I4,J4,K4,A4,Amat4,
     x                           L2,M2,N2,B2,Bmat2,
     x                           L1,M1,N1,B1,Bmat1,
     x                           L3,M3,N3,B3,Bmat3,
     x                           L4,M4,N4,B4,Bmat4,
     x                           nat,ngtg1,
     x                           pmass,cat,zan,
     x                           bcoef1,gamma1,
     x                           ans)

c                       call underflow(ans)

                        OMG3_213=OMG3_213+ans
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2
     x                      *Cof_ie3*Cof_je3
C---------------------OMG_213------------------------------------------)
C---------------------OMG_231------------------------------------------(
              call xcalc_GAM3_MD(I2,J2,K2,A2,Amat2,
     x                           I3,J3,K3,A3,Amat3,
     x                           I1,J1,K1,A1,Amat1,
     x                           I4,J4,K4,A4,Amat4,
     x                           L2,M2,N2,B2,Bmat2,
     x                           L3,M3,N3,B3,Bmat3,
     x                           L1,M1,N1,B1,Bmat1,
     x                           L4,M4,N4,B4,Bmat4,
     x                           nat,ngtg1,
     x                           pmass,cat,zan,
     x                           bcoef1,gamma1,
     x                           ans)

c                       call underflow(ans)

                        OMG3_231=OMG3_231+ans
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2
     x                      *Cof_ie3*Cof_je3
C---------------------OMG_231------------------------------------------)
C---------------------OMG_312------------------------------------------(
              call xcalc_GAM3_MD(I3,J3,K3,A3,Amat3,
     x                           I1,J1,K1,A1,Amat1,
     x                           I2,J2,K2,A2,Amat2,
     x                           I4,J4,K4,A4,Amat4,
     x                           L3,M3,N3,B3,Bmat3,
     x                           L1,M1,N1,B1,Bmat1,
     x                           L2,M2,N2,B2,Bmat2,
     x                           L4,M4,N4,B4,Bmat4,
     x                           nat,ngtg1,
     x                           pmass,cat,zan,
     x                           bcoef1,gamma1,
     x                           ans)

c                       call underflow(ans)

                        OMG3_312=OMG3_312+ans
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2
     x                      *Cof_ie3*Cof_je3
C---------------------OMG_312------------------------------------------)
C---------------------OMG_321------------------------------------------(
              call xcalc_GAM3_MD(I3,J3,K3,A3,Amat3,
     x                           I2,J2,K2,A2,Amat2,
     x                           I1,J1,K1,A1,Amat1,
     x                           I4,J4,K4,A4,Amat4,
     x                           L3,M3,N3,B3,Bmat3,
     x                           L2,M2,N2,B2,Bmat2,
     x                           L1,M1,N1,B1,Bmat1,
     x                           L4,M4,N4,B4,Bmat4,
     x                           nat,ngtg1,
     x                           pmass,cat,zan,
     x                           bcoef1,gamma1,
     x                           ans)

c                       call underflow(ans)

                        OMG3_321=OMG3_321+ans
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2
     x                      *Cof_ie3*Cof_je3
C---------------------OMG_321------------------------------------------)

           end do
          end do
         end do
        end do
       end do
      end do

C  Symmetrize 
      OMG3=(OMG3_123+OMG3_132+OMG3_213
     x     +OMG3_231+OMG3_312+OMG3_321)/six


      return
      end

