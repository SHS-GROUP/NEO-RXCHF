!=======================================================================
      subroutine DECOMPOSE_GAM3(Nchunks,nebf,npebf,npbf,ng3,ng3prm,
     x                          nat,ngtg1,cat,zan,bcoef1,gamma1,
     x                          DE,DP,S_TOTAL,
     x                          KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                          ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

!=======================================================================
      implicit none

! Input Variables
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)

      integer ng3    ! Number contracted gamma3 integrals
      integer ng3prm ! Number primitive gamma3 integrals
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions

      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas

      integer Nchunks

C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
C-------Basis Set Info-------)

      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)

      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
      integer KPESTR(nebf)  ! Map contracted index to primitive start
      integer KPEEND(nebf)  ! Map contracted index to primitive end


      double precision S_TOTAL 

! Variables Returned

! Local Variables
      character*20 FLNM
      integer IFIL
      double precision E_TOT 
      double precision E_gVEE  
      double precision E_gTEg  
      double precision E_gVECg 
      double precision E_gVEPg 
      double precision E_gVEEg1
      double precision E_gVEEg2

! Dummy variable:  pmass
      pmass=1.0d+00

      write(*,*)'IN DECOMPOSE_GAM3 about to call GAM3_DCMP_SET'
      call GAM3_DCMP_SET(nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
     x                   Nchunks,pmass,cat,zan,bcoef1,gamma1,
     x                   KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                   ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

! Evaluate gVEE contribution to GAM3:
      IFIL=920
      FLNM='GM3_gVEE.ufm'
      call GAM3_ELEMENT(FLNM,IFIL,nebf,npbf,ng3,
     x                  DE,DP,E_gVEE)

! Evaluate gTEg contribution to GAM3:
      IFIL=921
      FLNM='GM3_gTEg.ufm'
      call GAM3_ELEMENT(FLNM,IFIL,nebf,npbf,ng3,
     x                  DE,DP,E_gTEg)

! Evaluate gVECg contribution to GAM3:
      IFIL=922
      FLNM='GM3_gVECg.ufm'
      call GAM3_ELEMENT(FLNM,IFIL,nebf,npbf,ng3,
     x                  DE,DP,E_gVECg)

! Evaluate gVEPg contribution to GAM3:
      IFIL=923
      FLNM='GM3_gVEPg.ufm'
      call GAM3_ELEMENT(FLNM,IFIL,nebf,npbf,ng3,
     x                  DE,DP,E_gVEPg)

! Evaluate gVEEg1 contribution to GAM3:
      IFIL=924
      FLNM='GM3_gVEEg1.ufm'
      call GAM3_ELEMENT(FLNM,IFIL,nebf,npbf,ng3,
     x                  DE,DP,E_gVEEg1)

! Evaluate gVEEg2 contribution to GAM3:
      IFIL=925
      FLNM='GM3_gVEEg2.ufm'
      call GAM3_ELEMENT(FLNM,IFIL,nebf,npbf,ng3,
     x                  DE,DP,E_gVEEg2)

      E_TOT=E_gVEE  
     x  +E_gTEg  
     x  +E_gVECg 
     x  +E_gVEPg 
     x  +E_gVEEg1
     x  +E_gVEEg2

      write(*,9100)E_gVEE,E_gTEg,E_gVECg,E_gVEPg,E_gVEEg1,E_gVEEg2,
     xE_TOT

      E_gVEE   = E_gVEE   / S_TOTAL
      E_gTEg   = E_gTEg   / S_TOTAL 
      E_gVECg  = E_gVECg  / S_TOTAL 
      E_gVEPg  = E_gVEPg  / S_TOTAL 
      E_gVEEg1 = E_gVEEg1 / S_TOTAL 
      E_gVEEg2 = E_gVEEg2 / S_TOTAL 
      E_TOT    = E_TOT    / S_TOTAL 

      write(*,9200)E_gVEE,E_gTEg,E_gVECg,E_gVEPg,E_gVEEg1,E_gVEEg2,
     xE_TOT


 9100 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|          GAM3 ENERGETIC COMPONENTS          |',/,
     x        6X,'|           -NOT DIVIDED BY OVERLAP-          |',/,
     x        6X,'+---------------------------------------------+',/,
     x       12X,'     gVEE=',1X,F20.10/
     x       12X,'     gTEg=',1X,F20.10/
     x       12X,'    gVECg=',1X,F20.10/
     x       12X,'    gVEPg=',1X,F20.10/
     x       12X,'   gVEEg1=',1X,F20.10/
     x       12X,'   gVEEg2=',1X,F20.10/
     x       12X,'   ----------------------------',/,
     x       12X,'    TOTAL=',1X,F20.10/)

 9200 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|          GAM3 ENERGETIC COMPONENTS          |',/,
     x        6X,'|             -DIVIDED BY OVERLAP-            |',/,
     x        6X,'+---------------------------------------------+',/,
     x       12X,'     gVEE=',1X,F20.10/
     x       12X,'     gTEg=',1X,F20.10/
     x       12X,'    gVECg=',1X,F20.10/
     x       12X,'    gVEPg=',1X,F20.10/
     x       12X,'   gVEEg1=',1X,F20.10/
     x       12X,'   gVEEg2=',1X,F20.10/
     x       12X,'   ----------------------------',/,
     x       12X,'    TOTAL=',1X,F20.10/)


      return
      end
!=======================================================================
      subroutine GAM3_ELEMENT(FLNM,IFIL,nebf,npbf,ng3,
     x                        DE,DP,E_gam3)

!=======================================================================
      implicit none

! Input Variables
      character*20 FLNM
      integer IFIL

      double precision zero,half,two,three,four,eight
      parameter(zero=0.0d+00,half=0.5d+00,two=2.0d+00)
      parameter(three=3.0d+00,four=4.0d+00,eight=8.0d+00)

      integer nebf
      integer npbf
      integer ng3

      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision E_gam3

! Local Variables
      integer ie1
      integer je1
      integer ie2
      integer je2
      integer ie3
      integer je3
      integer ip
      integer jp

      integer ia_123
      integer ia_213
      integer ia_312
      integer ia_132
      integer ia_231
      integer ia_321

      double precision val_gam3
      double precision xfocke(nebf,nebf)
      double precision xfockp(npbf,npbf)

      double precision x123
      double precision x213
      double precision x312
      double precision x132
      double precision x231
      double precision x321


      open(IFIL,file=FLNM,form='unformatted',
     x status='unknown',access='direct',RECL=8)

      E_gam3=zero

      do ip=1,npbf
      do jp=1,npbf

            do ie1=1,nebf
            do je1=1,nebf

                  do ie2=1,nebf
                  do je2=1,nebf

                     do ie3=1,nebf
                     do je3=1,nebf
CCCC
C    AS Packed--> GAM_3(ip,jp,je1,ie1,je2,ie2,je3,ie3)
C                   = GAM_3PK(je3,ie3,je2,ie2,je1,ie1,jp,ip)

      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je1,ie2,je2,ie3,je3,ia_123)
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je2,ie2,je1,ie3,je3,ia_213)
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je3,ie2,je1,ie3,je2,ia_312)
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je1,ie2,je3,ie3,je2,ia_132)
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je2,ie2,je3,ie3,je1,ia_231)
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je3,ie2,je2,ie3,je1,ia_321)

CCWS-IO  start
                        read(IFIL,REC=ia_123) x123
                        read(IFIL,REC=ia_213) x213
                        read(IFIL,REC=ia_312) x312
                        read(IFIL,REC=ia_132) x132
                        read(IFIL,REC=ia_231) x231
                        read(IFIL,REC=ia_321) x321

                  val_gam3=eight*x123 -
     x                      four*x213 +
     x                       two*x312 -
     x                      four*x132 +
     x                       two*x231 -
     x                      four*x321
CCWS-IO  end
                  val_gam3=val_gam3/eight

            E_gam3=E_gam3+
     x            DE(ie1,je1)*DE(ie2,je2)*DE(ie3,je3)*DP(ip,jp)*val_gam3

                     end do
                     end do

                  end do
                  end do

            end do
            end do

      end do
      end do

CCWS-IO
      close(IFIL)


      return
      end
!=======================================================================
      subroutine GAM3_DCMP_SET(nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
     x                         Nchunks,pmass,cat,zan,bcoef1,gamma1,
     x                         KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                         ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!=======================================================================
      implicit none
C Input Variables
      integer ng3    ! Number contracted gamma3 integrals
      integer ng3prm ! Number primitive gamma3 integrals
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions

      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas

      integer Nchunks
C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
C-------Basis Set Info-------)

      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)

      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
      integer KPESTR(nebf)  ! Map contracted index to primitive start
      integer KPEEND(nebf)  ! Map contracted index to primitive end

! Local Variables
      logical LgVEE,LgTEg,LgVECg
      logical LgVEPg,LgVEEg1,LgVEEg2
      character*20 FLNM
      integer IFIL


!----------------------------------------------------------------------(
! Evaluate gVEE contribution to GAM3:
      IFIL=920
      FLNM='GM3_gVEE.ufm'
      LgVEE=.TRUE.
      LgTEg=.FALSE.
      LgVECg=.FALSE.
      LgVEPg=.FALSE.
      LgVEEg1=.FALSE.
      LgVEEg2=.FALSE.

      write(*,*)'IN GAM3_DCMP_SET about to call GAM3_OMP_DC'
      call GAM3_CONV_DC(LgVEE,LgTEg,LgVECg,
     x                  LgVEPg,LgVEEg1,LgVEEg2,
     x                  FLNM,IFIL,Nchunks,
     x                  nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
     x                  pmass,cat,zan,bcoef1,gamma1,
     x                  KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                  ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

!     call GAM3_OMP_DC(LgVEE,LgTEg,LgVECg,
!    x                 LgVEPg,LgVEEg1,LgVEEg2,
!    x                 FLNM,IFIL,
!    x                 nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
!    x                 pmass,cat,zan,bcoef1,gamma1,
!    x                 AMPEB2C,AGEBFCC,AGNBFCC,
!    x                 ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!----------------------------------------------------------------------)
!----------------------------------------------------------------------(
! Evaluate gTEg contribution to GAM3:
      IFIL=921
      FLNM='GM3_gTEg.ufm'
      LgVEE=.FALSE.
      LgTEg=.TRUE.
      LgVECg=.FALSE.
      LgVEPg=.FALSE.
      LgVEEg1=.FALSE.
      LgVEEg2=.FALSE.

      call GAM3_CONV_DC(LgVEE,LgTEg,LgVECg,
     x                  LgVEPg,LgVEEg1,LgVEEg2,
     x                  FLNM,IFIL,Nchunks,
     x                  nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
     x                  pmass,cat,zan,bcoef1,gamma1,
     x                  KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                  ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!     call GAM3_OMP_DC(LgVEE,LgTEg,LgVECg,
!    x                 LgVEPg,LgVEEg1,LgVEEg2,
!    x                 FLNM,IFIL,
!    x                 nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
!    x                 pmass,cat,zan,bcoef1,gamma1,
!    x                 AMPEB2C,AGEBFCC,AGNBFCC,
!    x                 ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!----------------------------------------------------------------------)
!----------------------------------------------------------------------(
! Evaluate gVECg contribution to GAM3:
      IFIL=922
      FLNM='GM3_gVECg.ufm'
      LgVEE=.FALSE.
      LgTEg=.FALSE.
      LgVECg=.TRUE.
      LgVEPg=.FALSE.
      LgVEEg1=.FALSE.
      LgVEEg2=.FALSE.

      call GAM3_CONV_DC(LgVEE,LgTEg,LgVECg,
     x                  LgVEPg,LgVEEg1,LgVEEg2,
     x                  FLNM,IFIL,Nchunks,
     x                  nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
     x                  pmass,cat,zan,bcoef1,gamma1,
     x                  KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                  ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!     call GAM3_OMP_DC(LgVEE,LgTEg,LgVECg,
!    x                 LgVEPg,LgVEEg1,LgVEEg2,
!    x                 FLNM,IFIL,
!    x                 nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
!    x                 pmass,cat,zan,bcoef1,gamma1,
!    x                 AMPEB2C,AGEBFCC,AGNBFCC,
!    x                 ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!----------------------------------------------------------------------)
!----------------------------------------------------------------------(
! Evaluate gVEPg contribution to GAM3:
      IFIL=923
      FLNM='GM3_gVEPg.ufm'
      LgVEE=.FALSE.
      LgTEg=.FALSE.
      LgVECg=.FALSE.
      LgVEPg=.TRUE.
      LgVEEg1=.FALSE.
      LgVEEg2=.FALSE.

      call GAM3_CONV_DC(LgVEE,LgTEg,LgVECg,
     x                  LgVEPg,LgVEEg1,LgVEEg2,
     x                  FLNM,IFIL,Nchunks,
     x                  nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
     x                  pmass,cat,zan,bcoef1,gamma1,
     x                  KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                  ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!     call GAM3_OMP_DC(LgVEE,LgTEg,LgVECg,
!    x                 LgVEPg,LgVEEg1,LgVEEg2,
!    x                 FLNM,IFIL,
!    x                 nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
!    x                 pmass,cat,zan,bcoef1,gamma1,
!    x                 AMPEB2C,AGEBFCC,AGNBFCC,
!    x                 ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!----------------------------------------------------------------------)
!----------------------------------------------------------------------(
! Evaluate gVEEg1 contribution to GAM3:
      IFIL=924
      FLNM='GM3_gVEEg1.ufm'
      LgVEE=.FALSE.
      LgTEg=.FALSE.
      LgVECg=.FALSE.
      LgVEPg=.FALSE.
      LgVEEg1=.TRUE.
      LgVEEg2=.FALSE.

      call GAM3_CONV_DC(LgVEE,LgTEg,LgVECg,
     x                  LgVEPg,LgVEEg1,LgVEEg2,
     x                  FLNM,IFIL,Nchunks,
     x                  nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
     x                  pmass,cat,zan,bcoef1,gamma1,
     x                  KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                  ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!     call GAM3_OMP_DC(LgVEE,LgTEg,LgVECg,
!    x                 LgVEPg,LgVEEg1,LgVEEg2,
!    x                 FLNM,IFIL,
!    x                 nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
!    x                 pmass,cat,zan,bcoef1,gamma1,
!    x                 AMPEB2C,AGEBFCC,AGNBFCC,
!    x                 ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!----------------------------------------------------------------------)
!----------------------------------------------------------------------(
! Evaluate gVEEg2 contribution to GAM3:
      IFIL=925
      FLNM='GM3_gVEEg2.ufm'
      LgVEE=.FALSE.
      LgTEg=.FALSE.
      LgVECg=.FALSE.
      LgVEPg=.FALSE.
      LgVEEg1=.FALSE.
      LgVEEg2=.TRUE.

      call GAM3_CONV_DC(LgVEE,LgTEg,LgVECg,
     x                  LgVEPg,LgVEEg1,LgVEEg2,
     x                  FLNM,IFIL,Nchunks,
     x                  nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
     x                  pmass,cat,zan,bcoef1,gamma1,
     x                  KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                  ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!     call GAM3_OMP_DC(LgVEE,LgTEg,LgVECg,
!    x                 LgVEPg,LgVEEg1,LgVEEg2,
!    x                 FLNM,IFIL,
!    x                 nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
!    x                 pmass,cat,zan,bcoef1,gamma1,
!    x                 AMPEB2C,AGEBFCC,AGNBFCC,
!    x                 ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
!----------------------------------------------------------------------)


      return
      end
C======================================================================
      subroutine GAM3_OMP_DC(LgVEE,LgTEg,LgVECg,
     x                       LgVEPg,LgVEEg1,LgVEEg2,
     x                       FLNM,IFIL,
     x                       nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
     x                       pmass,cat,zan,bcoef1,gamma1,
     x                       AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

C======================================================================
      implicit none
      include 'omp_lib.h'

C Input Variables
      logical LgVEE,LgTEg,LgVECg
      logical LgVEPg,LgVEEg1,LgVEEg2
      character*20 FLNM
      integer IFIL
      integer ng3    ! Number contracted gamma3 integrals
      integer ng3prm ! Number primitive gamma3 integrals
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions

      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas

C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
C-------Basis Set Info-------)

      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)

      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
C Local variables
      integer ia
      integer ia_123
      integer ia_132
      integer ia_213
      integer ia_231
      integer ia_312
      integer ia_321
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3
      integer iec1,jec1,iec2,jec2,iec3,jec3
      double precision Cof_ie1,Cof_je1
      double precision Cof_ie2,Cof_je2
      double precision Cof_ie3,Cof_je3
      double precision Cof_ip,Cof_jp

      double precision zero,half,six
      parameter(zero=0.0d+00,half=0.5d+00,six=6.0d+00)
      double precision ans

      double precision x123
      double precision x132
      double precision x213
      double precision x231
      double precision x312
      double precision x321
      double precision xxxx

C-------------------------------------
      double precision, allocatable :: GAM_3(:)
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
      integer istat

      double precision A1,Amat1(3) 
      double precision A2,Amat2(3) 
      double precision A3,Amat3(3) 
      double precision A4,Amat4(3) 
      double precision B1,Bmat1(3) 
      double precision B2,Bmat2(3) 
      double precision B3,Bmat3(3) 
      double precision B4,Bmat4(3) 
C--------------------------------)
C---OPENMP-RELATED-VARIABLES-----(
      integer id
      integer loopi,iLP
      double precision wtime
      integer loop_map(ng3prm,8)
C---OPENMP-RELATED-VARIABLES-----)

      write(*,*)'IN GAM3_OMP_DC about to allocate gam3'
      call flush(6)
      if(allocated(GAM_3)) deallocate(GAM_3)
      allocate( GAM_3(ng3),stat=istat )


C---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime()
C---OPENMP-TIMING------------------------------------------------------)

C-----------INITIALIZE-DATA-STRUCTURES-ON-MASTER-----------------------(
         write(*,*)
         write(*,*)'**************************************'
         write(*,*)'         GAM_3 DECOMPOSITION    '
         write(*,*)
         write(*,*)'METHOD = MD'
         write(*,*)'FILE NAME:  ',FLNM
         write(*,*)'npebf=',npebf
         write(*,*)'nebf=',nebf
         write(*,*)'npbf=',npbf
         write(*,*)'nprim_g3 =',ng3prm
         write(*,*)'ng3 =',ng3
         write(*,*)' Available processors: ',omp_get_num_procs()
         write(*,*)' Available threads     ',omp_get_max_threads()
         write(*,*)' Threads in use        ',omp_get_num_threads()
         write(*,*)'**************************************'
         write(*,*)

C Initialize 1-D arrays to hold the contracted integrals
         do ia=1,ng3
            GAM_3(ia)=zero
         end do

C Compress nested loops to 1-dimension
         Loopi=0
         do ip=1,npbf
         do jp=1,npbf
            do ie1=1,npebf
            do je1=1,npebf
               do ie2=1,npebf
               do je2=1,npebf
                  do ie3=1,npebf
                  do je3=1,npebf

                           Loopi=Loopi+1
                           loop_map(Loopi,1)=je3
                           loop_map(Loopi,2)=ie3
                           loop_map(Loopi,3)=je2
                           loop_map(Loopi,4)=ie2
                           loop_map(Loopi,5)=je1
                           loop_map(Loopi,6)=ie1
                           loop_map(Loopi,7)=jp
                           loop_map(Loopi,8)=ip

                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(LgVEE,LgTEg,LgVECg)
!$ompx shared(LgVEPg,LgVEEg1,LgVEEg2)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf)
!$ompx shared(ng3)
!$ompx shared(ng3prm)
!$ompx shared(GAM_3)
!$ompx private(iLp) 
!$ompx private(ia)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(ie3,je3) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(A1,I1,J1,K1,Amat1)
!$ompx private(A2,I2,J2,K2,Amat2)
!$ompx private(A3,I3,J3,K3,Amat3)
!$ompx private(A4,I4,J4,K4,Amat4)
!$ompx private(B1,L1,M1,N1,Bmat1)
!$ompx private(B2,L2,M2,N2,Bmat2)
!$ompx private(B3,L3,M3,N3,Bmat3)
!$ompx private(B4,L4,M4,N4,Bmat4)
!$ompx private(Cof_ie1,Cof_je1)
!$ompx private(Cof_ie2,Cof_je2)
!$ompx private(Cof_ie3,Cof_je3)
!$ompx private(Cof_ip,Cof_jp)
!$ompx private(ans)
!$ompx private(id)

      id= omp_get_thread_num()
      write(*,*)' Hello from process ',id
      if(id.eq.0) then
         write(*,*)'Threads in use', omp_get_num_threads()
      end if

!$omp do
      do iLP=1,ng3prm

C  Map loop indices
         je3=loop_map(iLP,1)
         ie3=loop_map(iLP,2)
         je2=loop_map(iLP,3)
         ie2=loop_map(iLP,4)
         je1=loop_map(iLP,5)
         ie1=loop_map(iLP,6)
         jp =loop_map(iLP,7)
         ip =loop_map(iLP,8)

C Get Basis set info:
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

         call xcalc_GAM3_DC(LgVEE,LgTEg,LgVECg,
     x                      LgVEPg,LgVEEg1,LgVEEg2,
     x                      I1,J1,K1,A1,Amat1,
     x                      I2,J2,K2,A2,Amat2,
     x                      I3,J3,K3,A3,Amat3,
     x                      I4,J4,K4,A4,Amat4,
     x                      L1,M1,N1,B1,Bmat1,
     x                      L2,M2,N2,B2,Bmat2,
     x                      L3,M3,N3,B3,Bmat3,
     x                      L4,M4,N4,B4,Bmat4,
     x                      nat,ngtg1,
     x                      pmass,cat,zan,
     x                      bcoef1,gamma1,
     x                      ans) 
!    x                      ans_gVEE, 
!    x                      ans_gTEg, 
!    x                      ans_gVECg,
!    x                      ans_gVEPg,
!    x                      ans_gVEEg1,
!    x                      ans_gVEEg2)

c                       call underflow(ans)

C  Map from primitive BF indices to contracted indices
                        iec1=AMPEB2C(ie1)
                        iec2=AMPEB2C(ie2)
                        iec3=AMPEB2C(ie3)
                        jec1=AMPEB2C(je1)
                        jec2=AMPEB2C(je2)
                        jec3=AMPEB2C(je3)
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
C  Map the 8-index contracted integral to 1-D:
                        call index_GAM_3PK(nebf,npbf,
     x                                     ip,jp,
     x                                     iec1,jec1,
     x                                     iec2,jec2,
     x                                     iec3,jec3,ia)

                        GAM_3(ia)=GAM_3(ia)+ans
     x                            *Cof_ip*Cof_jp
     x                            *Cof_ie1*Cof_je1
     x                            *Cof_ie2*Cof_je2
     x                            *Cof_ie3*Cof_je3
  
         end do
!$omp end do
!$omp end parallel      
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

C---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime() - wtime
      write(*,*)'TIME TO CALCULATE GAM_3 INTEGRALS: ',wtime
C---OPENMP-TIMING------------------------------------------------------)

CCWS-IO
         write(*,*)
         write(*,*)'FINISHED CALCULATING CONTRACTED GAM_3 INTEGRALS'
         write(*,*)'BEGINNING SYMMETRIZATION OF GAM_3 INTEGRALS'
         write(*,*)

         open(IFIL,file=FLNM,form='unformatted',
     x    status='unknown',access='direct',RECL=8)

      do ip=1,npbf
      do jp=1,npbf
         do ie1=1,nebf
         do je1=1,nebf
            do ie2=1,nebf
            do je2=1,nebf
               do ie3=1,nebf
               do je3=1,nebf

C  GAM_2 Symmetrization:
C  Determine packing indices for XGAM_3 integral matrices
C
C              As Packed-->       XGAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je1,ie2,je2,ie3,je3,ia)
                             ia_123=ia
                               x123=GAM_3(ia_123)

C              As Packed-->       XGAM_3(je2,ie2,je3,ie3,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je1,ie3,je3,ie2,je2,ia)
                             ia_132=ia
                               x132=GAM_3(ia_132)

C              As Packed-->       XGAM_3(je3,ie3,je1,ie1,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie2,je2,ie1,je1,ie3,je3,ia)
                             ia_213=ia
                               x213=GAM_3(ia_213)

C              As Packed-->       XGAM_3(je1,ie1,je3,ie3,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie2,je2,ie3,je3,ie1,je1,ia)
                             ia_231=ia
                               x231=GAM_3(ia_231)

C              As Packed-->       XGAM_3(je2,ie2,je1,ie1,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie3,je3,ie1,je1,ie2,je2,ia)
                             ia_312=ia
                               x312=GAM_3(ia_312)

C              As Packed-->       XGAM_3(je1,ie1,je2,ie2,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie3,je3,ie2,je2,ie1,je1,ia)
                             ia_321=ia
                               x321=GAM_3(ia_321)

                       xxxx=(x123+x132+x213+x231+x312+x321)/six

c                      write(*,*)'xxxx=',xxxx
cc                     GAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3)=xxxx 
c                      GAM_3PK(ia_123)=xxxx 
c                      call put_GAM3(ia_123,ng3,xxxx)
                       write(IFIL,REC=ia_123) xxxx

cc                     GAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2)=xxxx
c                      GAM_3PK(ia_132)=xxxx 
c                      call put_GAM3(ia_132,ng3,xxxx)
                       write(IFIL,REC=ia_132) xxxx

cc                     GAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3)=xxxx
c                      GAM_3PK(ia_213)=xxxx 
c                      call put_GAM3(ia_213,ng3,xxxx)
                       write(IFIL,REC=ia_213) xxxx

cc                     GAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1)=xxxx
c                      GAM_3PK(ia_231)=xxxx 
c                      call put_GAM3(ia_231,ng3,xxxx)
                       write(IFIL,REC=ia_231) xxxx

cc                     GAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2)=xxxx
c                      GAM_3PK(ia_312)=xxxx 
c                      call put_GAM3(ia_312,ng3,xxxx)
                       write(IFIL,REC=ia_312) xxxx

cc                     GAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1)=xxxx
c                      GAM_3PK(ia_321)=xxxx 
c                      call put_GAM3(ia_321,ng3,xxxx)
                       write(IFIL,REC=ia_321) xxxx


               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      close(IFIL)
      write(*,*)'ALL DONE WITH SYMMETRIZATION'

      if(allocated(GAM_3)) deallocate(GAM_3)


      return
      end
C======================================================================
      subroutine xcalc_GAM3_DC(LgVEE,LgTEg,LgVECg,
     x                         LgVEPg,LgVEEg1,LgVEEg2,
     x                         I1,J1,K1,A1,Amat1,
     x                         I2,J2,K2,A2,Amat2,
     x                         I3,J3,K3,A3,Amat3,
     x                         I4,J4,K4,A4,Amat4,
     x                         L1,M1,N1,B1,Bmat1,
     x                         L2,M2,N2,B2,Bmat2,
     x                         L3,M3,N3,B3,Bmat3,
     x                         L4,M4,N4,B4,Bmat4,
     x                         nat,ngtg1,
     x                         pmass,cat,zan,
     x                         bcoef1,gamma1,
     x                         ans) 
!    x                         ans_gVEE, 
!    x                         ans_gTEg, 
!    x                         ans_gVECg,
!    x                         ans_gVEPg,
!    x                         ans_gVEEg1,
!    x                         ans_gVEEg2)

C======================================================================
      implicit none

C Input Variables
      logical LgVEE,LgTEg,LgVECg
      logical LgVEPg,LgVEEg1,LgVEEg2
      integer nat
      integer ngtg1
      double precision pmass
      double precision zan(nat)
      double precision cat(3,nat)
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)

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

C Variables Returned
      double precision ans  
      double precision ans_gVEE  
      double precision ans_gTEg  
      double precision ans_gVECg 
      double precision ans_gVEPg 
      double precision ans_gVEEg1
      double precision ans_gVEEg2

C Local Variables
      integer iii    ! Index for looping over natoms
      integer ik,il  ! Indices for geminal loops
      integer iat

      double precision gamA14
      double precision gamA24
      double precision gamA34
      double precision gamB14
      double precision gamB24
      double precision gamB34
      double precision gamA
      double precision gamB

      double precision cmat(3)
      double precision znuc

      double precision xx,yy,zz
      double precision zero,half,one,two,four
      parameter(zero=0.0d+00,one=1.0d+00,two=2.0d+00,four=4.0d+00)
      parameter(half=0.5d+00)

      double precision gVEE
      double precision gTEg
      double precision gVECg
      double precision gVEPg
      double precision gVEEg1
      double precision gVEEg2

      double precision xgVEE
      double precision xgTEg
      double precision xgVECg
      double precision xgVEPg
      double precision xgVEEg1
      double precision xgVEEg2

      double precision xmass,coulomb_sign
      double precision xke,Vc
      double precision val_vec 
    

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     BASIS FUNCTIONS: ASSIGN CENTERS, EXPONENTS, ANGULAR MOM.
C

C     *****ie1 ::  electron 1 bra *****
C     *****ie2 ::  electron 2 bra *****      
C     *****ie3 ::  electron 3 bra *****      
C     *****ip  ::  proton bra     *****

C     *****je1 ::  electron 1 ket *****
C     *****je2 ::  electron 2 ket *****
C     *****je3 ::  electron 3 ket *****
C     *****jp  ::  proton ket     *****

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      ans_gVEE  =zero 
      ans_gTEg  =zero
      ans_gVECg =zero
      ans_gVEPg =zero
      ans_gVEEg1=zero
      ans_gVEEg2=zero

      gVEE=zero
      gTEg=zero
      gVECg=zero
      gVEPg=zero
      gVEEg1=zero
      gVEEg2=zero

      xgVEE=zero
      xgTEg=zero
      xgVECg=zero
      xgVEPg=zero
      xgVEEg1=zero
      xgVEEg2=zero

      if(LgVEE) then

         DO IK=1,NGTG1

            gamA=gamma1(ik)
            gamB=zero

C  xgVEE
C  --- g(1,p)VEE(2,3)---
            call G2_MD_xggs(I1,J1,K1,A1,Amat1,
     1                      I4,J4,K4,A4,Amat4,
     2                      L1,M1,N1,B1,Bmat1,
     3                      L4,M4,N4,B4,Bmat4,
     4                      gamA,gamB,xx)

c           call pgiovlap(I1,J1,K1,A1,Amat1,
c    x                    I4,J4,K4,A4,Amat4,
c    x                    L1,M1,N1,B1,Bmat1,
c    x                    L4,M4,N4,B4,Bmat4,
c    x                    gamA,gamB,xx)


            call gfvee(I2,J2,K2,A2,Amat2,
     x                 I3,J3,K3,A3,Amat3,
     x                 L2,M2,N2,B2,Bmat2,
     x                 L3,M3,N3,B3,Bmat3,
     x                 yy)

!           call underflow(xx)
!           call underflow(yy)

            xgVEE=xx*yy
            gVEE=gVEE+(bcoef1(ik)*xgVEE)

C End 1 gamma loop
         end do
C Begin 2 gamma loop
      end if ! end-if for LgVEE

      DO IK=1,NGTG1
         DO IL=1,NGTG1

C>>>>>>>>>>>>>>>>>>>>  xgHEg <<<<<<<<<<<<<<<<<<<<

            xmass=one
            coulomb_sign=-one

            if(LgTEg) then

               call gfke(I3,J3,K3,A3,Amat3,
     x                   L3,M3,N3,B3,Bmat3,
     x                   xmass,xke)

            end if ! end-if for LgTEg

            Vc = ZERO

            if(LgVECg) then

               do iat=1,nat
                    Cmat(1)=cat(1,iat)
                    Cmat(2)=cat(2,iat)
                    Cmat(3)=cat(3,iat)
                    call gfvec(I3,J3,K3,A3,Amat3,
     x                         L3,M3,N3,B3,Bmat3,
     x                         Cmat,val_vec)
                    Vc = Vc - (zan(iat)*val_vec)
                end do

             end if ! end-if for LgVECg

            if(LgTEg.or.LgVECg) then

               gamA = gamma1(ik)
               gamB = gamma1(il)

c              call G3ovlap(I1,J1,K1,A1,Amat1,
c    x                      I2,J2,K2,A2,Amat2,
c    x                      I4,J4,K4,A4,Amat4,
c    x                      L1,M1,N1,B1,Bmat1,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      L4,M4,N4,B4,Bmat4,
c    x                      ZERO,gamA,ZERO,
c    x                      ZERO,ZERO,gamB,yy)
cc   4                   gamA12,gamA13,gamA23,
cc   4                   gamB12,gamB13,gamB23,sval)
               call G3_MD_xggs(I1,J1,K1,A1,Amat1,
     x                         I2,J2,K2,A2,Amat2,
     x                         I4,J4,K4,A4,Amat4,
     x                         L1,M1,N1,B1,Bmat1,
     x                         L2,M2,N2,B2,Bmat2,
     x                         L4,M4,N4,B4,Bmat4,
     x                         ZERO,gamA,ZERO,
     x                         ZERO,ZERO,gamB,yy)

               xgTEg=yy*xke
               xgVECg=yy*Vc

            end if ! end-if for LgTEg or LgVECg

            if(LgVEPg) then

               coulomb_sign  = -ONE
               gamA14 = gamma1(ik)
               gamA24 = ZERO
               gamA34 = ZERO
               gamB14 = ZERO
               gamB24 = gamma1(il)
               gamB34 = ZERO

               call G4_MD_xgVepg(I1,J1,K1,A1,Amat1,
     *                           I2,J2,K2,A2,Amat2,
     *                           I3,J3,K3,A3,Amat3,
     *                           I4,J4,K4,A4,Amat4,
     *                           L1,M1,N1,B1,Bmat1,
     *                           L2,M2,N2,B2,Bmat2,
     *                           L3,M3,N3,B3,Bmat3,
     *                           L4,M4,N4,B4,Bmat4,
     *                           gamA14,gamB24,
     *                           xgVEPg)
c              call G4Vep_AUX_g14g24V34(I1,J1,K1,A1,Amat1,
c    *                                  I2,J2,K2,A2,Amat2,
c    *                                  I3,J3,K3,A3,Amat3,
c    *                                  I4,J4,K4,A4,Amat4,
c    *                                  L1,M1,N1,B1,Bmat1,
c    *                                  L2,M2,N2,B2,Bmat2,
c    *                                  L3,M3,N3,B3,Bmat3,
c    *                                  L4,M4,N4,B4,Bmat4,
c    *                                  gamA14,zero,
c    *                                  zero,gamB24,
c    *                                  xgVEPg)

!               call underflow(xx)
!               call underflow(yy)
!               call underflow(zz)
                xgVEPg=-xgVEPg
  
             end if ! end-if for LgVEPg


            if(LgVEEg1) then
C>>>>>>>>>>>>>>>>>>>>  xgVEEg1 <<<<<<<<<<<<<<<<<<<<
C  --- g(1,p)VEE(2,3)g(1,p)---
               gamA=gamma1(ik)
               gamB=gamma1(il)
               call G2_MD_xggs(I1,J1,K1,A1,Amat1,
     1                         I4,J4,K4,A4,Amat4,
     2                         L1,M1,N1,B1,Bmat1,
     3                         L4,M4,N4,B4,Bmat4,
     4                         gamA,gamB,xx)

c              call pgiovlap(I1,J1,K1,A1,Amat1,
c    x                       I4,J4,K4,A4,Amat4,
c    x                       L1,M1,N1,B1,Bmat1,
c    x                       L4,M4,N4,B4,Bmat4,
c    x                       gamA,gamB,xx)


               call gfvee(I2,J2,K2,A2,Amat2,
     x                    I3,J3,K3,A3,Amat3,
     x                    L2,M2,N2,B2,Bmat2,
     x                    L3,M3,N3,B3,Bmat3,
     x                    yy)


!              call underflow(xx)
!              call underflow(yy)

               xgVEEg1=xx*yy

            end if ! end-if for LgVEEg1

            if(LgVEEg2) then
C>>>>>>>>>>>>>>>>>>>>  xgVEEg2 <<<<<<<<<<<<<<<<<<<<
               gamA14=gamma1(ik)
               gamA24=zero
               gamA34=zero
               gamB14=zero
               gamB24=zero
               gamB34=gamma1(il)

               call G4_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                           I2,J2,K2,A2,Amat2,
     *                           I3,J3,K3,A3,Amat3,
     *                           I4,J4,K4,A4,Amat4,
     *                           L1,M1,N1,B1,Bmat1,
     *                           L2,M2,N2,B2,Bmat2,
     *                           L3,M3,N3,B3,Bmat3,
     *                           L4,M4,N4,B4,Bmat4,
     *                           gamA14,gamB34,
     *                           xgVEEg2)

c              call G4Vee_AUX_g14g34V12(I1,J1,K1,A1,Amat1,
c    *                                  I2,J2,K2,A2,Amat2,
c    *                                  I3,J3,K3,A3,Amat3,
c    *                                  I4,J4,K4,A4,Amat4,
c    *                                  L1,M1,N1,B1,Bmat1,
c    *                                  L2,M2,N2,B2,Bmat2,
c    *                                  L3,M3,N3,B3,Bmat3,
c    *                                  L4,M4,N4,B4,Bmat4,
c    *                                  gamA14,gamB34,
c    *                                  xgVEEg2)

            end if ! end-if for LgVEEg2

C  Sum terms against bcoeff
            gTEg  =gTEg  +(bcoef1(ik)*bcoef1(il)*xgTEg)
            gVECg =gVECg +(bcoef1(ik)*bcoef1(il)*xgVECg)
            gVEPg =gVEPg +(bcoef1(ik)*bcoef1(il)*xgVEPg)
            gVEEg1=gVEEg1+(bcoef1(ik)*bcoef1(il)*xgVEEg1)
            gVEEg2=gVEEg2+(bcoef1(ik)*bcoef1(il)*xgVEEg2)

C End 2 gamma loop
         end do
      end do


C Total integral build
!     ans=two*gVEE*half + gHEg + gVEEg1*half + four*gVEEg2*half

      ans_gVEE  =2.0d+00*gVEE
      ans_gTEg  =gTEg
      ans_gVECg =gVECg
      ans_gVEPg =gVEPg
      ans_gVEEg1=gVEEg1
      ans_gVEEg2=4.0d+00*gVEEg2

      ans= ans_gVEE  
     x   + ans_gTEg  
     x   + ans_gVECg 
     x   + ans_gVEPg 
     x   + ans_gVEEg1
     x   + ans_gVEEg2



      return
      end

!=======================================================================
!     subroutine GAM3_CONV(Nchunks,nebf,npebf,npbf,
!    x                     ng3,ng3prm,nat,ngtg1,
!    x                     pmass,cat,zan,bcoef1,gamma1,
!    x                     KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
!    x                     ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
      subroutine GAM3_CONV_DC(LgVEE,LgTEg,LgVECg,
     x                 LgVEPg,LgVEEg1,LgVEEg2,
     x                 FLNM,IFIL,Nchunks,
     x                 nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
     x                 pmass,cat,zan,bcoef1,gamma1,
     x                 KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                 ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

!=======================================================================
      implicit none
! Input Variables
      logical LgVEE,LgTEg,LgVECg,LgVEPg,LgVEEg1,LgVEEg2
      character*20 FLNM
      integer IFIL
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

! Variables Returned

! Local Variables
      integer istat,ichunk,istart,iend,ng3_seg
      integer Loopi,imas
      integer ip,jp,iec1,jec1,iec2,jec2,iec3,jec3
      integer,allocatable :: loop_map(:,:)
!     integer,allocatable :: gam_map(:)
      double precision,allocatable :: GAM_3(:)

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



      write(*,*)
      write(*,*)'**************************************'
      write(*,*)'    Computing GAM_3 Integrals    '
      write(*,*)
      write(*,*)'FILE NAME: ',FLNM
      write(*,*)'nebf      =',nebf
      write(*,*)'npbf      =',npbf
      write(*,*)'Total ng3 =',ng3
      write(*,*)'NChunks   =',nchunks
!     write(*,*)' Available processors: ',omp_get_num_procs()
!     write(*,*)' Available threads     ',omp_get_max_threads()
      write(*,*)'**************************************'
      write(*,*)


!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      do ichunk=1,Nchunks

         call loop_size(1,ng3,Nchunks,ichunk-1,istart,iend)
!        write(*,*)'after call loop size'
!        write(*,*)'NG3=',ng3
!        write(*,*)'Nchunks=',Nchunks
         write(*,*)'Chunk Number: ',ichunk
!        write(*,*)'istart=',istart
!        write(*,*)'iend=',iend

! Segment of ng3:
         ng3_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng3_seg,8),stat=istat )
!        write(*,*) 'allocate loop_map: ',istat

         if(allocated(GAM_3)) deallocate(GAM_3)
         allocate( GAM_3(ng3_seg),stat=istat )
!        write(*,*) 'allocate GAM_3: ',istat

!        if(allocated(gam_map)) deallocate(gam_map)
!        allocate( gam_map(ng3_seg),stat=istat )
!        write(*,*) 'allocate gam_map: ',istat

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

         call thread_gam3_DC(LgVEE,LgTEg,LgVECg,
     x                       LgVEPg,LgVEEg1,LgVEEg2,
     x                       istart,iend,ng3_seg,ng3,
     x                       nebf,npebf,npbf,nat,ngtg1,
     x                       pmass,cat,zan,bcoef1,gamma1,
     x                       loop_map,GAM_3,
     x                       KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

!        call thread_gam3_DC(istart,iend,ng3_seg,ng3,
!    x                         nebf,npebf,npbf,nat,ngtg1,
!    x                         pmass,cat,zan,bcoef1,gamma1,
!    x                         loop_map,GAM_3,
!    x                         KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
!    x                         ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
      if(allocated(GAM_3)) deallocate(GAM_3)
!     if(allocated(gam_map)) deallocate(gam_map)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!     return
!--------------------SYMMETRIZE----------------------------------------(
         write(*,*)
         write(*,*)'FINISHED CALCULATING CONTRACTED GAM_3 INTEGRALS'
         write(*,*)'BEGINNING SYMMETRIZATION OF GAM_3 INTEGRALS'
         write(*,*)

         open(IFIL,file=FLNM,form='unformatted',
     x    status='unknown',access='direct',RECL=8)

         open(905,file='GAM_3X.ufm',form='unformatted',
     x    status='unknown',access='direct',RECL=8)

      do ip=1,npbf
      do jp=1,npbf
         do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
            do jec2=1,nebf
               do iec3=1,nebf
               do jec3=1,nebf

C  GAM_3 Symmetrization:
C  Determine packing indices for XGAM_3 integral matrices

C              As Packed-->       XGAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) 
      call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
                             ia_123=ia
                             read(905,REC=ia_123) x123
                             ! x123=GAM_3(ia_123)

C              As Packed-->       XGAM_3(je2,ie2,je3,ie3,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2) 
      call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec1,jec1,iec3,jec3,iec2,jec2,ia)
                             ia_132=ia
                             read(905,REC=ia_132) x132
                             ! x132=GAM_3(ia_132)

C              As Packed-->       XGAM_3(je3,ie3,je1,ie1,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3) 
      call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec2,jec2,iec1,jec1,iec3,jec3,ia)
                             ia_213=ia
                             read(905,REC=ia_213) x213
                             ! x213=GAM_3(ia_213)

C              As Packed-->       XGAM_3(je1,ie1,je3,ie3,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1) 
      call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec2,jec2,iec3,jec3,iec1,jec1,ia)
                             ia_231=ia
                             read(905,REC=ia_231) x231
                             ! x231=GAM_3(ia_231)

C              As Packed-->       XGAM_3(je2,ie2,je1,ie1,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2) 
      call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec3,jec3,iec1,jec1,iec2,jec2,ia)
                             ia_312=ia
                             read(905,REC=ia_312) x312
                             ! x312=GAM_3(ia_312)

C              As Packed-->       XGAM_3(je1,ie1,je2,ie2,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1) 
      call index_GAM_3PK(nebf,npbf,
     x                   ip,jp,iec3,jec3,iec2,jec2,iec1,jec1,ia)
                             ia_321=ia
                             read(905,REC=ia_321) x321
                             ! x321=GAM_3(ia_321)

                       xxxx=(x123+x132+x213+x231+x312+x321)/six


c                      write(*,*)'xxxx=',xxxx
cc                     GAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3)=xxxx 
c                      GAM_3PK(ia_123)=xxxx 
c                      call put_GAM3(ia_123,ng3,xxxx)
                       write(IFIL,REC=ia_123) xxxx

cc                     GAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2)=xxxx
c                      GAM_3PK(ia_132)=xxxx 
c                      call put_GAM3(ia_132,ng3,xxxx)
                       write(IFIL,REC=ia_132) xxxx

cc                     GAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3)=xxxx
c                      GAM_3PK(ia_213)=xxxx 
c                      call put_GAM3(ia_213,ng3,xxxx)
                       write(IFIL,REC=ia_213) xxxx

cc                     GAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1)=xxxx
c                      GAM_3PK(ia_231)=xxxx 
c                      call put_GAM3(ia_231,ng3,xxxx)
                       write(IFIL,REC=ia_231) xxxx

cc                     GAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2)=xxxx
c                      GAM_3PK(ia_312)=xxxx 
c                      call put_GAM3(ia_312,ng3,xxxx)
                       write(IFIL,REC=ia_312) xxxx

cc                     GAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1)=xxxx
c                      GAM_3PK(ia_321)=xxxx 
c                      call put_GAM3(ia_321,ng3,xxxx)
                       write(IFIL,REC=ia_321) xxxx


               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      close(905)
      close(IFIL)
!--------------------SYMMETRIZE----------------------------------------)

      return
      end

C=======================================================================
      subroutine thread_gam3_DC(LgVEE,LgTEg,LgVECg,
     x                          LgVEPg,LgVEEg1,LgVEEg2,
     x                          istart,iend,ng3_seg,ng3,
     x                          nebf,npebf,npbf,nat,ngtg1,
     x                          pmass,cat,zan,bcoef1,gamma1,
     x                          loop_map,GAM_3,
     x                          KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                          ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

C=======================================================================
      implicit none
      include 'omp_lib.h'

C Input Variables
      logical LgVEE,LgTEg,LgVECg,LgVEPg,LgVEEg1,LgVEEg2
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
      double precision GAM_3(ng3_seg)

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
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(GAM_3)
!$ompx shared(LgVEE,LgTEg,LgVECg,LgVEPg,LgVEEg1,LgVEEg2)
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
!$ompx private(OMG3)
!$ompx private(id)

!     id= omp_get_thread_num()
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

         call contract_omega3_DC(LgVEE,LgTEg,LgVECg,
     x                             LgVEPg,LgVEEg1,LgVEEg2,
     x                             ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,
     x                             nebf,npebf,npbf,nat,ngtg1,
     x                             pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            OMG3)

         GAM_3(imap)=OMG3

      end do
!$omp end do
!$omp end parallel      
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

C---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime() - wtime
      write(*,*)'TIME TO CALCULATE GAM_3 INTEGRALS: ',wtime
!     write(*,*)'ISTART=',istart
!     write(*,*)'IEND  =',iend
C---OPENMP-TIMING------------------------------------------------------)
      write(*,*)'WRITING THEM TO DISK...'
      write(*,*)

      open(905,file='GAM_3X.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

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

C               As Packed--> XGAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
c                    XGAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) 
         call index_GAM_3PK(nebf,npbf,
     x                      ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,ia)
         write(905,REC=ia) GAM_3(imap)

      end do

      close(905)
!     close(805)


      return
      end

C======================================================================
      subroutine contract_omega3_DC(LgVEE,LgTEg,LgVECg,
     x                            LgVEPg,LgVEEg1,LgVEEg2,
     x                            ip,jp,iec1,jec1,iec2,jec2,iec3,jec3,
     x                            nebf,npebf,npbf,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            OMG3)

C======================================================================
      implicit none
c     include 'omp_lib.h'
c     include 'mpif.h'

C Input Variables
      logical LgVEE,LgTEg,LgVECg,LgVEPg,LgVEEg1,LgVEEg2
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

      OMG3=0.0d+00

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
              call xcalc_GAM3_DC(LgVEE,LgTEg,LgVECg,
     x                           LgVEPg,LgVEEg1,LgVEEg2,
     x                           I1,J1,K1,A1,Amat1,
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

!                       call underflow(ans)

                        OMG3=OMG3+ans
     x                      *Cof_ip*Cof_jp
     x                      *Cof_ie1*Cof_je1
     x                      *Cof_ie2*Cof_je2
     x                      *Cof_ie3*Cof_je3
C---------------------OMG_123------------------------------------------)

          end do
          end do
        end do
        end do
      end do
      end do


      return
      end

