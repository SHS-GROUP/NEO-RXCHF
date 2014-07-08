!======================================================================
      subroutine RXCHF_scf_MPI(nproc,rank,
     x                         nelec,NAE,NAE_alp,NAE_beta,
     x                         NBE,NEBFLT,NPBFLT,NUCST,
     x                         npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                         ngtg1,ng1,ng2,ng3,ng4,
     x                         NG2CHK,NG3CHK,NG4CHK,
     x                         read_CE,read_CP,
     x                         LG4DSCF,LG3DSCF,LG2DSCF,
     x                         LSOSCF,LDIIS,LSORXCHF,
     x                         OCBSETYP,LCMF,LADDEXCH,
     x                         ng2prm,ng3prm,nat,pmass,cat,zan,
     x                         bcoef1,gamma1,
     x                         KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                         ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                         LG2IC1,dimXCHF2,dimINT2,dimINT2ex,
     x                         XCHF_GAM2,INT_GAM2,INT_GAM2ex,XCHF_GAM2s,
     x                         LG3IC1,dimXCHF3,dimINT3,dimINT3ex,
     x                         XCHF_GAM3,INT_GAM3,
     x                         INT_GAM3ex1,INT_GAM3ex2,
     x                         LG4IC,dimXCHF4,dimINT4,
     x                         XCHF_GAM4,INT_GAM4)

!
! PERFORM A NUCLEAR-ELECTRONIC RESTRICTED XC HARTREE-FOCK CALCULATION
! FOR AN NAE-REGULAR ELECTRON NBE-SPECIAL ELECTRON ONE-PROTON SYSTEM
!
!     **DEFINITIONS:
!
!     *FOR REGULAR ELECTRONS:
!     NAE    ::  NUMBER OF REGULAR ELECTRONS
!     DAE    ::  NEW REGULAR ELECTRON DENSITY MATRIX
!     DAE0   ::  OLD REGULAR ELECTRON DENSITY MATRIX
!     VECAE  ::  REGULAR ELECTRON MOS
!     VECAE0 ::  OLD REGULAR ELECTRON MOS
!     AEE    ::  REGULAR ELECTRON ORBITAL EIGENVALUES
!     NAE_alp    ::  NUMBER OF ALPHA REGULAR ELECTRONS
!     DAE_alp    ::  NEW REGULAR ALPHA ELECTRON DENSITY MATRIX
!     DAE0_alp   ::  OLD REGULAR ALPHA ELECTRON DENSITY MATRIX
!     VECAE_alp  ::  REGULAR ALPHA ELECTRON MOS
!     VECAE0_alp ::  OLD REGULAR ALPHA ELECTRON MOS
!     AEE_alp    ::  REGULAR ALPHA ELECTRON ORBITAL EIGENVALUES
!     NAE_beta    ::  NUMBER OF BETA REGULAR ELECTRONS
!     DAE_beta    ::  NEW REGULAR BETA ELECTRON DENSITY MATRIX
!     DAE0_beta   ::  OLD REGULAR BETA ELECTRON DENSITY MATRIX
!     VECAE_beta  ::  REGULAR BETA ELECTRON MOS
!     VECAE0_beta ::  OLD REGULAR BETA ELECTRON MOS
!     AEE_beta    ::  REGULAR BETA ELECTRON ORBITAL EIGENVALUES

!     * FOR SPECIAL ELECTRONS:
!     NBE    ::  NUMBER OF SPECIAL ELECTRONS
!     DBE    ::  NEW SPECIAL ELECTRON DENSITY MATRIX
!     DBE0   ::  OLD SPECIAL ELECTRON DENSITY MATRIX
!     VECBE  ::  SPECIAL ELECTRON MOS
!     VECBE0 ::  OLD SPECIAL ELECTRON MOS
!     BEE    ::  SPECIAL ELECTRON ORBITAL EIGENVALUES
!
!     *FOR PROTONS:
!     DP    ::  NEW PROTON DENSITY MATRIX
!     DP0   ::  OLD PROTON DENSITY MATRIX
!     VECP  ::  PROTON MOS
!     EP    ::  PROTON ORBITAL EIGENVALUES
!
!======================================================================
      implicit none
      include 'mpif.h'

! Input Variables
      integer nproc,rank
      logical LNEOHF
      integer OCBSETYP ! =1 : OCBSE ; =2 : OCBSE2 ; =3 : OCBSE3
      logical LOCBSE   ! Use OCBSE scheme as is (restricted variational freedom for reg/sp elecs)
      logical LOCBSE2  ! Use modified OCBSE scheme (complete variational freedom for reg elecs)
      logical read_CE
      logical read_CP
      logical LGAM4
      logical LCMF
      logical LADDEXCH
      logical LG4DSCF
      logical LG3DSCF
      logical LG2IC1
      logical LG3IC1
      logical LG4IC
      logical LG2DSCF
      integer SZG2ICR
      integer SZG3IC1
      integer SZG4IC
      integer NG2CHK,NG3CHK,NG4CHK
      integer nelec
      integer NAE,NBE
      integer NAE_alp,NAE_beta
      integer NPRA,NPRB
      integer NPRA_alp, NPRA_beta
      integer NEBFLT
      integer NUCST
      integer nebf
      integer npbf
      integer nebf2
      integer npbf2
      integer ngee
      integer ng1             ! Number of GAM1 integrals
      integer ng2             ! Number of GAM2 integrals
      integer ng3             ! Number of GAM3 integrals
      integer ng4             ! Number of GAM3 integrals
!-----DIRECT-SCF-RELATED-----------------------------------------------(
      integer ngtg1
      integer npebf
      integer ng2prm,ng3prm,nat
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
!-----DIRECT-SCF-RELATED-----------------------------------------------)

      integer          dimXCHF2,dimXCHF3,dimXCHF4
      integer          dimINT2,dimINT3,dimINT4
      integer          dimINT2ex,dimINT3ex
      double precision XCHF_GAM2(dimXCHF2)         ! XCHF GAM2 integrals
      double precision INT_GAM2(dimINT2)           ! Interaction GAM2 integrals
      double precision INT_GAM2ex(dimINT2)         ! Exchange GAM2 integrals
      double precision XCHF_GAM2s(dimXCHF2)        ! XCHF GAM2s integrals
      double precision XCHF_GAM3(dimXCHF3)         ! XCHF GAM3 integrals
      double precision INT_GAM3(dimINT3)           ! Interaction GAM3 integrals
      double precision INT_GAM3ex1(dimINT3)        ! Exchange GAM3 integrals
      double precision INT_GAM3ex2(dimINT3)        ! Exchange GAM3 integrals
      double precision XCHF_GAM4(dimXCHF4)         ! XCHF GAM4 integrals
      double precision INT_GAM4(dimINT4)           ! Interaction GAM4 integrals

! Local variables
      double precision zero,one
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00) 

      double precision xxse(nebf,nebf)  ! Elec overlap matrix
      double precision xxsp(npbf,npbf)  ! Nuc overlap matrix
      double precision GAM_pcore(npbf2)
      double precision GAM_ecore(nebf2)
      double precision GAM_ep(ng1)
      double precision GAM_ee(ngee)

!     integer noccE ! Number of occupied elec orbs
!     integer noccP ! Number of occupied nuc orbs
      integer maxit,maxmicroit

      integer i,ielec
      integer j
      integer k
      integer l

      integer NPBFLT

      double precision TOLP
      double precision TOLE
      double precision diffE
      double precision diffAE
      double precision diffAE_alp 
      double precision diffAE_beta
      double precision diffBE
      double precision diffP

      double precision E_total

      double precision E_HF
      double precision E_HF_ecore
      double precision E_HF_ee

      double precision E_XCHF
      double precision E_XCHF_gam1
      double precision E_XCHF_gam2
      double precision E_XCHF_gam3
      double precision E_XCHF_gam4

      double precision E_int
      double precision E_int_OMG2
      double precision E_int_OMG3
      double precision E_int_OMG4

      double precision E_nuc

      double precision S_total
      double precision S_gam1
      double precision S_gam2

      double precision DAE(NEBF,NEBF)
      double precision DAE0(NEBF,NEBF)
      double precision DAE_alp(NEBF,NEBF)
      double precision DAE0_alp(NEBF,NEBF)
      double precision DAE_beta(NEBF,NEBF)
      double precision DAE0_beta(NEBF,NEBF)
      double precision DBE(NEBF,NEBF)
      double precision DBE0(NEBF,NEBF)
      double precision VECAE(NEBF,NEBF)
      double precision VECAE0(NEBF,NEBF)
      double precision VECAE_alp(NEBF,NEBF)
      double precision VECAE0_alp(NEBF,NEBF)
      double precision VECAE_beta(NEBF,NEBF)
      double precision VECAE0_beta(NEBF,NEBF)
      double precision VECBE(NEBF,NEBF)
      double precision VECBE0(NEBF,NEBF)
      double precision AEE(NEBF)
      double precision AEE_alp(NEBF)
      double precision AEE_beta(NEBF)
      double precision BEE(NEBF)
      double precision FAE(nebf,nebf)
      double precision XFAE(nebf,nebf)
      double precision FAE_alp(nebf,nebf)
      double precision XFAE_alp(nebf,nebf)
      double precision FAE_beta(nebf,nebf)
      double precision XFAE_beta(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision XFBE(nebf,nebf)

      double precision DP(NPBF,NPBF)
      double precision DP0(NPBF,NPBF)
      double precision VECP(NPBF,NPBF)
      double precision EP(NPBF)
      double precision FP(npbf,npbf)
      double precision XFP(npbf,npbf)

      double precision FAEint(nebf,nebf)
      double precision FBEint(nebf,nebf)
      double precision FPint(npbf,npbf)

      double precision SBE_XCHF(nebf,nebf)
      double precision SP_XCHF(npbf,npbf)

      double precision E_total_old
      double precision Delta_E_tot

      logical LDIFFE

!--------SOSCF-RELATED-VARIABLES------------(
      logical LSOSCF,LSOSCFA,LSOSCFB
      logical EIGAVL
      integer NA
      integer NA_alp
      integer NA_beta
      integer NB
      integer ITER
      integer ITSOA ! SOSCF iteration counter
      integer ITSOA_alp
      integer ITSOA_beta
      integer ITSOB ! SOSCF iteration counter
      integer L0,L1
      integer L0wA,L1wA
      integer L0wA_alp,L1wA_alp,L0wA_beta,L1wA_beta
      integer L0wB,L1wB
      integer NFT15
      integer NFT16
      integer NFT17
      double precision FLT(NEBFLT) !FLT: Lower triangle focke
      double precision, allocatable :: HSTARTA(:)
      double precision, allocatable :: GRADA(:)
      double precision, allocatable :: PGRADA(:)
      double precision, allocatable :: DISPLIA(:)
      double precision, allocatable :: DGRADA(:)  ! WRK1
      double precision, allocatable :: DISPLA(:)  ! WRK2
      double precision, allocatable :: UPDTA(:)   ! WRK3
      double precision, allocatable :: DISPLNA(:) ! WRK1+NPR
      double precision, allocatable :: DGRADIA(:) ! WRK2+NPR
      double precision, allocatable :: UPDTIA(:)  ! WRK3+NPR
      double precision, allocatable :: XA(:)
      double precision ORBGRDA
      double precision GA(nebf,nebf) !G(L0,L0)
      double precision, allocatable :: HSTARTA_alp(:)
      double precision, allocatable :: GRADA_alp(:)
      double precision, allocatable :: PGRADA_alp(:)
      double precision, allocatable :: DISPLIA_alp(:)
      double precision, allocatable :: DGRADA_alp(:)  ! WRK1
      double precision, allocatable :: DISPLA_alp(:)  ! WRK2
      double precision, allocatable :: UPDTA_alp(:)   ! WRK3
      double precision, allocatable :: DISPLNA_alp(:) ! WRK1+NPR
      double precision, allocatable :: DGRADIA_alp(:) ! WRK2+NPR
      double precision, allocatable :: UPDTIA_alp(:)  ! WRK3+NPR
      double precision, allocatable :: XA_alp(:)
      double precision ORBGRDA_alp
      double precision GA_alp(nebf,nebf) !G(L0,L0)
      double precision, allocatable :: HSTARTA_beta(:)
      double precision, allocatable :: GRADA_beta(:)
      double precision, allocatable :: PGRADA_beta(:)
      double precision, allocatable :: DISPLIA_beta(:)
      double precision, allocatable :: DGRADA_beta(:)  ! WRK1
      double precision, allocatable :: DISPLA_beta(:)  ! WRK2
      double precision, allocatable :: UPDTA_beta(:)   ! WRK3
      double precision, allocatable :: DISPLNA_beta(:) ! WRK1+NPR
      double precision, allocatable :: DGRADIA_beta(:) ! WRK2+NPR
      double precision, allocatable :: UPDTIA_beta(:)  ! WRK3+NPR
      double precision, allocatable :: XA_beta(:)
      double precision ORBGRDA_beta
      double precision GA_beta(nebf,nebf) !G(L0,L0)
      double precision WRK(nebf) !WRK(L0)
      double precision, allocatable :: HSTARTB(:)
      double precision, allocatable :: GRADB(:)
      double precision, allocatable :: PGRADB(:)
      double precision, allocatable :: DISPLIB(:)
      double precision, allocatable :: DGRADB(:)  ! WRK1
      double precision, allocatable :: DISPLB(:)  ! WRK2
      double precision, allocatable :: UPDTB(:)   ! WRK3
      double precision, allocatable :: DISPLNB(:) ! WRK1+NPR
      double precision, allocatable :: DGRADIB(:) ! WRK2+NPR
      double precision, allocatable :: UPDTIB(:)  ! WRK3+NPR
      double precision, allocatable :: XB(:)
      double precision GB(nebf,nebf) !G(L0,L0)
      double precision ORBGRDB
      double precision SMALL
      double precision SOGTOL ! ORBGRAD TOL to activate soscf
!cc   double precision CCC(nebf,nebf) !WRK(L0)
!cc   NPR=(L0-NA)*NA ! Line 2134 RHFCL ?NA is NUM ALPHA E?
!--------SOSCF-RELATED-VARIABLES------------)
C ARS( OCBSE variables
      integer nocca,noccb
      integer nwbfA,nwbfB                                          ! # bf in W basis
      integer nwbfltA,nwbfltB                                      ! one-dimensional size
      double precision, allocatable :: tempvecA(:,:),tempvecB(:,:) ! transformation matrix
      double precision, allocatable :: WA(:,:),WB(:,:)             ! transformation matrix
      double precision, allocatable :: wFAEw(:,:),wFBEw(:,:)       ! Fock matrix in W basis
      double precision, allocatable :: wvecAEw(:,:),wvecBEw(:,:)   ! eigenvectors in W basis from curr it
      double precision, allocatable :: wAEenw(:),wBEenw(:)         ! eigenvalues
      double precision, allocatable :: wFLTAw(:),wFLTBw(:)         ! flattened wFEw
      double precision, allocatable :: wGAw(:,:),wGBw(:,:)         ! exponential transformation
      double precision, allocatable :: wWRKAw(:),wWRKBw(:)         ! work array (L0w)
      integer nocca_alp,nocca_beta
      integer nwbfA_alp,NwbfA_beta                                          ! # bf in W basis
      integer nwbfltA_alp,nwbfltA_beta                                      ! one-dimensional size
      double precision, allocatable :: tempvecA_alp(:,:),
     x                                 tempvecA_beta(:,:) ! transformation matrix
      double precision, allocatable :: WA_alp(:,:),WA_beta(:,:)             ! transformation matrix
      double precision, allocatable :: wFAEw_alp(:,:),wFAEw_beta(:,:)       ! Fock matrix in W basis
      double precision, allocatable :: wvecAEw_alp(:,:),
     x                                 wvecAEw_beta(:,:)   ! eigenvectors in W basis from curr it
      double precision, allocatable :: wAEenw_alp(:),wAEenw_beta(:)         ! eigenvalues
      double precision, allocatable :: wFLTAw_alp(:),wFLTAw_beta(:)         ! flattened wFEw
      double precision, allocatable :: wGAw_alp(:,:),wGAw_beta(:,:)         ! exponential transformation
      double precision, allocatable :: wWRKAw_alp(:),wWRKAw_beta(:)         ! work array (L0w)
C )
C ARS( DIIS variables
      logical LDIIS
      integer itDIIS
      integer nstore
      parameter(nstore=4)
      double precision errDIIS,threshDIIS
      double precision, allocatable ::  errvec(:,:,:)   ! stored as (AO ind,MO ind,iter) where
      double precision vecDIIS(nebf,nebf,nstore)        ! MOs stored in order {occA,occB,virt}
C )
C ARS( SORXCHF variables
      logical LSORXCHF
      double precision, allocatable :: HESSA(:,:)       ! Hessian at current iteration
      double precision, allocatable :: HESSB(:,:)
      double precision, allocatable :: HESSA0(:,:)      ! Hessian at previous iteration
      double precision, allocatable :: HESSB0(:,:)
!     double precision, allocatable :: GRADA(:)         ! Gradient at current iteration
!     double precision, allocatable :: GRADB(:)
      double precision, allocatable :: GRADA0(:)        ! Gradient at previous iteration
      double precision, allocatable :: GRADB0(:)
!     double precision, allocatable :: DISPLIA(:)       ! Displ vec at current iteration
!     double precision, allocatable :: DISPLIB(:)
      double precision, allocatable :: DISPLIA0(:)      ! Displ vec at previous iteration
      double precision, allocatable :: DISPLIB0(:)
      double precision, allocatable :: WA0(:,:)         ! W bases at previous iteration
      double precision, allocatable :: WB0(:,:)
      double precision, allocatable :: HESSA_alp(:,:)       ! Hessian at current iteration
      double precision, allocatable :: HESSA_beta(:,:)
      double precision, allocatable :: HESSA0_alp(:,:)      ! Hessian at previous iteration
      double precision, allocatable :: HESSA0_beta(:,:)
!     double precision, allocatable :: GRADA(:)         ! Gradient at current iteration
!     double precision, allocatable :: GRADB(:)
      double precision, allocatable :: GRADA0_alp(:)        ! Gradient at previous iteration
      double precision, allocatable :: GRADA0_beta(:)
!     double precision, allocatable :: DISPLIA(:)       ! Displ vec at current iteration
!     double precision, allocatable :: DISPLIB(:)
      double precision, allocatable :: DISPLIA0_alp(:)      ! Displ vec at previous iteration
      double precision, allocatable :: DISPLIA0_beta(:)
      double precision, allocatable :: WA0_alp(:,:)         ! W bases at previous iteration
      double precision, allocatable :: WA0_beta(:,:)
C )C )
C ARS( testing variables
      logical LNOINT
      double precision FBEmo(nebf,nebf)
      double precision work1(nebf,nebf)
      double precision work2(nebf,nebf)
C )
#if MPI32
      integer*4 ierr
      integer*4 npbf_,npbf2_,nebf_,nebf2_,ng1_,ngee_
#else
      integer ierr
      integer npbf_,npbf2_,nebf_,nebf2_,ng1_,ngee_
#endif

!--------OUTPUT-FORMATTING---------------------------------------------(
 9000 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS ')

 9050 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS         ',
     * 'ORBGRAD_A ')

 9051 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS         ',
     * 'ORBGRAD_B ')

 9052 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS         ',
     * 'ORBGRAD_A        ORBGRAD_B ')

 9053 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS         ',
     * 'DIIS ERR  ')

 9054 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS         ',
     * 'ORBGRAD_A        ORBGRAD_B        DIIS ERR  ')
 9055 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA_alp DENS   ALPHA_beta DENS  BETA DENS        ',
     *  'QMP DENS         ',
     * 'ORBGRAD_A_alp    ORBGRAD_A_beta   ORBGRAD_B ')

 9100 FORMAT(1X,I3,F20.10,F17.10,3F17.10)

 9150 FORMAT(1X,I3,F20.10,F17.10,4F17.10)

 9151 FORMAT(1X,I3,F20.10,F17.10,5F17.10)

 9152 FORMAT(1X,I3,F20.10,F17.10,6F17.10)
 
 9153 FORMAT(1X,I3,F20.10,F17.10,7F17.10)

 9200 FORMAT(/1X,'FINAL NEORXCHF ENERGY IS',F20.10,' AFTER',I4,
     *           ' ITERATIONS')

 9300 FORMAT(/6X,'-----------------------------------------------',/,
     x        6X,'         NEORXCHF ENERGETIC COMPONENTS         ',/,
     x        6X,'-----------------------------------------------',/,
     x       12X,'    E_NUC=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,'  E_ECORE=',1X,F20.10/
     x       12X,'     E_EE=',1X,F20.10/
     x       12X,'     E_HF=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,'   E_GAM1=',1X,F20.10/
     x       12X,'   E_GAM2=',1X,F20.10/
     x       12X,'   E_GAM3=',1X,F20.10/
     x       12X,'   E_GAM4=',1X,F20.10/
     x       12X,'   E_XCHF=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,'   E_OMG2=',1X,F20.10/
     x       12X,'   E_OMG3=',1X,F20.10/
     x       12X,'   E_OMG4=',1X,F20.10/
     x       12X,'    E_int=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,'  S_TOTAL=',1X,F20.10/
     x       12X,'  E_TOTAL=',1X,F20.10/
     x        6X,'-----------------------------------------------',/)

 9400 FORMAT(/1X,'          INITIAL GUESS ENERGETICS:            ')

 9500 FORMAT(/6X,' ** BEGIN RXCHF SELF-CONSISTENT-FIELD CALCULATION **')

 9610 FORMAT(/1X,' REGULAR ELECTRONIC ORBITALS AND EIGENVALUES:       ')

 9620 FORMAT(/1X,' SPECIAL ELECTRONIC ORBITALS AND EIGENVALUES:       ')

 9700 FORMAT(/1X,'      QM PARTICLE ORBITALS AND EIGENVALUES:         ')

 9800 FORMAT(10X,15(1H-),'START SECOND ORDER SCF',15(1H-))

 9801 FORMAT(10X,15(1H-),' START DIIS PROCEDURE ',15(1H-))

 2001 FORMAT(/1X,'STARTING MICROITERATIONS FOR ITERATION',1X,I3)

 2000 FORMAT(1X,'CONVERGED ITERATION',1X,I3,1X,'IN',
     x       1X,I3,1X,'MICROITERATIONS',/)
                                           
!--------OUTPUT-FORMATTING---------------------------------------------)

#if MPI32
      npbf_=int(npbf,kind=4)
      npbf2_=int(npbf2,kind=4)
      nebf_=int(nebf,kind=4)
      nebf2_=int(nebf2,kind=4)
      ng1_=int(ng1,kind=4)
      ngee_=int(ngee,kind=4)
#else
      npbf_=npbf
      npbf2_=npbf2
      nebf_=nebf
      nebf2_=nebf2
      ng1_=ng1
      ngee_=ngee
#endif

      LOCBSE=.false.
      LOCBSE2=.false.
      if(OCBSETYP.eq.1) then
       LOCBSE=.true.
      else if (OCBSETYP.eq.2) then
       LOCBSE2=.true.
      end if
cKRB this is changed for UHF for specials
c      if (nae.gt.1) then
c       nocca=nae/2
c      else
       nocca=nae_alp
       nocca_alp =nae_alp
       nocca_beta=nae_beta
c      end if
      if (nbe.gt.1) then
       noccb=nbe/2
      else
       noccb=nbe
      end if

      if((LDIIS).and.(rank.eq.0)) then
       if(allocated(errvec)) deallocate(errvec)
       allocate(errvec(nebf,nebf,nstore))
       threshDIIS=0.5d+00
      end if

      if(LOCBSE.or.LOCBSE2) then

       nwbfB=nebf-nocca
       nwbfltB=nwbfB*(nwbfB+1)/2
       L0wB=nwbfB
       L1wB=nwbfB
       if (rank.eq.0) then
        if(allocated(tempvecB)) deallocate(tempvecB)
        allocate(tempvecB(nebf,nwbfB))
        if(allocated(WB)) deallocate(WB)
        allocate(WB(nebf,nwbfB))
        if(allocated(wFBEw)) deallocate(wFBEw)
        allocate(wFBEw(nwbfB,nwbfB))
        if(allocated(wvecBEw)) deallocate(wvecBEw)
        allocate(wvecBEw(nwbfB,nwbfB))
        if(allocated(wBEenw)) deallocate(wBEenw)
        allocate(wBEenw(nwbfB))
        if(allocated(wFLTBw)) deallocate(wFLTBw)
        allocate(wFLTBw(nwbfltB))
        if(allocated(wGBw)) deallocate(wGBw)
        allocate(wGBw(nwbfB,nwbfB))
        if(allocated(wWRKBw)) deallocate(wWRKBw)
        allocate(wWRKBw(nwbfB))
       end if

       if(LOCBSE2) then

        if(rank.eq.0) write(*,*) "Using LOCBSE2"
        if(LSORXCHF) then
         SOGTOL=0.40d+00
         npra=(nebf-nocca)*nocca       ! occ-vir pairs for regular elecs
         nprb=(nebf-noccb-nocca)*noccb ! occ-vir pairs for regular elecs
         NA=nocca
         NB=noccb
         L0=nebf
         L1=nebf
         if (rank.eq.0) then
          if(allocated(WB0))      deallocate(WB0)
          if(allocated(HSTARTA))  deallocate(HSTARTA)
          if(allocated(HSTARTB))  deallocate(HSTARTB)
          if(allocated(HESSA))    deallocate(HESSA)
          if(allocated(HESSB))    deallocate(HESSB)
          if(allocated(HESSA0))   deallocate(HESSA0)
          if(allocated(HESSB0))   deallocate(HESSB0)
          if(allocated(GRADA))    deallocate(GRADA)
          if(allocated(GRADB))    deallocate(GRADB)
          if(allocated(GRADA0))   deallocate(GRADA0)
          if(allocated(GRADB0))   deallocate(GRADB0)
          if(allocated(DISPLIA))  deallocate(DISPLIA)
          if(allocated(DISPLIB))  deallocate(DISPLIB)
          if(allocated(DISPLIA0)) deallocate(DISPLIA0)
          if(allocated(DISPLIB0)) deallocate(DISPLIB0)
          allocate(WB0(nebf,nwbfB))
          allocate(HSTARTA(NPRA))
          allocate(HSTARTB(NPRB))
          allocate(HESSA(NPRA,NPRA))
          allocate(HESSB(NPRB,NPRB))
          allocate(HESSA0(NPRA,NPRA))
          allocate(HESSB0(NPRB,NPRB))
          allocate(GRADA(NPRA))
          allocate(GRADB(NPRB))
          allocate(GRADA0(NPRA))
          allocate(GRADB0(NPRB))
          allocate(DISPLIA(NPRA))
          allocate(DISPLIB(NPRB))
          allocate(DISPLIA0(NPRA))
          allocate(DISPLIB0(NPRB))
         end if
        end if

       else
        if(rank.eq.0) write(*,*) "Using LOCBSE"
        nwbfA_alp=nebf-noccb
        nwbfltA_alp=nwbfA_alp*(nwbfA_alp+1)/2
        L0wA_alp=nwbfA_alp
        L1wA_alp=nwbfA_alp
        nwbfA_beta=nebf-noccb
        nwbfltA_beta=nwbfA_beta*(nwbfA_beta+1)/2
        L0wA_beta=nwbfA_beta
        L1wA_beta=nwbfA_beta
        nwbfA=nebf-noccb
        nwbfltA=nwbfA*(nwbfA+1)/2
        L0wA=nwbfA
        L1wA=nwbfA
        if (rank.eq.0) then
         if(allocated(tempvecA)) deallocate(tempvecA)
         allocate(tempvecA(nebf,nwbfA))
         if(allocated(WA)) deallocate(WA)
         allocate(WA(nebf,nwbfA))
         if(allocated(wFAEw)) deallocate(wFAEw)
         allocate(wFAEw(nwbfA,nwbfA))
         if(allocated(wvecAEw)) deallocate(wvecAEw)
         allocate(wvecAEw(nwbfA,nwbfA))
         if(allocated(wAEenw)) deallocate(wAEenw)
         allocate(wAEenw(nwbfA))
         if(allocated(wFLTAw)) deallocate(wFLTAw)
         allocate(wFLTAw(nwbfltA))
         if(allocated(wGAw)) deallocate(wGAw)
         allocate(wGAw(nwbfA,nwbfA))
         if(allocated(wWRKAw)) deallocate(wWRKAw)
         allocate(wWRKAw(nwbfA))

        if(allocated(tempvecA_alp)) deallocate(tempvecA_alp)
        allocate(tempvecA_alp(nebf,nwbfA_alp))
        if(allocated(WA_alp)) deallocate(WA_alp)
        allocate(WA_alp(nebf,nwbfA_alp))
        if(allocated(wFAEw_alp)) deallocate(wFAEw_alp)
        allocate(wFAEw_alp(nwbfA_alp,nwbfA_alp))
        if(allocated(wvecAEw_alp)) deallocate(wvecAEw_alp)
        allocate(wvecAEw_alp(nwbfA_alp,nwbfA_alp))
        if(allocated(wAEenw_alp)) deallocate(wAEenw_alp)
        allocate(wAEenw_alp(nwbfA_alp))
        if(allocated(wFLTAw_alp)) deallocate(wFLTAw_alp)
        allocate(wFLTAw_alp(nwbfltA_alp))
        if(allocated(wGAw_alp)) deallocate(wGAw_alp)
        allocate(wGAw_alp(nwbfA_alp,nwbfA_alp))
        if(allocated(wWRKAw_alp)) deallocate(wWRKAw_alp)
        allocate(wWRKAw_alp(nwbfA_alp))

        if(allocated(tempvecA_beta)) deallocate(tempvecA_beta)
        allocate(tempvecA_beta(nebf,nwbfA_beta))
        if(allocated(WA_beta)) deallocate(WA_beta)
        allocate(WA_beta(nebf,nwbfA_beta))
        if(allocated(wFAEw_beta)) deallocate(wFAEw_beta)
        allocate(wFAEw_beta(nwbfA_beta,nwbfA_beta))
        if(allocated(wvecAEw_beta)) deallocate(wvecAEw_beta)
        allocate(wvecAEw_beta(nwbfA,nwbfA_beta))
        if(allocated(wAEenw_beta)) deallocate(wAEenw_beta)
        allocate(wAEenw_beta(nwbfA_beta))
        if(allocated(wFLTAw_beta)) deallocate(wFLTAw_beta)
        allocate(wFLTAw_beta(nwbfltA_beta))
        if(allocated(wGAw_beta)) deallocate(wGAw_beta)
        allocate(wGAw_beta(nwbfA_beta,nwbfA_beta))
        if(allocated(wWRKAw_beta)) deallocate(wWRKAw_beta)
        allocate(wWRKAw_beta(nwbfA_beta))
        
        end if
        if(LSORXCHF) then
         SOGTOL=0.40d+00
         npra=(nebf-nocca-noccb)*nocca ! occ-vir pairs for regular elecs
         npra_alp=(nebf-nocca_alp-noccb)*nocca_alp ! occ-vir pairs for alpha regular elecs
         npra_beta=(nebf-nocca_beta-noccb)*nocca_beta ! occ-vir pairs for beta regular elecs
         nprb=(nebf-noccb-nocca)*noccb ! occ-vir pairs for special elecs
         NA=nocca
         NA_alp = nocca_alp
         NA_beta = nocca_beta
         NB=noccb
         if (rank.eq.0) then
          if(allocated(WA0))      deallocate(WA0)
          if(allocated(WB0))      deallocate(WB0)
          if(allocated(HSTARTA))  deallocate(HSTARTA)
          if(allocated(HSTARTB))  deallocate(HSTARTB)
          if(allocated(HESSA))    deallocate(HESSA)
          if(allocated(HESSB))    deallocate(HESSB)
          if(allocated(HESSA0))   deallocate(HESSA0)
          if(allocated(HESSB0))   deallocate(HESSB0)
          if(allocated(GRADA))    deallocate(GRADA)
          if(allocated(GRADB))    deallocate(GRADB)
          if(allocated(GRADA0))   deallocate(GRADA0)
          if(allocated(GRADB0))   deallocate(GRADB0)
          if(allocated(DISPLIA))  deallocate(DISPLIA)
          if(allocated(DISPLIB))  deallocate(DISPLIB)
          if(allocated(DISPLIA0)) deallocate(DISPLIA0)
          if(allocated(DISPLIB0)) deallocate(DISPLIB0)
          allocate(WA0(nebf,nwbfA))
          allocate(WB0(nebf,nwbfB))
          allocate(HSTARTA(NPRA))
          allocate(HSTARTB(NPRB))
          allocate(HESSA(NPRA,NPRA))
          allocate(HESSB(NPRB,NPRB))
          allocate(HESSA0(NPRA,NPRA))
          allocate(HESSB0(NPRB,NPRB))
          allocate(GRADA(NPRA))
          allocate(GRADB(NPRB))
          allocate(GRADA0(NPRA))
          allocate(GRADB0(NPRB))
          allocate(DISPLIA(NPRA))
          allocate(DISPLIB(NPRB))
          allocate(DISPLIA0(NPRA))
          allocate(DISPLIB0(NPRB))
          
          if(allocated(WA0_alp))      deallocate(WA0_alp)
          if(allocated(WA0_beta))      deallocate(WA0_beta)
          if(allocated(HSTARTA_alp))  deallocate(HSTARTA_alp)
          if(allocated(HSTARTA_beta))  deallocate(HSTARTA_beta)
          if(allocated(HESSA_alp))    deallocate(HESSA_alp)
          if(allocated(HESSA_beta))    deallocate(HESSA_beta)
          if(allocated(HESSA0_alp))   deallocate(HESSA0_alp)
          if(allocated(HESSA0_beta))   deallocate(HESSA0_beta)
          if(allocated(GRADA_alp))    deallocate(GRADA_alp)
          if(allocated(GRADA_beta))    deallocate(GRADA_beta)
          if(allocated(GRADA0_alp))   deallocate(GRADA0_alp)
          if(allocated(GRADA0_beta))   deallocate(GRADA0_beta)
          if(allocated(DISPLIA_alp))  deallocate(DISPLIA_alp)
          if(allocated(DISPLIA_beta))  deallocate(DISPLIA_beta)
          if(allocated(DISPLIA0_alp)) deallocate(DISPLIA0_alp)
          if(allocated(DISPLIA0_beta)) deallocate(DISPLIA0_beta)
          allocate(WA0_alp(nebf,nwbfA_alp))
          allocate(WA0_beta(nebf,nwbfA_beta))
          allocate(HSTARTA_alp(NPRA_alp))
          allocate(HSTARTA_beta(NPRA_beta))
          allocate(HESSA_alp(NPRA_alp,NPRA_alp))
          allocate(HESSA_beta(NPRA_beta,NPRA_beta))
          allocate(HESSA0_alp(NPRA_alp,NPRA_alp))
          allocate(HESSA0_beta(NPRA_beta,NPRA_beta))
          allocate(GRADA_alp(NPRA_alp))
          allocate(GRADA_beta(NPRA_beta))
          allocate(GRADA0_alp(NPRA_alp))
          allocate(GRADA0_beta(NPRA_beta))
          allocate(DISPLIA_alp(NPRA_alp))
          allocate(DISPLIA_beta(NPRA_beta))
          allocate(DISPLIA0_alp(NPRA_alp))
          allocate(DISPLIA0_beta(NPRA_beta))
         end if
        end if

       end if

      end if

      LGAM4=.true. ! Always calculate five-particle integrals
C ARS( no interaction
      LNOINT=.false.
C      LNOINT=.true.
      if ((LNOINT).and.(rank.eq.0)) then
       write(*,*)
       write(*,*) "******************"
       write(*,*) "  NO INTERACTION  "
       write(*,*) "******************"
       write(*,*)
      end if
C )

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------(
!      call class_nuc_rep(nat,zan,cat,E_nuc)
      if (rank.eq.0) then
       open(800,file='ENUCRP.dat',status='unknown')
       read(800,*) E_nuc
       close(800)
       write(*,*)'READ IN NUC REPULSION'
      end if
      call MPI_BCAST(E_nuc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------)

!--------------READ-INTEGRALS-NEEDED-FOR-NEO-HF------------------------(
      if (rank.eq.0) then
       call read_nuc_ovlap(npbf,xxsp)
       write(*,*)
       write(*,*)'READ IN NUC OVLAP'
       call read_elec_ovlap(nebf,xxse)
       write(*,*)'READ IN ELEC OVLAP'
       call read_GAM_ecore(nebf,nebf2,GAM_ecore)
       write(*,*)'READ IN GAM_ECORE'
       call read_GAM_pcore(npbf,npbf2,GAM_pcore)
       write(*,*)'READ IN GAM_PCORE'
       call read_GAM_ep(nebf,npbf,ng1,GAM_ep)
       write(*,*)'READ IN GAM_EP'
       call read_GAM_ee(nebf,ngee,GAM_ee)
       write(*,*)'READ IN GAM_EE'
       write(*,*)
      end if
      call MPI_BCAST(xxsp,npbf2_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xxse,nebf2_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(GAM_ecore,nebf2_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(GAM_pcore,npbf2_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(GAM_ep,ng1_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(GAM_ee,ngee_,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!--------------READ-INTEGRALS-NEEDED-FOR-NEO-HF------------------------)

!-------------INITIAL-GUESSES------------------------------------------(
      if(read_CE) then

       if (rank.eq.0) then
!       READ IN GUESS FOR E:
!       call read_elec_density(nebf,nelec,DE)
        call RXCHFmult_read_CAE(nebf,NAE,DAE,VECAE0)
        call RXCHFmult_read_CBE(nebf,NBE,DBE,VECBE0)
       end if

       call MPI_BCAST(DAE,nebf2_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(VECAE0,nebf2_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DBE,nebf2_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(VECBE0,nebf2_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      else
!       STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
        if(rank.eq.0) write(*,*)'ABOUT TO CALL guess_A_elec'
!       call guess_elec(nelec,nebf,xxse,GAM_ecore,DE)
        if ((LOCBSE).or.(LOCBSE2)) then
          call RXCHFmult_guess_elec_uhf(nae,nbe,nebf,xxse,GAM_ecore,
     x                          DAE,DBE,VECAE0,VECBE0,
     x                          nae_alp,nae_beta,DAE_alp,DAE_beta,
     x                           VECAE0_alp,VECAE0_beta)
          if(rank.eq.0) write(*,*)'BACK FROM guess_elec for OCBSE'
        else
         call RXCHFmult_guess_A_elec(NAE,nebf,xxse,GAM_ecore,DAE,VECAE0)
         call RXCHFmult_guess_A_elec(NBE,nebf,xxse,GAM_ecore,DBE,VECBE0)
          if(rank.eq.0) write(*,*)'BACK FROM guess_elec'
        end if
      end if

      if(read_CP) then
!        READ IN GUESS FOR N:

       if (rank.eq.0) then
        call read_nuc_density(npbf,1,NUCST,DP)
       end if

       call MPI_BCAST(DP,npbf2_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      else
!        STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
         if(rank.eq.0) write(*,*)'ABOUT TO CALL guess_prot'
!        call guess_prot(NUCST,npbf,nebf,xxsp,GAM_pcore,GAM_ep,DE,DP)
         call guess_prot2(NUCST,npbf,xxsp,GAM_pcore,DP)
         if(rank.eq.0) write(*,*)'BACK FROM guess_prot'
      end if

C ARS( debug: print out initial guess MOs here
      AEe=0.0d+00
      BEe=0.0d+00
      if ((LCMF).and.(rank.eq.0)) then
       write(*,*)
       write(*,*) "------------------"
       write(*,*) "INITIAL GUESS MOs:"
       write(*,*) "------------------"
       write(*,*)
       WRITE(*,9610)
       WRITE(*,*) "ALPHA ELECTRONS"
       call PREVNU(vecAE0_alp,AEE_alp,nebf,nebf,nebf)
       WRITE(*,9620)
       call PREVNU(vecBE0,BEE,nebf,nebf,nebf)
C vecP not defined yet
C       WRITE(*,9700)
C       call PREVNU(vecp,EP,npbf,npbf,npbf)
       write(*,*)
      end if
C )
!-------------INITIAL-GUESSES------------------------------------------)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------(
c@@@krb think we should not ever be here
      if (LSOSCF) then
         write(*,*) 'krb shouldnt ever be here, since should' 
         write(*,*) 'use SORXCHF'
         call exit
         SOGTOL=0.40d+00
         SMALL=1.0D-06
         L0=nebf
         L1=nebf
         LSOSCFA=.true.
         LSOSCFB=.true.
c@@@krb taking this out right now
         if(nae.eq.1) LSOSCFA=.FALSE.
         if(nbe.eq.1) LSOSCFB=.FALSE.
      else
         LSOSCFA=.false.
         LSOSCFB=.false.
      end if
      if(LSOSCFA) THEN
         NFT15=15
         if (rank.eq.0) then
          OPEN(NFT15, FILE='WORK15', STATUS='UNKNOWN',
     *         ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         end if
         NA=nocca
         if(LOCBSE) then
          npra=(nebf-nocca-noccb)*nocca ! occ-vir pairs for regular elecs
         else
          npra=(nebf-nocca)*nocca ! occ-vir pairs for regular elecs
         end if
! Allocate here (previously on stack)
         if(allocated(XA))      deallocate(XA)
         if(allocated(UPDTIA))  deallocate(UPDTIA)
         if(allocated(DGRADIA)) deallocate(DGRADIA)
         if(allocated(DISPLNA)) deallocate(DISPLNA)
         if(allocated(UPDTA))   deallocate(UPDTA)
         if(allocated(DISPLA))  deallocate(DISPLA)
         if(allocated(DGRADA))  deallocate(DGRADA)
         if(allocated(DISPLIA)) deallocate(DISPLIA)
         if(allocated(PGRADA))  deallocate(PGRADA)
         if(allocated(GRADA))   deallocate(GRADA)
         if(allocated(HSTARTA)) deallocate(HSTARTA)
         allocate(HSTARTA(NPRA))
         allocate(GRADA(NPRA))
         allocate(PGRADA(NPRA))
         allocate(DISPLIA(NPRA))
         allocate(DGRADA(NPRA))
         allocate(DISPLA(NPRA))
         allocate(UPDTA(NPRA))
         allocate(DISPLNA(NPRA))
         allocate(DGRADIA(NPRA))
         allocate(UPDTIA(NPRA))
         allocate(XA(NPRA))
      end if
      if(LSOSCFB) THEN
         NFT16=16
         if (rank.eq.0) then
          OPEN(NFT16, FILE='WORK16', STATUS='UNKNOWN',
     *         ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         end if
         NB=noccb
         if(LOCBSE.or.LOCBSE2) then
          nprb=(nebf-noccb-nocca)*noccb ! occ-vir pairs for special elecs
         else
          nprb=(nebf-noccb)*noccb ! occ-vir pairs for special elecs
         end if
! Allocate here (previously on stack)
         if(allocated(XB))      deallocate(XB)
         if(allocated(UPDTIB))  deallocate(UPDTIB)
         if(allocated(DGRADIB)) deallocate(DGRADIB)
         if(allocated(DISPLNB)) deallocate(DISPLNB)
         if(allocated(UPDTB))   deallocate(UPDTB)
         if(allocated(DISPLB))  deallocate(DISPLB)
         if(allocated(DGRADB))  deallocate(DGRADB)
         if(allocated(DISPLIB)) deallocate(DISPLIB)
         if(allocated(PGRADB))  deallocate(PGRADB)
         if(allocated(GRADB))   deallocate(GRADB)
         if(allocated(HSTARTB)) deallocate(HSTARTB)
         allocate(HSTARTB(NPRB))
         allocate(GRADB(NPRB))
         allocate(PGRADB(NPRB))
         allocate(DISPLIB(NPRB))
         allocate(DGRADB(NPRB))
         allocate(DISPLB(NPRB))
         allocate(UPDTB(NPRB))
         allocate(DISPLNB(NPRB))
         allocate(DGRADIB(NPRB))
         allocate(UPDTIB(NPRB))
         allocate(XB(NPRB))
      end if
!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------)

!     SET CONVERGENCE CRITERIA AND MAXIMUM ITERATIONS 
      TOLE = 1.0D-06
      TOLP = 1.0D-04
      maxit=30
      maxmicroit=100
      if(.not.(LSOSCF.or.LDIIS.or.LSORXCHF)) maxmicroit=200

!     ZERO OUT 'OLD' DENSITY MATRICES
      DAE0=0.0d+00
      DAE0_alp=0.0d+00
      DAE0_beta=0.0d+00
      DBE0=0.0d+00
      DP0=0.0d+00

!     BEGIN XCSCF ITERATIONS
      if(rank.eq.0) WRITE(*,9500)
      E_total_old=0.0d+00

      DO I=1,MAXIT

C Call UHF Fock build for NAE regular electrons
c@@@krb taking this out for uhf for regulars
c         call RXCHF_fock_hf_MPI(nproc,rank,
c     x                          LCMF,nebf,nebf2,NAE,ngee,
c     x                          DAE,GAM_ecore,GAM_ee,
c     x                          FAE,E_HF,E_HF_ecore,E_HF_ee)
         call RXCHF_fock_uhf_MPI(nproc,rank,
     x                          LCMF,nebf,nebf2,NAE,ngee,
     x                          DAE,GAM_ecore,GAM_ee,
     x							FAE,E_HF,E_HF_ecore,E_HF_ee,
     x                           DAE_alp,DAE_beta,FAE_alp,FAE_beta)

C Call XCHF Fock build for NBE special electrons and one QM particle
         call RXCHF_fock_xchf_MPI(nproc,rank,
     x                            LCMF,NBE,
     x                            nebf,nebf2,npbf,npbf2,
     x                            ngee,ng1,ng2,ng3,ng4,
     x                            dimXCHF2,dimXCHF3,dimXCHF4,
     x                            NG2CHK,NG3CHK,NG4CHK,
     x                            DBE,DP,
     x                            XCHF_GAM2,XCHF_GAM2s,
     x                            XCHF_GAM3,XCHF_GAM4,
     x                            FBE,FP,SBE_XCHF,SP_XCHF,
     x                            E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                            E_XCHF_gam3,E_XCHF_gam4,
     x                            S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
         call RXCHF_fock_int_MPI(nproc,rank,
     x                           LCMF,LADDEXCH,nelec,NAE,NBE,
     x                           nebf,nebf2,npbf,npbf2,
     x                           ng1,ng2,ng3,ng4,
     x                           dimINT2,dimINT3,dimINT4,
     x                           dimINT2ex,dimINT3ex,
     x                           NG2CHK,NG3CHK,NG4CHK,
     x                           DAE,DBE,DP,
     x                           INT_GAM2,INT_GAM3,INT_GAM4,
     x                           INT_GAM2ex,INT_GAM3ex1,INT_GAM3ex2,
     x                           S_total,S_gam2,SBE_XCHF,SP_XCHF,
     x                           FPint,FAEint,FBEint, 
     x                           E_int_OMG2,E_int_OMG3,E_int_OMG4,
     x                           E_int)

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C ARS( no interaction
      if(LNOINT) then
          E_int=0.0d+00
          E_HF=0.0d+00
      else
          call add2fock(npbf,FPint,FP)
          call add2fock(nebf,FAEint,FAE_alp)
          call add2fock(nebf,FAEint,FAE_beta)
          call add2fock(nebf,FBEint,FBE)
      end if
C )

          IF ((LCMF).and.(rank.eq.0)) then
           write(*,*)
           write(*,*) "FAE_alp:"
           call prt_lower_triangle(nebf,nebflt,FAE_alp)
           write(*,*)
           write(*,*) "FAE_beta:"
           call prt_lower_triangle(nebf,nebflt,FAE_beta)
           write(*,*)
           write(*,*) "FBE:"
           call prt_lower_triangle(nebf,nebflt,FBE)
           write(*,*)
           write(*,*) "FP:"
           call prt_lower_triangle(npbf,npbflt,FP)
           write(*,*)
          END IF  

          E_total=E_HF+E_XCHF+E_int+E_nuc

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------)
         if((I.eq.1).and.(rank.eq.0)) then
            WRITE(*,9400)

      WRITE(*,9300) E_nuc,E_HF_ecore,E_HF_ee,E_HF,
     x  E_XCHF_gam1,E_XCHF_gam2,E_XCHF_gam3,E_XCHF_gam4,E_XCHF,
     x  E_int_OMG2,E_int_OMG3,E_int_OMG4,E_int,
     x  S_total,E_total

         end if

!        Fockp diag
         if(rank.eq.0) call UROOTHAN(vecP,EP,xxsp,FP,npbf)
         call MPI_BCAST(VECP,npbf2_,MPI_DOUBLE_PRECISION,
     x                  0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(EP,npbf_,MPI_DOUBLE_PRECISION,
     x                  0,MPI_COMM_WORLD,ierr)
         
         call construct_DP(nucst,npbf,vecP,DP)

C ARS( reform elec Fock matrices

C ARS( microiterate
!        --> FIND LARGEST CHANGE IN P DENSITY
         CALL DENDIF(DP0,DP,NPBF,DIFFP)
         CALL COPYDEN(DP0,DP,NPBF)
C )

         if(LSOSCFA) ITSOA=0
         if(LSOSCFB) ITSOB=0
         if(LSORXCHF) then
          ITSOA=0
          ITSOA_alp=0
          ITSOA_beta=0
          ITSOB=0
          if(rank.eq.0) then
           if(allocated(HESSA0)) HESSA0=0.0d+00
           if(allocated(HESSA0_alp)) HESSA0_alp=0.0d+00
           if(allocated(HESSA0_beta)) HESSA0_beta=0.0d+00
           if(allocated(HESSB0)) HESSB0=0.0d+00
           if(allocated(GRADA0)) GRADA0=0.0d+00
           if(allocated(GRADA0_alp)) GRADA0_alp=0.0d+00
           if(allocated(GRADA0_beta)) GRADA0_beta=0.0d+00
           if(allocated(GRADB0)) GRADB0=0.0d+00
          end if
         end if
         ORBGRDA=0.0d+00
         ORBGRDA_alp=0.0d+00
         ORBGRDA_beta=0.0d+00
         ORBGRDB=0.0d+00
         if(allocated(PGRADA)) PGRADA=0.0D+00
         if(allocated(PGRADA_alp)) PGRADA_alp=0.0D+00
         if(allocated(PGRADA_beta)) PGRADA_beta=0.0D+00
         if(allocated(PGRADB)) PGRADB=0.0D+00

         if (rank.eq.0) then
          write(*,2001) I
          if((LSOSCFA).and.(LSOSCFB)) then 
           WRITE(*,9055)
          else if ((LSOSCFA).and.(.not.(LSOSCFB))) then
           WRITE(*,9050)
          else if ((LSOSCFB).and.(.not.(LSOSCFA))) then
           WRITE(*,9051)
          else if ((LSORXCHF).and.(LDIIS)) then
           WRITE(*,9054)
          else if (LSORXCHF) then
           WRITE(*,9055)
          else if (LDIIS) then
           WRITE(*,9053)
          else
           WRITE(*,9000)
          end if
         end if

         if((LDIIS).and.(rank.eq.0)) then
          errvec=zero
          vecDIIS=zero
          itDIIS=0
          errDIIS=zero
         end if

         do ielec=1,maxmicroit

       if((.not.((locbse).or.(locbse2))).or.(ielec.eq.1)) then
C Call UHF Fock build for NAE regular electrons
c         call RXCHF_fock_hf_MPI(nproc,rank,
c     x                          LCMF,nebf,nebf2,NAE,ngee,
c     x                          DAE,GAM_ecore,GAM_ee,
c     x                          FAE,E_HF,E_HF_ecore,E_HF_ee)
         call RXCHF_fock_uhf_MPI(nproc,rank,
     x                          LCMF,nebf,nebf2,NAE,ngee,
     x                          DAE,GAM_ecore,GAM_ee,
     x							FAE,E_HF,E_HF_ecore,E_HF_ee,
     x                           DAE_alp,DAE_beta,FAE_alp,FAE_beta)

C Call XCHF Fock build for NBE special electrons and one QM particle
         call RXCHF_fock_xchf_MPI(nproc,rank,
     x                            LCMF,NBE,
     x                            nebf,nebf2,npbf,npbf2,
     x                            ngee,ng1,ng2,ng3,ng4,
     x                            dimXCHF2,dimXCHF3,dimXCHF4,
     x                            NG2CHK,NG3CHK,NG4CHK,
     x                            DBE,DP,
     x                            XCHF_GAM2,XCHF_GAM2s,
     x                            XCHF_GAM3,XCHF_GAM4,
     x                            FBE,FP,SBE_XCHF,SP_XCHF,
     x                            E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                            E_XCHF_gam3,E_XCHF_gam4,
     x                            S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
         call RXCHF_fock_int_MPI(nproc,rank,
     x                           LCMF,LADDEXCH,nelec,NAE,NBE,
     x                           nebf,nebf2,npbf,npbf2,
     x                           ng1,ng2,ng3,ng4,
     x                           dimINT2,dimINT3,dimINT4,
     x                           dimINT2ex,dimINT3ex,
     x                           NG2CHK,NG3CHK,NG4CHK,
     x                           DAE,DBE,DP,
     x                           INT_GAM2,INT_GAM3,INT_GAM4,
     x                           INT_GAM2ex,INT_GAM3ex1,INT_GAM3ex2,
     x                           S_total,S_gam2,SBE_XCHF,SP_XCHF,
     x                           FPint,FAEint,FBEint, 
     x                           E_int_OMG2,E_int_OMG3,E_int_OMG4,
     x                           E_int)

           call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C ARS( no interaction
      if(LNOINT) then
            E_int=0.0d+00
            E_HF=0.0d+00
      else
            call add2fock(npbf,FPint,FP)
            call add2fock(nebf,FAEint,FAE_alp)
            call add2fock(nebf,FAEint,FAE_beta)
            call add2fock(nebf,FBEint,FBE)
      end if
C )
C ARS( microiteration
          E_total=E_HF+E_XCHF+E_int+E_nuc
C )
      end if
C )

         if (LOCBSE) then
! Do OCBSE procedure (restricted solutions for regular and special electrons)

          if (rank.eq.0) then

! Regular alpha electrons

! Transform FAE_alp (calculated at end of previous iteration) to new WA_alp basis
!  - WA_alp constructed from vecAE0_alp and vecBE0 (previous iteration MOs)
!  - vecAE_alp in AO basis from previous iteration is transformed to new WA_alp basis
!    (relevant for SOSCF only)
           do j=1,nebf-nocca_alp-noccb
             tempvecA_alp(:,j)=vecBE0(:,j+noccb)
           end do
           do j=1,nae_alp
             tempvecA_alp(:,j+nebf-nocca_alp-noccb)=vecAE0_alp(:,j)
           end do
           if(LCMF) then
            write(*,*) "regular tempvecA_alp"
            call checkovlap(nwbfA_alp,nebf,tempvecA_alp,xxse)
           end if

           call RXCHFmult_OCBSE_transF(nebf,nwbfA_alp,tempvecA_alp,
     x                                 FAE_alp,WA_alp,wFAEw_alp)

!-----------------------POSSIBLE-SOSCF-ALPHA----------------------------(
         if((LSORXCHF).or.(LSOSCFA)) THEN
           ITER=IELEC
           EIGAVL = ITER.GT.1
         end if
         IF((LSORXCHF.or.LSOSCFA) .AND. (EIGAVL)) THEN
           call pack_LT(nwbfA_alp,nwbfLTA_alp,wFAEw_alp,wFLTAw_alp)
           call RXCHFmult_OCBSE_transVt(nebf,nwbfA_alp,WA_alp,
     x                                  xxse,vecAE_alp,wvecAEw_alp)
           if(LCMF) then
            write(*,*) "regular transVt vecAE_alp"
            call checkovlap(nebf,nebf,vecAE_alp,xxse)
           end if
           call SOGRAD(GRADA_alp,wFLTAw_alp,wvecAEw_alp,wWRKAw_alp,
     x                      NPRA_alp,NA_alp,
     x                 L0wA_alp,L1wA_alp,NwBFLTA_alp,ORBGRDA_alp)
C      write(*,*) "regular grada:",grada
            IF(ORBGRDA_alp.LT.SOGTOL  .OR.  ITSOA_alp.GT.0) THEN
              IF(ITSOA_alp.EQ.0) THEN
                 WRITE(*,9800)
                 call SOHESS(HSTARTA_alp,wAEenw_alp,NPRA_alp,
     x              L0wA_alp,NA_alp,NA_alp)
c                 write(*,*) "regular hstarta:",hstarta
                 if(LSORXCHF) then
                  do j=1,npra_alp
                    HESSA0_alp(j,j)=HSTARTA_alp(j)
                    DISPLIA0_alp(j)=-HSTARTA_alp(j)*GRADA_alp(j)
                  end do
                 end if
              END IF
              ITSOA_alp = ITSOA_alp+1
              IF(LSORXCHF) then
               call RXCHF_SOCHGBAS(NPRA_alp,NA_alp,nwbfA_alp,nebf,
     x                             WA0_alp,WA_alp,xxse,
     x                      GRADA0_alp,HESSA0_alp,DISPLIA0_alp)
               call RXCHF_SONEWT(NPRA_alp,ITSOA_alp,HESSA0_alp,
     x                         HESSA_alp,GRADA0_alp,
     x                        GRADA_alp,DISPLIA0_alp,DISPLIA_alp)
c      write(*,*) "grad:",GRADA
c      write(*,*) "grad0:",GRADA0
c      write(*,*) "displia:",DISPLIA
              ELSE
               call SONEWT(HSTARTA_alp,GRADA_alp,PGRADA_alp,
     x                DISPLIA_alp,DGRADA_alp,
     x        DISPLA_alp,UPDTA_alp,DISPLNA_alp,DGRADIA_alp,UPDTIA_alp,
     x                     ORBGRDA_alp,NPRA_alp,ITSOA_alp,NFT15)
C                 write(*,*) "regular displia:",displia
              END IF
              call SOTRAN(DISPLIA_alp,wvecAEw_alp,wGAw_alp,
     x                 wWRKAw_alp,NPRA_alp,
     x                    L0wA_alp,L1wA_alp,NA_alp,NA_alp,ORBGRDA_alp)
              IF(LSORXCHF) THEN
                CALL DCOPY(NPRA_alp,GRADA_alp,1,GRADA0_alp,1)
                CALL DCOPY(NPRA_alp,DISPLIA_alp,1,DISPLIA0_alp,1)
                CALL COPYDEN(HESSA0_alp,HESSA_alp,NPRA_alp)
              ELSE
                CALL DCOPY(NPRA_alp,GRADA_alp,1,PGRADA_alp,1)
              END IF
              call RXCHFmult_OCBSE_transV(nebf,nwbfA_alp,WA_alp,
     x               wvecAEw_alp,wAEenw_alp, ! eigenvalues useless
     x                                    vecAE_alp,AEe_alp)
              if(LCMF) then
               write(*,*) "regular transV vecAE_alp"
               call checkovlap(nebf,nebf,vecAE_alp,xxse)
              end if
              call RXCHFmult_construct_DE_uhf(NAE_alp,nebf,
     x              vecAE_alp,DAE_alp)
              GO TO 250
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA----------------------------)

! No SOSCF
!  - Diagonalize Fock matrix in WA basis of this iteration
!  - Obtain updated vecAE in WA basis of this iteration
           call RXCHFmult_OCBSE_diag(nebf,nwbfA_alp,WA_alp,wFAEw_alp,
     x                      wvecAEw_alp,wAEenw_alp,vecAE_alp,AEe_alp)
           if(LCMF) then
            write(*,*) "regular diag vecAE_alp"
            call checkovlap(nebf,nebf,vecAE_alp,xxse)
           end if
           call RXCHFmult_construct_DE_uhf(NAE_alp,nebf,
     x         vecAE_alp,DAE_alp)

  250 CONTINUE

! Regular beta electrons

! Transform FAE_beta (calculated at end of previous iteration) to new WA_beta basis
!  - WA_beta constructed from vecAE0_beta and vecBE0 (previous iteration MOs)
!  - vecAE_beta in AO basis from previous iteration is transformed to new WA_beta basis
!    (relevant for SOSCF only)
           do j=1,nebf-nocca_alp-noccb
             tempvecA_beta(:,j)=vecBE0(:,j+noccb)
           end do
           do j=1,nae_alp
             tempvecA_beta(:,j+nebf-nocca_alp-noccb)=vecAE0_alp(:,j)
           end do
           
           if(LCMF) then
            write(*,*) "regular tempvecA_beta"
            call checkovlap(nwbfA_beta,nebf,tempvecA_beta,xxse)
           end if

           call RXCHFmult_OCBSE_transF(nebf,nwbfA_beta,tempvecA_beta,
     x                                 FAE_beta,WA_beta,wFAEw_beta)

!-----------------------POSSIBLE-SOSCF-ALPHA----------------------------(
         if((LSORXCHF).or.(LSOSCFA)) THEN
           ITER=IELEC
           EIGAVL = ITER.GT.1
         end if
         IF((LSORXCHF.or.LSOSCFA) .AND. (EIGAVL)) THEN
           call pack_LT(nwbfA_beta,nwbfLTA_beta,wFAEw_beta,wFLTAw_beta)
           call RXCHFmult_OCBSE_transVt(nebf,nwbfA_beta,WA_beta,
     x                                  xxse,vecAE_beta,wvecAEw_beta)
           if(LCMF) then
            write(*,*) "regular transVt vecAE_beta"
            call checkovlap(nebf,nebf,vecAE_beta,xxse)
           end if
           call SOGRAD(GRADA_beta,wFLTAw_beta,wvecAEw_beta,wWRKAw_beta,
     x                      NPRA_beta,NA_beta,
     x                 L0wA_beta,L1wA_beta,NwBFLTA_beta,ORBGRDA_beta)
C      write(*,*) "regular grada:",grada
            IF(ORBGRDA_beta.LT.SOGTOL  .OR.  ITSOA_beta.GT.0) THEN
              IF(ITSOA_beta.EQ.0) THEN
c                 WRITE(*,*) '********* krb this is for beta *********'
                 WRITE(*,9800)
                 call SOHESS(HSTARTA_beta,wAEenw_beta,NPRA_beta,
     x              L0wA_beta,NA_beta,NA_beta)
C                 write(*,*) "regular hstarta:",hstarta
                 if(LSORXCHF) then
                  do j=1,npra_beta
                    HESSA0_beta(j,j)=HSTARTA_beta(j)
                    DISPLIA0_beta(j)=-HSTARTA_beta(j)*GRADA_beta(j)
                  end do
                 end if
              END IF
              ITSOA_beta = ITSOA_beta+1
              IF(LSORXCHF) then
               call RXCHF_SOCHGBAS(NPRA_beta,NA_beta,nwbfA_beta,nebf,
     x                             WA0_beta,WA_beta,xxse,
     x                      GRADA0_beta,HESSA0_beta,DISPLIA0_beta)
               call RXCHF_SONEWT(NPRA_beta,ITSOA_beta,HESSA0_beta,
     x                         HESSA_beta,GRADA0_beta,
     x                        GRADA_beta,DISPLIA0_beta,DISPLIA_beta)
C      write(*,*) "grad:",GRADA
C      write(*,*) "grad0:",GRADA0
C      write(*,*) "displia:",DISPLIA
              ELSE
               call SONEWT(HSTARTA_beta,GRADA_beta,PGRADA_beta,
     x                DISPLIA_beta,DGRADA_beta,
     x     DISPLA_beta,UPDTA_beta,DISPLNA_beta,DGRADIA_beta,UPDTIA_beta,
     x                     ORBGRDA_beta,NPRA_beta,ITSOA_beta,NFT17)
C                 write(*,*) "regular displia:",displia
              END IF
              call SOTRAN(DISPLIA_beta,wvecAEw_beta,wGAw_beta,
     x                 wWRKAw_beta,NPRA_beta,
     x                L0wA_beta,L1wA_beta,NA_beta,NA_beta,ORBGRDA_beta)
              IF(LSORXCHF) THEN
                CALL DCOPY(NPRA_beta,GRADA_beta,1,GRADA0_beta,1)
                CALL DCOPY(NPRA_beta,DISPLIA_beta,1,DISPLIA0_beta,1)
                CALL COPYDEN(HESSA0_beta,HESSA_beta,NPRA_beta)
              ELSE
                CALL DCOPY(NPRA_beta,GRADA_beta,1,PGRADA_beta,1)
              END IF
              call RXCHFmult_OCBSE_transV(nebf,nwbfA_beta,WA_beta,
     x               wvecAEw_beta,wAEenw_beta, ! eigenvalues useless
     x                                    vecAE_beta,AEe_beta)
              if(LCMF) then
               write(*,*) "regular transV vecAE_beta"
               call checkovlap(nebf,nebf,vecAE_beta,xxse)
              end if
              call RXCHFmult_construct_DE_uhf(NAE_beta,nebf,
     x              vecAE_beta,DAE_beta)
              GO TO 300
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA----------------------------)

! No SOSCF
!  - Diagonalize Fock matrix in WA basis of this iteration
!  - Obtain updated vecAE in WA basis of this iteration
           call RXCHFmult_OCBSE_diag(nebf,nwbfA_beta,WA_beta,wFAEw_beta,
     x                               wvecAEw_beta,wAEenw_beta,
     x                        vecAE_beta,AEe_beta)
           if(LCMF) then
            write(*,*) "regular diag vecAE_beta"
            call checkovlap(nebf,nebf,vecAE_beta,xxse)
           end if
           call RXCHFmult_construct_DE_uhf(NAE_beta,nebf,
     x                      vecAE_beta,DAE_beta)

  300 CONTINUE

c add together alpha and beta density matrices
       dae = 0.0
       do k = 1,nebf
         do l = 1, nebf
            DAE(k,l) = DAE_alp(k,l)+DAE_beta(k,l)
         enddo
       enddo
     
! Special electrons

! Transform FBE (calculated at end of previous iteration) to new WB basis
!  - WB constructed from vecAE and vecBE0 (current and previous iteration MOs)
!  - vecBE in AO basis from previous iteration is transformed to new WB basis
!    (relevant for SOSCF only)

           do j=1,nebf-nocca_alp-noccb
             tempvecB(:,j)=vecAE_alp(:,j+nocca_alp)
           end do
           do j=1,noccb
             tempvecB(:,j+nebf-nocca_alp-noccb)=vecBE0(:,j)
           end do
           if(LCMF) then
            write(*,*) "special tempvecB"
            call checkovlap(nwbfB,nebf,tempvecB,xxse)
           end if

           call RXCHFmult_OCBSE_transF(nebf,nwbfB,tempvecB,
     x                                 FBE,WB,wFBEw)
!-----------------------POSSIBLE-SOSCF-BETA----------------------------(
         if((LSORXCHF).or.(LSOSCFB)) THEN
           ITER=IELEC
           EIGAVL = ITER.GT.1
         end if
         IF((LSORXCHF.or.LSOSCFB) .AND. (EIGAVL)) THEN
           call pack_LT(nwbfB,nwbfLTB,wFBEw,wFLTBw)
           call RXCHFmult_OCBSE_transVt(nebf,nwbfB,WB,
     x                                  xxse,vecBE,wvecBEw)
           if(LCMF) then
            write(*,*) "special transVt vecBE"
            call checkovlap(nebf,nebf,vecBE,xxse)
           end if
           call SOGRAD(GRADB,wFLTBw,wvecBEw,wWRKBw,NPRB,NB,
     x                 L0wB,L1wB,NwBFLTB,ORBGRDB)
C      write(*,*) "special gradb:",gradb
            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
              IF(ITSOB.EQ.0) THEN
                 WRITE(*,9800)
                 call SOHESS(HSTARTB,wBEenw,NPRB,L0wB,NB,NB)
C                 write(*,*) "special hstartb:",hstartb
                 if(LSORXCHF) then
                  do j=1,nprb
                    HESSB0(j,j)=HSTARTB(j)
                    DISPLIB0(j)=-HSTARTB(j)*GRADB(j)
                  end do
                 end if
              END IF
              ITSOB = ITSOB+1
              IF(LSORXCHF) then
               call RXCHF_SOCHGBAS(NPRB,NB,nwbfB,nebf,WB0,WB,xxse,
     x                             GRADB0,HESSB0,DISPLIB0)
               call RXCHF_SONEWT(NPRB,ITSOB,HESSB0,HESSB,GRADB0,
     x                           GRADB,DISPLIB0,DISPLIB)
C      write(*,*) "grad:",GRADB
C      write(*,*) "grad0:",GRADB0
C      write(*,*) "displia:",DISPLIB
              ELSE
               call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,
     x                     DISPLB,UPDTB,DISPLNB,DGRADIB,UPDTIB,
     x                     ORBGRDB,NPRB,ITSOB,NFT16)
C                 write(*,*) "special displib:",displib
              END IF
              call SOTRAN(DISPLIB,wvecBEw,wGBw,wWRKBw,NPRB,
     x                    L0wB,L1wB,NB,NB,ORBGRDB)
              IF(LSORXCHF) THEN
                CALL DCOPY(NPRB,GRADB,1,GRADB0,1)
                CALL DCOPY(NPRB,DISPLIB,1,DISPLIB0,1)
                CALL COPYDEN(HESSB0,HESSB,NPRB)
              ELSE
                CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
              END IF
              call RXCHFmult_OCBSE_transV(nebf,nwbfB,WB,wvecBEw,wBEenw, ! eigenvalues useless
     x                                    vecBE,BEe)
c            write(6,*) 'krb  vecbe after'
c                  call PREVNU(vecbe,AEE,nebf,nebf,nebf)
                  if(LCMF) then
               write(*,*) "special transV vecBE"
               call checkovlap(nebf,nebf,vecBE,xxse)
              end if
              call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)
              GO TO 350
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-BETA----------------------------)

! No SOSCF
!  - Diagonalize Fock matrix in WB basis of this iteration
!  - Obtain updated vecBE in WB basis of this iteration
           call RXCHFmult_OCBSE_diag(nebf,nwbfB,WB,wFBEw,
     x                               wvecBEw,wBEenw,vecBE,BEe)
c            write(6,*) 'krb  vecbe after2'
c                  call PREVNU(vecbe,AEE,nebf,nebf,nebf)
            if(LCMF) then
            write(*,*) "special diag vecBE"
            call checkovlap(nebf,nebf,vecBE,xxse)
           end if
           call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)

  350 CONTINUE


C --------------------------- POTENTIAL DIIS ---------------------------
      if(LDIIS) then

       CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
       CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)

       if(itDIIS.eq.0) then ! If DIIS not intiated, check if it should be
        if((DIFFAE.lt.threshDIIS).and.(DIFFBE.lt.threshDIIS)) then
         write(*,9801)
         itDIIS=itDIIS+1
         do j=1,nocca
           vecDIIS(:,j,nstore)=vecAE0(:,j)
         end do
         do j=1,noccb
           vecDIIS(:,j+nocca,nstore)=vecBE0(:,j)
         end do
         do j=1,nebf-nocca-noccb
           vecDIIS(:,j+nocca+noccb,nstore)=vecBE0(:,j+noccb)
         end do
        end if
       end if

       if(itDIIS.gt.0) then ! If DIIS already initiated, continue with DIIS
        call DIIS_driver(nebf,nocca,noccb,nstore,itDIIS,
     x                   vecAE,vecBE,xxse,
     x                   errvec,vecDIIS,errDIIS)
        itDIIS=itDIIS+1
        if(LCMF) then
         write(*,*) "regular DIIS vecAE"
         call checkovlap(nebf,nebf,vecAE,xxse)
         write(*,*) "special DIIS vecBE"
         call checkovlap(nebf,nebf,vecBE,xxse)
        end if
        call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)
        call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)
       end if

      end if
C ----------------------------------------------------------------------

           WA0_alp(:,:)=WA_alp(:,:)
           WA0_beta(:,:)=WA_beta(:,:)
           WB0(:,:)=WB(:,:)

          end if ! rank 0


          call MPI_BCAST(DAE,nebf2_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
     
          call MPI_BCAST(VECAE_alp,nebf2_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(AEe_alp,nebf_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(DAE_alp,nebf2_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(VECAE_beta,nebf2_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(AEe_beta,nebf_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(DAE_beta,nebf2_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(VECBE,nebf2_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(BEe,nebf_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(DBE,nebf2_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)

          CALL DENDIF(DAE0_alp,DAE_alp,NEBF,DIFFAE_alp)
          CALL DENDIF(DAE0_beta,DAE_beta,NEBF,DIFFAE_beta)
          CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)

          CALL COPYDEN(DAE0_alp,DAE_alp,NEBF)
          CALL COPYDEN(DAE0_beta,DAE_beta,NEBF)
          CALL COPYDEN(DBE0,DBE,NEBF)

          CALL COPYDEN(vecAE0_alp,vecAE_alp,NEBF)
          CALL COPYDEN(vecAE0_beta,vecAE_beta,NEBF)
          CALL COPYDEN(vecBE0,vecBE,NEBF)

! Calculate energy for this it and Fock matrices for next it

C Call UHF Fock build for NAE regular electrons
c         call RXCHF_fock_hf_MPI(nproc,rank,
c     x                          LCMF,nebf,nebf2,NAE,ngee,
c     x                          DAE,GAM_ecore,GAM_ee,
c     x                          FAE,E_HF,E_HF_ecore,E_HF_ee)
         call RXCHF_fock_uhf_MPI(nproc,rank,
     x                          LCMF,nebf,nebf2,NAE,ngee,
     x                          DAE,GAM_ecore,GAM_ee,
     x							FAE,E_HF,E_HF_ecore,E_HF_ee,
     x                           DAE_alp,DAE_beta,FAE_alp,FAE_beta)

C Call XCHF Fock build for NBE special electrons and one QM particle
         call RXCHF_fock_xchf_MPI(nproc,rank,
     x                            LCMF,NBE,
     x                            nebf,nebf2,npbf,npbf2,
     x                            ngee,ng1,ng2,ng3,ng4,
     x                            dimXCHF2,dimXCHF3,dimXCHF4,
     x                            NG2CHK,NG3CHK,NG4CHK,
     x                            DBE,DP,
     x                            XCHF_GAM2,XCHF_GAM2s,
     x                            XCHF_GAM3,XCHF_GAM4,
     x                            FBE,FP,SBE_XCHF,SP_XCHF,
     x                            E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                            E_XCHF_gam3,E_XCHF_gam4,
     x                            S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
         call RXCHF_fock_int_MPI(nproc,rank,
     x                           LCMF,LADDEXCH,nelec,NAE,NBE,
     x                           nebf,nebf2,npbf,npbf2,
     x                           ng1,ng2,ng3,ng4,
     x                           dimINT2,dimINT3,dimINT4,
     x                           dimINT2ex,dimINT3ex,
     x                           NG2CHK,NG3CHK,NG4CHK,
     x                           DAE,DBE,DP,
     x                           INT_GAM2,INT_GAM3,INT_GAM4,
     x                           INT_GAM2ex,INT_GAM3ex1,INT_GAM3ex2,
     x                           S_total,S_gam2,SBE_XCHF,SP_XCHF,
     x                           FPint,FAEint,FBEint, 
     x                           E_int_OMG2,E_int_OMG3,E_int_OMG4,
     x                           E_int)

           call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C ARS( no interaction
      if(LNOINT) then
            E_int=0.0d+00
            E_HF=0.0d+00
      else
            call add2fock(npbf,FPint,FP)
            call add2fock(nebf,FAEint,FAE_alp)
            call add2fock(nebf,FAEint,FAE_beta)
            call add2fock(nebf,FBEint,FBE)
      end if
C )

            IF ((LCMF).and.(rank.eq.0)) then
             write(*,*)
             write(*,*) "FAE_alp:"
             call prt_lower_triangle(nebf,nebflt,FAE_alp)
             write(*,*)
             write(*,*) "FAE_beta:"
             call prt_lower_triangle(nebf,nebflt,FAE_beta)
             write(*,*)
             write(*,*) "FBE:"
             call prt_lower_triangle(nebf,nebflt,FBE)
             write(*,*)
             write(*,*) "FP:"
             call prt_lower_triangle(npbf,npbflt,FP)
             write(*,*)
            END IF  

            E_total=E_HF+E_XCHF+E_int+E_nuc

         else if (LOCBSE2) then
! Do OCBSE2 procedure (restricted solutions for special electrons)

! Regular electrons
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------(
      if (rank.eq.0) then
         if((LSORXCHF).or.(LSOSCFA)) THEN
          ITER=IELEC
          EIGAVL = ITER.GT.1
         end if
         IF((LSORXCHF.or.LSOSCFA) .AND.  EIGAVL) THEN
           call pack_LT(nebf,nebfLT,FAE,FLT)
           call SOGRAD(GRADA,FLT,vecAE,WRK,NPRA,NA,L0,L1,NEBFLT,ORBGRDA)
!      write(*,*) "regular grada:",grada
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AEE,NPRA,L0,NA,NA)
!                 write(*,*) "regular hstarta:",hstarta
                 if(LSORXCHF) then
                  do j=1,npra
                    HESSA0(j,j)=HSTARTA(j)
                    DISPLIA0(j)=-HSTARTA(j)*GRADA(j)
                  end do
                 end if
              END IF
              ITSOA = ITSOA+1
              IF(LSORXCHF) then
               call RXCHF_SONEWT(NPRA,ITSOA,HESSA0,HESSA,GRADA0,
     x                           GRADA,DISPLIA0,DISPLIA)
C      write(*,*) "grad:",GRADA
C      write(*,*) "grad0:",GRADA0
C      write(*,*) "displia:",DISPLIA
              ELSE
C      write(*,*) "hess:",HSTARTA
C      write(*,*) "grad:",GRADA
C      write(*,*) "pgrad:",PGRADA
               call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,
     x                     DISPLA,UPDTA,DISPLNA,DGRADIA,UPDTIA,
     x                     ORBGRDA,NPRA,ITSOA,NFT15)
!                 write(*,*) "regular displia:",displia
              END IF
              call SOTRAN(DISPLIA,vecAE,GA,WRK,NPRA,L0,L1,NA,NA,ORBGRDA)
              IF(LSORXCHF) THEN
                CALL DCOPY(NPRA,GRADA,1,GRADA0,1)
                CALL DCOPY(NPRA,DISPLIA,1,DISPLIA0,1)
                CALL COPYDEN(HESSA0,HESSA,NPRA)
              ELSE
                CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              END IF
              call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)
              if(LCMF) then
               write(*,*) "regular SORXCHF vecAE"
               call checkovlap(nebf,nebf,vecAE,xxse)
              end if
              GO TO 950
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------)

!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAE,AEE,xxse,FAE,nebf)
         call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)
         if(LCMF) then
          write(*,*) "regular diag vecAE",ielec
          call checkovlap(nebf,nebf,vecAE,xxse)
         end if

  950 CONTINUE

C ARS( OCBSE/SOSCF

! Transform FBE (calculated at end of previous iteration) to new W basis
!  - W updated with new vecA from this iteration
!  - vecBE in AO basis from previous iteration is transformed to new W basis
!    (relevant for SOSCF only)
         call RXCHFmult_OCBSE_transF(nebf,nwbfB,vecAE(:,nocca+1:nebf),
     x                               FBE,WB,wFBEw)

!-----------------------POSSIBLE-SOSCF-BETA----------------------------(
         if((LSORXCHF).or.(LSOSCFB)) THEN
           ITER=IELEC
           EIGAVL = ITER.GT.1
         end if
         IF((LSORXCHF.or.LSOSCFB) .AND. (EIGAVL)) THEN
           call pack_LT(nwbfB,nwbfLTB,wFBEw,wFLTBw)
           call RXCHFmult_OCBSE_transVt(nebf,nwbfB,WB,
     x                                  xxse,vecBE,wvecBEw)
           if(LCMF) then
            write(*,*) "special transVt vecBE"
            call checkovlap(nebf,nebf,vecBE,xxse)
           end if
           call SOGRAD(GRADB,wFLTBw,wvecBEw,wWRKBw,NPRB,NB,
     x                 L0wB,L1wB,NwBFLTB,ORBGRDB)
!      write(*,*) "special gradb:",gradb
            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
              IF(ITSOB.EQ.0) THEN
                 WRITE(*,9800)
                 call SOHESS(HSTARTB,wBEenw,NPRB,L0wB,NB,NB)
!                 write(*,*) "special hstartb:",hstartb
                 if(LSORXCHF) then
                  do j=1,nprb
                    HESSB0(j,j)=HSTARTB(j)
                    DISPLIB0(j)=-HSTARTB(j)*GRADB(j)
                  end do
                 end if
              END IF
              ITSOB = ITSOB+1
              IF(LSORXCHF) then
               call RXCHF_SOCHGBAS(NPRB,NB,nwbfB,nebf,WB0,WB,xxse,
     x                             GRADB0,HESSB0,DISPLIB0)
               call RXCHF_SONEWT(NPRB,ITSOB,HESSB0,HESSB,GRADB0,
     x                           GRADB,DISPLIB0,DISPLIB)
C      write(*,*) "grad:",GRADB
C      write(*,*) "grad0:",GRADB0
C      write(*,*) "displia:",DISPLIB
              ELSE
               call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,
     x                     DISPLB,UPDTB,DISPLNB,DGRADIB,UPDTIB,
     x                     ORBGRDB,NPRB,ITSOB,NFT16)
C                 write(*,*) "special displib:",displib
              END IF
              call SOTRAN(DISPLIB,wvecBEw,wGBw,wWRKBw,NPRB,
     x                    L0wB,L1wB,NB,NB,ORBGRDB)
              IF(LSORXCHF) THEN
                CALL DCOPY(NPRB,GRADB,1,GRADB0,1)
                CALL DCOPY(NPRB,DISPLIB,1,DISPLIB0,1)
                CALL COPYDEN(HESSB0,HESSB,NPRB)
              ELSE
                CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
              END IF
              call RXCHFmult_OCBSE_transV(nebf,nwbfB,WB,wvecBEw,wBEenw, ! eigenvalues useless
     x                                    vecBE,BEe)
              if(LCMF) then
               write(*,*) "special transV vecBE"
               call checkovlap(nebf,nebf,vecBE,xxse)
              end if
              call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)
              GO TO 450
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-BETA----------------------------)

! No SOSCF
!  - Diagonalize Fock matrix in W basis of this iteration
!  - Obtain updated vecBE in W basis of this iteration
         call RXCHFmult_OCBSE_diag(nebf,nwbfB,WB,wFBEw,
     x                             wvecBEw,wBEenw,vecBE,BEe)
         if(LCMF) then
          write(*,*) "special diag vecBE"
          call checkovlap(nebf,nebf,vecBE,xxse)
         end if
         call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)

  450 CONTINUE

       WB0(:,:)=WB(:,:)

      end if ! rank 0

           call MPI_BCAST(VECAE,nebf2_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(AEe,nebf_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(DAE,nebf2_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(VECBE,nebf2_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(BEe,nebf_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(DBE,nebf2_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)

         CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
         CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)

         CALL COPYDEN(DAE0,DAE,NEBF)
         CALL COPYDEN(DBE0,DBE,NEBF)

C )

! Calculate energy for this it and Fock matrices for next it

C Call HF Fock build for NAE regular electrons
           call RXCHF_fock_hf_MPI(nproc,rank,
     x                            LCMF,nebf,nebf2,NAE,ngee,
     x                            DAE,GAM_ecore,GAM_ee,
     x                            FAE,E_HF,E_HF_ecore,E_HF_ee)

C Call XCHF Fock build for NBE special electrons and one QM particle
         call RXCHF_fock_xchf_MPI(nproc,rank,
     x                            LCMF,NBE,
     x                            nebf,nebf2,npbf,npbf2,
     x                            ngee,ng1,ng2,ng3,ng4,
     x                            dimXCHF2,dimXCHF3,dimXCHF4,
     x                            NG2CHK,NG3CHK,NG4CHK,
     x                            DBE,DP,
     x                            XCHF_GAM2,XCHF_GAM2s,
     x                            XCHF_GAM3,XCHF_GAM4,
     x                            FBE,FP,SBE_XCHF,SP_XCHF,
     x                            E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                            E_XCHF_gam3,E_XCHF_gam4,
     x                            S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
         call RXCHF_fock_int_MPI(nproc,rank,
     x                           LCMF,LADDEXCH,nelec,NAE,NBE,
     x                           nebf,nebf2,npbf,npbf2,
     x                           ng1,ng2,ng3,ng4,
     x                           dimINT2,dimINT3,dimINT4,
     x                           dimINT2ex,dimINT3ex,
     x                           NG2CHK,NG3CHK,NG4CHK,
     x                           DAE,DBE,DP,
     x                           INT_GAM2,INT_GAM3,INT_GAM4,
     x                           INT_GAM2ex,INT_GAM3ex1,INT_GAM3ex2,
     x                           S_total,S_gam2,SBE_XCHF,SP_XCHF,
     x                           FPint,FAEint,FBEint, 
     x                           E_int_OMG2,E_int_OMG3,E_int_OMG4,
     x                           E_int)

           call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C ARS( no interaction
      if(LNOINT) then
            E_int=0.0d+00
            E_HF=0.0d+00
      else
            call add2fock(npbf,FPint,FP)
            call add2fock(nebf,FAEint,FAE)
            call add2fock(nebf,FBEint,FBE)
      end if
C )

            IF ((LCMF).and.(rank.eq.0)) then
             write(*,*)
             write(*,*) "FAE:"
             call prt_lower_triangle(nebf,nebflt,FAE)
             write(*,*)
             write(*,*) "FBE:"
             call prt_lower_triangle(nebf,nebflt,FBE)
             write(*,*)
             write(*,*) "FP:"
             call prt_lower_triangle(npbf,npbflt,FP)
             write(*,*)
            END IF  

            E_total=E_HF+E_XCHF+E_int+E_nuc

         else

!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------(
      if (rank.eq.0) then
         if(LSOSCFA) THEN
          ITER=IELEC
          EIGAVL = ITER.GT.1
         end if
         IF(LSOSCFA .AND.  EIGAVL) THEN                ! first it. skip SOSCF (diag to get EE)
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAE,FLT)
          call SOGRAD(GRADA,FLT,vecAE,WRK,NPRA,NA,L0,L1,NEBFLT,ORBGRDA)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN   ! only on first SOSCF it. set up approx Hess
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AEE,NPRA,L0,NA,NA)
              END IF
              ITSOA = ITSOA+1
           call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,DISPLA,UPDTA,
     *                 DISPLNA,DGRADIA,UPDTIA,ORBGRDA,NPRA,ITSOA,NFT15)
            call SOTRAN(DISPLIA,vecAE,GA,WRK,NPRA,L0,L1,NA,NA,ORBGRDA)
             CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)
              GO TO 750  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------)

  700 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAE,AEE,xxse,FAE,nebf)
         call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)

      end if ! rank 0

  750 CONTINUE

           call MPI_BCAST(VECAE,nebf2_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(AEe,nebf_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(DAE,nebf2_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)

!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
         CALL COPYDEN(DAE0,DAE,NEBF)

!-----------------------POSSIBLE-SOSCF-BETA----------------------------(
      if (rank.eq.0) then
        if(LSOSCFB) THEN
         ITER=IELEC
         EIGAVL = ITER.GT.1
        end if
         IF(LSOSCFB .AND.  EIGAVL) THEN                ! first it. skip SOSCF (diag to get EE)
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FBE,FLT)
          call SOGRAD(GRADB,FLT,vecBE,WRK,NPRB,NB,L0,L1,NEBFLT,ORBGRDB)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
              IF(ITSOB.EQ.0) THEN   ! only on first SOSCF it. set up approx Hess
             WRITE(*,9800)
                 call SOHESS(HSTARTB,BEE,NPRB,L0,NB,NB)
              END IF
              ITSOB = ITSOB+1
           call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,DISPLB,UPDTB,
     *                 DISPLNB,DGRADIB,UPDTIB,ORBGRDB,NPRB,ITSOB,NFT16)
            call SOTRAN(DISPLIB,vecBE,GB,WRK,NPRB,L0,L1,NB,NB,ORBGRDB)
             CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
              call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)
              GO TO 850  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-BETA----------------------------)

  800 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DBE,vecBE,BEE,xxse,FBE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecBE,BEE,xxse,FBE,nebf)
         call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)

      end if ! rank 0

  850 CONTINUE

           call MPI_BCAST(VECBE,nebf2_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(BEe,nebf_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)
           call MPI_BCAST(DBE,nebf2_,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)

!        --> FIND LARGEST CHANGE IN Beta E DENSITY
         CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
         CALL COPYDEN(DBE0,DBE,NEBF)

         end if ! end if for not ocbse or ocbse2

C ARS( microiterate
C!        --> FIND LARGEST CHANGE IN P DENSITY
C         CALL DENDIF(DP0,DP,NPBF,DIFFP)
C         CALL COPYDEN(DP0,DP,NPBF)
C )

!        --> CALCULATE CHANGE IN TOTAL ENERGY
         Delta_E_tot=E_total-E_total_old
         E_total_old=E_total

!        --> PRINT SUMMARY OF THIS ITERATION
         if (rank.eq.0) then
          if((LSOSCFA).and.(LSOSCFB)) then
             WRITE(*,9153) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP,ORBGRDA,ORBGRDB
          else if ((LSOSCFA).and.(.not.(LSOSCFB))) then
             WRITE(*,9150) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP,ORBGRDA
          else if ((LSOSCFB).and.(.not.(LSOSCFA))) then
             WRITE(*,9150) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP,ORBGRDB
          else if ((LDIIS).and.(LSORXCHF)) then
             WRITE(*,9152) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP,ORBGRDA,ORBGRDB,errDIIS
          else if (LSORXCHF) then
             WRITE(*,9153) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE_alp,DIFFAE_beta,DIFFBE,DIFFP,
     x                      ORBGRDA_alp,ORBGRDA_beta,ORBGRDB
          else if (LDIIS) then
             WRITE(*,9150) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP,errDIIS
          else
             WRITE(*,9100) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP
          end if
C ARS( debug: print out MOs here
          if (LCMF) then
           WRITE(*,9610)
           WRITE(*,*) 'ALPHA ORBITALS'
           call PREVNU(vecAE_alp,AEE_alp,nebf,nebf,nebf)
           WRITE(*,9620)
           call PREVNU(vecBE,BEE,nebf,nebf,nebf)
           WRITE(*,9700)
           call PREVNU(vecp,EP,npbf,npbf,npbf)
          end if
C )
C ARS( check orthonormality of MOs
          if(LCMF) then
           write(*,*) "Regular ALPHA MO overlap matrix:"
           call checkovlap(nebf,nebf,vecAE_alp,xxse)
           write(*,*) "Regular BETA MO overlap matrix:"
           call checkovlap(nebf,nebf,vecAE_beta,xxse)
           write(*,*) "Special MO overlap matrix:"
           call checkovlap(nebf,nebf,vecBE,xxse)
          end if
C )
! Output the vectors for this iteration for restart if necessary:
          call write_MOs(860,nebf,VECAE_alp)
          call write_MOs(861,nebf,VECAE_beta)
          call write_MOs(862,nebf,VECBE)
          call write_MOs(853,npbf,VECP)
         end if ! rank 0

         LDIFFE=( (DIFFAE.LT.TOLE).and.(DIFFBE.LT.TOLE) )
         IF(LDIFFE) GOTO 200
         IF(IELEC.EQ.MAXMICROIT) GOTO 10

       END DO  ! microiterations

  200 CONTINUE
!      IF WE GET HERE - MICROITERATION CONVERGENCE ACHIEVED
       if(rank.eq.0) write(*,2000) i,ielec

       IF(DIFFP.LT.TOLP) GOTO 100
       IF(I.EQ.MAXIT) GOTO 10

      END DO   ! iterations
 
  10  CONTINUE
!     IF WE GET HERE SOMETHING WENT WRONG

      if (rank.eq.0) then
       if(LSOSCFA) close(NFT15)
       if(LSOSCFB) close(NFT16)

       WRITE(*,*)
       WRITE(*,*)'WARNING:  ITERATION LIMIT EXCEEDED'
       WRITE(*,*)
      end if

      E_total=zero

  100 CONTINUE
!     IF WE GET HERE WE ARE DONE - CONVERGENCE ACHIEVED

      if (rank.eq.0) then

       if(LSOSCFA) close(NFT15)
       if(LSOSCFB) close(NFT16)

!     PRINT FINAL ENERGY AND PUNCH THE ORBITALS
       WRITE(*,9200) E_total,I

       WRITE(*,9300) E_nuc,E_HF_ecore,E_HF_ee,E_HF,
     x   E_XCHF_gam1,E_XCHF_gam2,E_XCHF_gam3,E_XCHF_gam4,E_XCHF,
     x   E_int_OMG2,E_int_OMG3,E_int_OMG4,E_int,
     x   S_total,E_total

!  OUTPUT ELEC AND NUC EIGENVALUES AND EIGENVECTORS
       WRITE(*,9610)
       WRITE(*,*) '******ALPHA EIGENVALUES *****'
       call PREVNU(vecAE_alp,AEE_alp,nebf,nebf,nebf)
       WRITE(*,9610)
       WRITE(*,*) '******BETA EIGENVALUES *****'
       call PREVNU(vecAE_beta,AEE_beta,nebf,nebf,nebf)
       WRITE(*,9620)
       call PREVNU(vecBE,BEE,nebf,nebf,nebf)
       WRITE(*,9700)
       call PREVNU(vecp,EP,npbf,npbf,npbf)

! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------(
!     subroutine write_MOs(IFIL,nbf,VEC)
!     IFIL=852 :: FinalCE.dat
!     IFIL=853 :: FinalCP.dat
!     IFIL=860 :: FinalCAE.dat
!     IFIL=861 :: FinalCBE.dat
       call write_MOs(870,nebf,VECAE_alp)
       call write_MOs(871,nebf,VECAE_beta)       
       call write_MOs(861,nebf,VECBE)
       call write_MOs(853,npbf,VECP)
! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------)

      end if

      if((LOCBSE.or.LOCBSE2).and.(rank.eq.0)) then
       if(allocated(wWRKAw))   deallocate(wWRKAw)
       if(allocated(wGAw))     deallocate(wGAw)
       if(allocated(wFLTAw))   deallocate(wFLTAw)
       if(allocated(wAEenw))   deallocate(wAEenw)
       if(allocated(wvecAEw))  deallocate(wvecAEw)
       if(allocated(wFAEw))    deallocate(wFAEw)
       if(allocated(WA))       deallocate(WA)
       if(allocated(tempvecA)) deallocate(tempvecA)
       if(allocated(wWRKBw))   deallocate(wWRKBw)
       if(allocated(wGBw))     deallocate(wGBw)
       if(allocated(wFLTBw))   deallocate(wFLTBw)
       if(allocated(wBEenw))   deallocate(wBEenw)
       if(allocated(wvecBEw))  deallocate(wvecBEw)
       if(allocated(wFBEw))    deallocate(wFBEw)
       if(allocated(WB))       deallocate(WB)
       if(allocated(tempvecB)) deallocate(tempvecB)
       
       if(allocated(wWRKAw_alp))   deallocate(wWRKAw_alp)
       if(allocated(wGAw_alp))     deallocate(wGAw_alp)
       if(allocated(wFLTAw_alp))   deallocate(wFLTAw_alp)
       if(allocated(wAEenw_alp))   deallocate(wAEenw_alp)
       if(allocated(wvecAEw_alp))  deallocate(wvecAEw_alp)
       if(allocated(wFAEw_alp))    deallocate(wFAEw_alp)
       if(allocated(WA_alp))       deallocate(WA_alp)
       if(allocated(tempvecA_alp)) deallocate(tempvecA_alp)
       
       if(allocated(wWRKAw_beta))   deallocate(wWRKAw_beta)
       if(allocated(wGAw_beta))     deallocate(wGAw_beta)
       if(allocated(wFLTAw_beta))   deallocate(wFLTAw_beta)
       if(allocated(wAEenw_beta))   deallocate(wAEenw_beta)
       if(allocated(wvecAEw_beta))  deallocate(wvecAEw_beta)
       if(allocated(wFAEw_beta))    deallocate(wFAEw_beta)
       if(allocated(WA_beta))       deallocate(WA_beta)
       if(allocated(tempvecA_beta)) deallocate(tempvecA_beta)
       
      end if
      if((LDIIS).and.(rank.eq.0)) then
       if(allocated(errvec)) deallocate(errvec)
      end if
      if(LSORXCHF) then
       if(allocated(DISPLIB0)) deallocate(DISPLIB0)
       if(allocated(DISPLIA0)) deallocate(DISPLIA0)
!       if(allocated(DISPLIB))  deallocate(DISPLIB)
!       if(allocated(DISPLIA))  deallocate(DISPLIA)
       if(allocated(GRADB0))   deallocate(GRADB0)
       if(allocated(GRADA0))   deallocate(GRADA0)
!       if(allocated(GRADB))    deallocate(GRADB)
!       if(allocated(GRADA))    deallocate(GRADA)
       if(allocated(HESSB0))   deallocate(HESSB0)
       if(allocated(HESSA0))   deallocate(HESSA0)
       if(allocated(HESSB))    deallocate(HESSB)
       if(allocated(HESSA))    deallocate(HESSA)
       if(allocated(WB0))      deallocate(WB0)
       if(allocated(WA0))      deallocate(WA0)
      end if

      return
      end

