!======================================================================
      subroutine RXCHF_scf_MPI(nproc,rank,
     x                         nelec,NAE,NBE,NEBFLT,NPBFLT,NUCST,
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
!
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
      integer NPRA,NPRB
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
      double precision DBE(NEBF,NEBF)
      double precision DBE0(NEBF,NEBF)
      double precision VECAE(NEBF,NEBF)
      double precision VECAE0(NEBF,NEBF)
      double precision VECAEi(NEBF,NEBF)
      double precision VECBE(NEBF,NEBF)
      double precision VECBE0(NEBF,NEBF)
      double precision VECBEi(NEBF,NEBF)
      double precision AEE(NEBF)
      double precision BEE(NEBF)
      double precision FAE(nebf,nebf)
      double precision XFAE(nebf,nebf)
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
      integer NB
      integer ITER
      integer ITSOA ! SOSCF iteration counter
      integer ITSOB ! SOSCF iteration counter
      integer L0,L1
      integer L0wA,L1wA
      integer L0wB,L1wB
      integer NFT15
      integer NFT16
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
C )
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

 9100 FORMAT(1X,I3,F20.10,F17.10,3F17.10)

 9150 FORMAT(1X,I3,F20.10,F17.10,4F17.10)

 9151 FORMAT(1X,I3,F20.10,F17.10,5F17.10)

 9152 FORMAT(1X,I3,F20.10,F17.10,6F17.10)

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

      if (nae.gt.1) then
       nocca=nae/2
      else
       nocca=nae
      end if
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
        end if
        if(LSORXCHF) then
         SOGTOL=0.40d+00
         npra=(nebf-nocca-noccb)*nocca ! occ-vir pairs for regular elecs
         nprb=(nebf-noccb-nocca)*noccb ! occ-vir pairs for special elecs
         NA=nocca
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
          call RXCHFmult_guess_elec(nae,nbe,nebf,xxse,GAM_ecore,
     x                          DAE,DBE,VECAE0,VECBE0)
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
       call PREVNU(vecAE0,AEE,nebf,nebf,nebf)
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
      if (LSOSCF) then
         SOGTOL=0.40d+00
         SMALL=1.0D-06
         L0=nebf
         L1=nebf
         LSOSCFA=.true.
         LSOSCFB=.true.
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
      DBE0=0.0d+00
      DP0=0.0d+00

!     BEGIN XCSCF ITERATIONS
      if(rank.eq.0) WRITE(*,9500)
      E_total_old=0.0d+00

      DO I=1,MAXIT

C Call HF Fock build for NAE regular electrons
         call RXCHF_fock_hf_MPI(nproc,rank,
     x                          LCMF,nebf,nebf2,NAE,ngee,
     x                          DAE,GAM_ecore,GAM_ee,
     x                          FAE,E_HF,E_HF_ecore,E_HF_ee)

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

         call COPYDEN(vecAEi,vecAE0,nebf)
         call COPYDEN(vecBEi,vecBE0,nebf)

! MICROITERATE REGULAR ELECTRONS

         if(LSOSCFA) ITSOA=0
         if(LSOSCFB) ITSOB=0
         if(LSORXCHF) then
          ITSOA=0
          ITSOB=0
          if(rank.eq.0) then
           if(allocated(HESSA0)) HESSA0=0.0d+00
           if(allocated(HESSB0)) HESSB0=0.0d+00
           if(allocated(GRADA0)) GRADA0=0.0d+00
           if(allocated(GRADB0)) GRADB0=0.0d+00
          end if
         end if
         ORBGRDA=0.0d+00
         ORBGRDB=0.0d+00
         if(allocated(PGRADA)) PGRADA=0.0D+00
         if(allocated(PGRADB)) PGRADB=0.0D+00

         if (rank.eq.0) then
          write(*,2001) I
          if((LSOSCFA).and.(LSOSCFB)) then 
           WRITE(*,9052)
          else if ((LSOSCFA).and.(.not.(LSOSCFB))) then
           WRITE(*,9050)
          else if ((LSOSCFB).and.(.not.(LSOSCFA))) then
           WRITE(*,9051)
          else if ((LSORXCHF).and.(LDIIS)) then
           WRITE(*,9054)
          else if (LSORXCHF) then
           WRITE(*,9052)
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
C ARS( microiteration
          E_total=E_HF+E_XCHF+E_int+E_nuc
C )
      end if
C )

! Do OCBSE procedure (restricted solutions for regular and special electrons)

          if (rank.eq.0) then

! Regular electrons

! Transform FAE (calculated at end of previous iteration) to new WA basis
!  - WA constructed from vecAE0 and vecBE0 (previous iteration MOs)
!  - vecAE in AO basis from previous iteration is transformed to new WA basis
!    (relevant for SOSCF only)

           do j=1,nebf-nocca-noccb
             tempvecA(:,j)=vecBE0(:,j+noccb)
           end do
           do j=1,nocca
             tempvecA(:,j+nebf-nocca-noccb)=vecAEi(:,j)
           end do
           if(LCMF) then
            write(*,*) "regular tempvecA"
            call checkovlap(nwbfA,nebf,tempvecA,xxse)
           end if

           call RXCHFmult_OCBSE_transF(nebf,nwbfA,tempvecA,
     x                                 FAE,WA,wFAEw)

!-----------------------POSSIBLE-SOSCF-ALPHA----------------------------(
         if((LSORXCHF).or.(LSOSCFA)) THEN
           ITER=IELEC
           EIGAVL = ITER.GT.1
         end if
         IF((LSORXCHF.or.LSOSCFA) .AND. (EIGAVL)) THEN
           call pack_LT(nwbfA,nwbfLTA,wFAEw,wFLTAw)
           call RXCHFmult_OCBSE_transVt(nebf,nwbfA,WA,
     x                                  xxse,vecAE,wvecAEw)
           if(LCMF) then
            write(*,*) "regular transVt vecAE"
            call checkovlap(nebf,nebf,vecAE,xxse)
           end if
           call SOGRAD(GRADA,wFLTAw,wvecAEw,wWRKAw,NPRA,NA,
     x                 L0wA,L1wA,NwBFLTA,ORBGRDA)
C      write(*,*) "regular grada:",grada
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN
                 WRITE(*,9800)
                 call SOHESS(HSTARTA,wAEenw,NPRA,L0wA,NA,NA)
C                 write(*,*) "regular hstarta:",hstarta
                 if(LSORXCHF) then
                  do j=1,npra
                    HESSA0(j,j)=HSTARTA(j)
                    DISPLIA0(j)=-HSTARTA(j)*GRADA(j)
                  end do
                 end if
              END IF
              ITSOA = ITSOA+1
              IF(LSORXCHF) then
               call RXCHF_SOCHGBAS(NPRA,NA,nwbfA,nebf,WA0,WA,xxse,
     x                             GRADA0,HESSA0,DISPLIA0)
               call RXCHF_SONEWT(NPRA,ITSOA,HESSA0,HESSA,GRADA0,
     x                           GRADA,DISPLIA0,DISPLIA)
C      write(*,*) "grad:",GRADA
C      write(*,*) "grad0:",GRADA0
C      write(*,*) "displia:",DISPLIA
              ELSE
               call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,
     x                     DISPLA,UPDTA,DISPLNA,DGRADIA,UPDTIA,
     x                     ORBGRDA,NPRA,ITSOA,NFT15)
C                 write(*,*) "regular displia:",displia
              END IF
              call SOTRAN(DISPLIA,wvecAEw,wGAw,wWRKAw,NPRA,
     x                    L0wA,L1wA,NA,NA,ORBGRDA)
              IF(LSORXCHF) THEN
                CALL DCOPY(NPRA,GRADA,1,GRADA0,1)
                CALL DCOPY(NPRA,DISPLIA,1,DISPLIA0,1)
                CALL COPYDEN(HESSA0,HESSA,NPRA)
              ELSE
                CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              END IF
              call RXCHFmult_OCBSE_transV(nebf,nwbfA,WA,wvecAEw,wAEenw, ! eigenvalues useless
     x                                    vecAE,AEe)
              if(LCMF) then
               write(*,*) "regular transV vecAE"
               call checkovlap(nebf,nebf,vecAE,xxse)
              end if
              call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)
              GO TO 250
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA----------------------------)

! No SOSCF
!  - Diagonalize Fock matrix in WA basis of this iteration
!  - Obtain updated vecAE in WA basis of this iteration
           call RXCHFmult_OCBSE_diag(nebf,nwbfA,WA,wFAEw,
     x                               wvecAEw,wAEenw,vecAE,AEe)
           if(LCMF) then
            write(*,*) "regular diag vecAE"
            call checkovlap(nebf,nebf,vecAE,xxse)
           end if
           call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)

  250 CONTINUE

           WA0(:,:)=WA(:,:)

          end if ! rank 0


          call MPI_BCAST(VECAE,nebf2_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(AEe,nebf_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(DAE,nebf2_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)

          CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)

          CALL COPYDEN(DAE0,DAE,NEBF)

          CALL COPYDEN(vecAE0,vecAE,NEBF)

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
             npbflt=npbf*(npbf+1)/2
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

!        --> CALCULATE CHANGE IN TOTAL ENERGY
         Delta_E_tot=E_total-E_total_old
         E_total_old=E_total

!        --> PRINT SUMMARY OF THIS ITERATION
         if (rank.eq.0) then
          if((LSOSCFA).and.(LSOSCFB)) then
             WRITE(*,9151) IELEC,E_total,Delta_E_tot,
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
             WRITE(*,9151) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP,ORBGRDA,ORBGRDB
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
           call PREVNU(vecAE,AEE,nebf,nebf,nebf)
           WRITE(*,9620)
           call PREVNU(vecBE,BEE,nebf,nebf,nebf)
           WRITE(*,9700)
           call PREVNU(vecp,EP,npbf,npbf,npbf)
          end if
C )
! Output the vectors for this iteration for restart if necessary:
          call write_MOs(860,nebf,VECAE)
          call write_MOs(861,nebf,VECBE)
          call write_MOs(853,npbf,VECP)
         end if ! rank 0

         LDIFFE=( DIFFAE.LT.TOLE )
         IF(LDIFFE) GOTO 200
         IF(IELEC.EQ.MAXMICROIT) GOTO 10

       END DO  ! microiterations

  200 CONTINUE
!      IF WE GET HERE - MICROITERATION CONVERGENCE ACHIEVED
       if(rank.eq.0) write(*,2000) i,ielec

! MICROITERATE SPECIAL ELECTRONS

         if(LSOSCFA) ITSOA=0
         if(LSOSCFB) ITSOB=0
         if(LSORXCHF) then
          ITSOA=0
          ITSOB=0
          if(rank.eq.0) then
           if(allocated(HESSA0)) HESSA0=0.0d+00
           if(allocated(HESSB0)) HESSB0=0.0d+00
           if(allocated(GRADA0)) GRADA0=0.0d+00
           if(allocated(GRADB0)) GRADB0=0.0d+00
          end if
         end if
         ORBGRDA=0.0d+00
         ORBGRDB=0.0d+00
         if(allocated(PGRADA)) PGRADA=0.0D+00
         if(allocated(PGRADB)) PGRADB=0.0D+00

         if (rank.eq.0) then
          write(*,2001) I
          if((LSOSCFA).and.(LSOSCFB)) then 
           WRITE(*,9052)
          else if ((LSOSCFA).and.(.not.(LSOSCFB))) then
           WRITE(*,9050)
          else if ((LSOSCFB).and.(.not.(LSOSCFA))) then
           WRITE(*,9051)
          else if ((LSORXCHF).and.(LDIIS)) then
           WRITE(*,9054)
          else if (LSORXCHF) then
           WRITE(*,9052)
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
C ARS( microiteration
          E_total=E_HF+E_XCHF+E_int+E_nuc
C )
      end if
C )

! Do OCBSE procedure (restricted solutions for regular and special electrons)

          if (rank.eq.0) then

! Special electrons

! Transform FBE (calculated at end of previous iteration) to new WB basis
!  - WB constructed from vecAE and vecBE0 (current and previous iteration MOs)
!  - vecBE in AO basis from previous iteration is transformed to new WB basis
!    (relevant for SOSCF only)

           do j=1,nebf-nocca-noccb
             tempvecB(:,j)=vecAE(:,j+nocca)
           end do
           do j=1,noccb
             tempvecB(:,j+nebf-nocca-noccb)=vecBEi(:,j)
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
           if(LCMF) then
            write(*,*) "special diag vecBE"
            call checkovlap(nebf,nebf,vecBE,xxse)
           end if
           call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)

  350 CONTINUE

           WB0(:,:)=WB(:,:)

          end if ! rank 0


          call MPI_BCAST(VECBE,nebf2_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(BEe,nebf_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(DBE,nebf2_,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)

          CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)

          CALL COPYDEN(DBE0,DBE,NEBF)

          CALL COPYDEN(vecBE0,vecBE,NEBF)

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

!        --> CALCULATE CHANGE IN TOTAL ENERGY
         Delta_E_tot=E_total-E_total_old
         E_total_old=E_total

!        --> PRINT SUMMARY OF THIS ITERATION
         if (rank.eq.0) then
          if((LSOSCFA).and.(LSOSCFB)) then
             WRITE(*,9151) IELEC,E_total,Delta_E_tot,
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
             WRITE(*,9151) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP,ORBGRDA,ORBGRDB
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
           call PREVNU(vecAE,AEE,nebf,nebf,nebf)
           WRITE(*,9620)
           call PREVNU(vecBE,BEE,nebf,nebf,nebf)
           WRITE(*,9700)
           call PREVNU(vecp,EP,npbf,npbf,npbf)
          end if
C )
C ARS( check orthonormality of MOs
          if(LCMF) then
           write(*,*) "Regular MO overlap matrix:"
           call checkovlap(nebf,nebf,vecAE,xxse)
           write(*,*) "Special MO overlap matrix:"
           call checkovlap(nebf,nebf,vecBE,xxse)
          end if
C )
! Output the vectors for this iteration for restart if necessary:
          call write_MOs(860,nebf,VECAE)
          call write_MOs(861,nebf,VECBE)
          call write_MOs(853,npbf,VECP)
         end if ! rank 0

         LDIFFE=( DIFFBE.LT.TOLE )
         IF(LDIFFE) GOTO 300
         IF(IELEC.EQ.MAXMICROIT) GOTO 10

       END DO  ! microiterations

  300 CONTINUE
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
       call PREVNU(vecAE,AEE,nebf,nebf,nebf)
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
       call write_MOs(860,nebf,VECAE)
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

