!======================================================================
      subroutine xcuscf(nelec,NAE,NBE,NPRA,NPRB,NEBFLT,NUCST,
     x                  npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                  ngtg1,ng1,ng2,ng3,ng4,NG4CHK,NG3CHK,NG2CHK,
     x                  read_CE,read_CP,
     x                  LNEOHF,LGAM4,LCMF,LSOSCF,
     x                  ng2prm,ng3prm,nat,pmass,cat,zan,bcoef1,gamma1,
     x                  KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                  ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                  SZG2ICR,GM2ICR,GM2SICR,
     x                  SZG3IC1,GM3IC1,
     x                  SZG4IC,GM4ICR)

!
!    x                  GM2ICR,GM2SICR,GM3IC1,GM4ICR)
!
!     --- PERFORM A NUCLEAR-ELECTRONIC XC UNRESTRICTED HARTREE-FOCK
!     CALCULATION FOR A N-ELECTRON ONE-PROTON SYSTEM ---
!
!     **DEFINITIONS:
!
!     *FOR ELECTRONS:
!     DAE    ::  ALPHA NEW ELECTRON DENSITY MATRIX
!     DAE0   ::  ALPHA OLD ELECTRON DENSITY MATRIX
!     VECAE  ::  ALPHA ELECTRON MOS
!     AEE    ::  ALPHA ELECTRON ORBITAL EIGENVALUES
!
!     *FOR PROTONS:
!     DP    ::  NEW PROTON DENSITY MATRIX
!     DP0   ::  OLD PROTON DENSITY MATRIX
!     VECP  ::  PROTON MOS
!     EP    ::  PROTON ORBITAL EIGENVALUES
!
!======================================================================
      implicit none
! Input Variables
      logical LNEOHF
      logical read_CE
      logical read_CP
      logical LGAM4
      logical LCMF
!     logical LG3DSCF
!     logical LG2IC1
!     logical LG3IC1
!     logical LG4IC
!     logical LG2DSCF
      integer SZG2ICR
      integer SZG3IC1
      integer SZG4IC
      integer NG2CHK
      integer NG3CHK
      integer NG4CHK
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
      integer ng1
      integer ng2
      integer ng3
      integer ng4
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

      double precision GM2ICR(SZG2ICR)  ! POSSIBLE IN-CORE GAM2 INTEGRALS
      double precision GM2SICR(SZG2ICR) ! POSSIBLE IN-CORE GAM2 INTEGRALS
      double precision GM3IC1(SZG3IC1)  ! POSSIBLE IN-CORE GAM3 INTEGRALS
      double precision GM4ICR(SZG4IC)   ! POSSIBLE IN-CORE GAM4 INTEGRALS

!     double precision GM2ICR(ng2)  ! POSSIBLE IN-CORE GAM2 INTEGRALS
!     double precision GM2SICR(ng2) ! POSSIBLE IN-CORE GAM2 INTEGRALS
!     double precision GM3IC1(ng3)  ! POSSIBLE IN-CORE GAM3 INTEGRALS
!     double precision GM4ICR(ng4)   ! POSSIBLE IN-CORE GAM4 INTEGRALS

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
      integer maxit

      integer i
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

      double precision E_ecore
      double precision E_pcore
      double precision E_ee
      double precision E_ep
      double precision E_nuc
      double precision E_total

      double precision E_OMG1
      double precision E_OMG2
      double precision E_OMG3
      double precision E_OMG4
      double precision E_UHF

!     double precision S_gam1s
      double precision S_total
!     double precision S_OMG0
      double precision S_OMG1
      double precision S_OMG2

      double precision DAE(NEBF,NEBF)
      double precision DAE0(NEBF,NEBF)
      double precision DBE(NEBF,NEBF)
      double precision DBE0(NEBF,NEBF)
      double precision VECAE(NEBF,NEBF)
      double precision VECBE(NEBF,NEBF)
      double precision AEE(NEBF)
      double precision BEE(NEBF)
      double precision FAE(nebf,nebf)
      double precision XFAE(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision XFBE(nebf,nebf)
!     double precision feDUMMY(nebf,nebf)
!     double precision xeDUMMY(nebf,nebf)

      double precision DP(NPBF,NPBF)
      double precision DP0(NPBF,NPBF)
      double precision VECP(NPBF,NPBF)
      double precision EP(NPBF)
      double precision FP(npbf,npbf)
      double precision XFP(npbf,npbf)
!     double precision fpDUMMY(npbf,npbf)
!     double precision xpDUMMY(npbf,npbf)

      double precision E_total_old
      double precision Delta_E_tot

      logical LDIFFE

!--------SOSCF-RELATED-VARIABLES------------(
      logical LSOSCF
      logical EIGAVL
!     integer NA
      integer ITER
      integer ITSOA ! SOSCF iteration counter
      integer ITSOB ! SOSCF iteration counter
      integer L0,L1
      integer NFT15
      integer NFT16
      double precision FLT(NEBFLT) !FLT: Lower triangle focke
      double precision HSTARTA(NPRA)
      double precision HSTARTB(NPRB)
      double precision GRADA(NPRA)
      double precision GRADB(NPRB)
      double precision PGRADA(NPRA)
      double precision PGRADB(NPRB)
      double precision DISPLIA(NPRA)
      double precision DISPLIB(NPRB)
      double precision DGRADA(NPRA)  ! WRK1
      double precision DGRADB(NPRB)  ! WRK1
      double precision DISPLA(NPRA)  ! WRK2
      double precision DISPLB(NPRB)  ! WRK2
      double precision UPDTA(NPRA)   ! WRK3
      double precision UPDTB(NPRB)   ! WRK3
      double precision DISPLNA(NPRA) ! WRK1+NPR
      double precision DISPLNB(NPRB) ! WRK1+NPR
      double precision DGRADIA(NPRA) ! WRK2+NPR
      double precision DGRADIB(NPRB) ! WRK2+NPR
      double precision UPDTIA(NPRA)  ! WRK3+NPR
      double precision UPDTIB(NPRB)  ! WRK3+NPR
      double precision ORBGRDA
      double precision ORBGRDB
      double precision SMALL
      double precision SOGTOL ! ORBGRAD TOL to activate soscf
      double precision XA(NPRA)
      double precision XB(NPRB)
      double precision GA(nebf,nebf) !G(L0,L0)
      double precision GB(nebf,nebf) !G(L0,L0)
      double precision WRK(nebf) !WRK(L0)
!cc   double precision CCC(nebf,nebf) !WRK(L0)
!cc   NPR=(L0-NA)*NA ! Line 2134 RHFCL ?NA is NUM ALPHA E?
!--------SOSCF-RELATED-VARIABLES------------)

!--------OUTPUT-FORMATTING---------------------------------------------(
 9000 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS ')

 9050 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS         ',
     * 'ORBGRAD_A        ORBGRAD_B')

 9100 FORMAT(1X,I3,F20.10,F17.10,3F17.10)

 9150 FORMAT(1X,I3,F20.10,F17.10,5F17.10)

 9200 FORMAT(/1X,'FINAL NEOXCHF ENERGY IS',F20.10,' AFTER',I4,
     *           ' ITERATIONS')

 9300 FORMAT(/6X,'-----------------------------------------------',/,
     x        6X,'         NEOXCHF ENERGETIC COMPONENTS          ',/,
     x        6X,'-----------------------------------------------',/,
     x       12X,'    E_NUC=',1X,F20.10/
     x       12X,'  E_ECORE=',1X,F20.10/
     x       12X,'  E_PCORE=',1X,F20.10/
     x       12X,'     E_EP=',1X,F20.10/
     x       12X,'     E_EE=',1X,F20.10/
     x       12X,'   E_OMG1=',1X,F20.10/
     x       12X,'   E_OMG2=',1X,F20.10/
     x       12X,'   E_OMG3=',1X,F20.10/
     x       12X,'   E_OMG4=',1X,F20.10/
     x       12X,'  S_TOTAL=',1X,F20.10/
     x       12X,'  E_TOTAL=',1X,F20.10/)

 9400 FORMAT(/1X,'          INITIAL GUESS ENERGETICS:            ')

 9500 FORMAT(/6X,' ** BEGIN NEOXC SELF-CONSISTENT-FIELD CALCULATION **')

 9610 FORMAT(/1X,' ALPHA ELECTRONIC ORBITALS AND EIGENVALUES:         ')

 9620 FORMAT(/1X,'  BETA ELECTRONIC ORBITALS AND EIGENVALUES:         ')

 9700 FORMAT(/1X,'      QM PARTICLE ORBITALS AND EIGENVALUES:         ')

 9800 FORMAT(10X,15(1H-),'START SECOND ORDER SCF',15(1H-))
                                           
!--------OUTPUT-FORMATTING---------------------------------------------)


!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------(
!     call class_nuc_rep(nat,zan,cat,E_nuc)
!     write(*,*)
!     write(*,*)'IN SCF: E_nuc=',E_nuc
!     write(*,*)
      open(800,file='ENUCRP.dat',status='unknown')
      read(800,*) E_nuc
      close(800)
      write(*,*)'READ IN NUC REPULSION'
!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------)

!--------------READ-INTEGRALS-NEEDED-FOR-NEO-HF------------------------(
      nebf2=nebf*nebf
      npbf2=npbf*npbf
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

!     write(*,*)
!     write(*,*)'IN: xcuscf '
!     write(*,*)'AFTER READ IN GAM_EE='
!     write(*,*)GAM_ee
!     write(*,*)

!--------------READ-INTEGRALS-NEEDED-FOR-NEO-HF------------------------)

!-------------INITIAL-GUESSES------------------------------------------(
      if(read_CE) then
!        READ IN GUESS FOR E:
!        call read_elec_density(nebf,nelec,DE)
         call read_CAE(nebf,NAE,DAE)
         call read_CBE(nebf,NBE,DBE)
      else
!        STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
         write(*,*)'ABOUT TO CALL guess_A_elec'
!        call guess_elec(nelec,nebf,xxse,GAM_ecore,DE)
         call guess_A_elec(NAE,nebf,xxse,GAM_ecore,DAE)
         call guess_A_elec(NBE,nebf,xxse,GAM_ecore,DBE)
         write(*,*)'BACK FROM guess_elec'
         write(*,*)
      end if
      if(read_CP) then
!        READ IN GUESS FOR N:
         call read_nuc_density(npbf,1,NUCST,DP)
      else
!        STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
         write(*,*)'ABOUT TO CALL guess_prot'
!        call guess_prot(NUCST,npbf,nebf,xxsp,GAM_pcore,GAM_ep,DE,DP)
         call guess_prot2(NUCST,npbf,xxsp,GAM_pcore,DP)
         write(*,*)'BACK FROM guess_prot'
         write(*,*)
      end if
!-------------INITIAL-GUESSES------------------------------------------)

!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------(
!     SOSCF=.FALSE.
!     SOSCF=.TRUE.
      if(nelec.eq.1) then
         LSOSCF=.FALSE.
      end if

      if(LSOSCF) THEN
         NFT15=15
         OPEN(NFT15, FILE='WORK15', STATUS='UNKNOWN',
     *        ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         NFT16=16
         OPEN(NFT16, FILE='WORK16', STATUS='UNKNOWN',
     *        ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         SOGTOL=1.0d+00
         SMALL=1.0D-06
         ITSOA=0
         ITSOB=0
         L0=nebf
         L1=nebf
!        NA=nelec/2
      end if
!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------)

!
!     SET CONVERGENCE CRITERIA AND MAXIMUM ITERATIONS 
!
      TOLE = 1.0D-06
      TOLP = 1.0D-04
      maxit=100
!
!     ZERO OUT 'OLD' DENSITY MATRICES
!
      DAE0=0.0d+00
      DBE0=0.0d+00
      DP0=0.0d+00
!
!     BEGIN XCSCF ITERATIONS
      WRITE(*,9500)
!     WRITE(*,9000)
!
      E_total_old=0.0d+00
      ORBGRDA=0.0d+00
      ORBGRDB=0.0d+00

      DO I=1,MAXIT

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------(
!     write(*,*)
!     write(*,*)'IN: xcuscf '
!     write(*,*)'Before call to UHF_FOCK, GAM_EE='
!     write(*,*)GAM_ee
!     write(*,*)
         call UHF_FOCK(LCMF,LNEOHF,nelec,NAE,NBE,
     x                    nebf,nebf2,npbf,npbf2,ngee,
     x                    ng1,ng2,ng3,ng4,
     x                    SZG2ICR,SZG3IC1,SZG4IC,
     x                    NG2CHK,NG3CHK,NG4CHK,
     x                    DAE,DBE,DP,
     x                    GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                    GM2ICR,GM2sICR,GM3IC1,GM4ICR,
     x                    FP,FAE,FBE,
     x                    E_pcore,
     x                    E_ecore,
     x                    E_ep,
     x                    E_ee,
     x                    E_OMG1,
     x                    E_OMG2,
     x                    E_OMG3,
     x                    E_OMG4,
     x                    S_OMG1,
     x                    S_OMG2,
     x                    E_UHF,
     x                    E_total,
     x                    S_total)

         E_total=E_total+E_nuc

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------)
         if(I.eq.1) then
            WRITE(*,9400)
            WRITE(*,9300) E_nuc,E_ecore,E_pcore,E_ep,E_ee,
     x      E_OMG1,E_OMG2,E_OMG3,E_OMG4,S_total,E_total
            if(LSOSCF) then 
               WRITE(*,9050)
            else
               WRITE(*,9000)
            end if
         end if
!        Fockp diag
         call UROOTHAN(vecP,EP,xxsp,FP,npbf)
         call construct_DP(nucst,npbf,vecP,DP)

!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------(
         if(LSOSCF) THEN
         ITER=I
         EIGAVL = ITER.GT.1
         end if
         IF(LSOSCF .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAE,FLT)
          call SOGRAD(GRADA,FLT,vecAE,WRK,NPRA,NAE,L0,L1,NEBFLT,ORBGRDA)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AEE,NPRA,L0,NAE,NAE)
              END IF
              ITSOA = ITSOA+1
           call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,DISPLA,UPDTA,
     *                 DISPLNA,DGRADIA,UPDTIA,ORBGRDA,NPRA,ITSOA,NFT15)
            call SOTRAN(DISPLIA,vecAE,GA,WRK,NPRA,L0,L1,NAE,NAE,ORBGRDA)
             CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              call construct_DAE(NAE,nebf,vecAE,DAE)
              GO TO 750  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------)

  700 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAE,AEE,xxse,FAE,nebf)
         call construct_DAE(NAE,nebf,vecAE,DAE)

  750 CONTINUE
!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
         CALL COPYDEN(DAE0,DAE,NEBF)

!-----------------------POSSIBLE-SOSCF-BETA----------------------------(
!        if(LSOSCF) THEN
!        ITER=I
!        EIGAVL = ITER.GT.1
         IF(LSOSCF .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FBE,FLT)
          call SOGRAD(GRADB,FLT,vecBE,WRK,NPRB,NBE,L0,L1,NEBFLT,ORBGRDB)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
              IF(ITSOB.EQ.0) THEN
!             WRITE(*,9800)
                 call SOHESS(HSTARTB,BEE,NPRB,L0,NBE,NBE)
              END IF
              ITSOB = ITSOB+1
           call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,DISPLB,UPDTB,
     *                 DISPLNB,DGRADIB,UPDTIB,ORBGRDB,NPRB,ITSOB,NFT16)
            call SOTRAN(DISPLIB,vecBE,GB,WRK,NPRB,L0,L1,NBE,NBE,ORBGRDB)
             CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
              call construct_DAE(NBE,nebf,vecBE,DBE)
              GO TO 850  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-BETA----------------------------)

  800 CONTINUE
!        call ROOTHAN(DBE,vecBE,BEE,xxse,FBE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecBE,BEE,xxse,FBE,nebf)
         call construct_DAE(NBE,nebf,vecBE,DBE)

  850 CONTINUE
!        --> FIND LARGEST CHANGE IN Beta E DENSITY
         CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
         CALL COPYDEN(DBE0,DBE,NEBF)

!        --> FIND LARGEST CHANGE IN P DENSITY
         CALL DENDIF(DP0,DP,NPBF,DIFFP)
         CALL COPYDEN(DP0,DP,NPBF)

!        --> CALCULATE CHANGE IN TOTAL ENERGY
         Delta_E_tot=E_total-E_total_old
         E_total_old=E_total

!        --> PRINT SUMMARY OF THIS ITERATION
         if(LSOSCF) then
            WRITE(*,9150) I,E_total,Delta_E_tot,DIFFAE,DIFFBE,DIFFP,
     xORBGRDA,ORBGRDB
         else
            WRITE(*,9100) I,E_total,Delta_E_tot,DIFFAE,DIFFBE,DIFFP
         end if
! Output the vectors for this iteration for restart if necessary:
         call write_MOs(860,nebf,VECAE)
         call write_MOs(861,nebf,VECBE)
         call write_MOs(853,npbf,VECP)

         LDIFFE=( (DIFFAE.LT.TOLE).and.(DIFFBE.LT.TOLE) )
         if(npbf.gt.1) then
            IF((LDIFFE).and.(DIFFP.LT.TOLP)) GOTO 100
         else
            IF(LDIFFE) GOTO 100
         end if
         IF(I.EQ.MAXIT) GOTO 10
      END DO
 
  10  CONTINUE
!     IF WE GET HERE SOMETHING WENT WRONG

      if(LSOSCF) THEN
         close(NFT15)
         close(NFT16)
      end if

      WRITE(*,*)
      WRITE(*,*)'WARNING:  ITERATION LIMIT EXCEEDED'
      E_TOTAL=zero
      WRITE(*,*)
!     STOP
!
  100 CONTINUE
!     IF WE GET HERE WE ARE DONE - CONVERGENCE ACHIEVED

      if(LSOSCF) THEN
         close(NFT15)
         close(NFT16)
      end if

!     PRINT FINAL ENERGY AND PUNCH THE ORBITALS
      WRITE(*,9200) E_TOTAL,I
!
      WRITE(*,9300) E_nuc,E_ecore,E_pcore,E_ep,E_ee,
     xE_OMG1,E_OMG2,E_OMG3,E_OMG4,S_total,E_total

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
!

      RETURN
      END
!======================================================================
      SUBROUTINE guess_A_elec(NAE,nebf,xxse,GAM_ecore,DAE)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial electron guess density
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer NAE
      double precision xxse(nebf,nebf)
      double precision GAM_ecore(nebf,nebf)
! Variables Returned
      double precision DAE(nebf,nebf)
! Local variables
      double precision CAE(nebf,nebf)
      double precision EVF(nebf)


      call UROOTHAN(CAE,EVF,xxse,GAM_ecore,nebf)

      call construct_DAE(NAE,nebf,CAE,DAE)


      RETURN
      END
!======================================================================
      SUBROUTINE guess_prot2(NUCST,npbf,xxsp,GAM_pcore,DP)
!
!     Diagonalize the core proton Hamiltonian- 
!     to construct initial proton guess density
!======================================================================
      implicit none
! Input Variables
      integer npbf
      integer nucst
      double precision xxsp(npbf,npbf)
      double precision GAM_pcore(npbf,npbf)
! Variables Returned
      double precision DP(npbf,npbf)
! Local variables
      double precision fockp(npbf,npbf)
      double precision CP(npbf,npbf)
      double precision EVF(npbf)
      double precision E_pcore

      fockp=0.0d+00
      call E_from_GAM_pcore(npbf,npbf*npbf,
     * GAM_pcore,DP,fockp,E_pcore)

      call UROOTHAN(CP,EVF,xxsp,fockp,npbf)

      call construct_DP(nucst,npbf,CP,DP)


      RETURN
      END
C======================================================================
      SUBROUTINE UROOTHAN(C,EVF,S,F,NB)
C
C     SOLVES THE ROOTHAN EQUATIONS:
C     ==============================
C     ON INPUT:
C     S         ::  OVERLAP MATRIX
C     F         ::  FOCK MATRIX
C     NB        ::  NUMBER OF BASIS FUNCTIONS
C     ==============================
C     WORKING:
C     EV        ::  EIGENVALUES OF S MATRIX
C     EVECS     ::  EIGENVECTORS OF S MATRIX
C     EVECST    ::  TRANSPOSE OF EIGENVECTORS OF S MATRIX
C     X         ::  TRANSFORMATION MATRIX
C     XP        ::  TRANSPOSE OF TRANFORMATION MATRIX
C     FP        ::  TRANSFORMED FOCK MATRIX
C     CP        ::  EIGENVECTORS FROM FOCK MATRIX
C     FV1       ::  WORK SPACE
C     FV2       ::  WORK SPACE
C     FV3       ::  WORK SPACE
C     FV4       ::  WORK SPACE
C     ==============================
C     OUTPUT:
C     C         ::  TRANSFORMED EIGENVECTORS
C     EVF       ::  EIGENVALUES FROM FOCK MATRIX
C     ============================== 
C
C======================================================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00)
      DIMENSION S(NB,NB),F(NB,NB),FP(NB,NB),EVF(NB),X(NB,NB),XP(NB,NB),
     *          CP(NB,NB),C(NB,NB),D(NB,NB),EVS(NB),EVECS(NB,NB),
     *          FV1(NB),FV2(NB),FV3(NB,NB),FV4(NB,NB),DEVS(NB,NB),
     *          EVECST(NB,NB),xmat(NB,NB)
      LOGICAL DEBUG
C
C     CONSTRUCT ORTHONORMALIZING TRANSFORMATION MATRIX
C
C     ---> DIAGONALIZE OVERLAP MATRIX
C
C     DEBUG=.TRUE.
      DEBUG=.FALSE.
      CALL RS(NB,NB,S,EVS,2,EVECS,FV1,FV2,IERR)
      IF(DEBUG) THEN
         WRITE(*,*)
         WRITE(*,*)'---- OVERLAP MATRIX --> EIGENVALUES ----'
         WRITE(*,*)
         DO I=1,NB
               WRITE(*,*)'EVS(',I,')=',EVS(I)
         END DO
         WRITE(*,*)
         WRITE(*,*)'---- OVERLAP MATRIX --> EIGENVECTORS ----'
         WRITE(*,*)
         DO I=1,NB
            DO J=1,NB
               WRITE(*,*)'EVECS(',I,J,')=',EVECS(I,J)
            END DO
         END DO
      END IF
C           
C     ---- SYMMETRIC ORTHOGONALIZATION SCHEME ----
C
C     ---> FORM DIAGONAL EIGENVALUE MATRIX
C     ---> GET TRANSPOSE OF VECTORS FROM DIAGONALIZATION OF S
C
      DO I=1,NB
         DO J=1,NB
            DEVS(I,J)=ZERO
            EVECST(I,J)=EVECS(J,I)
         END DO
         DEVS(I,I)=1/SQRT(EVS(I))
      END DO
      IF(DEBUG) THEN
         WRITE(*,*)
         WRITE(*,*)'---- S^-(1/2) EIGENVALUE MATRIX ----'
         WRITE(*,*)
         DO I=1,NB
            DO J=1,NB
               WRITE(*,*)'DEVS(',I,J,')=',DEVS(I,J)
            END DO
         END DO
      END IF
C
C     ---> FORM TRANSFORMATION MATRIX 
C
C     MULTIPLY EVECS.DEVS.EVECST = X
C     
      CALL MATMULT(NB,NB,NB,NB,DEVS,EVECST,FV3)
      CALL MATMULT(NB,NB,NB,NB,EVECS,FV3,X)
      IF(DEBUG) THEN
         WRITE(*,*)
         WRITE(*,*)'---- X MATRIX ----'
         WRITE(*,*)
         DO I=1,NB
            DO J=1,NB
               WRITE(*,*)'X(',I,J,')=',X(I,J)
            END DO
         END DO
      END IF
C
C     FORM TRANSPOSE OF TRANSFORMATION MATRIX
C
      DO I=1,NB
         DO J=1,NB
            XP(I,J)=X(J,I)
         END DO
      END DO
C
C     TRANSFORM FOCK MATRIX 
C 
C     MULTIPLY XP.F.X = FP
C
      CALL MATMULT(NB,NB,NB,NB,F,X,FV4)
      CALL MATMULT(NB,NB,NB,NB,XP,FV4,FP)
C
      IF(DEBUG) THEN
         WRITE(*,*)
         WRITE(*,*) "FOCK BEFORE DIAGONALIZING:"
         CALL PREVNU(FP,EVS,NB,NB,NB)
      END IF
C
C     DIAGONALIZE FOCK MATRIX 
C
      CALL RS(NB,NB,FP,EVF,2,CP,FV1,FV2,IERR)
C
C     TRANSFORM MO COEFFICIENTS
C      C = X.CP  
C   
      CALL MATMULT(NB,NB,NB,NB,X,CP,C)
C
CDDD
c     write(*,*) 
c     write(*,*)'--------ORTHO-CHEC----------' 
c     do i=1,NB
c     do j=1,NB
c     xmat(i,j)=0.0d0
c     do k1=1,NB
c     do k2=1,NB
c     xmat(i,j)=xmat(i,j)+(C(k1,i)*S(k1,k2)*C(k2,j))
c     end do
c     end do
c     write(*,*) i,j,xmat(i,j)
c     end do
c     end do
c     write(*,*)'--------END OF ORTHO-CHEC----------' 
C
      RETURN
      END
