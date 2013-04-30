!======================================================================
      subroutine xcrxcuhf(nelec,NAalpE,NAbetE,NBE,NPRA,NPRB,NEBFLT,
     x                    NUCST,npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                    ngtg1,ng1,ng2,ng3,ng4,NG2CHK,NG3CHK,NG4CHK,
     x                    read_CE,read_CP,
     x                    LNEOHF,LGAM4,LCMF,LSOSCF,LOCBSE,
     x                    ng2prm,ng3prm,nat,pmass,cat,zan,bcoef1,gamma1,
     x                    KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                    ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                    SZG2ICR,GM2_1ICR,GM2_2ICR,GM2sICR,
     x                    LG3IC1,SZG3IC1,GM3_1IC1,GM3_2IC1,
     x                    LG4IC,SZG4IC,GM4ICR)

!
!     --- PERFORM A NUCLEAR-ELECTRONIC RESTRICTED XC HARTREE-FOCK
! CALCULATION FOR A N-REGULAR ELECTRON ONE ELECTRON-PROTON SYSTEM ---
!
!     **DEFINITIONS:
!
!     *FOR REGULAR ELECTRONS:
!     DAalpE    ::  NEW REGULAR ALPHA ELECTRON DENSITY MATRIX
!     DAalpE0   ::  OLD REGULAR ALPHA ELECTRON DENSITY MATRIX
!     VECAalpE  ::  REGULAR ALPHA ELECTRON MOS
!     VECAalpE0 ::  OLD REGULAR ALPHA ELECTRON MOS
!     DAbetE    ::  NEW REGULAR BETA ELECTRON DENSITY MATRIX
!     DAbetE0   ::  OLD REGULAR BETA ELECTRON DENSITY MATRIX
!     VECAbetE  ::  REGULAR BETA ELECTRON MOS
!     VECAbetE0 ::  OLD REGULAR BETA ELECTRON MOS
!     AalpEE    ::  REGULAR ALPHA ELECTRON ORBITAL EIGENVALUES
!     AbetEE    ::  REGULAR BETA ELECTRON ORBITAL EIGENVALUES
!
!     * FOR SPECIAL ELECTRON:
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
! Input Variables
      logical LNEOHF
      logical LOCBSE   ! Use OCBSE scheme as is (restricted variational freedom for reg/sp elecs)
      logical LOCBSE2  ! Use modified OCBSE scheme (complete variational freedom for reg elecs)
      logical read_CE
      logical read_CP
      logical LGAM4
      logical LCMF
!     logical LG3DSCF
!     logical LG2IC1
      logical LG3IC1
      logical LG4IC
!     logical LG2DSCF
      integer SZG2ICR
      integer SZG3IC1
      integer SZG4IC
      integer NG2CHK,NG3CHK,NG4CHK
      integer nelec
      integer NBE
      integer NAalpE,NAbetE
      integer NPRA
      integer NPRB
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

      double precision GM2_1ICR(SZG2ICR)
      double precision GM2_2ICR(SZG2ICR)
      double precision GM2sICR(SZG2ICR)
      double precision GM3_1IC1(SZG3IC1)
      double precision GM3_2IC1(SZG3IC1)
      double precision GM4ICR(SZG4IC)

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
      double precision diffAalpE
      double precision diffAbetE
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

      double precision DAalpE(NEBF,NEBF)
      double precision DAalpE0(NEBF,NEBF)
      double precision DAbetE(NEBF,NEBF)
      double precision DAbetE0(NEBF,NEBF)
      double precision DBE(NEBF,NEBF)
      double precision DBE0(NEBF,NEBF)
      double precision VECAalpE(NEBF,NEBF)
      double precision VECAalpE0(NEBF,NEBF)
      double precision VECAbetE(NEBF,NEBF)
      double precision VECAbetE0(NEBF,NEBF)
      double precision VECBE(NEBF,NEBF)
      double precision VECBE0(NEBF,NEBF)
      double precision AalpEE(NEBF)
      double precision AbetEE(NEBF)
      double precision BEE(NEBF)
      double precision FAalpE(nebf,nebf)
      double precision FAbetE(nebf,nebf)
      double precision FBE(nebf,nebf)
!     double precision feDUMMY(nebf,nebf)
!     double precision xeDUMMY(nebf,nebf)

      double precision DP(NPBF,NPBF)
      double precision DP0(NPBF,NPBF)
      double precision VECP(NPBF,NPBF)
      double precision EP(NPBF)
      double precision FP(npbf,npbf)
!     double precision fpDUMMY(npbf,npbf)
!     double precision xpDUMMY(npbf,npbf)

      double precision E_total_old
      double precision Delta_E_tot

      logical LDIFFE

!--------SOSCF-RELATED-VARIABLES------------(
      logical LSOSCF
      logical EIGAVL
      integer NA
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
     * 'ALPHA DENS       BETA DENS        SPEC DENS        ',
     * 'QMP DENS ')

 9050 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        SPEC DENS        ',
     x 'QMP DENS         ORBGRAD_A        ORBGRAD_B')

 9100 FORMAT(1X,I3,F20.10,F17.10,4F17.10)

 9150 FORMAT(1X,I3,F20.10,F17.10,6F17.10)

 9200 FORMAT(/1X,'FINAL NEORXCHF ENERGY IS',F20.10,' AFTER',I4,
     *           ' ITERATIONS')

 9300 FORMAT(/6X,'-----------------------------------------------',/,
     x        6X,'         NEORXCHF ENERGETIC COMPONENTS          ',/,
     x        6X,'-----------------------------------------------',/,
     x       12X,'    E_NUC=',1X,F20.10/
     x       12X,'   E_OMG1=',1X,F20.10/
     x       12X,'   E_OMG2=',1X,F20.10/
     x       12X,'   E_OMG3=',1X,F20.10/
     x       12X,'   E_OMG4=',1X,F20.10/
     x       12X,'   S_OMG1=',1X,F20.10/
     x       12X,'   S_OMG2=',1X,F20.10/
     x       12X,'  S_TOTAL=',1X,F20.10/
     x       12X,'  E_TOTAL=',1X,F20.10/)

 9400 FORMAT(/1X,'          INITIAL GUESS ENERGETICS:            ')

 9500 FORMAT(/6X,' ** BEGIN RXCHF SELF-CONSISTENT-FIELD CALCULATION **')

 9611 FORMAT(/1X,' REGULAR ALPHA ELECTRONIC ORBITALS AND EIGENVALUES: ')

 9612 FORMAT(/1X,' REGULAR BETA ELECTRONIC ORBITALS AND EIGENVALUES:  ')

 9620 FORMAT(/1X,' SPECIAL ELECTRONIC ORBITALS AND EIGENVALUES:       ')

 9700 FORMAT(/1X,'      QM PARTICLE ORBITALS AND EIGENVALUES:         ')

 9800 FORMAT(10X,15(1H-),'START SECOND ORDER SCF',15(1H-))
                                           
!--------OUTPUT-FORMATTING---------------------------------------------)

C User input LOCBSE
C      LOCBSE=.false.
C      LOCBSE2=.true. 
      LOCBSE2=LOCBSE
      if(LOCBSE2) then
       LOCBSE=.false.
       write(*,*) "Using LOCBSE2"
      end if
      if(LOCBSE) write(*,*) "Using LOCBSE"

!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------(
!      call class_nuc_rep(nat,zan,cat,E_nuc)
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
         call RXCUHF_read_CAE("guessCAalpE.inp",nebf,NAalpE,
     x                        DAalpE,VECAalpE0)
         call RXCUHF_read_CAE("guessCAbetE.inp",nebf,NAbetE,
     x                        DAbetE,VECAbetE0)
         call RXCHF_read_CBE(nebf,NBE,DBE,VECBE0)
      else
!        STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
         write(*,*)'ABOUT TO CALL guess_A_elec'
!        call guess_elec(nelec,nebf,xxse,GAM_ecore,DE)
         if ((LOCBSE).or.(LOCBSE2)) then
          call RXCUHF_guess_elec(NAalpE,NAbetE,NBE,nebf,xxse,
     x                           GAM_ecore,DAalpE,DAbetE,DBE,
     x                           VECAalpE0,VECAbetE0,VECBE0)
          write(*,*)'BACK FROM guess_elec for OCBSE'
         else
          call RXCUHF_guess_A_elec(NAalpE,nebf,xxse,GAM_ecore,
     x                             DAalpE,VECAalpE0)
          call RXCUHF_guess_A_elec(NAbetE,nebf,xxse,GAM_ecore,
     x                             DAbetE,VECAbetE0)
          call RXCHF_guess_B_elec(NBE,nebf,xxse,GAM_ecore,DBE,VECBE0)
          write(*,*)'BACK FROM guess_elec'
          write(*,*)
         end if
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
!     LSOSCF=.FALSE.
      if(((NAalpE.eq.1).and.(NAbetE.eq.1)).or.LOCBSE) then
         LSOSCF=.FALSE.
      end if

      if(LSOSCF) THEN
          NFT15=15
          OPEN(NFT15, FILE='WORK15', STATUS='UNKNOWN',
     *         ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
          ITSOA=0
          NFT16=16
          OPEN(NFT16, FILE='WORK16', STATUS='UNKNOWN',
     *         ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
          ITSOB=0
          SOGTOL=1.0d+00
          SMALL=1.0D-06
          L0=nebf
          L1=nebf
      end if
!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------)

!
!     SET CONVERGENCE CRITERIA AND MAXIMUM ITERATIONS 
!
      TOLE = 1.0D-06
      TOLP = 1.0D-04
      maxit=100
      if(LOCBSE) maxit=400
!
!     ZERO OUT 'OLD' DENSITY MATRICES
!
      DAalpE0=0.0d+00
      DAbetE0=0.0d+00
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
       if((.not.(locbse2)).or.(i.eq.1)) then
         call RXCUHF_FOCK(LCMF,LNEOHF,nelec,NBE,NAalpE,NAbetE,
     x                    nebf,nebf2,npbf,npbf2,ngee,
     x                    ng1,ng2,ng3,ng4,
     x                    SZG2ICR,SZG3IC1,SZG4IC,
     x                    NG2CHK,NG3CHK,NG4CHK,
     x                    DAalpE,DAbetE,DBE,DP,
     x                    GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                    GM2_1ICR,GM2_2ICR,GM2sICR,
     x                    GM3_1IC1,GM3_2IC1,
     x                    GM4ICR,
     x                    FP,FAalpE,FAbetE,FBE, 
     x                    E_OMG1,
     x                    E_OMG2,
     x                    E_OMG3,
     x                    E_OMG4,
     x                    S_OMG1,
     x                    S_OMG2,
     x                    E_total,
     x                    S_total)

         E_total=E_total+E_nuc
       end if

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------)
       if(I.eq.1) then
          WRITE(*,9400)
          WRITE(*,9300) E_nuc,E_OMG1,E_OMG2,E_OMG3,E_OMG4,
     x                  S_OMG1,S_OMG2,S_total,E_total
          if(LSOSCF) then 
             WRITE(*,9050)
          else
             WRITE(*,9000)
          end if
       end if
!        Fockp diag
       call UROOTHAN(vecP,EP,xxsp,FP,npbf)
       call construct_DP(nucst,npbf,vecP,DP)

       if (LOCBSE2) then

! Regular electrons
         if(LSOSCF) THEN
          ITER=I
          EIGAVL = ITER.GT.1
         end if

! SOSCF alpha
         IF(LSOSCF .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAalpE,FLT)
          call SOGRAD(GRADA,FLT,vecAalpE,WRK,NPRA,NAalpE,
     x                L0,L1,NEBFLT,ORBGRDA)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 900  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AalpEE,NPRA,L0,NAalpE,NAalpE)
              END IF
              ITSOA = ITSOA+1
           call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,DISPLA,UPDTA,
     *                 DISPLNA,DGRADIA,UPDTIA,ORBGRDA,NPRA,ITSOA,NFT15)
            call SOTRAN(DISPLIA,vecAalpE,GA,WRK,NPRA,
     x                  L0,L1,NAalpE,NAalpE,ORBGRDA)
             CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              call RXCUHF_construct_DAE(NAalpE,nebf,vecAalpE,DAalpE)
              GO TO 950  ! Use the new C's to form new density (change)
            END IF
         END IF

  900 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAalpE,AalpEE,xxse,FAalpE,nebf)
         call RXCUHF_construct_DAE(NAalpE,nebf,vecAalpE,DAalpE)

  950 CONTINUE
!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAalpE0,DAalpE,NEBF,DIFFAalpE)
         CALL COPYDEN(DAalpE0,DAalpE,NEBF)

! SOSCF beta
         IF(LSOSCF .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAbetE,FLT)
          call SOGRAD(GRADB,FLT,vecAbetE,WRK,NPRB,NAbetE,
     x                L0,L1,NEBFLT,ORBGRDB)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 910  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
              IF(ITSOB.EQ.0) THEN
!              WRITE(*,9800)
                 call SOHESS(HSTARTB,AbetEE,NPRB,L0,NAbetE,NAbetE)
              END IF
              ITSOB = ITSOB+1
           call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,DISPLB,UPDTB,
     *                 DISPLNB,DGRADIB,UPDTIB,ORBGRDB,NPRB,ITSOB,NFT16)
            call SOTRAN(DISPLIB,vecAbetE,GB,WRK,NPRB,
     x                  L0,L1,NAbetE,NAbetE,ORBGRDB)
             CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
              call RXCUHF_construct_DAE(NAbetE,nebf,vecAbetE,DAbetE)
              GO TO 960  ! Use the new C's to form new density (change)
            END IF
         END IF

  910 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAbetE,AbetEE,xxse,FAbetE,nebf)
         call RXCUHF_construct_DAE(NAbetE,nebf,vecAbetE,DAbetE)

  960 CONTINUE
!        --> FIND LARGEST CHANGE IN Beta E DENSITY
         CALL DENDIF(DAbetE0,DAbetE,NEBF,DIFFAbetE)
         CALL COPYDEN(DAbetE0,DAbetE,NEBF)

! Do OCBSE procedure for special electrons
         call RXCUHF_OCBSE2(nebf,npebf,NAalpE,NAbetE,NBE,
     x                      vecAalpE,vecAbetE,vecBE0,FBE,xxse,
     x                      elcam,elcbfc,elcex,
     x                      ampeb2c,agebfcc,kpestr,kpeend,
     x                      vecBE,BEe)

! Form special electronic density matrix and store stuff for next it
         call RXCHF_construct_DBE(NBE,nebf,vecBE,DBE)
         CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
         CALL COPYDEN(DBE0,DBE,NEBF)

         CALL COPYDEN(vecAalpE0,vecAalpE,NEBF)
         CALL COPYDEN(vecAbetE0,vecAbetE,NEBF)
         CALL COPYDEN(vecBE0,vecBE,NEBF)

! Calculate energy for this it and Fock matrices for next it
         call RXCUHF_FOCK(LCMF,LNEOHF,nelec,NBE,NAalpE,NAbetE,
     x                    nebf,nebf2,npbf,npbf2,ngee,
     x                    ng1,ng2,ng3,ng4,
     x                    SZG2ICR,SZG3IC1,SZG4IC,
     x                    NG2CHK,NG3CHK,NG4CHK,
     x                    DAalpE,DAbetE,DBE,DP,
     x                    GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                    GM2_1ICR,GM2_2ICR,GM2sICR,
     x                    GM3_1IC1,GM3_2IC1,
     x                    GM4ICR,
     x                    FP,FAalpE,FAbetE,FBE, 
     x                    E_OMG1,
     x                    E_OMG2,
     x                    E_OMG3,
     x                    E_OMG4,
     x                    S_OMG1,
     x                    S_OMG2,
     x                    E_total,
     x                    S_total)

         E_total=E_total+E_nuc

       else      ! if not ocbse2

! SOSCF alpha
         IF(LSOSCF .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAalpE,FLT)
          call SOGRAD(GRADA,FLT,vecAalpE,WRK,NPRA,NAalpE,
     x                L0,L1,NEBFLT,ORBGRDA)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AalpEE,NPRA,L0,NAalpE,NAalpE)
              END IF
              ITSOA = ITSOA+1
           call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,DISPLA,UPDTA,
     *                 DISPLNA,DGRADIA,UPDTIA,ORBGRDA,NPRA,ITSOA,NFT15)
            call SOTRAN(DISPLIA,vecAalpE,GA,WRK,NPRA,
     x                  L0,L1,NAalpE,NAalpE,ORBGRDA)
             CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              call RXCUHF_construct_DAE(NAalpE,nebf,vecAalpE,DAalpE)
              GO TO 750  ! Use the new C's to form new density (change)
            END IF
         END IF

  700 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAalpE,AalpEE,xxse,FAalpE,nebf)
         call RXCUHF_construct_DAE(NAalpE,nebf,vecAalpE,DAalpE)

  750 CONTINUE
!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAalpE0,DAalpE,NEBF,DIFFAalpE)
         CALL COPYDEN(DAalpE0,DAalpE,NEBF)

! SOSCF beta
         IF(LSOSCF .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAbetE,FLT)
          call SOGRAD(GRADB,FLT,vecAbetE,WRK,NPRB,NAbetE,
     x                L0,L1,NEBFLT,ORBGRDB)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 800  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
              IF(ITSOB.EQ.0) THEN
!              WRITE(*,9800)
                 call SOHESS(HSTARTB,AbetEE,NPRB,L0,NAbetE,NAbetE)
              END IF
              ITSOB = ITSOB+1
           call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,DISPLB,UPDTB,
     *                 DISPLNB,DGRADIB,UPDTIB,ORBGRDB,NPRB,ITSOB,NFT16)
            call SOTRAN(DISPLIB,vecAbetE,GB,WRK,NPRB,
     x                  L0,L1,NAbetE,NAbetE,ORBGRDB)
             CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
              call RXCUHF_construct_DAE(NAbetE,nebf,vecAbetE,DAbetE)
              GO TO 850  ! Use the new C's to form new density (change)
            END IF
         END IF

  800 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAbetE,AbetEE,xxse,FAbetE,nebf)
         call RXCUHF_construct_DAE(NAbetE,nebf,vecAbetE,DAbetE)

  850 CONTINUE
!        --> FIND LARGEST CHANGE IN Beta E DENSITY
         CALL DENDIF(DAbetE0,DAbetE,NEBF,DIFFAbetE)
         CALL COPYDEN(DAbetE0,DAbetE,NEBF)

! Do special electrons through direct diagonalization
         call UROOTHAN(vecBE,BEE,xxse,FBE,nebf)
         call RXCHF_construct_DBE(NBE,nebf,vecBE,DBE)

         CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
         CALL COPYDEN(DBE0,DBE,NEBF)

       end if ! end if for not ocbse2

!        --> FIND LARGEST CHANGE IN P DENSITY
       CALL DENDIF(DP0,DP,NPBF,DIFFP)
       CALL COPYDEN(DP0,DP,NPBF)

!        --> CALCULATE CHANGE IN TOTAL ENERGY
       Delta_E_tot=E_total-E_total_old
       E_total_old=E_total

!        --> PRINT SUMMARY OF THIS ITERATION
       if(LSOSCF) then
          WRITE(*,9150) I,E_total,Delta_E_tot,
     x                  DIFFAalpE,DIFFAbetE,DIFFBE,DIFFP,
     x                  ORBGRDA,ORBGRDB
       else
          WRITE(*,9100) I,E_total,Delta_E_tot,
     x                  DIFFAalpE,DIFFAbetE,DIFFBE,DIFFP
       end if
! Output the vectors for this iteration for restart if necessary:
       call write_MOs(870,nebf,VECAalpE)
       call write_MOs(871,nebf,VECAbetE)
       call write_MOs(861,nebf,VECBE)
       call write_MOs(853,npbf,VECP)

       LDIFFE=( (DIFFAalpE.LT.TOLE).and.
     x          (DIFFAbetE.LT.TOLE).and.
     x          (DIFFBE.LT.TOLE) )
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
      WRITE(*,9300) E_nuc,E_OMG1,E_OMG2,E_OMG3,E_OMG4,
     x              S_OMG1,S_OMG2,S_total,E_total

!  OUTPUT ELEC AND NUC EIGENVALUES AND EIGENVECTORS
      WRITE(*,9611)
      call PREVNU(vecAalpE,AalpEE,nebf,nebf,nebf)
      WRITE(*,9612)
      call PREVNU(vecAbetE,AbetEE,nebf,nebf,nebf)
      WRITE(*,9620)
      call PREVNU(vecBE,BEE,nebf,nebf,nebf)
      WRITE(*,9700)
      call PREVNU(vecp,EP,npbf,npbf,npbf)

! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------(
!     subroutine write_MOs(IFIL,nbf,VEC)
C     IFIL=853 :: FinalCP.dat
C     IFIL=861 :: FinalCBE.dat
C     IFIL=870 :: FinalCAalpE.dat
C     IFIL=871 :: FinalCAbetE.dat
      call write_MOs(870,nebf,VECAalpE)
      call write_MOs(871,nebf,VECAbetE)
      call write_MOs(861,nebf,VECBE)
      call write_MOs(853,npbf,VECP)
! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------)
!

      RETURN
      END
!======================================================================
      SUBROUTINE RXCUHF_guess_A_elec(NAE,nebf,xxse,GAM_ecore,DAE,CAE)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial regular electron guess density
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

      call RXCUHF_construct_DAE(NAE,nebf,CAE,DAE)


      RETURN
      END
!======================================================================
      subroutine RXCUHF_guess_elec(NAalpE,NAbetE,NBE,nebf,xxse,
     x                             GAM_ecore,DAalpE,DAbetE,DBE,
     x                             CAalpE,CAbetE,CBE)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial regular and special electronic guess density
!
!     Regular guess for max{NAalpE,NAbetE} regular electrons
!     Special guess determined from orthogonal space of min{NAalpE,NAbetE}
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer NAalpE,NAbetE
      integer NBE
      double precision xxse(nebf,nebf)
      double precision GAM_ecore(nebf,nebf)
! Variables Returned
      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision CAalpE(nebf,nebf)
      double precision CAbetE(nebf,nebf)
      double precision CBE(nebf,nebf)
! Local variables
      integer i,j
      integer nocca,noccb
      double precision C(nebf,nebf)
      double precision EVF(nebf)
      double precision zero
      parameter(zero=0.0d+00)

      DAalpE=zero
      DAbetE=zero
      DBE=zero
      CAalpE=zero
      CAbetE=zero
      CBE=zero
      noccb=1

      call UROOTHAN(C,EVF,xxse,GAM_ecore,nebf)

C Special electron has spin alpha
C   => It should be initialized to be orthogonal to alpha regular electrons

! Normal guess for beta regular electrons
      call RXCUHF_guess_A_elec(NAbetE,nebf,xxse,GAM_ecore,
     x                         DAbetE,CAbetE)

      nocca=NAalpE

! Modified guess for alpha regular electrons
! Store first nocca evectors as occ reg elec vectors
      do i=1,nocca
        do j=1,nebf
          CAalpE(j,i)=C(j,i)
        end do
      end do
      call RXCUHF_construct_DAE(NAalpE,nebf,CAalpE,DAalpE)

! Store nebf-nocca remaining evectors as occ and virt spec elec vectors
       do i=nocca+1,nebf
         do j=1,nebf
           CBE(j,i-nocca)=C(j,i)
         end do
       end do

      call RXCHF_construct_DBE(NBE,nebf,CBE,DBE)

      return
      end
!======================================================================
      subroutine RXCUHF_OCBSE2(nebf,npebf,NAalpE,NAbetE,NBE,
     x                         vecAalpE,vecAbetE,vecBE0,FBE,Selec,
     x                         elcam,elcbfc,elcex,
     x                         ampeb2c,agebfcc,kpestr,kpeend,
     x                         vecBE,BEen)
! 
!     Perform OCBSE procedure for special electrons only
!       vecBE0:    Special electron coefficients from previous iteration
!       vecAalpE:  Regular electron coefficients from current iteration
!       vecAbetE:  Regular electron coefficients from current iteration
!       vecBE:     Special electron coefficients from current iteration
!     Expands special electronic solution in orthogonal space to
!       NAalpE regular electrons (special electron has spin alpha)
!======================================================================
      implicit none
! Input Variables
      integer nebf,npebf
      integer NAalpe,NAbetE,NBE
      double precision vecAalpE(nebf,nebf),vecAbetE(nebf,nebf)
      double precision vecBE0(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
      integer ampeb2c(npebf)               ! Map prim index to contr index
      integer kpestr(nebf)                 ! Map contr index to prim start
      integer kpeend(nebf)                 ! Map contr index to prim end
      integer elcam(npebf,3)               ! Angular mom for electrons
      double precision agebfcc(npebf)      ! Map prim index to contr coeff
      double precision elcex(npebf)        ! Exponents: elec basis
      double precision elcbfc(npebf,3)     ! Basis centers: elec basis
! Variables Returned
      double precision vecBE(nebf,nebf)
      double precision BEen(nebf)
! Local variables
      integer nae
      integer nocca,noccb,nvirt 
      integer noccvirta,noccvirtb
      double precision vecAE(nebf,nebf)
      double precision zero
      parameter(zero=0.0d+00)

      noccb=1

C Special electron has spin alpha
C   => It should be constrained to be orthogonal to alpha regular electrons
      nae=NAalpE
      nocca=NAalpE
      vecAE=vecAalpE

      nvirt=nebf-nocca-noccb
      noccvirta=nebf-noccb
      noccvirtb=nebf-nocca

      vecBE=zero
      BEen=zero

      call RXCHF_OCBSE_driver2(nebf,npebf,nae,nbe,nocca,noccb,
     x                         nvirt,noccvirta,noccvirtb,
     x                         vecAE,vecBE0,FBE,Selec,
     x                         elcam,elcbfc,elcex,
     x                         ampeb2c,agebfcc,kpestr,kpeend,
     x                         vecBE,BEen)

      return
      end
