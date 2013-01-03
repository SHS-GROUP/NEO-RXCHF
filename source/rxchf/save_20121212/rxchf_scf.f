!======================================================================
      subroutine xcrxchf(nelec,NAE,NBE,NPRA,NEBFLT,NUCST,
     x                   npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                   ngtg1,ng1,ng2,ng3,ng4,NG2CHK,NG3CHK,NG4CHK,
     x                   read_CE,read_CP,
     x                   LNEOHF,LGAM4,LCMF,LSOSCF,LOCBSE,
     x                   ng2prm,ng3prm,nat,pmass,cat,zan,bcoef1,gamma1,
     x                   KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                   ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                   SZG2ICR,GM2_1ICR,GM2_2ICR,GM2sICR,
     x                   LG3IC1,SZG3IC1,GM3_1IC1,GM3_2IC1,
     x                   LG4IC,SZG4IC,GM4ICR)

!
!     --- PERFORM A NUCLEAR-ELECTRONIC RESTRICTED XC HARTREE-FOCK
! CALCULATION FOR A N-REGULAR ELECTRON ONE ELECTRON-PROTON SYSTEM ---
!
!     **DEFINITIONS:
!
!     *FOR REGULAR ELECTRONS:
!     DAE    ::  NEW REGULAR ELECTRON DENSITY MATRIX
!     DAE0   ::  OLD REGULAR ELECTRON DENSITY MATRIX
!     VECAE  ::  REGULAR ELECTRON MOS
!     VECAE0 ::  OLD REGULAR ELECTRON MOS
!     AEE    ::  REGULAR ELECTRON ORBITAL EIGENVALUES
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
      integer NAE,NBE
      integer NPRA
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
      double precision E_exch
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
      double precision VECAE0(NEBF,NEBF)
      double precision VECBE(NEBF,NEBF)
      double precision VECBE0(NEBF,NEBF)
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
      integer NA
      integer ITER
      integer ITSOA ! SOSCF iteration counter
      integer ITSOB ! SOSCF iteration counter
      integer L0,L1
      integer NFT15
      integer NFT16
      double precision FLT(NEBFLT) !FLT: Lower triangle focke
      double precision HSTARTA(NPRA)
C      double precision HSTARTB(NPRB)
      double precision GRADA(NPRA)
C      double precision GRADB(NPRB)
      double precision PGRADA(NPRA)
C      double precision PGRADB(NPRB)
      double precision DISPLIA(NPRA)
C      double precision DISPLIB(NPRB)
      double precision DGRADA(NPRA)  ! WRK1
C      double precision DGRADB(NPRB)  ! WRK1
      double precision DISPLA(NPRA)  ! WRK2
C      double precision DISPLB(NPRB)  ! WRK2
      double precision UPDTA(NPRA)   ! WRK3
C      double precision UPDTB(NPRB)   ! WRK3
      double precision DISPLNA(NPRA) ! WRK1+NPR
C      double precision DISPLNB(NPRB) ! WRK1+NPR
      double precision DGRADIA(NPRA) ! WRK2+NPR
C      double precision DGRADIB(NPRB) ! WRK2+NPR
      double precision UPDTIA(NPRA)  ! WRK3+NPR
C      double precision UPDTIB(NPRB)  ! WRK3+NPR
      double precision ORBGRDA
C      double precision ORBGRDB
      double precision SMALL
      double precision SOGTOL ! ORBGRAD TOL to activate soscf
      double precision XA(NPRA)
C      double precision XB(NPRB)
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
     * 'ORBGRAD_A ')

 9100 FORMAT(1X,I3,F20.10,F17.10,3F17.10)

 9150 FORMAT(1X,I3,F20.10,F17.10,4F17.10)

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

 9610 FORMAT(/1X,' REGULAR ELECTRONIC ORBITALS AND EIGENVALUES:       ')

 9620 FORMAT(/1X,' SPECIAL ELECTRONIC ORBITALS AND EIGENVALUES:       ')

 9700 FORMAT(/1X,'      QM PARTICLE ORBITALS AND EIGENVALUES:         ')

 9800 FORMAT(10X,15(1H-),'START SECOND ORDER SCF',15(1H-))
                                           
!--------OUTPUT-FORMATTING---------------------------------------------)

      write(*,*) "V_ee terms only"
      LOCBSE2=LOCBSE
C      LOCBSE2=.true. 
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
         call RXCHF_read_CAE(nebf,NAE,DAE,VECAE0)
         call RXCHF_read_CBE(nebf,NBE,DBE,VECBE0)
      else
!        STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
         write(*,*)'ABOUT TO CALL guess_A_elec'
!        call guess_elec(nelec,nebf,xxse,GAM_ecore,DE)
         if ((LOCBSE).or.(LOCBSE2)) then
          call RXCHF_guess_elec(nae,nbe,nebf,xxse,GAM_ecore,
     x                          DAE,DBE,VECAE0,VECBE0)
          write(*,*)'BACK FROM guess_elec for OCBSE'
         else
          call RXCHF_guess_A_elec(NAE,nebf,xxse,GAM_ecore,DAE,VECAE0)
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
      if((nae.eq.1).or.LOCBSE) then
         LSOSCF=.FALSE.
      end if
      if(LSOSCF) THEN
         NFT15=15
         OPEN(NFT15, FILE='WORK15', STATUS='UNKNOWN',
     *        ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         SOGTOL=1.0d+00
         SMALL=1.0D-06
         ITSOA=0
         L0=nebf
         L1=nebf
         NA=nae/2
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
C      ORBGRDB=0.0d+00

C ARS( gam2 testing
C      open(unit=98,form="unformatted",file="gam2-icr.ufm")
C      write(98) GM2ICR
C      close(98)
C )
      DO I=1,MAXIT

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------(
!     write(*,*)
!     write(*,*)'IN: xcuscf '
!     write(*,*)'Before call to UHF_FOCK, GAM_EE='
!     write(*,*)GAM_ee
!     write(*,*)
       if((.not.((locbse).or.(locbse2))).or.(i.eq.1)) then

C      write(*,*) LCMF,LNEOHF,nelec,NAE,NBE,
C     x                     nebf,nebf2,npbf,npbf2,ngee,
C     x                     ng1,ng2,ng3,ng4,
C     x                     SZG2ICR,SZG3IC1,SZG4IC,
C     x                     NG2CHK,NG3CHK,NG4CHK
C     x                     DAE,DBE,DP,
C     x                     GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
C     x                     GM2_1ICR,GM2_2ICR,GM2sICR,
C     x                     GM3_1IC1,GM3_2IC1,GM3_3IC1,
C     x                     GM4ICR,
           call RXCHF_FOCK(LCMF,LNEOHF,nelec,NAE,NBE,
     x                     nebf,nebf2,npbf,npbf2,ngee,
     x                     ng1,ng2,ng3,ng4,
     x                     SZG2ICR,SZG3IC1,SZG4IC,
     x                     NG2CHK,NG3CHK,NG4CHK,
     x                     DAE,DBE,DP,
     x                     GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                     GM2_1ICR,GM2_2ICR,GM2sICR,
     x                     GM3_1IC1,GM3_2IC1,
     x                     GM4ICR,
     x                     FP,FAE,FBE,
     x                     E_OMG1,
     x                     E_OMG2,
     x                     E_OMG3,
     x                     E_OMG4,
     x                     S_OMG1,
     x                     S_OMG2,
     x                     E_total,
     x                     S_total)

         E_total=E_total+E_nuc
       end if

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------)
         if(I.eq.1) then
            WRITE(*,9400)
            WRITE(*,9300) E_nuc,E_OMG1,E_OMG2,E_OMG3,E_OMG4,
     x                    S_OMG1,S_OMG2,S_total,E_total
            if(LSOSCF) then 
               WRITE(*,9050)
            else
               WRITE(*,9000)
            end if
         end if
!        Fockp diag
         call UROOTHAN(vecP,EP,xxsp,FP,npbf)
         call construct_DP(nucst,npbf,vecP,DP)

         if (LOCBSE) then

! Do OCBSE procedure
           call RXCHF_OCBSE(nebf,nae,nbe,vecAE0,vecBE0,FAE,FBE,xxse,
     x                       vecAE,vecBE,AEe,BEe)

! Form regular electronic density matrix and store stuff for next it
           call RXCHF_construct_DAE(NAE,nebf,vecAE,DAE)
           CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
           CALL COPYDEN(DAE0,DAE,NEBF)
           CALL COPYDEN(vecAE0,vecAE,NEBF)

! Form special electronic density matrix and store stuff for next it
           call RXCHF_construct_DBE(NBE,nebf,vecBE,DBE)
           CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
           CALL COPYDEN(DBE0,DBE,NEBF)
           CALL COPYDEN(vecBE0,vecBE,NEBF)

! Calculate energy for this it and Fock matrices for next it

             call RXCHF_FOCK(LCMF,LNEOHF,nelec,NAE,NBE,
     x                       nebf,nebf2,npbf,npbf2,ngee,
     x                       ng1,ng2,ng3,ng4,
     x                       SZG2ICR,SZG3IC1,SZG4IC,
     x                       NG2CHK,NG3CHK,NG4CHK,
     x                       DAE,DBE,DP,
     x                       GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                       GM2_1ICR,GM2_2ICR,GM2sICR,
     x                       GM3_1IC1,GM3_2IC1,
     x                       GM4ICR,
     x                       FP,FAE,FBE,
     x                       E_OMG1,
     x                       E_OMG2,
     x                       E_OMG3,
     x                       E_OMG4,
     x                       S_OMG1,
     x                       S_OMG2,
     x                       E_total,
     x                       S_total)

           E_total=E_total+E_nuc

         else if (LOCBSE2) then

! Regular electrons
         if(LSOSCF) THEN
          ITER=I
          EIGAVL = ITER.GT.1
         end if
         IF(LSOSCF .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAE,FLT)
          call SOGRAD(GRADA,FLT,vecAE,WRK,NPRA,NA,L0,L1,NEBFLT,ORBGRDA)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 900  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AEE,NPRA,L0,NA,NA)
              END IF
              ITSOA = ITSOA+1
           call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,DISPLA,UPDTA,
     *                 DISPLNA,DGRADIA,UPDTIA,ORBGRDA,NPRA,ITSOA,NFT15)
            call SOTRAN(DISPLIA,vecAE,GA,WRK,NPRA,L0,L1,NA,NA,ORBGRDA)
             CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              call RXCHF_construct_DAE(NAE,nebf,vecAE,DAE)
              GO TO 950  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------)

  900 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAE,AEE,xxse,FAE,nebf)
         call RXCHF_construct_DAE(NAE,nebf,vecAE,DAE)

  950 CONTINUE
!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
         CALL COPYDEN(DAE0,DAE,NEBF)

! Do OCBSE procedure for special electrons
           call RXCHF_OCBSE2(nebf,nae,nbe,vecAE,vecBE0,FBE,xxse,
C ARS( 07/01 testing
     x           elcam,elcbfc,ampeb2c,agebfcc,elcex,
C )
     x                       vecBE,BEe)

! Form special electronic density matrix and store stuff for next it
           call RXCHF_construct_DBE(NBE,nebf,vecBE,DBE)
           CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
           CALL COPYDEN(DBE0,DBE,NEBF)
           CALL COPYDEN(vecAE0,vecAE,NEBF)
           CALL COPYDEN(vecBE0,vecBE,NEBF)

! Calculate energy for this it and Fock matrices for next it

             call RXCHF_FOCK(LCMF,LNEOHF,nelec,NAE,NBE,
     x                       nebf,nebf2,npbf,npbf2,ngee,
     x                       ng1,ng2,ng3,ng4,
     x                       SZG2ICR,SZG3IC1,SZG4IC,
     x                       NG2CHK,NG3CHK,NG4CHK,
     x                       DAE,DBE,DP,
     x                       GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                       GM2_1ICR,GM2_2ICR,GM2sICR,
     x                       GM3_1IC1,GM3_2IC1,
     x                       GM4ICR,
     x                       FP,FAE,FBE,
     x                       E_OMG1,
     x                       E_OMG2,
     x                       E_OMG3,
     x                       E_OMG4,
     x                       S_OMG1,
     x                       S_OMG2,
     x                       E_total,
     x                       S_total)

           E_total=E_total+E_nuc

         else

!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------(
         if(LSOSCF) THEN
          ITER=I
          EIGAVL = ITER.GT.1
         end if
         IF(LSOSCF .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAE,FLT)
          call SOGRAD(GRADA,FLT,vecAE,WRK,NPRA,NA,L0,L1,NEBFLT,ORBGRDA)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AEE,NPRA,L0,NA,NA)
              END IF
              ITSOA = ITSOA+1
           call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,DISPLA,UPDTA,
     *                 DISPLNA,DGRADIA,UPDTIA,ORBGRDA,NPRA,ITSOA,NFT15)
            call SOTRAN(DISPLIA,vecAE,GA,WRK,NPRA,L0,L1,NA,NA,ORBGRDA)
             CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              call RXCHF_construct_DAE(NAE,nebf,vecAE,DAE)
              GO TO 750  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------)

  700 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAE,AEE,xxse,FAE,nebf)
         call RXCHF_construct_DAE(NAE,nebf,vecAE,DAE)

  750 CONTINUE
!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
         CALL COPYDEN(DAE0,DAE,NEBF)

!-----------------------POSSIBLE-SOSCF-BETA----------------------------(
!        if(LSOSCF) THEN
!        ITER=I
!        EIGAVL = ITER.GT.1
C         IF(LSOSCF .AND.  EIGAVL) THEN
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
C           call pack_LT(nebf,nebfLT,FBE,FLT)
C          call SOGRAD(GRADB,FLT,vecBE,WRK,NPRB,NBE,L0,L1,NEBFLT,ORBGRDB)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
C            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
C              IF(ITSOB.EQ.0) THEN
!             WRITE(*,9800)
C                 call SOHESS(HSTARTB,BEE,NPRB,L0,NBE,NBE)
C              END IF
C              ITSOB = ITSOB+1
C           call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,DISPLB,UPDTB,
C     *                 DISPLNB,DGRADIB,UPDTIB,ORBGRDB,NPRB,ITSOB,NFT16)
C            call SOTRAN(DISPLIB,vecBE,GB,WRK,NPRB,L0,L1,NBE,NBE,ORBGRDB)
C             CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
C              call construct_DAE(NBE,nebf,vecBE,DBE)
C              GO TO 850  ! Use the new C's to form new density (change)
C            END IF
C         END IF
!-----------------------POSSIBLE-SOSCF-BETA----------------------------)

  800 CONTINUE
!        call ROOTHAN(DBE,vecBE,BEE,xxse,FBE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecBE,BEE,xxse,FBE,nebf)
         call RXCHF_construct_DBE(NBE,nebf,vecBE,DBE)

  850 CONTINUE
!        --> FIND LARGEST CHANGE IN Beta E DENSITY
         CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
         CALL COPYDEN(DBE0,DBE,NEBF)

         end if ! end if for not ocbse or ocbse2

!        --> FIND LARGEST CHANGE IN P DENSITY
         CALL DENDIF(DP0,DP,NPBF,DIFFP)
         CALL COPYDEN(DP0,DP,NPBF)

!        --> CALCULATE CHANGE IN TOTAL ENERGY
         Delta_E_tot=E_total-E_total_old
         E_total_old=E_total

!        --> PRINT SUMMARY OF THIS ITERATION
         if(LSOSCF) then
            WRITE(*,9150) I,E_total,Delta_E_tot,DIFFAE,DIFFBE,DIFFP,
     x                    ORBGRDA
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
      end if

!     PRINT FINAL ENERGY AND PUNCH THE ORBITALS
      WRITE(*,9200) E_TOTAL,I
!
        WRITE(*,9300) E_nuc,E_OMG1,E_OMG2,E_OMG3,E_OMG4,
     x                S_OMG1,S_OMG2,S_total,E_total

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
      SUBROUTINE RXCHF_guess_A_elec(NAE,nebf,xxse,GAM_ecore,DAE,CAE)
 
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

      call RXCHF_construct_DAE(NAE,nebf,CAE,DAE)


      RETURN
      END
!======================================================================
      SUBROUTINE RXCHF_guess_B_elec(NBE,nebf,xxse,GAM_ecore,DBE,CBE)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial special electron guess density
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer NBE
      double precision xxse(nebf,nebf)
      double precision GAM_ecore(nebf,nebf)
! Variables Returned
      double precision DBE(nebf,nebf)
! Local variables
      double precision CBE(nebf,nebf)
      double precision EVF(nebf)


      call UROOTHAN(CBE,EVF,xxse,GAM_ecore,nebf)

      call RXCHF_construct_DBE(NBE,nebf,CBE,DBE)


      RETURN
      END
!======================================================================
      subroutine RXCHF_guess_elec(nae,nbe,nebf,xxse,GAM_ecore,
     x                            DAE,DBE,CAE,CBE)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial regular and special electronic guess density
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nae
      integer nbe
      double precision xxse(nebf,nebf)
      double precision GAM_ecore(nebf,nebf)
! Variables Returned
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision CAE(nebf,nebf)
      double precision CBE(nebf,nebf)
! Local variables
      integer i,j
      integer nocca,noccb
      double precision C(nebf,nebf)
      double precision EVF(nebf)
      double precision zero
      parameter(zero=0.0d+00)

      if (nae.gt.1) then
       nocca=nae/2
      else
       nocca=nae
      end if
      noccb=1

      DAE=zero
      DBE=zero
      CAE=zero
      CBE=zero

      call UROOTHAN(C,EVF,xxse,GAM_ecore,nebf)

! Store first nocca evectors as occ reg elec vectors
      do i=1,nocca
        do j=1,nebf
          CAE(j,i)=C(j,i)
        end do
      end do

! Store nebf-nocca remaining evectors as occ and virt spec elec vectors
      do i=nocca+1,nebf
        do j=1,nebf
          CBE(j,i-nocca)=C(j,i)
        end do
      end do

      call RXCHF_construct_DAE(NAE,nebf,CAE,DAE)
      call RXCHF_construct_DBE(NBE,nebf,CBE,DBE)

      return
      end
!======================================================================
      subroutine RXCHF_OCBSE(nebf,nae,nbe,vecAE0,vecBE0,FAE,FBE,Selec,
     x                       vecAE,vecBE,AEen,BEen)
! 
!     Perform OCBSE procedure
!       vecAE0: Regular electron coefficients from previous iteration
!       vecBE0: Special electron coefficients from previous iteration
!       vecAE:  Regular electron coefficients from current iteration
!       vecBE:  Special electron coefficients from current iteration
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nae,nbe
      double precision vecAE0(nebf,nebf),vecBE0(nebf,nebf)
      double precision FAE(nebf,nebf),FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
! Variables Returned
      double precision vecAE(nebf,nebf),vecBE(nebf,nebf)
      double precision AEen(nebf),BEen(nebf)
! Local variables
      integer nocca,noccb,nvirt 
      integer noccvirta,noccvirtb
      double precision zero
      parameter(zero=0.0d+00)

      if (nae.gt.1) then
       nocca=nae/2
      else
       nocca=nae
      end if
      noccb=1
      nvirt=nebf-nocca-noccb
      noccvirta=nebf-noccb
      noccvirtb=nebf-nocca

      vecAE=zero
      vecBE=zero
      AEen=zero
      BEen=zero

      call RXCHF_OCBSE_driver(nebf,nae,nbe,nocca,noccb,
     x                        nvirt,noccvirta,noccvirtb,
     x                        vecAE0,vecBE0,FAE,FBE,Selec,
     x                        vecAE,vecBE,AEen,BEen)

      return
      end
!======================================================================
      subroutine RXCHF_OCBSE_driver(nebf,nae,nbe,nocca,noccb,
     x                              nvirt,noccvirta,noccvirtb,
     x                              vecAE0,vecBE0,FAE,FBE,Selec,
     x                              vecAE,vecBE,AEen,BEen)
!
!     OCBSE Procedure:
!       - Construct transformation matrix for regular electrons
!       - Transform FAE and diagonalize to obtain vecAE
!       - Construct transformation matrix for special electrons
!       - Transform FBE and diagonalize to obtain vecBE
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nae,nbe
      double precision vecAE0(nebf,nebf),vecBE0(nebf,nebf)
      double precision FAE(nebf,nebf),FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
! Variables Returned
      double precision vecAE(nebf,nebf),vecBE(nebf,nebf)
      double precision AEen(nebf),BEen(nebf)
! Local variables
      integer i,j
      integer nocca,noccb,nvirt
      integer noccvirta,noccvirtb
      double precision WA(nebf,noccvirta),WB(nebf,noccvirtb)
      double precision WAtrans(noccvirta,nebf),WBtrans(noccvirtb,nebf)
      double precision wFAEw(noccvirta,noccvirta)
      double precision wFBEw(noccvirtb,noccvirtb)
      double precision wSAw(noccvirta,noccvirta)
      double precision wSBw(noccvirtb,noccvirtb)
      double precision AUXA(nebf,noccvirta)
      double precision AUXB(nebf,noccvirtb)
      double precision xAEen(noccvirta)
      double precision xBEen(noccvirtb)
      double precision xvecAE(noccvirta,noccvirta)
      double precision xvecBE(noccvirtb,noccvirtb)
      double precision blockvecAE(nebf,noccvirta)
      double precision blockvecBE(nebf,noccvirtb)
      double precision zero
      parameter(zero=0.0d+00)

      logical debug
      debug=.false.

! Initialize
      WA=zero
      WB=zero
      WAtrans=zero
      WBtrans=zero
      wFAEw=zero
      wFBEw=zero
      wSAw=zero
      wSBw=zero
      AUXA=zero
      AUXB=zero
      xAEen=zero
      xBEen=zero
      xvecAE=zero
      xvecBE=zero
      blockvecAE=zero
      blockvecBE=zero

C ARS(
      if (debug) then
      write(*,*) "MATRIX Previous vecAE:"
      call PREVNU(vecAE0,AEen,nebf,nebf,nebf)
      write(*,*) "MATRIX Previous vecBE:"
      call PREVNU(vecBE0,BEen,nebf,nebf,nebf)
      write(*,*) "MATRIX Previous S:"
      call PREVNU(Selec,BEen,nebf,nebf,nebf)
      end if
C )

!!!!!!!!!!!!!!!! Solve for regular electronic solution !!!!!!!!!!!!!!!!

! Form regular electronic transformation matrix
      do i=1,nvirt
        do j=1,nebf
          WA(j,i)=vecBE0(j,noccb+i)
        end do
      end do

      do i=1,nocca
        do j=1,nebf
          WA(j,i+nvirt)=vecAE0(j,i)
        end do
      end do

C ARS(
      if (debug) then
      write(*,*) "MATRIX WA:"
      call PREVNU(WA,xAEen,noccvirta,nebf,nebf)
      end if
C )

! Form transpose
      do i=1,noccvirta
        do j=1,nebf
          WAtrans(i,j)=WA(j,i)
        end do
      end do

! Transform FAE as Wtrans * FAE * W
      call RXCHF_matmult(nebf,nebf,nebf,noccvirta,
     x                   FAE,WA,AUXA)
      call RXCHF_matmult(noccvirta,nebf,nebf,noccvirta,
     x                   WAtrans,AUXA,wFAEw)

! Transform AO overlap matrix as Wtrans * S * W
      call RXCHF_matmult(nebf,nebf,nebf,noccvirta,
     x                   Selec,WA,AUXA)
      call RXCHF_matmult(noccvirta,nebf,nebf,noccvirta,
     x                   WAtrans,AUXA,wSAw)

C ARS(
      if (debug) then
      write(*,*) "MATRIX wSAw:"
      call PREVNU(wSAw,xAEen,noccvirta,noccvirta,noccvirta)
      end if
C )

! Diagonalize transformed FAE to obtain solutions in new basis
      call UROOTHAN(xvecAE,xAEen,wSAw,wFAEw,noccvirta)

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX xvecAE:"
      call PREVNU(xvecAE,xAEen,noccvirta,noccvirta,noccvirta)
      end if
C )

! Transform evectors to AO basis as vecAE = W * xvecAE
      call RXCHF_matmult(nebf,noccvirta,noccvirta,noccvirta,
     x                   WA,xvecAE,blockvecAE)

! Pass evectors to output variables (zeros for unfilled noccb part)
      do i=1,noccvirta
        AEen(i)=xAEen(i)
        do j=1,nebf
          vecAE(j,i)=blockvecAE(j,i)
        end do
      end do

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX New vecAE:"
      call PREVNU(vecAE,AEen,nebf,nebf,nebf)
      end if
C )

!!!!!!!!!!!!!!!! Solve for special electronic solution !!!!!!!!!!!!!!!!

! Form regular electronic transformation matrix
      do i=1,nvirt
        do j=1,nebf
          WB(j,i)=vecAE(j,nocca+i)
        end do
      end do

      do i=1,noccb
        do j=1,nebf
          WB(j,i+nvirt)=vecBE0(j,i)
        end do
      end do

C ARS(
      if (debug) then
      write(*,*) "MATRIX WB:"
      call PREVNU(WB,xBEen,noccvirtb,nebf,nebf)
      end if
C )

! Form transpose
      do i=1,noccvirtb
        do j=1,nebf
          WBtrans(i,j)=WB(j,i)
        end do
      end do

! Transform FBE as Wtrans * FBE * W
      call RXCHF_matmult(nebf,nebf,nebf,noccvirtb,
     x                   FBE,WB,AUXB)
      call RXCHF_matmult(noccvirtb,nebf,nebf,noccvirtb,
     x                   WBtrans,AUXB,wFBEw)

! Transform AO overlap matrix as Wtrans * S * W
      call RXCHF_matmult(nebf,nebf,nebf,noccvirtb,
     x                   Selec,WB,AUXB)
      call RXCHF_matmult(noccvirtb,nebf,nebf,noccvirtb,
     x                   WBtrans,AUXB,wSBw)

C ARS(
      if (debug) then
      write(*,*) "MATRIX wSBw:"
      call PREVNU(wSBw,xBEen,noccvirtb,noccvirtb,noccvirtb)
      end if
C )

! Diagonalize transformed FBE to obtain solutions in new basis
      call UROOTHAN(xvecBE,xBEen,wSBw,wFBEw,noccvirtb)

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX xvecBE:"
      call PREVNU(xvecBE,xBEen,noccvirtb,noccvirtb,noccvirtb)
      end if
C )

! Transform evectors to AO basis as vecBE = W * xvecBE
      call RXCHF_matmult(nebf,noccvirtb,noccvirtb,noccvirtb,
     x                   WB,xvecBE,blockvecBE)

! Pass evectors to output variables (zeros for unfilled noccb part)
      do i=1,noccvirtb
        BEen(i)=xBEen(i)
        do j=1,nebf
          vecBE(j,i)=blockvecBE(j,i)
        end do
      end do

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX New vecBE:"
      call PREVNU(vecBE,BEen,nebf,nebf,nebf)
      end if
C )

      return
      end
!======================================================================
      subroutine RXCHF_OCBSE2(nebf,nae,nbe,vecAE,vecBE0,FBE,Selec,
C ARS( 07/01 testing
     x           elcam,elcbfc,ampeb2c,agebfcc,elcex,
C )
     x                        vecBE,BEen)
! 
!     Perform OCBSE procedure for special electrons only
!       vecBE0: Special electron coefficients from previous iteration
!       vecAE:  Regular electron coefficients from current iteration
!       vecBE:  Special electron coefficients from current iteration
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nae,nbe
      double precision vecAE(nebf,nebf),vecBE0(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
! Variables Returned
      double precision vecBE(nebf,nebf)
      double precision BEen(nebf)
! Local variables
      integer nocca,noccb,nvirt 
      integer noccvirta,noccvirtb
      double precision zero
      parameter(zero=0.0d+00)
C ARS( 07/01 testing
      integer mo1,mo2,ie1,je1,iec1,jec1
      integer i1,j1,k1,l1,m1,n1,a1,b1
      double precision ans,ovlap,cof_ie1,cof_je1
      double precision Amat1(3), Bmat1(3)
      integer ELCAM(nebf,3)  ! Angular mom for electrons
      double precision ELCBFC(nebf,3) ! Basis centers: elec basis
      integer AMPEB2C(nebf) ! Map primitive index to contracted
      double precision AGEBFCC(nebf) ! Map prim index to contract coef
      double precision ELCEX(nebf) ! Exponents: elec basis
C )

      if (nae.gt.1) then
       nocca=nae/2
      else
       nocca=nae
      end if
      noccb=1
      nvirt=nebf-nocca-noccb
      noccvirta=nebf-noccb
      noccvirtb=nebf-nocca

      vecBE=zero
      BEen=zero

      call RXCHF_OCBSE_driver2(nebf,nae,nbe,nocca,noccb,
     x                         nvirt,noccvirta,noccvirtb,
     x                         vecAE,vecBE0,FBE,Selec,
C ARS( 07/01 testing
     x           elcam,elcbfc,ampeb2c,agebfcc,elcex,
C )
     x                         vecBE,BEen)

      return
      end
!======================================================================
      subroutine RXCHF_OCBSE_driver2(nebf,nae,nbe,nocca,noccb,
     x                               nvirt,noccvirta,noccvirtb,
     x                               vecAE,vecBE0,FBE,Selec,
C ARS( 07/01 testing
     x           elcam,elcbfc,ampeb2c,agebfcc,elcex,
C )
     x                               vecBE,BEen)
!
!     OCBSE Procedure:
!       - Construct transformation matrix for special electrons
!       - Transform FBE and diagonalize to obtain vecBE
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nae,nbe
      double precision vecAE(nebf,nebf),vecBE0(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
! Variables Returned
      double precision vecBE(nebf,nebf)
      double precision BEen(nebf)
! Local variables
      integer i,j
      integer nocca,noccb,nvirt
      integer noccvirta,noccvirtb
      double precision WB(nebf,noccvirtb)
      double precision WBtrans(noccvirtb,nebf)
      double precision wFBEw(noccvirtb,noccvirtb)
      double precision wSBw(noccvirtb,noccvirtb)
      double precision AUXB(nebf,noccvirtb)
      double precision xAEen(noccvirta)
      double precision xBEen(noccvirtb)
      double precision xvecBE(noccvirtb,noccvirtb)
      double precision blockvecBE(nebf,noccvirtb)
      double precision zero
      parameter(zero=0.0d+00)

C ARS( 07/01 testing
      integer mo1,mo2,ie1,je1,iec1,jec1,motodiscard
      integer i1,j1,k1,l1,m1,n1
      double precision ans,ovlap,cof_ie1,cof_je1,a1,b1,maxovlap
      double precision Amat1(3), Bmat1(3)
      integer ELCAM(nebf,3)  ! Angular mom for electrons
      double precision ELCBFC(nebf,3) ! Basis centers: elec basis
      integer AMPEB2C(nebf) ! Map primitive index to contracted
      double precision AGEBFCC(nebf) ! Map prim index to contract coef
      double precision ELCEX(nebf) ! Exponents: elec basis
      double precision ovlaparr(nebf),projorb(nebf) ! fix2
C )
      logical debug
      debug=.false.

! Initialize
      WB=zero
      WBtrans=zero
      wFBEw=zero
      wSBw=zero
      AUXB=zero
      xBEen=zero
      xvecBE=zero
      blockvecBE=zero

C ARS(
      if (debug) then
      write(*,*) "MATRIX vecAE:"
      call PREVNU(vecAE,BEen,nebf,nebf,nebf)
      write(*,*) "MATRIX Previous vecBE:"
      call PREVNU(vecBE0,BEen,nebf,nebf,nebf)
      write(*,*) "MATRIX Previous S:"
      call PREVNU(Selec,BEen,nebf,nebf,nebf)
      end if
C )

C ARS( 07/01 testing
      maxovlap=0.0d0
      do mo1=1,nebf
      mo2=1

        ovlap=0.0d0

        do ie1=1,nebf
        do je1=1,nebf

          iec1=AMPEB2C(ie1)
          jec1=AMPEB2C(je1)
          Cof_ie1=AGEBFCC(ie1)
          Cof_je1=AGEBFCC(je1)

          A1=ELCEX(ie1)
          I1=ELCAM(ie1,1)
          J1=ELCAM(ie1,2)
          K1=ELCAM(ie1,3)
          Amat1(1)=ELCBFC(ie1,1)
          Amat1(2)=ELCBFC(ie1,2)
          Amat1(3)=ELCBFC(ie1,3)

          B1=ELCEX(je1)
          L1=ELCAM(je1,1)
          M1=ELCAM(je1,2)
          N1=ELCAM(je1,3)
          Bmat1(1)=ELCBFC(je1,1)
          Bmat1(2)=ELCBFC(je1,2)
          Bmat1(3)=ELCBFC(je1,3)

          call gfovlap(I1,J1,K1,A1,Amat1,
     2                 L1,M1,N1,B1,Bmat1,
     3                 ans)

          ovlap=ovlap+cof_ie1*cof_je1*vecAE(ie1,mo1)*vecBE0(je1,mo2)*ans

        end do
        end do

        if (debug) then
         write(*,*) "Overlap (mo1,mo2,ovlap):"
         write(*,*) mo1,mo2,ovlap
        end if

        ovlaparr(mo1)=ovlap
        if ((mo1.gt.1).and.(dabs(ovlap).gt.maxovlap)) then
         maxovlap=dabs(ovlap)
         motodiscard=mo1
        end if

      end do

      if (debug) then
       write(*,*) "Discarding MO: ",motodiscard," which had overlap",
     x            maxovlap
      end if
C )

!!!!!!!!!!!!!!!! Solve for special electronic solution !!!!!!!!!!!!!!!!

! Form special electronic transformation matrix
      do i=1,nvirt
        do j=1,nebf
C ARS( 07/01 testing
C          WB(j,i)=vecAE(j,nocca+i)
          if (i.ne.(motodiscard-nocca)) then
           WB(j,i)=vecAE(j,nocca+i)
          else
           WB(j,i)=vecAE(j,nebf)
          end if
C )
        end do
      end do

C ARS( 07/01 testing
C FIX2
      do j=1,nebf
        projorb(j)=vecBE0(j,1)
        do i=1,nebf
          if (i.ne.motodiscard) then
           projorb(j)=projorb(j)-ovlaparr(i)*vecAE(j,i)
          end if
        end do
      end do

      ovlap=0.0d0
      do i=1,nebf
      do j=1,nebf
        ovlap=ovlap+projorb(i)*projorb(j)*Selec(i,j)
      end do
      end do
C )

      do i=1,noccb
        do j=1,nebf
C ARS( 07/01 testing
C          WB(j,i+nvirt)=vecBE0(j,i)
C
C *** ONLY WORKS FOR noccb=1 ***
C FIX1
C
C          WB(j,i+nvirt)=vecAE(j,motodiscard)
C
C FIX2
C
          WB(j,i+nvirt)=projorb(j)/dsqrt(ovlap)
C )
        end do
      end do

C ARS(
      if (debug) then
      write(*,*) "MATRIX WB:"
      call PREVNU(WB,xBEen,noccvirtb,nebf,nebf)
      end if
C )

! Form transpose
      do i=1,noccvirtb
        do j=1,nebf
          WBtrans(i,j)=WB(j,i)
        end do
      end do

! Transform FBE as Wtrans * FBE * W
      call RXCHF_matmult(nebf,nebf,nebf,noccvirtb,
     x                   FBE,WB,AUXB)
      call RXCHF_matmult(noccvirtb,nebf,nebf,noccvirtb,
     x                   WBtrans,AUXB,wFBEw)

! Transform AO overlap matrix as Wtrans * S * W
      call RXCHF_matmult(nebf,nebf,nebf,noccvirtb,
     x                   Selec,WB,AUXB)
      call RXCHF_matmult(noccvirtb,nebf,nebf,noccvirtb,
     x                   WBtrans,AUXB,wSBw)

C ARS(
      if (debug) then
      write(*,*) "MATRIX wSBw:"
      call PREVNU(wSBw,xBEen,noccvirtb,noccvirtb,noccvirtb)
      end if
C )

! Diagonalize transformed FBE to obtain solutions in new basis
      call UROOTHAN(xvecBE,xBEen,wSBw,wFBEw,noccvirtb)

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX xvecBE:"
      call PREVNU(xvecBE,xBEen,noccvirtb,noccvirtb,noccvirtb)
      end if
C )

! Transform evectors to AO basis as vecBE = W * xvecBE
      call RXCHF_matmult(nebf,noccvirtb,noccvirtb,noccvirtb,
     x                   WB,xvecBE,blockvecBE)

! Pass evectors to output variables (zeros for unfilled noccb part)
      do i=1,noccvirtb
        BEen(i)=xBEen(i)
        do j=1,nebf
          vecBE(j,i)=blockvecBE(j,i)
        end do
      end do

C ARS(
      if (debug) then
      WRITE(*,*) "MATRIX New vecBE:"
      call PREVNU(vecBE,BEen,nebf,nebf,nebf)
      end if
C )

      return
      end
!======================================================================
      subroutine RXCHF_matmult(rowsa,colsa,rowsb,colsb,A,B,C)
!
!     Computes C = A * B
!======================================================================
      implicit none
! Input Variables
      integer rowsa,colsa,rowsb,colsb
      double precision A(rowsa,colsa)
      double precision B(rowsb,colsb)
! Output Variables
      double precision C(rowsa,colsb)
! Local Variables
      integer i,j,k
      double precision zero
      parameter(zero=0.0d+00)

      C=zero

      if (colsa.ne.rowsb) then
       write(*,*) "ERROR IN RXCHF_MATMULT, QUITTING"
       STOP
      end if

      do i=1,colsb
        do j=1,rowsa
          do k=1,colsa
            C(j,i)=C(j,i)+A(j,k)*B(k,i)
          end do
        end do
      end do

      return
      end
