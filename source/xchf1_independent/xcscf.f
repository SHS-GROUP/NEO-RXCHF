!======================================================================
      subroutine xcscf(nelec,NPR,NEBFLT,NUCST,
     x                 npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                 ngtg1,ng1,ng2,ng3,ng4,NG4CHK,NG3CHK,NG2CHK,
     x                 read_CE,read_CP,
     x                 LNEOHF,LGAM4,LG4DSCF,LG3DSCF,LG2DSCF,LCMF,
     x                 ng2prm,ng3prm,nat,pmass,cat,zan,bcoef1,gamma1,
     x                 KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                 ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                 LG2IC1,SZG2ICR,GM2ICR,GM2SICR,
     x                 LG3IC1,SZG3IC1,GM3IC1,
     x                 LG4IC,SZG4IC,GM4ICR)
!    x                 LG4IC,SZG4IC,GM4ICR,LNONAD)
!
!     --- PERFORM A NUCLEAR-ELECTRONIC XC HARTREE-FOCK
!     CALCULATION FOR A N-ELECTRON ONE-PROTON SYSTEM ---
!
!     **DEFINITIONS:
!
!     *FOR ELECTRONS:
!     DE    ::  NEW ELECTRON DENSITY MATRIX
!     DE0   ::  OLD ELECTRON DENSITY MATRIX
!     VECE  ::  ELECTRON MOS
!     EE    ::  ELECTRON ORBITAL EIGENVALUES
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
!     logical LNONAD
      logical LNEOHF,LG4DSCF,read_CE,read_CP
      logical LGAM4
      logical LCMF
      logical LG3DSCF
      logical LG2IC1
      logical LG3IC1
      logical LG4IC
      logical LG2DSCF
      integer NG4CHK
      integer SZG2ICR
      integer SZG3IC1
      integer SZG4IC
      integer NG3CHK
      integer NG2CHK
      integer nelec
      integer NPR
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
!     double precision zan(nat)
!     double precision cat(3,nat)
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

! Local variables
      double precision zero,one
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00) 

      double precision xxse(nebf,nebf)  ! Elec overlap matrix
      double precision xxsp(npbf,npbf)  ! Nuc overlap matrix
      double precision GAM_pcore(npbf2)
      double precision GAM_ecore(nebf2)
      double precision GAM_ep(ng1)
      double precision GAM_ee(ngee)

      integer noccE ! Number of occupied elec orbs
      integer noccP ! Number of occupied nuc orbs
      integer maxit

      integer i
      integer j
      integer k
      integer l

      integer NPBFLT

      double precision TOLP
      double precision TOLE
      double precision diffE
      double precision diffP

      double precision E_ecore
      double precision E_pcore
      double precision E_ee
      double precision E_ep
      double precision E_nuc
      double precision E_total

      double precision E_gam1
      double precision E_gam2
      double precision E_gam3
      double precision E_gam4

      double precision S_gam1s
      double precision S_total
      double precision S_gam0
      double precision S_gam1
      double precision S_gam2

      double precision DE(NEBF,NEBF)
      double precision DE0(NEBF,NEBF)
      double precision VECE(NEBF,NEBF)
      double precision EE(NEBF)
      double precision focke(nebf,nebf)
      double precision xfocke(nebf,nebf)
      double precision feDUMMY(nebf,nebf)
      double precision xeDUMMY(nebf,nebf)

      double precision DP(NPBF,NPBF)
      double precision DP0(NPBF,NPBF)
      double precision VECP(NPBF,NPBF)
      double precision EP(NPBF)
      double precision fockp(npbf,npbf)
      double precision xfockp(npbf,npbf)
      double precision fpDUMMY(npbf,npbf)
      double precision xpDUMMY(npbf,npbf)

      double precision E_total_old
      double precision Delta_E_tot

!--------SOSCF-RELATED-VARIABLES------------(
      logical SOSCF
      logical EIGAVL
      integer NA
      integer ITER
      integer ITSO ! SOSCF iteration counter
      integer L0,L1
      integer NFT15
      double precision FLT(NEBFLT) !FLT: Lower triangle focke
!cc   double precision EE(nebf) ! EIG(L1)
      double precision HSTART(NPR)
      double precision GRAD(NPR)
      double precision PGRAD(NPR)
      double precision DISPLI(NPR)
      double precision DGRAD(NPR)  ! WRK1
      double precision DISPL(NPR)  ! WRK2
      double precision UPDT(NPR)   ! WRK3
      double precision DISPLN(NPR) ! WRK1+NPR
      double precision DGRADI(NPR) ! WRK2+NPR
      double precision UPDTI(NPR)  ! WRK3+NPR
      double precision ORBGRD
      double precision SMALL
      double precision SOGTOL ! ORBGRAD TOL to activate soscf
      double precision X(NPR)
!cc   double precision vecE(nebf,nebf) !C(L1,L0)
      double precision G(nebf,nebf) !G(L0,L0)
      double precision WRK(nebf) !WRK(L0)
!cc   double precision CCC(nebf,nebf) !WRK(L0)
!cc   NPR=(L0-NA)*NA ! Line 2134 RHFCL ?NA is NUM ALPHA E?
!--------SOSCF-RELATED-VARIABLES------------)

!--------OUTPUT-FORMATTING---------------------------------------------(
 9000 FORMAT(/' ITER     TOTAL ENERGY        E CHANGE     E DENSITY ',
     *        'CHANGE   P DENSITY CHANGE    ORB GRAD')

 9100 FORMAT(1X,I3,F20.10,F17.10,3F17.10)

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
     x       12X,'   E_GAM1=',1X,F20.10/
     x       12X,'   E_GAM2=',1X,F20.10/
     x       12X,'   E_GAM3=',1X,F20.10/
     x       12X,'   E_GAM4=',1X,F20.10/
     x       12X,'  S_TOTAL=',1X,F20.10/
     x       12X,'  E_TOTAL=',1X,F20.10/)

 9400 FORMAT(/1X,'          INITIAL GUESS ENERGETICS:            ')

 9500 FORMAT(/6X,' ** BEGIN NEOXC SELF-CONSISTENT-FIELD CALCULATION **')

 9600 FORMAT(/1X,'       ELECTRONIC ORBITALS AND EIGENVALUES:         ')

 9700 FORMAT(/1X,'         NUCLEAR ORBITALS AND EIGENVALUES:          ')

 9800 FORMAT(10X,15(1H-),'START SECOND ORDER SCF',15(1H-))
                                           
!--------OUTPUT-FORMATTING---------------------------------------------)

!-------------DETERMINE-OCCUPATION-NUMBERS-----------------------------(
!  Either electronic closed-shell or 1e/1p case...
      if(nelec.gt.1) then
         noccE=nelec/2
      else
         noccE=1
      endif
    
      noccP=1
!-------------DETERMINE-OCCUPATION-NUMBERS-----------------------------)

!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------(
! Currently E_NUC_REP is calculated when molecular integrals are done.
! This saves having to pass along nat,zan,cat to this routine.
!     call class_nuc_rep(nat,zan,cat,E_nuc)
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
!--------------READ-INTEGRALS-NEEDED-FOR-NEO-HF------------------------)

!-------------INITIAL-GUESSES------------------------------------------(
      if(read_CE) then
!        READ IN GUESS FOR E:
         call read_elec_density(nebf,nelec,DE)
      else
!        STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
         write(*,*)'ABOUT TO CALL guess_elec'
         call guess_elec(nelec,nebf,xxse,GAM_ecore,DE)
         write(*,*)'BACK FROM guess_elec'
      end if
      if(read_CP) then
!        READ IN GUESS FOR N:
         call read_nuc_density(npbf,1,NUCST,DP)
      else
!        STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
         write(*,*)'ABOUT TO CALL guess_prot'
         call guess_prot(NUCST,npbf,nebf,xxsp,GAM_pcore,GAM_ep,DE,DP)
         write(*,*)'BACK FROM guess_prot'
         write(*,*)
      end if
!     BOZO GUESS:
!     DO I=1,NEBF
!        DO J=1,NEBF
!           DE(I,J)=ZERO
!        END DO
!        DE(i,i)=one
!     END DO
!     DO K=1,NPBF
!        DO L=1,NPBF
!           DP(K,L)=ZERO
!        END DO
!        DP(k,k)=one
!     END DO
!-FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-FOR-INITIAL-GUESS-(
!        call xchf1_fock(LNEOHF,nebf,nebf2,npbf,npbf2,nelec,
!    x                   ngee,ng1,ng2,ng3,ng4,DE,DP,
!    x                   GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
!    x                   focke,fockp,Xfocke,Xfockp,
!    x                   E_total,E_nuc,E_ecore,E_pcore,E_ep,
!    x                   E_ee,E_gam1,E_gam2,E_gam3,E_gam4,
!    x                   S_total,S_gam1,S_gam2)
!     WRITE(*,9400)
!     WRITE(*,9300) E_nuc,E_ecore,E_pcore,E_ep,E_ee,
!    xE_gam1,E_gam2,E_gam3,E_gam4,S_total,E_total
!-FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-FOR-INITIAL-GUESS-)

!-------------INITIAL-GUESSES------------------------------------------)

!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------(
!c    SOSCF=.FALSE.
      SOSCF=.TRUE.
      if(nelec.eq.1) then
         SOSCF=.FALSE.
      end if

      if(SOSCF) THEN
         NFT15=15
         OPEN(NFT15, FILE='WORK15', STATUS='UNKNOWN',
     *        ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         SOGTOL=1.0d+00
         SMALL=1.0D-06
         ITSO=0
         L0=nebf
         L1=nebf
         NA=nelec/2
!c       write(*,*)'NPR=',npr
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
!     call zero_out(nebf,DE0)
!     call zero_out(npbf,DP0)
      DE0=0.0d+00
      DP0=0.0d+00
!
!     BEGIN XCSCF ITERATIONS
      WRITE(*,9500)
!     WRITE(*,9000)
!
      E_total_old=0.0d+00
      ORBGRD=0.0d+00

      DO I=1,MAXIT

!-NUC----------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------(
         call xchf1_fock(LNEOHF,LGAM4,LG4DSCF,LG4IC,
     x                   LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                   NG4CHK,NG3CHK,NG2CHK,
     x                   SZG4IC,SZG3IC1,SZG2ICR,
     x                   npebf,nebf,nebf2,npbf,npbf2,nelec,
     x                   ngee,ng1,ng2,ng3,ng4,DE,DP,
     x                   GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                   GM4ICR,GM3IC1,GM2ICR,GM2SICR,
     x                   ng2prm,ngtg1,ng3prm,
     x                   nat,pmass,cat,zan,bcoef1,gamma1,
     x                   KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                   ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                   feDUMMY,fockp,XeDUMMY,Xfockp,
     x                   E_total,E_nuc,E_ecore,E_pcore,E_ep,
     x                   E_ee,E_gam1,E_gam2,E_gam3,E_gam4,
     x                   S_total,S_gam1,S_gam2)

!     subroutine xchf1_fock(LNEOHF,LG4DSCF,LG3DSCF,NG4CHK,NG3CHK,
!    x                      npebf,nebf,nebf2,npbf,npbf2,nelec,
!    x                      ngee,ng1,ng2,ng3,ng4,DE,DP,
!    x                      GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
!    x                      ngtg1,ng3prm,nat,pmass,cat,zan,bcoef1,gamma1,
!    x                      KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
!    x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
!    x                      focke,fockp,Xfocke,Xfockp,
!    x                      E_total,E_nuc,E_ecore,E_pcore,E_ep,
!    x                      E_ee,E_gam1,E_gam2,E_gam3,E_gam4,
!    x                      S_total,S_gam1,S_gam2)

C-NUC----------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------)
         if(I.eq.1) then
            WRITE(*,9400)
            WRITE(*,9300) E_nuc,E_ecore,E_pcore,E_ep,E_ee,
     x      E_gam1,E_gam2,E_gam3,E_gam4,S_total,E_total
            WRITE(*,9000)
         end if
C        Fockp diag
         call ROOTHAN(DP,vecP,EP,xxsp,fockp,npbf,1,2,NUCST)
!CWS-04122011
!     write(*,*)
!     write(*,*)'NUCLEAR DENSITY MATRIX'
!     write(*,*)
!!    write(*,*) DP
!     NPBFLT=npbf*(npbf+1)/2
!     call prt_lower_triangle(npbf,NPBFLT,DP)
!     write(*,*)
!     WRITE(*,9700)
!     call PREVNU(vecp,EP,npbf,npbf,npbf)
!CWS-04122011

C-ELEC---------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------(
!        call xchf1_fock(LNEOHF,LG4DSCF,LG3DSCF,LG2DSCF,
!    x                   NG4CHK,NG3CHK,NG2CHK,
!    x                   npebf,nebf,nebf2,npbf,npbf2,nelec,
!    x                   ngee,ng1,ng2,ng3,ng4,DE,DP,
!    x                   GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,ng2prm,
!    x                     ngtg1,ng3prm,nat,pmass,cat,zan,bcoef1,gamma1,
!    x                     KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
!    x                     ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
!    x                   focke,fpDUMMY,Xfocke,XpDUMMY,
!    x                   E_total,E_nuc,E_ecore,E_pcore,E_ep,
!    x                   E_ee,E_gam1,E_gam2,E_gam3,E_gam4,
!    x                   S_total,S_gam1,S_gam2)

         call xchf1_fock(LNEOHF,LGAM4,LG4DSCF,LG4IC,
     x                   LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                   NG4CHK,NG3CHK,NG2CHK,
     x                   SZG4IC,SZG3IC1,SZG2ICR,
     x                   npebf,nebf,nebf2,npbf,npbf2,nelec,
     x                   ngee,ng1,ng2,ng3,ng4,DE,DP,
     x                   GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                   GM4ICR,GM3IC1,GM2ICR,GM2SICR,
     x                   ng2prm,ngtg1,ng3prm,
     x                   nat,pmass,cat,zan,bcoef1,gamma1,
     x                   KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                   ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                   focke,fpDUMMY,Xfocke,XpDUMMY,
     x                   E_total,E_nuc,E_ecore,E_pcore,E_ep,
     x                   E_ee,E_gam1,E_gam2,E_gam3,E_gam4,
     x                   S_total,S_gam1,S_gam2)

!     subroutine xchf1_fock(LNEOHF,LG4DSCF,LG3DSCF,NG4CHK,NG3CHK,
!    x                      npebf,nebf,nebf2,npbf,npbf2,nelec,
!    x                      ngee,ng1,ng2,ng3,ng4,DE,DP,
!    x                      GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
!    x                      ngtg1,ng3prm,nat,pmass,cat,zan,bcoef1,gamma1,
!    x                      KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
!    x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
!    x                      focke,fockp,Xfocke,Xfockp,
!    x                      E_total,E_nuc,E_ecore,E_pcore,E_ep,
!    x                      E_ee,E_gam1,E_gam2,E_gam3,E_gam4,
!    x                      S_total,S_gam1,S_gam2)

C-ELEC---------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------)

C-----------------------POSSIBLE-SOSCF---------------------------------(
c        if(SOSCF) THEN
         ITER=I
         EIGAVL = ITER.GT.1
         IF(SOSCF .AND.  EIGAVL) THEN
Cc          --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
            call pack_LT(nebf,nebfLT,focke,FLT)
            call SOGRAD(GRAD,FLT,vecE,WRK,NPR,NA,L0,L1,NEBFLT,ORBGRD)
cc          IF(ORBGRD.LT.SMALL) THEN
cc             DIFF = ZERO
cc             CVGING=.TRUE.
cc             GO TO 700  ! Check on convergence behavior
cc          END IF
            IF(ORBGRD.LT.SOGTOL  .OR.  ITSO.GT.0) THEN
               IF(ITSO.EQ.0) THEN
               WRITE(*,9800)
                  call SOHESS(HSTART,EE,NPR,L0,NA,NA)
               END IF
               ITSO = ITSO+1
               call SONEWT(HSTART,GRAD,PGRAD,DISPLI,DGRAD,DISPL,UPDT,
     *                     DISPLN,DGRADI,UPDTI,ORBGRD,NPR,ITSO,NFT15)
               call SOTRAN(DISPLI,vecE,G,WRK,NPR,L0,L1,NA,NA,ORBGRD)
               CALL DCOPY(NPR,GRAD,1,PGRAD,1)
               call make_density(nelec,nebf,vecE,DE)
               GO TO 800  ! Use the new C's to form new density (change)
            END IF
         END IF
C-----------------------POSSIBLE-SOSCF---------------------------------)

  700 CONTINUE
C        Focke diag
         call ROOTHAN(DE,vecE,EE,xxse,focke,nebf,nelec,1,NUCST)

  800 CONTINUE
C        --> FIND LARGEST CHANGE IN E DENSITY
         CALL DENDIF(DE0,DE,NEBF,DIFFE)
         CALL COPYDEN(DE0,DE,NEBF)

C        Fockp diag
c        call ROOTHAN(DP,vecP,EP,xxsp,fockp,npbf,1,2,NUCST)

C        --> FIND LARGEST CHANGE IN P DENSITY
         CALL DENDIF(DP0,DP,NPBF,DIFFP)
         CALL COPYDEN(DP0,DP,NPBF)

C        --> CALCULATE CHANGE IN TOTAL ENERGY
         Delta_E_tot=E_total-E_total_old
         E_total_old=E_total

C        --> PRINT SUMMARY OF THIS ITERATION
c        WRITE(*,*)I, E, DIFFE, DIFFP, EE(1), EP(1)
         WRITE(*,9100) I,E_total,Delta_E_tot,DIFFE,DIFFP,ORBGRD
! Output the vectors for this iteration for restart if necessary:
            call write_MOs(852,nebf,VECE)
            call write_MOs(853,npbf,VECP)

         if(npbf.gt.1) then
            IF((DIFFE.LT.TOLE).and.(DIFFP.LT.TOLP)) GOTO 100
         else
            IF(DIFFE.LT.TOLE) GOTO 100
         end if
         IF(I.EQ.MAXIT) GOTO 10
      END DO
 
  10  CONTINUE
C     IF WE GET HERE SOMETHING WENT WRONG

      if(SOSCF) THEN
         close(NFT15)
      end if

      WRITE(*,*)
!     WRITE(*,*)'ITERATION LIMIT EXCEEDED ... ABORTING'
      WRITE(*,*)'WARNING:  ITERATION LIMIT EXCEEDED'
      WRITE(*,*)
!     STOP
C
  100 CONTINUE
C     IF WE GET HERE WE ARE DONE - CONVERGENCE ACHIEVED

      if(SOSCF) THEN
         close(NFT15)
      end if

C     PRINT FINAL ENERGY AND PUNCH THE ORBITALS
      WRITE(*,9200) E_TOTAL,I
C
      WRITE(*,9300) E_nuc,E_ecore,E_pcore,E_ep,E_ee,
     xE_gam1,E_gam2,E_gam3,E_gam4,S_total,E_total

C  OUTPUT ELEC AND NUC EIGENVALUES AND EIGENVECTORS
      WRITE(*,9600)
      call PREVNU(vece,EE,nebf,nebf,nebf)
      WRITE(*,9700)
      call PREVNU(vecp,EP,npbf,npbf,npbf)

c        write(*,*)
c        write(*,*)'electron vectors'
c        write(*,*)
c        do l=1,nebf
c        do k=1,nebf
c        write(*,*)'Ce(',k,l,')=',vece(k,l)
c        end do
c        write(*,*)
c        end do
c        write(*,*)
c        write(*,*)'electron eigenvalues'
c        write(*,*)
c        do j=1,nebf
c        write(*,*)'ee(',j,')=',EE(J)
c        end do
c        write(*,*)
c        write(*,*)'proton vectors'
c        write(*,*)
c        do i=1,npbf
c        do l=1,npbf
c        write(*,*)'Cp(',l,i,')=',vecp(l,i)
c        end do
c        write(*,*)
c        end do
c        write(*,*)
c        write(*,*)'proton eigenvalues'
c        write(*,*)
c        do i=1,npbf
c        write(*,*)'ep(',i,')=',EP(I)
c        end do
C PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------(
c     subroutine write_MOs(IFIL,nbf,VEC)
C     IFIL=852 :: FinalCE.dat
C     IFIL=853 :: FinalCP.dat
      call write_MOs(852,nebf,VECE)
      call write_MOs(853,npbf,VECP)
C PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------)
C
!-----------NON-ADIABATIC-DIAGNOSTIC-----------------------------------(
c     if(LNONAD) then
c        call NONADD(LNEOHF,nebf,npebf,npbf,nocce,NUCST,ngtg1,
c    x               1.0d+00,PMASS,bcoef1,gamma1,
c    x               KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
c    x               ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
c    x               VECE,DE,VECP,DP)
c     end if
!-----------NON-ADIABATIC-DIAGNOSTIC-----------------------------------)

      RETURN
      END
C
C======================================================================
      SUBROUTINE guess_elec(nelec,nebf,xxse,GAM_ecore,DE)
C
C     Diagonalize the core electron Hamiltonian
C     to construct initial electron guess density
C======================================================================

      integer nebf
      integer nelec
      double precision xxse(nebf,nebf)
      double precision GAM_ecore(nebf,nebf)
      double precision DE(nebf,nebf)
C Local variables
      double precision C(nebf,nebf)
      double precision evf(nebf)

C n_particle=nelec
C i_particle=1
      call ROOTHAN(DE,C,EVF,xxse,GAM_ecore,nebf,nelec,1,1)

      RETURN
      END
C
C======================================================================
      SUBROUTINE guess_prot(NUCST,npbf,nebf,xxsp,GAM_pcore,GAM_ep,DE,DP)
C
C     Diagonalize the core proton Hamiltonian- 
C       ---has p and ep contributions---
C     to construct initial proton guess density
C======================================================================

      integer npbf
      integer nebf
      integer NUCST
      double precision xxsp(npbf,npbf)
      double precision GAM_pcore(npbf,npbf)
      double precision GAM_ep(npbf,npbf,nebf,nebf)
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
C Local variables
      double precision fockp(npbf,npbf)
      double precision C(npbf,npbf)
      double precision evf(npbf)
      double precision E_pcore

      fockp=0.0d+00
      call E_from_GAM_pcore(npbf,npbf*npbf,
     * GAM_pcore,DP,fockp,E_pcore)

c     subroutine E_from_GAM_pcore(npbf,npbf2,
c    * GAM_pcore,DP,fockp,E_pcore)


      call ROOTHAN(DP,C,EVF,xxsp,fockp,npbf,1,2,NUCST)

      RETURN
      END
C
C*MODULE GEMMOD  *DECK DENDIF
C======================================================================
      SUBROUTINE DENDIF(D0,D1,NB,DIFF)
C
C     ----- CALCULATE THE LARGEST ABSOLUTE CHANGE IN THE
C           DENSITY MATRIX -----
C======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION D0(NB,NB),D1(NB,NB)
C
      DIFF = 0.0D+00
      DO I=1,NB
         DO J=1,NB
            DIFIJ=ABS(D0(I,J)-D1(I,J))
            IF(DIFIJ.GT.DIFF) DIFF=DIFIJ
         END DO
      END DO
C
      RETURN
      END
C
C======================================================================
C*MODULE GEMMOD  *DECK COPYDEN
      SUBROUTINE COPYDEN(D0,D,N)
C
C     ---> COPY CURRENT DENSITY TO OLD DENSITY
C
C     D0 :: OLD DENSITY
C     D  :: CURRENT DENSITY
C     N  :: NUMBER OF BASIS FUNCTIONS
C
C======================================================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D0(N,N),D(N,N)

      DO I=1,N
         DO J=1,N
            D0(I,J)=D(I,J)
         END DO
      END DO
C
      RETURN
      END
C
C*MODULE GEMMOD  *DECK ROOTHAN
C======================================================================
      SUBROUTINE ROOTHAN(D,C,EVF,S,F,NB,n_particle,i_particle,NUCST)
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
C     D         ::  NEW DENSITY MATRIX
C     ============================== 
C
C     n_particle :: number of quantum particles of type i_particle
C     i_particle :: 1 for electron(s) 2 for quantum nucleus
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
C     DIAGONALIZE FOCK MATRIX 
C
      CALL RS(NB,NB,FP,EVF,2,CP,FV1,FV2,IERR)
C
C     TRANSFORM MO COEFFICIENTS
C      C = X.CP  
C   
      CALL MATMULT(NB,NB,NB,NB,X,CP,C)

C-----FORM-DENSITY-MATRIX----------------------------------------------(
      if(i_particle.eq.1) then
         if(n_particle.gt.1) then
            nocc=n_particle/2
            prefac=2.0d+00
         else
            nocc=1
            prefac=1.0d+00
         end if
         kstart=1
         klast=nocc
      end if

      if(i_particle.eq.2) then
         nocc=1
         prefac=1.0d+00
         kstart=NUCST
         klast=NUCST
      end if

      DO I=1,NB
         DO J=1,NB

          D(i,j)=0.0d+00
c         do k=1,nocc
          do k=kstart,klast
c           D(I,J)=C(I,1)*C(J,1)
            D(I,J) = D(i,j) + prefac*C(I,k)*C(J,k)
          end do

         END DO
      END DO
C-----FORM-DENSITY-MATRIX----------------------------------------------)

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
C
C=======================================================================
      subroutine invC(nebf,vecE,CCC)
C=======================================================================
      implicit none
      integer nebf
      double precision vecE(nebf,nebf),CCC(nebf,nebf)
      integer i,j

      do i=1,nebf
         do j=1,nebf
            CCC(i,j)=vecE(j,i)
         end do
      end do


      return
      end
