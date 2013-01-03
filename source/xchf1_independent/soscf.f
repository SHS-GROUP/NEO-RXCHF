C*MODULE SCFLIB  *DECK SOGRAD
c     SUBROUTINE SOGRAD(GRAD,F,C,WRK,NPR,NA,L0,L1,ORBGRD)
      SUBROUTINE SOGRAD(GRAD,F,C,WRK,NPR,NA,L0,L1,NEBFLT,ORBGRD)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
c     DIMENSION GRAD(NPR),F(*),C(L1,L0),WRK(L1)
      DIMENSION GRAD(NPR),F(NEBFLT),C(L1,L0),WRK(L1)
c     dimension C2(3,3),wrk2(3) 
C
      PARAMETER (ZERO=0.0D+00, MXSEQ=150)
C
c     LOGICAL GOPARR,DSKWRK,MASWRK,PARR
C
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c     COMMON /SCFOPT/ CONVHF,MAXIT,MCONV,NPUNCH,NPREO(4)
C
C     -----CALCULATE ORBITAL GRADIENT <OCC|FMO|VIRT> -----
C     WHICH IS C-DAGGER * F * C  TO FORM GRADIENT VECTOR -GRAD-
C          F      - FOCK MATRIX IN AO BASIS (TRIANGULAR)
C          C      - MATRIX OF MO COEFFICIENTS
C          L0     - NUMBER OF LINEARLY INDEPENDENT BASIS FUNCS (MO'S)
C          L1     - NUMBER OF BASIS FUNCTIONS (AO'S)
C          NA     - NUMBER OF OCCUPIED ORBITALS
C          NPR    - NUMBER OF ROTATION PARAMETERS
C          ORBGRD - MAXIMUM ORBITAL GRADIENT COMPONENT
C
CCWS(
C    For our purposes we hardwire convhf to the same value used in 
C    our xcscf.  CONVHF is the density tolerance.
      convhf=1.0d-10
CCWS)
      SMALL = 1.0D-03 * CONVHF
      NVIR = L0-NA
      CALL VCLR(WRK,1,L1)
      KG=0
c     PARR = GOPARR  .AND.  NA.GT.MXSEQ
C
C          MULTIPLY NA (OCCUPIED) FIRST ROWS OF C-DAGGER BY -F-
C
      DO 300 JO = 1,NA
C
C     ----- GO PARALLEL! -----
C
c        IF(PARR) THEN
c           IF(MOD(JO,NPROC).NE.ME) GO TO 300
c        END IF
C
         IK = 0
         DO 150 I = 1,L1
            IM1 = I-1
            DUM = ZERO
            CDUM = C(I,JO)
            IF (IM1.GT.0) THEN
               DO 100 K = 1,IM1
                  IK = IK+1
                  WRK(K) = WRK(K)+F(IK)*CDUM
                  DUM = DUM+F(IK)*C(K,JO)
  100          CONTINUE
            END IF
            IK = IK+1
            WRK(I) = DUM+F(IK)*CDUM
  150    CONTINUE
C
C           MULTIPLY THESE PRODUCT ROWS BY NVIR COLUMNS IN -C-,
C                   CORRESPONDING ONLY TO VIRTUAL ORBITALS
C
CCWS DEBUG(
c           write(*,*)'WRK=',WRK
c           write(*,*)'C=',C
c           write(*,*)'L1=',L1
c           write(*,*)'L0=',L0
c           write(*,*)'NA=',NA
c           write(*,*)'NVIR=',NVIR
c           write(*,*)'C(1,3)=',C(1,3)
c           write(*,*)'TEST1=',-DDOT(L1,WRK,1,C(1,3),1)
c           wrk2(1)=0.809571798828913
c           wrk2(2)=0.809571798828913
c           wrk2(3)=0.936265144208501
c           C2(1,3)=-0.867309539539584
c           write(*,*)'wrk2:',wrk2
c           write(*,*)'C2(1,3)=',C2(1,3)
c           TEST2=-DDOT(L1,WRK2,1,C2(1,3),1)
c           write(*,*)'TEST2=',TEST2 
CCWS DEBUG)
         KG = (JO-1)*NVIR
         DO 200 JV = NA+1,L0
            KG=KG+1
CCWS DEBUG(
c           write(*,*)'NVIR=',nvir
c           write(*,*)'WRK=',WRK
c           write(*,*)'KG=  JV=',KG,JV
cc          write(*,*)'C(1,JV)=',C(1,JV)
cc          write(*,*)'C(1,1)=',C(1,1)
cc          write(*,*)'C(1,2)=',C(1,2)
c           write(*,*)'C(1,3)=',C(1,3)
cc          write(*,*)'C(2,1)=',C(2,1)
cc          write(*,*)'C(2,2)=',C(2,2)
c           write(*,*)'C(2,3)=',C(2,3)
cc          write(*,*)'C(3,1)=',C(3,1)
cc          write(*,*)'C(3,2)=',C(3,2)
c           write(*,*)'C(3,3)=',C(3,3)
CCWS DEBUG)

            GRAD(KG) = -DDOT(L1,WRK,1,C(1,JV),1)
CCWS DEBUG(
c           write(*,*)'C(1,JV)=',C(1,JV)
c           write(*,*)'GRAD KG=',GRAD(KG)
CCWS DEBUG)
c           IF (ABS(GRAD(KG)).LT.SMALL) GRAD(KG)=ZERO
  200    CONTINUE
  300 CONTINUE
C
C     ----- GO PARALLEL! -----
C     COLLECT COLUMNS FROM EACH NODE TO MASTER, THEN BROADCAST TO ALL
C
c     IF(PARR) THEN
c        DO 350 I=1,NA
c           IFROM = MOD(I,NPROC)
c           IF(IFROM.EQ.MASTER) GO TO 350
c           II= (I-1)*NVIR+1
c           IF(MASWRK) THEN
c              CALL DDI_RECV(GRAD(II),8*NVIR,IFROM)
c           ELSE
c              IF(IFROM.EQ.ME) CALL DDI_SEND(GRAD(II),8*NVIR,MASTER)
c           END IF
c 350    CONTINUE
c        CALL DDI_BCAST(1011,'F',GRAD,NA*NVIR,MASTER)
c     END IF
C
CCWS DEBUG(
c     write(*,*)'IN SOGRAD:'
c     write(*,*)'GRAD VECTOR:'
c     write(*,*)GRAD
c     write(*,*)
CCWS DEBUG)
      IMAX   = IDAMAX(NPR,GRAD,1)
      ORBGRD = ABS(GRAD(IMAX))
      RETURN
      END
C
C*MODULE SCFLIB  *DECK SOHESS
      SUBROUTINE SOHESS(HSTART,EIG,NPR,L1,NA,NB)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION HSTART(NPR),EIG(L1)
C
      PARAMETER (ONE=1.0D+00, PT05=0.5D-02, TWENTY=ONE/PT05)
C
C     ----- SET APPROXIMATE INVERSE HESSIAN -----
C     DIAGONAL ELEMENTS OF INITIAL ORBITAL HESSIAN -HSTART-
C     ARE SET USING ONLY THE ENERGY EIGENVALUES -EIG-
C
C         L1  = NUMBER OF BASIS FUNCTIONS
C         NA  = NUMBER OF ALPHA MO'S
C         NB  = NUMBER OF BETA MO'S
C         NPR = NUMBER OF PAIR ROTATIONS
C
      K=0
      DO 110 J=1,NA
         IF (J.LE.NB) THEN
            N=NB+1
         ELSE
            N=NA+1
         END IF
         DO 100 I=N,L1
            K=K+1
            DE = EIG(I)-EIG(J)
            IF(DE.GT.PT05) THEN
               HSTART(K)=ONE/DE
            ELSE
               HSTART(K)=TWENTY
            END IF
  100    CONTINUE
  110 CONTINUE
      RETURN
      END
C*MODULE SCFLIB  *DECK SONEWT
      SUBROUTINE SONEWT(HSTART,GRAD,PGRAD,DISPLI,DGRAD,DISPL,UPDT,
     *                  DISPLN,DGRADI,UPDTI,ORBGRD,NPR,ITSO,NFT15)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION HSTART(NPR),GRAD(NPR),PGRAD(NPR),DISPLI(NPR),
     *          DGRAD(NPR),DISPL(NPR),UPDT(NPR),DISPLN(NPR),
     *          DGRADI(NPR),UPDTI(NPR)
C
c     LOGICAL GOPARR,DSKWRK,MASWRK,SVDSKW,GPSAVE
C
c     COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c     COMMON /SCFOPT/ CONVHF,MAXIT,MCONV,NPUNCH,NPREO(4)
C
      PARAMETER (ONE=1.0D+00)
C
C     ----- THIS ROUTINE PERFORMS INVERSE HESSIAN MATRIX UPDATES -----
C                USING THE DIRECT UPDATE PROCEDURE GIVEN IN
C           T.H.FISCHER, J.ALMLOF.  J.PHYS.CHEM.  1992,96,9768-74
C
C        HSTART - APPROXIMATE INITIAL DIAGONAL INVERSE HESSIAN
C        GRAD   - GRADIENT VECTOR AT CURRENT POINT
C        PGRAD  - GRADIENT FROM PREVIOUS POINT
C        DISPLI - ON EXIT, THE DISPLACEMENT VECTOR
C        ALL OTHER ARGUMENTS ARE SCRATCH STORAGE.
C
CCWS(
C    For our purposes we hardwire convhf to the same value used in 
C    our xcscf.  CONVHF is the density tolerance.
      convhf=1.0d-05
CCWS)
      SMALL = 1.0D-03 * CONVHF
      TOOBIG = 1.0D+00
      BIGROT = 0.1D+00
C
C     THIS ROUTINE CONTAINS WORK PROPORTIONAL TO N**2, AND SOME I/O,
C     AND THUS IS MOST EFFICIENTLY DONE BY THE MASTER NODE ONLY.
C
CCWS(
c     SVDSKW = DSKWRK
c     GPSAVE = GOPARR
c     DSKWRK = .FALSE.
c     GOPARR = .FALSE.
c     IF(.NOT.MASWRK) GO TO 500
CCWS)
C
C     INITIALIZE DISPLACEMENT VECTOR
C
      DO 100 I=1,NPR
         DISPLN(I)=HSTART(I)*GRAD(I)
  100 CONTINUE
C
      IF (ITSO.EQ.1) GO TO 400
C
      IF (ORBGRD.LT.SMALL) THEN
         ITSO = ITSO - 1
         GO TO 400
      END IF
C
      DO 120 I=1,NPR
         DGRADI(I)= GRAD(I)-PGRAD(I)
         UPDTI(I) = HSTART(I)*DGRADI(I)
  120 CONTINUE
C
      IF (ITSO.EQ.2) GO TO 300
C
C     CALCULATE VECTOR RECURSIVELY
C
      CALL SEQREW(NFT15)
C
      DO 200 K=1,ITSO-2
C
         CALL SQREAD(NFT15,DISPL,NPR)
         CALL SQREAD(NFT15,DGRAD,NPR)
         CALL SQREAD(NFT15, UPDT,NPR)
C
         S1 = DDOT(NPR, DISPL,1, DGRAD,1)
         S2 = DDOT(NPR, DGRAD,1,  UPDT,1)
         S3 = DDOT(NPR, DISPL,1,  GRAD,1)
         S4 = DDOT(NPR,  UPDT,1,  GRAD,1)
         S5 = DDOT(NPR, DISPL,1,DGRADI,1)
         S6 = DDOT(NPR,  UPDT,1,DGRADI,1)
C
         S1=ONE/S1
         S2=ONE/S2
         T=ONE+S1/S2
         T2=S1*S3
         T4=S1*S5
         T1=T*T2-S1*S4
         T3=T*T4-S1*S6
C
         CALL DAXPY(NPR, T1,DISPL,1,DISPLN,1)
         CALL DAXPY(NPR,-T2, UPDT,1,DISPLN,1)
         CALL DAXPY(NPR, T3,DISPL,1, UPDTI,1)
         CALL DAXPY(NPR,-T4, UPDT,1, UPDTI,1)
C
  200 CONTINUE
C
C     UPDATE VECTOR
C
  300 CONTINUE
      S1 = DDOT(NPR,DISPLI,1,DGRADI,1)
      S2 = DDOT(NPR,DGRADI,1, UPDTI,1)
      S3 = DDOT(NPR,DISPLI,1,  GRAD,1)
      S4 = DDOT(NPR, UPDTI,1,  GRAD,1)
C
      S1=ONE/S1
      S2=ONE/S2
      T=ONE+S1/S2
      T2=S1*S3
      T1=T*T2-S1*S4
C
      CALL DAXPY(NPR, T1,DISPLI,1,DISPLN,1)
      CALL DAXPY(NPR,-T2, UPDTI,1,DISPLN,1)
C
      CALL SQWRIT(NFT15,DISPLI,NPR)
      CALL SQWRIT(NFT15,DGRADI,NPR)
      CALL SQWRIT(NFT15, UPDTI,NPR)
C
C       KURT'S 1/2009 IDEA TO FIX PROBLEMS IF G-INVERSE ISN'T POS.DEF.,
C       NAMELY, IF STEP IS IN THE WRONG DIRECTION, SIMPLY ZERO IT OUT.
C       USING THIS PROBABLY CHANGES THE CONVERGENCE OF MANY TESTS,
C       INCLUDING (ADVERSELY, PROBABLY) EXAM33 AND EXAM36.  IT SHOULD
C       BE THOROUGHLY TESTED BEFORE ADAPTING IT.  THE SECOND LINE MAY
C       VERY WELL AVOID SUCH PROBLEMS?
C
C---  ZERO=0.0D+00
C---  DO I=1,NPR
C---     IF(DISPLN(I)*GRAD(I).LT.ZERO) DISPLN(I) = ZERO
CXXX     IF(DISPLN(I)*GRAD(I).LT.ZERO) DISPLN(I) = -DISPLN(I)/100.0D+00
C---  ENDDO
C
  400 CONTINUE
      DO I=1,NPR
         DISPLI(I) = -DISPLN(I)
      ENDDO
C
C        SCALE DISPLACEMENT SO THAT SQCDF DOESN'T EXCEED 0.1
C
      SQCDF = SQRT(DDOT(NPR,DISPLI,1,DISPLI,1)/NPR)
      IF(SQCDF.GT.TOOBIG  .AND.  ITSO.GT.5) THEN
C         KURT'S IDEA FROM JANUARY 2009 IS TO RESET THE HESSIAN,
C         BY TURNING OFF SOSCF ITERATION COUNTER, IN HOPES THAT
C         THE RUN MIGHT CURE ITSELF.  BOMBING IS VERY FATAL.
C-BOMB-         IF(MASWRK) WRITE(IW,9010) SQCDF
C-BOMB-         ISTAT=1
C-BOMB-         CALL DDI_BCAST(1016,'I',ISTAT,1,MASTER)
C-BOMB-         IF(NFG.EQ.0) THEN
C-BOMB-           CALL ABRT
C-BOMB-         ELSE
C-BOMB-C          ABORTING HERE RESULTS IN DEADLOCKS FOR GDDI
C-BOMB-           ORBGRD=-123
C-BOMB-           DSKWRK = SVDSKW
C-BOMB-           GOPARR = GPSAVE
C-BOMB-           RETURN
C-BOMB-         ENDIF
CCWS(
c        IF(MASWRK) WRITE(IW,9011) SQCDF
         WRITE(*,9011) SQCDF
CCWS)
         ITSO = 0
      END IF
C
C        FOR THE FIRST TWO ITERATIONS, BE MORE CONSERVATIVE
C          (THIS ALSO IS KURT'S IDEA, IT ADDS AN ITERATION TO EXAM09
C
C-KRG-  IF (SQCDF.GT.BIGROT/3.0D+00 .AND. ITSO.EQ.1) THEN
C-KRG-     IF(MASWRK) WRITE(IW,9020) SQCDF
C-KRG-     SCAL=SQRT(BIGROT/3.0D+00/SQCDF)
C-KRG-     CALL DSCAL(NPR,SCAL,DISPLI,1)
C-KRG-  END IF
C-KRG-  IF (SQCDF.GT.BIGROT/2.0D+00 .AND. ITSO.EQ.2) THEN
C-KRG-     IF(MASWRK) WRITE(IW,9020) SQCDF
C-KRG-     SCAL=SQRT(BIGROT/2.0D+00/SQCDF)
C-KRG-     CALL DSCAL(NPR,SCAL,DISPLI,1)
C-KRG-  END IF
C-KRG-  IF(SQCDF.GT.BIGROT          .AND. ITSO.GT.2) THEN
C-KRG-     IF(MASWRK) WRITE(IW,9020) SQCDF
C-KRG-     SCAL=SQRT(BIGROT/SQCDF)
C-KRG-     CALL DSCAL(NPR,SCAL,DISPLI,1)
C-KRG-  END IF
C
C         THE CODE IS GOING TO CALL -SOTRAN- AFTER -SONEWT-
C         SO ONE LAST SCALING WILL BE NEEDED, EVEN IF ITSO IS RESET.
      IF(SQCDF.GT.BIGROT) THEN
CCWS(
c        IF(MASWRK  .AND.  ITSO.GT.0) WRITE(IW,9020) SQCDF
         IF(ITSO.GT.0) WRITE(*,9020) SQCDF
CCWS)
         SCAL=SQRT(BIGROT/SQCDF)
         CALL DSCAL(NPR,SCAL,DISPLI,1)
      END IF
C
C-BOMB-      ISTAT=0
C-BOMB-      CALL DDI_BCAST(1016,'I',ISTAT,1,MASTER)
C
C        GIVE DISPLACEMENT TO ANY OTHER NODES
C
  500 CONTINUE
CCWS(
c     DSKWRK = SVDSKW
c     GOPARR = GPSAVE
CCWS)
C
C       OTHER NODES MUST CHECK TO SEE IF THE SOSCF RAN OK ON MASTER
C
C-BOMB-      IF(.NOT.MASWRK) CALL DDI_BCAST(1016,'I',ISTAT,1,MASTER)
C-BOMB-      IF(ISTAT.EQ.1) THEN
C-BOMB-         IF(NFG.EQ.0) THEN
C-BOMB-           CALL ABRT
C-BOMB-         ELSE
C-BOMB-           ORBGRD=-123
C-BOMB-           RETURN
C-BOMB-         ENDIF
C-BOMB-      END IF
C
c     IF(GOPARR) THEN
c        CALL DDI_BCAST(1017,'I',ITSO  ,1  ,MASTER)
c        CALL DDI_BCAST(1015,'F',DISPLI,NPR,MASTER)
c     END IF
      RETURN
C
C9010 FORMAT(1X,'SOSCF ENCOUNTERS A SERIOUS PROBLEM IN -SONEWT-'/
C    *       1X,'THE ROTATION ANGLE VECTOR HAS A HUGE NORM, SQCDF=',
C    *          1P,E12.3/1X,'REEXAMINE STARTING VECTORS, ',
C    *             'APPROPRIATENESS OF YOUR SCFTYP, ETC...')
 9011 FORMAT(1X,'*** RESETTTING SOSCF, UPON ENCOUNTERING A HUGE',
     *          ' TOTAL ROTATION=',1P,E12.3)
 9020 FORMAT(1X,'SOSCF IS SCALING ROTATION ANGLE MATRIX, SQCDF=',F12.6)
      END
C*MODULE SCFLIB  *DECK SOTRAN
      SUBROUTINE SOTRAN(X,C,G,WRK,NPR,L0,L1,NA,NB,ORBGRD)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION X(NPR),C(L1,L0),G(L0,L0),WRK(L0)
C
      PARAMETER (ZERO=0.0D+00, HALF=0.5D+00, ONE=1.0D+00, MXSEQ=150)
c     LOGICAL GOPARR,DSKWRK,MASWRK,PARR
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ----- APPROXIMATE SECOND ORDER ORBITAL TRANSFORMATION -----
C          X      = ROTATION PARAMETERS
C          G      = TRANSFORMATION MATRIX ( G = EXP(X) )
C          C      = MATRIX OF MO COEFFICIENTS TO BE ROTATED
C          L0     = NUMBER OF LINEARLY INDEPENDENT BASIS FUNCS (MO'S)
C          L1     = NUMBER OF BASIS FUNCTIONS (AO'S)
C          NA     = NUMBER OF ALPHA MO'S
C          NB     = NUMBER OF BETA MO'S
C          NPR    = NUMBER OF ROTATION PARAMETERS ( NPR=NOCC*NVIR)
C
C     ORIGINAL ROUTINE WRITTEN BY GALINA CHABAN AT ISU, 1995
C     PARALLELIZATION IS BY KAZUYA ISHIMURA AT IMS, 2004
C
      CALL VCLR(G,1,L0*L0)
CCWS(
c     PARR = GOPARR  .AND.  L0.GT.MXSEQ
CCWS)
C
C     ----- FORM MATRIX G = I + X + ( X*X/2 ) -----
C
      K=0
      DO 110 J=1,NA
         IF (J.LE.NB) THEN
            N=NB+1
         ELSE
            N=NA+1
         END IF
         DO 100 I=N,L0
            K=K+1
            G(J,I)= X(K)
            G(I,J)=-X(K)
  100    CONTINUE
  110 CONTINUE
C
C        TESTS SHOWED THE FULL SECOND ORDER EXPANSION IS OF LITTLE USE
C        SO WE ALWAYS SKIP IT, BUT SAVE THE CODE FOR FUTURE REFERENCE.
C        (SECOND ORDER EXPANSION CODE WORKS ONLY FOR NA=NB)
C
      IF(NA.NE.NB) GO TO 300
      IF(ORBGRD.GT.1.0D-20) GO TO 300
C
      NOCC=NA
      NVIR=L0-NOCC
C
      DO 230 I=1,NOCC
         DO 220 J=I,NOCC
            DUM=ZERO
            DO 210 K=1,NVIR
               LI=NVIR*(I-1)+K
               LJ=NVIR*(J-1)+K
               DUM=DUM+X(LI)*X(LJ)
  210       CONTINUE
            DUM = -DUM*HALF
            G(I,J)= G(I,J) + DUM
            G(J,I)= G(J,I) + DUM
  220    CONTINUE
  230 CONTINUE
C
      DO 260 I=1,NVIR
         DO 250 J=I,NVIR
            DUM=ZERO
            DO 240 K=1,NOCC
               LI=NVIR*(K-1)+I
               LJ=LI-I+J
               DUM=DUM+X(LI)*X(LJ)
  240       CONTINUE
            DUM = -DUM*HALF
            I1=NOCC+I
            J1=NOCC+J
            G(I1,J1)= G(I1,J1) + DUM
            G(J1,I1)= G(J1,I1) + DUM
  250    CONTINUE
  260 CONTINUE
C
  300 CONTINUE
      DO 310 I=1,L0
         G(I,I) = G(I,I) + ONE
  310 CONTINUE
C
C     ----- ORTHONORMALIZE THE TRANSFORMATION MATRIX -----
C
CCWS(
c     IF(PARR) THEN
c        DO 520 I=1,L0
c           DUM = DDOT(L0,G(1,I),1,G(1,I),1)
c           CALL DSCAL(L0,ONE/SQRT(DUM),G(1,I),1)
c           IF(I.EQ.L0) GO TO 520
c           IP1=I+1
c           DO 510 J=IP1,L0
c              IF(MOD(J,NPROC).NE.ME) GO TO 510
c              DUM = DDOT(L0,G(1,I),1,G(1,J),1)
c              CALL DAXPY(L0,-DUM,G(1,I),1,G(1,J),1)
c 510       CONTINUE
c           CALL DDI_BCAST(1011,'F',G(1,I+1),L0,MOD(I+1,NPROC))
c 520    CONTINUE
c     ELSE
         DO 540 I=1,L0
            DUM = DDOT(L0,G(1,I),1,G(1,I),1)
            CALL DSCAL(L0,ONE/SQRT(DUM),G(1,I),1)
            IF(I.EQ.L0) GO TO 540
            IP1=I+1
            DO 530 J=IP1,L0
               DUM = DDOT(L0,G(1,I),1,G(1,J),1)
               CALL DAXPY(L0,-DUM,G(1,I),1,G(1,J),1)
  530       CONTINUE
  540    CONTINUE
c     END IF
CCWS)
C
C     ----- ROTATE THE ORBITALS -C- BY TRANSFORMATION -G- -----
C
CCWS(
c     IF(PARR) THEN
c        DO 650 I=1,L1
c           DO 610 K=1,L0
c              WRK(K)=C(I,K)
c 610       CONTINUE
c           DO 630 J=1,L0
c              IF(MOD(J,NPROC).NE.ME) GO TO 630
c              DUM=ZERO
c              DO 620 K=1,L0
c                 DUM=DUM+WRK(K)*G(K,J)
c 620          CONTINUE
c              C(I,J)=DUM
c 630       CONTINUE
c 650    CONTINUE
C
c        DO 660 J=1,L0
c           JFROM = MOD(J,NPROC)
c           IF(JFROM.EQ.MASTER) GO TO 660
c           IF(MASWRK) THEN
c              CALL DDI_RECV(C(1,J),8*L1,JFROM)
c           ELSE
c              IF(JFROM.EQ.ME) CALL DDI_SEND(C(1,J),8*L1,MASTER)
c           END IF
c 660    CONTINUE
c        CALL DDI_BCAST(1012,'F',C,L1*L0,MASTER)
c     ELSE
         DO 750 I=1,L1
            DO 710 K=1,L0
               WRK(K)=C(I,K)
  710       CONTINUE
            DO 730 J=1,L0
               DUM=ZERO
               DO 720 K=1,L0
                  DUM=DUM+WRK(K)*G(K,J)
  720          CONTINUE
               C(I,J)=DUM
  730       CONTINUE
  750    CONTINUE
c     END IF
CCWS)
      RETURN
      END
C
C*MODULE MTHLIB  *DECK VCLR
      SUBROUTINE VCLR(A,INCA,N)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION A(*)
C
      PARAMETER (ZERO=0.0D+00)
C
C     ----- ZERO OUT VECTOR -A-, USING INCREMENT -INCA- -----
C
      IF (INCA .NE. 1) GO TO 200
      DO 110 L=1,N
         A(L) = ZERO
  110 CONTINUE
      RETURN
C
  200 CONTINUE
      LA=1-INCA
      DO 210 L=1,N
         LA=LA+INCA
         A(LA) = ZERO
  210 CONTINUE
      RETURN
      END
C*MODULE IOLIB   *DECK SQWRIT
      SUBROUTINE SQWRIT(LFILE,REGION,LENGTH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c     LOGICAL GOPARR,DSKWRK,MASWRK
      DIMENSION REGION(LENGTH)
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     WRITE AN ARRAY -REGION- OF LENGTH -LENGTH- TO UNIT -LFILE-
C     SEE ALSO -SQREAD- AND -SEQADV- WHICH ARE PARTNERS TO SQWRIT.
C
C         THE FLAG -DSKWRK- DISTINGUISHES THE PARALLEL I/O STRATEGY:
C           IF DSKWRK=.TRUE.,  ALL PROCESSES WRITE THEIR OWN FILES.
C           IF DSKWRK=.FALSE., ONLY THE MASTER WILL WRITE THE FILE.
C
CCWS(
c     IF (DSKWRK.OR.MASWRK) WRITE(LFILE,ERR=300) REGION
      WRITE(LFILE,ERR=300) REGION
CCWS)
      RETURN
C
  300 CONTINUE
CCWS(
c     WRITE(6,9020) ME,LFILE
      WRITE(*,9020) 0,LFILE
c     CALL ABRT
      CALL ABORT
CCWS)
      STOP
C
 9020 FORMAT(1X,'SQWRIT: NODE',I4,
     *          ' ENCOUNTERED I/O ERROR WRITING UNIT',I4)
      END
C*MODULE IOLIB   *DECK SQREAD
      SUBROUTINE SQREAD(LFILE,REGION,LENGTH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION REGION(LENGTH)
c     LOGICAL GOPARR,DSKWRK,MASWRK
c     COMMON /IOFILE/ IR,IW,IP,IIS,IPK,IDAFX,NAV,IODAX(950)
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     READ -LENGTH- FLOATING POINT WORDS INTO -REGION- FROM -LFILE-
C     A POSSIBLE END-OF-FILE RESULTS IN RETURNING -LENGTH- AS ZERO!
C     NOTE THAT DUE TO ITS POSSIBLE RESET TO 0, -LENGTH- SHOULD
C     BE A INTEGER VARIABLE, RATHER THAN AN INTEGER CONSTANT.
C     SEE ALSO -SQWRIT- AND -SEQADV- WHICH ARE PARTNERS TO SQREAD.
C
C         THE FLAG -DSKWRK- DISTINGUISHES THE PARALLEL I/O STRATEGY:
C           IF DSKWRK=.TRUE.,  ALL PROCESSES READ THEIR OWN FILES.
C           IF DSKWRK=.FALSE., ONLY THE MASTER WILL READ THE FILE,
C           BUT THE DATA IS THE BROADCAST TO ALL OTHER PROCESSES.
C
CCWS(
c     IF (DSKWRK.OR.MASWRK) READ(LFILE, END=200, ERR=300) REGION
       READ(LFILE, END=200, ERR=300) REGION
CCWS)
C
C         IF RUNNING IN PARALLEL, AND THE FILE EXISTS ONLY
C         ON THE MASTER NODE (DSKWRK=.FALSE.), THEN THE DATA
C         SHOULD BE BROADCAST FROM THE MASTER TO ALL OTHER NODES.
C
CCWS(
c     IF (GOPARR.AND.(.NOT.DSKWRK)) THEN
c        CALL DDI_BCAST(230,'F',REGION,LENGTH,MASTER)
c     END IF
CCWS)
      RETURN
C
C                  END OF FILE
C        THIS IS HANDLED BY RETURNING ZERO LENGTH READ, SO THE CALLER
C        CAN DETERMINE IF THIS IS REALLY AN ERROR, OR WAS EXPECTED.
C
  200 CONTINUE
      LENGTH=0
      RETURN
C
C                  ERROR READING FILE, PULL THE PLUG ON THE JOB
C
  300 CONTINUE
CCWS(
c     WRITE(IW,9000) LFILE,ME,LENGTH
      WRITE(*,9000) LFILE,0,LENGTH
c     CALL ABRT
      CALL ABORT
CCWS)
C
      RETURN
 9000 FORMAT(1X,'ERROR READING FILE',I4,' ON NODE',I5,' LENGTH=',I10)
      END
C*MODULE IOLIB   *DECK SEQREW
      SUBROUTINE SEQREW(IFILE)
C
c     LOGICAL GOPARR,DSKWRK,MASWRK
C
C     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C   ----- REWIND FILE IFILE ----
C
CCWS(
C     IF (DSKWRK.OR.MASWRK) REWIND (UNIT=IFILE, ERR=300)
      REWIND (UNIT=IFILE, ERR=300)
CCWS)
  300 CONTINUE
      RETURN
      END
C=======================================================================
      subroutine abort
C=======================================================================
      implicit none
      integer ierr

      write(*,*)'THERE WAS A PROBLEM.  STOPPING ALL PROCESSES'
      call mpi_final(ierr) 
      stop

      return
      end

C=======================================================================
      subroutine pack_LT(N,NLT,X,XLT)
C=======================================================================
      implicit none
      integer N,NLT
      double precision X(N,N),XLT(NLT)
      integer i,j,ia

      ia=0
      do i=1,N
         do j=1,i
            ia=ia+1
c           XLT(ia)=X(i,j)
            XLT(ia)=X(j,i)
         end do
      end do

      return
      end
C=======================================================================
      subroutine make_density(n_particle,nebf,C,DE)
C=======================================================================
      implicit none
      integer n_particle,nocc,nebf
      double precision DE(nebf,nebf),C(nebf,nebf)
      double precision prefac
      integer ie1,je1,k,kstart,klast

      if(n_particle.lt.2) then
         write(*,*)'IN make_density:'
         write(*,*)'n_particle must be gt 1 and even'
         call abort
      end if
C-----FORM-DENSITY-MATRIX----------------------------------------------(
c     if(i_particle.eq.1) then
c     if(n_particle.gt.1) then
         nocc=n_particle/2
         prefac=2.0d+00
c     else
c        nocc=1
c        prefac=1.0d+00
c     end if
      kstart=1
      klast=nocc
c     end if

c     if(i_particle.eq.2) then
c        nocc=1
c        prefac=1.0d+00
c        kstart=NUCST
c        klast=NUCST
c     end if

      DO ie1=1,nebf
         DO je1=1,nebf

            DE(ie1,je1)=0.0d+00

            do k=kstart,klast
               DE(ie1,je1) = DE(ie1,je1) + prefac*C(ie1,k)*C(je1,k)
            end do

         END DO
      END DO
C-----FORM-DENSITY-MATRIX----------------------------------------------)

      return
      end
C=======================================================================
      subroutine abrt
C=======================================================================
      implicit none
      integer ierr

      write(*,*)'THERE WAS A PROBLEM.  STOPPING ALL PROCESSES'
      call mpi_final(ierr) 
      stop


      return
      end

