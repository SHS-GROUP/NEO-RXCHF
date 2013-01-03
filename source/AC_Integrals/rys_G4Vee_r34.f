C  Taken from driver_xcneo.f
C*************************************************************************
      SUBROUTINE rys_G4Vee_r34(NQUAD_coul,NQUAD_ovlap,
     1                I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     1                I4,J4,K4,alp4,Amat4,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     3                L4,M4,N4,beta4,Bmat4,
     4                gamA14,gamA24,gamA34,
     4                gamB14,gamB24,gamB34,
     4                ans)
C*************************************************************************
C ans = <GA(1)GA(2)GA(3)GB(4)|g(1,4)g(2,4)g(3,4)/(r3-r4)|GB(1)GB(2)GB(3)GB(4)>
C
      implicit none
C input
      integer  NDIM
      parameter (NDIM=3)
      integer NQUAD_coul
      integer NQUAD_ovlap
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  I4,J4,K4
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3
      integer  L4,M4,N4
      double precision  alp1
      double precision  alp2
      double precision  alp3
      double precision  alp4
      double precision  beta1
      double precision  beta2
      double precision  beta3
      double precision  beta4
      double precision  Amat1(NDIM), Bmat1(NDIM)
      double precision  Amat2(NDIM), Bmat2(NDIM)
      double precision  Amat3(NDIM), Bmat3(NDIM)
      double precision  Amat4(NDIM), Bmat4(NDIM)
      double precision  gamA14,gamA24,gamA34
      double precision  gamB14,gamB24,gamB34
C output
      double precision ans
C local
      integer iloop
      double precision t
      double precision exact
      double precision tstep
      double precision xsum
      double precision xval
      double precision PI
      double precision  xgamB34

C     --quadrature--
      double precision x
      double precision u
      double precision rho
      double precision coef
      double precision xjaco
      double precision xn(NQUAD_coul)
      double precision wt(NQUAD_coul)
      double precision sx1,sx2,sx34
      integer isum
      isum = I1+J1+K1+
     1       I2+J2+K2+
     1       I3+J3+K3+
     2       L1+M1+N1+
     3       L2+M2+N2+
     3       L3+M3+N3+1
C------------------------------------------------------------------------------
      if(NQUAD_coul .eq. 0) then
         call gfovlap(I1,J1,K1,alp1,Amat1,
     2                L1,M1,N1,beta1,Bmat1,sx1)
         call gfovlap(I2,J2,K2,alp2,Amat2,
     2                L2,M2,N2,beta2,Bmat2,sx2)
         call gfvee(I3,J3,K3,alp3,Amat3,
     1              I4,J4,K4,alp4,Amat4,
     2              L3,M3,N3,beta3,Bmat3,
     3              L4,M4,N4,beta4,Bmat4,
     4              sx34)
         ans = sx1 * sx2 * sx34
         RETURN
      end if
C------------------------------------------------------------------------------
C initialize
      exact   =  ans
      PI      = dacos(-1.0d0)
C------------------------------------------------------------------------------
C     --get quadrature--
      call g5k15_kron(NQUAD_coul,(alp4+beta4)*0.50d0,xn,wt)
      x = (1.0d0/(alp1+beta1))
     1  + (1.0d0/(alp2+beta2))
     2  + (1.0d0/(alp3+beta3))
     3  + (1.0d0/(alp4+beta4))
      rho  = 1.0d0/x
      coef = dsqrt(rho/PI)
C------------------------------------------------------------------------------
      xsum = 0.0d0
      do iloop=1,(10*isum)*NQUAD_coul

         call XG4ovlap_typ2(NQUAD_ovlap,
     1                 I1,J1,K1,alp1,Amat1,
     1                 I2,J2,K2,alp2,Amat2,
     1                 I3,J3,K3,alp3,Amat3,
     1                 I4,J4,K4,alp4,Amat4,
     2                 L1,M1,N1,beta1,Bmat1,
     3                 L2,M2,N2,beta2,Bmat2,
     3                 L3,M3,N3,beta3,Bmat3,
     3                 L4,M4,N4,beta4,Bmat4,
     4                 gamA14,gamA24,gamA34,
     4                 gamB14,gamB24,xgamB34,xval)

      end do

      xsum = 0.0d0
      do iloop=1,NQUAD_coul
         u = xn(iloop)*dsqrt( rho/(1.0d0-(xn(iloop)*xn(iloop))) )
         xjaco = 1.0d0/(1.0d0-(xn(iloop)*xn(iloop)))
         xjaco = xjaco*dsqrt(xjaco)
         xgamB34 = gamB34 + (u*u)

         call XG4ovlap_typ2(NQUAD_ovlap,
     1                 I1,J1,K1,alp1,Amat1,
     1                 I2,J2,K2,alp2,Amat2,
     1                 I3,J3,K3,alp3,Amat3,
     1                 I4,J4,K4,alp4,Amat4,
     2                 L1,M1,N1,beta1,Bmat1,
     3                 L2,M2,N2,beta2,Bmat2,
     3                 L3,M3,N3,beta3,Bmat3,
     3                 L4,M4,N4,beta4,Bmat4,
     4                 gamA14,gamA24,gamA34,
     4                 gamB14,gamB24,xgamB34,xval)

         xsum = xsum + (wt(iloop)*xval*xjaco)

      end do
      ans = coef * xsum

      END

      SUBROUTINE g5k15_kron(nq,p,x,w)
      IMPLICIT NONE
      INTEGER nq
      REAL*8 x(nq)
      REAL*8 w(nq)

      INTEGER*4   n, i, j, m

      REAL*8      EPS, x1, x2, ck, cp, theta

      PARAMETER   (EPS=3.d-14)
      PARAMETER   (THETA=1.565807874005600d0)

      REAL*8      p,p1, p2, p3, pp, xl, xm, z, z1
      REAL*8 A
      REAL*8 ADIFF
      REAL*8 ADIFF1
      REAL*8 BIG
      REAL*8 CEPS
      REAL*8 CEPSF
      REAL*8 CEPST
      REAL*8 CRIT
      REAL*8 DIFF
      REAL*8 EFACT
      REAL*8 EPMACH
      REAL*8 ERROR
      REAL*8 EST
      REAL*8 EST1
      REAL*8 EST2
      REAL*8 ESTST(30)
      REAL*8 FACERR
      REAL*8 FIFTH
      REAL*8 FUN
      REAL*8 FX1
      REAL*8 FX2
      REAL*8 FX3
      REAL*8 FX3ST(30)
      REAL*8 FX4
      REAL*8 FX5
      REAL*8 FX5ST(30)
      INTEGER LEV
      INTEGER LEVTAG
      INTEGER NIM
      INTEGER NO
      INTEGER NOM
      INTEGER NUM
      REAL*8 PREDIF(30)
      REAL*8 QCEPS
      REAL*8 RUM
      REAL*8 SIM
      REAL*8 SQUANK
      REAL*8 SUM
      REAL*8 THIRD
      REAL*8 X3
      REAL*8 X3ST(30)
      REAL*8 X4
      REAL*8 X5
      REAL*8 X5ST(30)
      REAL*8 XZERO

CCWS
      error=0.0d0
      big=0.0d0
      a=0.0d0
      fx1=0.0d0
      fx2=0.0d0
      fx3=0.0d0
      fx4=0.0d0
      fx5=0.0d0
CCWS
      n=nq
      x1 = -1.0d0
      x2 =  1.0d0
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)

      DO 12 i=1,m
         z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
    1       CONTINUE
         p1=1.d0
         p2=0.d0
         DO 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
   11       CONTINUE
      pp=n*(z*p1-p2)/(z*z-1.d0)
      z1=z
      cp = cos(theta)
      ck = 1.0d0-dexp(-cp*p)

          z=z1-p1/pp
      IF(DABS(z-z1).gt.EPS) GOTO 1
      x(i)=xm-xl*z
      x(n+1-i)=xm+xl*z
      w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
      w(n+1-i)=w(i)
   12    CONTINUE
      do i=1,nq
         w(i) = w(i)+ck
      end do
      EPMACH = 0.0000000000075D+00
      SUM = 0.0D+00
      SIM = 0.0D+00
      CEPSF = 180.0D+00 * ERROR / ( BIG - A )
      CEPS = CEPSF
      ADIFF = 0.0D+00
      LEVTAG = -1
      FACERR = 1.0D+00
      XZERO = A
      EFACT = 0.0D+00
      NIM = 1
      LEV = 0
      X1 = A
      X5 = BIG
      X3 = 0.5D+00 * ( A + BIG )
C     FX1 = FUN ( X1 )
C     FX3 = FUN ( X3 )
C     FX5 = FUN ( X3 )
      NO = 3
      EST = FX1 + FX5 + 4.0D+00 * FX3
      IF ( CEPSF ) 295, 205, 295
205   LEVTAG = 0
      FACERR = 15.0D+00
      CEPS = EPMACH * ABS ( FX1 )
      IF ( FX1 ) 295, 210, 295
210   CEPS = EPMACH * ABS ( FX3 )
      LEVTAG = 3
      IF ( FX3 ) 295, 215, 295
215   CEPS = EPMACH * ABS ( FX5 )
      IF ( FX5 ) 295, 220, 295
220   CEPS = EPMACH
295   QCEPS = 0.25D+00 * CEPS
300   CONTINUE
      X2 = 0.5D+00 * ( X1 + X3 )
      X4 = 0.5D+00 * ( X3 + X5 )
C     FX2 = FUN ( X2 )
C     FX4 = FUN ( X4 )
      NO = NO + 2
      EST1 = FX1 + 4.0D+00 * FX2 + FX3
      EST2 = FX3 + 4.0D+00 * FX4 + FX5
      ADIFF1 = ADIFF
      DIFF = EST + EST - EST1 - EST2
      IF ( LEV - 30 ) 305, 800, 800
305   ADIFF = ABS ( DIFF )
      CRIT = ADIFF - CEPS
      IF ( CRIT ) 700, 700, 400
400   CONTINUE
      IF ( ADIFF1 - ADIFF ) 410, 410, 500
410   IF ( LEV - 5 ) 500, 415, 415
415   EFACT = EFACT + CEPS * ( X1 - XZERO ) * FACERR
      XZERO = X1
      FACERR = 15.0D+00
      IF ( ADIFF - 2.0D+00 * CEPS ) 420, 420, 425
420   CEPS = ADIFF
      LEVTAG = 0
      GO TO 780
425   IF ( ADIFF1 - ADIFF ) 435, 430, 435
430   CEPS = ADIFF
      GO TO 445
435   CEPS = 2.0D+00 * CEPS
      IF ( LEVTAG - 3 ) 440, 445, 445
440   LEVTAG = 2
445   QCEPS = 0.25D+00 * CEPS
500   CONTINUE
      NIM = 2 * NIM
      LEV = LEV + 1
      ESTST(LEV) = EST2
      X3ST(LEV) = X4
      X5ST(LEV) = X5
      FX3ST(LEV) = FX4
      FX5ST(LEV) = FX5
      PREDIF(LEV) = ADIFF
      X5 = X3
      X3 = X2
      FX5 = FX3
      FX3 = FX2
      EST = EST1
      GO TO 300
700   CONTINUE
      IF ( LEV ) 400, 400, 705
C
705   IF ( LEVTAG ) 800, 710, 710
710   CEPST = 15.0D+00 * CEPS
      IF ( CRIT ) 715, 800, 800
715   IF ( LEVTAG - 2 ) 720, 740, 750
720   IF ( ADIFF ) 800, 800, 725
725   IF ( ADIFF - QCEPS ) 730, 800, 800
730   IF ( ADIFF - CEPSF ) 770, 770, 735
735   LEVTAG = 0
      CEPS = ADIFF
      EFACT = EFACT + CEPST * ( X1 - XZERO )
      XZERO = X1
      GO TO 445
740   LEVTAG = 0
      IF ( ADIFF ) 765, 765, 725
C
C  LEVTAG = 3.
C
750   LEVTAG = 0
      IF ( ADIFF ) 775, 775, 730
765   CEPS = ADIFF1
      GO TO 775
770   LEVTAG = -1
      FACERR = 1.0D+00
      CEPS = CEPSF
775   EFACT = EFACT + CEPST * ( X1 - XZERO )
      XZERO = X1
780   CONTINUE
      QCEPS = 0.25D+00 * CEPS
800   CONTINUE
      SUM = SUM + ( EST1 + EST2 ) * ( X5 - X1 )
      IF ( LEVTAG ) 805, 810, 810
805   SIM = SIM + DIFF * ( X5 - X1 )
810   CONTINUE
905   NUM = NIM / 2
      NOM = NIM - 2 * NUM
      IF ( NOM ) 910, 915, 910
910   NIM = NUM
      LEV = LEV - 1
      GO TO 905
915   NIM = NIM + 1
      IF ( LEV ) 1100, 1100, 1000
C
1000  CONTINUE
      X1 = X5
      FX1 = FX5
      X3 = X3ST(LEV)
      X5 = X5ST(LEV)
      FX3 = FX3ST(LEV)
      FX5 = FX5ST(LEV)
      EST = ESTST(LEV)
      ADIFF = PREDIF(LEV)
      GO TO 300
1100  CONTINUE

      EFACT = EFACT + CEPS * ( BIG - XZERO ) * FACERR
      RUM = EFACT / 180.0D+00
      THIRD = SUM / 12.0D+00
      FIFTH = -SIM / 180.0D+00
      SQUANK = THIRD + FIFTH
      RETURN

      END


