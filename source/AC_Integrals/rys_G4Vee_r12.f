C*************************************************************************
      SUBROUTINE rys_G4Vee_r12(NQUAD_coul,NQUAD_ovlap,
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
C ans = <GA(1)GA(2)GA(3)GB(4)|g(1,4)g(2,4)g(3,4)/(r1-r2)|GB(1)GB(2)GB(3)GB(4)>
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
      double precision  gamA12
      double precision  gamB12
      double precision  xgamB12

C     --quadrature--
      double precision x
      double precision u
      double precision rho
      double precision coef
      double precision xjaco
      double precision xn(NQUAD_coul)
      double precision wt(NQUAD_coul)
      double precision sx12,sx3,sx4
      integer isum
      isum = I1+J1+K1+
     1       I2+J2+K2+
     1       I3+J3+K3+
     2       L1+M1+N1+
     3       L2+M2+N2+
     3       L3+M3+N3+1
C------------------------------------------------------------------------------
C initialize
      exact   =  ans
C------------------------------------------------------------------------------
      if(NQUAD_coul .eq. 0) then
         call gfvee(I1,J1,K1,alp1,Amat1,
     1              I2,J2,K2,alp2,Amat2,
     2              L1,M1,N1,beta1,Bmat1,
     3              L2,M2,N2,beta2,Bmat2,
     4              sx12)

         call gfovlap(I3,J3,K3,alp3,Amat3,
     2                L3,M3,N3,beta3,Bmat3,sx3)

         call gfovlap(I4,J4,K4,alp4,Amat4,
     2                L4,M4,N4,beta4,Bmat4,sx4)
         ans = sx12 * sx3 * sx4
         RETURN
      end if
C------------------------------------------------------------------------------
C initialize
      exact   =  ans
      PI      = dacos(-1.0d0)
      gamA12  = 0.0d0
      gamB12  = 0.0d0
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

         call XG4ovlap_typ1(NQUAD_ovlap,
     1                 I1,J1,K1,alp1,Amat1,
     1                 I2,J2,K2,alp2,Amat2,
     1                 I3,J3,K3,alp3,Amat3,
     1                 I4,J4,K4,alp4,Amat4,
     2                 L1,M1,N1,beta1,Bmat1,
     3                 L2,M2,N2,beta2,Bmat2,
     3                 L3,M3,N3,beta3,Bmat3,
     3                 L4,M4,N4,beta4,Bmat4,
     4                 gamA12,gamA14,gamA24,gamA34,
     4                 xgamB12,gamB14,gamB24,gamB34,xval)

      end do

      xsum = 0.0d0
      do iloop=1,NQUAD_coul
         u = xn(iloop)*dsqrt( rho/(1.0d0-(xn(iloop)*xn(iloop))) )
         xjaco = 1.0d0/(1.0d0-(xn(iloop)*xn(iloop)))
         xjaco = xjaco*dsqrt(xjaco)
         xgamB12 = gamB12+(u*u)

C        call rys_G4ovlap(I1,J1,K1,alp1,Amat1,
C    1                I2,J2,K2,alp2,Amat2,
C    1                I3,J3,K3,alp3,Amat3,
C    1                I4,J4,K4,alp4,Amat4,
C    2                L1,M1,N1,beta1,Bmat1,
C    3                L2,M2,N2,beta2,Bmat2,
C    3                L3,M3,N3,beta3,Bmat3,
C    4                L4,M4,N4,beta4,Bmat4,
C    4                gamA12,gamA13,gamA14,gamA23,gamA24,gamA34,
C    4                xgamB12,gamB13,gamB14,gamB23,gamB24,gamB34,
C    4                xval)
         call XG4ovlap_typ1(NQUAD_ovlap,
     1                 I1,J1,K1,alp1,Amat1,
     1                 I2,J2,K2,alp2,Amat2,
     1                 I3,J3,K3,alp3,Amat3,
     1                 I4,J4,K4,alp4,Amat4,
     2                 L1,M1,N1,beta1,Bmat1,
     3                 L2,M2,N2,beta2,Bmat2,
     3                 L3,M3,N3,beta3,Bmat3,
     3                 L4,M4,N4,beta4,Bmat4,
     4                 gamA12,gamA14,gamA24,gamA34,
     4                 xgamB12,gamB14,gamB24,gamB34,xval)

         xsum = xsum + (wt(iloop)*xval*xjaco)

      end do
      ans = coef * xsum

      END

