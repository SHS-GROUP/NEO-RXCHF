C*************************************************************************
      SUBROUTINE interface_G3Vec(NQUAD_coul,
     1                I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     4                gamA12,gamA13,gamA23,
     4                gamB12,gamB13,gamB23,Cmat1,ans)
C*************************************************************************
C ans = <GA(1)GA(2)GA(3)|g(1,2)g(1,3)g(2,3)/(r1-r2)|GB(1)GB(2)GB(3)>
C
      implicit none
C input
      integer  NDIM
      parameter (NDIM=3)
      integer NQUAD_coul
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3
      double precision  alp1
      double precision  alp2
      double precision  alp3
      double precision  beta1
      double precision  beta2
      double precision  beta3
      double precision  Amat1(NDIM), Bmat1(NDIM)
      double precision  Amat2(NDIM), Bmat2(NDIM)
      double precision  Amat3(NDIM), Bmat3(NDIM)
      double precision  Cmat1(NDIM)
      double precision  gamA12,gamA13,gamA23
      double precision  gamB12,gamB13,gamB23
C output
      double precision ans

      integer JVAL
      double precision sx1,sx2,sx3


      JVAL = I1 + J1 + K1
     1     + I2 + J2 + K2
     2     + I3 + J3 + K3
     3     + L1 + M1 + N1
     3     + L2 + M2 + N2
     3     + L3 + M3 + N3

      if(NQUAD_coul .ne. 0) then
         if(JVAL .eq. 0) then
            call G3Vec(I1,J1,K1,alp1,Amat1,
     1              I2,J2,K2,alp2,Amat2,
     1              I3,J3,K3,alp3,Amat3,
     2              L1,M1,N1,beta1,Bmat1,
     3              L2,M2,N2,beta2,Bmat2,
     3              L3,M3,N3,beta3,Bmat3,
     4              gamA12,gamA13,gamA23,
     4              gamB12,gamB13,gamB23,Cmat1,ans)
         else
            call rys_G3Vec(NQUAD_coul,
     1                I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     4                gamA12,gamA13,gamA23,
     4                gamB12,gamB13,gamB23,Cmat1,ans)
         end if
      else
c        call gfvec(I1,J1,K1,alp1,Amat1,
c    2              L1,M1,N1,beta1,Bmat1,
c    4              Cmat1,sx1)

c        call gfovlap(I2,J2,K2,alp2,Amat2,
c    2                L2,M2,N2,beta2,Bmat2,sx2)

c        call gfovlap(I3,J3,K3,alp3,Amat3,
c    2                L3,M3,N3,beta3,Bmat3,sx3)
c        ans = sx1 * sx2 * sx3
      end if



      END


C*************************************************************************
      SUBROUTINE interface_G3Vee(NQUAD_coul,
     1                I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     4                gamA12,gamA13,gamA23,
     4                gamB12,gamB13,gamB23,ans)
C*************************************************************************
C ans = <GA(1)GA(2)GA(3)|g(1,2)g(1,3)g(2,3)/(r1-r2)|GB(1)GB(2)GB(3)>
C
      implicit none
C input
      integer  NDIM
      parameter (NDIM=3)
      integer NQUAD_coul
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3
      double precision  alp1
      double precision  alp2
      double precision  alp3
      double precision  beta1
      double precision  beta2
      double precision  beta3
      double precision  Amat1(NDIM), Bmat1(NDIM)
      double precision  Amat2(NDIM), Bmat2(NDIM)
      double precision  Amat3(NDIM), Bmat3(NDIM)
      double precision  gamA12,gamA13,gamA23
      double precision  gamB12,gamB13,gamB23
C output
      double precision ans

      integer JVAL
      double precision sx12,sx3

      JVAL = I1 + J1 + K1
     1     + I2 + J2 + K2
     2     + I3 + J3 + K3
     3     + L1 + M1 + N1
     3     + L2 + M2 + N2
     3     + L3 + M3 + N3

      if(NQUAD_coul .ne. 0) then
         if(JVAL .eq. 0) then
            call G3Vee(I1,J1,K1,alp1,Amat1,
     1              I2,J2,K2,alp2,Amat2,
     1              I3,J3,K3,alp3,Amat3,
     2              L1,M1,N1,beta1,Bmat1,
     3              L2,M2,N2,beta2,Bmat2,
     3              L3,M3,N3,beta3,Bmat3,
     4              gamA12,gamA13,gamA23,
     4              gamB12,gamB13,gamB23,ans)
         else
            call rys_G3Vee(NQUAD_coul,
     1              I1,J1,K1,alp1,Amat1,
     1              I2,J2,K2,alp2,Amat2,
     1              I3,J3,K3,alp3,Amat3,
     2              L1,M1,N1,beta1,Bmat1,
     3              L2,M2,N2,beta2,Bmat2,
     3              L3,M3,N3,beta3,Bmat3,
     4              gamA12,gamA13,gamA23,
     4              gamB12,gamB13,gamB23,ans)
         end if
      else
c        call gfvee(I1,J1,K1,alp1,Amat1,
c    1              I2,J2,K2,alp2,Amat2,
c    2              L1,M1,N1,beta1,Bmat1,
c    2              L2,M2,N2,beta2,Bmat2,sx12)

c        call gfovlap(I3,J3,K3,alp3,Amat3,
c    2                L3,M3,N3,beta3,Bmat3,sx3)
c        ans = sx12 * sx3
      end if


      END

C*************************************************************************
      SUBROUTINE rys_G3Vec(NQUAD_coul,
     1                I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     4                gamA12,gamA13,gamA23,
     4                gamB12,gamB13,gamB23,Cmat1,ans)
C*************************************************************************
C ans = <GA(1)GA(2)GA(3)|g(1,2)g(1,3)g(2,3)/(r1-r2)|GB(1)GB(2)GB(3)>
C
      implicit none
C input
      integer  NDIM
      parameter (NDIM=3)
      integer NQUAD_coul
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3
      double precision  alp1
      double precision  alp2
      double precision  alp3
      double precision  beta1
      double precision  beta2
      double precision  beta3
      double precision  Amat1(NDIM), Bmat1(NDIM)
      double precision  Amat2(NDIM), Bmat2(NDIM)
      double precision  Amat3(NDIM), Bmat3(NDIM)
      double precision  Cmat1(NDIM)
      double precision  gamA12,gamA13,gamA23
      double precision  gamB12,gamB13,gamB23
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
      double precision  xgamA12,xgamA13,xgamA23
      double precision  xgamB12,xgamB13,xgamB23

C     --quadrature--
      double precision x
      double precision u
      double precision rho
      double precision coef
      double precision xjaco
      double precision xn(NQUAD_coul)
      double precision wt(NQUAD_coul)
      double precision cx1
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
      xgamA12 =  gamA12
      xgamA13 =  gamA13
      xgamA23 =  gamA23

      xgamB12 = gamB12
      xgamB13 = gamB13
      xgamB23 = gamB23
      PI      = dacos(-1.0d0)
C------------------------------------------------------------------------------
C     --get quadrature--
      call g5k15_kron(NQUAD_coul,(alp3+beta3)*0.50d0,xn,wt)
      x = (1.0d0/(alp1+beta1))
     1  + (1.0d0/(alp2+beta2))
     2  + (1.0d0/(alp3+beta3))
      rho  = 1.0d0/x
      coef = dsqrt(rho/PI)
C------------------------------------------------------------------------------
C     write(*,*) '--INITIAL--',pf
C------------------------------------------------------------------------------
C           --g(re2,rp1)VEE(re1,re2)--         
      xsum = 0.0d0
      do iloop=1,(10*isum)*NQUAD_coul

         call XG3ovlap_typ2(I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     3                cx1,Cmat1,
     4                xgamA12,xgamA13,xgamA23,
     4                xgamB12,xgamB13,xgamB23,xval)

      end do

      xsum = 0.0d0
      do iloop=1,NQUAD_coul
         u = xn(iloop)*dsqrt( rho/(1.0d0-(xn(iloop)*xn(iloop))) )
         xjaco = 1.0d0/(1.0d0-(xn(iloop)*xn(iloop)))
         xjaco = xjaco*dsqrt(xjaco)
         cx1   = (u*u)+(PI*PI)

         call XG3ovlap_typ2(I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     3                cx1,Cmat1,
     4                xgamA12,xgamA13,xgamA23,
     4                xgamB12,xgamB13,xgamB23,xval)

         xsum = xsum + (wt(iloop)*xval*xjaco)

      end do
      ans = coef * xsum
C     write(*,*) '--ANS--',ans,exact

C     STOP '****END OF NUMER***'
      END

C*************************************************************************
      SUBROUTINE rys_G3Vee(NQUAD_coul,
     1                I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     4                gamA12,gamA13,gamA23,
     4                gamB12,gamB13,gamB23,ans)
C*************************************************************************
C ans = <GA(1)GA(2)GA(3)|g(1,2)g(1,3)g(2,3)/(r1-r2)|GB(1)GB(2)GB(3)>
C    
      implicit none
C input
      integer  NDIM
      parameter (NDIM=3)
      integer NQUAD_coul
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3
      double precision  alp1
      double precision  alp2
      double precision  alp3
      double precision  beta1
      double precision  beta2
      double precision  beta3
      double precision  Amat1(NDIM), Bmat1(NDIM)
      double precision  Amat2(NDIM), Bmat2(NDIM)
      double precision  Amat3(NDIM), Bmat3(NDIM)
      double precision  gamA12,gamA13,gamA23
      double precision  gamB12,gamB13,gamB23
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
      double precision  xgamA12,xgamA13,xgamA23
      double precision  xgamB12,xgamB13,xgamB23

C     --quadrature--
      double precision x
      double precision u
      double precision rho
      double precision coef
      double precision xjaco
      double precision xn(NQUAD_coul)
      double precision wt(NQUAD_coul)

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
      xgamA12 =  gamA12
      xgamA13 =  gamA13
      xgamA23 =  gamA23

      xgamB12 = gamB12
      xgamB13 = gamB13
      xgamB23 = gamB23
      PI      = dacos(-1.0d0)
C------------------------------------------------------------------------------
C     --get quadrature--
      call g5k15_kron(NQUAD_coul,(alp3+beta3)*0.50d0,xn,wt)
      x = (1.0d0/(alp1+beta1))
     1  + (1.0d0/(alp2+beta2))
     2  + (1.0d0/(alp3+beta3))
      rho  = 1.0d0/x
      coef = dsqrt(rho/PI)
C------------------------------------------------------------------------------
C     write(*,*) '--INITIAL--',pf
C------------------------------------------------------------------------------
C           --g(re2,rp1)VEE(re1,re2)--         

      do iloop=1,(10*isum)*NQUAD_coul

C        call rys_G3ovlap(I1,J1,K1,alp1,Amat1,
C        call XG3ovlap_typ1(I1,J1,K1,alp1,Amat1,
         call G3ovlap(I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     4                xgamA12,xgamA13,xgamA23,
     4                xgamB12,xgamB13,xgamB23,xval)

      end do

      xsum = 0.0d0
      do iloop=1,NQUAD_coul
         u = xn(iloop)*dsqrt( rho/(1.0d0-(xn(iloop)*xn(iloop))) )
         xjaco = 1.0d0/(1.0d0-(xn(iloop)*xn(iloop)))
         xjaco = xjaco*dsqrt(xjaco)
         xgamB12 = gamB12+(u*u)+(PI*PI)

C        call rys_G3ovlap(I1,J1,K1,alp1,Amat1,
C        call XG3ovlap_typ1(I1,J1,K1,alp1,Amat1,
         call G3ovlap(I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     4                xgamA12,xgamA13,xgamA23,
     4                xgamB12,xgamB13,xgamB23,xval)

         xsum = xsum + (wt(iloop)*xval*xjaco)

      end do


      ans = coef * xsum
C     write(*,*) '--ANS--',ans,exact

C     STOP '****END OF NUMER***'
      END


