C23456789
C*****************************************************************
      SUBROUTINE XG3ovlap_typ1(MAXQ,
     1                 I1,J1,K1,alp1,Amat1,
     1                 I2,J2,K2,alp2,Amat2,
     1                 I3,J3,K3,alp3,Amat3,
     2                 L1,M1,N1,beta1,Bmat1,
     3                 L2,M2,N2,beta2,Bmat2,
     3                 L3,M3,N3,beta3,Bmat3,
     4                 gamA12,gamA13,gamA23,
     4                 gamB12,gamB13,gamB23,ans)
C*****************************************************************
C
      implicit none
C Input variables
      integer NDIM
      parameter (NDIM=3)
      integer MAXQ
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      double precision alp1,Amat1(NDIM)
      double precision alp2,Amat2(NDIM)
      double precision alp3,Amat3(NDIM)
      double precision beta1,Bmat1(NDIM)
      double precision beta2,Bmat2(NDIM)
      double precision beta3,Bmat3(NDIM)
      double precision gamA12,gamA13,gamA23
      double precision gamB12,gamB13,gamB23
C Output variables
      double precision ans

C Local variables
      integer NQX,NQY,NQZ
      double precision sx,sy,sz

C----------------------------------------------------------------
C initialize quadrature
C
      NQX = I1+I2+I3
      NQY = J1+J2+J3
      NQZ = K1+K2+K3
      NQX = max(NQX,5)
      NQY = max(NQY,5)
      NQZ = max(NQZ,5)
      if(MAXQ .gt. 0) then
         NQX = MAXQ
         NQY = MAXQ
         NQZ = MAXQ
      end if
C----------------------------------------------------------------
      call G3ovlap_typ1_1D(NQX,
     1                     I1,alp1,Amat1(1),
     1                     I2,alp2,Amat2(1),
     1                     I3,alp3,Amat3(1),
     2                     L1,beta1,Bmat1(1),
     3                     L2,beta2,Bmat2(1),
     4                     L3,beta3,Bmat3(1),
     5                     gamA12,gamA13,gamA23,
     5                     gamB12,gamB13,gamB23,
     6                     sx)

      call G3ovlap_typ1_1D(NQY,
     1                     J1,alp1,Amat1(2),
     1                     J2,alp2,Amat2(2),
     1                     J3,alp3,Amat3(2),
     2                     M1,beta1,Bmat1(2),
     3                     M2,beta2,Bmat2(2),
     4                     M3,beta3,Bmat3(2),
     5                     gamA12,gamA13,gamA23,
     5                     gamB12,gamB13,gamB23,
     6                     sy)


      call G3ovlap_typ1_1D(NQZ,
     1                     K1,alp1,Amat1(3),
     1                     K2,alp2,Amat2(3),
     1                     K3,alp3,Amat3(3),
     2                     N1,beta1,Bmat1(3),
     3                     N2,beta2,Bmat2(3),
     4                     N3,beta3,Bmat3(3),
     5                     gamA12,gamA13,gamA23,
     5                     gamB12,gamB13,gamB23,
     6                     sz)

      ans = sx*sy*sz


      END

C*****************************************************************
      SUBROUTINE G3ovlap_typ1_1D(NQ,
     1                      I1,alp1,Amat1,
     1                      I2,alp2,Amat2,
     1                      I3,alp3,Amat3,
     2                      L1,beta1,Bmat1,
     3                      L2,beta2,Bmat2,
     4                      L3,beta3,Bmat3,
     5                      gamA12,gamA13,gamA23,
     5                      gamB12,gamB13,gamB23,
     6                      ans)
C*****************************************************************
C Calculates the following integral
C <GA(1)GA(2)GA(3)|g(1,2)g(1,3)g(2,3)|GB(1)GB(2)GB(3)>
C
      implicit none
C Input variables
      integer NQ
      integer I1
      integer I2
      integer I3
      integer I4
      integer L1
      integer L2
      integer L3
      integer L4
      double precision alp1,Amat1
      double precision alp2,Amat2
      double precision alp3,Amat3
      double precision alp4,Amat4
      double precision beta1,Bmat1
      double precision beta2,Bmat2
      double precision beta3,Bmat3
      double precision beta4,Bmat4
      double precision gamA12,gamA13,gamA23
      double precision gamB12,gamB13,gamB23
C Output variables
      double precision ans
C Local variables
      double precision gam12,gam13,gam23
C     --quadrature--
      integer i
      double precision xval12
      double precision xval3
      double precision xval4
      double precision expo3 
      double precision poly3 
      double precision expo12
      double precision poly12
      double precision xsum
      double precision p3,Pmat3,coef3
      double precision c1,c2,c3
      double precision tx1,tx2,tx3
      double precision xn(NQ),wt(NQ)
      double precision sx12,sx3,sx4,xneo
C----------------------------------------------------------------
C analytical
C     call gfovlap_1D(I3,alp3,Amat3,
C    2                L3,beta3,Bmat3,sx3)

C     call pgiovlap_1D(I1,alp1,Amat1,
C    1                 I2,alp2,Amat2,
C    2                 L1,beta1,Bmat1,
C    3                 L2,beta2,Bmat2,
C    4                 gamA12,gamB12,sx12)

C     xneo = sx12*sx3
C     if(dabs(xneo) .eq. 0.0d0) RETURN
C----------------------------------------------------------------
   
C----------------------------------------------------------------
C initialize
      gam12 = gamA12 + gamB12
      gam13 = gamA13 + gamB13 
      gam23 = gamA23 + gamB23 
C----------------------------------------------------------------
C Gaussian prod of GA(3) and GB(3)
C GP(3) = GA(3)*GB(3)
      p3 = alp3 + beta3
      Pmat3 = ((alp3*Amat3)+(beta3*Bmat3))/p3
      coef3 = alp3*beta3*(Amat3-Bmat3)*(Amat3-Bmat3)/p3
      coef3 = dexp(-coef3)
C     call gp(alp3,beta3,Amat3,Bmat3,coef3,p3,Pmat3) 
C----------------------------------------------------------------
      call hermite_quad(NQ,1.0d0,p3,Pmat3,xn,wt)
      xsum = 0.0d0
      do i=1,NQ 
         call multiC_G2ovlap_1D(I1,alp1,Amat1,
     1                          I2,alp2,Amat2,
     2                          L1,beta1,Bmat1,
     3                          L2,beta2,Bmat2,
     3                          gam13,xn(i),
     3                          gam23,xn(i),
     5                          gam12,xval12)

         xval3 = ((xn(i)-Amat3)**I3)*((xn(i)-Bmat3)**L3)

         xsum = xsum + (wt(i)*xval12*xval3)
C        xsum = xsum + (wt(i)*xval12*xval3/xneo)
      end do
      ans = coef3 * xsum 
C     ans = coef3 * xsum * xneo

      END
C***************************************************************************
      SUBROUTINE gp(alp,beta,A,B,coef,p,cent)
C***************************************************************************
      implicit none
C     --input--
      double precision alp
      double precision beta
      double precision A
      double precision B
C     --output--
      double precision coef
      double precision p
      double precision cent

      p = alp + beta
      cent = ((alp*A) + (beta*B))/p

      coef = alp*beta*(A-B)*(A-B)/p
      coef = dexp(-coef)
      

      END




