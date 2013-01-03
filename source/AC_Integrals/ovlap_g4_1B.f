C23456789
C*****************************************************************
      SUBROUTINE XG4ovlap_typ2(MAXQ,
     1                 I1,J1,K1,alp1,Amat1,
     1                 I2,J2,K2,alp2,Amat2,
     1                 I3,J3,K3,alp3,Amat3,
     1                 I4,J4,K4,alp4,Amat4,
     2                 L1,M1,N1,beta1,Bmat1,
     3                 L2,M2,N2,beta2,Bmat2,
     3                 L3,M3,N3,beta3,Bmat3,
     3                 L4,M4,N4,beta4,Bmat4,
     4                 gamA14,gamA24,gamA34,
     4                 gamB14,gamB24,gamB34,ans)
C*****************************************************************
C
      implicit none
C Input variables
      integer MAXQ
      integer NDIM
      parameter (NDIM=3)
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4
      double precision alp1,Amat1(NDIM)
      double precision alp2,Amat2(NDIM)
      double precision alp3,Amat3(NDIM)
      double precision alp4,Amat4(NDIM)
      double precision beta1,Bmat1(NDIM)
      double precision beta2,Bmat2(NDIM)
      double precision beta3,Bmat3(NDIM)
      double precision beta4,Bmat4(NDIM)
      double precision gamA14,gamA24,gamA34
      double precision gamB14,gamB24,gamB34
C Output variables
      double precision ans

C Local variables
      integer NQX,NQY,NQZ
      double precision sx,sy,sz

C----------------------------------------------------------------
C initialize quadrature
C
      NQX = I1+I2+I3+I4
      NQY = J1+J2+J3+J4
      NQZ = K1+K2+K3+K4
      NQX = max(NQX,5)
      NQY = max(NQY,5)
      NQZ = max(NQZ,5)
      if(MAXQ .gt. 0) then
         NQX = MAXQ
         NQY = MAXQ
         NQZ = MAXQ
      end if
C----------------------------------------------------------------
      call G4ovlap_typ2_1D(NQX,
     1                     I1,alp1,Amat1(1),
     1                     I2,alp2,Amat2(1),
     1                     I3,alp3,Amat3(1),
     1                     I4,alp4,Amat4(1),
     2                     L1,beta1,Bmat1(1),
     3                     L2,beta2,Bmat2(1),
     4                     L3,beta3,Bmat3(1),
     4                     L4,beta4,Bmat4(1),
     5                     gamA14,gamA24,gamA34,
     5                     gamB14,gamB24,gamB34,
     6                     sx)

      call G4ovlap_typ2_1D(NQY,
     1                     J1,alp1,Amat1(2),
     1                     J2,alp2,Amat2(2),
     1                     J3,alp3,Amat3(2),
     1                     J4,alp4,Amat4(2),
     2                     M1,beta1,Bmat1(2),
     3                     M2,beta2,Bmat2(2),
     4                     M3,beta3,Bmat3(2),
     4                     M4,beta4,Bmat4(2),
     5                     gamA14,gamA24,gamA34,
     5                     gamB14,gamB24,gamB34,
     6                     sy)


      call G4ovlap_typ2_1D(NQZ,
     1                     K1,alp1,Amat1(3),
     1                     K2,alp2,Amat2(3),
     1                     K3,alp3,Amat3(3),
     1                     K4,alp4,Amat4(3),
     2                     N1,beta1,Bmat1(3),
     3                     N2,beta2,Bmat2(3),
     4                     N3,beta3,Bmat3(3),
     4                     N4,beta4,Bmat4(3),
     5                     gamA14,gamA24,gamA34,
     5                     gamB14,gamB24,gamB34,
     6                     sz)

      ans = sx*sy*sz


      END

C*****************************************************************
      SUBROUTINE G4ovlap_typ2_1D(NQ,
     1                      I1,alp1,Amat1,
     1                      I2,alp2,Amat2,
     1                      I3,alp3,Amat3,
     1                      I4,alp4,Amat4,
     2                      L1,beta1,Bmat1,
     3                      L2,beta2,Bmat2,
     4                      L3,beta3,Bmat3,
     4                      L4,beta4,Bmat4,
     5                      gamA14,gamA24,gamA34,
     5                      gamB14,gamB24,gamB34,
     6                      ans)
C*****************************************************************
C Calculates the following integral
C <GA(1)GA(2)GA(3)GA(4)|g(1,2)g(1,4)g(2,4)g(3,4)|GB(1)GB(2)GB(3)GB(4)>
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
      double precision gamA14,gamA24,gamA34
      double precision gamB14,gamB24,gamB34
C Output variables
      double precision ans
C Local variables
      double precision gam12,gam14,gam24,gam34
C     --quadrature--
      integer i
      double precision xval1
      double precision xval2
      double precision xval3
      double precision xval4
      double precision xsum
      double precision p4,Pmat4,coef4
      double precision xn(NQ),wt(NQ)
      double precision sx1,sx2,sx3,sx4,xneo
C----------------------------------------------------------------
C analytical
C     call gfovlap_1D(I1,alp1,Amat1,
C    2                L1,beta1,Bmat1,sx1)

C     call gfovlap_1D(I2,alp2,Amat2,
C    2                L2,beta2,Bmat2,sx2)

C     call gfovlap_1D(I3,alp3,Amat3,
C    2                L3,beta3,Bmat3,sx3)

C     call gfovlap_1D(I4,alp4,Amat4,
C    2                L4,beta4,Bmat4,sx4)

C     xneo = sx1*sx2*sx3*sx4

C     if(dabs(xneo) .ne. 0.0d0) then
C----------------------------------------------------------------
C initialize
         gam14 = gamA14 + gamB14 
         gam24 = gamA24 + gamB24 
         gam34 = gamA34 + gamB34 
C----------------------------------------------------------------
C Gaussian prod
         p4 = alp4 + beta4
         Pmat4 = ((alp4*Amat4)+(beta4*Bmat4))/p4
         coef4 = alp4*beta4*(Amat4-Bmat4)*(Amat4-Bmat4)/p4
         coef4 = dexp(-coef4)
C----------------------------------------------------------------
         call hermite_quad(NQ,1.0d0,p4,Pmat4,xn,wt)
         xsum = 0.0d0
         do i=1,NQ 
            call multiC_G1ovlap_1D(I1,alp1,Amat1,
     2                          L1,beta1,Bmat1,
     3                          gam14,xn(i),xval1)

            call multiC_G1ovlap_1D(I2,alp2,Amat2,
     2                          L2,beta2,Bmat2,
     3                          gam24,xn(i),xval2)

            call multiC_G1ovlap_1D(I3,alp3,Amat3,
     2                          L3,beta3,Bmat3,
     3                          gam34,xn(i),xval3)

            xval4 = ((xn(i)-Amat4)**I4)*((xn(i)-Bmat4)**L4)

            xsum = xsum + (wt(i)*xval1*xval2*xval3*xval4)
C           xsum = xsum + (wt(i)*xval1*xval2*xval3*xval4/xneo)
         end do
         ans = coef4 * xsum 
C        ans = coef4 * xsum * xneo
C     else
C        ans = xneo
C     end if


      END
