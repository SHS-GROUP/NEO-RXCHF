C*****************************************************************
      SUBROUTINE G3ke_g1Tpg2(I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     4                gamA12,gamA13,gamA23,
     4                gamB12,gamB13,gamB23,xmass,ansKE)
C*****************************************************************
C Geminal 3-particle KE  matrix element
C <GA(1)GA(2)GA(3)g(1,3)|T(3)|GB(1)GB(2)GB(3)g(2,3)>
C
C The Laplacian is calculated with respect to the coordinates of 
C particle 3
C
C------------Definitions----------
C GA and Gb are Cartesian Gaussian functions
C GA(1) : (x-Ax)^i
C GA(2) : (y-Ay)^j
C GA(3) : (z-Az)^k
C GA(4) : alpha
C GA(5) : Ax
C GA(6) : Ay
C GA(7) : Az
C
C GA(r1) = (x-Ax)^i (y-Ay)^j (z-Az)^k Exp[-alpha*{(x-Ax)^2+(y-Ay)^2+(z-Az)^2}]
C GB(r1) = (x-Bx)^l (y-By)^m (z-Bz)^n Exp[-beta *{(x-Bx)^2+(y-By)^2+(z-Bz)^2}]
C 
C gamA: Geminal coupling between GA(r1) and GA(r2) == Exp[-GamA*r12*r12]
C gamB: Geminal coupling between GB(r1) and GB(r2) == Exp[-GamB*r12*r12]
C
      implicit none
C Input variables
      integer NDIM
      parameter (NDIM=3)
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
      double precision xmass
C Output variables
      double precision ansKE
      double precision sval
C Local variables
      double precision DELTA
      parameter (DELTA = 1.0d-4)

      double precision Sx,Sy,Sz
      double precision Tx,Ty,Tz
      double precision ddx,ddy,ddz
C--------------Overlaps------------------------
      call G3ovlap_1D(I1,alp1,Amat1(1),
     1                I2,alp2,Amat2(1),
     1                I3,alp3,Amat3(1),
     2                L1,beta1,Bmat1(1),
     3                L2,beta2,Bmat2(1),
     4                L3,beta3,Bmat3(1),
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,Sx)

      call G3ovlap_1D(J1,alp1,Amat1(2),
     1                J2,alp2,Amat2(2),
     1                J3,alp3,Amat3(2),
     2                M1,beta1,Bmat1(2),
     3                M2,beta2,Bmat2(2),
     4                M3,beta3,Bmat3(2),
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,Sy)

      call G3ovlap_1D(K1,alp1,Amat1(3),
     1                K2,alp2,Amat2(3),
     1                K3,alp3,Amat3(3),
     2                N1,beta1,Bmat1(3),
     3                N2,beta2,Bmat2(3),
     4                N3,beta3,Bmat3(3),
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,Sz)
C--------------Laplacian--------------------
      call G3Lap_1D_g1Tpg2(I1,alp1,Amat1(1),
     1                I2,alp2,Amat2(1),
     1                I3,alp3,Amat3(1),
     2                L1,beta1,Bmat1(1),
     3                L2,beta2,Bmat2(1),
     4                L3,beta3,Bmat3(1),
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,DELTA,ddx)

      call G3Lap_1D_g1Tpg2(J1,alp1,Amat1(2),
     1                J2,alp2,Amat2(2),
     1                J3,alp3,Amat3(2),
     2                M1,beta1,Bmat1(2),
     3                M2,beta2,Bmat2(2),
     4                M3,beta3,Bmat3(2),
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,DELTA,ddy)

      call G3Lap_1D_g1Tpg2(K1,alp1,Amat1(3),
     1                K2,alp2,Amat2(3),
     1                K3,alp3,Amat3(3),
     2                N1,beta1,Bmat1(3),
     3                N2,beta2,Bmat2(3),
     4                N3,beta3,Bmat3(3),
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,DELTA,ddz)
C--------------kinetic energy--------------------
      Tx = -ddx*0.50d0/xmass
      Ty = -ddy*0.50d0/xmass
      Tz = -ddz*0.50d0/xmass

      ansKE = (Tx*Sy*Sz)
     1    + (Sx*Ty*Sz) 
     2    + (Sx*Sy*Tz) 

      sval = Sx * Sy * Sz
      END

C*****************************************************************
      SUBROUTINE G3Lap_1D_g1Tpg2(I1,alp1,Amat1,
     1                    I2,alp2,Amat2,
     1                    I3,alp3,Amat3,
     2                    L1,beta1,Bmat1,
     3                    L2,beta2,Bmat2,
     4                    L3,beta3,Bmat3,
     5                    gamA12,gamA13,gamA23,
     6                    gamB12,gamB13,gamB23,DELTA,ans)
C*****************************************************************
C Solve
C <GA(1)GA(2)GA(3)g(1,3)|d^2/dx3^2|GB(1)GB(2)GB(3)g(2,3)>
C
C GA,GB: are primitive Gaussian basis functions
C g(1,3) = geminal function between particle  1 and 3
C
C particles 1 and 2 are electrons and 3 is proton
C
      implicit none
C Input variables
      integer I1
      integer I2
      integer I3
      integer L1
      integer L2
      integer L3
      double precision alp1,Amat1
      double precision alp2,Amat2
      double precision alp3,Amat3
      double precision beta1,Bmat1
      double precision beta2,Bmat2
      double precision beta3,Bmat3
      double precision gamA12,gamA13,gamA23
      double precision gamB12,gamB13,gamB23
      double precision DELTA
C Output variables
      double precision ans
C Local variables
      double precision ZERO
      parameter (ZERO = 0.0d0)

      double precision S1
      double precision S2
      double precision S3
      double precision coef
      double precision AM2
      double precision BM2
      double precision BM3
C===================================================================================
C overlap integral
C S1(-Delta) = <GA(1)GA(2)GA(3)g(1,3)|GB(1)GB(2)GB(3-DELTA)g(2,3-DELTA)>
      coef = -1.0d0
      AM2  = Amat2 - (coef*Delta)
      BM2  = Bmat2 - (coef*Delta)
      BM3  = Bmat3 - (coef*Delta)

      call G3ovlap_1D(I1,alp1,Amat1,
     1                I2,alp2,AM2,
     1                I3,alp3,Amat3,
     2                L1,beta1,Bmat1,
     3                L2,beta2,BM2,
     4                L3,beta3,BM3,
     5                ZERO,gamA13,ZERO,
     6                ZERO,ZERO,gamB23,S1)
C----------------------------------------------------------------
C S2(0) = <GA(1)GA(2)GA(3)g(1,3)|GB(1)GB(2)GB(3)g(2,3)>
      coef = 0.0d0
      AM2  = Amat2 - (coef*Delta)
      BM2  = Bmat2 - (coef*Delta)
      BM3  = Bmat3 - (coef*Delta)

      call G3ovlap_1D(I1,alp1,Amat1,
     1                I2,alp2,AM2,
     1                I3,alp3,Amat3,
     2                L1,beta1,Bmat1,
     3                L2,beta2,BM2,
     4                L3,beta3,BM3,
     5                ZERO,gamA13,ZERO,
     6                ZERO,ZERO,gamB23,S2)
C----------------------------------------------------------------
C S3(Delta) = <GA(1)GA(2)GA(3)g(1,3)|GB(1)GB(2)GB(3+DELTA)g(2,3+DELTA)>
      coef = 1.0d0
      AM2  = Amat2 - (coef*Delta)
      BM2  = Bmat2 - (coef*Delta)
      BM3  = Bmat3 - (coef*Delta)

      call G3ovlap_1D(I1,alp1,Amat1,
     1                I2,alp2,AM2,
     1                I3,alp3,Amat3,
     2                L1,beta1,Bmat1,
     3                L2,beta2,BM2,
     4                L3,beta3,BM3,
     5                ZERO,gamA13,ZERO,
     6                ZERO,ZERO,gamB23,S3)
C----------------------------------------------------------------
C finite difference
      ans = (S1-(2.0d0*S2)+S3)/(DELTA*DELTA)


      END
