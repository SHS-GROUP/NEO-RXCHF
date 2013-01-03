C*****************************************************************
      SUBROUTINE G3ke_gnT1g2(I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     4                gamA12,gamA13,gamA23,
     4                gamB12,gamB13,gamB23,xmass,ansKE)
C*****************************************************************
C Geminal 3-particle KE  matrix element
C <GA(1)GA(2)GA(3)g(1,2)g(1,3)g(2,3)|T(1)|GB(1)GB(2)GB(3)g(2,3)>
C
C The Laplacian is calculated with respect to the coordinates of 
C particle 1
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
      call G3Lap_1D_gnT1g2(I1,alp1,Amat1(1),
     1                I2,alp2,Amat2(1),
     1                I3,alp3,Amat3(1),
     2                L1,beta1,Bmat1(1),
     3                L2,beta2,Bmat2(1),
     4                L3,beta3,Bmat3(1),
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,ddx)

      call G3Lap_1D_gnT1g2(J1,alp1,Amat1(2),
     1                J2,alp2,Amat2(2),
     1                J3,alp3,Amat3(2),
     2                M1,beta1,Bmat1(2),
     3                M2,beta2,Bmat2(2),
     4                M3,beta3,Bmat3(2),
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,ddy)

      call G3Lap_1D_gnT1g2(K1,alp1,Amat1(3),
     1                K2,alp2,Amat2(3),
     1                K3,alp3,Amat3(3),
     2                N1,beta1,Bmat1(3),
     3                N2,beta2,Bmat2(3),
     4                N3,beta3,Bmat3(3),
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,ddz)
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
      SUBROUTINE G3Lap_1D_gnT1g2(I1,alp1,Amat1,
     1                    I2,alp2,Amat2,
     1                    I3,alp3,Amat3,
     2                    L1,beta1,Bmat1,
     3                    L2,beta2,Bmat2,
     4                    L3,beta3,Bmat3,
     5                    gamA12,gamA13,gamA23,
     6                    gamB12,gamB13,gamB23,ans)
C*****************************************************************
C Calculate 1D Laplacian element with respect to x1:
C <GA(1)GA(2)GA(3)g(1,2)g(1,3)g(2,3)|D^2/Dx1^2|GB(1)GB(2)GB(3)g(2,3)>
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
C Output variables
      double precision ans
C Local variables
      double precision ZERO
      parameter (ZERO = 0.0d0)

      double precision c1,c2
      double precision d1,d2,d3
      double precision p1,p2,p3
      double precision q1,q2
      double precision t1,t2
      double precision W1,W2,W3
C--------------Initialize-----------------------
      d1 = 0.0d0
      d2 = 0.0d0
      d3 = 0.0d0
      c1 = 0.0d0
      c2 = 0.0d0
      p1 = 0.0d0
      p2 = 0.0d0
      p3 = 0.0d0
      q1 = 0.0d0
      q2 = 0.0d0
      t1 = 0.0d0
      t2 = 0.0d0
C--------------Calculate ceoff-----------------------
C get the expansion coeff
      if(L1-1 .ge. 0) d1 = dble(L1*(L1-1))
      d2 = -2.0d0*beta1*dble((2*L1)+1)
      d3 = 4.0d0*beta1*beta1
      c1 = dble(L1)
      c2 = -2.0d0*beta1
C--------------Calculate W1-----------------------
      if(L1-2 .ge. 0) then
C        call gfovlap_1D(I1,alp1,Amat1,
C    2                 L1-2,beta1,Bmat1,
C    4                 p1)

      call G3ovlap_1D(I1,alp1,Amat1,
     1                I2,alp2,Amat2,
     1                I3,alp3,Amat3,
     2                L1-2,beta1,Bmat1,
     3                L2,beta2,Bmat2,
     4                L3,beta3,Bmat3,
     5                gamA12,gamA13,gamA23,
     6                ZERO,ZERO,gamB23,p1)
      end if


C     call gfovlap_1D(I1,alp1,Amat1,
C    2                 L1,beta1,Bmat1,
C    4                 p2)

      call G3ovlap_1D(I1,alp1,Amat1,
     1                I2,alp2,Amat2,
     1                I3,alp3,Amat3,
     2                L1,beta1,Bmat1,
     3                L2,beta2,Bmat2,
     4                L3,beta3,Bmat3,
     5                gamA12,gamA13,gamA23,
     6                ZERO,ZERO,gamB23,p2)


C     call gfovlap_1D(I1,alp1,Amat1,
C    2                 L1+2,beta1,Bmat1,
C    4                 p3)


      call G3ovlap_1D(I1,alp1,Amat1,
     1                I2,alp2,Amat2,
     1                I3,alp3,Amat3,
     2                L1+2,beta1,Bmat1,
     3                L2,beta2,Bmat2,
     4                L3,beta3,Bmat3,
     5                gamA12,gamA13,gamA23,
     6                ZERO,ZERO,gamB23,p3)


      W1=(d1*p1)+(d2*p2)+(d3*p3)

      ans = W1

      END
