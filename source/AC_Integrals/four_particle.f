C23456789
C*****************************************************************
      SUBROUTINE G4ovlap_typ1(I1,J1,K1,alp1,Amat1,
     1                 I2,J2,K2,alp2,Amat2,
     1                 I4,J4,K4,alp4,Amat4,
     2                 L1,M1,N1,beta1,Bmat1,
     3                 L2,M2,N2,beta2,Bmat2,
     3                 L4,M4,N4,beta4,Bmat4,
     4                 gamA14,gamA24,
     4                 gamB14,gamB24,ans)
C*****************************************************************
C Gaussian Integral for Vee (PGIVee)
C
C Calculates the e-e integral using
C Gaussian type geminal (GTG)
C and Cartesian Gaussian functions
C
C  ans = <GA(r1)GA(r2)GA(r4)|
C         gamma(r1,r4)*gamma(r2,r4)
C        |GB(r1)GB(r2)GB(r4)>
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
C NOTE: The sign of the integral is positive
C       To use it for Coulomb attraction operator, multiply ans by -1.0d0
C
      implicit none
C Input variables
      integer NDIM
      parameter (NDIM=3)
      integer I1,J1,K1
      integer I2,J2,K2
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L4,M4,N4
      double precision alp1,Amat1(NDIM)
      double precision alp2,Amat2(NDIM)
      double precision alp4,Amat4(NDIM)
      double precision beta1,Bmat1(NDIM)
      double precision beta2,Bmat2(NDIM)
      double precision beta4,Bmat4(NDIM)
      double precision gamA14,gamA24
      double precision gamB14,gamB24
C Output variables
      double precision ans

C Local
      integer i
      double precision gam14
      double precision gam24
      double precision xk
      double precision coef
      double precision ovlap
      double precision beta5,Bmat5(NDIM)

C---------------------------------------------------------- 
C define new ket vector of GB(r4)
C GB'(r4) =  <GA(r1)GA(r2)|g(r1,r4)g(r2,r4)|GB(r1)GB(r2)>*GB(r4)
C and store it in beta5,Bmat5
C
C     --x,y,z component--
      gam14 = gamA14 + gamB14
      gam24 = gamA24 + gamB24
      coef = 1.0d0
      do i=1,NDIM
         call basis_update(alp1,beta1,alp2,beta2,beta4,
     1                  Amat1(i),Bmat1(i),Amat2(i),Bmat2(i),Bmat4(i),
     1                  gam14,gam24,xk,beta5,Bmat5(i))
         coef = coef * xk
      end do
C----------------------------------------------------------------------     
C  ans = coef * <GA(r4)|GB'(r4)>
      call gfovlap(I4,J4,K4,alp4,Amat4,
     2             L4,M4,N4,beta5,Bmat5,ovlap)

      ans = coef * ovlap



      END



C*****************************************************************
      SUBROUTINE G4vee_typ1(I1,J1,K1,alp1,Amat1,
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
C Gaussian Integral for Vee (PGIVee)
C
C Calculates the e-e integral using
C Gaussian type geminal (GTG)
C and Cartesian Gaussian functions
C
C  ans = <GA(r1)GA(r2)GA(r3)GA(r4)|
C         gamma(r1,r4)*gamma(r2,r4)/(r3-r4)
C        |GB(r1)GB(r2)GB(r3)GB(r4)>
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
C NOTE: The sign of the integral is positive
C       To use it for Coulomb attraction operator, multiply ans by -1.0d0
C
      implicit none
C Input variables
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

C Local
      integer i
      double precision gam14
      double precision gam24
      double precision xk
      double precision coef
      double precision vee
      double precision beta5,Bmat5(NDIM)

C---------------------------------------------------------- 
C define new ket vector of GB(r4)
C GB'(r4) =  <GA(r1)GA(r2)|g(r1,r4)g(r2,r4)|GB(r1)GB(r2)>*GB(r4)
C and store it in beta5,Bmat5
C
C     --x,y,z component--
      gam14 = gamA14 + gamB14
      gam24 = gamA24 + gamB24
      coef = 1.0d0
      do i=1,NDIM
         call basis_update(alp1,beta1,alp2,beta2,beta4,
     1                  Amat1(i),Bmat1(i),Amat2(i),Bmat2(i),Bmat4(i),
     1                  gam14,gam24,xk,beta5,Bmat5(i))
         coef = coef * xk
      end do
C----------------------------------------------------------------------     
C  ans = coef * <GA(r3)GA(r4)|1/(r3-r4)|GB(r3)GB'(r4)>
      call gfvee(I3,J3,K3,alp3,Amat3,
     1           I4,J4,K4,alp4,Amat4,
     2           L3,M3,N3,beta3,Bmat3,
     3           L4,M4,N4,beta5,Bmat5,
     4           vee)
      ans = coef * vee

      END

C*****************************************************************
      SUBROUTINE G4vee_typ2(I1,J1,K1,alp1,Amat1,
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
C Calculates the e-e integral using
C Gaussian type geminal (GTG)
C and Cartesian Gaussian functions
C
C  ans = <GA(r1)GA(r2)GA(r3)GA(r4)|
C         gamma(r1,r4)*gamma(r3,r4)/(r1-r2)
C        |GB(r1)GB(r2)GB(r3)GB(r4)>
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
C NOTE: The sign of the integral is positive
C       To use it for Coulomb attraction operator, multiply ans by -1.0d0
C
      implicit none
C Input variables
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

C Local
      integer i
      double precision gam14
      double precision gam24
      double precision gam34
      double precision xk
      double precision coef
      double precision vee
      double precision beta5,Bmat5(NDIM)

C---------------------------------------------------------- 
C define new ket vector of GB(r4)
C GB'(r1) =  <GA(r3)GA(r4)|g(r1,r4)g(r3,r4)|GB(r3)GB(r4)>*GB(r1)
C and store it in beta5,Bmat5
C
C     --x,y,z component--
      gam14 = gamA14 + gamB14
      gam34 = gamA34 + gamB34
      coef = 1.0d0
      do i=1,NDIM
         call basis_update(alp3,beta3,alp4,beta4,beta1,
     1                  Amat3(i),Bmat3(i),Amat4(i),Bmat4(i),Bmat1(i),
     1                  gam14,gam34,xk,beta5,Bmat5(i))
         coef = coef * xk
      end do
C----------------------------------------------------------------------     
C  ans = coef * <GA(r1)GA(r2)|1/(r1-r2)|GB'(r1)GB(r2)>
      call gfvee(I1,J1,K1,alp1,Amat1,
     1           I2,J2,K2,alp2,Amat2,
     2           L1,M1,N1,beta5,Bmat5,
     3           L2,M2,N2,beta2,Bmat2,
     4           vee)
      ans = coef * vee

      END


C******************************************************************
      SUBROUTINE basis_update(a1,b1,a2,b2,b3,AX1,BX1,AX2,BX2,BX3,
     1                        gam1,gam2,coef,p4,PX4)
C******************************************************************
C  phi_1A(x1) = Exp[-a1*(x1-AX1)^2]
C  phi_1B(x1) = Exp[-b1*(x1-BX1)^2]
C
C  phi_2A(x2) = Exp[-a2*(x2-AX1)^2]
C  phi_2B(x2) = Exp[-b2*(x2-BX1)^2]
C
C  g1 = Exp[-gam1*(x1-x4)^2]
C  g2 = Exp[-gam2*(x2-x4)^2]
C
C  ans = coef * Exp[-p4*(x4-PX4)^2]
C      = <phi_1A|g1|phi_1B> <phi_2A|g2|phi_2B>
C
      implicit none
      double precision a1,b1,AX1,BX1
      double precision a2,b2,AX2,BX2
      double precision b3,BX3
      double precision gam1,gam2
      double precision coef,p4,PX4

      double precision xk1,xk2,xk3,xk4
      double precision p1,PX1
      double precision p2,PX2
      double precision p3,PX3
      double precision xjunk1,xjunk2,xjunk3,xjunk4

C     --integrate over g(x1,x4)*phi(x1)*phi'(x1)--
C      = xk1 Exp[-p1*(x4-PX1)^2]
      call three_term_gauss(a1,b1,gam1,AX1,BX1,xk1,p1,PX1)

C     --integrate over g(x2,x4)*phi(x2)*phi'(x2)--
C     = xk2 Exp[-p2*(x4-PX2)^2]
      call three_term_gauss(a2,b2,gam2,AX2,BX2,xk2,p2,PX2)

C     --product of p1 and p2 gaussian--
C     = xk3 Exp[-p3(x4-PX3)^2]
      call gauss_g2_prod(p1,p2,PX1,PX2,p3,xjunk3,PX3,xjunk2,xk3,xjunk1)

C     --answer--
C     --product p3 * with original basis--
C     Exp[-p3*(x4-PX3)^2] Exp[-b3*(x4-BX3)^2]
C     = xk4 Exp[-p4*(x4-PX4)^2]
C
      call gauss_g2_prod(p3,b3,PX3,BX3,p4,xjunk1,PX4,xjunk2,xk4,xjunk3)
      coef = xk1 * xk2 * xk3 * xk4


C     write(*,*) 'a1,AX1',a1,AX1
C     write(*,*) 'b1,BX1',b1,BX1
C     write(*,*) 'a2,AX2',a2,AX2
C     write(*,*) 'b2,BX2',b2,BX2
C     write(*,*) 'b2,BX2',b2,BX2
C     write(*,*) 'b3,BX3',b3,BX3
C     write(*,*) 'gam1',gam1
C     write(*,*) 'gam2',gam2
C     write(*,*) '----ANSWER----'
C     write(*,*) 'xk1,p1',xk1,p1
C     write(*,*) 'xk2,p2',xk2,p2
C     write(*,*) 'xk3,p3',xk3,p3
C     write(*,*) 'coef',coef
C     write(*,*) 'p4,PX4',p4,PX4

C     STOP '--AFTER BASIS UPDATE--'

      
      END


C******************************************************************
      SUBROUTINE three_term_gauss(a1,a2,a3,AX1,AX2,xk4,a4,AX4)
C******************************************************************
C ans = <p1(x1)|g(x1,x4)|p2(x1)>
C ans = xk4 * Exp[-a4(x4-AX4)^2]
C
      implicit none
      double precision a1,a2,a3
      double precision AX1,AX2
      double precision xk4,a4,AX4

      double precision p1,q1,PX1,QX1,xk1
      double precision p2,q2,xk2
      double precision xjunk1

C     --combine AX1 and AX2--
C     g1*g2 = k1*Exp[-p(x-PX)^2] = k1*p1
      call gauss_g2_prod(a1,a2,AX1,AX2,p1,q1,PX1,QX1,xk1,xjunk1)

C     --combine PX1 and AX3--
C     p1*g3 = Exp[-q*(AX3-PX1)^2] Exp[-p2*(x-PX2)^]
C     call gauss_g2_prod(p1,a3,PX1,AX3,p2,q,xjunk4,xjunk1,xjunk2,xjunk3)
      p2 = p1 + a3
      q2 = p1*a3/p2

C     --integrate over x--
      xk2 = dsqrt(dacos(-1.0d0)/p2)

C     --final pre-factor--
      xk4 = xk1 * xk2
      a4  = q2
      AX4 = PX1
      
      END 

C*************************************************************
      SUBROUTINE gauss_g2_prod(a,b,AX,BX,p,q,PX,QX,xKab,xint)
C*************************************************************
C Calculate 1-D overlap density
C Exp[-a(x-AX)^2] Exp[-b(x-BX) = Exp[-q*QX*QX] Exp[-p(x-PX)^2]
C                              = xKab          Exp[-p(x-PX)^2]
C integration = xKab*sqrt(PI/p)
C
      implicit none
C input variables
      double precision a,b,AX,BX
C output variables
      double precision p,q
      double precision PX,QX
      double precision xKab
      double precision xint
C Gaussian product
      p    = a + b
      if( p .ne. 0.0d0) then
         PX   = ((a*AX)+(b*BX))/p
      else
         PX   = 0.0d0
      end if

      if( p .ne. 0.0d0) then
         q    = (a*b)/p
      else
         q    = 0.0d0
      end if

      Qx   = Ax-Bx
      xKab = exp(-q*Qx*Qx)


      if( p .ne. 0.0d0) then
         xint = dsqrt(dacos(-1.0d0)/p) * xKab
      else
         xint = 0.0d0
      end if

      END

C*************************************************************
      SUBROUTINE gauss_gn_prod(N,a,AX,q,QX,p,PX)
C*************************************************************
      implicit none
      integer N
      double precision a(N)
      double precision AX(N)
      double precision q(N,N)
      double precision QX(N,N)
      double precision p
      double precision PX

      integer i,j

C     --calculate alpha for Gaussian--
      p = 0.0d0
      do i=1,N
         p = p + a(i)
      end do

C     --calculate center for Gaussian--
      PX = 0.0d0
      do i=1,N
         PX = PX + (a(i)*A(i))
      end do
      PX = PX/p

C     --calculate alpha for the const--
      do i=1,N
      do j=1,N
         q(i,j) = a(i)*a(j)/p
         q(i,i) = 0.0d0
         q(j,j) = 0.0d0
      end do
      end do

C     --calculate center for the const--
      do i=1,N
      do j=1,N
         QX(i,j) = (AX(i)-AX(j))**2
         QX(i,i) = 0.0d0
         QX(j,j) = 0.0d0
      end do 
      end do 
      
      END

