C23456789
C---------------------------------------------------------------------
C START_OF_THREE_PARTICLE_INTEGRALS
C---------------------------------------------------------------------
C Routines to calculate all the matrix elements invovling
C Cartesian Gaussian functions.
C
C References:
C [1]. "Molecular integrals over Gaussian-type geminal basis functions"
C       by Joakim Persson and Peter R. Taylor
C       Theoretical Chemistry Accounts, Vol. 97, 240-250, 1997
C
C*****************************************************************
      SUBROUTINE G3ovlap_1D(I1,alp1,Amat1,
     1                      I2,alp2,Amat2,
     1                      I3,alp3,Amat3,
     2                      L1,beta1,Bmat1,
     3                      L2,beta2,Bmat2,
     4                      L3,beta3,Bmat3,
     5                      gamA12,gamA13,gamA23,
     6                      gamB12,gamB13,gamB23,ans)
C*****************************************************************
C Primitive Gaussian Integral for overlap in 1D (PGIovlap)
C
C Calculates the e-e integral using
C Gaussian type geminal (GTG)
C and Cartesian Gaussian functions
C
C  ans = <GA(r1)GA(r2)|Exp[-gamma*r12^2]|GB(r1)GB(r2)>
C
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
      double precision gam12,gam13,gam23
      double precision p1,q1
      double precision p2,q2
      double precision p3,q3
      double precision xKab1
      double precision xKab2
      double precision xKab3
      double precision QX1,PX1
      double precision QX2,PX2
      double precision QX3,PX3

      integer it1,iu1,iv1
      integer it2,iu2,iv2
      integer it3,iu3,iv3
      double precision E1,E2,E3,xsum,herm
      double precision ex1(0:I1+L1)
      double precision ex2(0:I2+L2)
      double precision ex3(0:I3+L3)

C------------------calc gamma--------------------
      gam12=gamA12+gamB12
      gam13=gamA13+gamB13
      gam23=gamA23+gamB23
C------------------overlap--------------------
C get overlap and hermite expansion coefficients
C xyz for electron 1
C
      call omega1D(I1,L1,alp1,beta1,Amat1,Bmat1,p1,q1,PX1,
     1             QX1,xKab1,ex1)

C xyz for electron 2
      call omega1D(I2,L2,alp2,beta2,Amat2,Bmat2,p2,q2,PX2,
     1             QX2,xKab2,ex2)

C xyz for electron 3
      call omega1D(I3,L3,alp3,beta3,Amat3,Bmat3,p3,q3,PX3,
     1             QX3,xKab3,ex3)
C------------------hermite-loop--------------------
C loop over hermite integrals
      xsum = 0.0d0
      do it1=0,I1+L1
         do it2=0,I2+L2
            do it3=0,I3+L3
               call G3hint0_1D(it1,it2,it3,p1,p2,p3,
     1                         gam12,gam13,gam23,
     2                         PX1,PX2,PX3,herm)
               E1  =ex1(it1)
               E2  =ex2(it2)
               E3  =ex3(it3)
               xsum=xsum+(E1*E2*E3*herm)
            end do
         end do
      end do
      ans = xsum

      END

C*****************************************************************
      SUBROUTINE G3vee(I1,J1,K1,alp1,Amat1,
     1                 I2,J2,K2,alp2,Amat2,
     1                 I3,J3,K3,alp3,Amat3,
     2                 L1,M1,N1,beta1,Bmat1,
     3                 L2,M2,N2,beta2,Bmat2,
     3                 L3,M3,N3,beta3,Bmat3,
     4                 gamA12,gamA13,gamA23,
     4                 gamB12,gamB13,gamB23,ans)
C*****************************************************************
C Primitive Gaussian Integral for Vee (PGIVee)
C
C Calculates the e-e integral using
C Gaussian type geminal (GTG)
C and Cartesian Gaussian functions
C
C  ans = <GA(r1)GA(r2)|Exp[-gamma*r12^2]/r12|GB(r1)GB(r2)>
C
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
      double precision gam12,gam13,gam23
      double precision p1,p2,p3,q1,q2,q3
      double precision xKab1(NDIM),xKab2(NDIM),xKab3(NDIM)
      double precision Qmat1(NDIM),Pmat1(NDIM)
      double precision Qmat2(NDIM),Pmat2(NDIM)
      double precision Qmat3(NDIM),Pmat3(NDIM)

      integer it1,iu1,iv1
      integer it2,iu2,iv2
      integer it3,iu3,iv3
      double precision E1,E2,E3,xsum,herm
      double precision ex1(0:I1+L1),ey1(0:J1+M1),ez1(0:K1+N1)
      double precision ex2(0:I2+L2),ey2(0:J2+M2),ez2(0:K2+N2)
      double precision ex3(0:I3+L3),ey3(0:J3+M3),ez3(0:K3+N3)

C------------------calc gamma--------------------
      gam12 = gamA12 + gamB12
      gam13 = gamA13 + gamB13
      gam23 = gamA23 + gamB23
C------------------overlap--------------------
C get overlap and hermite expansion coefficients
C xyz for electron 1
C
      call omega1D(I1,L1,alp1,beta1,Amat1(1),Bmat1(1),p1,q1,Pmat1(1),
     1             Qmat1(1),xKab1(1),ex1)
      call omega1D(J1,M1,alp1,beta1,Amat1(2),Bmat1(2),p1,q1,Pmat1(2),
     1             Qmat1(2),xKab1(2),ey1)
      call omega1D(K1,N1,alp1,beta1,Amat1(3),Bmat1(3),p1,q1,Pmat1(3),
     1             Qmat1(3),xKab1(3),ez1)

C xyz for electron 2
      call omega1D(I2,L2,alp2,beta2,Amat2(1),Bmat2(1),p2,q2,Pmat2(1),
     1             Qmat2(1),xKab2(1),ex2)
      call omega1D(J2,M2,alp2,beta2,Amat2(2),Bmat2(2),p2,q2,Pmat2(2),
     1             Qmat2(2),xKab2(2),ey2)
      call omega1D(K2,N2,alp2,beta2,Amat2(3),Bmat2(3),p2,q2,Pmat2(3),
     1             Qmat2(3),xKab2(3),ez2)

C xyz for electron 3
      call omega1D(I3,L3,alp3,beta3,Amat3(1),Bmat3(1),p3,q3,Pmat3(1),
     1             Qmat3(1),xKab3(1),ex3)
      call omega1D(J3,M3,alp3,beta3,Amat3(2),Bmat3(2),p3,q3,Pmat3(2),
     1             Qmat3(2),xKab3(2),ey3)
      call omega1D(K3,N3,alp3,beta3,Amat3(3),Bmat3(3),p3,q3,Pmat3(3),
     1             Qmat3(3),xKab3(3),ez3)
C------------------hermite-loop--------------------
C loop over hermite integrals
      xsum = 0.0d0
      do it1=0,I1+L1
      do iu1=0,J1+M1
      do iv1=0,K1+N1
         do it2=0,I2+L2
         do iu2=0,J2+M2
         do iv2=0,K2+N2
            do it3=0,I3+L3
            do iu3=0,J3+M3
            do iv3=0,K3+N3
         call G3hint2(it1,iu1,iv1,it2,iu2,iv2,it3,iu3,iv3,
     1                p1,p2,p3,gam12,gam13,gam23,Pmat1,Pmat2,Pmat3,
     2                herm)
         E1  =ex1(it1)*ey1(iu1)*ez1(iv1)
         E2  =ex2(it2)*ey2(iu2)*ez2(iv2)
         E3  =ex3(it3)*ey3(iu3)*ez3(iv3)
         xsum=xsum+(E1*E2*E3*herm)
            end do
            end do
            end do
         end do
         end do
         end do
      end do
      end do
      end do
      ans = xsum

      END

C*****************************************************************
      SUBROUTINE G3vec(I1,J1,K1,alp1,Amat1,
     1                 I2,J2,K2,alp2,Amat2,
     1                 I3,J3,K3,alp3,Amat3,
     2                 L1,M1,N1,beta1,Bmat1,
     3                 L2,M2,N2,beta2,Bmat2,
     3                 L3,M3,N3,beta3,Bmat3,
     4                 gamA12,gamA13,gamA23,
     4                 gamB12,gamB13,gamB23,Cmat,ans)
C*****************************************************************
C Primitive Gaussian Integral for Vee (PGIVee)
C
C Calculates the e-e integral using
C Gaussian type geminal (GTG)
C and Cartesian Gaussian functions
C
C  ans = <GA(r1)GA(r2)|Exp[-gamma*r12^2]/r12|GB(r1)GB(r2)>
C
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
      double precision Cmat(NDIM)
C Output variables
      double precision ans
C Local variables
      double precision gam12,gam13,gam23
      double precision p1,p2,p3,q1,q2,q3
      double precision xKab1(NDIM),xKab2(NDIM),xKab3(NDIM)
      double precision Qmat1(NDIM),Pmat1(NDIM)
      double precision Qmat2(NDIM),Pmat2(NDIM)
      double precision Qmat3(NDIM),Pmat3(NDIM)

      integer it1,iu1,iv1
      integer it2,iu2,iv2
      integer it3,iu3,iv3
      double precision E1,E2,E3,xsum,herm
      double precision ex1(0:I1+L1),ey1(0:J1+M1),ez1(0:K1+N1)
      double precision ex2(0:I2+L2),ey2(0:J2+M2),ez2(0:K2+N2)
      double precision ex3(0:I3+L3),ey3(0:J3+M3),ez3(0:K3+N3)

C------------------calc gamma--------------------
      gam12 = gamA12 + gamB12
      gam13 = gamA13 + gamB13
      gam23 = gamA23 + gamB23
C------------------overlap--------------------
C get overlap and hermite expansion coefficients
C xyz for electron 1
C
      call omega1D(I1,L1,alp1,beta1,Amat1(1),Bmat1(1),p1,q1,Pmat1(1),
     1             Qmat1(1),xKab1(1),ex1)
      call omega1D(J1,M1,alp1,beta1,Amat1(2),Bmat1(2),p1,q1,Pmat1(2),
     1             Qmat1(2),xKab1(2),ey1)
      call omega1D(K1,N1,alp1,beta1,Amat1(3),Bmat1(3),p1,q1,Pmat1(3),
     1             Qmat1(3),xKab1(3),ez1)

C xyz for electron 2
      call omega1D(I2,L2,alp2,beta2,Amat2(1),Bmat2(1),p2,q2,Pmat2(1),
     1             Qmat2(1),xKab2(1),ex2)
      call omega1D(J2,M2,alp2,beta2,Amat2(2),Bmat2(2),p2,q2,Pmat2(2),
     1             Qmat2(2),xKab2(2),ey2)
      call omega1D(K2,N2,alp2,beta2,Amat2(3),Bmat2(3),p2,q2,Pmat2(3),
     1             Qmat2(3),xKab2(3),ez2)

C xyz for electron 3
      call omega1D(I3,L3,alp3,beta3,Amat3(1),Bmat3(1),p3,q3,Pmat3(1),
     1             Qmat3(1),xKab3(1),ex3)
      call omega1D(J3,M3,alp3,beta3,Amat3(2),Bmat3(2),p3,q3,Pmat3(2),
     1             Qmat3(2),xKab3(2),ey3)
      call omega1D(K3,N3,alp3,beta3,Amat3(3),Bmat3(3),p3,q3,Pmat3(3),
     1             Qmat3(3),xKab3(3),ez3)
C------------------hermite-loop--------------------
C loop over hermite integrals
      xsum = 0.0d0
      do it1=0,I1+L1
      do iu1=0,J1+M1
      do iv1=0,K1+N1
         do it2=0,I2+L2
         do iu2=0,J2+M2
         do iv2=0,K2+N2
            do it3=0,I3+L3
            do iu3=0,J3+M3
            do iv3=0,K3+N3
         call G3hint1(it1,iu1,iv1,it2,iu2,iv2,it3,iu3,iv3,
     1                p1,p2,p3,gam12,gam13,gam23,Pmat1,Pmat2,Pmat3,Cmat,
     2                herm)
         E1  =ex1(it1)*ey1(iu1)*ez1(iv1)
         E2  =ex2(it2)*ey2(iu2)*ez2(iv2)
         E3  =ex3(it3)*ey3(iu3)*ez3(iv3)
         xsum=xsum+(E1*E2*E3*herm)
            end do
            end do
            end do
         end do
         end do
         end do
      end do
      end do
      end do
      ans = xsum

      END

C************************************************************************
      SUBROUTINE G3hint1(it1,iu1,iv1,it2,iu2,iv2,it3,iu3,iv3,
     1                   p1,p2,p3,gam12,gam13,gam23,Pmat1,Pmat2,Pmat3,
     1                   Cmat,ans)
C************************************************************************
C Three particle one dimensional Hermite Gaussian integral 
C
C Input varaibles
      implicit none
      integer NDIM
      parameter (NDIM=3)
      integer it1,iu1,iv1
      integer it2,iu2,iv2
      integer it3,iu3,iv3
      double precision p1
      double precision p2
      double precision p3
      double precision gam12
      double precision gam13
      double precision gam23
      double precision Pmat1(NDIM)
      double precision Pmat2(NDIM)
      double precision Pmat3(NDIM)
      double precision Cmat(NDIM)
C Output varaibles
      double precision ans
C FunCall
      double precision binoc
C SubCall
C     gdel_x
C Local varaibles
      double precision PI
      parameter (PI=3.141592650d0)
      integer i
      integer i1,i2,i3
      integer j1,j2,j3
      integer k1,k2,k3
      double precision b1,b2,b3
      double precision covlp
      double precision c01D,c03D,c12,c13,c23
      double precision eta11,eta13,eta22,eta23

      double precision xsum,val
      double precision cov
      double precision d11,d12,d13
      double precision d21,d22,d23
      double precision d31,d32,d33
      double precision beta1,beta2,beta3

      double precision u,u1,u2,u3,v
      double precision xk1,xk2,xk3
      double precision t1,Tmat1(NDIM)
      double precision t2,Tmat2(NDIM)
      double precision s1,Smat1(NDIM)
      double precision s2,Smat2(NDIM)
      double precision gam12_new
      double precision h0

      double precision GX,GY,GZ,Rnlm
      double precision bt1,bt2,bt3
      double precision bu1,bu2,bu3
      double precision bv1,bv2,bv3
      double precision f1,f2,f3
      double precision g1,g2,g3
      double precision sign

      double precision tmp1,tmp2
C-----------Integrate over x3---------------
      call G3remove(p3,gam13,gam23,c01D,c03D,c12,c13,c23)
C
C Integration over r3 gives one coefficient 
C c03D
C and three Gaussian functions:
C Exp[-c12(r1-r2)**2]
C Exp[-c13(r1-P3)**2]
C Exp[-c23(r2-P3)**2]
C
C update gam12 to:
C gam12 = gam12 + c12
C and carry out the following Gauss product
C
C Exp[-p1(r1-P1)**2]*Exp[-c13(r1-P3)**2] = K1 * Exp[-s1(r1-S1)]
C Exp[-p2(r2-P2)**2]*Exp[-c23(r2-P3)**2] = K2 * Exp[-s2(r2-S2)]
C
C (r1C^-1) integration of
C Integrate[dx1,dx2] (r1C^-1) Exp[-s1(r1-S1)]*Exp[-s2(r2-S2)]*Exp[-gam12*r12**2]  = cov * K3 *F(v*RCQ)
C
C The final expression is: c03D * cov * K1 * K2 * K3
C
C cov = (PI/(s2+gam12))**(3/2) (2*PI/v)
C K1  = Exp[-t1*T1**2]      = Exp[-t1(P1-P3)**2]
C K2  = Exp[-t2*T2**2]      = Exp[-t2(P2-P3)**2]
C K3  = Exp[-u*(S1-S2)**2] = Exp[-u(e1*P1 + e2*P2 + e3*P3)**2]
C u   = s1*s2*gam12/(s1*s2 + s1*gam12 + s2*gam12)
C v   = (s1*s2 + s1*gam12 + s2*gam12)/(s2 + gam12)
C Q   = f1*S1 + f2*S2
C
C f1  = (s1*s2 + s1*gam12)/(s1*s2 + s1*gam12 + s2*gam12)
C f2  = s2*gam12)/(s1*s2 + s1*gam12 + s2*gam12)
C
C S1 = (p1/s1)*P1 + (c13/s1)*P3
C S2 = (p2/s2)*P2 + (c23/s2)*P3
C S1-S2 = (p1/s1)*P1 -(p2/s2)*P2 + [ (c13/s1)-(c23/s2) ]*P3
C
C TX1 = PX1-PX3
C TX2 = PX2-PX3
C
C Q   = f1*S1 + f2*S2
C Q   = f1*(p1/s1)*P1 + f1*(c13/s1)*P3 +  f2*(p2/s2)*P2 + f2*(c23/s2)*P3
C Q   = f1*(p1/s1)*P1 + f2*(p2/s2)*P2  + [f1*(c13/s1)   + f2*(c23/s2)]*P3
C Q   = g1*P1 + g2*P2 + g3*P3
C
C
C
C
      gam12_new = gam12 + c12

      do i=1,NDIM
       call G3gprod_1D(p1,c13,Pmat3(i),Pmat3(i),t1,s1,Tmat1(i),Smat1(i))
       call G3gprod_1D(p2,c23,Pmat2(i),Pmat3(i),t2,s2,Tmat2(i),Smat2(i))
      end do

      v   = ((s1*s2)+(s1*gam12_new)+(s2*gam12_new))/(s2+gam12_new)
      tmp1= sqrt(PI/(s2+gam12_new))**3
      cov = tmp1*2.0d0*PI/v
      u   = s1*s2*gam12_new/( (s1*s2)+(s1*gam12_new)+(s2*gam12_new) )

      tmp2= (s1*s2)+(s1*gam12_new)+(s2*gam12_new)
      f1  = ((s1*s2)+(s1*gam12_new))/tmp2
      f2  = s2*gam12_new/tmp2
      g1  = f1*(p1/s1)
      g2  = f2*(p2/s2)
      g3  = (f1*c13/s1) + (f2*c23/s2)
C
C Using geneal expression:
C K1 = Exp[-beta1*(d11*PX1+d12*PX2+d13*PX3)**2]
C K2 = Exp[-beta2*(d21*PX1+d22*PX2+d23*PX3)**2]
C K3 = Exp[-beta3*(d31*PX1+d32*PX2+d33*PX3)**2]
C
C d11= 1      d12=  0       d13= -1
C d21= 0      d22=  1       d23= -1
C d31= e1     d22=  e2      d23= e3
C-----------prepare for deriv--------------------
      beta1=t1    ; beta2=t2     ; beta3=u
      d11  =1.0d0 ; d12  =0.0d0  ; d13  =-1.0d0
      d21  =0.0d0 ; d22  =1.0d0  ; d23  =-1.0d0
      d31  =p1/s1 ; d32  =-p2/s2 ; d33  =(c13/s1)-(c23/s2)
C-----------Derivative wrt PX1 PX2 and PX3--------
      xsum = 0.0d0
      do i1=0,it1
        bt1=binoc(it1,i1)
        do i2=0,it2
          bt2=binoc(it2,i2)
          do i3=0,it3
            bt3=binoc(it3,i3)
            call G3KKK_1D(it1-i1,it2-i2,it3-i3,
     1                    beta1,beta2,beta3,
     2                    d11,d12,d13,
     3                    d21,d22,d23,
     4                    d31,d32,d33,
     5                    Pmat1(1),Pmat2(1),Pmat3(1),
     6                    GX)
            do j1=0,iu1
              bu1=binoc(iu1,j1)
              do j2=0,iu2
                bu2=binoc(iu2,j2)
                do j3=0,iu3
                  bu3=binoc(iu3,j3)
                  call G3KKK_1D(iu1-j1,iu2-j2,iu3-j3,
     1                    beta1,beta2,beta3,
     2                    d11,d12,d13,
     3                    d21,d22,d23,
     4                    d31,d32,d33,
     5                    Pmat1(2),Pmat2(2),Pmat3(2),
     6                    GY)
                  do k1=0,iv1
                    bv1=binoc(iv1,k1)
                    do k2=0,iv2
                      bv2=binoc(iv2,k2)
                      do k3=0,iv3
                        bv3=binoc(iv3,k3)
                        call G3KKK_1D(iv1-k1,iv2-k2,iv3-k3,
     1                    beta1,beta2,beta3,
     2                    d11,d12,d13,
     3                    d21,d22,d23,
     4                    d31,d32,d33,
     5                    Pmat1(3),Pmat2(3),Pmat3(3),
     6                    GZ)

                        sign = (-1.0d0)**(i1+j1+k1+i2+j2+k2+i3+j3+k3)
                        call G3Rnlm(NDIM,i1,j1,k1,i2,j2,k2,i3,j3,k3,
     1                    v,g1,g2,g3,
     2                    Pmat1,Pmat2,Pmat3,Cmat,Rnlm)

                          val=bt1*bt2*bt3*
     1                        bu1*bu2*bu3*
     2                        bv1*bv2*bv3*
     3                        GX*GX*GZ*
     4                        Rnlm * sign
                          xsum=xsum+val
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      ans = cov*c03D*val


      END

C************************************************************************
      SUBROUTINE G3hint2(it1,iu1,iv1,it2,iu2,iv2,it3,iu3,iv3,
     1                   p1,p2,p3,gam12,gam13,gam23,Pmat1,Pmat2,Pmat3,
     1                   ans)
C************************************************************************
C Three particle one dimensional Hermite Gaussian integral 
C
C Input varaibles
      implicit none
      integer NDIM
      parameter (NDIM=3)
      integer it1,iu1,iv1
      integer it2,iu2,iv2
      integer it3,iu3,iv3
      double precision p1
      double precision p2
      double precision p3
      double precision gam12
      double precision gam13
      double precision gam23
      double precision Pmat1(NDIM)
      double precision Pmat2(NDIM)
      double precision Pmat3(NDIM)
C Output varaibles
      double precision ans
C FunCall
      double precision binoc
C SubCall
C     gdel_x
C Local varaibles
      double precision PI
      parameter (PI=3.141592650d0)
      integer i
      integer i1,i2,i3
      integer j1,j2,j3
      integer k1,k2,k3
      double precision b1,b2,b3
      double precision covlp
      double precision c01D,c03D,c12,c13,c23
      double precision eta11,eta13,eta22,eta23

      double precision xsum,val
      double precision cov
      double precision d11,d12,d13
      double precision d21,d22,d23
      double precision d31,d32,d33
      double precision beta1,beta2,beta3

      double precision u,u1,u2,u3,v
      double precision xk1,xk2,xk3
      double precision t1,Tmat1(NDIM)
      double precision t2,Tmat2(NDIM)
      double precision s1,Smat1(NDIM)
      double precision s2,Smat2(NDIM)
      double precision gam12_new
      double precision h0

      double precision GX,GY,GZ,Rnlm
      double precision bt1,bt2,bt3
      double precision bu1,bu2,bu3
      double precision bv1,bv2,bv3

      double precision tmp1,tmp2
      double precision ZEROMAT(NDIM)
C-----------Initialize---------------
      do i=1,NDIM
         ZEROMAT(i) = 0.0d0
      end do
C-----------Integrate over x3---------------
      call G3remove(p3,gam13,gam23,c01D,c03D,c12,c13,c23)
C
C Integration over x3 gives one coefficient 
C c01D
C and three Gaussian functions:
C Exp[-c12(x1-x2)**2]
C Exp[-c13(x1-PX3)**2]
C Exp[-c23(x2-PX3)**2]
C
C update gam12 to:
C gam12 = gam12 + c12
C and carry out the following Gauss product
C
C Exp[-p1(x1-PX1)**2]*Exp[-c13(x1-PX3)**2] = K1 * Exp[-s1(x1-SX1)]
C Exp[-p2(x2-PX2)**2]*Exp[-c23(x2-PX3)**2] = K2 * Exp[-s2(x2-SX2)]
C
C Overlap integration of
C Integrate[dx1,dx2] Exp[-s1(x1-SX1)]*Exp[-s2(x2-SX2)]*Exp[-gam12*x12**2]  = cov * K3
C
C The final expression is: c01D * cov * K1 * K2 * K3
C
C cov = PI/(s1*s2 + s1*gam12 + s2*gam12)
C K1  = Exp[-t1*TX1**2]      = Exp[-t1(PX1-PX3)**2]
C K2  = Exp[-t2*TX2**2]      = Exp[-t2(PX2-PX3)**2]
C K3  = Exp[-u*(SX1-SX2)**2] = Exp[-u(e1*PX1 + e2*PX2 + e3*PX3)**2]
C u   = s1*s2*gam12/(s1*s2 + s1*gam12 + s2*gam12)
C
C SX1 = (p1/s1)*PX1 + (c13/s1)*PX3
C SX2 = (p2/s2)*PX2 + (c23/s2)*PX3
C SX1-SX2 = (p1/s1)*PX1 -(p2/s2)*PX2 + [ (c13/s1)-(c23/s2) ]*PX3
C
C TX1 = PX1-PX3
C TX2 = PX2-PX3
C
C
      gam12_new = gam12 + c12

      do i=1,NDIM
       call G3gprod_1D(p1,c13,Pmat3(i),Pmat3(i),t1,s1,Tmat1(i),Smat1(i))
       call G3gprod_1D(p2,c23,Pmat2(i),Pmat3(i),t2,s2,Tmat2(i),Smat2(i))
      end do

      tmp1= ((s1*s2)+(s1*gam12_new)+(s2*gam12_new))*sqrt(s1+s2)
      cov = 2.0d0*sqrt(PI**5)/tmp1

      u   = s1*s2*gam12_new/( (s1*s2)+(s1*gam12_new)+(s2*gam12_new) )

      tmp2= ((s1*s2)+(s1*gam12_new)+(s2*gam12_new))*(s1+s2)
      v   = s1*s1*s2*s2/tmp2  
C
C Using geneal expression:
C K1 = Exp[-beta1*(d11*PX1+d12*PX2+d13*PX3)**2]
C K2 = Exp[-beta2*(d21*PX1+d22*PX2+d23*PX3)**2]
C K3 = Exp[-beta3*(d31*PX1+d32*PX2+d33*PX3)**2]
C
C d11= 1      d12=  0       d13= -1
C d21= 0      d22=  1       d23= -1
C d31= e1     d22=  e2      d23= e3
C-----------prepare for deriv--------------------
      beta1=t1    ; beta2=t2     ; beta3=u
      d11  =1.0d0 ; d12  =0.0d0  ; d13  =-1.0d0
      d21  =0.0d0 ; d22  =1.0d0  ; d23  =-1.0d0
      d31  =p1/s1 ; d32  =-p2/s2 ; d33  =(c13/s1)-(c23/s2)
C-----------Derivative wrt PX1 PX2 and PX3--------
      xsum = 0.0d0
      do i1=0,it1
        bt1=binoc(it1,i1)
        do i2=0,it2
          bt2=binoc(it2,i2)
          do i3=0,it3
            bt3=binoc(it3,i3)
            call G3KKK_1D(it1-i1,it2-i2,it3-i3,
     1                    beta1,beta2,beta3,
     2                    d11,d12,d13,
     3                    d21,d22,d23,
     4                    d31,d32,d33,
     5                    Pmat1(1),Pmat2(1),Pmat3(1),
     6                    GX)
            do j1=0,iu1
              bu1=binoc(iu1,j1)
              do j2=0,iu2
                bu2=binoc(iu2,j2)
                do j3=0,iu3
                  bu3=binoc(iu3,j3)
                  call G3KKK_1D(iu1-j1,iu2-j2,iu3-j3,
     1                    beta1,beta2,beta3,
     2                    d11,d12,d13,
     3                    d21,d22,d23,
     4                    d31,d32,d33,
     5                    Pmat1(2),Pmat2(2),Pmat3(2),
     6                    GY)
                  do k1=0,iv1
                    bv1=binoc(iv1,k1)
                    do k2=0,iv2
                      bv2=binoc(iv2,k2)
                      do k3=0,iv3
                        bv3=binoc(iv3,k3)
                        call G3KKK_1D(iv1-k1,iv2-k2,iv3-k3,
     1                    beta1,beta2,beta3,
     2                    d11,d12,d13,
     3                    d21,d22,d23,
     4                    d31,d32,d33,
     5                    Pmat1(3),Pmat2(3),Pmat3(3),
     6                    GZ)

                        call G3Rnlm(NDIM,i1,j1,k1,i2,j2,k2,i3,j3,k3,
     1                    v,d31,d32,d33,
     2                    Pmat1,Pmat2,Pmat3,ZEROMAT,Rnlm)

                          val=bt1*bt2*bt3*
     1                        bu1*bu2*bu3*
     2                        bv1*bv2*bv3*
     3                        GX*GX*GZ*
     4                        Rnlm
                          xsum=xsum+val
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      ans = cov*c03D*val


      END

C************************************************************************
      SUBROUTINE G3hint0_1D(it1,it2,it3,p1,p2,p3,gam12,gam13,gam23,
     1                      PX1,PX2,PX3,ans)
C************************************************************************
C Three particle one dimensional Hermite Gaussian integral 
C
C Input varaibles
      implicit none
      integer it1
      integer it2
      integer it3
      double precision p1
      double precision p2
      double precision p3
      double precision gam12
      double precision gam13
      double precision gam23
      double precision PX1
      double precision PX2
      double precision PX3
C Output varaibles
      double precision ans
C FunCall
      double precision binoc
C SubCall
C     gdel_x
C Local varaibles
      double precision PI
      parameter (PI=3.141592650d0)
      integer i1,i2,i3
      double precision b1,b2,b3
      double precision covlp
      double precision c01D,c03D,c12,c13,c23
      double precision eta11,eta13,eta22,eta23

      double precision xsum,val
      double precision cov
      double precision d11,d12,d13
      double precision d21,d22,d23
      double precision d31,d32,d33
      double precision beta1,beta2,beta3

      double precision u,u1,u2,u3
      double precision xk1,xk2,xk3
      double precision t1,TX1
      double precision t2,TX2
      double precision s1,SX1
      double precision s2,SX2
      double precision gam12_new
      double precision h0
C-----------Integrate over x3---------------
      call G3remove(p3,gam13,gam23,c01D,c03D,c12,c13,c23)
C
C Integration over x3 gives one coefficient 
C c01D
C and three Gaussian functions:
C Exp[-c12(x1-x2)**2]
C Exp[-c13(x1-PX3)**2]
C Exp[-c23(x2-PX3)**2]
C
C update gam12 to:
C gam12 = gam12 + c12
C and carry out the following Gauss product
C
C Exp[-p1(x1-PX1)**2]*Exp[-c13(x1-PX3)**2] = K1 * Exp[-s1(x1-SX1)]
C Exp[-p2(x2-PX2)**2]*Exp[-c23(x2-PX3)**2] = K2 * Exp[-s2(x2-SX2)]
C
C Overlap integration of
C Integrate[dx1,dx2] Exp[-s1(x1-SX1)]*Exp[-s2(x2-SX2)]*Exp[-gam12*x12**2]  = cov * K3
C
C The final expression is: c01D * cov * K1 * K2 * K3
C
C cov = PI/(s1*s2 + s1*gam12 + s2*gam12)
C K1  = Exp[-t1*TX1**2]      = Exp[-t1(PX1-PX3)**2]
C K2  = Exp[-t2*TX2**2]      = Exp[-t2(PX2-PX3)**2]
C K3  = Exp[-u*(SX1-SX2)**2] = Exp[-u(e1*PX1 + e2*PX2 + e3*PX3)**2]
C u   = s1*s2*gam12/(s1*s2 + s1*gam12 + s2*gam12)
C
C SX1 = (p1/s1)*PX1 + (c13/s1)*PX3
C SX2 = (p2/s2)*PX2 + (c23/s2)*PX3
C SX1-SX2 = (p1/s1)*PX1 -(p2/s2)*PX2 + [ (c13/s1)-(c23/s2) ]*PX3
C
C TX1 = PX1-PX3
C TX2 = PX2-PX3
C
C
      gam12_new = gam12 + c12

      call G3gprod_1D(p1,c13,PX1,PX3,t1,s1,TX1,SX1)
      call G3gprod_1D(p2,c23,PX2,PX3,t2,s2,TX2,SX2)

      cov = PI/sqrt( (s1*s2)+(s1*gam12_new)+(s2*gam12_new) )
      u   = s1*s2*gam12_new/( (s1*s2)+(s1*gam12_new)+(s2*gam12_new) )

      xk1 = exp(-t1*TX1*TX1)
      xk2 = exp(-t2*TX2*TX2)
      xk3 = exp(-u*(SX1-SX2)*(SX1-SX2))

      h0 = c01D*cov*xk1*xk2*xk3
C
C Using geneal expression:
C K1 = Exp[-beta1*(d11*PX1+d12*PX2+d13*PX3)**2]
C K2 = Exp[-beta2*(d21*PX1+d22*PX2+d23*PX3)**2]
C K3 = Exp[-beta3*(d31*PX1+d32*PX2+d33*PX3)**2]
C
C d11= 1      d12=  0       d13= -1
C d21= 0      d22=  1       d23= -1
C d31= e1     d22=  e2      d23= e3
C-----------prepare for deriv--------------------
      beta1=t1    ; beta2=t2     ; beta3=u
      d11  =1.0d0 ; d12  =0.0d0  ; d13  =-1.0d0
      d21  =0.0d0 ; d22  =1.0d0  ; d23  =-1.0d0
      d31  =p1/s1 ; d32  =-p2/s2 ; d33  =(c13/s1)-(c23/s2)

      u1=(d11*PX1)+(d12*PX2)+(d13*PX3)
      u2=(d21*PX1)+(d22*PX2)+(d23*PX3)
      u3=(d31*PX1)+(d32*PX2)+(d33*PX3)

C-----------Derivative wrt PX1 PX2 and PX3--------
      call G3KKK_1D(it1,it2,it3,
     1              beta1,beta2,beta3,
     2              d11,d12,d13,
     3              d21,d22,d23,
     4              d31,d32,d33,
     5              PX1,PX2,PX3,
     6              val)

      ans = cov*c01D*val


      END



C**************************************************************************
      SUBROUTINE G3remove(p3,gam13,gam23,c01D,c03D,c12,c13,c23)
C**************************************************************************
C Integrates over particle 3 using Boys 2nd Theorem:
C
C FX = Integrate[-Inf,+Inf] dr3
C      exp[-p3*(r3-P3)**2] exp[-gam13*r13**2]  exp[-gam23*r23**2]
C
C Using Boys thm. the result is:
C
C FX = c03D*exp[-c12*r12**2]*exp[-c13*(r1-P3)**2]*exp[-c23*(r2-P3)**2]
C
C NOTE: If the integration is carried in 1D (i.e. over dx3 instead of dr3)
C       then the use c01D instead of c03D
C
C
      implicit none
C Input
      double precision p3,gam12,gam13,gam23
C Output
      double precision  c01D
      double precision  c03D
      double precision  c12
      double precision  c13
      double precision  c23
C Local
      double precision PI
      parameter (PI=3.141592650d0)
      double precision xsum 

      xsum = p3+gam13+gam23
      c01D = sqrt(PI/xsum)
      c03D = c01D**3

      c12=gam13*gam23/xsum
      c13=p3*gam13/xsum
      c23=p3*gam23/xsum

      END
C end of removeR3

C**************************************************************************
      SUBROUTINE G3KKK_1D(N1,N2,N3,
     1                    alp1,alp2,alp3,
     2                    a1,a2,a3,
     3                    b1,b2,b3,
     4                    c1,c2,c3,
     5                    PX1,PX2,PX3,
     6                    ans)
C**************************************************************************
C Calculate the derivative of the 3 Gaussian
C K1 = Exp[-alp1*(a1*PX1 + a2*PX2 + a3*PX3)**2]
C K2 = Exp[-alp2*(b1*PX1 + b2*PX2 + b3*PX3)**2]
C K3 = Exp[-alp3*(c1*PX1 + c2*PX2 + c3*PX3)**2]
C
C ans = (d/dPX1)^N1 (d/dPX2)^N2 (d/dPX3)^N3 (K1*K2*K3)
C
      implicit none
C Input
      integer N1,N2,N3
      double precision alp1,alp2,alp3
      double precision a1,a2,a3
      double precision b1,b2,b3
      double precision c1,c2,c3
      double precision PX1,PX2,PX3
C Output
      double precision ans
C FunCall
      double precision binoc
C SubCall
C     G3K_1D
C Local
      integer i1,i2,i3
      integer j1,j2,j3
      double precision xsum
      double precision d1,d2,d3
      double precision e1,e2,e3
      double precision v1,v2,v3

      xsum = 0.0d0
      do i1=0,N1
        d1=binoc(N1,i1)
        do i2=0,N2
          d2=binoc(N2,i2)
          do i3=0,N3
           d3=binoc(N3,i3)
           call G3K_1D(N1-i1,N2-i2,N3-i3,alp3,c1,c2,c3,PX1,PX2,PX3,v3)
           do j1=0,i1
             e1=binoc(j1,i1)
             do j2=0,i2
               e2=binoc(j2,i2)
               do j3=0,i3
                 e3=binoc(j3,i3)
                 call G3K_1D(j1,j2,j3,alp1,a1,a2,a3,PX1,PX2,PX3,v1)
                 call G3K_1D(i1-j1,i2-j2,i3-j3,alp2,b1,b2,b3,
     1                       PX1,PX2,PX3,v2)
                 xsum = xsum +(d1*d2*d3*e1*e2*e3*v1*v2*v3)
               end do
             end do
           end do
          end do
        end do
      end do
      ans = xsum

      END

C**************************************************************************
      SUBROUTINE G3K_1D(i1,i2,i3,alp,a,b,c,PX1,PX2,PX3,ans)
C**************************************************************************
C Calculates the following derivative
C
C ans = (d/dPX1)^i1 (d/dPX2)^i2 (d/dPX3)^i3 K
C
C K = Exp[ -alp(a*PX1 + b*PX2 + c*PX3)^2 ]
C
      implicit none
C Input
      integer i1
      integer i2
      integer i3
      double precision alp
      double precision a
      double precision b
      double precision c
      double precision PX1
      double precision PX2
      double precision PX3
C Output
      double precision ans
C FunCall
C     --no calls--     
C SubCall
C     gdel_x
C Local
      double precision u
      double precision gd
      double precision coeff
C
C u = a*PX1 + b*PX2 + c*PX3
C
C d/dPX1 = a d/du
C d/dPX2 = b d/du
C d/dPX2 = c d/du
C
      u =  (a*PX1) + (b*PX2) + (c*PX3)
      coeff = (a**i1)*(b**i2)*(c**i3)
      call gdel_x(i1+i2+i3,u,alp,0.0d0,gd)
      ans = coeff * gd

      END

C**************************************************************************
      SUBROUTINE G3gprod_1D(a,b,AX,BX,q,p,QX,PX)
C**************************************************************************
C Does Gaussian product rule
C Exp[-a*(x-AX)**2]*Exp[-b(x-BX)**2] = Exp[-q*QX**2]*Exp[-p(x-PX)**2]
C
      implicit none
C Input
      double precision a,b
      double precision AX,BX
C Output
      double precision q
      double precision QX
      double precision p
      double precision PX

      p = a+b
      q = a*b/(a+b)
      QX=AX-BX
      PX=( (a*AX)+(b*BX) )/(a+b)

      END 



C**************************************************************************
      SUBROUTINE XG3KKK_1D(N1,N2,N3,
     1                    alp1,alp2,alp3,
     2                    a1,a2,a3,
     3                    b1,b2,b3,
     4                    c1,c2,c3,
     5                    PX1,PX2,PX3,
     6                    ans)
C**************************************************************************
C Calculate the derivative of the 3 Gaussian
C K1 = Exp[-alp1*(a1*PX1 + a2*PX2 + a3*PX3)**2]
C K2 = Exp[-alp2*(b1*PX1 + b2*PX2 + b3*PX3)**2]
C K3 = Exp[-alp3*(c1*PX1 + c2*PX2 + c3*PX3)**2]
C
C ans = (d/dPX1)^N1 (d/dPX2)^N2 (d/dPX3)^N3 (K1*K2*K3)
C
      implicit none
C Input
      integer N1,N2,N3
      double precision alp1,alp2,alp3
      double precision a1,a2,a3
      double precision b1,b2,b3
      double precision c1,c2,c3
      double precision PX1,PX2,PX3
C Output
      double precision ans
C FunCall
      double precision binoc
C SubCall
C     G3K_1D
C Local
      integer i1,i2,i3
      integer j1,j2,j3
      double precision xsum
      double precision d1,d2,d3
      double precision e1,e2,e3
      double precision va,vb,vc
      double precision xk1,xk2,xk3
      double precision gk1,gk2,gk3

      va=(a1*PX1) + (a2*PX2) + (a3*PX3)
      vb=(b1*PX1) + (b2*PX2) + (b3*PX3)
      vc=(c1*PX1) + (c2*PX2) + (c3*PX3)

      xk1=dexp(-alp1*va*va)
      xk2=dexp(-alp2*vb*vb)
      xk3=dexp(-alp3*vc*vc)

      if(N1 .eq. 1) then
         gk1=xk1*(-2.0d0*alp1*va)*a1
         gk2=xk2*(-2.0d0*alp2*vb)*b1
         gk3=xk3*(-2.0d0*alp3*vc)*c1
         ans=(gk1*xk2*xk3)+(xk1*gk2*xk3)+(xk1*xk2*gk3)
      else
         gk1=xk1
         gk2=xk2
         gk3=xk3
         ans=xk1*xk2*xk3
      end if

      call G3K_1D(N1,N2,N3,alp1,a1,a2,a3,PX1,PX2,PX3,gk1)
      call G3K_1D(N1,N2,N3,alp2,b1,b2,b3,PX1,PX2,PX3,gk2)
      call G3K_1D(N1,N2,N3,alp3,c1,c2,c3,PX1,PX2,PX3,gk3)
      ans=(gk1*xk2*xk3)+(xk1*gk2*xk3)+(xk1*xk2*gk3)

      


      END

C**************************************************************************
      SUBROUTINE G3Rnlm(NDIM,
     1                  i1,j1,k1,
     2                  i2,j2,k2,
     3                  i3,j3,k3,
     4                  alp,d1,d2,d3,
     5                  Pmat1,Pmat2,Pmat3,Cmat,ans)
C**************************************************************************
C Derivative of Boys functions for 3 particle integrals
C (d/dPX1)^i1 (d/dPY1)^j1 (d/dPZ1)^k1
C (d/dPX2)^i2 (d/dPY2)^j2 (d/dPZ2)^k2
C (d/dPX3)^i3 (d/dPY3)^j1 (d/dPZ3)^k3   F(alp*S^2)
C
C S = d1*Pmat1 + d2*Pmat2 + d3*Pmat3 - Cmat
C
C
C
      implicit none
C Input
      integer NDIM
      integer i1,j1,k1
      integer i2,j2,k2
      integer i3,j3,k3
      double precision alp
      double precision d1,d2,d3
      double precision Pmat1(NDIM),Pmat2(NDIM),Pmat3(NDIM)
      double precision Cmat(NDIM)
C Output
      double precision ans
      double precision SX,SY,SZ
      double precision N1,N2,N3
      double precision Rnlm

C--------components of S-------------
      SX=(d1*Pmat1(1))+(d2*Pmat2(1))+(d3*Pmat3(1))-Cmat(1)
      SY=(d1*Pmat1(2))+(d2*Pmat2(2))+(d3*Pmat3(2))-Cmat(2)
      SZ=(d1*Pmat1(3))+(d2*Pmat2(3))+(d3*Pmat3(3))-Cmat(3)
C
C (d/dPX1) = (dSX/dPX1) * (d/dSX)
C          = d1 * (d/dSX)
C 
C (d/dPX1)^i1 (d/dPY1)^j1 (d/dPZ1)^k1 = d1**(i1+j1+k1) (d/dSX)**(i1+j1+k1)
C
C
      N1=i1+j1+k1
      N2=i2+j2+k2
      N3=i3+j3+k3

      call auxR2(N1,N2,N3,alp,SX,SY,SZ,Rnlm)
      ans = Rnlm*(d1**N1)*(d2**N2)*(d3**N3)

      END

C*****************************************************************
      SUBROUTINE G3Lap_1D(I1,alp1,Amat1,
     1                    I2,alp2,Amat2,
     1                    I3,alp3,Amat3,
     2                    L1,beta1,Bmat1,
     3                    L2,beta2,Bmat2,
     4                    L3,beta3,Bmat3,
     5                    gamA12,gamA13,gamA23,
     6                    gamB12,gamB13,gamB23,ans)
C*****************************************************************
C Primitive Gaussian Integral for overlap in 1D (PGIovlap)
C
C Calculates the e-e integral using
C Gaussian type geminal (GTG)
C and Cartesian Gaussian functions
C
C  ans = <GA(r1)GA(r2)|Exp[-gamma*r12^2]|GB(r1)GB(r2)>
C
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
      integer i,j,k
      integer n1,n2,n3
      double precision W1,W2,W3,W4,W5,W6
      double precision S2A,S2B,S3A,S3B,S4,S5,S6
      double precision t1,t2,t3
      double precision d1,d2,d3
      double precision c1,c2
      double precision gam12,gam13
      double precision XAB,XAC
      double precision S1(0:L1+2,0:L2+2,0:L3+2)
C-------------Initialize-------------------
      i = L1
      j = L2
      k = L3
      gam12=gamB12
      gam13=gamB13
      XAB = Bmat1-Bmat2
      XAC = Bmat1-Bmat3
C--------------Calculate ceoff-----------------------
C get the expansion coeff
      d1=0.0d0
      if(L1-1 .ge. 0) d1 = dble(L1*(L1-1))
      d2 = -2.0d0*beta1*dble((2*L1)+1)
      d3 = 4.0d0*beta1*beta1
      c1 = dble(L1)
      c2 = -2.0d0*beta1
C-------------Get overlap S1-------------------
      do n1=i-2,i+2
         if(n1 .lt. 0) cycle
         do n2=j-2,j+2
            if(n2 .lt. 0) cycle
            do n3=k-2,k+2
               if(n3 .lt. 0) cycle
               S1(n1,n2,n3) = 0.0d0
               call G3S1(I1,alp1,Amat1,
     1                I2,alp2,Amat2,
     1                I3,alp3,Amat3,
     2                n1,beta1,Bmat1,
     3                n2,beta2,Bmat2,
     4                n3,beta3,Bmat3,
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,S1(n1,n2,n3))
            end do
         end do
      end do
C-------------Get overlap S2-------------------
      t1 = 0.0d0
      t2 = 0.0d0
      if(i-1 .ge. 0) then
          t1 = S1(i-1,j+1,k)
          t2 = S1(i-1,j,k)
      end if
      
      S2A=S1(i,j,k)  -t1+(XAB*t2)
      S2B=S1(i+2,j,k)-S1(i+1,j+1,k)+(XAB*S1(i+1,j,k))

C-------------Get overlap S3-------------------
      t1 = 0.0d0
      t2 = 0.0d0
      if(i-1 .ge. 0) then
          t1 = S1(i-1,j,k+1)
          t2 = S1(i-1,j,k)
      end if

      S3A=S1(i,j,k)  -t1+(XAC*t2)
      S3B=S1(i+2,j,k)-S1(i+1,j,k+1)+(XAC*S1(i+1,j,k))

C-------------Get overlap S4-------------------
      S4=0.0d0
      S4=          S1(i+2,j,k) 
     1  -          S1(i+1,j,k+1) 
     2  -( XAC*    S1(i+1,j,k)   )
     3  -(         S1(i+1,j+1,k) )
     4  +(         S1(i,j+1,k+1) )
     5  -( XAC*    S1(i,j+1,k)   )
     6  +( XAB*    S1(i+1,j,k)   )
     7  -( XAB*    S1(i,j+1,k)   )
     8  +( XAB*XAC*S1(i,j,k)     )
C-------------Get overlap S5-------------------

      S5=0.0d0
      S5=            S1(i+2,j,k)+
     1               S1(i,j+2,k)+
     2   (XAB*XAB*   S1(i,j,k)) + 
     3   (-2.0d0*    S1(i+1,j+1,k)) +
     4   (-2.0d0*XAB*S1(i,j+1,k)) +
     4   (2.0d0*XAB* S1(i+1,j,k)) 
      
C-------------Get overlap S6-------------------
      S6=0.0d0
      S6=            S1(i+2,j,k)+
     1               S1(i,j,k+2)+
     2   (XAC*XAC*   S1(i,j,k)) + 
     3   (-2.0d0*    S1(i+1,j,k+1)) +
     4   (-2.0d0*XAC*S1(i,j,k+1)) +
     4   (2.0d0*XAC* S1(i+1,j,k)) 
      
C-------------Calculate W1-------------------
      W1 = 0.0d0
      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      if(i-2 .ge. 0) t1=d1*S1(i-2,j,k)
      t2 = d2*S1(i,j,k)
      t3 = d3*S1(i+2,j,k)
      W1 = t1 + t2 + t3
C-------------Calculate W2-------------------
      W2 = 0.0d0
      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      t1 = -2.0d0*gam12*S1(i,j,k)
      t2 = 4.0d0*gam12*gam12*S5
      W2 = t1 + t2
C-------------Calculate W3-------------------
      W3 = 0.0d0
      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      t1 = -2.0d0*gam13*S1(i,j,k)
      t2 = 4.0d0*gam13*gam13*S6
      W3 = t1 + t2
C-------------Calculate W4-------------------
      W4 = 0.0d0
      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      if(i-1 .ge. 0) t1 = -2.0d0*gam12*c1*S2A
      t2 = -2.0d0*gam12*c2*S2B
      W4 = t1 + t2
C-------------Calculate W5-------------------
      W5 = 0.0d0
      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      if(i-1 .ge. 0) t1 = -2.0d0*gam13*c1*S3A
      t2 = -2.0d0*gam13*c2*S3B
      W5 = t1 + t2
C-------------Calculate W6-------------------
      W6 = 0.0d0
      t1 = 0.0d0
      t2 = 0.0d0
      t3 = 0.0d0
      W6=4.0d0*gam12*gam13*S4
C-------------final ans-------------------
      ans = W1+W2+W3+(2.0d0*(W4+W5+W6))
      
      
      END
      
C*****************************************************************
      SUBROUTINE G3S1(I1,alp1,Amat1,
     1                I2,alp2,Amat2,
     1                I3,alp3,Amat3,
     2                L1,beta1,Bmat1,
     3                L2,beta2,Bmat2,
     4                L3,beta3,Bmat3,
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,ans)
C*****************************************************************
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
      ans = 0.0d0
C     if( (I1.ge.0).and.(I2.ge.0).and.(I3.ge.0).and.
C    1 (L1.ge.0).and.(L2.ge.0).and.(L3.ge.0) ) then

            call G3ovlap_1D(I1,alp1,Amat1,
     1                I2,alp2,Amat2,
     1                I3,alp3,Amat3,
     2                L1,beta1,Bmat1,
     3                L2,beta2,Bmat2,
     4                L3,beta3,Bmat3,
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,ans)

C     end if

      END

C*****************************************************************
      SUBROUTINE G3ke(I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     4                gamA12,gamA13,gamA23,
     4                gamB12,gamB13,gamB23,xmass,ansKE)
C*****************************************************************
C Geminal 3-particle KE  matrix element
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
      call G3Lap_1D(I1,alp1,Amat1(1),
     1                I2,alp2,Amat2(1),
     1                I3,alp3,Amat3(1),
     2                L1,beta1,Bmat1(1),
     3                L2,beta2,Bmat2(1),
     4                L3,beta3,Bmat3(1),
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,ddx)

      call G3Lap_1D(J1,alp1,Amat1(2),
     1                J2,alp2,Amat2(2),
     1                J3,alp3,Amat3(2),
     2                M1,beta1,Bmat1(2),
     3                M2,beta2,Bmat2(2),
     4                M3,beta3,Bmat3(2),
     5                gamA12,gamA13,gamA23,
     6                gamB12,gamB13,gamB23,ddy)

      call G3Lap_1D(K1,alp1,Amat1(3),
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
CCWS_DEBUG
c     write(*,*) 'From G3ke Sx    = ',Sx
c     write(*,*) 'From G3ke Sy    = ',Sy
c     write(*,*) 'From G3ke Sz    = ',Sz
c     write(*,*) 'From G3ke ddx   = ',ddx
c     write(*,*) 'From G3ke ddy   = ',ddy
c     write(*,*) 'From G3ke ddz   = ',ddz
c     write(*,*) 'From G3ke Tx    = ',Tx
c     write(*,*) 'From G3ke Ty    = ',Ty
c     write(*,*) 'From G3ke Tz    = ',Tz
c     write(*,*) 'From G3ke ansKE = ',ansKE
CCWS_DEBUG
      END

C*****************************************************************
      SUBROUTINE G3ovlap(I1,J1,K1,alp1,Amat1,
     1                I2,J2,K2,alp2,Amat2,
     1                I3,J3,K3,alp3,Amat3,
     2                L1,M1,N1,beta1,Bmat1,
     3                L2,M2,N2,beta2,Bmat2,
     3                L3,M3,N3,beta3,Bmat3,
     4                gamA12,gamA13,gamA23,
     4                gamB12,gamB13,gamB23,sval)
C*****************************************************************
C Geminal 3-particle KE  matrix element
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
C Output variables
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

CCWS_DEBUG
c     write(*,*) 'From G3ovlap_1D Sx = ',Sx
c     write(*,*) 'From G3ovlap_1D Sy = ',Sy
c     write(*,*) 'From G3ovlap_1D Sz = ',Sz
CCWS_DEBUG
      sval = Sx * Sy * Sz
      END
