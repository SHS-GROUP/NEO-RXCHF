C*************************************************************************
      SUBROUTINE G4Vee_r34(NQUAD_coul,NQUAD_ovlap,
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
C works only for S functions
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
      integer i
      double precision g14
      double precision g24
      double precision vee
      double precision c5,p5,Pmat5(NDIM)
      double precision c6,p6,Pmat6(NDIM)
      double precision c7,p7,Pmat7(NDIM)
      double precision c8,p8,Pmat8(NDIM)

      double precision xtmp
      double precision xjunk1
      double precision xjunk2
      double precision xjunk3
C------------------------------------------------------------------------------
C ans = <GA(1)GA(2)GA(3)GB(4)|g(1,4)g(2,4)g(3,4)/(r3-r4)|GB(1)GB(2)GB(3)GB(4)>
C f1(r_4)  = <GA(1)|g(1,4)|GB(1)> = c5 * EXP[-p5(r_4-Pmat5)^2]
C f2(r_4)  = <GA(2)|g(2,4)|GB(2)> = c6 * EXP[-p6(r_4-Pmat6)^2]
C f3(r_4)  = f1 * f2              = c7 * EXP[-p7(r_4-Pmat7)^2]
C
C------------------------------------------------------------------------------
C integrate over r_1 = c5 * EXP[-p5(r_4-PX6)^2]
C
      g14 = gamA14 + gamB14
      c5  = 1.0d0
      do i=1,NDIM
         call three_term_gauss(alp1,beta1,g14,
     1                         Amat1(i),Bmat1(i),xtmp,p5,Pmat5(i))
         c5 = c5 * xtmp
      end do
C------------------------------------------------------------------------------
C integrate over r_2 = c6 * EXP[-p6(r_4-PX6)^2]
C
      g24 = gamA24 + gamB24
      c6  = 1.0d0
      do i=1,NDIM
         call three_term_gauss(alp2,beta2,g24,
     1                         Amat2(i),Bmat2(i),xtmp,p6,Pmat6(i))
         c6 = c6 * xtmp
      end do
C------------------------------------------------------------------------------
C product of  EXP[-p5(r_4-PX6)^2] * EXP[-p6(r_4-PX6)^2]
C  = c7 * EXP[-p7(r_4-PX7)^2]
C 
      c7  = 1.0d0
      do i=1,NDIM
         call gauss_g2_prod(p5,p6,Pmat5(i),Pmat6(i),
     1                      p7,xjunk1,Pmat7(i),xjunk2,xtmp,xjunk3)
         c7 = c7 * xtmp
      end do
C------------------------------------------------------------------------------
C product of  EXP[-p7(r_4-PX7)^2] with GB(4)
C
      c8 = 1.0d0
      do i=1,NDIM
         call gauss_g2_prod(beta4,p7,Bmat4(i),Pmat7(i),
     1                      p8,xjunk1,Pmat8(i),xjunk2,xtmp,xjunk3)
         c8 = c8 * xtmp
      end do
C------------------------------------------------------------------------------
      call pgivee(I3,J3,K3,alp3,Amat3,
     1            I4,J4,K4,alp4,Amat4,
     2            L3,M3,N3,beta3,Bmat3,
     3            L4,M4,N4,p8,Pmat8,
     4            gamA34,gamB34,vee)
C------------------------------------------------------------------------------
      ans = c5*c6*c7*c8*vee
      call underflow(ans)


      END

