C*************************************************************************
      SUBROUTINE interface_G4Vee_r34(NQUAD_coul,NQUAD_ovlap,
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
      double precision xjunk

      integer JVAL
C------------------------------------------------------------------------------
C ans = <GA(1)GA(2)GA(3)GB(4)|g(1,4)g(2,4)g(3,4)/(r3-r4)|GB(1)GB(2)GB(3)GB(4)>
C------------------------------------------------------------------------------
      JVAL = I1+J1+K1
     1     + I2+J2+K2
     1     + I3+J3+K3
     1     + I4+J4+K4
     2     + L1+M1+N1
     3     + L2+M2+N2
     3     + L3+M3+N3
     3     + L4+M4+N4

      if(JVAL .eq. 0) then
         call G4Vee_r34(NQUAD_coul,NQUAD_ovlap,
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
      else
         call rys_G4Vee_r34(NQUAD_coul,NQUAD_ovlap,
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
      end if


      END

