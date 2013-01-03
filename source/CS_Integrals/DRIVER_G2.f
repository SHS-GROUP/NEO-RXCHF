C Drivers for calculating select 2-particle (G2)
C integrals used in NEO-XC calculations.
C The routines here are for unit testing the integral codes.
C
C Routines in this file:
C G2_MD_DRIVER
C G2_AC_DRIVER
C G2_RI_DRIVER
C
C======================================================================
      subroutine G2_MD_DRIVER(I1,J1,K1,A1,Amat1,
     x                        I2,J2,K2,A2,Amat2,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L2,M2,N2,B2,Bmat2,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        xggs,xggVEC,xggVEP)

C  This is a driver used for unit testing G2 integrals.
C  It calls G2 McMurchie-Davidson routines to evaluate
C     xggs, xggVEC, and XggVEP integrals (returned).
C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2

      double precision A1
      double precision A2
      double precision B1
      double precision B2

      double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2
      double precision Cmat(3)
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Bmat1(3)
      double precision Bmat2(3)

C Variables Returned
      double precision xggs   
      double precision xggVEC
      double precision xggVEP
      
C Local variables


c        call cws_gam1_xggs(I1,J1,K1,A1,Amat1,
         call G2_MD_xggs(I1,J1,K1,A1,Amat1,
     1                   I2,J2,K2,A2,Amat2,
     2                   L1,M1,N1,B1,Bmat1,
     3                   L2,M2,N2,B2,Bmat2,
     4                   gamma1,gamma2,xggs)


c        call cws_gam1_xggvee(I1,J1,K1,A1,Amat1,
         call G2_MD_xggvee(I1,J1,K1,A1,Amat1,
     1                     I2,J2,K2,A2,Amat2,
     2                     L1,M1,N1,B1,Bmat1,
     3                     L2,M2,N2,B2,Bmat2,
     4                     gamma1,gamma2,xggVEP)
         xggVEP=-1.0d+00*xggVEP


c        call cws_gam1_xggvec(I1,J1,K1,A1,Amat1,
         call G2_MD_xggvec(I1,J1,K1,A1,Amat1,
     1                     I2,J2,K2,A2,Amat2,
     2                     L1,M1,N1,B1,Bmat1,
     3                     L2,M2,N2,B2,Bmat2,
     4                     gamma1,gamma2,Cmat,ZNUC,
     5                     xggVEC)
         xggVEC=-1.0d+00*ZNUC*xggVEC


      return
      end

C======================================================================
      subroutine G2_AC_DRIVER(I1,J1,K1,A1,Amat1,
     x                        I2,J2,K2,A2,Amat2,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L2,M2,N2,B2,Bmat2,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        xggs,xggVEC,xggVEP)

C  This is a driver used for unit testing G2 integrals.
C  It calls G2 McMurchie-Davidson routines of Ari to evaluate
C     xggs, xggVEC, and XggVEP integrals (returned).
C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2

      double precision A1
      double precision A2
      double precision B1
      double precision B2

      double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2
      double precision Cmat(3)
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Bmat1(3)
      double precision Bmat2(3)

C Variables Returned
      double precision xggs   
      double precision xggVEC
      double precision xggVEP
      
C Local variables


C     xggs
C     --g^A(e1,p1) g^B(e1,p1)--
      call pgiovlap(I1,J1,K1,A1,Amat1,
     1              I2,J2,K2,A2,Amat2,
     2              L1,M1,N1,B1,Bmat1,
     3              L2,M2,N2,B2,Bmat2,
     4              gamma1,gamma2,xggs)


C     xgVEPg
C     --g(e1,p1) V^ep(e1,p1)--
      CALL PGIVEE(I1,J1,K1,A1,Amat1,
     1            I2,J2,K2,A2,Amat2,
     2            L1,M1,N1,B1,Bmat1,
     3            L2,M2,N2,B2,Bmat2,
     4            gamma1,gamma2,xggVEP)
      xggVEP=-1.0d+00*xggVEP


C     xgVECg 
C     -- g^A(e1,p1) V^eC(e1) g^B(e1,p1) --
      CALL PGIVEC(I1,J1,K1,A1,Amat1,
     1            I2,J2,K2,A2,Amat2,
     2            L1,M1,N1,B1,Bmat1,
     3            L2,M2,N2,B2,Bmat2,
     4            gamma1,gamma2,CMAT,xggVEC)
      xggVEC=-1.0d+00*ZNUC*xggVEC



      return
      end

C======================================================================
      subroutine G2_RI_DRIVER(I1,J1,K1,A1,Amat1,
     x                        I2,J2,K2,A2,Amat2,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L2,M2,N2,B2,Bmat2,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        xggVEC,xggVEP)

C  This is a driver used for unit testing G2 integrals.
C  It calls G2 RI routines to evaluate
C     xggVEC, and XggVEP integrals (returned).
C======================================================================
      implicit none

C Input variables
c     logical orth_abs

      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2

      double precision A1
      double precision A2
      double precision B1
      double precision B2

      double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2
      double precision Cmat(3)
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Bmat1(3)
      double precision Bmat2(3)

C Variables Returned
      double precision xggVEC
      double precision xggVEP

c     double precision xgVEC
c     double precision xgVPC
c     double precision xgVEP
c     double precision xgVEP_DBG
c     double precision xgVECg
c     double precision xgVPCg
c     double precision xgVEPg

C Local variables
      double precision zero
      parameter(zero=0.0d+00)


C Orthonormalize the electron and proton ABS
c     if(orth_abs) call OrthoABS


c     call RI_G2_xgVEC(I1,J1,K1,A1,Amat1,
c    *                 I2,J2,K2,A2,Amat2,
c    *                 L1,M1,N1,B1,Bmat1,
c    *                 L2,M2,N2,B2,Bmat2,
c    *                 gamma1,zero,
c    *                 Cmat,ZNUC,
c    *                 xgVEC)

      call RI_G2_xgVEC(I1,J1,K1,A1,Amat1,
     *                 I2,J2,K2,A2,Amat2,
     *                 L1,M1,N1,B1,Bmat1,
     *                 L2,M2,N2,B2,Bmat2,
     *                 gamma1,gamma2,
     *                 Cmat,ZNUC,
     *                 xggVEC)

c     call RI_G2_xgVPC(I1,J1,K1,A1,Amat1,
c    *                 I2,J2,K2,A2,Amat2,
c    *                 L1,M1,N1,B1,Bmat1,
c    *                 L2,M2,N2,B2,Bmat2,
c    *                 gamma1,zero,
c    *                 Cmat,ZNUC,
c    *                 xgVPC)

c     call RI_G2_xgVPC(I1,J1,K1,A1,Amat1,
c    *                 I2,J2,K2,A2,Amat2,
c    *                 L1,M1,N1,B1,Bmat1,
c    *                 L2,M2,N2,B2,Bmat2,
c    *                 gamma1,gamma2,
c    *                 Cmat,ZNUC,
c    *                 xgVPCg)

c     call RI_G2_xgVEP(I1,J1,K1,A1,Amat1,
c    *                 I2,J2,K2,A2,Amat2,
c    *                 L1,M1,N1,B1,Bmat1,
c    *                 L2,M2,N2,B2,Bmat2,
c    *                 gamma1,zero,
c    *                 xgVEP)

c     call RI_G2_xgVEP_DEBUG(I1,J1,K1,A1,Amat1,
c     call RI_G2_xgVEP_alt(I1,J1,K1,A1,Amat1,
c    *                     I2,J2,K2,A2,Amat2,
c    *                     L1,M1,N1,B1,Bmat1,
c    *                     L2,M2,N2,B2,Bmat2,
c    *                     gamma1,zero,
c    *                     xgVEP_DBG)
c     xgVEP_DBG=0.0d+00 

      call RI_G2_xgVEPg(I1,J1,K1,A1,Amat1,
     *                  I2,J2,K2,A2,Amat2,
     *                  L1,M1,N1,B1,Bmat1,
     *                  L2,M2,N2,B2,Bmat2,
     *                  gamma1,gamma2,
     *                  xggVEP)



      return
      end
