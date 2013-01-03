C Drivers for calculating select 3-particle (G3)
C integrals used in NEO-XC calculations.
C The routines here are for unit testing the integral codes.
C
C Routines in this file:
C G3_MD_DRIVER
C G3_AC_DRIVER
C G3_RI_DRIVER
C G3_AE_DRIVER
C
C======================================================================
      subroutine G3_MD_DRIVER(I1,J1,K1,A1,Amat1,
     x                        I2,J2,K2,A2,Amat2,
     x                        I3,J3,K3,A3,Amat3,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L2,M2,N2,B2,Bmat2,
     x                        L3,M3,N3,B3,Bmat3,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        xgsg,xgTE,xgTEg1,xgTEg3,
     x                        xgTEg2,xgTPg,
     x                        xgVEC,xgVPCg,
     x                        xgVECg1,xgVECg2,xgVEP,
     x                        xgVEE,xgVEPg1,xgVEPg2,
     x                        xgVEEg1,xgVEEg2)

C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3

      double precision A1
      double precision A2
      double precision A3
      double precision B1
      double precision B2
      double precision B3

      double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2
      double precision Cmat(3)
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)

C Variables Returned
      double precision xgsg   
      double precision xgTE   
      double precision xgTEg1   
      double precision xgTEg2   
      double precision xgTEg3   
      double precision xgTPg   
      double precision xgVEC
      double precision xgVPCg   
      double precision xgVECg1
      double precision xgVECg2
      double precision xgVEP
      double precision xgVEE
      double precision xgVEPg1
      double precision xgVEPg2
      double precision xgVEEg1
      double precision xgVEEg2

C Local variables
      double precision xmass
      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23
      double precision coul_sign
      double precision zero
      parameter(zero=0.0d+00)


         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)

c        gamA12=gamma1
c        gamB13=gamma2

c        call cws_gam2_xggs(I3,J3,K3,A3,Amat3,
c        call G3_MD_xggs(I3,J3,K3,A3,Amat3,
c    *                   I2,J2,K2,A2,Amat2,
c    *                   I1,J1,K1,A1,Amat1,
c    *                   L3,M3,N3,B3,Bmat3,
c    *                   L2,M2,N2,B2,Bmat2,
c    *                   L1,M1,N1,B1,Bmat1,
c    *                   gamA12,gamA13,gamA23,
c    *                   gamB12,gamB13,gamB23,
c    *                   xgsg)

         gamA13=gamma1
         gamB23=gamma2 
         

         call xG3_MD_xggs(I1,J1,K1,A1,Amat1,
     *                    I2,J2,K2,A2,Amat2,
     *                    I3,J3,K3,A3,Amat3,
     *                    L1,M1,N1,B1,Bmat1,
     *                    L2,M2,N2,B2,Bmat2,
     *                     L3,M3,N3,B3,Bmat3,
     *                    gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23,
     *                    xgsg)
      write(*,*)'XGSG   =',xgsg

         call G3_MD_xggs(I1,J1,K1,A1,Amat1,
     *                   I2,J2,K2,A2,Amat2,
     *                   I3,J3,K3,A3,Amat3,
     *                   L1,M1,N1,B1,Bmat1,
     *                   L2,M2,N2,B2,Bmat2,
     *                   L3,M3,N3,B3,Bmat3,
     *                   gamA12,gamA13,gamA23,
     *                   gamB12,gamB13,gamB23,
     *                   xgsg)

C----- Kinetic Energy-----
         xgTE=zero
         xgTEg1=zero
         xgTEg2=zero
         xgTEg3=zero
         xgTPg=zero

C        xgTE  (2,3)
C        --g(e2,p1) T^e(e1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         xmass=1.0d+00
c     call G3_MD_TEg013(I1,J1,K1,A1,Amat1,
c    *                  I2,J2,K2,A2,Amat2,
c    *                  I3,J3,K3,A3,Amat3,
c    *                  L1,M1,N1,B1,Bmat1,
c    *                  L2,M2,N2,B2,Bmat2,
c    *                  L3,M3,N3,B3,Bmat3,
c    *                  gamA12,gamA13,gamA23,
c    *                  gamB12,gamB13,gamB23,
c    *                  xgTE)
      call G3_MD_KE(I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              I3,J3,K3,A3,Amat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              L3,M3,N3,B3,Bmat3,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              xmass,xgTE)

C        xgTEg1  
C        --gA(e1,p1) T^e(e1) gB(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA13=gamma1
         gamB23=gamma2
         xmass=1.0d+00
c     call G3_MD_TEg013(I1,J1,K1,A1,Amat1,
c    *                  I2,J2,K2,A2,Amat2,
c    *                  I3,J3,K3,A3,Amat3,
c    *                  L1,M1,N1,B1,Bmat1,
c    *                  L2,M2,N2,B2,Bmat2,
c    *                  L3,M3,N3,B3,Bmat3,
c    *                  gamA12,gamA13,gamA23,
c    *                  gamB12,gamB13,gamB23,
c    *                  xgTEg1)
      call G3_MD_KE(I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              I3,J3,K3,A3,Amat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              L3,M3,N3,B3,Bmat3,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              xmass,xgTEg1)

C        xgTEg3  
C        --gA(e2,p1) T^e(e1) gB(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         gamB23=gamma2
         xmass=1.0d+00
c     call G3_MD_TEg013(I1,J1,K1,A1,Amat1,
c    *                  I2,J2,K2,A2,Amat2,
c    *                  I3,J3,K3,A3,Amat3,
c    *                  L1,M1,N1,B1,Bmat1,
c    *                  L2,M2,N2,B2,Bmat2,
c    *                  L3,M3,N3,B3,Bmat3,
c    *                  gamA12,gamA13,gamA23,
c    *                  gamB12,gamB13,gamB23,
c    *                  xgTEg3)
      call G3_MD_KE(I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              I3,J3,K3,A3,Amat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              L3,M3,N3,B3,Bmat3,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              xmass,xgTEg3)

C        xgTEg2  
C        --gA(e2,p1) T^e(e1) gB(e1,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         gamB13=gamma2
         xmass=1.0d+00
      call G3_MD_KE(I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              I3,J3,K3,A3,Amat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              L3,M3,N3,B3,Bmat3,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              xmass,xgTEg2)

C        xgTPg  
C        --gA(e1,p1) T^p(p) gB(e2,p1)--
C        --gA(1,3) T^p(3) gB(2,3)--
C        --gA(3,1) T^p(1) gB(2,1)--
C        --gA(2,1) T^p(1) gB(3,1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA12=gamma1
         gamB13=gamma2
         xmass=1836.0d+00
      call G3_MD_KE(I3,J3,K3,A3,Amat3,
     x              I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              L3,M3,N3,B3,Bmat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              xmass,xgTPg)


C        xgVEC  (2,3)
C        --g(e2,p1) V^{eC}(e1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         coul_sign=-1.0d+00
c        call G3_CWS_xgVEC(I1,J1,K1,A1,Amat1,
         call G3_MD_xgVxCg(I1,J1,K1,A1,Amat1,
     *                     I2,J2,K2,A2,Amat2,
     *                     I3,J3,K3,A3,Amat3,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L2,M2,N2,B2,Bmat2,
     *                     L3,M3,N3,B3,Bmat3,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     Cmat,ZNUC,
     *                     xgVEC)

          xgVEC=xgVEC*coul_sign*ZNUC
C
C        xgVPCg (1,2)  (1,3)
C        --g(e1,p1) V^{PC}(p1) g(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA12=gamma1
         gamB13=gamma2
         coul_sign=1.0d+00
c        call cws_gam2_xgVxCg(I3,J3,K3,A3,Amat3,
         call G3_MD_xgVxCg(I3,J3,K3,A3,Amat3,
     *                     I1,J1,K1,A1,Amat1,
     *                     I2,J2,K2,A2,Amat2,
     *                     L3,M3,N3,B3,Bmat3,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L2,M2,N2,B2,Bmat2,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     Cmat,ZNUC,
     *                     xgVPCg)
          xgVPCg=xgVPCg*coul_sign*ZNUC

C
C        xgVECg1 (1,3)  (2,3)
C        --g(e1,p1) V^{eC}(e1) g(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA13=gamma1
         gamB23=gamma2
         coul_sign=-1.0d+00
c        call cws_gam2_xgVxCg(I1,J1,K1,A1,Amat1,
         call G3_MD_xgVxCg(I1,J1,K1,A1,Amat1,
     *                     I2,J2,K2,A2,Amat2,
     *                     I3,J3,K3,A3,Amat3,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L2,M2,N2,B2,Bmat2,
     *                     L3,M3,N3,B3,Bmat3,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     Cmat,ZNUC,
     *                     xgVECg1)

          xgVECg1=xgVECg1*coul_sign*ZNUC

C        xgVECg2 (2,3) (2,3)
C        --g(e2,p1) V^{eC}(e1) g(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         gamB23=gamma2
         coul_sign=-1.0d+00
c        call cws_G3_xgVEC(I1,J1,K1,A1,Amat1,
c        call G3_CWS_xgVEC(I1,J1,K1,A1,Amat1,
         call G3_MD_xgVxCg(I1,J1,K1,A1,Amat1,
     *                     I2,J2,K2,A2,Amat2,
     *                     I3,J3,K3,A3,Amat3,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L2,M2,N2,B2,Bmat2,
     *                     L3,M3,N3,B3,Bmat3,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     Cmat,ZNUC,
     *                     xgVECg2)

          xgVECg2=xgVECg2*coul_sign*ZNUC

C        xgVEP  (2,3)
C        --g(e2,p1) V^{ep}(e1,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         coul_sign=-1.0d+00
c        call cws_gam2_xgVeeg(I1,J1,K1,A1,Amat1,
         call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                     I3,J3,K3,A3,Amat3,
     *                     I2,J2,K2,A2,Amat2,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L3,M3,N3,B3,Bmat3,
     *                     L2,M2,N2,B2,Bmat2,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     xgVEP)
         xgVEP=xgVEP*coul_sign

C        xgVEE  (2,3)
C        --g(e2,p1) V^{ee}(e1,e2)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
c        call cws_gam2_xgVeeg(I1,J1,K1,A1,Amat1,
         call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                     I2,J2,K2,A2,Amat2,
     *                     I3,J3,K3,A3,Amat3,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L2,M2,N2,B2,Bmat2,
     *                     L3,M3,N3,B3,Bmat3,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     xgVEE)

C        xgVEEg2  (1,3) (2,3)
C        --g(e1,p1) V^{ee}(e1,e2) g(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA13=gamma1
         gamB23=gamma2
c        call cws_gam2_xgVeeg(I1,J1,K1,A1,Amat1,
         call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                     I2,J2,K2,A2,Amat2,
     *                     I3,J3,K3,A3,Amat3,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L2,M2,N2,B2,Bmat2,
     *                     L3,M3,N3,B3,Bmat3,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     xgVEEg2)

C        xgVEEg1  (1,3) (1,3)
C        --g(e1,p1) V^{ee}(e1,e2) g(e1,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA13=gamma1
         gamB13=gamma2
c        call cws_gam2_xgVeeg(I1,J1,K1,A1,Amat1,
         call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                     I2,J2,K2,A2,Amat2,
     *                     I3,J3,K3,A3,Amat3,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L2,M2,N2,B2,Bmat2,
     *                     L3,M3,N3,B3,Bmat3,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     xgVeeg1)

C        xgVEPg1  (1,2) (2,3)
C        --g(e1,p1) V^{ep}(e1,p1) g(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA12=gamma1
         gamB23=gamma2
         coul_sign=-1.0d+00
c        call cws_gam2_xgVeeg(I1,J1,K1,A1,Amat1,
         call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                     I3,J3,K3,A3,Amat3,
     *                     I2,J2,K2,A2,Amat2,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L3,M3,N3,B3,Bmat3,
     *                     L2,M2,N2,B2,Bmat2,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     xgVEPg1)
         xgVEPg1=xgVEPg1*coul_sign

C        xgVEPg2  (2,3) (2,3)
C        --g(e2,p1) V^{ep}(e1,p1) g(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         gamB23=gamma2
         coul_sign=-1.0d+00
         call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                     I3,J3,K3,A3,Amat3,
     *                     I2,J2,K2,A2,Amat2,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L3,M3,N3,B3,Bmat3,
     *                     L2,M2,N2,B2,Bmat2,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     xgVEPg2)
         xgVEPg2=xgVEPg2*coul_sign


      return
      end

C======================================================================
      subroutine G3_RI_DRIVER(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamma1,gamma2,
     *                        ZNUC,Cmat,
     x                        xgVEC,xgVPCg,
     x                        xgVECg1,xgVECg2,xgVEP,
     x                        xgVEE,xgVEPg1,xgVEPg2,
     x                        xgVEEg1,xgVEEg2)

C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3

      double precision A1
      double precision A2
      double precision A3
      double precision B1
      double precision B2
      double precision B3

      double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2
      double precision Cmat(3)
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)

C Variables Returned
      double precision xgVEC
      double precision xgVPCg   
      double precision xgVECg1
      double precision xgVECg2
      double precision xgVEP
      double precision xgVEE
      double precision xgVEPg1
      double precision xgVEPg2
      double precision xgVEEg1
      double precision xgVEEg2

C Local variables
      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23

      double precision coul_sign


C     *RI* xgVEC  (2,3)
C     --g(e2,p1) V^{eC}(e1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
      gamA23=gamma1
c     call cws_G3_xgVEC(I1,J1,K1,A1,Amat1,
      call G3_CWS_xgVEC(I1,J1,K1,A1,Amat1,
     *                  I2,J2,K2,A2,Amat2,
     *                  I3,J3,K3,A3,Amat3,
     *                  L1,M1,N1,B1,Bmat1,
     *                  L2,M2,N2,B2,Bmat2,
     *                  L3,M3,N3,B3,Bmat3,
     *                  gamA12,gamA13,gamA23,
     *                  gamB12,gamB13,gamB23,
     *                  Cmat,ZNUC,
     *                  xgVEC)

      coul_sign=-1.0d+00
      xgVEC=xgVEC*coul_sign*ZNUC

C
C     xgVPCg (1,3)  (2,3)
C     --g(e1,p1) V^{PC}(p1) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA13=gamma1
      gamB23=gamma2
      call RI_G3_xgVPCg(I1,J1,K1,A1,Amat1,
     *                  I2,J2,K2,A2,Amat2,
     *                  I3,J3,K3,A3,Amat3,
     *                  L1,M1,N1,B1,Bmat1,
     *                  L2,M2,N2,B2,Bmat2,
     *                  L3,M3,N3,B3,Bmat3,
     *                  gamA12,gamA13,gamA23,
     *                  gamB12,gamB13,gamB23,
     *                  Cmat,ZNUC,
     *                  xgVPCg)

C
C     *RI* xgVECg1 (1,3)  (2,3)
C     --g(e1,p1) V^{eC}(e1) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA13=gamma1
      gamB23=gamma2
      call RI_G3_xgVECg1(I1,J1,K1,A1,Amat1,
     *                   I2,J2,K2,A2,Amat2,
     *                   I3,J3,K3,A3,Amat3,
     *                   L1,M1,N1,B1,Bmat1,
     *                   L2,M2,N2,B2,Bmat2,
     *                   L3,M3,N3,B3,Bmat3,
     *                   gamA12,gamA13,gamA23,
     *                   gamB12,gamB13,gamB23,
     *                   Cmat,ZNUC,
     *                   xgVECg1)

c     coul_sign=-1.0d+00
c     xgVECg1=xgVECg1*coul_sign

C
C     *RI* xgVECg2 (2,3)  (2,3)
C     --g(e2,p1) V^{eC}(e1) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA23=gamma1
      gamB23=gamma2
c     call cws_G3_xgVEC(I1,J1,K1,A1,Amat1,
      call G3_CWS_xgVEC(I1,J1,K1,A1,Amat1,
     *                  I2,J2,K2,A2,Amat2,
     *                  I3,J3,K3,A3,Amat3,
     *                  L1,M1,N1,B1,Bmat1,
     *                  L2,M2,N2,B2,Bmat2,
     *                  L3,M3,N3,B3,Bmat3,
     *                  gamA12,gamA13,gamA23,
     *                  gamB12,gamB13,gamB23,
     *                  Cmat,ZNUC,
     *                  xgVECg2)

      coul_sign=-1.0d+00
      xgVECg2=xgVECg2*coul_sign*ZNUC


C>>>>>>>>>>>>>>>>>>>> VEP/VEE <<<<<<<<<<<<<<<<<<<

C    *RI* xgVEP  (2,3)
C     --g(e2,p1) V^{ep}(e1,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA23=gamma1
      call RI_G3_xgVEP(I1,J1,K1,A1,Amat1,
     *                 I2,J2,K2,A2,Amat2,
     *                 I3,J3,K3,A3,Amat3,
     *                 L1,M1,N1,B1,Bmat1,
     *                 L2,M2,N2,B2,Bmat2,
     *                 L3,M3,N3,B3,Bmat3,
     *                 gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23,
     *                 xgVEP)

c     coul_sign=-1.0d+00
c     xgVEP=xgVEP*coul_sign

C     *RI* xgVEE  (2,3)
C     --g(e2,p1) V^{ee}(e1,e2)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA23=gamma1
      call RI_G3_xgVEE(I1,J1,K1,A1,Amat1,
     *                 I2,J2,K2,A2,Amat2,
     *                 I3,J3,K3,A3,Amat3,
     *                 L1,M1,N1,B1,Bmat1,
     *                 L2,M2,N2,B2,Bmat2,
     *                 L3,M3,N3,B3,Bmat3,
     *                 gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23,
     *                 xgVEE)

C     *RI* xgVEPg1  (1,3) (2,3)
C     --g(e1,p1) V^{ep}(e1,p1) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA13=gamma1
      gamB23=gamma2
      call RI_G3_xgVEPg1(I1,J1,K1,A1,Amat1,
     *                   I2,J2,K2,A2,Amat2,
     *                   I3,J3,K3,A3,Amat3,
     *                   L1,M1,N1,B1,Bmat1,
     *                   L2,M2,N2,B2,Bmat2,
     *                   L3,M3,N3,B3,Bmat3,
     *                   gamA12,gamA13,gamA23,
     *                   gamB12,gamB13,gamB23,
     *                   xgVEPg1)

c     coul_sign=-1.0d+00
c     xgVEPg1=xgVEPg1*coul_sign

C     *RI* xgVEPg2  (2,3) (2,3)
C     --g(e2,p1) V^{ep}(e1,p1) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA23=gamma1
      gamB23=gamma2
      call RI_G3_xgVEPg2(I1,J1,K1,A1,Amat1,
     *                   I2,J2,K2,A2,Amat2,
     *                   I3,J3,K3,A3,Amat3,
     *                   L1,M1,N1,B1,Bmat1,
     *                   L2,M2,N2,B2,Bmat2,
     *                   L3,M3,N3,B3,Bmat3,
     *                   gamA12,gamA13,gamA23,
     *                   gamB12,gamB13,gamB23,
     *                   xgVEPg2)

c     coul_sign=-1.0d+00
c     xgVEPg2=xgVEPg2*coul_sign

C     *RI* xgVEEg1  (1,3) (1,3)
C     --g(e1,p1) V^{ee}(e1,e2) g(e1,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA13=gamma1
      gamB13=gamma2
      call RI_G3_xgVEEg1(I1,J1,K1,A1,Amat1,
     *                   I2,J2,K2,A2,Amat2,
     *                   I3,J3,K3,A3,Amat3,
     *                   L1,M1,N1,B1,Bmat1,
     *                   L2,M2,N2,B2,Bmat2,
     *                   L3,M3,N3,B3,Bmat3,
     *                   gamA12,gamA13,gamA23,
     *                   gamB12,gamB13,gamB23,
     *                   xgVEEg1)


C     *RI* xgVEEg2  (1,3) (2,3)
C     --g(e1,p1) V^{ee}(e1,e2) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA13=gamma1
      gamB23=gamma2
      call RI_G3_xgVEEg2(I1,J1,K1,A1,Amat1,
     *                   I2,J2,K2,A2,Amat2,
     *                   I3,J3,K3,A3,Amat3,
     *                   L1,M1,N1,B1,Bmat1,
     *                   L2,M2,N2,B2,Bmat2,
     *                   L3,M3,N3,B3,Bmat3,
     *                   gamA12,gamA13,gamA23,
     *                   gamB12,gamB13,gamB23,
     *                   xgVEEg2)


      return
      end


C======================================================================
      subroutine G3_AC_DRIVER(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamma1,gamma2,
     *                        ZNUC,Cmat,
     x                        xgsg,xgTE,xgTEg1,xgTEg3,
     x                        xgTEg2,xgTPg,
     x                        xgVEC,xgVPCg,
     x                        xgVECg1,xgVECg2,xgVEP,
     x                        xgVEE,xgVEPg1,xgVEPg2,
     x                        xgVEEg1,xgVEEg2)

C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3

      double precision A1
      double precision A2
      double precision A3
      double precision B1
      double precision B2
      double precision B3

      double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2
      double precision Cmat(3)
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)

C Variables Returned
      double precision xgsg   
      double precision xgTE   
      double precision xgTEg1   
      double precision xgTEg2   
      double precision xgTEg3   
      double precision xgTPg   
      double precision xgVEC
      double precision xgVEE
      double precision xgVEP
      double precision xgVPCg
      double precision xgVECg1
      double precision xgVECg2
      double precision xgVEPg1
      double precision xgVEPg2
      double precision xgVEEg1
      double precision xgVEEg2

C Local variables
      integer NQUAD_coul
      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23
      double precision xmass
      double precision vec_ans


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NQUAD_coul=5
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)

         gamA13=gamma1
         gamB23=gamma2 
         
         call G3ovlap(I1,J1,K1,A1,Amat1,
     *                I2,J2,K2,A2,Amat2,
     *                I3,J3,K3,A3,Amat3,
     *                L1,M1,N1,B1,Bmat1,
     *                L2,M2,N2,B2,Bmat2,
     *                L3,M3,N3,B3,Bmat3,
     *                gamA12,gamA13,gamA23,
     *                gamB12,gamB13,gamB23,
     *                xgsg)

C     xgTE  (2,3)
C     --g(e2,p1) T^e(e1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA23=GAMMA1
      xmass=1.0d+00
      call G3ke_gnT1g2(I1,J1,K1,A1,Amat1,
     *                 I2,J2,K2,A2,Amat2,
     *                 I3,J3,K3,A3,Amat3,
     *                 L1,M1,N1,B1,Bmat1,
     *                 L2,M2,N2,B2,Bmat2,
     *                 L3,M3,N3,B3,Bmat3,
     *                 gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23,
     *                 xmass,xgTE)

C     xgVEC  (2,3)
C     --g(e2,p1) V^{eC}(e1)--
      call interface_G3Vec(NQUAD_coul,
     *              I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              I3,J3,K3,A3,Amat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              L3,M3,N3,B3,Bmat3,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              Cmat,vec_ans)

      call underflow(vec_ans)
      xgVEC=-ZNUC*vec_ans
 
C        xgVEP  (2,3)
C        --g(e2,p1) V^{ep}(e1,p1)--
      call interface_G3Vee(NQUAD_coul,
     *              I1,J1,K1,A1,Amat1,
     *              I3,J3,K3,A3,Amat3,
     *              I2,J2,K2,A2,Amat2,
     *              L1,M1,N1,B1,Bmat1,
     *              L3,M3,N3,B3,Bmat3,
     *              L2,M2,N2,B2,Bmat2,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              xgVEP)
      call underflow(xgVEP)
      xgVEP = -xgVEP

C        xgVEE  (2,3)
C        --g(e2,p1) V^{ee}(e1,e2)--
      call interface_G3Vee(NQUAD_coul,
     *              I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              I3,J3,K3,A3,Amat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              L3,M3,N3,B3,Bmat3,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              xgVEE)
      call underflow(xgVEE)

C  Sum up 1 gamma loop terms
c        gTE= gTE +bcoef(ik)*xgTE
c        gVEC=gVEC+bcoef(ik)*xgVEC
c        gVEP=gVEP+bcoef(ik)*xgVEP
c        gVEE=gVEE+bcoef(ik)*xgVEE

C --- end of 1 gamma loop ---
C
C --- 2 gamma loop ---
C
C%%%%%%%%%%--2 Gamma Kinetic Energy Integrals--%%%%%%%%%%

C             xgTPg (1,3) (2,3)
C-------------g(e1,p1) T^p(p1) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA13=GAMMA1
      gamB23=GAMMA2
c     xmass=pmass
      xmass=1836.0d+00
      call G3ke_g1Tpg2(I1,J1,K1,A1,Amat1,
     *                 I2,J2,K2,A2,Amat2,
     *                 I3,J3,K3,A3,Amat3,
     *                 L1,M1,N1,B1,Bmat1,
     *                 L2,M2,N2,B2,Bmat2,
     *                 L3,M3,N3,B3,Bmat3,
     *                 gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23,
     *                 xmass,xgTPg)

C     xgTEg1 (1,3) (2,3)
C-----g(e1,p1) T^e(e1) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     1              gamB12,gamB13,gamB23)
      gamA13=GAMMA1
      gamB23=GAMMA2
      xmass=1.d0
      call G3ke_gnT1g2(I1,J1,K1,A1,Amat1,
     1                 I2,J2,K2,A2,Amat2,
     1                 I3,J3,K3,A3,Amat3,
     2                 L1,M1,N1,B1,Bmat1,
     3                 L2,M2,N2,B2,Bmat2,
     3                 L3,M3,N3,B3,Bmat3,
     4                 gamA12,gamA13,gamA23,
     4                 gamB12,gamB13,gamB23,
     *                 xmass,xgTEg1)

C      xgTEg2 (2,3) (1,3)
C------g(e2,p1) T^e(e1) g(e1,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA23=GAMMA1
      gamB13=GAMMA2
      xmass=1.d0
      call G3ke_g2T1g1(I1,J1,K1,A1,Amat1,
     *                 I2,J2,K2,A2,Amat2,
     *                 I3,J3,K3,A3,Amat3,
     *                 L1,M1,N1,B1,Bmat1,
     *                 L2,M2,N2,B2,Bmat2,
     *                 L3,M3,N3,B3,Bmat3,
     *                 gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23,
     *                 xmass,xgTEg2)

C       xgTEg3 (2,3) (2,3)
C-------g(e2,p1)T^e(e1)g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA23=GAMMA1
      gamB23=GAMMA2
      xmass=1.d0
      call G3ke_gnT1g2(I1,J1,K1,A1,Amat1,
     *                 I2,J2,K2,A2,Amat2,
     *                 I3,J3,K3,A3,Amat3,
     *                 L1,M1,N1,B1,Bmat1,
     *                 L2,M2,N2,B2,Bmat2,
     *                 L3,M3,N3,B3,Bmat3,
     *                 gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23,
     *                 xmass,xgTEg3)


C%%%%%%%%%%--End of 2 Gamma Kinetic Energy Integrals--%%%%%%%%%%
C
C%%%%%%%%%%--2 Gamma VEC Integrals--%%%%%%%%%%

C
C     xgVPCg (1,2)  (1,3)
C     --g(e1,p1) V^{PC}(p1) g(e2,p1)--
C
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA12=GAMMA1
      gamB13=GAMMA2
C
      call interface_G3Vec(NQUAD_coul,
     *              I3,J3,K3,A3,Amat3,
     *              I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              L3,M3,N3,B3,Bmat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              Cmat,vec_ans)

      call underflow(vec_ans)
      xgVPCg=ZNUC*vec_ans

C     xgVECg1 (1,3)  (2,3)
C     --g(e1,p1) V^{eC}(e1) g(e2,p1)--
C
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA13=GAMMA1
      gamB23=GAMMA2

      call interface_G3Vec(NQUAD_coul,
     *              I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              I3,J3,K3,A3,Amat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              L3,M3,N3,B3,Bmat3,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              Cmat,vec_ans)

      call underflow(vec_ans)
      xgVECg1=-ZNUC*vec_ans

C     xgVECg2 (2,3)  (2,3)
C     --g(e2,p1) V^{eC}(e1) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA23=GAMMA1
      gamB23=GAMMA2

      call interface_G3Vec(NQUAD_coul,
     *              I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              I3,J3,K3,A3,Amat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              L3,M3,N3,B3,Bmat3,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              Cmat,vec_ans)
      call underflow(vec_ans)

      xgVECg2=-ZNUC*vec_ans
C End of loop over NAT (classical nuclei)

C%%%%%%%%%%--End of 2 Gamma VEC Integrals--%%%%%%%%%%

C%%%%%%%%%%--2 Gamma VEP Integrals--%%%%%%%%%%

C           xgVEPg1  (1,2) (2,3)
C           --g(e1,p1) V^{ep}(e1,p1) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA12=GAMMA1
      gamB23=GAMMA2

      call interface_G3Vee(NQUAD_coul,
     *              I1,J1,K1,A1,Amat1,
     *              I3,J3,K3,A3,Amat3,
     *              I2,J2,K2,A2,Amat2,
     *              L1,M1,N1,B1,Bmat1,
     *              L3,M3,N3,B3,Bmat3,
     *              L2,M2,N2,B2,Bmat2,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              xgVEPg1)

      call underflow(xgVEPg1)
      xgVEPg1=-xgVEPg1

C     xgVEPg2  (2,3) (2,3)
C     --g(e2,p1) V^{ep}(e1,p1) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA23=GAMMA1
      gamB23=GAMMA2

      call interface_G3Vee(NQUAD_coul,
     *              I1,J1,K1,A1,Amat1,
     *              I3,J3,K3,A3,Amat3,
     *              I2,J2,K2,A2,Amat2,
     *              L1,M1,N1,B1,Bmat1,
     *              L3,M3,N3,B3,Bmat3,
     *              L2,M2,N2,B2,Bmat2,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              xgVEPg2)

      call underflow(xgVEPg2)
      xgVEPg2=-xgVEPg2

C%%%%%%%%%%--End of 2 Gamma VEP Integrals--%%%%%%%%%%

C%%%%%%%%%%--2 Gamma VEE Integrals--%%%%%%%%%%

C     xgVEEg1  (1,3) (1,3)
C     --g(e1,p1) V^{ee}(e1,e2) g(e1,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23)
      gamA13=GAMMA1
      gamB13=GAMMA2

      call interface_G3Vee(NQUAD_coul,
     *              I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              I3,J3,K3,A3,Amat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              L3,M3,N3,B3,Bmat3,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              xgVEEg1)
      call underflow(xgVEEg1)

C     xgVEEg2  (1,3) (2,3)
C     --g(e1,p1) V^{ee}(e1,e2) g(e2,p1)--
      call zero_gam(gamA12,gamA13,gamA23,
     1              gamB12,gamB13,gamB23)
      gamA13=GAMMA1
      gamB23=GAMMA2

      call interface_G3Vee(NQUAD_coul,
     *              I1,J1,K1,A1,Amat1,
     *              I2,J2,K2,A2,Amat2,
     *              I3,J3,K3,A3,Amat3,
     *              L1,M1,N1,B1,Bmat1,
     *              L2,M2,N2,B2,Bmat2,
     *              L3,M3,N3,B3,Bmat3,
     *              gamA12,gamA13,gamA23,
     *              gamB12,gamB13,gamB23,
     *              xgVEEg2)
      call underflow(xgVEEg2)


      return
      end

C======================================================================
      subroutine G3_AE_DRIVER(I1,J1,K1,A1,Amat1,
     x                        I2,J2,K2,A2,Amat2,
     x                        I3,J3,K3,A3,Amat3,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L2,M2,N2,B2,Bmat2,
     x                        L3,M3,N3,B3,Bmat3,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        xgVEC,xgVPCg,
     x                        xgVECg1,xgVECg2,xgVEP,
     x                        xgVEE,xgVEPg1,xgVEPg2,
     x                        xgVEEg1)

C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3

      double precision A1
      double precision A2
      double precision A3
      double precision B1
      double precision B2
      double precision B3

      double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2
      double precision Cmat(3)
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)

C Variables Returned
      double precision xgVEC
      double precision xgVPCg   
      double precision xgVECg1
      double precision xgVECg2
      double precision xgVEP
      double precision xgVEE
      double precision xgVEPg1
      double precision xgVEPg2
      double precision xgVEEg1

C Local variables
      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23
      double precision coul_sign
      double precision zero
      parameter(zero=0.0d+00) 

      xgVEC  =zero
      xgVPCg =zero
      xgVECg1=zero
      xgVECg2=zero
      xgVEP  =zero
      xgVEE  =zero
      xgVEPg1=zero
      xgVEPg2=zero
      xgVEEg1=zero



C        xgVEC  (2,3)
C        --g(e2,p1) V^{eC}(e1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         call G3VeC_AUX_g13g23V1(I1,J1,K1,A1,Amat1,
     *                           I2,J2,K2,A2,Amat2,
     *                           I3,J3,K3,A3,Amat3,
     *                           L1,M1,N1,B1,Bmat1,
     *                           L2,M2,N2,B2,Bmat2,
     *                           L3,M3,N3,B3,Bmat3,
     *                           gamA13,gamA23,
     *                           gamB13,gamB23,
     *                           Cmat,ZNUC,xgVEC)
         xgVEC=-1.0d+00*xgVEC*ZNUC

C        xgVECg1 (1,3)  (2,3)
C        --g(e1,p1) V^{eC}(e1) g(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA13=gamma1
         gamB23=gamma2
         call G3VeC_AUX_g13g23V1(I1,J1,K1,A1,Amat1,
     *                           I2,J2,K2,A2,Amat2,
     *                           I3,J3,K3,A3,Amat3,
     *                           L1,M1,N1,B1,Bmat1,
     *                           L2,M2,N2,B2,Bmat2,
     *                           L3,M3,N3,B3,Bmat3,
     *                           gamA13,gamA23,
     *                           gamB13,gamB23,
     *                           Cmat,ZNUC,xgVECg1)
         xgVECg1=-1.0d+00*xgVECg1*ZNUC

C        xgVECg2 (2,3) (2,3)
C        --g(e2,p1) V^{eC}(e1) g(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         gamB23=gamma2
         call G3VeC_AUX_g13g23V1(I1,J1,K1,A1,Amat1,
     *                           I2,J2,K2,A2,Amat2,
     *                           I3,J3,K3,A3,Amat3,
     *                           L1,M1,N1,B1,Bmat1,
     *                           L2,M2,N2,B2,Bmat2,
     *                           L3,M3,N3,B3,Bmat3,
     *                           gamA13,gamA23,
     *                           gamB13,gamB23,
     *                           Cmat,ZNUC,xgVECg2)
         xgVECg2=-1.0d+00*xgVECg2*ZNUC

C
C        xgVPCg (1,2)  (1,3)
C        --g(e1,p1) V^{PC}(p1) g(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA13=gamma1
         gamB23=gamma2
         call G3VpC_AUX_g13g23V3(I1,J1,K1,A1,Amat1,
     *                           I2,J2,K2,A2,Amat2,
     *                           I3,J3,K3,A3,Amat3,
     *                           L1,M1,N1,B1,Bmat1,
     *                           L2,M2,N2,B2,Bmat2,
     *                           L3,M3,N3,B3,Bmat3,
     *                           gamA13,gamA23,
     *                           gamB13,gamB23,
     *                           Cmat,ZNUC,xgVPCg)
         xgVPCg=xgVPCg*ZNUC


C        xgVEP  (2,3)
C        --g(e2,p1) V^{ep}(e1,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         call G3Vep_AUX_g13g23V13(I1,J1,K1,A1,Amat1,
     *                            I2,J2,K2,A2,Amat2,
     *                            I3,J3,K3,A3,Amat3,
     *                            L1,M1,N1,B1,Bmat1,
     *                            L2,M2,N2,B2,Bmat2,
     *                            L3,M3,N3,B3,Bmat3,
     *                            gamA13,gamA23,
     *                            gamB13,gamB23,
     *                            xgVEP)
         xgVEP=-1.0d+00*xgVEP


C        xgVEPg1  (1,2) (2,3)
C        --g(e1,p1) V^{ep}(e1,p1) g(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA13=gamma1
         gamB23=gamma2
         call G3Vep_AUX_g13g23V13(I1,J1,K1,A1,Amat1,
     *                            I2,J2,K2,A2,Amat2,
     *                            I3,J3,K3,A3,Amat3,
     *                            L1,M1,N1,B1,Bmat1,
     *                            L2,M2,N2,B2,Bmat2,
     *                            L3,M3,N3,B3,Bmat3,
     *                            gamA13,gamA23,
     *                            gamB13,gamB23,
     *                            xgVEPg1)
         xgVEPg1=-1.0d+00*xgVEPg1

C        xgVEPg2  (2,3) (2,3)
C        --g(e2,p1) V^{ep}(e1,p1) g(e2,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         gamB23=gamma2
         call G3Vep_AUX_g13g23V13(I1,J1,K1,A1,Amat1,
     *                            I2,J2,K2,A2,Amat2,
     *                            I3,J3,K3,A3,Amat3,
     *                            L1,M1,N1,B1,Bmat1,
     *                            L2,M2,N2,B2,Bmat2,
     *                            L3,M3,N3,B3,Bmat3,
     *                            gamA13,gamA23,
     *                            gamB13,gamB23,
     *                            xgVEPg2)
         xgVEPg2=-1.0d+00*xgVEPg2


C        xgVEE  (2,3)
C        --g(e2,p1) V^{ee}(e1,e2)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=gamma1
         call G3Vee_AUX_g23V12(I1,J1,K1,A1,Amat1,
     *                         I2,J2,K2,A2,Amat2,
     *                         I3,J3,K3,A3,Amat3,
     *                         L1,M1,N1,B1,Bmat1,
     *                         L2,M2,N2,B2,Bmat2,
     *                         L3,M3,N3,B3,Bmat3,
     *                         gamA13,gamA23,
     *                         gamB13,gamB23,
     *                         xgVEE)


C        xgVEEg1  (1,3) (1,3)
C        --g(e1,p1) V^{ee}(e1,e2) g(e1,p1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA13=gamma1
         gamB13=gamma2
         call G3Vee_AUX_g13V12(I1,J1,K1,A1,Amat1,
     *                         I2,J2,K2,A2,Amat2,
     *                         I3,J3,K3,A3,Amat3,
     *                         L1,M1,N1,B1,Bmat1,
     *                         L2,M2,N2,B2,Bmat2,
     *                         L3,M3,N3,B3,Bmat3,
     *                         gamA13,gamA23,
     *                         gamB13,gamB23,
     *                         xgVEEg1)



      return
      end

