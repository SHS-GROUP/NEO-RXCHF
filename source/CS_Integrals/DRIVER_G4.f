C Drivers for calculating select 4-particle (G4)
C integrals used in NEO-XC calculations.
C The routines here are for unit testing the integral codes.
C
C Routines in this file:
C G4_MD_DRIVER
C G4_AC_DRIVER
C G4_RI_DRIVER
C G4_AE_DRIVER
C G4_sonly_DRIVER
C
C======================================================================
      subroutine G4_MD_DRIVER(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        I4,J4,K4,A4,Amat4,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gamma1,gamma2,
     x                        xgVEPg,xgVEEg2)

C
C  Test the following 2 integrals:
C xgVEEg:
C <ie1 ie2 ie3 ip | g(1,p) Vee(1,2) g(3,p) | je1 je2 je3 jp >
C xgVEPg:
C <ie1 ie2 ie3 ip | g(1,p) Vep(3,p) g(2,p) | je1 je2 je3 jp >
C
C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1
      double precision A2
      double precision A3
      double precision A4
      double precision B1
      double precision B2
      double precision B3
      double precision B4

      double precision GAMMA1
      double precision GAMMA2
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Amat4(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)
      double precision Bmat4(3)

C Variables Returned
      double precision xgVEPg,xgVEEg2

C Local variables
c     double precision gamA14
c     double precision gamA24
c     double precision gamA34
c     double precision gamB14
c     double precision gamB24
c     double precision gamB34

      double precision zero
      parameter(zero=0.0d+00)

      double precision gam14
      double precision gam34
      double precision gam24


      gam14=gamma1
      gam34=gamma2
   
C McMurchie-Davidson
      call G4_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                  I2,J2,K2,A2,Amat2,
     *                  I3,J3,K3,A3,Amat3,
     *                  I4,J4,K4,A4,Amat4,
     *                  L1,M1,N1,B1,Bmat1,
     *                  L2,M2,N2,B2,Bmat2,
     *                  L3,M3,N3,B3,Bmat3,
     *                  L4,M4,N4,B4,Bmat4,
     *                  gam14,gam34,
     *                  xgVeeg2)

      gam14=gamma1
      gam24=gamma2
      call G4_MD_xgVepg(I1,J1,K1,A1,Amat1,
     *                  I2,J2,K2,A2,Amat2,
     *                  I3,J3,K3,A3,Amat3,
     *                  I4,J4,K4,A4,Amat4,
     *                  L1,M1,N1,B1,Bmat1,
     *                  L2,M2,N2,B2,Bmat2,
     *                  L3,M3,N3,B3,Bmat3,
     *                  L4,M4,N4,B4,Bmat4,
     *                  gam14,gam24,
     *                  xgVepg)
      xgVepg=-1.0d+00*xgVepg

      return
      end

C======================================================================
      subroutine G4_sonly_DRIVER(I1,J1,K1,A1,Amat1,
     *                           I2,J2,K2,A2,Amat2,
     *                           I3,J3,K3,A3,Amat3,
     *                           I4,J4,K4,A4,Amat4,
     *                           L1,M1,N1,B1,Bmat1,
     *                           L2,M2,N2,B2,Bmat2,
     *                           L3,M3,N3,B3,Bmat3,
     *                           L4,M4,N4,B4,Bmat4,
     *                           gamma1,gamma2,
     x                           xgVeeg_SONLY1,
     x                           xgVeeg_SONLY2,
     x                           xgVepg_SONLY1,
     x                           xgVepg_SONLY2)

C
C  Simple analytical routines for zero total angular momentum
C
C  Test the following 2 integrals:
C xgVEEg:
C <ie1 ie2 ie3 ip | g(1,p) Vee(1,2) g(3,p) | je1 je2 je3 jp >
C xgVEPg:
C <ie1 ie2 ie3 ip | g(1,p) Vep(3,p) g(2,p) | je1 je2 je3 jp >
C
C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1
      double precision A2
      double precision A3
      double precision A4
      double precision B1
      double precision B2
      double precision B3
      double precision B4

      double precision GAMMA1
      double precision GAMMA2
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Amat4(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)
      double precision Bmat4(3)

C Variables Returned
      double precision xgVeeg_SONLY1 
      double precision xgVeeg_SONLY2 
      double precision xgVepg_SONLY1 
      double precision xgVepg_SONLY2 

C Local variables
c     double precision gamA14
c     double precision gamA24
c     double precision gamA34
c     double precision gamB14
c     double precision gamB24
c     double precision gamB34
      double precision zero
      parameter(zero=0.0d+00)
      double precision gam14
      double precision gam34
      double precision gam24



      gam14=gamma1
      gam34=gamma2
      gam24=gamma2

      call cws_G4_xgVeeg_SONLY1(I1,J1,K1,A1,Amat1,
     *                          I2,J2,K2,A2,Amat2,
     *                          I3,J3,K3,A3,Amat3,
     *                          I4,J4,K4,A4,Amat4,
     *                          L1,M1,N1,B1,Bmat1,
     *                          L2,M2,N2,B2,Bmat2,
     *                          L3,M3,N3,B3,Bmat3,
     *                          L4,M4,N4,B4,Bmat4,
     *                          gam14,gam34,
     *                          xgVeeg_SONLY1)

      call cws_G4_xgVeeg_SONLY2(I1,J1,K1,A1,Amat1,
     *                          I2,J2,K2,A2,Amat2,
     *                          I3,J3,K3,A3,Amat3,
     *                          I4,J4,K4,A4,Amat4,
     *                          L1,M1,N1,B1,Bmat1,
     *                          L2,M2,N2,B2,Bmat2,
     *                          L3,M3,N3,B3,Bmat3,
     *                          L4,M4,N4,B4,Bmat4,
     *                          gam14,gam34,
     *                          xgVeeg_SONLY2)

      call cws_G4_xgVepg_SONLY1(I1,J1,K1,A1,Amat1,
     *                          I2,J2,K2,A2,Amat2,
     *                          I3,J3,K3,A3,Amat3,
     *                          I4,J4,K4,A4,Amat4,
     *                          L1,M1,N1,B1,Bmat1,
     *                          L2,M2,N2,B2,Bmat2,
     *                          L3,M3,N3,B3,Bmat3,
     *                          L4,M4,N4,B4,Bmat4,
     *                          gam14,gam24,
     *                          xgVepg_SONLY1)
      xgVepg_SONLY1=-1.0d+00*xgVepg_SONLY1

      call cws_G4_xgVepg_SONLY2(I1,J1,K1,A1,Amat1,
     *                          I2,J2,K2,A2,Amat2,
     *                          I3,J3,K3,A3,Amat3,
     *                          I4,J4,K4,A4,Amat4,
     *                          L1,M1,N1,B1,Bmat1,
     *                          L2,M2,N2,B2,Bmat2,
     *                          L3,M3,N3,B3,Bmat3,
     *                          L4,M4,N4,B4,Bmat4,
     *                          gam14,gam24,
     *                          xgVepg_SONLY2)
      xgVepg_SONLY2=-1.0d+00*xgVepg_SONLY2


      return
      end

C======================================================================
      subroutine G4_AE_DRIVER(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        I4,J4,K4,A4,Amat4,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gamma1,gamma2,
     x                        xgVEEg2,xgVEPg)

C
C  Test the following 2 integrals:
C xgVEEg:
C <ie1 ie2 ie3 ip | g(1,p) Vee(1,2) g(3,p) | je1 je2 je3 jp >
C xgVEPg:
C <ie1 ie2 ie3 ip | g(1,p) Vep(3,p) g(2,p) | je1 je2 je3 jp >
C
C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1
      double precision A2
      double precision A3
      double precision A4
      double precision B1
      double precision B2
      double precision B3
      double precision B4

      double precision GAMMA1
      double precision GAMMA2
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Amat4(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)
      double precision Bmat4(3)

C Local variables
      double precision gamA14
      double precision gamA24
      double precision gamA34
      double precision gamB14
      double precision gamB24
      double precision gamB34

      double precision xgVEPg
      double precision xgVEEg2

      double precision zero
      parameter(zero=0.0d+00)


      gamA14=gamma1
      gamA24=zero
      gamA34=zero
      gamB14=zero
      gamB24=gamma2
      gamB34=gamma2


C Evaluates the following Vep integral:
C xgVEPg  ::  gA(1,4) V(3,4) gB(2,4)
C Note:  Proton is assumed to be particle 4

      call G4Vep_AUX_g14g24V34(I1,J1,K1,A1,Amat1,
     *                         I2,J2,K2,A2,Amat2,
     *                         I3,J3,K3,A3,Amat3,
     *                         I4,J4,K4,A4,Amat4,
     *                         L1,M1,N1,B1,Bmat1,
     *                         L2,M2,N2,B2,Bmat2,
     *                         L3,M3,N3,B3,Bmat3,
     *                         L4,M4,N4,B4,Bmat4,
     *                         gamA14,gamA24,
     *                         gamB14,gamB24,
     *                         xgVEPg)
      xgVEPg=-1.0d+00*xgVEPg

C Evaluates the following Vep integral:
C xgVEPg  ::  gA(1,4) V(1,2) gB(3,4)
C Note:  Proton is assumed to be particle 4

      call G4Vee_AUX_g14g34V12(I1,J1,K1,A1,Amat1,
     *                         I2,J2,K2,A2,Amat2,
     *                         I3,J3,K3,A3,Amat3,
     *                         I4,J4,K4,A4,Amat4,
     *                         L1,M1,N1,B1,Bmat1,
     *                         L2,M2,N2,B2,Bmat2,
     *                         L3,M3,N3,B3,Bmat3,
     *                         L4,M4,N4,B4,Bmat4,
     *                         gamA14,gamB34,
     *                         xgVeeg2)



      return
      end

C======================================================================
      subroutine G4_AC_DRIVER(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        I4,J4,K4,A4,Amat4,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gamma1,gamma2,
     *                        xgVEPg,xgVEEg2)

C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1
      double precision A2
      double precision A3
      double precision A4
      double precision B1
      double precision B2
      double precision B3
      double precision B4

      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Amat4(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)
      double precision Bmat4(3)

      double precision gamma1
      double precision gamma2

C Variables Returned
      double precision xgVEPg
      double precision xgVEEg2

C Local variables
      integer NQUAD_coul
      integer NQUAD_ovlap

      double precision gamA14
      double precision gamA24
      double precision gamA34
      double precision gamB14
      double precision gamB24
      double precision gamB34

      double precision coulomb_sign
      double precision one
      parameter (one=1.0d+00)
      double precision zero
      parameter (zero=0.0d+00)


      NQUAD_coul=5
      NQUAD_ovlap=5

      coulomb_sign  = -ONE
      gamA14 = gamma1
      gamA24 = ZERO
      gamA34 = ZERO
      gamB14 = ZERO
      gamB24 = gamma2
      gamB34 = ZERO
c     write(*,*)'>>>call interface_G4vee_r34'
      call interface_G4vee_r34(NQUAD_coul,NQUAD_ovlap,
     x                              I1,J1,K1,A1,Amat1,
     x                              I2,J2,K2,A2,Amat2,
     x                              I3,J3,K3,A3,Amat3,
     x                              I4,J4,K4,A4,Amat4,
     x                              L1,M1,N1,B1,Bmat1,
     x                              L2,M2,N2,B2,Bmat2,
     x                              L3,M3,N3,B3,Bmat3,
     x                              L4,M4,N4,B4,Bmat4,
     x                           gamA14,gamA24,gamA34,
     x                           gamB14,gamB24,gamB34,xgVEPg)

      call underflow(xgVEPg)
      xgVEPg=xgVEPg*coulomb_sign



C>>>>>>>>>>>>>>>>>>>>  xgVEEg2 <<<<<<<<<<<<<<<<<<<<
C ans = <GA(1)GA(2)GA(3)GB(4)|g(1,4)g(2,4)g(3,4)/(r1-r2)|GB(1)GB(2)GB(3)GB(4)>
C ans = <ie1 ie2 ie3 ip|g(1,p)g(2,p)g(3,p)/(r1-r2)|je1 je2 je3 jp>
c     write(*,*)'doing xgVEEg2'
      gamA14=gamma1
      gamA24=zero
      gamA34=zero
      gamB14=zero
      gamB24=zero
      gamB34=gamma2
c     write(*,*)'>>>call rys_G4Vee_r12'
      call rys_G4Vee_r12(NQUAD_coul,NQUAD_ovlap,
     x                        I1,J1,K1,A1,Amat1,
     x                        I2,J2,K2,A2,Amat2,
     x                        I3,J3,K3,A3,Amat3,
     x                        I4,J4,K4,A4,Amat4,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L2,M2,N2,B2,Bmat2,
     x                        L3,M3,N3,B3,Bmat3,
     x                        L4,M4,N4,B4,Bmat4,
     x                        gamA14,gamA24,gamA34,
     x                        gamB14,gamB24,gamB34,
     x                        xgVEEg2)


      return
      end

C======================================================================
c     subroutine RI_G4(orth_abs,
      subroutine G4_RI_DRIVER(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        I4,J4,K4,A4,Amat4,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gamma1,gamma2,
     *                        xgVEPg,xgVEEg2)

C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1
      double precision A2
      double precision A3
      double precision A4
      double precision B1
      double precision B2
      double precision B3
      double precision B4

      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Amat4(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)
      double precision Bmat4(3)

      double precision GAMMA1
      double precision GAMMA2

C Variables Returned
      double precision xgVEPg
      double precision xgVEEg2

C Local variables
      double precision gamA14
      double precision gamA24
      double precision gamA34
      double precision gamB14
      double precision gamB24
      double precision gamB34

      double precision coulomb_sign

      double precision one
      parameter (one=1.0d+00)

      double precision zero
      parameter (zero=0.0d+00)


      gamA14 = gamma1
      gamA24 = ZERO
      gamA34 = ZERO
      gamB14 = ZERO
      gamB24 = gamma2
      gamB34 = ZERO
      call RI_G4_xgVEPg(I1,J1,K1,A1,Amat1,
     x                  I2,J2,K2,A2,Amat2,
     x                  I3,J3,K3,A3,Amat3,
     x                  I4,J4,K4,A4,Amat4,
     x                  L1,M1,N1,B1,Bmat1,
     x                  L2,M2,N2,B2,Bmat2,
     x                  L3,M3,N3,B3,Bmat3,
     x                  L4,M4,N4,B4,Bmat4,
     x                  gamA14,gamA24,gamA34,
     x                  gamB14,gamB24,gamB34,
     x                  xgVEPg)

      call underflow(xgVEPg)


C>>>>>>>>>>>>>>>>>>>>  xgVEEg2 <<<<<<<<<<<<<<<<<<<<
C ans = <GA(1)GA(2)GA(3)GB(4)|g(1,4)g(2,4)g(3,4)/(r1-r2)|GB(1)GB(2)GB(3)GB(4)>
C ans = <ie1 ie2 ie3 ip|g(1,p)g(2,p)g(3,p)/(r1-r2)|je1 je2 je3 jp>
      gamA14=gamma1
      gamA24=zero
      gamA34=zero
      gamB14=zero
      gamB24=zero
      gamB34=gamma2
      call RI_G4_xgVEEg2(I1,J1,K1,A1,Amat1,
     x                   I2,J2,K2,A2,Amat2,
     x                   I3,J3,K3,A3,Amat3,
     x                   I4,J4,K4,A4,Amat4,
     x                   L1,M1,N1,B1,Bmat1,
     x                   L2,M2,N2,B2,Bmat2,
     x                   L3,M3,N3,B3,Bmat3,
     x                   L4,M4,N4,B4,Bmat4,
     x                   gamA14,gamA24,gamA34,
     x                   gamB14,gamB24,gamB34,
     x                   xgVEEg2)

      call underflow(xgVEEg2)


      return
      end

