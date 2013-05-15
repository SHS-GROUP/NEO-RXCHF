C======================================================================
      subroutine RXCHFmult_xcalc_GAM3_MD(I1,J1,K1,A1,Amat1,
     x                                   I2,J2,K2,A2,Amat2,
     x                                   I3,J3,K3,A3,Amat3,
     x                                   I4,J4,K4,A4,Amat4,
     x                                   L1,M1,N1,B1,Bmat1,
     x                                   L2,M2,N2,B2,Bmat2,
     x                                   L3,M3,N3,B3,Bmat3,
     x                                   L4,M4,N4,B4,Bmat4,
     x                                   nat,ngtg1,
     x                                   pmass,cat,zan,
     x                                   bcoef1,gamma1,
     x                                   ans1,ans2)

C Adapted ../gam_3_OMP.f to account for INT_GAM3 terms separately
C======================================================================
      implicit none

C Input Variables
      integer nat
      integer ngtg1
      double precision pmass
      double precision zan(nat)
      double precision cat(3,nat)
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)

      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1,Amat1(3) 
      double precision A2,Amat2(3) 
      double precision A3,Amat3(3) 
      double precision A4,Amat4(3) 
      double precision B1,Bmat1(3) 
      double precision B2,Bmat2(3) 
      double precision B3,Bmat3(3) 
      double precision B4,Bmat4(3) 

C Variables Returned
      double precision ans1   ! XCHF OMG3 contribution
      double precision ans2   ! INT OMG3 contribution

C Local Variables
      integer iii    ! Index for looping over natoms
      integer ik,il  ! Indices for geminal loops
      integer iat

      double precision gamA14
      double precision gamA24
      double precision gamA34
      double precision gamB14
      double precision gamB24
      double precision gamB34
      double precision gamA
      double precision gamB

      double precision cmat(3)
      double precision znuc

      double precision gVEE
      double precision xgVEE
      double precision xx,yy,zz
      double precision zero,half,one,two,four
      parameter(zero=0.0d+00,one=1.0d+00,two=2.0d+00,four=4.0d+00)
      parameter(half=0.5d+00)

      double precision gHEg
      double precision gVEEg1
      double precision gVEEg2
      double precision gVEPg1  ! INT_GAM3 only
      double precision xgHEg
      double precision xgVEEg1 
      double precision xgVEEg2 
      double precision xgVEPg1 ! INT_GAM3 only

      double precision xmass,coulomb_sign
      double precision xke,Vc
      double precision val_vec 
    
      integer NQUAD_coul
      integer NQUAD_ovlap


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     BASIS FUNCTIONS: ASSIGN CENTERS, EXPONENTS, ANGULAR MOM.
C

C     *****ie1 ::  electron 1 bra *****
C     *****ie2 ::  electron 2 bra *****      
C     *****ie3 ::  electron 3 bra *****      
C     *****ip  ::  proton bra     *****

C     *****je1 ::  electron 1 ket *****
C     *****je2 ::  electron 2 ket *****
C     *****je3 ::  electron 3 ket *****
C     *****jp  ::  proton ket     *****

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C      gVEE=zero
C
C      DO IK=1,NGTG1
C
C         gamA=gamma1(ik)
C         gamB=zero
C
CC  xgVEE
CC  --- g(1,p)VEE(2,3)---
C         call G2_MD_xggs(I1,J1,K1,A1,Amat1,
C     1                   I4,J4,K4,A4,Amat4,
C     2                   L1,M1,N1,B1,Bmat1,
C     3                   L4,M4,N4,B4,Bmat4,
C     4                   gamA,gamB,xx)
C
Cc        call pgiovlap(I1,J1,K1,A1,Amat1,
Cc    x                 I4,J4,K4,A4,Amat4,
Cc    x                 L1,M1,N1,B1,Bmat1,
Cc    x                 L4,M4,N4,B4,Bmat4,
Cc    x                 gamA,gamB,xx)
C
C
C         call gfvee(I2,J2,K2,A2,Amat2,
C     x              I3,J3,K3,A3,Amat3,
C     x              L2,M2,N2,B2,Bmat2,
C     x              L3,M3,N3,B3,Bmat3,
C     x              yy)
C
C         call underflow(xx)
C         call underflow(yy)
C
C         xgVEE=xx*yy
C         gVEE=gVEE+(bcoef1(ik)*xgVEE)
C
CC End 1 gamma loop
C      end do
C Begin 2 gamma loop

      gHEg=zero
      gVEEg1=zero
      gVEEg2=zero
      gVEPg1=zero ! INT_GAM3 only

      DO IK=1,NGTG1
         DO IL=1,NGTG1

C>>>>>>>>>>>>>>>>>>>>  xgHEg <<<<<<<<<<<<<<<<<<<<
c           ndim=3
c           natom=nat
c           xmass=one
c           coulomb_sign=-one
c           call o3_Hcore_val(NDIM,NATOM,
c    x                        xmass,zan,cmat,
c    x                        coulomb_sign,
c    x                        I3,J3,K3,A3,Amat3,
c    x                        L3,M3,N3,B3,Bmat3,
c    x                        xx)
CCCCCC-rather than call o3_Hcore_val evaluate 
C everything right here...

            xmass=one
            coulomb_sign=-one

            call gfke(I3,J3,K3,A3,Amat3,
     x                L3,M3,N3,B3,Bmat3,
     x                xmass,xke)
            Vc = ZERO
            do iat=1,nat
                 Cmat(1)=cat(1,iat)
                 Cmat(2)=cat(2,iat)
                 Cmat(3)=cat(3,iat)
                 call gfvec(I3,J3,K3,A3,Amat3,
     x                      L3,M3,N3,B3,Bmat3,
     x                      Cmat,val_vec)
                 Vc = Vc + (zan(iat)*val_vec)
             end do
c            hcore = xke + (Vc*coulomb_sign)
             xx = xke + (Vc*coulomb_sign)
CCCCCC

c           gamA14 = gamma1(ik)
c           gamA24 = ZERO
c           gamA34 = ZERO
c           gamB14 = ZERO
c           gamB24 = gamma1(il)
c           gamB34 = ZERO

            gamA = gamma1(ik)
            gamB = gamma1(il)

c           call G4ovlap_typ1(I1,J1,K1,A1,Amat1,
c           call G3ovlap(I1,J1,K1,A1,Amat1,
c    x                   I2,J2,K2,A2,Amat2,
c    x                   I4,J4,K4,A4,Amat4,
c    x                   L1,M1,N1,B1,Bmat1,
c    x                   L2,M2,N2,B2,Bmat2,
c    x                   L4,M4,N4,B4,Bmat4,
c    x                   ZERO,gamA,ZERO,
c    x                   ZERO,ZERO,gamB,yy)
cc   4                gamA12,gamA13,gamA23,
cc   4                gamB12,gamB13,gamB23,sval)
            call G3_MD_xggs(I1,J1,K1,A1,Amat1,
     x                      I2,J2,K2,A2,Amat2,
     x                      I4,J4,K4,A4,Amat4,
     x                      L1,M1,N1,B1,Bmat1,
     x                      L2,M2,N2,B2,Bmat2,
     x                      L4,M4,N4,B4,Bmat4,
     x                      ZERO,gamA,ZERO,
     x                      ZERO,ZERO,gamB,yy)

C RXCHFmult(  reorder indicies for INT_GAM3 - should not affect XCHF_GAM3
C    index 1: regular electron
C    index 2: special electron 1
C    index 3: special electron 2
C    index 4: proton
C )
C           --g(e2,p1) V^{ep}(e1,p1) g(e3,p1)--
            coulomb_sign  = -ONE
            gamA14 = gamma1(ik)
            gamA24 = ZERO
            gamA34 = ZERO
            gamB14 = ZERO
            gamB24 = gamma1(il)
            gamB34 = ZERO
!           call rys_G4vee_r34(NQUAD_coul,NQUAD_ovlap,
c           call interface_G4vee_r34(NQUAD_coul,NQUAD_ovlap,
c    x                               I1,J1,K1,A1,Amat1,
c    x                               I2,J2,K2,A2,Amat2,
c    x                               I3,J3,K3,A3,Amat3,
c    x                               I4,J4,K4,A4,Amat4,
c    x                               L1,M1,N1,B1,Bmat1,
c    x                               L2,M2,N2,B2,Bmat2,
c    x                               L3,M3,N3,B3,Bmat3,
c    x                               L4,M4,N4,B4,Bmat4,
c    x                               gamA14,gamA24,gamA34,
c    x                               gamB14,gamB24,gamB34,zz)
            call G4_MD_xgVepg(I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        I1,J1,K1,A1,Amat1,
     *                        I4,J4,K4,A4,Amat4,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gamA14,gamB24,
     *                        zz)
c           call G4Vep_AUX_g14g24V34(I1,J1,K1,A1,Amat1,
c    *                               I2,J2,K2,A2,Amat2,
c    *                               I3,J3,K3,A3,Amat3,
c    *                               I4,J4,K4,A4,Amat4,
c    *                               L1,M1,N1,B1,Bmat1,
c    *                               L2,M2,N2,B2,Bmat2,
c    *                               L3,M3,N3,B3,Bmat3,
c    *                               L4,M4,N4,B4,Bmat4,
c    *                               gamA14,zero,
c    *                               zero,gamB24,
c    *                               zz)

             call underflow(xx)
             call underflow(yy)
             call underflow(zz)

             xgHEg = ((xx*yy) + (zz*coulomb_sign) )
           xgVEPg1 = zz*coulomb_sign


C>>>>>>>>>>>>>>>>>>>>  xgVEEg1 <<<<<<<<<<<<<<<<<<<<
C RXCHFmult( reorder indicies for INT_GAM3 - should not affect XCHF_GAM3
C    index 1: regular electron
C    index 2: special electron 1
C    index 3: special electron 2
C    index 4: proton
C )
C  --- g(2,p)VEE(1,3)g(2,p)---
            gamA=gamma1(ik)
            gamB=gamma1(il)
            call G2_MD_xggs(I2,J2,K2,A2,Amat2,
     1                      I4,J4,K4,A4,Amat4,
     2                      L2,M2,N2,B2,Bmat2,
     3                      L4,M4,N4,B4,Bmat4,
     4                      gamA,gamB,xx)

c           call pgiovlap(I1,J1,K1,A1,Amat1,
c    x                    I4,J4,K4,A4,Amat4,
c    x                    L1,M1,N1,B1,Bmat1,
c    x                    L4,M4,N4,B4,Bmat4,
c    x                    gamA,gamB,xx)


            call gfvee(I1,J1,K1,A1,Amat1,
     x                 I3,J3,K3,A3,Amat3,
     x                 L1,M1,N1,B1,Bmat1,
     x                 L3,M3,N3,B3,Bmat3,
     x                 yy)


            call underflow(xx)
            call underflow(yy)

            xgVEEg1=xx*yy

C>>>>>>>>>>>>>>>>>>>>  xgVEEg2 <<<<<<<<<<<<<<<<<<<<
C ans = <GA(1)GA(2)GA(3)GB(4)|g(1,4)g(2,4)g(3,4)/(r1-r2)|GB(1)GB(2)GB(3)GB(4)>
C ans = <ie1 ie2 ie3 ip|g(1,p)g(2,p)g(3,p)/(r1-r2)|je1 je2 je3 jp>
            gamA14=gamma1(ik)
            gamA24=zero
            gamA34=zero
            gamB14=zero
            gamB24=zero
            gamB34=gamma1(il)
CCWS- 11-08-2010(:  There is a problem with the rys_G4Vee_r12
C                   AC routine, so call a CS routine to 
C                   evaluate the xgVEEg2 integral:
c           NQUAD_coul=5
c           NQUAD_ovlap=5
c           call rys_G4Vee_r12(NQUAD_coul,NQUAD_ovlap,
c    x                              I1,J1,K1,A1,Amat1,
c    x                              I2,J2,K2,A2,Amat2,
c    x                              I3,J3,K3,A3,Amat3,
c    x                              I4,J4,K4,A4,Amat4,
c    x                              L1,M1,N1,B1,Bmat1,
c    x                              L2,M2,N2,B2,Bmat2,
c    x                              L3,M3,N3,B3,Bmat3,
c    x                              L4,M4,N4,B4,Bmat4,
c    x                              gamA14,gamA24,gamA34,
c    x                              gamB14,gamB24,gamB34,
c    x                              xgVEEg2)
C RXCHFmult(  reorder indicies for INT_GAM2 - should not affect XCHF_GAM2
C    index 1: regular electron
C    index 2: special electron 1
C    index 3: special electron 2
C    index 4: proton
C  --- g(2,p)VEE(1,2)g(3,p)---
            call G4_MD_xgVeeg(I2,J2,K2,A2,Amat2,
     *                        I1,J1,K1,A1,Amat1,
     *                        I3,J3,K3,A3,Amat3,
     *                        I4,J4,K4,A4,Amat4,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gamA14,gamB34,
     *                        xgVEEg2)
CCWS- 11-08-2010 )
c           call G4Vee_AUX_g14g34V12(I1,J1,K1,A1,Amat1,
c    *                               I2,J2,K2,A2,Amat2,
c    *                               I3,J3,K3,A3,Amat3,
c    *                               I4,J4,K4,A4,Amat4,
c    *                               L1,M1,N1,B1,Bmat1,
c    *                               L2,M2,N2,B2,Bmat2,
c    *                               L3,M3,N3,B3,Bmat3,
c    *                               L4,M4,N4,B4,Bmat4,
c    *                               gamA14,gamB34,
c    *                               xgVEEg2)


C  Sum terms against bcoeff
            gHEg  =gHEg  +(bcoef1(ik)*bcoef1(il)*xgHEg)
            gVEEg1=gVEEg1+(bcoef1(ik)*bcoef1(il)*xgVEEg1)
            gVEEg2=gVEEg2+(bcoef1(ik)*bcoef1(il)*xgVEEg2)
            gVEPg1=gVEPg1+(bcoef1(ik)*bcoef1(il)*xgVEPg1) ! INT_GAM3 only

C End 2 gamma loop
         end do
      end do


C Total integral build
      ans1 = gHEg + gVEEg1*half + four*gVEEg2*half
      ans2 = gVEPg1 + gVEEg1 + two*gVEEg2
CCWS-debug
c                    write(*,*)'========='
c                    write(*,*)'ie1=',ie1
c                    write(*,*)'je1=',je1
c                    write(*,*)'ie2=',ie2
c                    write(*,*)'je2=',je2
c                    write(*,*)'ie3=',ie3
c                    write(*,*)'je3=',je3
c                    write(*,*)'gHEg   =',gHEg
c                    write(*,*)'gVEEg1 =',gVEEg1
c                    write(*,*)'gVEEg2 =',gVEEg2
c                    write(*,*)'ans    =',ans
c                    write(*,*)'========='
CCWS-debug



      return
      end

