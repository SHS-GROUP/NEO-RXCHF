C======================================================================
      subroutine RXCHFmult_xcalc_GAM2_MD(I1,J1,K1,A1,Amat1,
     x                                   I2,J2,K2,A2,Amat2,
     x                                   I3,J3,K3,A3,Amat3,
     x                                   L1,M1,N1,B1,Bmat1,
     x                                   L2,M2,N2,B2,Bmat2,
     x                                   L3,M3,N3,B3,Bmat3,
     x                                   nat,ngtg1,
     x                                   pmass,cat,zan,
     x                                   bcoef1,gamma1,
     x                                   ansE1,ansE2,ansS)

C Adapted ../gam_2_OMP.f to account for INT_GAM2 terms separately
C======================================================================
      implicit none

C Input Variables
c     integer ie1,ie2,ip
c     integer je1,je2,jp
      integer nat
      integer ngtg1
      double precision pmass
      double precision zan(nat)
      double precision cat(3,nat)
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)

C Variables Returned
      double precision ansE1   ! XCHF OMG2 contribution
      double precision ansE2   ! INT OMG2 contribution
      double precision ansS    ! XCHF OMG2s contribution

C Local Variables
      integer iii    ! Index for looping over natoms
      integer ik,il  ! Indices for geminal loops

      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3

      double precision A1,Amat1(3) 
      double precision A2,Amat2(3) 
      double precision A3,Amat3(3) 
      double precision B1,Bmat1(3) 
      double precision B2,Bmat2(3) 
      double precision B3,Bmat3(3) 

      double precision cmat(3)
      double precision znuc
      double precision vec_ans
      double precision xmass 

      double precision gamA12,gamA13,gamA23 
      double precision gamB12,gamB13,gamB23 

      double precision zero,half,two,four
      parameter(zero=0.0d+00,half=0.5d+00,two=2.0d+00,four=4.0d+00)

      double precision xgTE      
      double precision xgVEC     
      double precision xgVEE     
      double precision xgVEP     

      double precision gTE      
      double precision gVEC     
      double precision gVEE     
      double precision gVEP     

      double precision xgTPg      
      double precision xgTEg1     
      double precision xgTEg2     
      double precision xgTEg3     
      double precision xgVPCg      
      double precision xgVECg1     
      double precision xgVECg2     
      double precision xgVEPg1     
      double precision xgVEPg2     
      double precision xgVEEg1     
      double precision xgVEEg2     
      double precision xgsg     

      double precision gTPg      
      double precision gTEg1     
      double precision gTEg2     
      double precision gTEg3     
      double precision gVPCg      
      double precision gVECg1     
      double precision gVECg2     
      double precision gVEPg1     
      double precision gVEPg2     
      double precision gVEEg1     
      double precision gVEEg2     
      double precision gsg     

C******
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     BASIS FUNCTIONS: ASSIGN CENTERS, EXPONENTS, ANGULAR MOM.
C

C     *****ie1 ::  electron 1 bra *****
C     *****ie2 ::  electron 2 bra *****      
C     *****ip  ::  proton bra     *****

C     *****je1 ::  electron 1 ket *****
C     *****je2 ::  electron 2 ket *****
C     *****jp  ::  proton ket     *****

c     call get_BF(1,ie1,I1,J1,K1,A1,Amat1)
c     call get_BF(1,ie2,I2,J2,K2,A2,Amat2)
c     call get_BF(2,ip,I3,J3,K3,A3,Amat3)

c     call get_BF(1,je1,L1,M1,N1,B1,Bmat1)
c     call get_BF(1,je2,L2,M2,N2,B2,Bmat2)
c     call get_BF(2,jp,L3,M3,N3,B3,Bmat3)


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C      gTE =zero
C      gVEC=zero
C      gVEE=zero
C      gVEP=zero
C
C      DO IK=1,NGTG1
C
CC        xgTE  (2,3)
CC        --g(e2,p1) T^e(e1)--
C         call zero_gam(gamA12,gamA13,gamA23,
C     *                 gamB12,gamB13,gamB23)
C         gamA23=GAMMA1(IK)
C         xmass=1.0d+00
C         call G3_MD_KE(I1,J1,K1,A1,Amat1,
C     *                 I2,J2,K2,A2,Amat2,
C     *                 I3,J3,K3,A3,Amat3,
C     *                 L1,M1,N1,B1,Bmat1,
C     *                 L2,M2,N2,B2,Bmat2,
C     *                 L3,M3,N3,B3,Bmat3,
C     *                 gamA12,gamA13,gamA23,
C     *                 gamB12,gamB13,gamB23,
C     *                 xmass,xgTE)
C
CC        xgVEC  (2,3)
CC        --g(e2,p1) V^{eC}(e1)--
C         xgVEC=zero
C         DO III=1,NAT
C            Cmat(1)=cat(1,III)
C            Cmat(2)=cat(2,III)
C            Cmat(3)=cat(3,III)
C            ZNUC=ZAN(III)
C 
C            call G3_MD_xgVxCg(I1,J1,K1,A1,Amat1,
C     *                        I2,J2,K2,A2,Amat2,
C     *                        I3,J3,K3,A3,Amat3,
C     *                        L1,M1,N1,B1,Bmat1,
C     *                        L2,M2,N2,B2,Bmat2,
C     *                        L3,M3,N3,B3,Bmat3,
C     *                        gamA12,gamA13,gamA23,
C     *                        gamB12,gamB13,gamB23,
C     *                        Cmat,ZNUC,
C     *                        vec_ans)
C
C            call underflow(vec_ans)
C            xgVEC=xgVEC-ZNUC*vec_ans
C 
C         END DO
C
CC        xgVEP  (2,3)
CC        --g(e2,p1) V^{ep}(e1,p1)--
C         call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
C     *                     I3,J3,K3,A3,Amat3,
C     *                     I2,J2,K2,A2,Amat2,
C     *                     L1,M1,N1,B1,Bmat1,
C     *                     L3,M3,N3,B3,Bmat3,
C     *                     L2,M2,N2,B2,Bmat2,
C     *                     gamA12,gamA13,gamA23,
C     *                     gamB12,gamB13,gamB23,
C     *                     xgVEP)
C
C         call underflow(xgVEP)
C         xgVEP = -xgVEP
C
CC        xgVEE  (2,3)
CC        --g(e2,p1) V^{ee}(e1,e2)--
C         call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
C     *                     I2,J2,K2,A2,Amat2,
C     *                     I3,J3,K3,A3,Amat3,
C     *                     L1,M1,N1,B1,Bmat1,
C     *                     L2,M2,N2,B2,Bmat2,
C     *                     L3,M3,N3,B3,Bmat3,
C     *                     gamA12,gamA13,gamA23,
C     *                     gamB12,gamB13,gamB23,
C     *                     xgVEE)
C         call underflow(xgVEE)
C
C
CC  Sum up 1 gamma loop terms
C         gTE= gTE +bcoef1(ik)*xgTE
C         gVEC=gVEC+bcoef1(ik)*xgVEC
C         gVEP=gVEP+bcoef1(ik)*xgVEP
C         gVEE=gVEE+bcoef1(ik)*xgVEE
C
C      end do
CC --- end of 1 gamma loop ---
C
C --- 2 gamma loop ---
C
      gTPg =zero
      gTEg1=zero
      gTEg2=zero
      gTEg3=zero

      gVPCg =zero
      gVECg1=zero
      gVECg2=zero

      gVEPg1=zero
      gVEPg2=zero
      gVEEg1=zero
      gVEEg2=zero

      gsg = zero

      DO IK=1,NGTG1
         DO IL=1,NGTG1

C%%%%%%%%%%--2 Gamma Kinetic Energy Integrals--%%%%%%%%%%

C           xgTPg  
C           --gA(e1,p1) T^p(p) gB(e2,p1)--
C           --gA(1,3) T^p(3) gB(2,3)--
C           --gA(3,1) T^p(1) gB(2,1)--
C           --gA(2,1) T^p(1) gB(3,1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
            gamA12=gamma1(IK)
            gamB13=gamma1(IL)
            xmass=pmass
            call G3_MD_KE(I3,J3,K3,A3,Amat3,
     x                    I1,J1,K1,A1,Amat1,
     *                    I2,J2,K2,A2,Amat2,
     *                    L3,M3,N3,B3,Bmat3,
     *                    L1,M1,N1,B1,Bmat1,
     *                    L2,M2,N2,B2,Bmat2,
     *                    gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23,
     *                    xmass,xgTPg)

C             xgTEg1 (1,3) (2,3)
C-------------g(e1,p1) T^e(e1) g(e2,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     1                    gamB12,gamB13,gamB23)
            gamA13=GAMMA1(IK)
            gamB23=GAMMA1(IL)
            xmass=1.0d+00
            call G3_MD_KE(I1,J1,K1,A1,Amat1,
     *                    I2,J2,K2,A2,Amat2,
     *                    I3,J3,K3,A3,Amat3,
     *                    L1,M1,N1,B1,Bmat1,
     *                    L2,M2,N2,B2,Bmat2,
     *                    L3,M3,N3,B3,Bmat3,
     *                    gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23,
     *                    xmass,xgTEg1)


C             xgTEg2 (2,3) (1,3)
C-------------g(e2,p1) T^e(e1) g(e1,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23)
            gamA23=GAMMA1(IK)
            gamB13=GAMMA1(IL)
            xmass=1.0d+00
            call G3_MD_KE(I1,J1,K1,A1,Amat1,
     *                    I2,J2,K2,A2,Amat2,
     *                    I3,J3,K3,A3,Amat3,
     *                    L1,M1,N1,B1,Bmat1,
     *                    L2,M2,N2,B2,Bmat2,
     *                    L3,M3,N3,B3,Bmat3,
     *                    gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23,
     *                    xmass,xgTEg2)

C             xgTEg3 (2,3) (2,3)
C-------------g(e2,p1)T^e(e1)g(e2,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23)
            gamA23=GAMMA1(IK)
            gamB23=GAMMA1(IL)
            xmass=1.0d+00
            call G3_MD_KE(I1,J1,K1,A1,Amat1,
     *                    I2,J2,K2,A2,Amat2,
     *                    I3,J3,K3,A3,Amat3,
     *                    L1,M1,N1,B1,Bmat1,
     *                    L2,M2,N2,B2,Bmat2,
     *                    L3,M3,N3,B3,Bmat3,
     *                    gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23,
     *                    xmass,xgTEg3)


C%%%%%%%%%%--End of 2 Gamma Kinetic Energy Integrals--%%%%%%%%%%
C
C%%%%%%%%%%--2 Gamma VEC Integrals--%%%%%%%%%%

            xgVPCg =zero
            xgVECg1=zero
            xgVECg2=zero

            DO III=1,NAT
               Cmat(1)=cat(1,III)
               Cmat(2)=cat(2,III)
               Cmat(3)=cat(3,III)
               ZNUC=ZAN(III)
C
C              xgVPCg (1,2)  (1,3)
C              --g(e1,p1) V^{PC}(p1) g(e2,p1)--
C
               call zero_gam(gamA12,gamA13,gamA23,
     *                       gamB12,gamB13,gamB23)
               gamA12=GAMMA1(IK)
               gamB13=GAMMA1(IL)

               call G3_MD_xgVxCg(I3,J3,K3,A3,Amat3,
     x                           I1,J1,K1,A1,Amat1,
     x                           I2,J2,K2,A2,Amat2,
     x                           L3,M3,N3,B3,Bmat3,
     x                           L1,M1,N1,B1,Bmat1,
     x                           L2,M2,N2,B2,Bmat2,
     x                           gamA12,gamA13,gamA23,
     x                           gamB12,gamB13,gamB23,
     x                           Cmat,ZNUC,
     x                           vec_ans)

               call underflow(vec_ans)
               xgVPCg=xgVPCg+ZNUC*vec_ans

C              xgVECg1 (1,3)  (2,3)
C              --g(e1,p1) V^{eC}(e1) g(e2,p1)--
C
               call zero_gam(gamA12,gamA13,gamA23,
     *                       gamB12,gamB13,gamB23)
               gamA13=GAMMA1(IK)
               gamB23=GAMMA1(IL)

               call G3_MD_xgVxCg(I1,J1,K1,A1,Amat1,
     *                           I2,J2,K2,A2,Amat2,
     *                           I3,J3,K3,A3,Amat3,
     *                           L1,M1,N1,B1,Bmat1,
     *                           L2,M2,N2,B2,Bmat2,
     *                           L3,M3,N3,B3,Bmat3,
     *                           gamA12,gamA13,gamA23,
     *                           gamB12,gamB13,gamB23,
     *                           Cmat,ZNUC,
     *                           vec_ans)
               call underflow(vec_ans)
               xgVECg1=xgVECg1-ZNUC*vec_ans

C              xgVECg2 (2,3)  (2,3)
C              --g(e2,p1) V^{eC}(e1) g(e2,p1)--
               call zero_gam(gamA12,gamA13,gamA23,
     *                       gamB12,gamB13,gamB23)
               gamA23=GAMMA1(IK)
               gamB23=GAMMA1(IL)

               call G3_MD_xgVxCg(I1,J1,K1,A1,Amat1,
     *                           I2,J2,K2,A2,Amat2,
     *                           I3,J3,K3,A3,Amat3,
     *                           L1,M1,N1,B1,Bmat1,
     *                           L2,M2,N2,B2,Bmat2,
     *                           L3,M3,N3,B3,Bmat3,
     *                           gamA12,gamA13,gamA23,
     *                           gamB12,gamB13,gamB23,
     *                           Cmat,ZNUC,
     *                           vec_ans)
               call underflow(vec_ans)

               xgVECg2=xgVECg2-ZNUC*vec_ans
C End of loop over NAT (classical nuclei)
            END DO

C%%%%%%%%%%--End of 2 Gamma VEC Integrals--%%%%%%%%%%

C%%%%%%%%%%--2 Gamma VEP Integrals--%%%%%%%%%%

C           xgVEPg1  (1,2) (2,3)
C           --g(e1,p1) V^{ep}(e1,p1) g(e2,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23)
            gamA12=GAMMA1(IK)
            gamB23=GAMMA1(IL)

            call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                        I3,J3,K3,A3,Amat3,
     *                        I2,J2,K2,A2,Amat2,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L2,M2,N2,B2,Bmat2,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        xgVEPg1)
            call underflow(xgVEPg1)
            xgVEPg1=-xgVEPg1

C RXCHFmult( reorder indicies for INT_GAM2 - should not affect XCHF_GAM2
C    index 1: regular electron
C    index 2: special electron
C    index 3: proton
C )
C           --g(e2,p1) V^{ep}(e1,p1) g(e2,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23)
            gamA23=GAMMA1(IK)
            gamB23=GAMMA1(IL)

            call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                        I3,J3,K3,A3,Amat3,
     *                        I2,J2,K2,A2,Amat2,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L2,M2,N2,B2,Bmat2,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        xgVEPg2)
            call underflow(xgVEPg2)
            xgVEPg2=-xgVEPg2

C%%%%%%%%%%--End of 2 Gamma VEP Integrals--%%%%%%%%%%

C%%%%%%%%%%--2 Gamma VEE Integrals--%%%%%%%%%%

C RXCHFmult(  reorder indicies for INT_GAM2 - should not affect XCHF_GAM2
C    index 1: regular electron
C    index 2: special electron
C    index 3: proton
C )
C           --g(e2,p1) V^{ee}(e1,e2) g(e2,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23)
            gamA23=GAMMA1(IK)
            gamB23=GAMMA1(IL)

            call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        xgVEEg1)
            call underflow(xgVEEg1)

C           xgVEEg2  (1,3) (2,3)
C           --g(e1,p1) V^{ee}(e1,e2) g(e2,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     1                    gamB12,gamB13,gamB23)
            gamA13=GAMMA1(IK)
            gamB23=GAMMA1(IL)

            call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        xgVEEg2)
            call underflow(xgVEEg2)

c           evaluate overlap
c           --g(re1,rp1)g(re2,rp1)--         
            call zero_gam(gamA12,gamA13,gamA23, 
     *            gamB12,gamB13,gamB23)
            gamA13 = gamma1(ik)
            gamB23 = gamma1(il)
c           call G3ovlap(I1,J1,K1,A1,Amat1,
c    *                   I2,J2,K2,A2,Amat2,
c    *                   I3,J3,K3,A3,Amat3,
c    *                   L1,M1,N1,B1,Bmat1,
c    *                   L2,M2,N2,B2,Bmat2,
c    *                   L3,M3,N3,B3,Bmat3,
c    *                   gamA12,gamA13,gamA23,
c    *                   gamB12,gamB13,gamB23,
c    *                   xgsg)
            call G3_MD_xggs(I1,J1,K1,A1,Amat1,
     *                      I2,J2,K2,A2,Amat2,
     *                      I3,J3,K3,A3,Amat3,
     *                      L1,M1,N1,B1,Bmat1,
     *                      L2,M2,N2,B2,Bmat2,
     *                      L3,M3,N3,B3,Bmat3,
     *                      gamA12,gamA13,gamA23,
     *                      gamB12,gamB13,gamB23,
     *                      xgsg)


C  Sum up 2 gamma loop terms
            gTPg =gTPg +bcoef1(ik)*bcoef1(iL)*xgTPg
            gTEg1=gTEg1+bcoef1(ik)*bcoef1(iL)*xgTEg1
            gTEg2=gTEg2+bcoef1(ik)*bcoef1(iL)*xgTEg2
            gTEg3=gTEg3+bcoef1(ik)*bcoef1(iL)*xgTEg3

            gVPCg =gVPCg +bcoef1(ik)*bcoef1(iL)*xgVPCg
            gVECg1=gVECg1+bcoef1(ik)*bcoef1(iL)*xgVECg1
            gVECg2=gVECg2+bcoef1(ik)*bcoef1(iL)*xgVECg2

            gVEPg1=gVEPg1+bcoef1(ik)*bcoef1(iL)*xgVEPg1
            gVEPg2=gVEPg2+bcoef1(ik)*bcoef1(iL)*xgVEPg2

            gVEEg1=gVEEg1+bcoef1(ik)*bcoef1(iL)*xgVEEg1
            gVEEg2=gVEEg2+bcoef1(ik)*bcoef1(iL)*xgVEEg2

            gsg=gsg+bcoef1(ik)*bcoef1(iL)*xgsg

C%%%%%%%%%%--End of 2 Gamma VEE Integrals--%%%%%%%%%%
      end do
         end do
C --- end of 2 gamma loop ---

C sum of terms: 
      ansE1 = gTPg + gVPCg
     x      + gTEg1 + gTEg3
     x      + two * (gVEPg1 + gVECg1)
     x      + gTEg2 + gVEPg2 + gVECg2
     x      + two*gVEEg1*half
     x      + two*gVEEg2*half

      ansE2 = gVEPg2 + gVEEg1

      ansS=gsg

      return
      end

