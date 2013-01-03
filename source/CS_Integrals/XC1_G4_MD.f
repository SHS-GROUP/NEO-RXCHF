C
C These routines calculate the G4 integrals for s-functions
C only.  They are for testing/debugging purposes:
C cws_G4_xgVeeg_SONLY1  Applies Boys Lemma 2 times, calls 2-particle
C cws_G4_xgVeeg_SONLY2  Applies Boys Lemma 1 time, calls 3-particle
C cws_G4_xgVepg_SONLY1  Applies Boys Lemma 2 times, calls 2-particle
C cws_G4_xgVepg_SONLY2  Applies Boys Lemma 1 time, calls 3-particle
C
C These routines calculate the G4 integrals for any angular momentum:
C G4_MD_xgVeeg
C G4_MD_xgVepg
C
C=======================================================================
      subroutine cws_G4_xgVeeg_SONLY1(I1,J1,K1,A1,Amat1,
     *                                I2,J2,K2,A2,Amat2,
     *                                I3,J3,K3,A3,Amat3,
     *                                I4,J4,K4,A4,Amat4,
     *                                L1,M1,N1,B1,Bmat1,
     *                                L2,M2,N2,B2,Bmat2,
     *                                L3,M3,N3,B3,Bmat3,
     *                                L4,M4,N4,B4,Bmat4,
     *                                gam14,gam34,
     *                                xgVeeg)


C=======================================================================
      implicit none

C Input Variables
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

      double precision gam14
      double precision gam34

C Variables Returned
      double precision xgVeeg

C Local Variables
      double precision KAB1(3)
      double precision KAB2(3)
      double precision KAB3(3)
      double precision KAB4(3)

      double precision P1
      double precision P2
      double precision P3
      double precision P4

      double precision Pmat1(3)
      double precision Pmat2(3)
      double precision Pmat3(3)
      double precision Pmat4(3)

      double precision qAB1
      double precision qAB2
      double precision qAB3
      double precision qAB4

      double precision Qmat1(3)
      double precision Qmat2(3)
      double precision Qmat3(3)
      double precision Qmat4(3)

      double precision c0
      double precision c1
      double precision c3
      double precision c13

      double precision KP1P4(3)
      double precision KP3P4(3)
      double precision qPC1
      double precision qPC3
      double precision S1
      double precision S3
      double precision Smat1(3)
      double precision Smat3(3)
      double precision SQmat1(3)
      double precision SQmat3(3)

      double precision CB0
      double precision CB1
      double precision PI
      parameter(PI=3.14159265358979d+00)
      double precision half
      parameter(half=0.5d+00)

      double precision KS1S3(3)
      double precision qSCB1
      double precision W1
      double precision Wmat13(3)
      double precision WQmat13(3)

      double precision C4
      double precision F0
      double precision alpha
      double precision expT
      double precision R2WP
      double precision RWP(3)
      double precision Cof1
      double precision Cof2
      double precision Cof3
      double precision Cof4
CCWS CHECK(
      double precision xans,halfW1,halfP2,CHK,CHKANS
CCWS CHECK)


C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
      call Gauss_prod(A3,B3,Amat3,Bmat3,KAB3,qAB3,P3,Pmat3,Qmat3)
      call Gauss_prod(A4,B4,Amat4,Bmat4,KAB4,qAB4,P4,Pmat4,Qmat4)
C  Apply Boys' Lemma 2 to particle 4:
C     (
C IN: |dx4 Exp(-P4(x4-Pmat4)**2) Exp(-g14(x1-x4)**2) Exp(-g34(x3-x4)**2)
C     )
C OUT: C0*Exp(-C1(x1-Pmat4)**2 Exp(-C3(x3-Pmat4)**2 Exp(-C13(x1-x3)**2
c     C1=(P4*gam14)/xsum
c     C3=(P4*gam34)/xsum
c     C13=(gam14*gam34)/xsum

      call SFBoys2(P4,gam14,gam34,Pmat4,
     x             C0,C1,C3,C13)

C  Should now have the following integral:
C  (
C  | Exp(-P1(x1-Pmat1)**2) Exp(-P2(x2-Pmat2)**2) Exp(-P3(x3-Pmat3)**2)  
C  | Exp(-C1(x1-Pmat4)**2)                       Exp(-C3(x3-Pmat4)**2)  
C  | Exp(-C13(x1-x3)**2) 1/|x1-x2| 
C  )

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(P1,C1,Pmat1,Pmat4,KP1P4,qPC1,S1,Smat1,SQmat1)
      call Gauss_prod(P3,C3,Pmat3,Pmat4,KP3P4,qPC3,S3,Smat3,SQmat3)

C  Should now have the following integral:
C  (        
C  | Exp(-S1(x1-Smat1)**2) Exp(-P2(x2-Pmat2)**2 Exp(-S3(x3-Smat3)**2)  
C  | Exp(-c13(x1-x3)**2) 1/|x1-x2| 
C  )

C  Apply Boys' Lemma 2 to particle 3:
C     (
C IN: |dx3 Exp(-S3(x3-Smat3)**2) Exp(-c13(x1-x3)**2) 
C     )
C OUT: CB0*Exp(-CB1(x1-Smat3)**2 
c     CB0=( PI / (s3 + c13) )**(3/2)
c     CB1=( S3 C13 ) / ( S3 + C13 )

      CB0=(PI/(S3+C13))*sqrt(PI/(S3+C13))
      CB1=(S3*C13)/(S3+C13)

C  Should now have the following 2-particle integral:
C  (        
C  | Exp(-S1(x1-Smat1)**2) Exp(-P2(x2-Pmat2)**2 Exp(-CB1(x1-Smat3)**2)  
C  | 1/|x1-x2| 
C  )

C  Calculate overlap distributions.  
C  Form products of Gaussian functions for particle 1:
      call Gauss_prod(S1,CB1,Smat1,Smat3,KS1S3,qSCB1,W1,Wmat13,WQmat13)

C  Should now have the following 2-particle integral:
C  (        
C  | Exp(-W1(x1-Wmat13)**2) Exp(-P2(x2-Pmat2)**2  
C  | 1/|x1-x2| 
C  )

CCWS CHECK(
      halfW1=0.5d+00*W1
      halfP2=0.5d+00*P2
      call gfvee(0,0,0,halfW1,Wmat13,
     1           0,0,0,halfP2,Pmat2,
     2           0,0,0,halfW1,Wmat13,
     3           0,0,0,halfP2,Pmat2,
     4           xans)

      CHKANS=KAB1(1)*KAB2(1)*KAB3(1)*KAB4(1)
     x      *KAB1(2)*KAB2(2)*KAB3(2)*KAB4(2)
     x      *KAB1(3)*KAB2(3)*KAB3(3)*KAB4(3)
     x      *KP1P4(1)*KP3P4(1)*KS1S3(1)
     x      *KP1P4(2)*KP3P4(2)*KS1S3(2)
     x      *KP1P4(3)*KP3P4(3)*KS1S3(3)
     x      *C0*CB0*xans


      CHK   =(Wmat13(1)-Pmat2(1))**2 
     x      +(Wmat13(2)-Pmat2(2))**2 
     x      +(Wmat13(3)-Pmat2(3))**2
      write(*,*)'CHK R2WP  =',CHK
      write(*,*)'CHKANS INT=',CHKANS
CCWS CHECK)
C  Evaluate integral over spherical Gaussians:

      Cof1=P1/(S1+CB1)
      Cof2=-1.0d+00
      Cof3=(CB1/(S1+CB1))*(P3/(P3+C3))
      Cof4=(C1/(S1+CB1))+(CB1/(S1+CB1))*(C3/(P3+C3))

      RWP(1)=Cof1*Pmat1(1)+Cof2*Pmat2(1)+Cof3*Pmat3(1)+Cof4*Pmat4(1)
      RWP(2)=Cof1*Pmat1(2)+Cof2*Pmat2(2)+Cof3*Pmat3(2)+Cof4*Pmat4(2)
      RWP(3)=Cof1*Pmat1(3)+Cof2*Pmat2(3)+Cof3*Pmat3(3)+Cof4*Pmat4(3)

      R2WP=RWP(1)**2+RWP(2)**2+RWP(3)**2
      write(*,*)'    R2WP  =',R2WP

      C4=(2.0d+00*PI*PI*sqrt(PI)) / (W1*P2*sqrt(W1+P2))
      alpha=(W1*P2)/(W1+P2)
      expT=alpha*R2WP

CCCC
c     F0=half*sqrt(pi/expT)*erf(sqrt(expT))
      call RTUV(0,0,0,0,
     x          expT,alpha,RWP(1),RWP(2),RWP(3),F0)

CCCC

      xgVeeg=KAB1(1)*KAB2(1)*KAB3(1)*KAB4(1)
     x      *KAB1(2)*KAB2(2)*KAB3(2)*KAB4(2)
     x      *KAB1(3)*KAB2(3)*KAB3(3)*KAB4(3)
     x      *KP1P4(1)*KP3P4(1)*KS1S3(1)
     x      *KP1P4(2)*KP3P4(2)*KS1S3(2)
     x      *KP1P4(3)*KP3P4(3)*KS1S3(3)
     x      *C0*CB0*C4*F0

c     write(*,*)
c     write(*,*)'---SO1---'
c     write(*,*)'KAB1=',KAB1(1)*KAB1(2)*KAB1(3)
c     write(*,*)'KAB2=',KAB2(1)*KAB2(2)*KAB2(3)
c     write(*,*)'KAB3=',KAB3(1)*KAB3(2)*KAB3(3)
c     write(*,*)'KP1P4=',KP1P4(1)*KP1P4(2)*KP1P4(3)
c     write(*,*)'KP3P4=',KP3P4(1)*KP3P4(2)*KP3P4(3)
c     write(*,*)'KS1S3=',KS1S3(1)*KS1S3(2)*KS1S3(3)
c     write(*,*)'C0=',C0
c     write(*,*)'CB0=',CB0
c     write(*,*)'C4=',C4
c     write(*,*)'alpha=',alpha
c     write(*,*)'R2WP=',R2WP
c     write(*,*)'expT=',expT
c     write(*,*)'F0=',F0
c     write(*,*)'---SO1---'
c     write(*,*)


      return
      end

C=======================================================================
      subroutine cws_G4_xgVeeg_SONLY2(I1,J1,K1,A1,Amat1,
     *                                I2,J2,K2,A2,Amat2,
     *                                I3,J3,K3,A3,Amat3,
     *                                I4,J4,K4,A4,Amat4,
     *                                L1,M1,N1,B1,Bmat1,
     *                                L2,M2,N2,B2,Bmat2,
     *                                L3,M3,N3,B3,Bmat3,
     *                                L4,M4,N4,B4,Bmat4,
     *                                gam14,gam34,
     *                                xgVeeg)


C=======================================================================
      implicit none

C Input Variables
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

      double precision gam14
      double precision gam34

C Variables Returned
      double precision xgVeeg

C Local Variables
      double precision KAB1(3)
      double precision KAB2(3)
      double precision KAB3(3)
      double precision KAB4(3)

      double precision P1
      double precision P2
      double precision P3
      double precision P4

      double precision Pmat1(3)
      double precision Pmat2(3)
      double precision Pmat3(3)
      double precision Pmat4(3)

      double precision qAB1
      double precision qAB2
      double precision qAB3
      double precision qAB4

      double precision Qmat1(3)
      double precision Qmat2(3)
      double precision Qmat3(3)
      double precision Qmat4(3)

      double precision c0
      double precision c1
      double precision c3
      double precision c13

      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23
      double precision xans 

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
      call Gauss_prod(A3,B3,Amat3,Bmat3,KAB3,qAB3,P3,Pmat3,Qmat3)
      call Gauss_prod(A4,B4,Amat4,Bmat4,KAB4,qAB4,P4,Pmat4,Qmat4)
C  Apply Boys' Lemma 2 to particle 4:
C     (
C IN: |dx4 Exp(-P4(x4-Pmat4)**2) Exp(-g14(x1-x4)**2) Exp(-g34(x3-x4)**2)
C     )
C OUT: C0*Exp(-C1(x1-Pmat4)**2 Exp(-C3(x3-Pmat4)**2 Exp(-C13(x1-x3)**2
c     C1=(P4*gam14)/xsum
c     C3=(P4*gam34)/xsum
c     C13=(gam14*gam34)/xsum

      call SFBoys2(P4,gam14,gam34,Pmat4,
     x             C0,C1,C3,C13)

C  Should now have the following integral:
C  (
C  | Exp(-P1(x1-Pmat1)**2) Exp(-P2(x2-Pmat2)**2) Exp(-P3(x3-Pmat3)**2)  
C  | Exp(-C1(x1-Pmat4)**2)                       Exp(-C3(x3-Pmat4)**2)  
C  | Exp(-C13(x1-x3)**2) 1/|x1-x2| 
C  )

      gamA12=0.0d+00 
      gamA13=C13 
      gamA23=0.0d+00 
      gamB12=0.0d+00 
      gamB13=0.0d+00 
      gamB23=0.0d+00 

c     call cws_gam2_xgVeeg(0,0,0,P1,Pmat1,
      call G3_MD_xgVeeg(0,0,0,P1,Pmat1,
     *                  0,0,0,A2,Amat2,
     *                  0,0,0,P3,Pmat3,
     *                  0,0,0,C1,Pmat4,
     *                  0,0,0,B2,Bmat2,
     *                  0,0,0,C3,Pmat4,
     *                  gamA12,gamA13,gamA23,
     *                  gamB12,gamB13,gamB23,
     *                  xans)

      xgVeeg=KAB1(1)*KAB3(1)*KAB4(1)
     x      *KAB1(2)*KAB3(2)*KAB4(2)
     x      *KAB1(3)*KAB3(3)*KAB4(3)
     x      *C0*xans


      return
      end

C=======================================================================
      subroutine cws_G4_xgVepg_SONLY1(I1,J1,K1,A1,Amat1,
     *                                I2,J2,K2,A2,Amat2,
     *                                I3,J3,K3,A3,Amat3,
     *                                I4,J4,K4,A4,Amat4,
     *                                L1,M1,N1,B1,Bmat1,
     *                                L2,M2,N2,B2,Bmat2,
     *                                L3,M3,N3,B3,Bmat3,
     *                                L4,M4,N4,B4,Bmat4,
     *                                gam14,gam24,
     *                                xgVepg)


C=======================================================================
      implicit none

C Input Variables
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

      double precision gam14
      double precision gam24

C Variables Returned
      double precision xgVepg

C Local Variables
      double precision KAB1(3)
      double precision KAB2(3)
      double precision KAB3(3)
      double precision KAB4(3)

      double precision P1
      double precision P2
      double precision P3
      double precision P4

      double precision Pmat1(3)
      double precision Pmat2(3)
      double precision Pmat3(3)
      double precision Pmat4(3)

      double precision qAB1
      double precision qAB2
      double precision qAB3
      double precision qAB4

      double precision Qmat1(3)
      double precision Qmat2(3)
      double precision Qmat3(3)
      double precision Qmat4(3)

      double precision c0
      double precision c1
      double precision c3
      double precision c13

      double precision KP1P4(3)
      double precision qPC1
      double precision qPC3
      double precision S4
      double precision S3
      double precision Smat4(3)
      double precision Smat3(3)
      double precision SQmat4(3)
      double precision SQmat3(3)

      double precision CB0
      double precision CB1
      double precision PI
      parameter(PI=3.14159265358979d+00)
      double precision half
      parameter(half=0.5d+00)
      double precision two,three
      parameter(two=2.0d+00,three=3.0d+00)

      double precision KS4P2(3)
      double precision qSCB1
      double precision T4
      double precision Tmat4(3)
      double precision TQmat4(3)

      double precision C4
      double precision F0
      double precision alpha
      double precision expT
      double precision R2PT
      double precision RPT(3)
c     double precision Cof1
c     double precision Cof2
c     double precision Cof3
c     double precision Cof4



C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
      call Gauss_prod(A3,B3,Amat3,Bmat3,KAB3,qAB3,P3,Pmat3,Qmat3)
      call Gauss_prod(A4,B4,Amat4,Bmat4,KAB4,qAB4,P4,Pmat4,Qmat4)
C  Apply Boys' Lemma 2 to particle 1:
C     (
C IN: |dx1 Exp(-P1(x1-Pmat1)**2) Exp(-g14(x1-x4)**2) 
C     )
C OUT: C0*Exp(-C1(x4-Pmat1)**2
C      C0=(PI/(P1+gam14))**(3/2) 
c      C1=(P1*gam14)/(P1+gam14)

       C0=(PI/(P1+gam14))**(three/two)
       C1=(P1*gam14)/(P1+gam14)
      
C  Should now have the following integral:
C  (
C  | Exp(-P2(x2-Pmat2)**2) Exp(-P3(x3-Pmat3)**2)  
C  | Exp(-C1(x4-Pmat1)**2) Exp(-P4(x4-Pmat4)**2)  
C  | Exp(-gam24(x2-x4)**2) 1/|x3-x4| 
C  )

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(C1,P4,Pmat1,Pmat4,KP1P4,qPC1,S4,Smat4,SQmat4)

C  Should now have the following integral:
C  (        
C  | Exp(-P2(x2-Pmat2)**2) Exp(-P3(x3-Pmat3)**2 Exp(-S4(x4-Smat4)**2)  
C  | Exp(-gam14(x2-x4)**2) 1/|x3-x4| 
C  )

C  Apply Boys' Lemma 2 to particle 2:
C     (
C IN: |dx2 Exp(-P2(x2-Pmat2)**2) Exp(-gam24(x2-x4)**2) 
C     )
C OUT: CB0*Exp(-CB1(x1-Smat3)**2 
c      CB0=( PI / (P2 + gam24) )**(3/2)
c      CB1=( P2 gam24 ) / ( P2 + gam24 )

      CB0=(PI/(P2+gam24))*sqrt(PI/(P2+gam24))
      CB1=(P2*gam24)/(P2+gam24)

C  Should now have the following 2-particle integral:
C  (        
C  | Exp(-P3(x3-Pmat3)**2) Exp(-S4(x4-Smat4)**2 Exp(-CB1(x4-Pmat2)**2)  
C  | 1/|x3-x4| 
C  )

C  Form products of Gaussian functions for particle 4:
      call Gauss_prod(S4,CB1,Smat4,Pmat2,KS4P2,qSCB1,T4,Tmat4,TQmat4)

C  Should now have the following 2-particle integral:
C  (        
C  | Exp(-P3(x3-Pmat3)**2) Exp(-T4(x4-Tmat4)**2  
C  | 1/|x3-x4| 
C  )

C  Evaluate integral over spherical Gaussians:

c     Cof1=P1/(S1+CB1)
c     Cof2=-1.0d+00
c     Cof3=(CB1/(S1+CB1))*(P3/(P3+C3))
c     Cof4=(C1/(S1+CB1))+(CB1/(S1+CB1))*(C3/(P3+C3))

c     RWP(1)=Cof1*Pmat1(1)+Cof2*Pmat2(1)+Cof3*Pmat3(1)+Cof4*Pmat4(1)
c     RWP(2)=Cof1*Pmat1(2)+Cof2*Pmat2(2)+Cof3*Pmat3(2)+Cof4*Pmat4(2)
c     RWP(3)=Cof1*Pmat1(3)+Cof2*Pmat2(3)+Cof3*Pmat3(3)+Cof4*Pmat4(3)

c     R2WP=RWP(1)**2+RWP(2)**2+RWP(3)**2

      RPT(1)=Pmat3(1)-Tmat4(1)
      RPT(2)=Pmat3(2)-Tmat4(2)
      RPT(3)=Pmat3(3)-Tmat4(3)

      R2PT=RPT(1)**2 + RPT(2)**2 + RPT(3)**2

      C4=(2.0d+00*PI*PI*sqrt(PI)) / (P3*T4*sqrt(P3+T4))
      alpha=(P3*T4)/(P3+T4)
      expT=alpha*R2PT

CCCC
c     F0=half*sqrt(pi/expT)*erf(sqrt(expT))
      call RTUV(0,0,0,0,
     x          expT,alpha,RPT(1),RPT(2),RPT(3),F0)

c      write(*,*)'VEP SO1 dfBoys=',F0
CCCC

      xgVepg=KAB1(1)*KAB2(1)*KAB3(1)*KAB4(1)
     x      *KAB1(2)*KAB2(2)*KAB3(2)*KAB4(2)
     x      *KAB1(3)*KAB2(3)*KAB3(3)*KAB4(3)
     x      *KP1P4(1)*KS4P2(1)
     x      *KP1P4(2)*KS4P2(2)
     x      *KP1P4(3)*KS4P2(3)
     x      *C0*CB0*C4*F0

c     write(*,*)
c     write(*,*)'---G4Vep SO1---'
c     write(*,*)'KAB1=',KAB1(1)*KAB1(2)*KAB1(3)
c     write(*,*)'KAB2=',KAB2(1)*KAB2(2)*KAB2(3)
c     write(*,*)'KAB3=',KAB3(1)*KAB3(2)*KAB3(3)
c     write(*,*)'KP1P4=',KP1P4(1)*KP1P4(2)*KP1P4(3)
c     write(*,*)'KS4P2=',KS4P2(1)*KS4P2(2)*KS4P2(3)
c     write(*,*)'C0=',C0
c     write(*,*)'CB0=',CB0
c     write(*,*)'C4=',C4
c     write(*,*)'alpha=',alpha
c     write(*,*)'R2PT=',R2PT
c     write(*,*)'expT=',expT
c     write(*,*)'F0=',F0
c     write(*,*)'---G4Vep SO1---'
c     write(*,*)


      return
      end

C=======================================================================
      subroutine cws_G4_xgVepg_SONLY2(I1,J1,K1,A1,Amat1,
     *                                I2,J2,K2,A2,Amat2,
     *                                I3,J3,K3,A3,Amat3,
     *                                I4,J4,K4,A4,Amat4,
     *                                L1,M1,N1,B1,Bmat1,
     *                                L2,M2,N2,B2,Bmat2,
     *                                L3,M3,N3,B3,Bmat3,
     *                                L4,M4,N4,B4,Bmat4,
     *                                gam14,gam24,
     *                                xgVepg)


C=======================================================================
      implicit none

C Input Variables
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

      double precision gam14
      double precision gam24

C Variables Returned
      double precision xgVepg

C Local Variables
      double precision KAB1(3)
      double precision KAB2(3)
      double precision KAB3(3)
      double precision KAB4(3)

      double precision P1
      double precision P2
      double precision P3
      double precision P4

      double precision Pmat1(3)
      double precision Pmat2(3)
      double precision Pmat3(3)
      double precision Pmat4(3)

      double precision qAB1
      double precision qAB2
      double precision qAB3
      double precision qAB4

      double precision Qmat1(3)
      double precision Qmat2(3)
      double precision Qmat3(3)
      double precision Qmat4(3)

      double precision c0
      double precision c1
      double precision c3
      double precision c13

      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23
      double precision xans 

      double precision PI
      parameter(PI=3.14159265358979d+00)
      double precision two,three
      parameter(two=2.0d+00,three=3.0d+00)

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
      call Gauss_prod(A3,B3,Amat3,Bmat3,KAB3,qAB3,P3,Pmat3,Qmat3)
      call Gauss_prod(A4,B4,Amat4,Bmat4,KAB4,qAB4,P4,Pmat4,Qmat4)

C  Apply Boys' Lemma 2 to particle 1:
C     (
C IN: |dx1 Exp(-P1(x1-Pmat1)**2) Exp(-g14(x1-x4)**2) 
C     )
C OUT: C0*Exp(-C1(x4-Pmat1)**2
C      C0=(PI/(P1+gam14))**(3/2) 
c      C1=(P1*gam14)/(P1+gam14)

       C0=(PI/(P1+gam14))**(three/two)
       C1=(P1*gam14)/(P1+gam14)
      
C  Should now have the following integral:
C  (
C  | Exp(-P2(x2-Pmat2)**2) Exp(-P3(x3-Pmat3)**2)  
C  | Exp(-C1(x4-Pmat1)**2) Exp(-P4(x4-Pmat4)**2)  
C  | Exp(-gam24(x2-x4)**2) 1/|x3-x4| 
C  )


C>>OLD
C  Apply Boys' Lemma 2 to particle 4:
C     (
C IN: |dx4 Exp(-P4(x4-Pmat4)**2) Exp(-g14(x1-x4)**2) Exp(-g34(x3-x4)**2)
C     )
C OUT: C0*Exp(-C1(x1-Pmat4)**2 Exp(-C3(x3-Pmat4)**2 Exp(-C13(x1-x3)**2
c     C1=(P4*gam14)/xsum
c     C3=(P4*gam34)/xsum
c     C13=(gam14*gam34)/xsum

c     call SFBoys2(P4,gam14,gam34,Pmat4,
c    x             C0,C1,C3,C13)

C  Should now have the following integral:
C  (
C  | Exp(-P1(x1-Pmat1)**2) Exp(-P2(x2-Pmat2)**2) Exp(-P3(x3-Pmat3)**2)  
C  | Exp(-C1(x1-Pmat4)**2)                       Exp(-C3(x3-Pmat4)**2)  
C  | Exp(-C13(x1-x3)**2) 1/|x1-x2| 
C  )
C>>OLD

      gamA12=0.0d+00 
      gamA13=0.0d+00 
      gamA23=gam24 
      gamB12=0.0d+00 
      gamB13=0.0d+00 
      gamB23=0.0d+00 

c     call cws_gam2_xgVeeg(0,0,0,A3,Amat3,
      call G3_MD_xgVeeg(0,0,0,A3,Amat3,
     *                  0,0,0,C1,Pmat1,
     *                  0,0,0,A2,Amat2,
     *                  0,0,0,B3,Bmat3,
     *                  0,0,0,P4,Pmat4,
     *                  0,0,0,B2,Bmat2,
     *                  gamA12,gamA13,gamA23,
     *                  gamB12,gamB13,gamB23,
     *                  xans)

      xgVepg=KAB1(1)*KAB4(1)
     x      *KAB1(2)*KAB4(2)
     x      *KAB1(3)*KAB4(3)
     x      *C0*xans


      return
      end

C=======================================================================
      subroutine G4_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        I4,J4,K4,A4,Amat4,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gam14,gam34,
     *                        xgVeeg)


C=======================================================================
      implicit none

C Input Variables
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

      double precision gam14
      double precision gam34

C Variables Returned
      double precision xgVeeg

C Local Variables
      integer t1,u1,v1
      integer t2,u2,v2
      integer t3,u3,v3
      integer t4,u4,v4

      double precision xi1,e1,f1,g1,h1 
      double precision xi2,e2,f2,g2,h2 
      double precision xi3,e3,f3,g3,h3 

      double precision KAB1(3)
      double precision KAB2(3)
      double precision KAB3(3)
      double precision KAB4(3)

      double precision P1
      double precision P2
      double precision P3
      double precision P4

      double precision Pmat1(3)
      double precision Pmat2(3)
      double precision Pmat3(3)
      double precision Pmat4(3)

      double precision qAB1
      double precision qAB2
      double precision qAB3
      double precision qAB4

      double precision Qmat1(3)
      double precision Qmat2(3)
      double precision Qmat3(3)
      double precision Qmat4(3)

      double precision Et1(0:I1+L1)
      double precision Eu1(0:J1+M1)
      double precision Ev1(0:K1+N1)

      double precision Et2(0:I2+L2)
      double precision Eu2(0:J2+M2)
      double precision Ev2(0:K2+N2)

      double precision Et3(0:I3+L3)
      double precision Eu3(0:J3+M3)
      double precision Ev3(0:K3+N3)

      double precision Et4(0:I4+L4)
      double precision Eu4(0:J4+M4)
      double precision Ev4(0:K4+N4)

      double precision c0
      double precision c1
      double precision c3
      double precision c13

      double precision KP1P4(3)
      double precision KP3P4(3)
      double precision qPC1
      double precision qPC3
      double precision S1
      double precision S3
      double precision Smat1(3)
      double precision Smat3(3)
      double precision SQmat1(3)
      double precision SQmat3(3)

      double precision CB0
      double precision CB1
      double precision PI
      parameter(PI=3.14159265358979d+00)

      double precision KS1S3(3)
      double precision qSCB1
      double precision W1
      double precision Wmat13(3)
      double precision WQmat13(3)

      double precision C4
      double precision alpha
      double precision expT
      double precision R2WP 
      double precision RWP(3) 
      double precision Cof1
      double precision Cof2
      double precision Cof3
      double precision Cof4

      double precision Dgf
      double precision ans



C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
      call Gauss_prod(A3,B3,Amat3,Bmat3,KAB3,qAB3,P3,Pmat3,Qmat3)
      call Gauss_prod(A4,B4,Amat4,Bmat4,KAB4,qAB4,P4,Pmat4,Qmat4)

C  Calculate Hermite expansion coefficients:
      call ghec(I1,L1,KAB1(1),A1,B1,P1,qAB1,Qmat1(1),Et1)
      call ghec(J1,M1,KAB1(2),A1,B1,P1,qAB1,Qmat1(2),Eu1)
      call ghec(K1,N1,KAB1(3),A1,B1,P1,qAB1,Qmat1(3),Ev1)

      call ghec(I2,L2,KAB2(1),A2,B2,P2,qAB2,Qmat2(1),Et2)
      call ghec(J2,M2,KAB2(2),A2,B2,P2,qAB2,Qmat2(2),Eu2)
      call ghec(K2,N2,KAB2(3),A2,B2,P2,qAB2,Qmat2(3),Ev2)

      call ghec(I3,L3,KAB3(1),A3,B3,P3,qAB3,Qmat3(1),Et3)
      call ghec(J3,M3,KAB3(2),A3,B3,P3,qAB3,Qmat3(2),Eu3)
      call ghec(K3,N3,KAB3(3),A3,B3,P3,qAB3,Qmat3(3),Ev3)

      call ghec(I4,L4,KAB4(1),A4,B4,P4,qAB4,Qmat4(1),Et4)
      call ghec(J4,M4,KAB4(2),A4,B4,P4,qAB4,Qmat4(2),Eu4)
      call ghec(K4,N4,KAB4(3),A4,B4,P4,qAB4,Qmat4(3),Ev4)

C Evaluate integral over hermite gaussians

      xgVeeg=0.0d+00

      do t1=0,I1+L1
      do u1=0,J1+M1
      do v1=0,K1+N1

         do t2=0,I2+L2
         do u2=0,J2+M2
         do v2=0,K2+N2

            do t3=0,I3+L3
            do u3=0,J3+M3
            do v3=0,K3+N3

               do t4=0,I4+L4
               do u4=0,J4+M4
               do v4=0,K4+N4

C  Apply Boys' Lemma 2 to particle 4:
C     (
C IN: |dx4 Exp(-P4(x4-Pmat4)**2) Exp(-g14(x1-x4)**2) Exp(-g34(x3-x4)**2)
C     )
C OUT: C0*Exp(-C1(x1-Pmat4)**2 Exp(-C3(x3-Pmat4)**2 Exp(-C13(x1-x3)**2
c     C1=(P4*gam14)/xsum
c     C3=(P4*gam34)/xsum
c     C13=(gam14*gam34)/xsum

                  call SFBoys2(P4,gam14,gam34,Pmat4,
     x                         C0,C1,C3,C13)

C  Should now have the following integral:
C  (
C  | Exp(-P1(x1-Pmat1)**2) Exp(-P2(x2-Pmat2)**2) Exp(-P3(x3-Pmat3)**2)  
C  | Exp(-C1(x1-Pmat4)**2)                       Exp(-C3(x3-Pmat4)**2)  
C  | Exp(-C13(x1-x3)**2) 1/|x1-x2| 
C  )

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(P1,C1,Pmat1,Pmat4,KP1P4,qPC1,S1,Smat1,SQmat1)
      call Gauss_prod(P3,C3,Pmat3,Pmat4,KP3P4,qPC3,S3,Smat3,SQmat3)

C  Should now have the following integral:
C  (        
C  | Exp(-S1(x1-Smat1)**2) Exp(-P2(x2-Pmat2)**2 Exp(-S3(x3-Smat3)**2)  
C  | Exp(-c13(x1-x3)**2) 1/|x1-x2| 
C  )

C  Apply Boys' Lemma 2 to particle 3:
C     (
C IN: |dx3 Exp(-S3(x3-Smat3)**2) Exp(-c13(x1-x3)**2) 
C     )
C OUT: CB0*Exp(-CB1(x1-Smat3)**2 
c     CB0=( PI / (s3 + c13) )**(3/2)
c     CB1=( S3 C13 ) / ( S3 + C13 )

                 CB0=(PI/(S3+C13))*sqrt(PI/(S3+C13))
                 CB1=(S3*C13)/(S3+C13)
      
C  Should now have the following 2-particle integral:
C  (        
C  | Exp(-S1(x1-Smat1)**2) Exp(-P2(x2-Pmat2)**2 Exp(-CB1(x1-Smat3)**2)  
C  | 1/|x1-x2| 
C  )

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(S1,CB1,Smat1,Smat3,KS1S3,qSCB1,W1,Wmat13,WQmat13)

C  Should now have the following 2-particle integral:
C  (        
C  | Exp(-W1(x1-Wmat13)**2) Exp(-P2(x2-Pmat2)**2  
C  | 1/|x1-x2| 
C  )

C  Evaluate integral over spherical Gaussians:

      Cof1=P1/(S1+CB1)
      Cof2=-1.0d+00
      Cof3=(CB1/(S1+CB1))*(P3/(P3+C3))
      Cof4=(C1/(S1+CB1))+(CB1/(S1+CB1))*(C3/(P3+C3))

      RWP(1)=Cof1*Pmat1(1)+Cof2*Pmat2(1)+Cof3*Pmat3(1)+Cof4*Pmat4(1)
      RWP(2)=Cof1*Pmat1(2)+Cof2*Pmat2(2)+Cof3*Pmat3(2)+Cof4*Pmat4(2)
      RWP(3)=Cof1*Pmat1(3)+Cof2*Pmat2(3)+Cof3*Pmat3(3)+Cof4*Pmat4(3)

      R2WP=RWP(1)**2+RWP(2)**2+RWP(3)**2

      C4=(2.0d+00*PI*PI*sqrt(PI)) / (W1*P2*sqrt(W1+P2))
      alpha=(W1*P2)/(W1+P2)
      expT=alpha*R2WP

C  Organize terms for differentiation of f_KKK:
C Prepare the following KKK Gaussian functions:
C KP1P4 <=> GK1 = exp(-xi1*( e1*x1 + f1*x2 +g1*x3 + h1*x4 )
C KP3P4 <=> GK2 = exp(-xi2*( e2*x1 + f2*x2 +g2*x3 + h2*x4 )
C KS1S3 <=> GK3 = exp(-xi3*( e3*x1 + f3*x2 +g3*x3 + h3*x4 )

      call G4_xgVeeg_prep_fKKK(P1,P3,C1,C3,S1,CB1,
     x                         xi1,e1,f1,g1,h1,
     x                         xi2,e2,f2,g2,h2,
     x                         xi3,e3,f3,g3,h3)

C  Differentiation of integral over spherical Gaussians:
      call G4_xgVee_dhermite(t1,t2,t3,t4,
     x                       u1,u2,u3,u4,
     x                       v1,v2,v3,v4,
     x                       xi1,xi2,xi3,
     x                       e1,f1,g1,h1,
     x                       e2,f2,g2,h2,
     x                       e3,f3,g3,h3,
     x                       Cof1,Cof2,Cof3,Cof4,
     x                       expT,alpha,
     x                       Pmat1,Pmat2,Pmat3,Pmat4,Dgf)


c     ans=C0*C3*CB0*C4*
      ans=C0*CB0*C4*
     x    Et1(t1)*Et2(t2)*Et3(t3)*Et4(t4) * 
     x    Eu1(u1)*Eu2(u2)*Eu3(u3)*Eu4(u4) *
     x    Ev1(v1)*Ev2(v2)*Ev3(v3)*Ev4(v4) *
     x    Dgf

      xgVeeg=xgVeeg+ans

c          write(*,*)'C0=',C0
c          write(*,*)'CB0=',CB0
c          write(*,*)'C3=',C3
c          write(*,*)'C4=',C4
c          write(*,*)'Et1 t1=',Et1(t1)
c          write(*,*)'Et2 t2=',Et2(t2)
c          write(*,*)'Et3 t3=',Et3(t3)
c          write(*,*)'Et4 t4=',Et4(t4)

c          write(*,*)'Eu1 u1=',Eu1(u1)
c          write(*,*)'Eu2 u2=',Eu2(u2)
c          write(*,*)'Eu3 u3=',Eu3(u3)
c          write(*,*)'Eu4 u4=',Eu4(u4)

c          write(*,*)'Ev1 v1=',Ev1(v1)
c          write(*,*)'Ev2 v2=',Ev2(v2)
c          write(*,*)'Ev3 v3=',Ev3(v3)
c          write(*,*)'Ev4 v4=',Ev4(v4)

c          write(*,*)'Dgf=',Dgf
c          write(*,*)


               end do
               end do
               end do

            end do
            end do
            end do

         end do
         end do
         end do

      end do
      end do
      end do


      return
      end

C=======================================================================
      subroutine G4_xgVeeg_prep_fKKK(P1,P3,C1,C3,S1,CB1,
     x                               xi1,e1,f1,g1,h1,
     x                               xi2,e2,f2,g2,h2,
     x                               xi3,e3,f3,g3,h3)

C=======================================================================
      implicit none

C Input Variables
      double precision P1
      double precision C1
      double precision P3
      double precision C3
      double precision S1
      double precision CB1

C Variables Returned
      double precision xi1,e1,f1,g1,h1
      double precision xi2,e2,f2,g2,h2
      double precision xi3,e3,f3,g3,h3

C Local Variables

      xi1=P1*C1/(P1+C1)
      xi2=P3*C3/(P3+C3)
      xi3=S1*CB1/(S1+CB1)

      e1=1.0d+00
      f1=0.0d+00
      g1=0.0d+00
      h1=-1.0d+00

      e2=0.0d+00
      f2=0.0d+00
      g2=1.0d+00
      h2=-1.0d+00

      e3=P1/(P1+C1)
      f3=0.0d+00
      g3=-P3/(P3+C3)
      h3=C1/(P1+C1) - C3/(P3+C3)


      return
      end

C=======================================================================
      subroutine G4_xgVee_dhermite(t1,t2,t3,t4,
     x                             u1,u2,u3,u4,
     x                             v1,v2,v3,v4,
     x                             xi1,xi2,xi3,
     x                             e1,f1,g1,h1,
     x                             e2,f2,g2,h2,
     x                             e3,f3,g3,h3,
     x                             Cof1,Cof2,Cof3,Cof4,
     x                             expT,alpha,
     x                             Pmat1,Pmat2,Pmat3,Pmat4,Dgf)

C=======================================================================
      implicit none

C Input variables
      integer t1,t2,t3,t4
      integer u1,u2,u3,u4
      integer v1,v2,v3,v4

      double precision xi1,xi2,xi3
      double precision e1,f1,g1,h1
      double precision e2,f2,g2,h2
      double precision e3,f3,g3,h3
      double precision expT
      double precision alpha
      double precision Pmat1(3),Pmat2(3),Pmat3(3),Pmat4(3)

C Variables returned
      double precision Dgf

C Local Variables
      integer xb1,xb2,xb3,xb4
      integer yb1,yb2,yb3,yb4
      integer zb1,zb2,zb3,zb4

      double precision bin_tx1
      double precision bin_tx2
      double precision bin_tx3
      double precision bin_tx4

      double precision bin_uy1
      double precision bin_uy2
      double precision bin_uy3
      double precision bin_uy4

      double precision bin_vz1
      double precision bin_vz2
      double precision bin_vz3
      double precision bin_vz4

      double precision dfKKK_x
      double precision dfKKK_y
      double precision dfKKK_z

c     double precision v4
      double precision dfBoys 
      double precision ans 

      double precision Cof1
      double precision Cof2
      double precision Cof3
      double precision Cof4


      Dgf=0.0d+00

      do xb1=0,t1
         call cbinom(t1,xb1,bin_tx1)
      do xb2=0,t2
         call cbinom(t2,xb2,bin_tx2)
      do xb3=0,t3
         call cbinom(t3,xb3,bin_tx3)
      do xb4=0,t4
         call cbinom(t4,xb4,bin_tx4)

         call G4xggs_DGK123(t1-xb1,t2-xb2,t3-xb3,t4-xb4,
     x                      xi1,xi2,xi3,
     x                      e1,f1,g1,h1,
     x                      e2,f2,g2,h2,
     x                      e3,f3,g3,h3,
     x                      Pmat1(1),Pmat2(1),Pmat3(1),Pmat4(1),dfKKK_x)


         do yb1=0,u1
            call cbinom(u1,yb1,bin_uy1)
         do yb2=0,u2
            call cbinom(u2,yb2,bin_uy2)
         do yb3=0,u3
            call cbinom(u3,yb3,bin_uy3)
         do yb4=0,u4
            call cbinom(u4,yb4,bin_uy4)

            call G4xggs_DGK123(u1-yb1,u2-yb2,u3-yb3,u4-yb4,
     x                         xi1,xi2,xi3,
     x                         e1,f1,g1,h1,
     x                         e2,f2,g2,h2,
     x                         e3,f3,g3,h3,
     x                      Pmat1(2),Pmat2(2),Pmat3(2),Pmat4(2),dfKKK_y)



            do zb1=0,v1
               call cbinom(v1,zb1,bin_vz1)
            do zb2=0,v2
               call cbinom(v2,zb2,bin_vz2)
            do zb3=0,v3
               call cbinom(v3,zb3,bin_vz3)
            do zb4=0,v4
               call cbinom(v4,zb4,bin_vz4)

               call G4xggs_DGK123(v1-zb1,v2-zb2,v3-zb3,v4-zb4,
     x                            xi1,xi2,xi3,
     x                            e1,f1,g1,h1,
     x                            e2,f2,g2,h2,
     x                            e3,f3,g3,h3,
     x                      Pmat1(3),Pmat2(3),Pmat3(3),Pmat4(3),dfKKK_z)

c              call gam3_RTUV(xb1,yb1,zb1,
               call G4_RTUV(xb1,yb1,zb1,
     x                      xb2,yb2,zb2,
     x                      xb3,yb3,zb3,
     x                      xb4,yb4,zb4,
     x                      expT,alpha,Cof1,Cof2,Cof3,Cof4,
     x                      Pmat1,Pmat2,Pmat3,Pmat4,dfBoys)

               ans=bin_tx1*bin_tx2*bin_tx3*bin_tx4
     x            *bin_uy1*bin_uy2*bin_uy3*bin_uy4
     x            *bin_vz1*bin_vz2*bin_vz3*bin_vz4
     x            *dfKKK_x*dfKKK_y*dfKKK_z
     x            *dfBoys
               
               Dgf=Dgf+ans



            end do
            end do
            end do
            end do

         end do
         end do
         end do
         end do

      end do
      end do
      end do
      end do





      return
      end

C=======================================================================
      subroutine G4xggs_DGK123(xb1,xb2,xb3,xb4,
     x                         xi1,xi2,xi3,
     x                         e1,f1,g1,h1,
     x                         e2,f2,g2,h2,
     x                         e3,f3,g3,h3,
     x                         x1,x2,x3,x4,dGK123)

C
C Calculate derivative P_x1,x2,x3,x4 (GK1 GK2 GK3)
C 
C = (d/P1x)^xb1 (d/dP2x)^xb2 (d/dP3x)^xb3 (d/dP4x)^xb4 GK1 GK2 GK3
C
C=======================================================================
      implicit none

C Input variables
      integer xb1,xb2,xb3,xb4

      double precision xi1,xi2,xi3
      double precision e1,f1,g1,h1
      double precision e2,f2,g2,h2
      double precision e3,f3,g3,h3
      double precision x1,x2,x3,x4

C Variables returned
      double precision dGK123

C Local variables
      integer i1,i2,i3,i4
      integer j1,j2,j3,j4

      double precision bin_i1
      double precision bin_i2
      double precision bin_i3
      double precision bin_i4

      double precision bin_j1
      double precision bin_j2
      double precision bin_j3
      double precision bin_j4

      double precision dGK1,dGK2,dGK3
      double precision ans



      dGK123=0.0d+00

      do i1=0,xb1
         call cbinom(xb1,i1,bin_i1)
      do i2=0,xb2
         call cbinom(xb2,i2,bin_i2)
      do i3=0,xb3
         call cbinom(xb3,i3,bin_i3)
      do i4=0,xb4
         call cbinom(xb4,i4,bin_i4)

         call DKGx_4(xb1-i1,xb2-i2,xb3-i3,xb4-i4,
     x               xi1,e1,f1,g1,h1,x1,x2,x3,x4,dGK1)

             do j1=0,i1
                call cbinom(i1,j1,bin_j1)
             do j2=0,i2
                call cbinom(i2,j2,bin_j2)
             do j3=0,i3
                call cbinom(i3,j3,bin_j3)
             do j4=0,i4
                call cbinom(i4,j4,bin_j4)

                call DKGx_4(i1-j1,i2-j2,i3-j3,i4-j4,
     x                      xi2,e2,f2,g2,h2,x1,x2,x3,x4,dGK2)

                call DKGx_4(j1,j2,j3,j4,
     x                      xi3,e3,f3,g3,h3,x1,x2,x3,x4,dGK3)

                ans=bin_i1*bin_i2*bin_i3*bin_i4
     x             *bin_j1*bin_j2*bin_j3*bin_j4
     x             *dGK1*dGK2*dGK3

                dGK123=dGK123+ans

             end do
             end do
             end do
             end do

      end do
      end do
      end do
      end do


      return
      end

C=======================================================================
      subroutine DKGx_4(n1,n2,n3,n4,xi,e,f,g,h,x1,x2,x3,x4,ans)
C
C Calculate derivative:
C 
C   (d/dx1)^n1 (d/dx2)^n2 (d/dx3)^n3 (d/dx4)^n4
C    [exp( -xi*(e*x1 + f*x2 + g*x3 + h*x4) ) ] 
C
C=======================================================================
      implicit none
C Input variables
      integer n1,n2,n3,n4

      double precision xi
      double precision e,f,g,h
      double precision x1,x2,x3,x4

C Variables returned
      double precision ans

C Local variables
      integer ntotal

      double precision u
      double precision c1
      double precision dGx


c     u =  (a*PX1) + (b*PX2) + (c*PX3)
c     coeff = (a**i1)*(b**i2)*(c**i3)
c     call gdel_x(i1+i2+i3,u,alp,0.0d0,gd)
c     ans = coeff * gd

      ntotal=n1+n2+n3+n4
      u = e*x1 + f*x2 + g*x3 + h*x4
      c1=(e**n1)*(f**n2)*(g**n3)*(h**n4)
      call gdel_x(ntotal,u,xi,0.0d0,dGx)
      ans=c1*dGx



      return
      end


C=======================================================================
c     subroutine gam3_RTUV(xb1,yb1,zb1,
      subroutine G4_RTUV(xb1,yb1,zb1,
     x                   xb2,yb2,zb2,
     x                   xb3,yb3,zb3,
     x                   xb4,yb4,zb4,
     x                   expT,alpha,e,f,g,h,
     x                   Pmat1,Pmat2,Pmat3,Pmat4,dfBoys)

C=======================================================================
      implicit none

C Input variables
      integer xb1,yb1,zb1
      integer xb2,yb2,zb2
      integer xb3,yb3,zb3
      integer xb4,yb4,zb4

      double precision expT
      double precision alpha
      double precision e,f,g,h
      double precision Pmat1(3),Pmat2(3),Pmat3(3),Pmat4(3)

C Variables returned
      double precision dfBoys

C Local variables
      integer N,L,M
      integer nLim

      double precision c1
      double precision c2
      double precision c3
      double precision c4
      double precision ans
      double precision XWP
      double precision YWP
      double precision ZWP


      XWP=(e*Pmat1(1))+(f*Pmat2(1))+(g*Pmat3(1))+(h*Pmat4(1))
      YWP=(e*Pmat1(2))+(f*Pmat2(2))+(g*Pmat3(2))+(h*Pmat4(2))
      ZWP=(e*Pmat1(3))+(f*Pmat2(3))+(g*Pmat3(3))+(h*Pmat4(3))

      N=xb1+xb2+xb3+xb4
      L=yb1+yb2+yb3+yb4
      M=zb1+zb2+zb3+zb4

      c1=e**(xb1+yb1+zb1)
      c2=f**(xb2+yb2+zb2)
      c3=g**(xb3+yb3+zb3)
      c4=h**(xb4+yb4+zb4)

      nLim=N+L+M

CCWS-DEBUG
c     zwp=-7.165109034267912D-002
CCWS-DEBUG
      call RTUV(N,L,M,nLim,
     x          expT,alpha,XWP,YWP,ZWP,ans)

      dfBoys=c1*c2*c3*c4*ans

CCWS-DEBUG
c     write(*,*)'----G4MD----'
c     write(*,*)'N=',N
c     write(*,*)'L=',L
c     write(*,*)'M=',M
c     write(*,*)'XWP=',XWP
c     write(*,*)'YWP=',YWP
c     write(*,*)'ZWP=',ZWP
c     write(*,*)'RNLM=',ans
c     write(*,*)'c1=',c1
c     write(*,*)'c2=',c2
c     write(*,*)'c3=',c3
c     write(*,*)'c4=',c4
c     write(*,*)'alpha=',alpha
c     write(*,*)'----G4MD----'
CCWS-DEBUG

      


      return
      end


C=======================================================================
      subroutine G4_MD_xgVepg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        I4,J4,K4,A4,Amat4,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gam14,gam24,
     *                        xgVepg)

C=======================================================================
      implicit none

C Input Variables
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

      double precision gam14
      double precision gam24

C Variables Returned
      double precision xgVepg

C Local Variables
      integer t1,u1,v1
      integer t2,u2,v2
      integer t3,u3,v3
      integer t4,u4,v4

      double precision xi1,e1,f1,g1,h1 
      double precision xi2,e2,f2,g2,h2 
c     double precision xi3,e3,f3,g3,h3 

      double precision KAB1(3)
      double precision KAB2(3)
      double precision KAB3(3)
      double precision KAB4(3)

      double precision P1
      double precision P2
      double precision P3
      double precision P4

      double precision Pmat1(3)
      double precision Pmat2(3)
      double precision Pmat3(3)
      double precision Pmat4(3)

      double precision qAB1
      double precision qAB2
      double precision qAB3
      double precision qAB4

      double precision Qmat1(3)
      double precision Qmat2(3)
      double precision Qmat3(3)
      double precision Qmat4(3)

      double precision Et1(0:I1+L1)
      double precision Eu1(0:J1+M1)
      double precision Ev1(0:K1+N1)

      double precision Et2(0:I2+L2)
      double precision Eu2(0:J2+M2)
      double precision Ev2(0:K2+N2)

      double precision Et3(0:I3+L3)
      double precision Eu3(0:J3+M3)
      double precision Ev3(0:K3+N3)

      double precision Et4(0:I4+L4)
      double precision Eu4(0:J4+M4)
      double precision Ev4(0:K4+N4)

      double precision C0
      double precision C1
c     double precision c3
c     double precision c13

      double precision KP1P4(3)
      double precision qPC1
      double precision qPC3

      double precision S4
      double precision Smat4(3)
      double precision SQmat4(3)

c     double precision S1
c     double precision Smat1(3)
c     double precision SQmat1(3)
c     double precision S3
c     double precision Smat3(3)
c     double precision SQmat3(3)

      double precision CB0
      double precision CB1
      double precision PI
      parameter(PI=3.14159265358979d+00)

c     double precision KS4P2(3)
c     double precision qSCB1
c     double precision W1
c     double precision Wmat13(3)
c     double precision WQmat13(3)

      double precision KS4P2(3)
      double precision qSCB1
      double precision Tx4
      double precision Tmat4(3)
      double precision TQmat4(3)

      double precision C4
      double precision alpha
      double precision expT
      double precision R2TP 
      double precision RTP(3) 
      double precision xcof
      double precision Cof1
      double precision Cof2
      double precision Cof3
      double precision Cof4

      double precision Dgf
      double precision ans

CCWS-debug
c     double precision kp1p4_ans
c     double precision kp1p4X
c     double precision kp1p4Y
c     double precision kp1p4Z
c     double precision ks4p2_ans
c     double precision ks4p2X
c     double precision ks4p2Y
c     double precision ks4p2Z
CCWS-debug


C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
      call Gauss_prod(A3,B3,Amat3,Bmat3,KAB3,qAB3,P3,Pmat3,Qmat3)
      call Gauss_prod(A4,B4,Amat4,Bmat4,KAB4,qAB4,P4,Pmat4,Qmat4)

C  Calculate Hermite expansion coefficients:
      call ghec(I1,L1,KAB1(1),A1,B1,P1,qAB1,Qmat1(1),Et1)
      call ghec(J1,M1,KAB1(2),A1,B1,P1,qAB1,Qmat1(2),Eu1)
      call ghec(K1,N1,KAB1(3),A1,B1,P1,qAB1,Qmat1(3),Ev1)

      call ghec(I2,L2,KAB2(1),A2,B2,P2,qAB2,Qmat2(1),Et2)
      call ghec(J2,M2,KAB2(2),A2,B2,P2,qAB2,Qmat2(2),Eu2)
      call ghec(K2,N2,KAB2(3),A2,B2,P2,qAB2,Qmat2(3),Ev2)

      call ghec(I3,L3,KAB3(1),A3,B3,P3,qAB3,Qmat3(1),Et3)
      call ghec(J3,M3,KAB3(2),A3,B3,P3,qAB3,Qmat3(2),Eu3)
      call ghec(K3,N3,KAB3(3),A3,B3,P3,qAB3,Qmat3(3),Ev3)

      call ghec(I4,L4,KAB4(1),A4,B4,P4,qAB4,Qmat4(1),Et4)
      call ghec(J4,M4,KAB4(2),A4,B4,P4,qAB4,Qmat4(2),Eu4)
      call ghec(K4,N4,KAB4(3),A4,B4,P4,qAB4,Qmat4(3),Ev4)

C Evaluate integral over hermite gaussians

      xgVepg=0.0d+00

      do t1=0,I1+L1
      do u1=0,J1+M1
      do v1=0,K1+N1

         do t2=0,I2+L2
         do u2=0,J2+M2
         do v2=0,K2+N2

            do t3=0,I3+L3
            do u3=0,J3+M3
            do v3=0,K3+N3

               do t4=0,I4+L4
               do u4=0,J4+M4
               do v4=0,K4+N4

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C  Apply Boys' Lemma 2 to particle 1:
C     (
C IN: |dx1 Exp(-P1(x1-Pmat1)**2) Exp(-g14(x1-x4)**2) 
C     )
C OUT: C0*Exp(-C1(x4-Pmat1)**2
C      C0=(PI/(P1+gam14))**(3/2) 
c      C1=(P1*gam14)/(P1+gam14)

c      C0=(PI/(P1+gam14))**(three/two)
       C0=(PI/(P1+gam14))*sqrt(PI/(P1+gam14))
       C1=(P1*gam14)/(P1+gam14)
      
C  Should now have the following integral:
C  (
C  | Exp(-P2(x2-Pmat2)**2) Exp(-P3(x3-Pmat3)**2)  
C  | Exp(-C1(x4-Pmat1)**2) Exp(-P4(x4-Pmat4)**2)  
C  | Exp(-gam24(x2-x4)**2) 1/|x3-x4| 
C  )

C  Form products of Gaussian functions for particle 4:
      call Gauss_prod(C1,P4,Pmat1,Pmat4,KP1P4,qPC1,S4,Smat4,SQmat4)

C  Should now have the following integral:
C  (        
C  | Exp(-P2(x2-Pmat2)**2) Exp(-P3(x3-Pmat3)**2 Exp(-S4(x4-Smat4)**2)  
C  | Exp(-gam14(x2-x4)**2) 1/|x3-x4| 
C  )

C  Apply Boys' Lemma 2 to particle 2:
C     (
C IN: |dx2 Exp(-P2(x2-Pmat2)**2) Exp(-gam24(x2-x4)**2) 
C     )
C OUT: CB0*Exp(-CB1(x1-Smat3)**2 
c      CB0=( PI / (P2 + gam24) )**(3/2)
c      CB1=( P2 gam24 ) / ( P2 + gam24 )

      CB0=(PI/(P2+gam24))*sqrt(PI/(P2+gam24))
      CB1=(P2*gam24)/(P2+gam24)

C  Should now have the following 2-particle integral:
C  (        
C  | Exp(-P3(x3-Pmat3)**2) Exp(-S4(x4-Smat4)**2 Exp(-CB1(x4-Pmat2)**2)  
C  | 1/|x3-x4| 
C  )

C  Form products of Gaussian functions for particle 4:
      call Gauss_prod(S4,CB1,Smat4,Pmat2,KS4P2,qSCB1,Tx4,Tmat4,TQmat4)

C  Should now have the following 2-particle integral:
C  (        
C  | Exp(-P3(x3-Pmat3)**2) Exp(-Tx4(x4-Tmat4)**2  
C  | 1/|x3-x4| 
C  )

C  Evaluate integral over spherical Gaussians:

      xcof=(S4+CB1)*(C1+P4)
c     Cof1=C1/( (S4+CB1) * (C1+P4) )
      Cof1=S4*C1/xcof
      Cof2=CB1/(S4+CB1)
      Cof3=-1.0d+00
c     Cof4=P4/( (S4+CB1) * (C1+P4) )
      Cof4=S4*P4/xcof

CCWS-DEBUG
c      Cof1=-Cof1
c      Cof2=-Cof2
c      Cof3=-Cof3
c      Cof4=-Cof4
CCWS-DEBUG
CCWS-DEBUG(
c     RTP(1)=Cof1*Pmat1(1)+Cof2*Pmat2(1)+Cof4*Pmat4(1)
c     RTP(2)=Cof1*Pmat1(2)+Cof2*Pmat2(2)+Cof4*Pmat4(2)
c     RTP(3)=Cof1*Pmat1(3)+Cof2*Pmat2(3)+Cof4*Pmat4(3)
c     write(*,*)'COEF TMAT4:'
c     write(*,*)'TMAT4X=',RTP(1)
c     write(*,*)'TMAT4Y=',RTP(2)
c     write(*,*)'TMAT4Z=',RTP(3)
c     write(*,*)'ACTUAL TMAT4:'
c     write(*,*)'TMAT4X=',Tmat4(1)
c     write(*,*)'TMAT4Y=',Tmat4(2)
c     write(*,*)'TMAT4Z=',Tmat4(3)
CCWS-DEBUG)
CCWS-DEBUG(
c     RTP(1)=Cof1*Pmat1(1)+Cof2*Pmat2(1)+Cof3*Pmat3(1)+Cof4*Pmat4(1)
c     RTP(2)=Cof1*Pmat1(2)+Cof2*Pmat2(2)+Cof3*Pmat3(2)+Cof4*Pmat4(2)
c     RTP(3)=Cof1*Pmat1(3)+Cof2*Pmat2(3)+Cof3*Pmat3(3)+Cof4*Pmat4(3)
c     write(*,*)'COEF RTP:'
c     write(*,*)'RTP1=',RTP(1)
c     write(*,*)'RTP2=',RTP(2)
c     write(*,*)'RTP3=',RTP(3)
CCWS-DEBUG)

      RTP(1)=Tmat4(1)-Pmat3(1)
      RTP(2)=Tmat4(2)-Pmat3(2)
      RTP(3)=Tmat4(3)-Pmat3(3)

CCWS-DEBUG
c     write(*,*)'STRAIGHT RTP:'
c     write(*,*)'RTP1=',RTP(1)
c     write(*,*)'RTP2=',RTP(2)
c     write(*,*)'RTP3=',RTP(3)
CCWS-DEBUG

      R2TP=RTP(1)**2 + RTP(2)**2 + RTP(3)**2

      C4=(2.0d+00*PI*PI*sqrt(PI)) / (P3*Tx4*sqrt(P3+Tx4))
      alpha=(P3*Tx4)/(P3+Tx4)
      expT=alpha*R2TP

C  Organize terms for differentiation of f_KK:
C Prepare the following KKK Gaussian functions:
C KP1P4 <=> GK1 = exp(-xi1*( e1*x1 + f1*x2 +g1*x3 + h1*x4 )
C KP3P4 <=> GK2 = exp(-xi2*( e2*x1 + f2*x2 +g2*x3 + h2*x4 )

      xi1 = (P4*C1)/(P4+C1)
      e1 = 1.0d+00
      f1 = 0.0d+00
      g1 = 0.0d+00
      h1 = -1.0d+00

      xi2 = (S4*CB1)/(S4+CB1)
      e2 = C1/(C1+P4)
      f2 = -1.0d+00
      g2 = 0.0d+00
      h2 = P4/(C1+P4)
CCWS-debug
c     kp1p4X=exp(-xi1*(e1*Pmat1(1)
c    x                +f1*Pmat2(1)
c    x                +g1*Pmat3(1)
c    x                +h1*Pmat4(1))**2)
c     kp1p4Y=exp(-xi1*(e1*Pmat1(2)
c    x                +f1*Pmat2(2)
c    x                +g1*Pmat3(2)
c    x                +h1*Pmat4(2))**2)
c     kp1p4Z=exp(-xi1*(e1*Pmat1(3)
c    x                +f1*Pmat2(3)
c    x                +g1*Pmat3(3)
c    x                +h1*Pmat4(3))**2)
c     kp1p4_ans=kp1p4X*kp1p4Y*kp1p4Z
c     write(*,*)'MD:  KP1P4=',kp1p4_ans
c
c     ks4p2X=exp(-xi2*(e2*Pmat1(1)
c    x                +f2*Pmat2(1)
c    x                +g2*Pmat3(1)
c    x                +h2*Pmat4(1))**2)
c     ks4p2Y=exp(-xi2*(e2*Pmat1(2)
c    x                +f2*Pmat2(2)
c    x                +g2*Pmat3(2)
c    x                +h2*Pmat4(2))**2)
c     ks4p2Z=exp(-xi2*(e2*Pmat1(3)
c    x                +f2*Pmat2(3)
c    x                +g2*Pmat3(3)
c    x                +h2*Pmat4(3))**2)
c     ks4p2_ans=ks4p2X*ks4p2Y*ks4p2Z
c     write(*,*)'MD:  KS4P2=',ks4p2_ans
CCWS-debug

C  Differentiation of integral over spherical Gaussians:
      call G4_xgVep_dhermite(t1,t2,t3,t4,
     x                       u1,u2,u3,u4,
     x                       v1,v2,v3,v4,
     x                       xi1,xi2,
     x                       e1,f1,g1,h1,
     x                       e2,f2,g2,h2,
     x                       Cof1,Cof2,Cof3,Cof4,
     x                       expT,alpha,
     x                       Pmat1,Pmat2,Pmat3,Pmat4,Dgf)

      ans=C0*CB0*C4*
     x    Et1(t1)*Et2(t2)*Et3(t3)*Et4(t4) * 
     x    Eu1(u1)*Eu2(u2)*Eu3(u3)*Eu4(u4) *
     x    Ev1(v1)*Ev2(v2)*Ev3(v3)*Ev4(v4) *
     x    Dgf

      xgVepg=xgVepg+ans

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c          write(*,*)'C0=',C0
c          write(*,*)'CB0=',CB0
c          write(*,*)'C3=',C3
c          write(*,*)'C4=',C4
c          write(*,*)'Et1 t1=',Et1(t1)
c          write(*,*)'Et2 t2=',Et2(t2)
c          write(*,*)'Et3 t3=',Et3(t3)
c          write(*,*)'Et4 t4=',Et4(t4)

c          write(*,*)'Eu1 u1=',Eu1(u1)
c          write(*,*)'Eu2 u2=',Eu2(u2)
c          write(*,*)'Eu3 u3=',Eu3(u3)
c          write(*,*)'Eu4 u4=',Eu4(u4)

c          write(*,*)'Ev1 v1=',Ev1(v1)
c          write(*,*)'Ev2 v2=',Ev2(v2)
c          write(*,*)'Ev3 v3=',Ev3(v3)
c          write(*,*)'Ev4 v4=',Ev4(v4)

c          write(*,*)'Dgf=',Dgf
c          write(*,*)


               end do
               end do
               end do

            end do
            end do
            end do

         end do
         end do
         end do

      end do
      end do
      end do


      return
      end

C=======================================================================
      subroutine G4_xgVep_dhermite(t1,t2,t3,t4,
     x                             u1,u2,u3,u4,
     x                             v1,v2,v3,v4,
     x                             xi1,xi2,
     x                             e1,f1,g1,h1,
     x                             e2,f2,g2,h2,
     x                             Cof1,Cof2,Cof3,Cof4,
     x                             expT,alpha,
     x                             Pmat1,Pmat2,Pmat3,Pmat4,Dgf)

C=======================================================================
      implicit none

C Input variables
      integer t1,t2,t3,t4
      integer u1,u2,u3,u4
      integer v1,v2,v3,v4

      double precision xi1,xi2
      double precision e1,f1,g1,h1
      double precision e2,f2,g2,h2
      double precision expT
      double precision alpha
      double precision Pmat1(3),Pmat2(3),Pmat3(3),Pmat4(3)

      double precision Cof1
      double precision Cof2
      double precision Cof3
      double precision Cof4

C Variables returned
      double precision Dgf

C Local Variables
      integer xb1,xb2,xb3,xb4
      integer yb1,yb2,yb3,yb4
      integer zb1,zb2,zb3,zb4

      double precision bin_tx1
      double precision bin_tx2
      double precision bin_tx3
      double precision bin_tx4

      double precision bin_uy1
      double precision bin_uy2
      double precision bin_uy3
      double precision bin_uy4

      double precision bin_vz1
      double precision bin_vz2
      double precision bin_vz3
      double precision bin_vz4

      double precision dfKKK_x
      double precision dfKKK_y
      double precision dfKKK_z

c     double precision v4
      double precision dfBoys 
      double precision ans 


      Dgf=0.0d+00

      do xb1=0,t1
         call cbinom(t1,xb1,bin_tx1)
      do xb2=0,t2
         call cbinom(t2,xb2,bin_tx2)
      do xb3=0,t3
         call cbinom(t3,xb3,bin_tx3)
      do xb4=0,t4
         call cbinom(t4,xb4,bin_tx4)

         call G4xggs_DGK12(t1-xb1,t2-xb2,t3-xb3,t4-xb4,
     x                     xi1,xi2,
     x                     e1,f1,g1,h1,
     x                     e2,f2,g2,h2,
     x                     Pmat1(1),Pmat2(1),Pmat3(1),Pmat4(1),dfKKK_x)


         do yb1=0,u1
            call cbinom(u1,yb1,bin_uy1)
         do yb2=0,u2
            call cbinom(u2,yb2,bin_uy2)
         do yb3=0,u3
            call cbinom(u3,yb3,bin_uy3)
         do yb4=0,u4
            call cbinom(u4,yb4,bin_uy4)

            call G4xggs_DGK12(u1-yb1,u2-yb2,u3-yb3,u4-yb4,
     x                        xi1,xi2,
     x                        e1,f1,g1,h1,
     x                        e2,f2,g2,h2,
     x                      Pmat1(2),Pmat2(2),Pmat3(2),Pmat4(2),dfKKK_y)



            do zb1=0,v1
               call cbinom(v1,zb1,bin_vz1)
            do zb2=0,v2
               call cbinom(v2,zb2,bin_vz2)
            do zb3=0,v3
               call cbinom(v3,zb3,bin_vz3)
            do zb4=0,v4
               call cbinom(v4,zb4,bin_vz4)

               call G4xggs_DGK12(v1-zb1,v2-zb2,v3-zb3,v4-zb4,
     x                           xi1,xi2,
     x                           e1,f1,g1,h1,
     x                           e2,f2,g2,h2,
     x                      Pmat1(3),Pmat2(3),Pmat3(3),Pmat4(3),dfKKK_z)
CCWS-debug
c        dfKKK_z=0.944869744244220*0.862004514134267
CCWS-debug

c              call gam3_RTUV(xb1,yb1,zb1,
CCWS-debug
c         write(*,*)'Call G4Rtuv'
c         write(*,*)'v4=',v4
c         write(*,*)'zb1=',zb1
c         write(*,*)'zb2=',zb2
c         write(*,*)'zb3=',zb3
c         write(*,*)'zb4=',zb4
CCWS-debug
               call G4_RTUV(xb1,yb1,zb1,
     x                      xb2,yb2,zb2,
     x                      xb3,yb3,zb3,
     x                      xb4,yb4,zb4,
     x                      expT,alpha,Cof1,Cof2,Cof3,Cof4,
     x                      Pmat1,Pmat2,Pmat3,Pmat4,dfBoys)

               ans=bin_tx1*bin_tx2*bin_tx3*bin_tx4
     x            *bin_uy1*bin_uy2*bin_uy3*bin_uy4
     x            *bin_vz1*bin_vz2*bin_vz3*bin_vz4
     x            *dfKKK_x*dfKKK_y*dfKKK_z
     x            *dfBoys

CCWS--debug
c     dfKKK_z=
c     ans=0.944869744244220*0.862004514134267*0.87999629806098
               
c      write(*,*)'VEP MD dfBoys=',dfBoys
c      write(*,*)'VEP MD dfKKK_x=',dfKKK_x
c      write(*,*)'VEP MD dfKKK_y=',dfKKK_y
c      write(*,*)'VEP MD dfKKK_z=',dfKKK_z
               Dgf=Dgf+ans


            end do
            end do
            end do
            end do

         end do
         end do
         end do
         end do

      end do
      end do
      end do
      end do

CCWS--debug
c     Dgf=0.944869744244220*0.862004514134267*0.87999629806098

      return
      end

C=======================================================================
      subroutine G4xggs_DGK12(xb1,xb2,xb3,xb4,
     x                        xi1,xi2,
     x                        e1,f1,g1,h1,
     x                        e2,f2,g2,h2,
     x                        x1,x2,x3,x4,dGK12)

C Used for G4Vep integral calculation
C Derivative of product of two GK terms
C
C Calculate derivative P_x1,x2,x3,x4 (GK1 GK2)
C 
C = (d/P1x)^xb1 (d/dP2x)^xb2 (d/dP3x)^xb3 (d/dP4x)^xb4 GK1 GK2 
C
C=======================================================================
      implicit none

C Input variables
      integer xb1,xb2,xb3,xb4

      double precision xi1,xi2
      double precision e1,f1,g1,h1
      double precision e2,f2,g2,h2
      double precision x1,x2,x3,x4

C Variables returned
      double precision dGK12

C Local variables
      integer i1,i2,i3,i4

      double precision bin_i1
      double precision bin_i2
      double precision bin_i3
      double precision bin_i4

      double precision dGK1,dGK2
      double precision ans


      dGK12=0.0d+00

      do i1=0,xb1
         call cbinom(xb1,i1,bin_i1)
      do i2=0,xb2
         call cbinom(xb2,i2,bin_i2)
      do i3=0,xb3
         call cbinom(xb3,i3,bin_i3)
      do i4=0,xb4
         call cbinom(xb4,i4,bin_i4)

         call DKGx_4(xb1-i1,xb2-i2,xb3-i3,xb4-i4,
     x               xi1,e1,f1,g1,h1,x1,x2,x3,x4,dGK1)

         call DKGx_4(i1,i2,i3,i4,
     x               xi2,e2,f2,g2,h2,x1,x2,x3,x4,dGK2)

CCWS-debug
c        write(*,*)'dGK1=',dGK1
c        write(*,*)'dGK2=',dGK2
CCWS-debug
         ans=bin_i1*bin_i2*bin_i3*bin_i4
     x      *dGK1*dGK2
         dGK12=dGK12+ans

      end do
      end do
      end do
      end do


      return
      end

