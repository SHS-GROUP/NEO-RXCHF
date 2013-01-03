C======================================================================
      subroutine G3_MD_xggs(I1,J1,K1,A1,Amat1,
     *                       I2,J2,K2,A2,Amat2,
     *                       I3,J3,K3,A3,Amat3,
     *                       L1,M1,N1,B1,Bmat1,
     *                       L2,M2,N2,B2,Bmat2,
     *                       L3,M3,N3,B3,Bmat3,
     *                       gamA12,gamA13,gamA23,
     *                       gamB12,gamB13,gamB23,
     *                       xggs)

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

      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23

      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)

C Variables returned
      double precision xggs

C Local Variables
      double precision Ix,Iy,Iz


C Evaluate integrals over hermite gaussians
      call G3_MD_xggs1D(I1,A1,Amat1(1),
     *                  I2,A2,Amat2(1),
     *                  I3,A3,Amat3(1),
     *                  L1,B1,Bmat1(1),
     *                  L2,B2,Bmat2(1),
     *                  L3,B3,Bmat3(1),
     *                  gamA12,gamA13,gamA23,
     *                  gamB12,gamB13,gamB23,
     *                  Ix)

      call G3_MD_xggs1D(J1,A1,Amat1(2),
     *                  J2,A2,Amat2(2),
     *                  J3,A3,Amat3(2),
     *                  M1,B1,Bmat1(2),
     *                  M2,B2,Bmat2(2),
     *                  M3,B3,Bmat3(2),
     *                  gamA12,gamA13,gamA23,
     *                  gamB12,gamB13,gamB23,
     *                  Iy)

      call G3_MD_xggs1D(K1,A1,Amat1(3),
     *                  K2,A2,Amat2(3),
     *                  K3,A3,Amat3(3),
     *                  N1,B1,Bmat1(3),
     *                  N2,B2,Bmat2(3),
     *                  N3,B3,Bmat3(3),
     *                  gamA12,gamA13,gamA23,
     *                  gamB12,gamB13,gamB23,
     *                  Iz)


      xggs=Ix*Iy*Iz


      return
      end

C======================================================================
      subroutine G3_MD_xggs1D(I1,A1,Amat1,
     *                        I2,A2,Amat2,
     *                        I3,A3,Amat3,
     *                        L1,B1,Bmat1,
     *                        L2,B2,Bmat2,
     *                        L3,B3,Bmat3,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        ans)

C======================================================================
      implicit none

C     Input Variables
      integer I1,I2,I3
      integer L1,L2,L3
      double precision A1,A2,A3
      double precision B1,B2,B3
      double precision Amat1,Amat2,Amat3
      double precision Bmat1,Bmat2,Bmat3
      double precision gamA12,gamA13,gamA23
      double precision gamB12,gamB13,gamB23

C     Variables Returned
      double precision ans

C     Local Variables
      integer t1,t2,t3
      double precision zero
      double precision one
      double precision two
      double precision three
      double precision pi
      parameter(one=1.0d+00,zero=0.0d+00)
      parameter(two=2.0d+00,three=3.0d+00)
      parameter(pi=3.14159265358979d+00)
      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
      double precision Et1(0:I1+L1)
      double precision Et2(0:I2+L2)
      double precision Et3(0:I3+L3)
      double precision dfKKK
      double precision xans

      double precision gam13
      double precision gam23
      double precision gam12

      double precision KAB1
      double precision KAB2
      double precision KAB3

      double precision P1
      double precision P2
      double precision P3

      double precision Pmat1
      double precision Pmat2
      double precision Pmat3

      double precision qAB1
      double precision qAB2
      double precision qAB3

      double precision Qmat1
      double precision Qmat2
      double precision Qmat3

      double precision C0
      double precision C1
      double precision C2
      double precision C12
      double precision xsum
      double precision KP1P3
      double precision KP2P3
      double precision qPC1
      double precision qPC2
      double precision S1
      double precision S2
      double precision Smat1
      double precision Smat2
      double precision SQmat1
      double precision SQmat2

      double precision tempvar1
      double precision tempvar2
      double precision C3


      gam13=gamA13+gamB13
      gam23=gamA23+gamB23
      gam12=gamA12+gamB12

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod_1D(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod_1D(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
      call Gauss_prod_1D(A3,B3,Amat3,Bmat3,KAB3,qAB3,P3,Pmat3,Qmat3)

C  Calculate Hermite expansion coefficients:
      call ghec(I1,L1,KAB1,A1,B1,P1,qAB1,Qmat1,Et1)
      call ghec(I2,L2,KAB2,A2,B2,P2,qAB2,Qmat2,Et2)
      call ghec(I3,L3,KAB3,A3,B3,P3,qAB3,Qmat3,Et3)

C  Apply Boys' Lemma 2 to particle 3:
C     (
C IN: /dx3 Exp(-P3(x3-Pmat3)**2) Exp(-g13(x3-x1)**2) Exp(-g23(x3-x2)**2)
C     )
C OUT: C0*Exp(-C1(x1-Pmat3)**2 Exp(-C2(x2-Pmat3)**2 Exp(-C12(x1-x2)**2
c     C1=(P3*gam13)/xsum
c     C2=(P3*gam23)/xsum
c     C12=(gam13*gam23)/xsum

c     call SFBoys2(P3,gam13,gam23,Pmat3,
c    x             C0,C1,C2,C12)
      xsum=P3+gam13+gam23

c     C0=(pi/xsum)**(three/two)
      C0=sqrt(pi/xsum)

      C1=(P3*gam13)/xsum
      C2=(P3*gam23)/xsum
      C12=(gam13*gam23)/xsum


C  Update gamma(1,2):
      gam12=gam12+C12

CCWS-DEBUG(
c     write(*,*)'gam12=',gam12
c     write(*,*)'C1=',C1
c     write(*,*)'C12=',C12
CCWS-DEBUG)
C  Should now have the following integral:
C  (
C  | Exp(-C1(x1-Pmat3)**2) Exp(-C2(x2-Pmat3)**2) Exp(-gam12(x1-x2)**2)
C  | Exp(-P1(x1-Pmat1)**2) Exp(-P2(x2-Pmat2)**2)  
C  )

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod_1D(P1,C1,Pmat1,Pmat3,KP1P3,qPC1,S1,Smat1,SQmat1)
      call Gauss_prod_1D(P2,C2,Pmat2,Pmat3,KP2P3,qPC2,S2,Smat2,SQmat2)

C  Should now have the following integral:
C  (
C  | Exp(-S1(x1-Smat1)**2) Exp(-S2(x2-Smat2)**2) Exp(-gam12(x1-x2)**2)
C  )

C  Evaluate integral over spherical Gaussians:
      tempvar1=(S1*S2+S1*gam12+S2*gam12)
      tempvar2=Pi*Pi/tempvar1
c     C3=tempvar2*sqrt(tempvar2)
      C3=sqrt(tempvar2)

C  Organize terms for differentiation of f_KKK:
C Prepare the following KKK Gaussian functions:
C KP1P3 <=> GK1 = exp(-xi1*( f1*x1 + g1*x2 + h1*x3 )
C KP2P3 <=> GK2 = exp(-xi2*( f2*x1 + g2*x2 + h2*x3 )
C KSS   <=> GK3 = exp(-xi3*( f3*x1 + g3*x2 + h3*x3 )

      xi1=p1*c1/(p1+c1)
      f1=one
      g1=zero
      h1=-one

      xi2=p2*c2/(p2+c2)
      f2=zero
      g2=one
      h2=-one

      xi3=s1*s2*gam12/(s1*s2+s1*gam12+s2*gam12)
      f3=p1/(p1+c1)
      g3=-p2/(p2+c2)
      h3=c1/(p1+c1) - c2/(p2+c2)


      ans=0.0d+00

      do t1=0,I1+L1
         do t2=0,I2+L2
            do t3=0,I3+L3

               call G3xggs_DGK123(t1,t2,t3,
     x                            xi1,xi2,xi3,
     x                            f1,g1,h1,
     x                            f2,g2,h2,
     x                            f3,g3,h3,
     x                            Pmat1,Pmat2,Pmat3,dfKKK)

               xans=Et1(t1)*Et2(t2)*Et3(t3)*dfKKK

               ans=ans+xans

            end do
         end do
      end do

      ans=ans*c0*c3

      return
      end

C======================================================================
c     subroutine cws_gam2_xggs(I1,J1,K1,A1,Amat1,
      subroutine XG3_MD_xggs(I1,J1,K1,A1,Amat1,
     *                      I2,J2,K2,A2,Amat2,
     *                      I3,J3,K3,A3,Amat3,
     *                      L1,M1,N1,B1,Bmat1,
     *                      L2,M2,N2,B2,Bmat2,
     *                      L3,M3,N3,B3,Bmat3,
     *                      gamA12,gamA13,gamA23,
     *                      gamB12,gamB13,gamB23,
     *                      xggs)

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

      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23

      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)

C Variables returned
      double precision xggs

C Local variables
      integer t1,u1,v1
      integer t2,u2,v2
      integer t3,u3,v3
c     integer Nt,Nu,Nv

      double precision zero
      double precision half
      double precision one
      double precision two
      double precision three
      double precision five
      double precision pi
      parameter(two=2.0d+00,three=3.0d+00,pi=3.14159265358979d+00)
      parameter(half=0.5d+00,one=1.0d+00,five=5.0d+00,zero=0.0d+00)

      double precision gam13
      double precision gam23
      double precision gam12
      double precision gam12_save

      double precision KAB1(3)
      double precision KAB2(3)
      double precision KAB3(3)

      double precision P1
      double precision P2
      double precision P3

      double precision Pmat1(3)
      double precision Pmat2(3)
      double precision Pmat3(3)

      double precision qAB1
      double precision qAB2
      double precision qAB3

      double precision Qmat1(3)
      double precision Qmat2(3)
      double precision Qmat3(3)

      double precision ansg1
      double precision C0
      double precision C1
      double precision C2
      double precision C12

      double precision Et1(0:I1+L1)
      double precision Eu1(0:J1+M1)
      double precision Ev1(0:K1+N1)

      double precision Et2(0:I2+L2)
      double precision Eu2(0:J2+M2)
      double precision Ev2(0:K2+N2)

      double precision Et3(0:I3+L3)
      double precision Eu3(0:J3+M3)
      double precision Ev3(0:K3+N3)

      double precision KP1P3
      double precision KP2P3
      double precision qPC1
      double precision qPC2
      double precision S1
      double precision S2
      double precision Smat1(3)
      double precision Smat2(3)
      double precision SQmat1(3)
      double precision SQmat2(3)

      double precision tempvar1
      double precision tempvar2
      double precision C3
      double precision KSS
c     double precision alpha
      double precision R2_S1S2
c     double precision expT

      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
      double precision ans
      double precision DfKKK



      gam13=gamA13+gamB13
      gam23=gamA23+gamB23
      gam12=gamA12+gamB12

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
      call Gauss_prod(A3,B3,Amat3,Bmat3,KAB3,qAB3,P3,Pmat3,Qmat3)

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

C  Apply Boys' Lemma 2 to particle 3:
C     (
C IN: |dx3 Exp(-P3(x3-Pmat3)**2) Exp(-g13(x3-x1)**2) Exp(-g23(x3-x2)**2)
C     )
C OUT: C0*Exp(-C1(x1-Pmat3)**2 Exp(-C2(x2-Pmat3)**2 Exp(-C12(x1-x2)**2
c     C1=(P3*gam13)/xsum
c     C2=(P3*gam23)/xsum
c     C12=(gam13*gam23)/xsum

      call SFBoys2(P3,gam13,gam23,Pmat3,
     x             C0,C1,C2,C12)

C  Update gamma(1,2):
c     gam12=gam12_save
      gam12=gam12+C12

C  Should now have the following integral:
C  (
C  | Exp(-C1(x1-Pmat3)**2) Exp(-C2(x2-Pmat3)**2) Exp(-gam12(x1-x2)**2)
C  | Exp(-P1(x1-Pmat1)**2) Exp(-P2(x2-Pmat2)**2)  
C  )

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(P1,C1,Pmat1,Pmat3,KP1P3,qPC1,S1,Smat1,SQmat1)
      call Gauss_prod(P2,C2,Pmat2,Pmat3,KP2P3,qPC2,S2,Smat2,SQmat2)

C  Should now have the following integral:
C  (
C  | Exp(-S1(x1-Smat1)**2) Exp(-S2(x2-Smat2)**2) Exp(-gam12(x1-x2)**2)
C  | dx1dx2 
C  )

C  Evaluate integral over spherical Gaussians:
      tempvar1=(S1*S2+S1*gam12+S2*gam12)
      C3=(pi*pi)/tempvar1
      C3=C3*sqrt(C3)

      R2_S1S2=(Smat1(1)-Smat2(1))**2 +
     x        (Smat1(2)-Smat2(2))**2 +
     x        (Smat1(3)-Smat2(3))**2 

      KSS=exp(((-S1*S2*gam12)/tempvar1)*R2_S1S2)

C  Organize terms for differentiation:
C Prepare the following KKK Gaussian functions:
C KP1P3 <=> GK1 = exp(-xi1*( f1*x1 + g1*x2 + h1*x3 )
C KP2P3 <=> GK2 = exp(-xi2*( f2*x1 + g2*x2 + h2*x3 )
C KSS   <=> GK3 = exp(-xi3*( f3*x1 + g3*x2 + h3*x3 )
      call G3xggs_dfKKK_prep(p1,c1,p2,c2,s1,s2,gam12,
     x                       xi1,xi2,xi3,
     x                       f1,g1,h1,
     x                       f2,g2,h2,
     x                       f3,g3,h3)
C Evaluate integral over hermite gaussians

      xggs=0.0d+00
c     gam12_save=gam12

      do t1=0,I1+L1
      do u1=0,J1+M1
      do v1=0,K1+N1

         do t2=0,I2+L2
         do u2=0,J2+M2
         do v2=0,K2+N2

            do t3=0,I3+L3
            do u3=0,J3+M3
            do v3=0,K3+N3


C  Differentiation of integral over spherical Gaussians:               
      call G3xggs_dfKKK(t1,t2,t3,
     x                  u1,u2,u3,
     x                  v1,v2,v3,
     x                  xi1,xi2,xi3,
     x                  f1,g1,h1,
     x                  f2,g2,h2,
     x                  f3,g3,h3,
     x                  Pmat1,Pmat2,Pmat3,DfKKK)

      ans=C0*C3*
     x   Et1(t1)*Et2(t2)*Et3(t3) *
     x   Eu1(u1)*Eu2(u2)*Eu3(u3) *
     x   Ev1(v1)*Ev2(v2)*Ev3(v3) *
     x   DfKKK

      xggs=xggs+ans
 


            end do
            end do
            end do

         end do
         end do
         end do

      end do
      end do
      end do
CCCCC---CHOP
C  Now have:
C KAB1 Exp(-P1(x1-Pmat1)**2) 
C KAB2 Exp(-P2(x2-Pmat2)**2) 
C KAB3 Exp(-P3(x3-Pmat3)**2)


c     write(*,*)
c     write(*,*)
c     write(*,*)
c     write(*,*)'KAB1=',KAB1
c     write(*,*)'KAB2=',KAB2
c     write(*,*)'KAB3=',KAB3
c     write(*,*)'P1=',P1
c     write(*,*)'P2=',P2
c     write(*,*)'P3=',P3
c     write(*,*)'C0=',C0
c     write(*,*)'C1=',C1
c     write(*,*)'C2=',C2
c     write(*,*)'C12=',C12
c     write(*,*)'gam12=',gam12
c     write(*,*)'gam13=',gam13
c     write(*,*)'gam23=',gam23
c     write(*,*)'ZNUC=',ZNUC
c     write(*,*)'ansg1=',ansg1
c     write(*,*)'xgVxCg=',xgVxCg



      return
      end

C=======================================================================
      subroutine G3xggs_dfKKK_prep(p1,c1,p2,c2,s1,s2,gam12,
     x                             xi1,xi2,xi3,
     x                             f1,g1,h1,
     x                             f2,g2,h2,
     x                             f3,g3,h3)
c    x                             Pmat1,Pmat2,Pmat3,
c    x                             Rmat)

C Prepare the following KKK Gaussian functions:
C KP1P3 <=> GK1 = exp(-xi1*( f1*x1 + g1*x2 + h1*x3 )
C KP2P3 <=> GK2 = exp(-xi2*( f2*x1 + g2*x2 + h2*x3 )
C KSS   <=> GK3 = exp(-xi3*( f3*x1 + g3*x2 + h3*x3 )
C=======================================================================
      implicit none
C Input variables
      double precision p1,c1
      double precision p2,c2
      double precision s1
      double precision s2
      double precision gam12
c     double precision Pmat1(3),Pmat2(3),Pmat3(3)
C Variables returned
      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
c     double precision Rmat(3)
C Local Variables

C     Gaussian K...: GK...
C     GK1 <=> KP1P3:
      xi1=(p1*c1)/(p1+c1)
      f1=1.0d+00
      g1=0.0d+00
      h1=-1.0d+00
      
C     GK2 <=> KP2P3:
      xi2=(p2*c2)/(p2+c2)
      f2=0.0d+00
      g2=1.0d+00
      h2=-1.0d+00
      
C     GK3 <=> KSS:
      xi3=(s1*s2*gam12)/(s1*s2+s1*gam12+s2*gam12)
      f3=p1/(p1+c1)
      g3=-p2/(p2+c2)
      h3=(c1/(p1+c1))-(c2/(p2+c2))
     
C Rmat needed for Boys func derivatives
c     double precision Pmat1(3),Pmat2(3),Pmat3(3)
c     double precision Rmat(3)
c     Rmat(1)=(f3*Pmat1(1))+(g3*Pmat2(1))+(h3*Pmat3(1)) 
c     Rmat(2)=(f3*Pmat1(2))+(g3*Pmat2(2))+(h3*Pmat3(2)) 
c     Rmat(3)=(f3*Pmat1(3))+(g3*Pmat2(3))+(h3*Pmat3(3)) 


      return
      end

C=======================================================================
c     subroutine gam2_vee_dhermite(t1,t2,t3,
      subroutine G3xggs_dfKKK(t1,t2,t3,
     x                        u1,u2,u3,
     x                        v1,v2,v3,
     x                        xi1,xi2,xi3,
     x                        f1,g1,h1,
     x                        f2,g2,h2,
     x                        f3,g3,h3,
     x                        Pmat1,Pmat2,Pmat3,DfKKK)

C=======================================================================
      implicit none

C Input Variables
      integer t1,u1,v1
      integer t2,u2,v2
      integer t3,u3,v3
      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
      double precision Pmat1(3)
      double precision Pmat2(3)
      double precision Pmat3(3)
C Variables Returned
      double precision DfKKK 
C Local Variables
      double precision DfKKK_x 
      double precision DfKKK_y 
      double precision DfKKK_z 



      call G3xggs_DGK123(t1,t2,t3,
     x                   xi1,xi2,xi3,
     x                   f1,g1,h1,
     x                   f2,g2,h2,
     x                   f3,g3,h3,
     x                   Pmat1(1),Pmat2(1),Pmat3(1),dfKKK_x)

      call G3xggs_DGK123(u1,u2,u3,
     x                   xi1,xi2,xi3,
     x                   f1,g1,h1,
     x                   f2,g2,h2,
     x                   f3,g3,h3,
     x                   Pmat1(2),Pmat2(2),Pmat3(2),dfKKK_y)

      call G3xggs_DGK123(v1,v2,v3,
     x                   xi1,xi2,xi3,
     x                   f1,g1,h1,
     x                   f2,g2,h2,
     x                   f3,g3,h3,
     x                   Pmat1(3),Pmat2(3),Pmat3(3),dfKKK_z)

      DfKKK=DfKKK_x*DfKKK_y*DfKKK_z


      return
      end


C======================================================================
c     subroutine cws_gam2_xgVxCg(I1,J1,K1,A1,Amat1,
      subroutine G3_MD_xgVxCg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        Cmat,ZNUC,
     *                        xgVxCg)

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

      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23

      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)

      double precision Cmat(3)
      double precision ZNUC

C Variables returned
      double precision xgVxCg

C Local variables
      integer t1,u1,v1
      integer t2,u2,v2
      integer t3,u3,v3
c     integer Nt,Nu,Nv

      double precision zero
      double precision half
      double precision one
      double precision two
      double precision three
      double precision five
      double precision pi
      parameter(two=2.0d+00,three=3.0d+00,pi=3.14159265358979d+00)
      parameter(half=0.5d+00,one=1.0d+00,five=5.0d+00,zero=0.0d+00)

      double precision gam13
      double precision gam23
      double precision gam12
      double precision gam12_save

      double precision KAB1(3)
      double precision KAB2(3)
      double precision KAB3(3)

      double precision P1
      double precision P2
      double precision P3

      double precision Pmat1(3)
      double precision Pmat2(3)
      double precision Pmat3(3)

      double precision qAB1
      double precision qAB2
      double precision qAB3

      double precision Qmat1(3)
      double precision Qmat2(3)
      double precision Qmat3(3)

      double precision ansg1
      double precision C0
      double precision C1
      double precision C2
      double precision C12

      double precision Et1(0:I1+L1)
      double precision Eu1(0:J1+M1)
      double precision Ev1(0:K1+N1)

      double precision Et2(0:I2+L2)
      double precision Eu2(0:J2+M2)
      double precision Ev2(0:K2+N2)

      double precision Et3(0:I3+L3)
      double precision Eu3(0:J3+M3)
      double precision Ev3(0:K3+N3)

      double precision KP1P3(3)
      double precision KP2P3(3)
      double precision qPC1
      double precision qPC2
      double precision S1
      double precision S2
      double precision Smat(3)
      double precision Smat1(3)
      double precision Smat2(3)
      double precision SQmat1(3)
      double precision SQmat2(3)

      double precision tempvar1
      double precision tempvar2
      double precision tempvar3
      double precision C3
      double precision KSS
      double precision alpha
      double precision R2_S1S2
      double precision expT

      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
      double precision ans
      double precision Dgf

      double precision RCS_f,RCS_g,RCS_h
      double precision R2_CS


      gam13=gamA13+gamB13
      gam23=gamA23+gamB23
      gam12=gamA12+gamB12

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
      call Gauss_prod(A3,B3,Amat3,Bmat3,KAB3,qAB3,P3,Pmat3,Qmat3)

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

C Evaluate integral over hermite gaussians

      xgVxCg=0.0d+00

      gam12_save=gam12

      do t1=0,I1+L1
      do u1=0,J1+M1
      do v1=0,K1+N1

         do t2=0,I2+L2
         do u2=0,J2+M2
         do v2=0,K2+N2

            do t3=0,I3+L3
            do u3=0,J3+M3
            do v3=0,K3+N3

C  Apply Boys' Lemma 2 to particle 3:
C     (
C IN: |dx3 Exp(-P3(x3-Pmat3)**2) Exp(-g13(x3-x1)**2) Exp(-g23(x3-x2)**2)
C     )
C OUT: C0*Exp(-C1(x1-Pmat3)**2 Exp(-C2(x2-Pmat3)**2 Exp(-C12(x1-x2)**2
c     C1=(P3*gam13)/xsum
c     C2=(P3*gam23)/xsum
c     C12=(gam13*gam23)/xsum

      call SFBoys2(P3,gam13,gam23,Pmat3,
     x             C0,C1,C2,C12)

C  Update gamma(1,2):
      gam12=gam12_save
      gam12=gam12+C12

CCWS-DEBUG(
c     write(*,*)'gam12=',gam12
c     write(*,*)'C1=',C1
c     write(*,*)'C12=',C12
CCWS-DEBUG)
C  Should now have the following integral:
C  (
C  | Exp(-C1(x1-Pmat3)**2) Exp(-C2(x2-Pmat3)**2) Exp(-gam12(x1-x2)**2)
C  | Exp(-P1(x1-Pmat1)**2) Exp(-P2(x2-Pmat2)**2) 1/|C-x1| 
C  )

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(P1,C1,Pmat1,Pmat3,KP1P3,qPC1,S1,Smat1,SQmat1)
      call Gauss_prod(P2,C2,Pmat2,Pmat3,KP2P3,qPC2,S2,Smat2,SQmat2)

C  Should now have the following integral:
C  (
C  | Exp(-S1(x1-Smat1)**2) Exp(-S2(x2-Smat2)**2) Exp(-gam12(x1-x2)**2)
C  | 1/|C-x1| 
C  )

C  Evaluate integral over spherical Gaussians:
      tempvar1=(s2+gam12)
      tempvar2=(S1*S2+S1*gam12+S2*gam12)
      tempvar3=tempvar2/tempvar1

      R2_S1S2=(Smat1(1)-Smat2(1))**2 +
     x        (Smat1(2)-Smat2(2))**2 +
     x        (Smat1(3)-Smat2(3))**2 

      C3=sqrt((pi/tempvar1)**3)*(2.0d+00*pi/tempvar3)

      KSS=exp(((-S1*S2*gam12)/tempvar2)*R2_S1S2)

      Smat(1)=((s1*s2+s1*gam12)/tempvar2)*Smat1(1)
     x             +((s2*gam12)/tempvar2)*Smat2(1)
      Smat(2)=((s1*s2+s1*gam12)/tempvar2)*Smat1(2)
     x             +((s2*gam12)/tempvar2)*Smat2(2)
      Smat(3)=((s1*s2+s1*gam12)/tempvar2)*Smat1(3)
     x             +((s2*gam12)/tempvar2)*Smat2(3)

c     R2_CS=(Smat(1)-Cmat(1))**2 +
c    x      (Smat(2)-Cmat(2))**2 +
c    x      (Smat(3)-Cmat(3))**2 
      R2_CS=(Cmat(1)-Smat(1))**2 +
     x      (Cmat(2)-Smat(2))**2 +
     x      (Cmat(3)-Smat(3))**2 

      alpha=tempvar3

      expT=alpha*R2_CS

     
C  Boys RCS f,g,h:
      RCS_f=(s1*s2+s1*gam12)/(s1*s2+s1*gam12+s2*gam12)*(p1/(p1+c1))
      RCS_g=(s2*gam12)/(s1*s2+s1*gam12+s2*gam12)*(p2/(p2+c2))
c     RCS_h=(s1*s2+s1*gam12)/(s1*s2+s1*gam12+s2*gam12)*(p1/(p1+c1))
c    x     +(s2*gam12)/(s1*s2+s1*gam12+s2*gam12)*(p2/(p2+c2))
      RCS_h=(s1*s2+s1*gam12)/(s1*s2+s1*gam12+s2*gam12)*(c1/(p1+c1))
     x     +(s2*gam12)/(s1*s2+s1*gam12+s2*gam12)*(c2/(p2+c2))

CCWS-Debug
c     R2_CS=(Cmat(1)-(RCS_f*Pmat1(1)+RCS_g*Pmat2(2)+RCS_h*Pmat3(1)))**2
c    x     +(Cmat(2)-(RCS_f*Pmat1(2)+RCS_g*Pmat2(2)+RCS_h*Pmat3(2)))**2
c    x     +(Cmat(3)-(RCS_f*Pmat1(3)+RCS_g*Pmat2(3)+RCS_h*Pmat3(3)))**2
c     alpha=tempvar3
c     expT=alpha*R2_CS
CCWS-Debug

C  Organize terms for differentiation of f_KKK:
C Prepare the following KKK Gaussian functions:
C KP1P3 <=> GK1 = exp(-xi1*( f1*x1 + g1*x2 + h1*x3 )
C KP2P3 <=> GK2 = exp(-xi2*( f2*x1 + g2*x2 + h2*x3 )
C KSS   <=> GK3 = exp(-xi3*( f3*x1 + g3*x2 + h3*x3 )
      call vec_dhermite_prep(p1,c1,p2,c2,s1,s2,gam12,
     x                       xi1,xi2,xi3,
     x                       f1,g1,h1,
     x                       f2,g2,h2,
     x                       f3,g3,h3)

C  Differentiation of integral over spherical Gaussians:               
c     call gam2_vec_dhermite(t1,t2,t3,

CCWS-DEBUG
c     write(*,*)''
c     write(*,*)'>>> V3=',V3,'<<<'
c     write(*,*)'    T1=',T1
c     write(*,*)'    U1=',U1
c     write(*,*)'    V1=',V1
c     write(*,*)'    T2=',T2
c     write(*,*)'    U2=',U2
c     write(*,*)'    V2=',V2
c     write(*,*)'    T3=',T3
c     write(*,*)'    U3=',U3
c     write(*,*)'    V3=',V3
c     write(*,*)'alpha=',alpha
c     write(*,*)'RCS_f=',RCS_f
c     write(*,*)'RCS_g=',RCS_g
c     write(*,*)'RCS_h=',RCS_h
c     write(*,*)'xi1=',xi1
c     write(*,*)'xi2=',xi2
c     write(*,*)'xi3=',xi3
c     write(*,*)'f1=',f1
c     write(*,*)'g1=',g1
c     write(*,*)'h1=',h1
c     write(*,*)'f2=',f2
c     write(*,*)'g2=',g2
c     write(*,*)'h2=',h2
c     write(*,*)'f3=',f3
c     write(*,*)'g3=',g3
c     write(*,*)'h3=',h3
c     write(*,*)''
CCWS-DEBUG

      call g3_vec_dhermite(t1,t2,t3,
     x                     u1,u2,u3,
     x                     v1,v2,v3,
     x                     xi1,xi2,xi3,
     x                     f1,g1,h1,
     x                     f2,g2,h2,
     x                     f3,g3,h3,
     x                     RCS_f,RCS_g,RCS_h,
     x                     expT,alpha,
     x                     Pmat1,Pmat2,Pmat3,Cmat,Dgf)

      ans=C0*C3*
     x   Et1(t1)*Et2(t2)*Et3(t3) *
     x   Eu1(u1)*Eu2(u2)*Eu3(u3) *
     x   Ev1(v1)*Ev2(v2)*Ev3(v3) *
     x   Dgf

      xgVxCg=xgVxCg+ans
 


            end do
            end do
            end do

         end do
         end do
         end do

      end do
      end do
      end do
CCWS-DEBUG
c     write(*,*)'end -of- integral'
CCWS-DEBUG
CCCCC---CHOP
C  Now have:
C KAB1 Exp(-P1(x1-Pmat1)**2) 
C KAB2 Exp(-P2(x2-Pmat2)**2) 
C KAB3 Exp(-P3(x3-Pmat3)**2)


c     write(*,*)
c     write(*,*)
c     write(*,*)
c     write(*,*)'KAB1=',KAB1
c     write(*,*)'KAB2=',KAB2
c     write(*,*)'KAB3=',KAB3
c     write(*,*)'P1=',P1
c     write(*,*)'P2=',P2
c     write(*,*)'P3=',P3
c     write(*,*)'C0=',C0
c     write(*,*)'C1=',C1
c     write(*,*)'C2=',C2
c     write(*,*)'C12=',C12
c     write(*,*)'gam12=',gam12
c     write(*,*)'gam13=',gam13
c     write(*,*)'gam23=',gam23
c     write(*,*)'ZNUC=',ZNUC
c     write(*,*)'ansg1=',ansg1
c     write(*,*)'xgVxCg=',xgVxCg



      return
      end

C=======================================================================
      subroutine vec_dhermite_prep(p1,c1,p2,c2,s1,s2,gam12,
     x                             xi1,xi2,xi3,
     x                             f1,g1,h1,
     x                             f2,g2,h2,
     x                             f3,g3,h3)
c    x                             Pmat1,Pmat2,Pmat3,
c    x                             Rmat)

C Prepare the following KKK Gaussian functions:
C KP1P3 <=> GK1 = exp(-xi1*( f1*x1 + g1*x2 + h1*x3 )
C KP2P3 <=> GK2 = exp(-xi2*( f2*x1 + g2*x2 + h2*x3 )
C KSS   <=> GK3 = exp(-xi3*( f3*x1 + g3*x2 + h3*x3 )
C=======================================================================
      implicit none
C Input variables
      double precision p1,c1
      double precision p2,c2
      double precision s1
      double precision s2
      double precision gam12
c     double precision Pmat1(3),Pmat2(3),Pmat3(3)
C Variables returned
      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
c     double precision Rmat(3)
C Local Variables

C     Gaussian K...: GK...
C     GK1 <=> KP1P3:
      xi1=(p1*c1)/(p1+c1)
      f1=1.0d+00
      g1=0.0d+00
      h1=-1.0d+00
      
C     GK2 <=> KP2P3:
      xi2=(p2*c2)/(p2+c2)
      f2=0.0d+00
      g2=1.0d+00
      h2=-1.0d+00
      
C     GK3 <=> KSS:
      xi3=(s1*s2*gam12)/(s1*s2+s1*gam12+s2*gam12)
      f3=p1/(p1+c1)
      g3=-p2/(p2+c2)
      h3=(c1/(p1+c1))-(c2/(p2+c2))

C Rmat needed for Boys func derivatives
c     double precision Pmat1(3),Pmat2(3),Pmat3(3)
c     double precision Rmat(3)
c     Rmat(1)=(f3*Pmat1(1))+(g3*Pmat2(1))+(h3*Pmat3(1)) 
c     Rmat(2)=(f3*Pmat1(2))+(g3*Pmat2(2))+(h3*Pmat3(2)) 
c     Rmat(3)=(f3*Pmat1(3))+(g3*Pmat2(3))+(h3*Pmat3(3)) 


      return
      end

C=======================================================================
c     subroutine gam2_vec_dhermite(t1,t2,t3,
      subroutine g3_vec_dhermite(t1,t2,t3,
     x                           u1,u2,u3,
     x                           v1,v2,v3,
     x                           xi1,xi2,xi3,
     x                           f1,g1,h1,
     x                           f2,g2,h2,
     x                           f3,g3,h3,
     x                           RCS_f,RCS_g,RCS_h,
     x                           expT,alpha,
     x                           Pmat1,Pmat2,Pmat3,Cmat,Dgf)

C=======================================================================
      implicit none

C Input variables
      integer t1,t2,t3
      integer u1,u2,u3
      integer v1,v2,v3

      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
      double precision RCS_f,RCS_g,RCS_h
      double precision expT
      double precision alpha
      double precision Pmat1(3),Pmat2(3),Pmat3(3) 
      double precision Cmat(3) 

C Variables returned
      double precision Dgf

C Local Variables
      integer xb1,xb2,xb3
      integer yb1,yb2,yb3
      integer zb1,zb2,zb3

      double precision bin_tx1 
      double precision bin_tx2
      double precision bin_tx3
      double precision bin_uy1
      double precision bin_uy2
      double precision bin_uy3
      double precision bin_vz1
      double precision bin_vz2
      double precision bin_vz3
      double precision dfKKK_x 
      double precision dfKKK_y 
      double precision dfKKK_z 
      double precision dfBoys
      double precision ans 
      

      Dgf=0.0d+00

      do xb1=0,t1
         call cbinom(t1,xb1,bin_tx1)
      do xb2=0,t2
         call cbinom(t2,xb2,bin_tx2)
      do xb3=0,t3
         call cbinom(t3,xb3,bin_tx3)

c        call DKG123(t1-xb1,t2-xb2,t3-xb3,
         call G3xggs_DGK123(t1-xb1,t2-xb2,t3-xb3,
     x                      xi1,xi2,xi3,
     x                      f1,g1,h1,
     x                      f2,g2,h2,
     x                      f3,g3,h3,
     x                      Pmat1(1),Pmat2(1),Pmat3(1),dfKKK_x)

         do yb1=0,u1
            call cbinom(u1,yb1,bin_uy1)
         do yb2=0,u2
            call cbinom(u2,yb2,bin_uy2)
         do yb3=0,u3
            call cbinom(u3,yb3,bin_uy3)

c           call DKG123(u1-yb1,u2-yb2,u3-yb3,
            call G3xggs_DGK123(u1-yb1,u2-yb2,u3-yb3,
     x                         xi1,xi2,xi3,
     x                         f1,g1,h1,
     x                         f2,g2,h2,
     x                         f3,g3,h3,
     x                         Pmat1(2),Pmat2(2),Pmat3(2),dfKKK_y)


            do zb1=0,v1
               call cbinom(v1,zb1,bin_vz1)
            do zb2=0,v2
               call cbinom(v2,zb2,bin_vz2)
            do zb3=0,v3
               call cbinom(v3,zb3,bin_vz3)

c              call DKG123(v1-zb1,v2-zb2,v3-zb3,
               call G3xggs_DGK123(v1-zb1,v2-zb2,v3-zb3,
     x                            xi1,xi2,xi3,
     x                            f1,g1,h1,
     x                            f2,g2,h2,
     x                            f3,g3,h3,
     x                            Pmat1(3),Pmat2(3),Pmat3(3),dfKKK_z)

C
c              call gam2_VxC_RTUV(xb1,yb1,zb1,
               call g3_VxC_RTUV(xb1,yb1,zb1,
     x                            xb2,yb2,zb2,
     x                            xb3,yb3,zb3,
     x                            expT,alpha,RCS_f,RCS_g,RCS_h,
     x                            Pmat1,Pmat2,Pmat3,Cmat,dfBoys)
C
CCWS-DEBUG
c              write(*,*)'zb1=',zb1
c              write(*,*)'bin_vz1=',bin_vz1
c              write(*,*)'dfKKK_z=',dfKKK_z
c              write(*,*)'dfBoys=',dfBoys
CCWS-DEBUG

               ans=bin_tx1*bin_tx2*bin_tx3
     x            *bin_uy1*bin_uy2*bin_uy3
     x            *bin_vz1*bin_vz2*bin_vz3
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

CCWS-DEBUG
c              write(*,*)'Dgf=',Dgf
CCWS-DEBUG

      return
      end


C=======================================================================
c     subroutine gam2_VxC_RTUV(xb1,yb1,zb1,
      subroutine g3_VxC_RTUV(xb1,yb1,zb1,
     x                       xb2,yb2,zb2,
     x                       xb3,yb3,zb3,
     x                       expT,alpha,f,g,h,
     x                       Pmat1,Pmat2,Pmat3,Cmat,dfBoys)
C
C
C=======================================================================
      implicit none

C Input variables
      integer xb1,yb1,zb1 
      integer xb2,yb2,zb2 
      integer xb3,yb3,zb3 

      double precision expT
      double precision alpha
      double precision f,g,h
      double precision Pmat1(3),Pmat2(3),Pmat3(3)
      double precision Cmat(3)

C Variables returned
      double precision dfBoys

C Local variables
      integer NRx,NRy,NRz
      integer N,L,M
      integer nLim 

      double precision c1
      double precision ans
      double precision Rmat(3)


C Rmat needed for Boys func derivatives
      Rmat(1)=(f*Pmat1(1))+(g*Pmat2(1))+(h*Pmat3(1))-Cmat(1) 
      Rmat(2)=(f*Pmat1(2))+(g*Pmat2(2))+(h*Pmat3(2))-Cmat(2)
      Rmat(3)=(f*Pmat1(3))+(g*Pmat2(3))+(h*Pmat3(3))-Cmat(3)
c     Rmat(1)=Cmat(1)-(f*Pmat1(1))+(g*Pmat2(1))+(h*Pmat3(1)) 
c     Rmat(2)=Cmat(2)-(f*Pmat1(2))+(g*Pmat2(2))+(h*Pmat3(2))
c     Rmat(3)=Cmat(3)-(f*Pmat1(3))+(g*Pmat2(3))+(h*Pmat3(3))

      NRx=xb1+yb1+zb1
      NRy=xb2+yb2+zb2
      NRz=xb3+yb3+zb3

      N=xb1+xb2+xb3
      L=yb1+yb2+yb3
      M=zb1+zb2+zb3

      c1=(f**NRx)*(g**NRy)*(h**NRz)
c     c1=-1.0d+00*(-f**NRx)*(g**NRy)*(h**NRz)
CCWS-DEBUG
c     write(*,*)'f^NRx=',f**NRx
c     write(*,*)'c1=',c1
c     write(*,*)' M=',M
CCWS-DEBUG

c     nLim=NRx+NRy+NRz
      nLim=N+L+M

c     call RTUV(NRx,NRy,NRz,nLim,
c    x          expT,alpha,Rmat(1),Rmat(2),Rmat(3),ans)
      call RTUV(N,L,M,nLim,
     x          expT,alpha,Rmat(1),Rmat(2),Rmat(3),ans)

      dfBoys=c1*ans



      return
      end


C======================================================================
c     subroutine cws_gam2_xgVeeg(I1,J1,K1,A1,Amat1,
      subroutine G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        xgVeeg)

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

      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23

      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)

C Variables returned
      double precision xgVeeg

C Local variables
      integer t1,u1,v1
      integer t2,u2,v2
      integer t3,u3,v3
c     integer Nt,Nu,Nv

      double precision zero
      double precision half
      double precision one
      double precision two
      double precision three
      double precision five
      double precision pi
      parameter(two=2.0d+00,three=3.0d+00,pi=3.14159265358979d+00)
      parameter(half=0.5d+00,one=1.0d+00,five=5.0d+00,zero=0.0d+00)

      double precision gam13
      double precision gam23
      double precision gam12
      double precision gam12_save

      double precision KAB1(3)
      double precision KAB2(3)
      double precision KAB3(3)

      double precision P1
      double precision P2
      double precision P3

      double precision Pmat1(3)
      double precision Pmat2(3)
      double precision Pmat3(3)

      double precision qAB1
      double precision qAB2
      double precision qAB3

      double precision Qmat1(3)
      double precision Qmat2(3)
      double precision Qmat3(3)

      double precision ansg1
      double precision C0
      double precision C1
      double precision C2
      double precision C12

      double precision Et1(0:I1+L1)
      double precision Eu1(0:J1+M1)
      double precision Ev1(0:K1+N1)

      double precision Et2(0:I2+L2)
      double precision Eu2(0:J2+M2)
      double precision Ev2(0:K2+N2)

      double precision Et3(0:I3+L3)
      double precision Eu3(0:J3+M3)
      double precision Ev3(0:K3+N3)

      double precision KP1P3(3)
      double precision KP2P3(3)
      double precision qPC1
      double precision qPC2
      double precision S1
      double precision S2
      double precision Smat1(3)
      double precision Smat2(3)
      double precision SQmat1(3)
      double precision SQmat2(3)

      double precision tempvar1
      double precision tempvar2
      double precision tempvar3
      double precision C3
      double precision KSS
      double precision alpha
      double precision R2_S1S2
      double precision expT

      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
      double precision ans
      double precision Dgf



      gam13=gamA13+gamB13
      gam23=gamA23+gamB23
      gam12=gamA12+gamB12

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
      call Gauss_prod(A3,B3,Amat3,Bmat3,KAB3,qAB3,P3,Pmat3,Qmat3)

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

C Evaluate integral over hermite gaussians

      xgVeeg=0.0d+00
      gam12_save=gam12

      do t1=0,I1+L1
      do u1=0,J1+M1
      do v1=0,K1+N1

         do t2=0,I2+L2
         do u2=0,J2+M2
         do v2=0,K2+N2

            do t3=0,I3+L3
            do u3=0,J3+M3
            do v3=0,K3+N3

C  Apply Boys' Lemma 2 to particle 3:
C     (
C IN: |dx3 Exp(-P3(x3-Pmat3)**2) Exp(-g13(x3-x1)**2) Exp(-g23(x3-x2)**2)
C     )
C OUT: C0*Exp(-C1(x1-Pmat3)**2 Exp(-C2(x2-Pmat3)**2 Exp(-C12(x1-x2)**2
c     C1=(P3*gam13)/xsum
c     C2=(P3*gam23)/xsum
c     C12=(gam13*gam23)/xsum

      call SFBoys2(P3,gam13,gam23,Pmat3,
     x             C0,C1,C2,C12)

C  Update gamma(1,2):
      gam12=gam12_save
      gam12=gam12+C12

C  Should now have the following integral:
C  (
C  | Exp(-C1(x1-Pmat3)**2) Exp(-C2(x2-Pmat3)**2) Exp(-gam12(x1-x2)**2)
C  | Exp(-P1(x1-Pmat1)**2) Exp(-P2(x2-Pmat2)**2) 1/|x1-x2| 
C  )

C  Calculate overlap distributions.  
C  Form products of Gaussian functions:
      call Gauss_prod(P1,C1,Pmat1,Pmat3,KP1P3,qPC1,S1,Smat1,SQmat1)
      call Gauss_prod(P2,C2,Pmat2,Pmat3,KP2P3,qPC2,S2,Smat2,SQmat2)

C  Should now have the following integral:
C  (
C  | Exp(-S1(x1-Smat1)**2) Exp(-S2(x2-Smat2)**2) Exp(-gam12(x1-x2)**2)
C  | 1/|x1-x2| 
C  )

C  Evaluate integral over spherical Gaussians:
      tempvar1=(S1*S2+S1*gam12+S2*gam12)*sqrt(S1+S2)
      tempvar2=(S1*S2+S1*gam12+S2*gam12)
      tempvar3=(S1*S2+S1*gam12+S2*gam12)*(S1+S2)

      R2_S1S2=(Smat1(1)-Smat2(1))**2 +
     x        (Smat1(2)-Smat2(2))**2 +
     x        (Smat1(3)-Smat2(3))**2 

      C3=(two*sqrt(pi**5))/tempvar1

      KSS=exp(((-S1*S2*gam12)/tempvar2)*R2_S1S2)

      alpha=(S1*S1*S2*S2)/(tempvar3)

      expT=alpha*R2_S1S2

C  Organize terms for differentiation of f_KKK:
C Prepare the following KKK Gaussian functions:
C KP1P3 <=> GK1 = exp(-xi1*( f1*x1 + g1*x2 + h1*x3 )
C KP2P3 <=> GK2 = exp(-xi2*( f2*x1 + g2*x2 + h2*x3 )
C KSS   <=> GK3 = exp(-xi3*( f3*x1 + g3*x2 + h3*x3 )
      call vee_dhermite_prep(p1,c1,p2,c2,s1,s2,gam12,
     x                       xi1,xi2,xi3,
     x                       f1,g1,h1,
     x                       f2,g2,h2,
     x                       f3,g3,h3)

C  Differentiation of integral over spherical Gaussians:               
c     call gam2_vee_dhermite(t1,t2,t3,
      call g3_vee_dhermite(t1,t2,t3,
     x                     u1,u2,u3,
     x                     v1,v2,v3,
     x                     xi1,xi2,xi3,
     x                     f1,g1,h1,
     x                     f2,g2,h2,
     x                     f3,g3,h3,
     x                     expT,alpha,
     x                     Pmat1,Pmat2,Pmat3,Dgf)

      ans=C0*C3*
     x   Et1(t1)*Et2(t2)*Et3(t3) *
     x   Eu1(u1)*Eu2(u2)*Eu3(u3) *
     x   Ev1(v1)*Ev2(v2)*Ev3(v3) *
     x   Dgf

      xgVeeg=xgVeeg+ans
 


            end do
            end do
            end do

         end do
         end do
         end do

      end do
      end do
      end do
CCCCC---CHOP
C  Now have:
C KAB1 Exp(-P1(x1-Pmat1)**2) 
C KAB2 Exp(-P2(x2-Pmat2)**2) 
C KAB3 Exp(-P3(x3-Pmat3)**2)


c     write(*,*)
c     write(*,*)
c     write(*,*)
c     write(*,*)'KAB1=',KAB1
c     write(*,*)'KAB2=',KAB2
c     write(*,*)'KAB3=',KAB3
c     write(*,*)'P1=',P1
c     write(*,*)'P2=',P2
c     write(*,*)'P3=',P3
c     write(*,*)'C0=',C0
c     write(*,*)'C1=',C1
c     write(*,*)'C2=',C2
c     write(*,*)'C12=',C12
c     write(*,*)'gam12=',gam12
c     write(*,*)'gam13=',gam13
c     write(*,*)'gam23=',gam23
c     write(*,*)'ZNUC=',ZNUC
c     write(*,*)'ansg1=',ansg1
c     write(*,*)'xgVxCg=',xgVxCg



      return
      end

C=======================================================================
      subroutine vee_dhermite_prep(p1,c1,p2,c2,s1,s2,gam12,
     x                             xi1,xi2,xi3,
     x                             f1,g1,h1,
     x                             f2,g2,h2,
     x                             f3,g3,h3)
c    x                             Pmat1,Pmat2,Pmat3,
c    x                             Rmat)

C Prepare the following KKK Gaussian functions:
C KP1P3 <=> GK1 = exp(-xi1*( f1*x1 + g1*x2 + h1*x3 )
C KP2P3 <=> GK2 = exp(-xi2*( f2*x1 + g2*x2 + h2*x3 )
C KSS   <=> GK3 = exp(-xi3*( f3*x1 + g3*x2 + h3*x3 )
C=======================================================================
      implicit none
C Input variables
      double precision p1,c1
      double precision p2,c2
      double precision s1
      double precision s2
      double precision gam12
c     double precision Pmat1(3),Pmat2(3),Pmat3(3)
C Variables returned
      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
c     double precision Rmat(3)
C Local Variables

C     Gaussian K...: GK...
C     GK1 <=> KP1P3:
      xi1=(p1*c1)/(p1+c1)
      f1=1.0d+00
      g1=0.0d+00
      h1=-1.0d+00
      
C     GK2 <=> KP2P3:
      xi2=(p2*c2)/(p2+c2)
      f2=0.0d+00
      g2=1.0d+00
      h2=-1.0d+00
      
C     GK3 <=> KSS:
      xi3=(s1*s2*gam12)/(s1*s2+s1*gam12+s2*gam12)
      f3=p1/(p1+c1)
      g3=-p2/(p2+c2)
      h3=(c1/(p1+c1))-(c2/(p2+c2))
     
C Rmat needed for Boys func derivatives
c     double precision Pmat1(3),Pmat2(3),Pmat3(3)
c     double precision Rmat(3)
c     Rmat(1)=(f3*Pmat1(1))+(g3*Pmat2(1))+(h3*Pmat3(1)) 
c     Rmat(2)=(f3*Pmat1(2))+(g3*Pmat2(2))+(h3*Pmat3(2)) 
c     Rmat(3)=(f3*Pmat1(3))+(g3*Pmat2(3))+(h3*Pmat3(3)) 


      return
      end

C=======================================================================
c     subroutine gam2_vee_dhermite(t1,t2,t3,
      subroutine g3_vee_dhermite(t1,t2,t3,
     x                           u1,u2,u3,
     x                           v1,v2,v3,
     x                           xi1,xi2,xi3,
     x                           f1,g1,h1,
     x                           f2,g2,h2,
     x                           f3,g3,h3,
     x                           expT,alpha,
     x                           Pmat1,Pmat2,Pmat3,Dgf)

C=======================================================================
      implicit none

C Input variables
      integer t1,t2,t3
      integer u1,u2,u3
      integer v1,v2,v3

      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
      double precision expT
      double precision alpha
      double precision Pmat1(3),Pmat2(3),Pmat3(3) 

C Variables returned
      double precision Dgf

C Local Variables
      integer xb1,xb2,xb3
      integer yb1,yb2,yb3
      integer zb1,zb2,zb3

      double precision bin_tx1 
      double precision bin_tx2
      double precision bin_tx3
      double precision bin_uy1
      double precision bin_uy2
      double precision bin_uy3
      double precision bin_vz1
      double precision bin_vz2
      double precision bin_vz3
      double precision dfKKK_x 
      double precision dfKKK_y 
      double precision dfKKK_z 
      double precision dfBoys
      double precision ans 
      

      Dgf=0.0d+00

      do xb1=0,t1
         call cbinom(t1,xb1,bin_tx1)
      do xb2=0,t2
         call cbinom(t2,xb2,bin_tx2)
      do xb3=0,t3
         call cbinom(t3,xb3,bin_tx3)

c        call DKG123(t1-xb1,t2-xb2,t3-xb3,
         call G3xggs_DGK123(t1-xb1,t2-xb2,t3-xb3,
     x                      xi1,xi2,xi3,
     x                      f1,g1,h1,
     x                      f2,g2,h2,
     x                      f3,g3,h3,
     x                      Pmat1(1),Pmat2(1),Pmat3(1),dfKKK_x)

         do yb1=0,u1
            call cbinom(u1,yb1,bin_uy1)
         do yb2=0,u2
            call cbinom(u2,yb2,bin_uy2)
         do yb3=0,u3
            call cbinom(u3,yb3,bin_uy3)

c           call DKG123(u1-yb1,u2-yb2,u3-yb3,
            call G3xggs_DGK123(u1-yb1,u2-yb2,u3-yb3,
     x                         xi1,xi2,xi3,
     x                         f1,g1,h1,
     x                         f2,g2,h2,
     x                         f3,g3,h3,
     x                         Pmat1(2),Pmat2(2),Pmat3(2),dfKKK_y)


            do zb1=0,v1
               call cbinom(v1,zb1,bin_vz1)
            do zb2=0,v2
               call cbinom(v2,zb2,bin_vz2)
            do zb3=0,v3
               call cbinom(v3,zb3,bin_vz3)

c              call DKG123(v1-zb1,v2-zb2,v3-zb3,
               call G3xggs_DGK123(v1-zb1,v2-zb2,v3-zb3,
     x                            xi1,xi2,xi3,
     x                            f1,g1,h1,
     x                            f2,g2,h2,
     x                            f3,g3,h3,
     x                            Pmat1(3),Pmat2(3),Pmat3(3),dfKKK_z)

c     subroutine gam2_RTUV(xb1,yb1,zb1,
c    x                     xb2,yb2,zb2,
c    x                     xb3,yb3,zb3,
c    x                     expT,alpha,f,g,h,
c    x                     Pmat1,Pmat2,Pmat3,dF0)
C
c              call gam2_RTUV(xb1,yb1,zb1,
               call g3_RTUV(xb1,yb1,zb1,
     x                      xb2,yb2,zb2,
     x                      xb3,yb3,zb3,
     x                      expT,alpha,f3,g3,h3,
     x                      Pmat1,Pmat2,Pmat3,dfBoys)
 

               ans=bin_tx1*bin_tx2*bin_tx3
     x            *bin_uy1*bin_uy2*bin_uy3
     x            *bin_vz1*bin_vz2*bin_vz3
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


      return
      end

C=======================================================================
      subroutine DKG123(xb1,xb2,xb3,
     x                  xi1,xi2,xi3,
     x                  f1,g1,h1,
     x                  f2,g2,h2,
     x                  f3,g3,h3,
     x                  x1,x2,x3,dG1G2G3)
C
C Calculate derivative P_x1,x2,x3 (GK1 GK2 GK3)
C 
C = (d/P1x)^xb1 (d/dP2x)^xb2 (d/dP3x)^xb3 GK1 GK2 GK3
C
C=======================================================================
      implicit none

C Input variables
      integer xb1,xb2,xb3

      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
      double precision x1,x2,x3

C Variables returned
      double precision dG1G2G3

C Local variables
      integer i1,i2,i3
      integer j1,j2,j3

      double precision bin_i1
      double precision bin_i2
      double precision bin_i3
      double precision bin_j1
      double precision bin_j2
      double precision bin_j3
      double precision dG1,dG2,dG3
      double precision ans



      dG1G2G3=0.0d+00

      do i1=0,xb1
         call cbinom(xb1,i1,bin_i1)

         do i2=0,xb2
            call cbinom(xb2,i2,bin_i2)

            do i3=0,xb3
               call cbinom(xb3,i3,bin_i3)
               call DKGx(xb1-i1,xb2-i2,xb3-i3,
     x                   xi1,f1,g1,h1,x1,x2,x3,dG1)


               do j1=0,i1
                  call cbinom(i1,j1,bin_j1)

                  do j2=0,i2
                     call cbinom(i2,j2,bin_j2)

                     do j3=0,i3
                        call cbinom(i3,j3,bin_j3)

                        call DKGx(i1-j1,i2-j2,i3-j3,
     x                            xi2,f2,g2,h2,x1,x2,x3,dG2)

                        call DKGx(j1,j2,j3,
     x                            xi3,f3,g3,h3,x1,x2,x3,dG3)

                        ans=bin_i1*bin_i2*bin_i3
     x                     *bin_j1*bin_j2*bin_j3
     x                     *dG1*dG2*dG3

                        dG1G2G3=dG1G2G3+ans

                     end do
                  end do
               end do

            end do
         end do
      end do



      return
      end


C=======================================================================
c     subroutine gam2_RTUV(xb1,yb1,zb1,
      subroutine g3_RTUV(xb1,yb1,zb1,
     x                   xb2,yb2,zb2,
     x                   xb3,yb3,zb3,
     x                   expT,alpha,f,g,h,
     x                   Pmat1,Pmat2,Pmat3,dfBoys)
C
C
C=======================================================================
      implicit none

C Input variables
      integer xb1,yb1,zb1 
      integer xb2,yb2,zb2 
      integer xb3,yb3,zb3 

      double precision expT
      double precision alpha
      double precision f,g,h
      double precision Pmat1(3),Pmat2(3),Pmat3(3)

C Variables returned
      double precision dfBoys

C Local variables
      integer NRx,NRy,NRz
      integer N,L,M
      integer nLim 

      double precision c1
      double precision ans
      double precision Rmat(3)


C Rmat needed for Boys func derivatives
      Rmat(1)=(f*Pmat1(1))+(g*Pmat2(1))+(h*Pmat3(1)) 
      Rmat(2)=(f*Pmat1(2))+(g*Pmat2(2))+(h*Pmat3(2)) 
      Rmat(3)=(f*Pmat1(3))+(g*Pmat2(3))+(h*Pmat3(3)) 

      NRx=xb1+yb1+zb1
      NRy=xb2+yb2+zb2
      NRz=xb3+yb3+zb3

      N=xb1+xb2+xb3
      L=yb1+yb2+yb3
      M=zb1+zb2+zb3

      c1=(f**NRx)*(g**NRy)*(h**NRz)

c     nLim=NRx+NRy+NRz
      nLim=N+L+M

      call RTUV(N,L,M,nLim,
     x          expT,alpha,Rmat(1),Rmat(2),Rmat(3),ans)

      dfBoys=c1*ans

C Rmat needed for Boys func derivatives
c     Rmat(1)=(f*Pmat1(1))+(g*Pmat2(1))+(h*Pmat3(1)) 
c     Rmat(2)=(f*Pmat1(2))+(g*Pmat2(2))+(h*Pmat3(2)) 
c     Rmat(3)=(f*Pmat1(3))+(g*Pmat2(3))+(h*Pmat3(3)) 

c original
c     NRx=xb1+yb1+zb1
c     NRy=xb2+yb2+zb2
c     NRz=xb3+yb3+zb3

c     c1=(f**NRx)*(g**NRy)*(h**NRz)

c     nLim=NRx+NRy+NRz

c     call RTUV(NRx,NRy,NRz,nLim,
c    x          expT,alpha,Rmat(1),Rmat(2),Rmat(3),ans)

c     dfBoys=c1*ans



      return
      end

C=======================================================================
      subroutine G3xggs_DGK123(xb1,xb2,xb3,
     x                         xi1,xi2,xi3,
     x                         f1,g1,h1,
     x                         f2,g2,h2,
     x                         f3,g3,h3,
     x                         x1,x2,x3,dGK123)
C
C Calculate derivative P_x1,x2,x3 (GK1 GK2 GK3)
C 
C = (d/P1x)^xb1 (d/dP2x)^xb2 (d/dP3x)^xb3 GK1 GK2 GK3
C
C=======================================================================
      implicit none

C Input variables
      integer xb1,xb2,xb3

      double precision xi1,xi2,xi3
      double precision f1,g1,h1
      double precision f2,g2,h2
      double precision f3,g3,h3
      double precision x1,x2,x3

C Variables returned
      double precision dGK123

C Local variables
      integer i1,i2,i3
      integer j1,j2,j3

      double precision bin_i1
      double precision bin_i2
      double precision bin_i3
      double precision bin_j1
      double precision bin_j2
      double precision bin_j3
      double precision dGK1,dGK2,dGK3
      double precision ans



      dGK123=0.0d+00

      do i1=0,xb1
         call cbinom(xb1,i1,bin_i1)

         do i2=0,xb2
            call cbinom(xb2,i2,bin_i2)

            do i3=0,xb3
               call cbinom(xb3,i3,bin_i3)
               call DKGx(xb1-i1,xb2-i2,xb3-i3,
     x                   xi1,f1,g1,h1,x1,x2,x3,dGK1)


               do j1=0,i1
                  call cbinom(i1,j1,bin_j1)

                  do j2=0,i2
                     call cbinom(i2,j2,bin_j2)

                     do j3=0,i3
                        call cbinom(i3,j3,bin_j3)

                        call DKGx(i1-j1,i2-j2,i3-j3,
     x                            xi2,f2,g2,h2,x1,x2,x3,dGK2)

                        call DKGx(j1,j2,j3,
     x                            xi3,f3,g3,h3,x1,x2,x3,dGK3)

                        ans=bin_i1*bin_i2*bin_i3
     x                     *bin_j1*bin_j2*bin_j3
     x                     *dGK1*dGK2*dGK3

                        dGK123=dGK123+ans

                     end do
                  end do
               end do

            end do
         end do
      end do



      return
      end

C=======================================================================
      subroutine DKGx(n1,n2,n3,xi,f,g,h,x1,x2,x3,ans)
C
C Calculate derivative:
C 
C   (d/dx1)^n1 (d/dx2)^n2 (d/dx3)^n3 [exp( -xi*(f*x1 + g*x2 + h*x3) ) ] 
C
C=======================================================================
      implicit none
C Input variables
      integer n1,n2,n3

      double precision xi
      double precision f,g,h
      double precision x1,x2,x3

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

      ntotal=n1+n2+n3
      u = f*x1 + g*x2 + h*x3
      c1=(f**n1)*(g**n2)*(h**n3)
      call gdel_x(ntotal,u,xi,0.0d0,dGx) 
      ans=c1*dGx



      return
      end


C======================================================================
c     subroutine cws_G3_xgVEC(I1,J1,K1,A1,Amat1,
      subroutine G3_CWS_xgVEC(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        Cmat,ZNUC,
     *                        xgVEC)

C Calculates: xgVEC    g^A(2,3) VeC(1)
C             xgVECg2  g^A(2,3) VeC(1) g^B(2,3)
C
C Uses natural decomposition between geminal and operator
C to calculate integral as a product of overlap and VeC.
C
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
      double precision Cmat(3)

      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23

      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)

C Variables returned
      double precision xgVEC

C Local variables
      double precision zero
      parameter(zero=0.0d+00)
      double precision xvec
      double precision xggs
      double precision ans
      double precision gamma1
      double precision gamma2

      gamma1=gamA23+gamA12+gamA13
      gamma2=gamB23+gamB12+gamB13
c     gamma2=zero

      xgvec=zero
      xvec=zero
      xggs=zero

C  Call to Ari routine for evaluating VeC 
      call gfvec(I1,J1,K1,A1,Amat1,
     2           L1,M1,N1,B1,Bmat1,
     3           Cmat,xvec)



c     call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
      call G2_MD_xggs(I2,J2,K2,A2,Amat2,
     1                I3,J3,K3,A3,Amat3,
     2                L2,M2,N2,B2,Bmat2,
     3                L3,M3,N3,B3,Bmat3,
     4                gamma1,gamma2,xggs)

c     call pgiovlap(I2,J2,K2,A2,Amat2,
c    1              I3,J3,K3,A3,Amat3,
c    3              L2,M2,N2,B2,Bmat2,
c    2              L3,M3,N3,B3,Bmat3,
c    4              gamma1,gamma2,xggs)

c     write(*,*)'CHECK'
c     write(*,*)'K1=',K1
c     write(*,*)'K2=',K2
c     write(*,*)'K3=',K3
c     write(*,*)'N1=',N1
c     write(*,*)'N2=',N2
c     write(*,*)'N3=',N3
c     write(*,*)'xvec=',xvec
c     write(*,*)'xggs=',xggs
      xgVEC=xvec*xggs
c     xgVEC=ZNUC*xvec*xggs
c     xgVEC=-1.0d+00*xgVEC

c     write(*,*)
c     write(*,*)'>>> xvec  =',xvec
c     write(*,*)'>>> xggs  =',xggs
c     write(*,*)'>>> xgvec =',xgvec
c     write(*,*)


      return
      end

C======================================================================
      subroutine G3_MD_TEg013(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        xgTEg)

C Evaluates Kinetic energy for the following integrals:
C xgTE    gA(2,3) T(1)
C xgTEg1  gA(1,3) T(1) gB(2,3)
C xgTEg3  gA(2,3) T(1) gB(2,3)
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

      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23

      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)

C Variables returned
      double precision xgTEg

C Local Variables
      integer LX1
      integer MX1
      integer NX1
      double precision ans
      double precision Sx1,Sy1,Sz1
      double precision Sx2,Sy2,Sz2
      double precision Sx3,Sy3,Sz3
      double precision Sx,Sy,Sz
      double precision TSx,TSy,TSz
      double precision zero,one,two,four
      parameter(zero=0.0d+00,one=1.0d+00,two=2.0d+00,four=4.0d+00)


C-------------------------- Sx ---------------------------
      call G3_MD_xggs1D(I1,a1,Amat1(1),
     x                  I2,a2,Amat2(1),
     x                  I3,a3,Amat3(1),
     x                  L1,b1,Bmat1(1),
     x                  L2,b2,Bmat2(1),
     x                  L3,b3,Bmat3(1),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sx)
      
C-------------------------- Sy ---------------------------
      call G3_MD_xggs1D(J1,a1,Amat1(2),
     x                  J2,a2,Amat2(2),
     x                  J3,a3,Amat3(2),
     x                  M1,b1,Bmat1(2),
     x                  M2,b2,Bmat2(2),
     x                  M3,b3,Bmat3(2),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sy)
      
C-------------------------- Sz ---------------------------
      call G3_MD_xggs1D(K1,a1,Amat1(3),
     x                  K2,a2,Amat2(3),
     x                  K3,a3,Amat3(3),
     x                  N1,b1,Bmat1(3),
     x                  N2,b2,Bmat2(3),
     x                  N3,b3,Bmat3(3),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sz)
      
C-------------------------- Tx ---------------------------
      Sx1=zero
      if(L1.ge.2) then
         LX1=L1-2
         call G3_MD_xggs1D(I1,a1,Amat1(1),
     x                     I2,a2,Amat2(1),
     x                     I3,a3,Amat3(1),
     x                     LX1,b1,Bmat1(1),
     x                     L2,b2,Bmat2(1),
     x                     L3,b3,Bmat3(1),
     x                     gamA12,gamA13,gamA23,
     x                     gamB12,gamB13,gamB23,
     x                     Sx1)
      end if
      
      LX1=L1
      call G3_MD_xggs1D(I1,a1,Amat1(1),
     x                  I2,a2,Amat2(1),
     x                  I3,a3,Amat3(1),
     x                  LX1,b1,Bmat1(1),
     x                  L2,b2,Bmat2(1),
     x                  L3,b3,Bmat3(1),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sx2)
      
      LX1=L1+2
      call G3_MD_xggs1D(I1,a1,Amat1(1),
     x                  I2,a2,Amat2(1),
     x                  I3,a3,Amat3(1),
     x                  LX1,b1,Bmat1(1),
     x                  L2,b2,Bmat2(1),
     x                  L3,b3,Bmat3(1),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sx3)
      
C-------------------------- Ty ---------------------------
      Sy1=zero
      if(M1.ge.2) then
         MX1=M1-2
         call G3_MD_xggs1D(J1,a1,Amat1(2),
     x                     J2,a2,Amat2(2),
     x                     J3,a3,Amat3(2),
     x                     MX1,b1,Bmat1(2),
     x                     M2,b2,Bmat2(2),
     x                     M3,b3,Bmat3(2),
     x                     gamA12,gamA13,gamA23,
     x                     gamB12,gamB13,gamB23,
     x                     Sy1)
      end if
      
      MX1=M1
      call G3_MD_xggs1D(J1,a1,Amat1(2),
     x                  J2,a2,Amat2(2),
     x                  J3,a3,Amat3(2),
     x                  MX1,b1,Bmat1(2),
     x                  M2,b2,Bmat2(2),
     x                  M3,b3,Bmat3(2),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sy2)
      
      MX1=M1+2
      call G3_MD_xggs1D(J1,a1,Amat1(2),
     x                  J2,a2,Amat2(2),
     x                  J3,a3,Amat3(2),
     x                  MX1,b1,Bmat1(2),
     x                  M2,b2,Bmat2(2),
     x                  M3,b3,Bmat3(2),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sy3)
      
C-------------------------- Tz ---------------------------
      Sz1=zero
      if(N1.ge.2) then
         NX1=N1-2
         call G3_MD_xggs1D(K1,a1,Amat1(3),
     x                     K2,a2,Amat2(3),
     x                     K3,a3,Amat3(3),
     x                     NX1,b1,Bmat1(3),
     x                     N2,b2,Bmat2(3),
     x                     N3,b3,Bmat3(3),
     x                     gamA12,gamA13,gamA23,
     x                     gamB12,gamB13,gamB23,
     x                     Sz1)
      end if
      
      NX1=N1
      call G3_MD_xggs1D(K1,a1,Amat1(3),
     x                  K2,a2,Amat2(3),
     x                  K3,a3,Amat3(3),
     x                  NX1,b1,Bmat1(3),
     x                  N2,b2,Bmat2(3),
     x                  N3,b3,Bmat3(3),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sz2)
      
      NX1=N1+2
      call G3_MD_xggs1D(K1,a1,Amat1(3),
     x                  K2,a2,Amat2(3),
     x                  K3,a3,Amat3(3),
     x                  NX1,b1,Bmat1(3),
     x                  N2,b2,Bmat2(3),
     x                  N3,b3,Bmat3(3),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sz3)
      
      

c     TSx=dble(L1)*(dble(L1)-one)*Sx1
c    x   +(-two*b1)*(two*dble(L1)+one)*Sx2
c    x   +(-two*b1)**2 *Sx3

c     TSy=dble(M1)*(dble(M1)-one)*Sy1
c    x   +(-two*b1)*(two*dble(M1)+one)*Sy2
c    x   +(-two*b1)**2 *Sy3

c     TSz=dble(N1)*(dble(N1)-one)*Sz1
c    x   +(-two*b1)*(two*dble(N1)+one)*Sz2
c    x   +(-two*b1)**2 *Sz3

CCCC
      TSx=dble(L1*(L1)-1)*Sx1
     x   +(-two*b1)*dble(2*L1+1)*Sx2
     x   + four*b1*b1*Sx3

      TSy=dble(M1*(M1)-1)*Sy1
     x   +(-two*b1)*dble(2*M1+1)*Sy2
     x   + four*b1*b1*Sy3

      TSz=dble(N1*(N1)-1)*Sz1
     x   +(-two*b1)*dble(2*N1+1)*Sz2
     x   + four*b1*b1*Sz3


      xgTEg=(TSx * Sy * Sz
     x       +Sx *TSy * Sz
     x       +Sx * Sy *TSz)

      xgTEg=-0.5d+00*xgTEg


      return
      end

C======================================================================
      subroutine G3_MD_KE(I1,J1,K1,A1,Amat1,
     *                    I2,J2,K2,A2,Amat2,
     *                    I3,J3,K3,A3,Amat3,
     *                    L1,M1,N1,B1,Bmat1,
     *                    L2,M2,N2,B2,Bmat2,
     *                    L3,M3,N3,B3,Bmat3,
     *                    gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23,
     *                    xmass,xgTxg)

C Evaluates Kinetic energy for the following integrals:
C xgTE    gA(2,3) T(1)
C xgTEg1  gA(1,3) T(1) gB(2,3)
C xgTEg3  gA(2,3) T(1) gB(2,3)
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

      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23

      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)

      double precision xmass

C Variables returned
      double precision xgTxg

C Local Variables
      double precision Sx,Sy,Sz
      double precision Tx,Ty,Tz


C-------------------------- Sx ---------------------------
      call G3_MD_xggs1D(I1,a1,Amat1(1),
     x                  I2,a2,Amat2(1),
     x                  I3,a3,Amat3(1),
     x                  L1,b1,Bmat1(1),
     x                  L2,b2,Bmat2(1),
     x                  L3,b3,Bmat3(1),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sx)
      
C-------------------------- Sy ---------------------------
      call G3_MD_xggs1D(J1,a1,Amat1(2),
     x                  J2,a2,Amat2(2),
     x                  J3,a3,Amat3(2),
     x                  M1,b1,Bmat1(2),
     x                  M2,b2,Bmat2(2),
     x                  M3,b3,Bmat3(2),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sy)
      
C-------------------------- Sz ---------------------------
      call G3_MD_xggs1D(K1,a1,Amat1(3),
     x                  K2,a2,Amat2(3),
     x                  K3,a3,Amat3(3),
     x                  N1,b1,Bmat1(3),
     x                  N2,b2,Bmat2(3),
     x                  N3,b3,Bmat3(3),
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  Sz)
      
C-------------------------- Tx ---------------------------
      call G3_MD_T1d(I1,a1,Amat1(1),
     x               I2,a2,Amat2(1),
     x               I3,a3,Amat3(1),
     x               L1,b1,Bmat1(1),
     x               L2,b2,Bmat2(1),
     x               L3,b3,Bmat3(1),
     x               gamA12,gamA13,gamA23,
     x               gamB12,gamB13,gamB23,
     x               Tx)

C-------------------------- Ty ---------------------------
      call G3_MD_T1d(J1,a1,Amat1(2),
     x               J2,a2,Amat2(2),
     x               J3,a3,Amat3(2),
     x               M1,b1,Bmat1(2),
     x               M2,b2,Bmat2(2),
     x               M3,b3,Bmat3(2),
     x               gamA12,gamA13,gamA23,
     x               gamB12,gamB13,gamB23,
     x               Ty)

C-------------------------- Tz ---------------------------
      call G3_MD_T1d(K1,a1,Amat1(3),
     x               K2,a2,Amat2(3),
     x               K3,a3,Amat3(3),
     x               N1,b1,Bmat1(3),
     x               N2,b2,Bmat2(3),
     x               N3,b3,Bmat3(3),
     x               gamA12,gamA13,gamA23,
     x               gamB12,gamB13,gamB23,
     x               Tz)

      xgTxg=( Tx * Sy * Sz
     x       +Sx * Ty * Sz
     x       +Sx * Sy * Tz)

      xgTxg=-0.5d+00*xgTxg/xmass



      return
      end

C======================================================================
      subroutine G3_MD_T1d(I1,a1,Amat1,
     x                     I2,a2,Amat2,
     x                     I3,a3,Amat3,
     x                     L1,b1,Bmat1,
     x                     L2,b2,Bmat2,
     x                     L3,b3,Bmat3,
     x                     gamA12,gamA13,gamA23,
     x                     gamB12,gamB13,gamB23,
     x                     ans)

C======================================================================
      implicit none

C     Input Variables
      integer I1,I2,I3
      integer L1,L2,L3
      double precision A1,A2,A3
      double precision B1,B2,B3
      double precision Amat1,Amat2,Amat3
      double precision Bmat1,Bmat2,Bmat3
      double precision gamA12,gamA13,gamA23
      double precision gamB12,gamB13,gamB23

C     Variables Returned
      double precision ans

C     Local Variables
      integer ii
      integer XL1,XL3
      double precision x(10)
      double precision XB1B3
      double precision zero,one,two,three,four,eight
      parameter(zero=0.0d+00)
      parameter(one=1.0d+00)
      parameter(two=2.0d+00)
      parameter(three=3.0d+00)
      parameter(four=4.0d+00)
      parameter(eight=8.0d+00)


      do ii=1,10
         x(ii)=zero
      end do

      XB1B3=Bmat1-Bmat3

      if(L1.gt.0) then
         if(L1.gt.1) then

            XL1=L1-2
            call G3_MD_xggs1D(I1,a1,Amat1,
     x                        I2,a2,Amat2,
     x                        I3,a3,Amat3,
     x                        XL1,b1,Bmat1,
     x                        L2,b2,Bmat2,
     x                        L3,b3,Bmat3,
     x                        gamA12,gamA13,gamA23,
     x                        gamB12,gamB13,gamB23,
     x                        x(1))

         end if

         XL1=L1-1
         call G3_MD_xggs1D(I1,a1,Amat1,
     x                     I2,a2,Amat2,
     x                     I3,a3,Amat3,
     x                     XL1,b1,Bmat1,
     x                     L2,b2,Bmat2,
     x                     L3,b3,Bmat3,
     x                     gamA12,gamA13,gamA23,
     x                     gamB12,gamB13,gamB23,
     x                     x(2))

         XL1=L1-1
         XL3=L3+1
         call G3_MD_xggs1D(I1,a1,Amat1,
     x                     I2,a2,Amat2,
     x                     I3,a3,Amat3,
     x                     XL1,b1,Bmat1,
     x                     L2,b2,Bmat2,
     x                     XL3,b3,Bmat3,
     x                     gamA12,gamA13,gamA23,
     x                     gamB12,gamB13,gamB23,
     x                     x(3))
      end if

      XL1=L1
      call G3_MD_xggs1D(I1,a1,Amat1,
     x                  I2,a2,Amat2,
     x                  I3,a3,Amat3,
     x                  XL1,b1,Bmat1,
     x                  L2,b2,Bmat2,
     x                  L3,b3,Bmat3,
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  x(4))

      XL1=L1
      XL3=L3+1
      call G3_MD_xggs1D(I1,a1,Amat1,
     x                  I2,a2,Amat2,
     x                  I3,a3,Amat3,
     x                  XL1,b1,Bmat1,
     x                  L2,b2,Bmat2,
     x                  XL3,b3,Bmat3,
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  x(5))

      XL1=L1
      XL3=L3+2
      call G3_MD_xggs1D(I1,a1,Amat1,
     x                  I2,a2,Amat2,
     x                  I3,a3,Amat3,
     x                  XL1,b1,Bmat1,
     x                  L2,b2,Bmat2,
     x                  XL3,b3,Bmat3,
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  x(6))

      XL1=L1+1
      call G3_MD_xggs1D(I1,a1,Amat1,
     x                  I2,a2,Amat2,
     x                  I3,a3,Amat3,
     x                  XL1,b1,Bmat1,
     x                  L2,b2,Bmat2,
     x                  L3,b3,Bmat3,
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  x(7))

      XL1=L1+1
      XL3=L3+1
      call G3_MD_xggs1D(I1,a1,Amat1,
     x                  I2,a2,Amat2,
     x                  I3,a3,Amat3,
     x                  XL1,b1,Bmat1,
     x                  L2,b2,Bmat2,
     x                  XL3,b3,Bmat3,
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  x(8))

      XL1=L1+2
      call G3_MD_xggs1D(I1,a1,Amat1,
     x                  I2,a2,Amat2,
     x                  I3,a3,Amat3,
     x                  XL1,b1,Bmat1,
     x                  L2,b2,Bmat2,
     x                  L3,b3,Bmat3,
     x                  gamA12,gamA13,gamA23,
     x                  gamB12,gamB13,gamB23,
     x                  x(9))


c     ans = dble(L1*(L1-1))                          * x(1)
c    x    - four*gamB13*XB1B3*dble(L1)               * x(2)
c    x    + four*gamB13*dble(L1)                     * x(3)
c    x    - two*b1*(two*dble(L1)+one)                * x(4)
c    x    - four*gamB13*dble(L1)                     * x(4)
c    x    + four*gamB13*XB1B3*XB1B3                  * x(4)
c    x    - two*gamB13                               * x(4)
c    x    - eight*gamB13*XB1B3                       * x(5)
c    x    + four*gamB13                              * x(6)
c    x    + eight*gamB13*XB1B3*(b1+one)              * x(7)
c    x    - eight*gamB13*(b1+one)                    * x(8)
c    x    + (four*b1*b1+eight*b1*gamB13+four*gamB13) * x(9)

      ans = dble(L1*(L1-1))                          * x(1)
     x    - four*gamB13*XB1B3*dble(L1)               * x(2)
     x    + four*gamB13*dble(L1)                     * x(3)
     x    - two*b1*(two*dble(L1)+one)                * x(4)
     x    - four*gamB13*dble(L1)                     * x(4)
     x    + four*gamB13*gamB13*XB1B3*XB1B3           * x(4)
     x    - two*gamB13                               * x(4)
     x    - eight*gamB13*gamB13*XB1B3                * x(5)
     x    + four*gamB13*gamB13                       * x(6)
     x    + eight*b1*gamB13*XB1B3                    * x(7)
     x    + eight*gamB13*gamB13*XB1B3                * x(7)
     x    - eight*b1*gamB13                          * x(8)
     x    - eight*gamB13*gamB13                      * x(8)
     x    + (four*b1*b1+eight*b1*gamB13+four*gamB13*gamB13) * x(9)




      return
      end

