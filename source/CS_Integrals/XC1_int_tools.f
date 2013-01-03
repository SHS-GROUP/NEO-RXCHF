C=======================================================================
      subroutine Gauss_prod(A,B,Amat,Bmat,KAB,qAB,P,Pmat,Qmat)

C     IN: Exp(-A(x1-Amat)**2) Exp(-B(x1-Bmat)**2)
C    OUT: KAB Exp(-P(x1-Pmat)**2)
C         qAB Qmat
C=======================================================================
      implicit none

C Input variables
      double precision A
      double precision B
      double precision Amat(3)
      double precision Bmat(3)

C Variables returned
      double precision KAB(3)
      double precision qAB
      double precision P
      double precision Pmat(3)
      double precision Qmat(3)

C Local Variables
c     double precision Q2AB


      qAB=A*B/(A+B)
      Qmat(1)=Amat(1)-Bmat(1)
      Qmat(2)=Amat(2)-Bmat(2)
      Qmat(3)=Amat(3)-Bmat(3)
c     Q2AB=Qmat(1)**2 + Qmat(2)**2 + Qmat(3)**2

c     KAB=exp(-qAB*Q2AB)
      KAB(1)=exp(-qAB*Qmat(1)*Qmat(1))
      KAB(2)=exp(-qAB*Qmat(2)*Qmat(2))
      KAB(3)=exp(-qAB*Qmat(3)*Qmat(3))

      P=A+B

      Pmat(1)=( A*Amat(1)+B*Bmat(1) ) / (A+B)
      Pmat(2)=( A*Amat(2)+B*Bmat(2) ) / (A+B)
      Pmat(3)=( A*Amat(3)+B*Bmat(3) ) / (A+B)
      

      return
      end

C=======================================================================
      subroutine Gauss_prod_1D(A,B,Amat,Bmat,KAB,qAB,P,Pmat,Qmat)

C     IN: Exp(-A(x1-Amat)**2) Exp(-B(x1-Bmat)**2)
C    OUT: KAB Exp(-P(x1-Pmat)**2)
C         qAB Qmat
C     1-dimensional version of Gauss_prod routine
C=======================================================================
      implicit none

C Input variables
      double precision A
      double precision B
      double precision Amat
      double precision Bmat

C Variables returned
      double precision KAB
      double precision qAB
      double precision P
      double precision Pmat
      double precision Qmat

C Local Variables
c     double precision Q2AB


      qAB=A*B/(A+B)
      Qmat=Amat-Bmat
c     Q2AB=Qmat**2 + Qmat**2 + Qmat**2

c     KAB=exp(-qAB*Q2AB)
      KAB=exp(-qAB*Qmat*Qmat)

      P=A+B

      Pmat=( A*Amat+B*Bmat ) / (A+B)
      

      return
      end

C
C=======================================================================
      subroutine Gauss_prod_S(A,B,Amat,Bmat,KAB,P,Pmat)

C     IN: Exp(-A(x1-Amat)**2) Exp(-B(x1-Bmat)**2)
C    OUT: KAB Exp(-P(x1-Pmat)**2)
C=======================================================================
      implicit none

C Input variables
      double precision A
      double precision B
      double precision Amat(3)
      double precision Bmat(3)

C Variables returned
      double precision KAB
      double precision P
      double precision Pmat(3)

C Local Variables
      double precision qAB
      double precision Qmat(3)
      double precision Q2AB


      qAB=A*B/(A+B)
      Qmat(1)=Amat(1)-Bmat(1)
      Qmat(2)=Amat(2)-Bmat(2)
      Qmat(3)=Amat(3)-Bmat(3)
      Q2AB=Qmat(1)**2 + Qmat(2)**2 + Qmat(3)**2

      KAB=exp(-qAB*Q2AB)

      P=A+B

      Pmat(1)=( A*Amat(1)+B*Bmat(1) ) / (A+B)
      Pmat(2)=( A*Amat(2)+B*Bmat(2) ) / (A+B)
      Pmat(3)=( A*Amat(3)+B*Bmat(3) ) / (A+B)
      

      return
      end


C=======================================================================
      subroutine ghec(iLim,jLim,KAB,a,b,p,q,Qx,EOUT)

C Generate Hermite Expansion Coef.
C Evaluates Eq. 16 & 17 in Persson & Taylor paper.
C=======================================================================
      implicit none

C Input variables
      integer iLim
      integer jLim
      double precision KAB
      double precision q
      double precision Qx
      double precision a
      double precision b
      double precision p

C Variables returned
      double precision EOUT(0:iLim+jLim) 

C Local variables
      integer i 
      integer j 
      integer t 
      double precision E(0:iLim,0:jLim,0:iLim+jLim+1) 
      double precision x1
      double precision x2
      double precision x3
      double precision zero
      parameter(zero=0.0d+00) 
      double precision two
      parameter(two=2.0d+00) 
      double precision one
      parameter(one=1.0d+00) 

      x1=one/(two*p)
      x2=-q*Qx/a
      x3=zero

      do i=0,iLim
         do j=0,jLim
            do t=0,iLim+jLim+1
               E(i,j,t)=zero
            end do
         end do
      end do

C     ==E(0,0,0)==
      E(0,0,0)=KAB

C     ==E(i,0,t)==
c     j=0
      do i=1,iLim
         E(i,0,0)=x2*E(i-1,0,0)+E(i-1,0,1)
         do t=1,i
            x3=dble(t+1)
            E(i,0,t)=x1*E(i-1,0,t-1) +
     x               x2*E(i-1,0,t)   +
     x               x3*E(i-1,0,t+1)
         end do
      end do

C     ==E(i,j,t)==
      x2=q*Qx/b
      do j=1,jLim
         do i=0,iLim
            E(i,j,0)=x2*E(i,j-1,0)+E(i,j-1,1)
            do t=1,i+j
               x3=dble(t+1)
               E(i,j,t)=x1*E(i,j-1,t-1) +
     x                  x2*E(i,j-1,t)   +
     x                  x3*E(i,j-1,t+1)
            end do
         end do
      end do
      
C  All we need:
      EOUT(0:iLim+jLim)=E(iLim,jLim,0:iLim+jLim)


C Check Expansion Coefs:
c     do i=0,iLim
c        do j=0,jLim
c           do t=0,iLim+jLim
c              write(*,*)'i=',i,'j=',j,'t=',t,'E=',E(i,j,t)
c           end do
c        end do
c     end do


      return
      end
C
C=======================================================================
      subroutine SFBoys2(P3,gam13,gam23,Pmat3,
     x                   C0,C1,C2,C12)

C     (
C IN: |dx3 Exp(-P3(x3-Pmat3)**2) Exp(-g13(x3-x1)**2) Exp(-g23(x3-x2)**2)
C     )
C OUT: C0*Exp(-C1(x1-Pmat3)**2 Exp(-C2(x2-Pmat3)**2 Exp(-C12(x1-x2)**2
C=======================================================================
      implicit none

C Input variables
      double precision P3
      double precision gam13
      double precision gam23
      double precision Pmat3(3)

C Variables returned
      double precision C0
      double precision C1
      double precision C2
      double precision C12

C Local variables
      double precision xsum
      double precision two
      double precision three
      double precision pi
      parameter(two=2.0d+00)
      parameter(three=3.0d+00)
      parameter(pi=3.14159265358979d+00)

      xsum=P3+gam13+gam23

      C0=(pi/xsum)**(three/two)

      C1=(P3*gam13)/xsum
      C2=(P3*gam23)/xsum
      C12=(gam13*gam23)/xsum

      return
      end
C
C=======================================================================
      subroutine SFBoys2_4p(P4,gam14,gam24,gam34,Pmat4,
     x                      C0,C1,C2,C3,C12,C13,C23)

C  Apply Boys' Lemma 2 to particle 4:
C     (
C IN: |dx4 Exp(-P4(x4-Pmat4)**2) Exp(-g14(x4-x1)**2) Exp(-g24(x4-x2)**2)
C     |    Exp(-g34(x4-x3)**2)
C     )
C
C OUT: C0*
C         Exp(-C1(x1-Pmat4)**2
C     *   Exp(-C2(x2-Pmat4)**2
C     *   Exp(-C3(x3-Pmat4)**2
C     *   Exp(-C12(x1-x2)**2
C     *   Exp(-C13(x1-x3)**2
C     *   Exp(-C23(x2-x3)**2
C
c     C1=(P4*gam14)/xsum     
c     C2=(P4*gam24)/xsum     
c     C3=(P4*gam34)/xsum     
c     C12=(gam14*gam24)/xsum 
c     C13=(gam14*gam34)/xsum 
c     C23=(gam24*gam34)/xsum 
C
C=======================================================================
      implicit none

C Input variables
      double precision P4
      double precision gam14
      double precision gam24
      double precision gam34
      double precision Pmat4(3)

C Variables returned
      double precision C0
      double precision C1
      double precision C2
      double precision C3
      double precision C12
      double precision C13
      double precision C23

C Local variables
      double precision xsum
      double precision two
      double precision three
      double precision pi
      parameter(two=2.0d+00)
      parameter(three=3.0d+00)
      parameter(pi=3.14159265358979d+00)

      xsum=P4+gam14+gam24+gam34

      C0=(pi/xsum)**(three/two)

      C1=(P4*gam14)/xsum     
      C2=(P4*gam24)/xsum     
      C3=(P4*gam34)/xsum     
      C12=(gam14*gam24)/xsum 
      C13=(gam14*gam34)/xsum 
      C23=(gam24*gam34)/xsum 


      return
      end

C=======================================================================
      subroutine RTUV(tLim,uLim,vLim,nLim,
     x                expt,alpha,XPC,YPC,ZPC,RTUV_OUT)
C
C=======================================================================
      implicit none

C Input variables
      integer tLim
      integer uLim
      integer vLim
      integer nLim
      double precision expt 
      double precision alpha 
      double precision XPC 
      double precision YPC 
      double precision ZPC 

C Variables Returned
      double precision RTUV_OUT

C Local Variables
      integer t,u,v,n
      double precision Fmx
      double precision zero
      parameter(zero=0.0d+00)
      double precision two
      parameter(two=2.0d+00)
c     double precision F(0:nLim)
      double precision R(0:tLim+1,0:uLim+1,0:vLim+1,0:nLim+1)

C     Initialization
      do t=0,tLim
         do u=0,uLim
            do v=0,vLim
               do n=0,nLim+1
                  R(t,u,v,n)=zero
               end do
            end do
         end do
      end do

C     ===R(0,0,0,n)===
      do n=0,nLim
         call iboys2(n,expt,Fmx)
         R(0,0,0,n)=(-two * alpha)**dble(n) * Fmx
      end do  
c     call gammaF(F,expt,nLim)
c     R(0,0,0,0:nLim)=F(0:nLim)

C     ===R(1,0,0,n)===
      t=1
      do n=0,nLim
         R(t,0,0,n)=XPC*R(t-1,0,0,n+1)
      end do  

C     ===R(t>=2,0,0,n)===
      do t=2,tLim
         do n=0,nLim
            R(t,0,0,n)=dble(t-1)*R(t-2,0,0,n+1) + XPC*R(t-1,0,0,n+1)
         end do  
      end do  

CCCCCCCCCCCCCC

C     ===R(t,u=1,0,n)===
      u=1
c     v=0
      do t=0,tLim
         do n=0,nLim
            R(t,u,0,n)=YPC*R(t,u-1,0,n+1)
         end do  
      end do  

C     ===R(t,u>=2,0,n)===
      do u=2,uLim
         do t=0,tLim
            do n=0,nLim
               R(t,u,0,n)=dble(u-1)*R(t,u-2,0,n+1) + YPC*R(t,u-1,0,n+1)
            end do  
         end do  
      end do
 
CCCCCCCCCCCCCC

C     ===R(t,u,v=1,n)===
      v=1
      do t=0,tLim
         do u=0,uLim
            do n=0,nLim
               R(t,u,v,n)=ZPC*R(t,u,v-1,n+1)
            end do  
         end do  
      end do  

C     ===R(t,u,v>=2,n)===
      do v=2,vLim
         do t=0,tLim
            do u=0,uLim
               do n=0,nLim
                R(t,u,v,n)=dble(v-1)*R(t,u,v-2,n+1) + ZPC*R(t,u,v-1,n+1)
               end do  
            end do  
         end do
      end do
 
CCCCCCCCCCCCCC

C  All we need:
      RTUV_OUT=R(tLim,uLim,vLim,0)

C Print output for testing:
c     do t=0,tLim
c        do u=0,uLim
c           do v=0,vLim
c              do n=0,nLim
c               write(*,*)'R(',t,u,v,')=',RTUV_OUT(t,u,v)
c              end do  
c           end do  
c        end do
c     end do


      return
      end

C=======================================================================
      subroutine cbinom(n,k,coef)
C
C  / n \
C  |   | = n! / ( k! (n-k)! )
C  \ k /
C=======================================================================

C Input variables
      integer n
      integer k
C Variables returned
      double precision coef
C Local variables
      integer i
      integer nk
      double precision nfac
      double precision kfac
      double precision nkfac
      double precision one
      parameter(one=1.0d+00)


C Initialization
      nfac=one
      kfac=one
      nkfac=one


      if(n.eq.0.or.n.eq.1) then
         nfac=one
      else
         do i=n,1,-1
            nfac=nfac*dble(i)
         end do
      end if


      if(k.eq.0.or.k.eq.1) then
         kfac=one
      else
         do i=k,1,-1
            kfac=kfac*dble(i)
         end do
      end if


      nk=n-k
      if(nk.eq.0.or.nk.eq.1) then
         nkfac=one
      else
         do i=nk,1,-1
            nkfac=nkfac*dble(i)
         end do
      end if

c           write(*,*)'nk=',nk
c           write(*,*)'nfac=',nfac
c           write(*,*)'kfac=',kfac
c           write(*,*)'nkfac=',nkfac

      coef= nfac / ( kfac * nkfac )



      return
      end

C=======================================================================
c     program test_cbinom

C=======================================================================
c     integer n
c     integer k
c     double precision coef
      
c     n=4
c     k=2

c     call cbinom(n,k,coef)

c     write(*,*)'coef=',coef
c     


c     end program test_cbinom

C=======================================================================
      subroutine evaluate_CHI(xb,yb,zb,a,ans)

C=======================================================================
      implicit none

C Input Variables
      integer xb,yb,zb
      double precision a

C Variables Returned
      double precision ans

C Local Variables
      integer I1,J1,K1
      integer L1,M1,N1
      double precision zero,two
      parameter(zero=0.0d+00,two=2.0d+00)
      double precision a1,b1
      double precision Cmat(3)

      ans=zero

      I1=xb
      J1=yb
      K1=zb

      L1=0
      M1=0
      N1=0

      a1=a/two
      b1=a1

      Cmat(1)=zero
      Cmat(2)=zero
      Cmat(3)=zero

      call gfovlap(I1,J1,K1,a1,Cmat,
     x             L1,M1,N1,b1,Cmat,ans)



      return
      end

C=======================================================================
      subroutine contract_BF(a,b,Amat,Bmat,p,Pmat,
     x                       KAB,XPA,YPA,ZPA,XPB,YPB,ZPB)

C=======================================================================
      implicit none

C Input Variables
      double precision Amat(3)
      double precision Bmat(3)
      double precision a,b

C Variables returned
      double precision p
      double precision Pmat(3)
      double precision KAB
      double precision XPA,YPA,ZPA
      double precision XPB,YPB,ZPB

C Local Variables
      double precision R2AB

      p=a+b

      R2AB=( Amat(1)-Bmat(1) )**2
     x    +( Amat(2)-Bmat(2) )**2
     x    +( Amat(3)-Bmat(3) )**2


      Pmat(1)=(a*Amat(1) + b*Bmat(1))/p
      Pmat(2)=(a*Amat(2) + b*Bmat(2))/p
      Pmat(3)=(a*Amat(3) + b*Bmat(3))/p

      XPA=Pmat(1)-Amat(1)
      YPA=Pmat(2)-Amat(2)
      ZPA=Pmat(3)-Amat(3)

      XPB=Pmat(1)-Bmat(1)
      YPB=Pmat(2)-Bmat(2)
      ZPB=Pmat(3)-Bmat(3)

      KAB=exp(-a*b/p * R2AB )

      return
      end


