C*****************************************************************
      SUBROUTINE gfvec(I1,J1,K1,alp1,Amat1,
     2                  L1,M1,N1,beta1,Bmat1,
     4                  Cmat,ans)
C*****************************************************************
C Primitive Gaussian Integral for Vee (PGIVee)
C
C Calculates the e-e integral using
C Gaussian type geminal (GTG)
C and Cartesian Gaussian functions
C
C  ans = <GA(r1)GA(r2)|Exp[-gamma*r12^2]/r12|GB(r1)GB(r2)>
C
C
C------------Definitions----------
C GA and Gb are Cartesian Gaussian functions
C GA(1) : (x-Ax)^i
C GA(2) : (y-Ay)^j
C GA(3) : (z-Az)^k
C GA(4) : alpha
C GA(5) : Ax
C GA(6) : Ay
C GA(7) : Az
C
C GA(r1) = (x-Ax)^i (y-Ay)^j (z-Az)^k Exp[-alpha*{(x-Ax)^2+(y-Ay)^2+(z-Az)^2}]
C GB(r1) = (x-Bx)^l (y-By)^m (z-Bz)^n Exp[-beta *{(x-Bx)^2+(y-By)^2+(z-Bz)^2}]
C 
C gamA: Geminal coupling between GA(r1) and GA(r2) == Exp[-GamA*r12*r12]
C gamB: Geminal coupling between GB(r1) and GB(r2) == Exp[-GamB*r12*r12]
C
C NOTE: The sign of the integral is positive
C       To use it for Coulomb attraction operator, multiply ans by -1.0d0
C
      implicit none
C Input variables
      integer NDIM
      parameter (NDIM=3)
      integer I1,J1,K1
      integer L1,M1,N1
      double precision alp1,Amat1(NDIM)
      double precision beta1,Bmat1(NDIM)
      double precision Cmat(NDIM)
C Output variables
      double precision ans
C Local variables
      double precision p1,q1
      double precision xKab1(NDIM)
      double precision Qmat1(NDIM),Pmat1(NDIM)

      integer it1,iu1,iv1
      double precision E1,xsum,herm
      double precision ex1(0:I1+L1),ey1(0:J1+M1),ez1(0:K1+N1)
      double precision sx,sy,sz
C------------------overlap--------------------
C get overlap and hermite expansion coefficients
C xyz for electron 1
C
      call omega1D(I1,L1,alp1,beta1,Amat1(1),Bmat1(1),p1,q1,Pmat1(1),
     1             Qmat1(1),xKab1(1),ex1)
      call omega1D(J1,M1,alp1,beta1,Amat1(2),Bmat1(2),p1,q1,Pmat1(2),
     1             Qmat1(2),xKab1(2),ey1)
      call omega1D(K1,N1,alp1,beta1,Amat1(3),Bmat1(3),p1,q1,Pmat1(3),
     1             Qmat1(3),xKab1(3),ez1)

C------------------hermite-loop--------------------
C loop over hermite integrals
      xsum = 0.0d0
      do it1=0,I1+L1
      do iu1=0,J1+M1
      do iv1=0,K1+N1
         E1  =ex1(it1)*ey1(iu1)*ez1(iv1)
C        call hint0_1D(it1,p1,Pmat1(1),sx)
C        call hint0_1D(iu1,p1,Pmat1(2),sy)
C        call hint0_1D(iv1,p1,Pmat1(3),sz)
C        herm = sx*sy*sz
CCWS-DEBUG
c        write(*,*)'t1=',it1
c        write(*,*)'u1=',iu1
c        write(*,*)'v1=',iv1
CCWS-DEBUG
         call hint1(NDIM,it1,iu1,iv1,p1,Pmat1,Cmat,herm)
         xsum=xsum+(E1*herm)
      end do
      end do
      end do
      ans = xsum

CCWS-DEBUG
c     write(*,*)'End of AC int'
CCWS-DEBUG
      END


C end of subroutine hint1

C************************************************************************
      SUBROUTINE hint1(NDIM,it,iu,iv,p,Pmat,Cmat,ans)
C************************************************************************
      implicit none
C input
      integer NDIM
      integer it,iu,iv
      double precision p
      double precision Pmat(NDIM),Cmat(NDIM)
C output
      double precision ans
C local
      integer i
      double precision PI
      double precision Rnlm
      double precision PC(NDIM)

      PI = acos(-1.0d0)
      do i=1,NDIM
         PC(i) = Pmat(i)-Cmat(i)
      end do
      call auxR2(it,iu,iv,p,PC(1),PC(2),PC(3),Rnlm)
CCWS-DEBUG
c     write(*,*)'AC:  RNLM=',Rnlm
CCWS-DEBUG

      ans = 2.0d0*PI*Rnlm/p

      END
