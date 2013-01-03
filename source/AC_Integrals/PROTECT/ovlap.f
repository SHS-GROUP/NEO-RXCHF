C*****************************************************************
      SUBROUTINE gfovlap(I1,J1,K1,alp1,Amat1,
     2                   L1,M1,N1,beta1,Bmat1,ans)
C*****************************************************************
      implicit none
C Input variables
      integer NDIM
      parameter (NDIM=3)
      integer I1,J1,K1
      integer L1,M1,N1
      double precision alp1,Amat1(NDIM)
      double precision beta1,Bmat1(NDIM)
C Output variables
      double precision ans
C Local variables
      double precision Sx,Sy,Sz

      call gfovlap_1D(I1,alp1,Amat1(1),L1,beta1,Bmat1(1),Sx)
      call gfovlap_1D(J1,alp1,Amat1(2),M1,beta1,Bmat1(2),Sy)
      call gfovlap_1D(K1,alp1,Amat1(3),N1,beta1,Bmat1(3),Sz)
      ans = Sx*Sy*Sz

      END
C*****************************************************************
      SUBROUTINE gfovlap_1D(I1,alp1,Amat1,
     2                  L1,beta1,Bmat1,ans)
C*****************************************************************
C
      implicit none
C Input variables
      integer I1
      integer L1
      double precision alp1,Amat1
      double precision beta1,Bmat1
C Output variables
      double precision ans
C Local variables
      double precision p1,q1
      double precision xKab1
      double precision Qmat1,Pmat1

      integer it1,iu1,iv1
      double precision E1,xsum,herm
      double precision ex1(0:I1+L1)

      double precision PI
C     parameter (PI=3.141592650d0)

      PI = dacos(-1.0d0)
C------------------overlap--------------------
C get overlap and hermite expansion coefficients
C xyz for electron 1
C
      call omega1D(I1,L1,alp1,beta1,Amat1,Bmat1,p1,q1,Pmat1,
     1             Qmat1,xKab1,ex1)

C------------------hermite-loop--------------------
C loop over hermite integrals
C      xsum = 0.0d0
C      do it1=0,I1+L1
C         call hint0_1D(it1,p1,Pmat1,herm)
C         E1  =ex1(it1)
C         xsum=xsum+(E1*herm)
C      end do
C      ans = xsum

      ans = ex1(0)*sqrt(PI/p1)

      END

C************************************************************************
      SUBROUTINE hint0_1D(it1,p1,Pmat1,ans)
C************************************************************************
C Input varaibles
      implicit none
      integer it1
      double precision p1
      double precision Pmat1
C Output varaibles
      double precision ans
C Local varaibles
      double precision PI
C     parameter (PI=3.141592650d0)

      PI = dacos(-1.0d0)
      ans = 0.0d0
      if(it1 .eq. 0) ans = sqrt(PI/p1)

      END
