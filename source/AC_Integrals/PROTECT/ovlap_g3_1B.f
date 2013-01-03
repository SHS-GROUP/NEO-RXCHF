C23456789
C*****************************************************************
      SUBROUTINE XG3ovlap_typ2(I1,J1,K1,alp1,Amat1,
     1                 I2,J2,K2,alp2,Amat2,
     1                 I3,J3,K3,alp3,Amat3,
     2                 L1,M1,N1,beta1,Bmat1,
     3                 L2,M2,N2,beta2,Bmat2,
     3                 L3,M3,N3,beta3,Bmat3,
     3                 cx1,Cmat1,
     4                 gamA12,gamA13,gamA23,
     4                 gamB12,gamB13,gamB23,ans)
C*****************************************************************
C
      implicit none
C Input variables
      integer NDIM
      parameter (NDIM=3)
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      double precision alp1,Amat1(NDIM)
      double precision alp2,Amat2(NDIM)
      double precision alp3,Amat3(NDIM)
      double precision beta1,Bmat1(NDIM)
      double precision beta2,Bmat2(NDIM)
      double precision beta3,Bmat3(NDIM)
      double precision cx1,Cmat1(NDIM)
      double precision gamA12,gamA13,gamA23
      double precision gamB12,gamB13,gamB23
C Output variables
      double precision ans

C Local variables
      double precision sx,sy,sz

C----------------------------------------------------------------
      call G3ovlap_typ2_1D(I1,alp1,Amat1(1),
     1                     I2,alp2,Amat2(1),
     1                     I3,alp3,Amat3(1),
     2                     L1,beta1,Bmat1(1),
     3                     L2,beta2,Bmat2(1),
     4                     L3,beta3,Bmat3(1),
     4                     cx1,Cmat1(1),
     5                     gamA12,gamA13,gamA23,
     5                     gamB12,gamB13,gamB23,
     6                     sx)

      call G3ovlap_typ2_1D(J1,alp1,Amat1(2),
     1                     J2,alp2,Amat2(2),
     1                     J3,alp3,Amat3(2),
     2                     M1,beta1,Bmat1(2),
     3                     M2,beta2,Bmat2(2),
     4                     M3,beta3,Bmat3(2),
     4                     cx1,Cmat1(2),
     5                     gamA12,gamA13,gamA23,
     5                     gamB12,gamB13,gamB23,
     6                     sy)


      call G3ovlap_typ2_1D(K1,alp1,Amat1(3),
     1                     K2,alp2,Amat2(3),
     1                     K3,alp3,Amat3(3),
     2                     N1,beta1,Bmat1(3),
     3                     N2,beta2,Bmat2(3),
     4                     N3,beta3,Bmat3(3),
     4                     cx1,Cmat1(3),
     5                     gamA12,gamA13,gamA23,
     5                     gamB12,gamB13,gamB23,
     6                     sz)

      ans = sx*sy*sz


      END

C*****************************************************************
      SUBROUTINE G3ovlap_typ2_1D(I1,alp1,Amat1,
     1                      I2,alp2,Amat2,
     1                      I3,alp3,Amat3,
     2                      L1,beta1,Bmat1,
     3                      L2,beta2,Bmat2,
     4                      L3,beta3,Bmat3,
     4                      cx1,Cmat1,
     5                      gamA12,gamA13,gamA23,
     5                      gamB12,gamB13,gamB23,
     6                      ans)
C*****************************************************************
C Calculates the following integral
C <GA(1)GA(2)GA(3)|g(1,2)g(1,3)g(2,3)|GB(1)GB(2)GB(3)>
C
      implicit none
C Input variables
      integer I1
      integer I2
      integer I3
      integer I4
      integer L1
      integer L2
      integer L3
      integer L4
      double precision alp1,Amat1
      double precision alp2,Amat2
      double precision alp3,Amat3
      double precision alp4,Amat4
      double precision beta1,Bmat1
      double precision beta2,Bmat2
      double precision beta3,Bmat3
      double precision beta4,Bmat4
      double precision cx1,Cmat1
      double precision gamA12,gamA13,gamA23
      double precision gamB12,gamB13,gamB23
C Output variables
      double precision ans
C Functions
      double precision binoc
C Local variables
C     --quadrature--
      integer i,j
      double precision p1,Pmat1,coef1
      double precision y,z
      double precision xsum
C----------------------------------------------------------------
C Gaussian product
      call gp(beta1,cx1,Bmat1,Cmat1,coef1,p1,Pmat1)
C----------------------------------------------------------------
       xsum = 0.0d0
       do j=0,L1
            z = binoc(L1,j) * ( (Pmat1-Bmat1)**(L1-j) )
            call G3ovlap_1D(I1,alp1,Amat1,
     1                      I2,alp2,Amat2,
     1                      I3,alp3,Amat3,
     2                      j,p1,Pmat1,
     3                      L2,beta2,Bmat2,
     4                      L3,beta3,Bmat3,
     5                      gamA12,gamA13,gamA23,
     6                      gamB12,gamB13,gamB23,y)

            xsum = xsum + (y*z)
       end do
       ans = xsum * coef1

      END


