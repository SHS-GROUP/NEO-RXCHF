C*****************************************************************
      SUBROUTINE gfke(I1,J1,K1,alp1,Amat1,
     2                L1,M1,N1,beta1,Bmat1,xmass1,ans)
C*****************************************************************
C
      implicit none
C Input variables
      integer NDIM
      parameter (NDIM=3)
      integer I1,J1,K1
      integer L1,M1,N1
      double precision alp1,Amat1(NDIM)
      double precision beta1,Bmat1(NDIM)
      double precision xmass1
C Output variables
      double precision ans
C Local variables
      double precision Sx,Sy,Sz
      double precision Tx,Ty,Tz
      double precision ddx,ddy,ddz

C--------------Overlaps------------------------
      call gfovlap_1D(I1,alp1,Amat1(1),
     2                 L1,beta1,Bmat1(1),
     4                 Sx)

      call gfovlap_1D(J1,alp1,Amat1(2),
     2                 M1,beta1,Bmat1(2),
     4                 Sy)

      call gfovlap_1D(K1,alp1,Amat1(3),
     2                 N1,beta1,Bmat1(3),
     4                 Sz)
C--------------Laplacian--------------------
C Get KE along x,y,z for electron 1
      call Lapx1(I1,alp1,Amat1(1),
     2           L1,beta1,Bmat1(1),
     4           ddx)

      call Lapx1(J1,alp1,Amat1(2),
     2           M1,beta1,Bmat1(2),
     4           ddy)

      call Lapx1(K1,alp1,Amat1(3),
     2           N1,beta1,Bmat1(3),
     4           ddz)

C--------------kinetic energy--------------------
      Tx = -ddx*0.50d0/xmass1
      Ty = -ddy*0.50d0/xmass1
      Tz = -ddz*0.50d0/xmass1

      ans = (Tx*Sy*Sz)
     1    + (Sx*Ty*Sz) 
     2    + (Sx*Sy*Tz) 

      
      END


C*****************************************************************
      SUBROUTINE Lapx1(I1,alp1,Amat1,
     2                 L1,beta1,Bmat1,
     4                 ans)
C*****************************************************************
C Calculate 1D Laplacian element with respect to x1:
C  <I1,I2,exp(-gamA*x12*x12)| D^2/Dx1^2 | L1,L2,exp(-gamB*x12^2)>
C
C  NOTE: The Laplacian operates on only on x1 in
C        G_L1(x1) G_L2(x2) exp(-gamB*x12*x12)
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
      double precision c1,c2
      double precision d1,d2,d3
      double precision p1,p2,p3
      double precision q1,q2
      double precision t1,t2
      double precision W1,W2,W3
C--------------Initialize-----------------------
      d1 = 0.0d0
      d2 = 0.0d0
      d3 = 0.0d0
      c1 = 0.0d0
      c2 = 0.0d0
      p1 = 0.0d0
      p2 = 0.0d0
      p3 = 0.0d0
      q1 = 0.0d0
      q2 = 0.0d0
      t1 = 0.0d0
      t2 = 0.0d0
C--------------Calculate ceoff-----------------------
C get the expansion coeff
      if(L1-1 .ge. 0) d1 = dble(L1*(L1-1))
      d2 = -2.0d0*beta1*dble((2*L1)+1)
      d3 = 4.0d0*beta1*beta1
      c1 = dble(L1)
      c2 = -2.0d0*beta1
C--------------Calculate W1-----------------------
      if(L1-2 .ge. 0) then
         call gfovlap_1D(I1,alp1,Amat1,
     2                 L1-2,beta1,Bmat1,
     4                 p1)
      end if

      call gfovlap_1D(I1,alp1,Amat1,
     2                 L1,beta1,Bmat1,
     4                 p2)

      call gfovlap_1D(I1,alp1,Amat1,
     2                 L1+2,beta1,Bmat1,
     4                 p3)

      W1=(d1*p1)+(d2*p2)+(d3*p3)

      ans = W1

      END
