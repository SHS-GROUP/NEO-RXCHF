
C*****************************************************************
      SUBROUTINE gfvee(I1,J1,K1,alp1,Amat1,
     1                  I2,J2,K2,alp2,Amat2,
     2                  L1,M1,N1,beta1,Bmat1,
     3                  L2,M2,N2,beta2,Bmat2,
     4                  ans)
C*****************************************************************
C
      implicit none
C Input variables
      integer NDIM
      parameter (NDIM=3)
      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2
      double precision alp1,Amat1(NDIM)
      double precision alp2,Amat2(NDIM)
      double precision beta1,Bmat1(NDIM)
      double precision beta2,Bmat2(NDIM)
C Output variables
      double precision ans
C Local variables
      double precision p1,p2,q1,q2
      double precision xKab1(NDIM),xKab2(NDIM)
      double precision Qmat1(NDIM),Pmat1(NDIM)
      double precision Qmat2(NDIM),Pmat2(NDIM)

      integer it1,iu1,iv1
      integer it2,iu2,iv2
      double precision E1,E2,xsum,herm
      double precision ex1(0:I1+L1),ey1(0:J1+M1),ez1(0:K1+N1)
      double precision ex2(0:I2+L2),ey2(0:J2+M2),ez2(0:K2+N2)

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

C xyz for electron 2
      call omega1D(I2,L2,alp2,beta2,Amat2(1),Bmat2(1),p2,q2,Pmat2(1),
     1             Qmat2(1),xKab2(1),ex2)
      call omega1D(J2,M2,alp2,beta2,Amat2(2),Bmat2(2),p2,q2,Pmat2(2),
     1             Qmat2(2),xKab2(2),ey2)
      call omega1D(K2,N2,alp2,beta2,Amat2(3),Bmat2(3),p2,q2,Pmat2(3),
     1             Qmat2(3),xKab2(3),ez2)
C------------------hermite-loop--------------------
C loop over hermite integrals
      xsum = 0.0d0
      do it1=0,I1+L1
      do iu1=0,J1+M1
      do iv1=0,K1+N1
         do it2=0,I2+L2
         do iu2=0,J2+M2
         do iv2=0,K2+N2
         call hint2(it1,iu1,iv1,it2,iu2,iv2,
     1                 p1,p2,Pmat1,Pmat2,
     2                 herm)
         E1  =ex1(it1)*ey1(iu1)*ez1(iv1)
         E2  =ex2(it2)*ey2(iu2)*ez2(iv2)
         xsum=xsum+(E1*E2*herm)
         end do
         end do
         end do
      end do
      end do
      end do
      ans = xsum

      END

C************************************************************************
      SUBROUTINE hint2(it1,iu1,iv1,it2,iu2,iv2,
     1                 p1,p2,Pmat1,Pmat2,
     2                 ans)
C************************************************************************
C Input varaibles
      implicit none
      integer NDIM
      parameter(NDIM=3)
      integer it1
      integer iu1
      integer iv1
      integer it2
      integer iu2
      integer iv2
      double precision p1
      double precision p2
      double precision Pmat1(NDIM)
      double precision Pmat2(NDIM)
C Output varaibles
      double precision ans

C Function called
      double precision binoc

C Local varaibles
      double precision PI
C     parameter (PI=3.141592650d0)
      integer NMAX
      integer LMAX
      integer MMAX
      double precision coe1
      double precision coe2
      double precision alpha
      double precision beta
      double precision Rpp
      double precision Rnlm
      double precision vpp(NDIM)
      double precision gamma
C----------- Initialize----------------
      PI   = dacos(-1.0d0)
      NMAX = it1+it2
      LMAX = iu1+iu2
      MMAX = iv1+iv2

C prefactors
      gamma = 0.0d0
      coe1 = (-1.0d0)**(it2+iu2+iv2)
      call get_coefEE(NDIM,PI,p1,p2,gamma,Pmat1,Pmat2,
     1                coe2,alpha,beta,Rpp,vpp)

      call auxR2(NMAX,LMAX,MMAX,beta,vpp(1),vpp(2),vpp(3),Rnlm)
CCWS-DEBUG
c     write(*,*)'----pgvee----'
c     write(*,*)'N=',NMAX
c     write(*,*)'L=',LMAX
c     write(*,*)'M=',MMAX
c     write(*,*)'XWP=',vpp(1)
c     write(*,*)'YWP=',vpp(2)
c     write(*,*)'ZWP=',vpp(3)
c     write(*,*)'RNLM=',rnlm
c     write(*,*)'coe1=',coe1
c     write(*,*)'coe2=',coe2
c     write(*,*)'alpha=',beta
c     write(*,*)'----pgvee----'
CCWS-DEBUG

      ans = coe1*coe2*Rnlm

      END



C************************************************************************
      SUBROUTINE get_coefEE(NDIM,PI,p1,p2,gamma,Pmat1,Pmat2,
     1                      coe2,alpha,beta,R,vec)
C************************************************************************
C get the pre-factor for e-e repulsion integrals
C Input variables
      implicit none
      integer NDIM
      double precision PI
      double precision p1
      double precision p2
      double precision gamma
      double precision Pmat1(NDIM)
      double precision Pmat2(NDIM)
C Output variables
      double precision coe2
      double precision alpha
      double precision beta
      double precision R
      double precision vec(NDIM)
C Local variables
      integer i
      double precision c1,c2,c3

C coefficient
      c1 = ( 2.0d0*sqrt(PI**5) ) 
      c2 = ((p1*p2)+(p1*gamma)+(p2*gamma)) *
     1     sqrt(p1+p2) 
      coe2 = c1/c2
C alpha
      c1 = p1*p2*gamma
      c2 = (p1*p2)+(p1*gamma)+(p2*gamma)
      alpha = c1/c2
C beta
      c1 = (p1*p2)**2
      c2 = (p1*p2)+(p1*gamma)+(p2*gamma)
      c3 = p1+p2
      beta = c1/(c2*c3)
C distance
      R = 0.0d0
      do i=1,NDIM
         vec(i) = Pmat1(i)-Pmat2(i)
         R = R+(vec(i)**2)
      end do
      R = sqrt(R)
         
      END
