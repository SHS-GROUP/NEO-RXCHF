C************************************************************************
C END MODULE auxR
C************************************************************************
C23456789
C************************************************************************
C MODULE gauss
C************************************************************************
C Routine to calculate Gaussian recursion relation 
C
C************************************************************************
      SUBROUTINE gdel(N,L,M,alp,Amat,rmat,gdnlm)
C************************************************************************
C computes derivative of the a Gaussian function with respect to 
C its center A
C           d       d       d
C gdnlm = (---)^N (---)^L (---)^M Exp[-alp*(r-A)^2]
C          dAx     dAy     dAz
C
      implicit none
      integer NDIM
      parameter(NDIM=3)
      integer N
      integer L
      integer M
      double precision alp
      double precision Amat(NDIM)
      double precision rmat(NDIM)
      double precision gdnlm
      double precision gxval

      integer IZERO
      parameter(IZERO=0)
      integer i
      integer nmat(NDIM)
      double precision gval
      double precision ans

      nmat(1) = N
      nmat(2) = L
      nmat(3) = M

      ans = 1.0d0
      do i=1,NDIM
         call gdel_Ax(nmat(i),rmat(i),alp,Amat(i),gval)
         ans = ans * gval
      end do
      gdnlm = ans

      END
C end of subroutine hint1




C************************************************************************
      FUNCTION binoc(n,k)
C************************************************************************
      implicit none
      integer n
      integer k
      double precision binoc

      integer i
      integer ia,ib,ic
   
      
      ia = 1
      do i=1,n
         ia = ia * i
      end do

      ib = 1
      do i=1,k
         ib = ib * i
      end do

      ic = 1
      do i=1,n-k
         ic = ic * i
      end do

      binoc = dble(ia)/dble(ib*ic)

         
      END

C end of subroutine get_coefEE

C************************************************************************
      SUBROUTINE pgi_get_coefEC(NDIM,PI,p1,p2,gamma,Pmat1,Pmat2,
     1                      Cmat,coe2,alpha,beta,
     2                      Rpp,vpp,Rcs,vcs)
C************************************************************************
C get coeff for electron-proton attraction integral
C Input variables
      implicit none
      integer NDIM
      double precision PI
      double precision p1
      double precision p2
      double precision gamma
      double precision Cmat(NDIM)
      double precision Pmat1(NDIM)
      double precision Pmat2(NDIM)
C Output variables
      double precision coe2
      double precision alpha
      double precision beta
      double precision Rpp
      double precision vpp(NDIM)
      double precision Rcs
      double precision vcs(NDIM)
C Local variables
      integer i
      double precision s,Smat(NDIM)
      double precision c1,c2,c3

C coeff
      s =((p1*p2)+(p1*gamma)+(p2*gamma))/(p2+gamma)
      c1=(sqrt(PI/(p2+gamma)))**3
      coe2 = c1*(2.0d0*PI/s)

C alpha
      alpha=(p1*p2*gamma)/((p1*p2)+(p1*gamma)+(p2*gamma))
C beta
      beta = s

C Smat
      c1=(p1*p2)+(p1*gamma)
      c2=p2*gamma
      c3=(p1*p2)+(p1*gamma)+(p2*gamma)

      do i=1,NDIM
         Smat(i)=(c1*Pmat1(i)/c3)+(c2*Pmat2(i)/c3)
      end do

C distances
      Rpp = 0.0d0
      Rcs = 0.0d0
      do i=1,NDIM
         vpp(i) = Pmat1(i)-Pmat2(i)
         vcs(i) = Cmat(i) -Smat(i)
         Rpp    = Rpp+(vpp(i)**2)
         Rcs    = Rcs+(vcs(i)**2)
      end do
      Rpp = sqrt(Rpp)
      Rcs = sqrt(Rcs)
      

      END 
C************************************************************************
      SUBROUTINE pgi_get_coefEE(NDIM,PI,p1,p2,gamma,Pmat1,Pmat2,
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

C************************************************************************
      SUBROUTINE pgi_hint0_1D(it1,it2,p1,p2,gamma,Pmat1,Pmat2,
     2                 ans)
C************************************************************************
C Input varaibles
      implicit none
      integer it1
      integer it2
      double precision p1
      double precision p2
      double precision gamma
      double precision Pmat1
      double precision Pmat2
C Output varaibles
      double precision ans
C Local varaibles
      double precision PI
      parameter (PI=3.141592650d0)
      integer i
      double precision coe1
      double precision s,q,Rpp,gd0

      coe1 = (-1.0d0)**it2

      s=(p1*p2)+(p1*gamma)+(p2*gamma) 
      q=(p1*p2*gamma)/s
      Rpp=0.0d0
      Rpp=((Pmat1-Pmat2)**2)
      Rpp =sqrt(Rpp)

      call gdel_x(it1+it2,Pmat1,q,Pmat2,gd0)
      ans = coe1*sqrt(PI*PI/s)*gd0
CCWS
c     write(*,*)'gd0=',gd0

      END
C end of subroutine hint1

C************************************************************************
      SUBROUTINE pgi_hint1_new(it1,iu1,iv1,it2,iu2,iv2,
     1                 p1,p2,gamma,Pmat1,Pmat2,
     2                 Cmat,ans)
C************************************************************************
C Hermite integration involving e-p term
C
C ans = Lambda(it1,iu1,iv1) (1/r12) Lambda(it2,iu2,iv2) Exp[-gamma*r12^2]
C
C      
C ans = Coe1 * Coe2 * Exp[-alpha*Rpp^2] * F0(beta*Rpp^2)
C
C where: 
C 
C Coe1= (-1)^(it2+iu2+iv2)
C
C                2*PI^5/2
C Coe2=------------------------------ 
C      (p1*p2 + p1*gamma + p2*gamma)
C
C
C                p1*p2*gamma
C alpha=------------------------------
C        p1*p2 + p1*gamma + p2*gamma
C
C
C        (p1*p2 + p1*gamma + p2*gamma)
C beta=------------------------------
C                (p2+gamma)
C
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
      double precision gamma
      double precision Pmat1(NDIM)
      double precision Pmat2(NDIM)
      double precision Cmat(NDIM)
C Output varaibles
      double precision ans

C Function called
      double precision binoc

C Local varaibles
      double precision PI
      parameter (PI=3.141592650d0)
      integer NMAX
      integer LMAX
      integer MMAX
      integer i,n,l,m
      double precision c1
      double precision c2
      double precision c3
      double precision coe1
      double precision coe2
      double precision alpha
      double precision beta
      double precision Rpp
      double precision Rcs
      double precision Gnlm
      double precision Rnlm
      double precision xsum
      double precision vpp(NDIM)
      double precision vcs(NDIM)

      integer n1,n2,n3
      integer m1,m2,m3
      integer k1,k2,k3
      double precision d1,d2,d3
      double precision s1,s2,s3
      double precision eta1,eta2
      double precision val
C----------- Initialize----------------
      NMAX = it1+it2
      LMAX = iu1+iu2
      MMAX = iv1+iv2

CCWS
c        write(*,*)
c        write(*,*)'t1=',it1
c        write(*,*)'u1=',iu1
c        write(*,*)'v1=',iv1
c        write(*,*)'t2=',it2
c        write(*,*)'u2=',iu2
c        write(*,*)'v2=',iv2


C prefactor in Eq. 27 and Eq. 79-83 of Ref. 1
      coe1 = 1.0d0

      call pgi_get_coefEC(NDIM,PI,p1,p2,gamma,Pmat1,Pmat2,
     1                Cmat,coe2,alpha,beta,
     2                Rpp,vpp,Rcs,vcs)

      eta1=( (p1*p2)+(gamma*p1) )/( (p1*p2)+(p1*gamma)+(p2*gamma) ) 
      eta2=( gamma*p2 )/( (p1*p2)+(p1*gamma)+(p2*gamma) ) 
C----------- calc derivative----------------
C Eq. 27
      xsum = 0.0d0
      do n1=0,it1
        do n2=0,iu1
          do n3=0,iv1
            do m1=0,it2
              do m2=0,iu2
                do m3=0,iv2
                  c1=binoc(it1,n1) 
                  c2=binoc(iu1,n2) 
                  c3=binoc(iv1,n3) 
                  d1=binoc(it2,m1) 
                  d2=binoc(iu2,m2) 
                  d3=binoc(iv2,m3) 
                  s1=(-1.0d0)**(it2+iu2+iv2-m1-m2-m3)
                  s2=(-eta1)**(n1+n2+n3)
                  s3=(-eta2)**(m1+m2+m3)
                  k1=it1+it2-n1-m1
                  k2=iu1+iu2-n2-m2
                  k3=iv1+iv2-n3-m3

                  call auxR2(n1+m1,n2+m2,n3+m3,beta,
     1                      vcs(1),vcs(2),vcs(3),Rnlm)
                  call gdel(k1,k2,k3,alpha,Pmat1,Pmat2,Gnlm)
                  val=c1*c2*c3*
     1               d1*d2*d3*
     2               s1*s2*s3*
     3               Rnlm*Gnlm
                  xsum=xsum + val
CCWS
c      write(*,*)'a1=',c1
c      write(*,*)'a2=',c2
c      write(*,*)'a3=',c3
c      write(*,*)'b1=',d1
c      write(*,*)'b2=',d2
c      write(*,*)'b3=',d3
c      write(*,*)'c1=',s1
c      write(*,*)'c2=',s2
c      write(*,*)'c3=',s3
c      write(*,*)'DC1=',eta1
c      write(*,*)'DC2=',eta2
c      write(*,*)'dgxyz=',Gnlm
c      write(*,*)'dfxyz=',Rnlm

                end do
              end do
            end do
          end do
        end do
      end do
      




C     do n=0,NMAX
C       c1=binoc(NMAX,n) 
C       do l=0,LMAX
C         c2=binoc(LMAX,l) 
C         do m=0,MMAX
C           c3=binoc(MMAX,m) 
C           call auxR2(n,l,m,beta,vcs(1),vcs(2),vcs(3),Rnlm)
C           call gdel(NMAX-n,LMAX-l,MMAX-m,alpha,Pmat2,Pmat1,Gnlm)
C           xsum=xsum+(c1*c2*c3*Rnlm*Gnlm)
C         end do
C       end do
C     end do

      ans = coe1*coe2*xsum

      END

C*****************************************************************
C END MODULE primitive_integral
C*****************************************************************
C23456789
C************************************************************************
C MODULE hermite
C************************************************************************
C Routine to calculate hermite integrations  
C
C************************************************************************
      SUBROUTINE pgi_hint2(it1,iu1,iv1,it2,iu2,iv2,
     1                 p1,p2,gamma,Pmat1,Pmat2,
     2                 ans)
C************************************************************************
C Hermite integration involving e-e repulsion term
C
C ans = Lambda(it1,iu1,iv1) (1/r12) Lambda(it2,iu2,iv2) Exp[-gamma*r12^2]
C
C      
C ans = Coe1 * Coe2 * Exp[-alpha*Rpp^2] * F0(beta*Rpp^2)
C
C where: 
C 
C Coe1= (-1)^(it2+iu2+iv2)
C
C                2*PI^5/2
C Coe2=------------------------------ 
C      (p1*p2 + p1*gamma + p2*gamma)
C
C
C                p1*p2*gamma
C alpha=------------------------------
C        p1*p2 + p1*gamma + p2*gamma
C
C
C                (p1*p2)^2
C beta=------------------------------
C        (p1*p2 + p1*gamma + p2*gamma)*(p1+p2)
C
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
      double precision gamma
      double precision Pmat1(NDIM)
      double precision Pmat2(NDIM)
C Output varaibles
      double precision ans

C Function called
      double precision binoc

C Local varaibles
      double precision PI
      parameter (PI=3.141592650d0)
      integer NMAX
      integer LMAX
      integer MMAX
      integer i,n,l,m
      double precision c1
      double precision c2
      double precision c3
      double precision coe1
      double precision coe2
      double precision alpha
      double precision beta
      double precision Rpp
      double precision Gnlm
      double precision Rnlm
      double precision xsum
      double precision vpp(NDIM)
      double precision val
C----------- Initialize----------------
      NMAX = it1+it2
      LMAX = iu1+iu2
      MMAX = iv1+iv2

C prefactor in Eq. 27
      coe1 = (-1.0d0)**(it2+iu2+iv2)
      call pgi_get_coefEE(NDIM,PI,p1,p2,gamma,Pmat1,Pmat2,
     1                coe2,alpha,beta,Rpp,vpp)

C get distance
C     Rpp = 0.0d0
C     do i=1,NDIM
C        vpp(i) = Pmat1(i)-Pmat2(i)
C        Rpp = Rpp+(vpp(i)**2)
C     end do
C     Rpp = sqrt(Rpp)
C----------- calc derivative----------------
C Eq. 27
      xsum = 0.0d0
      do n=0,NMAX
        c1=binoc(NMAX,n) 
        do l=0,LMAX
          c2=binoc(LMAX,l) 
          do m=0,MMAX
            c3=binoc(MMAX,m) 
            call auxR2(n,l,m,beta,vpp(1),vpp(2),vpp(3),Rnlm)
C DDD           call XauxR2(n,l,m,beta,vpp(1),vpp(2),vpp(3),Rnlm)
            call gdel(NMAX-n,LMAX-l,MMAX-m,alpha,Pmat1,Pmat2,Gnlm)
            val = c1*c2*c3*Rnlm*Gnlm
            xsum=xsum+val
          end do
        end do
      end do

C DDD
C     WRITE(*,*) 'NMAX LMAX MMAX:',NMAX,LMAX,MMAX
C      call auxR2(NMAX,LMAX,MMAX,beta,vpp(1),vpp(2),vpp(3),xsum)
C     WRITE(*,*) 'MY_RNLM:',xsum
C     call XauxR2(NMAX,LMAX,MMAX,beta,vpp(1),vpp(2),vpp(3),xsum)
C     WRITE(*,*) 'RNLM:',xsum
C     WRITE(*,*) 

      ans = coe1*coe2*xsum


      END


C*****************************************************************
      SUBROUTINE pgike1(I1,J1,K1,alp1,Amat1,
     1                  I2,J2,K2,alp2,Amat2,
     2                  L1,M1,N1,beta1,Bmat1,
     3                  L2,M2,N2,beta2,Bmat2,
     4                  gamA,gamB,xmass1,ans,sval)
C*****************************************************************
C Kinetic energy of electron 1
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
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2
      double precision alp1,Amat1(NDIM)
      double precision alp2,Amat2(NDIM)
      double precision beta1,Bmat1(NDIM)
      double precision beta2,Bmat2(NDIM)
      double precision gamA,gamB
      double precision xmass1
C Output variables
      double precision ans
      double precision sval
C Local variables
      double precision Sx,Sy,Sz
      double precision Tx,Ty,Tz
      double precision ddx,ddy,ddz

C--------------Overlaps------------------------
      call pgiovlap_1D(I1,alp1,Amat1(1),
     1                 I2,alp2,Amat2(1),
     2                 L1,beta1,Bmat1(1),
     3                 L2,beta2,Bmat2(1),
     4                 gamA,gamB,Sx)

      call pgiovlap_1D(J1,alp1,Amat1(2),
     1                 J2,alp2,Amat2(2),
     2                 M1,beta1,Bmat1(2),
     3                 M2,beta2,Bmat2(2),
     4                 gamA,gamB,Sy)

      call pgiovlap_1D(K1,alp1,Amat1(3),
     1                 K2,alp2,Amat2(3),
     2                 N1,beta1,Bmat1(3),
     3                 N2,beta2,Bmat2(3),
     4                 gamA,gamB,Sz)
C--------------Laplacian--------------------
C Get KE along x,y,z for electron 1
      call pgi_Lapx1(I1,alp1,Amat1(1),
     1           I2,alp2,Amat2(1),
     2           L1,beta1,Bmat1(1),
     3           L2,beta2,Bmat2(1),
     4           gamA,gamB,ddx)

      call pgi_Lapx1(J1,alp1,Amat1(2),
     1           J2,alp2,Amat2(2),
     2           M1,beta1,Bmat1(2),
     3           M2,beta2,Bmat2(2),
     4           gamA,gamB,ddy)

      call pgi_Lapx1(K1,alp1,Amat1(3),
     1           K2,alp2,Amat2(3),
     2           N1,beta1,Bmat1(3),
     3           N2,beta2,Bmat2(3),
     4           gamA,gamB,ddz)

C--------------kinetic energy--------------------
      Tx = -ddx*0.50d0/xmass1
      Ty = -ddy*0.50d0/xmass1
      Tz = -ddz*0.50d0/xmass1

      ans = (Tx*Sy*Sz)
     1    + (Sx*Ty*Sz) 
     2    + (Sx*Sy*Tz) 
C--------------sval--------------------
      sval = Sx * Sy * Sz

      
      END


C*****************************************************************
      SUBROUTINE pgi_Lapx1(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1,beta1,Bmat1,
     3                 L2,beta2,Bmat2,
     4                 gamA,gamB,ans)
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
      integer I2
      integer L1
      integer L2
      double precision alp1,Amat1
      double precision alp2,Amat2
      double precision beta1,Bmat1
      double precision beta2,Bmat2
      double precision gamA,gamB
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
         call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1-2,beta1,Bmat1,
     3                 L2,beta2,Bmat2,
     4                 gamA,gamB,p1)
      end if

      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1,beta1,Bmat1,
     3                 L2,beta2,Bmat2,
     4                 gamA,gamB,p2)

      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1+2,beta1,Bmat1,
     3                 L2,beta2,Bmat2,
     4                 gamA,gamB,p3)

      W1=(d1*p1)+(d2*p2)+(d3*p3)
C
C--------------Calculate W2-----------------------
      if(L1-1 .ge. 0) then
         call pgi_calcS2(I1,alp1,Amat1,
     1            I2,alp2,Amat2,
     2            L1-1,beta1,Bmat1,
     3            L2,beta2,Bmat2,
     4            gamA,gamB,q1)
      end if

      call pgi_calcS2(I1,alp1,Amat1,
     1            I2,alp2,Amat2,
     2            L1+1,beta1,Bmat1,
     3            L2,beta2,Bmat2,
     4            gamA,gamB,q2)

      W2=-2.0d0*gamB*( (c1*q1)+(c2*q2) )
C
C--------------Calculate W3-----------------------
      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1,beta1,Bmat1,
     3                 L2,beta2,Bmat2,
     4                 gamA,gamB,t1)

      call pgi_calcS3(I1,alp1,Amat1,
     1            I2,alp2,Amat2,
     2            L1,beta1,Bmat1,
     3            L2,beta2,Bmat2,
     4            gamA,gamB,t2)

      W3=(-2.0d0*gamB*t1)+(4.0d0*gamB*gamB*t2)
C
C--------------Calculate KE-----------------------
      ans = W1+(2.0d0*W2)+W3

      END

C*****************************************************************
      SUBROUTINE pgi_calcS2(I1,alp1,Amat1,
     1                  I2,alp2,Amat2,
     2                  L1,beta1,Bmat1,
     3                  L2,beta2,Bmat2,
     4                  gamA,gamB,ans)
C*****************************************************************
C Calculate the following three overlap integrals:
C  F   = Exp[-gamB*x12]
C  S1  = <G_I1(x1) G_I2(x2) Exp[-gamA*x12] | G_L1(x1) G_L2(x2) F>
C  S2  = <G_I1(x1) G_I2(x2) Exp[-gamA*x12] | G_L1(x1) G_L2(x2) F'>
C  S3  = <G_I1(x1) G_I2(x2) Exp[-gamA*x12] | G_L1(x1) G_L2(x2) F''>
C
C  F'  = -2*gamB*x12*F
C  F'' = -2*gamB*F + (2*gamB*x12)^2 * F
C
C  x12 is expanded as:
C
C  x12   = x1B - x2B + (B1-B2)
C  x12^2 = x1B^2 + x2B^2 + (B1-B2)^2
C        + 2*x1B*(B1-B2) - 2*x2B*(B1-B2) - 2*x1B*x2B
C 
C
      implicit none
C Input variables
      integer I1
      integer I2
      integer L1
      integer L2
      double precision alp1,Amat1
      double precision alp2,Amat2
      double precision beta1,Bmat1
      double precision beta2,Bmat2
      double precision gamA,gamB
C Output variables
      double precision ans
      double precision p1,p2,p3

      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1+1,beta1,Bmat1,
     3                 L2,beta2,Bmat2,
     4                 gamA,gamB,p1)

      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1,beta1,Bmat1,
     3                 L2+1,beta2,Bmat2,
     4                 gamA,gamB,p2)

      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1,beta1,Bmat1,
     3                 L2,beta2,Bmat2,
     4                 gamA,gamB,p3)

      ans = p1-p2+((Bmat1-Bmat2)*p3)

      END


C*****************************************************************
      SUBROUTINE pgi_calcS3(I1,alp1,Amat1,
     1                  I2,alp2,Amat2,
     2                  L1,beta1,Bmat1,
     3                  L2,beta2,Bmat2,
     4                  gamA,gamB,ans)
C*****************************************************************
C Calculate the following three overlap integrals:
C  F   = Exp[-gamB*x12]
C  S1  = <G_I1(x1) G_I2(x2) Exp[-gamA*x12] | G_L1(x1) G_L2(x2) F>
C  S2  = <G_I1(x1) G_I2(x2) Exp[-gamA*x12] | G_L1(x1) G_L2(x2) F'>
C  S3  = <G_I1(x1) G_I2(x2) Exp[-gamA*x12] | G_L1(x1) G_L2(x2) F''>
C
C  F'  = -2*gamB*x12*F
C  F'' = -2*gamB*F + (2*gamB*x12)^2 * F
C
C  x12 is expanded as:
C
C  x12   = x1B - x2B + (B1-B2)
C  x12^2 = x1B^2 + x2B^2 + (B1-B2)^2
C        + 2*x1B*(B1-B2) - 2*x2B*(B1-B2) - 2*x1B*x2B
C 
C
      implicit none
C Input variables
      integer I1
      integer I2
      integer L1
      integer L2
      double precision alp1,Amat1
      double precision alp2,Amat2
      double precision beta1,Bmat1
      double precision beta2,Bmat2
      double precision gamA,gamB
C Output variables
      double precision X12
      double precision ans
      double precision p1,p2,p3
      double precision p4,p5,p6

      X12 = Bmat1-Bmat2

      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1+2,beta1,Bmat1,
     3                 L2,beta2,Bmat2,
     4                 gamA,gamB,p1)

      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1,beta1,Bmat1,
     3                 L2+2,beta2,Bmat2,
     4                 gamA,gamB,p2)

      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1,beta1,Bmat1,
     3                 L2,beta2,Bmat2,
     4                 gamA,gamB,p3)

      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1+1,beta1,Bmat1,
     3                 L2+1,beta2,Bmat2,
     4                 gamA,gamB,p4)

      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1,beta1,Bmat1,
     3                 L2+1,beta2,Bmat2,
     4                 gamA,gamB,p5)

      call pgiovlap_1D(I1,alp1,Amat1,
     1                 I2,alp2,Amat2,
     2                 L1+1,beta1,Bmat1,
     3                 L2,beta2,Bmat2,
     4                 gamA,gamB,p6)

      ans = p1+p2+(X12*X12*p3)
     1    - (2.0d0*p4)
     2    - (2.0d0*X12*p5)
     3    + (2.0d0*X12*p6)

      END


C*****************************************************************
      SUBROUTINE pgiovlap_1D(I1,alp1,Amat1,
     1                  I2,alp2,Amat2,
     2                  L1,beta1,Bmat1,
     3                  L2,beta2,Bmat2,
     4                  gamA,gamB,ans)
C*****************************************************************
C Primitive Gaussian Integral for overlap in 1D (PGIovlap)
C
C Calculates the e-e integral using
C Gaussian type geminal (GTG)
C and Cartesian Gaussian functions
C
C  ans = <GA(r1)GA(r2)|Exp[-gamma*r12^2]|GB(r1)GB(r2)>
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
C
      implicit none
C Input variables
      integer I1
      integer I2
      integer L1
      integer L2
      double precision alp1,Amat1
      double precision alp2,Amat2
      double precision beta1,Bmat1
      double precision beta2,Bmat2
      double precision gamA,gamB
C Output variables
      double precision ans
C Local variables
      double precision gamma
      double precision p1,p2,q1,q2
      double precision xKab1,xKab2
      double precision Qmat1,Pmat1
      double precision Qmat2,Pmat2

      integer it1,iu1,iv1
      integer it2,iu2,iv2
      double precision E1,E2,xsum,herm
      double precision ex1(0:I1+L1)
      double precision ex2(0:I2+L2)

C------------------calc gamma--------------------
      gamma = gamA + gamB
C------------------overlap--------------------
C get overlap and hermite expansion coefficients
C xyz for electron 1
C
      call omega1D(I1,L1,alp1,beta1,Amat1,Bmat1,p1,q1,Pmat1,
     1             Qmat1,xKab1,ex1)

C xyz for electron 2
      call omega1D(I2,L2,alp2,beta2,Amat2,Bmat2,p2,q2,Pmat2,
     1             Qmat2,xKab2,ex2)
CCWS
c     write(*,*)
c     write(*,*)'EX1',Ex1
c     write(*,*)'EX2',Ex2

C------------------hermite-loop--------------------
C loop over hermite integrals
CCWS
c        write(*,*)'ovlap1d'
      xsum = 0.0d0
      do it1=0,I1+L1
         do it2=0,I2+L2
           call pgi_hint0_1D(it1,it2,p1,p2,gamma,Pmat1,Pmat2,
     2              herm)
         E1  =ex1(it1)
         E2  =ex2(it2)
CCWS
c        write(*,*)'E1=',e1
c        write(*,*)'E2=',e2
         xsum=xsum+(E1*E2*herm)
         end do
      end do
      ans = xsum

      END
C*****************************************************************
      SUBROUTINE pgiovlap(I1,J1,K1,alp1,Amat1,
     1                  I2,J2,K2,alp2,Amat2,
     2                  L1,M1,N1,beta1,Bmat1,
     3                  L2,M2,N2,beta2,Bmat2,
     4                  gamA,gamB,sval)
C*****************************************************************
C Kinetic energy of electron 1
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
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2
      double precision alp1,Amat1(NDIM)
      double precision alp2,Amat2(NDIM)
      double precision beta1,Bmat1(NDIM)
      double precision beta2,Bmat2(NDIM)
      double precision gamA,gamB
C Output variables
      double precision sval
C Local variables
      double precision Sx,Sy,Sz
      double precision Tx,Ty,Tz
      double precision ddx,ddy,ddz

C--------------Overlaps------------------------
      call pgiovlap_1D(I1,alp1,Amat1(1),
     1                 I2,alp2,Amat2(1),
     2                 L1,beta1,Bmat1(1),
     3                 L2,beta2,Bmat2(1),
     4                 gamA,gamB,Sx)

      call pgiovlap_1D(J1,alp1,Amat1(2),
     1                 J2,alp2,Amat2(2),
     2                 M1,beta1,Bmat1(2),
     3                 M2,beta2,Bmat2(2),
     4                 gamA,gamB,Sy)

      call pgiovlap_1D(K1,alp1,Amat1(3),
     1                 K2,alp2,Amat2(3),
     2                 N1,beta1,Bmat1(3),
     3                 N2,beta2,Bmat2(3),
     4                 gamA,gamB,Sz)
C--------------sval--------------------
      sval = Sx * Sy * Sz

      END


C*****************************************************************
      SUBROUTINE pgivec(I1,J1,K1,alp1,Amat1,
     1                  I2,J2,K2,alp2,Amat2,
     2                  L1,M1,N1,beta1,Bmat1,
     3                  L2,M2,N2,beta2,Bmat2,
     4                  gamA,gamB,Cmat,ans)
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
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2
      double precision alp1,Amat1(NDIM)
      double precision alp2,Amat2(NDIM)
      double precision beta1,Bmat1(NDIM)
      double precision beta2,Bmat2(NDIM)
      double precision gamA,gamB
      double precision Cmat(NDIM)
C Output variables
      double precision ans
C Local variables
      double precision gamma
      double precision p1,p2,q1,q2
      double precision xKab1(NDIM),xKab2(NDIM)
      double precision Qmat1(NDIM),Pmat1(NDIM)
      double precision Qmat2(NDIM),Pmat2(NDIM)

      integer it1,iu1,iv1
      integer it2,iu2,iv2
      double precision E1,E2,xsum,herm
      double precision ex1(0:I1+L1),ey1(0:J1+M1),ez1(0:K1+N1)
      double precision ex2(0:I2+L2),ey2(0:J2+M2),ez2(0:K2+N2)

C------------------calc gamma--------------------
      gamma = gamA + gamB
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
CCWS
c     write(*,*)
c     write(*,*)'gam=',gamma
c     write(*,*)'I1+L1=',I1+L1
c     write(*,*)'J1+M1=',J1+M1
c     write(*,*)'K1+N1=',K1+N1
c     write(*,*)'I2+L2=',I2+L2
c     write(*,*)'J2+M2=',J2+M2
c     write(*,*)'K2+N2=',K2+N2

      do it1=0,I1+L1
      do iu1=0,J1+M1
      do iv1=0,K1+N1
         do it2=0,I2+L2
         do iu2=0,J2+M2
         do iv2=0,K2+N2
C DDD    call hint1(it1,iu1,iv1,it2,iu2,iv2,
         call pgi_hint1_new(it1,iu1,iv1,it2,iu2,iv2,
     1                 p1,p2,gamma,Pmat1,Pmat2,Cmat,
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



C*****************************************************************
      SUBROUTINE pgivee(I1,J1,K1,alp1,Amat1,
     1                  I2,J2,K2,alp2,Amat2,
     2                  L1,M1,N1,beta1,Bmat1,
     3                  L2,M2,N2,beta2,Bmat2,
     4                  gamA,gamB,ans)
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
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2
      double precision alp1,Amat1(NDIM)
      double precision alp2,Amat2(NDIM)
      double precision beta1,Bmat1(NDIM)
      double precision beta2,Bmat2(NDIM)
      double precision gamA,gamB
C Output variables
      double precision ans
C Local variables
      double precision gamma
      double precision p1,p2,q1,q2
      double precision xKab1(NDIM),xKab2(NDIM)
      double precision Qmat1(NDIM),Pmat1(NDIM)
      double precision Qmat2(NDIM),Pmat2(NDIM)

      integer it1,iu1,iv1
      integer it2,iu2,iv2
      double precision E1,E2,xsum,herm
      double precision ex1(0:I1+L1),ey1(0:J1+M1),ez1(0:K1+N1)
      double precision ex2(0:I2+L2),ey2(0:J2+M2),ez2(0:K2+N2)

C------------------calc gamma--------------------
      gamma = gamA + gamB
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
         call pgi_hint2(it1,iu1,iv1,it2,iu2,iv2,
     1                 p1,p2,gamma,Pmat1,Pmat2,
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
