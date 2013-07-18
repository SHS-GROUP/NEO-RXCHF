C23456789
C*************************************************************
      SUBROUTINE omega1D(I,J,a,b,AX,BX,p,q,PX,QX,xKab,Eijt)
C*************************************************************
C Calculate 1-D overlap density
C
      implicit none
C input variables
      integer I,J
      double precision a,b,AX,BX
C output variables
      double precision p,q
      double precision PX,QX
      double precision xKab
      double precision Eijt(0:I+J) 
C Gaussian product
      p    = a + b
      PX   = ((a*AX)+(b*BX))/(a+b)
      q    = (a*b)/(a+b)
      Qx   = Ax-Bx
      xKab = exp(-q*Qx*Qx)

C     WRITE(*,*) 'a:',a
C     WRITE(*,*) 'b:',b
C     WRITE(*,*) 'AX:',AX
C     WRITE(*,*) 'BX:',BX
C     WRITE(*,*) 'p:',p
C     WRITE(*,*) 'q:',q
C     WRITE(*,*) 'PX:',PX
C     WRITE(*,*) 'QX:',QX

      call hecod(I,J,a,b,p,q,Qx,xKab,Eijt)

      END

C*****************************************************
      SUBROUTINE hecod(IMAX,JMAX,a,b,p,q,Qx,E0,Eijt)
C*****************************************************
C Hermite Expansion Coefficients for Overlap Density (hecod)
C *       *         *                *       *
C E(0,0,0)   = E0
C E(i+1,j,k) = C1*E(i,j,k-1) + C2*E(i,j,k) + (K+1)*E(i,j,k+1)
C E(i,j+1,k) = D1*E(i,j,k-1) + D2*E(i,j,k) + (K+1)*E(i,j,k+1)
C
C Actual implimentation:
C E(i,j,k) = C1*E(i-1,j,k-1) + C2*E(i-1,j,k) + (K+1)*E(i-1,j,k+1)
C E(i,j,k) = D1*E(i,j-1,k-1) + D2*E(i,j-1,k) + (K+1)*E(i,j-1,k+1)
C
C
      implicit none
C input variables
      integer IMAX
      integer JMAX
      double precision a,b
      double precision p,q,Qx
      double precision E0
C output variables
      double precision Eijt(0:IMAX+JMAX)
C local variables
      integer i,j,k
      double precision c1,c2,c3
      double precision d1,d2,d3
      double precision E(0:IMAX,0:JMAX,0:IMAX+JMAX)
C----------initialize-------------------
      do i=0,IMAX
      do j=0,JMAX
      do k=0,IMAX+JMAX
         E(i,j,k) = 0.0d0
      end do
      end do
      end do

      c1 = 1.0d0/(2.0d0*p)
      c2 = -q*Qx/a
      c3 = 0.0d0
      d1 = 1.0d0/(2.0d0*p)
      d2 = q*Qx/b
      d3 = 0.0d0
C----------compute E(IMAX,0,K)-------------------
      E(0,0,0) = E0
      j = 0
      do i=1,IMAX
         E(i,j,0)=(c2*E(i-1,j,0))+E(i-1,0,1)
         do k=1,i-1
            c3=dble(k+1)
            E(i,j,k)=(c1*E(i-1,j,k-1)) +
     1               (c2*E(i-1,j,k)  ) +
     2                (c3*E(i-1,j,k+1)) 
         end do
         E(i,j,i)=(c1*E(i-1,j,i-1)) +
     1            (c2*E(i-1,j,i)  ) 
      end do
C----------compute E(IMAX,JMAX,K)-------------------
      do j=1,JMAX
         do i=0,IMAX
            E(i,j,0)=(d2*E(i,j-1,0)) +
     1               E(i,j-1,1) 
            do k=1,I+J-1
               d3=dble(k+1)
               E(i,j,k)=(d1*E(i,j-1,k-1)) +
     1                  (d2*E(i,j-1,k)  ) +
     2                  (d3*E(i,j-1,k+1)) 
            end do
            E(i,j,i+j)=(d1*E(i,j-1,I+J-1)) +
     1                 (d2*E(i,j-1,I+J)  ) 
         end do
      end do
C----------prepare output-------------------
      do k=0,IMAX+JMAX
         Eijt(k)=E(IMAX,JMAX,k)
      end do

      END

C************************************************************************
      SUBROUTINE gdel_Ax(N,x,alp,Ax,gd0)
C************************************************************************
C This routine computes the derivative of Gaussian func
C           d     
C gd0   = (---)^N  Exp[-alp*(x-Ax)^2]
C          dAx    
C
C using Hermite functions:
C H[N,y] = (-1^N) Exp[y^2] (d^N/dy^N) Exp[-y^2]
C
C (d^N/dy^N) Exp[-y^2] = (-1^N) H[N,y] Exp[-y^2]
C 
C y = alp*(x-Ax)
C
C (d/dAx) = (dy/dAx) (d/dy)
C
      implicit none
      integer N
      double precision x
      double precision alp
      double precision Ax
      double precision gd0
      double precision HermiteH

      double precision y
      double precision dy
      double precision xj

C get derivative wrt y
      y=sqrt(alp)*(x-Ax)
      dy=((-1.0d0)**N)*HermiteH(N,y)*exp(-y*y)
C Jacobian factor
      xj=(-sqrt(alp))**N

      gd0 = xj * dy

      END



   
C************************************************************************
      SUBROUTINE gdel_x(N,x,alp,Ax,gd0)
C************************************************************************
C This routine computes the derivative of Gaussian func
C           d     
C gd0   = (---)^N  Exp[-alp*(x-Ax)^2]
C          dx    
C
C using the relation:
C  
C   d                 d      
C (---)^N  = (-1)^N (----)^N
C  dx                 dA   
C
C
      implicit none
      integer N
      double precision x
      double precision alp
      double precision Ax
      double precision gd0
      double precision gAx

      call gdel_Ax(N,x,alp,Ax,gAx)
      gd0 = ((-1.0d0)**N)*gAx

      END
C************************************************************************
      FUNCTION HermiteH(n,x)
C************************************************************************
      implicit none
      integer n
      double precision x
      double precision HermiteH

      integer i
      double precision ans
      double precision A(n+2)

      A(1) = 2.0d0*x
      A(2) = (4.0d0*x*x)-2.0d0

      do i=2,n
         A(i+1)=(2.0d0*x*A(i))-(2.0d0*dble(i)*A(i-1))
      end do

      if( n .eq. 0) then
         ans = 1.0d0
      else
         ans = A(n)
      end if
      HermiteH = ans

      END


C************************************************************************
      SUBROUTINE auxR2(NMAX,LMAX,MMAX,alp,a,b,c,Rnlm)
C************************************************************************
C  Using eq. 4.4 to 4.8
C
      implicit none
C  I/O variables
      integer NMAX
      integer LMAX
      integer MMAX
      double precision alp
      double precision a
      double precision b
      double precision c
      double precision Rnlm
C  local variables
      integer i,j,m,l,n
      integer jp
      integer JMAX
      integer im
      integer il
      integer in
      double precision T
      double precision Fj
      double precision Rm(MMAX+1,NMAX+LMAX+MMAX+1)
      double precision Rl(LMAX+1,NMAX+LMAX+MMAX+1)
      double precision Rn(NMAX+1,NMAX+LMAX+MMAX+1)
C initialize
      JMAX = NMAX + MMAX + LMAX 
      T    = alp*((a*a)+(b*b)+(c*c))
      do jp=1,JMAX+1
         do im=1,MMAX+1
            Rm(im,jp) = 0.0d0
         end do
         do il=1,LMAX+1
            Rl(il,jp) = 0.0d0
         end do
         do in=1,NMAX+1
            Rn(in,jp) = 0.0d0
         end do
      end do
C
C------------M-LOOP-------------------      
      do im=1,MMAX+1
         m = im - 1
         if(im .ge. 3) then
            do j=1,JMAX
               Rm(im,j) = (c*Rm(im-1,j+1))+(dble(m-1)*Rm(im-2,j+1)) 
            end do
         elseif(im .eq. 2) then
            do j=1,JMAX
               Rm(im,j) = c*Rm(im-1,j+1)
            end do
         elseif(im .eq. 1) then
            do jp=1,JMAX+1
               j = jp - 1
               call iboys(j,T,Fj)
               Rm(im,jp) = ((-2.0d0*alp)**j)*Fj
            end do
         end if
      end do
C
C------------L-LOOP-------------------      
      do il=1,LMAX+1
         l = il - 1
         if(il .ge. 3) then
            do j=1,JMAX
               Rl(il,j) = (b*Rl(il-1,j+1))+(dble(l-1)*Rl(il-2,j+1)) 
            end do
         elseif(il .eq. 2) then
            do j=1,JMAX
               Rl(il,j) = b*Rl(il-1,j+1)
            end do
         elseif(il .eq. 1) then
            do jp=1,JMAX+1
               Rl(il,jp) = Rm(MMAX+1,jp)
            end do
         end if
      end do
C
C------------N-LOOP-------------------      
      do in=1,NMAX+1
         n = in - 1
         if(in .ge. 3) then
            do j=1,JMAX
               Rn(in,j) = (a*Rn(in-1,j+1))+(dble(n-1)*Rn(in-2,j+1)) 
            end do
         elseif(in .eq. 2) then
            do j=1,JMAX
               Rn(in,j) = a*Rn(in-1,j+1)
            end do
         elseif(in .eq. 1) then
            do jp=1,JMAX+1
               Rn(in,jp) = Rl(LMAX+1,jp)
            end do
         end if
      end do

      Rnlm = Rn(NMAX+1,1)


      END 

C################gammaf#################################################
      SUBROUTINE XGAMMAF(F,T,M)
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
C--------------------------------------------------+
C     COMPUTE INCOMPLETE XGAMMA FUNCTION            |
C                                                  |
C               /1  2m (-Tx^2)                     |
C     Fm(T)  =  |   x  e       dx                  |
C               /0                                 |
C                                                  |
C--------------------------------------------------+
      PARAMETER (TMAX=20, ACCY=1.D-12, SQRPI=0.88622692545275801365D0)
      DIMENSION F(0:M)
      PT5PM = 0.5D0+M
      IF (T.GT.TMAX) THEN
        F(0) = SQRPI/SQRT(T)
        DO 10 L=1,M
  10      F(L) = F(L-1)*(1D0*L-0.5D0)/T
      ELSE
        TI2   = T*2D0
        EXPTI = EXP(-T)
        TERM  = 0.5D0/PT5PM
        SUMK  = TERM
        DO 25 K=1,1000
          TERM = TERM*T/(PT5PM+K)
          SUMK = SUMK + TERM
  25      IF (TERM.LT.ACCY) GOTO 35
  35    F(M) = EXPTI*SUMK
        DO 50 L=M-1,0,-1
  50      F(L) = (F(L+1)*TI2 + EXPTI) / (2D0*L + 1D0)
      ENDIF
      RETURN
      END

C****************************************
      SUBROUTINE iboys2(m,x,Fmx)
C****************************************
C calculates F_m(x) where m is integer
C
      implicit none
      integer m
      double precision x
      double precision Fmx

      double precision F(0:M)

      call XGAMMAF(F,x,M)
      Fmx = F(M)

      END
C23456789
C**********************
C MODULE XGAMMA
C**********************
C computes incomplete gamma function
C using routines from numerical routines
C in fortran 77 
C
C********************************************
		SUBROUTINE gaminc(a,x,ans)
C********************************************
C Calculates incomplete gamma function
C
		implicit none
		double precision a
		double precision x
		double precision ans

		double precision gamser
		double precision gln

      call gser(gamser,a,x,gln)
		ans = gamser * exp(gln)
		
		END

C*****************************************************
      FUNCTION gammln(xx)
C*****************************************************
C Returns the value of ln[Gamma(xx)] for xx > 0
C
      DOUBLE PRECISION gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END

C********************************************
      SUBROUTINE gser(gamser,a,x,gln)
C********************************************
C Return P(a,x) and ln[Gamma(a)] in gamser 
C and gln resp. 
C incomplete gamma(a,x)  = P(a,x) * Gamma(a)
C
      INTEGER ITMAX
      DOUBLE PRECISION a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
      INTEGER n
      DOUBLE PRECISION ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END

C************************************************************************
C END MODULE gauss
C************************************************************************
C23456789
C****************************************
C MODULE BOYS
C****************************************
C computes boys function using downward
C recursion formula.
C Reference:
C "On evaluation of Boys function using downward recursion ralation"
C by  B. A. Mamedov
C Journal of Mathematical Chemistry, Vol. 36, No. 3, July 2004
C
C****************************************
      SUBROUTINE iboys(m,x,Fmx)
C****************************************
C calculates F_m(x) where m is integer
C
      implicit none
      integer m
      double precision x
      double precision Fmx

      integer ISTEP
      parameter (ISTEP=10)

      integer max
      double precision F0
      double precision gam
      double precision xinv
      double precision xexpo
      double precision a
      double precision c1
      double precision fmat(m)

      integer i
      double precision ans
C-----------------------------------------------
C DDD
      call iboys2(m,x,Fmx)
      RETURN
C-----------------------------------------------
C initialize
      max = m + ISTEP
      a = 0.50d0
      xinv = 0.0d0
      if(x .ne. 0.0d0) xinv = 1.0d0/x
      xexpo = exp(-x)
      call gaminc(a,x,gam)
      if( x .gt. 0.0d0) then
         F0 = 0.5d0 * gam / sqrt(x)
      else
         F0 = 1.0d0
      end if
      ans = F0
C-----------------------------------------------
      if(m .gt. 0) then
        if(x .gt. 0.0d0) then
           fmat(1) = 0.50d0*xinv*(F0-xexpo)
           do i=2,m
              c1 = dble((2*i)-1)
              fmat(i) = 0.50d0*xinv*((c1*fmat(i-1))-xexpo)
           end do
           ans = fmat(m)
        else
           ans = 1.0d0/dble( (2*m)+1 )
        end if
      end if

      fmx = ans


      END

C************************************
C END MODULE GAMMA
C************************************
