C23456789
C***************************************************************
      SUBROUTINE hermite_quad(N,a,b,c,xn,wt)
C***************************************************************
      implicit none
      integer N
      double precision a,b,c
      double precision xn(N),wt(N)

      integer i
C----------------------------------------------------------------
C This subroutine converts the standard gauss-hermite weights
C and nodes to the generalized weights and nodes.
C The variables "a","b" and "c" are defined as :
C 
C Integrate[-Inf,+Inf] a*Exp[-b*(x-c)**2] f(x)
C
C The standard gauss-hermite quadrature is defined as :
C     a = 1   b = 1  c = 0
C 
C Variable List :  
C     N    : Order of quadrature
C     xn : node of generalized gauss-hermite quad
C     wt: weight of generalized gauss-hermite quad
C
C NOTE : The nodes and weights of standard gauss-hermite
C        quadrature can be obtained from the subroutine "gauher"
C----------------------------------------------------------------

      if(N .le. 0) STOP 'N should be greater than zero in hermite_quad'
C     if(a .eq. 0.0d0) STOP 'A should be greater than zero
C1   in hermite_quad'
      if(b .le. 0) STOP 'B should be greater than zero in hermite_quad'

            call gauher(xn,wt,N)
            do i=1,N
                  wt(i) = wt(i)*a/dsqrt(b)
                  xn(i)  = (xn(i)/dsqrt(b)) + c
            end do
      return
      END

C***************************************************************
      SUBROUTINE gauher(x,w,n)
C***************************************************************
      implicit none
      INTEGER n,MAXIT
      DOUBLE PRECISION w(n),x(n)
      DOUBLE PRECISION EPS,PIM4
      PARAMETER (EPS=3.0D-14,PIM4=.7511255444649425D0,MAXIT=10)
      INTEGER i,its,j,m
      DOUBLE PRECISION p1,p2,p3,pp,z,z1
      m=(n+1)/2
      do 13 i=1,m
        if(i.eq.1)then
          z=dsqrt(dble(2*n+1))-1.855750d0*dble(2*n+1)**(-0.16667)
        else if(i.eq.2)then
          z=z-1.140d0*dble(n)**0.4260d0/z
        else if (i.eq.3)then
          z=1.860d0*z-0.860d0*x(1)
        else if (i.eq.4)then
          z=1.910d0*z-0.910d0*x(2)
        else
          z=2.0d0*z-x(i-2)
        endif
        do 12 its=1,MAXIT
          p1=PIM4
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=z*dsqrt(2.0d0/dble(j))*p2-dsqrt(dble(j-1)/dble(j))*p3
11        continue
          pp=dsqrt(2.0d0*dble(n))*p2
          z1=z
          z=z1-p1/pp
          if(dabs(z-z1).le.EPS)goto 1
12      continue
        pause 'too many iterations in gauher'
1       x(i)=z
        x(n+1-i)=-z
        w(i)=2.0d0/(pp*pp)
        w(n+1-i)=w(i)
13    continue
      return
      END


