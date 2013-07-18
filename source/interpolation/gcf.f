        subroutine gcf(gammcf,a,x,gln)
        implicit none
        integer, parameter :: itmax=100
        double precision :: a,gammcf,gln,x
        double precision, parameter :: eps=1.0d-12,fpmin=1.0d-30,one=1.0d0
        double precision, parameter :: two=2.0d0
        integer :: i
        double precision :: an,b,c,d,del,h,gammln

        gln=gammln(a)

        b = x + one -a
        c = one / fpmin
        d = one / b
        h = d

        do i=1, itmax
           an=-i*(i-a)
           b= b + two
          
           d = an*d + b
           if(abs(d) .lt. fpmin)d=fpmin
                    
           c = b + an/c
           if(abs(c) .lt. fpmin)c=fpmin

           d=one/d
           del = d*c
           h = h*del
           if(abs(del-one) .lt.eps) goto 100
        enddo

        write(*,*) 'a is too large, itmax is too small in gcf',a,itmax
  100 gammcf=exp(-x+a*log(x)-gln)*h
      return
      end subroutine

        
