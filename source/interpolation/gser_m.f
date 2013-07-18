        subroutine gser_m(gamser,a,x)
        implicit none
        integer :: itmax,n
        double precision :: a,gamser,x,eps,zero
        double precision :: ap,del,tsum,gammln,one
        parameter(itmax=100, eps=1.0d-12,zero=0.0d0,one=1.0d0)

!        gln=gammln(a)

        if(x .le.zero)then
          if (x. lt. zero) then
          write(*,*) ' GSER: X < 0',x
          stop
          endif
        gamser=zero
        return
        endif

        ap=a
        tsum=one/a
        del=tsum
        do n=1,itmax
          ap = ap + one
          del = del * x / ap
          tsum = tsum + del
          if(abs(del) .lt. abs(tsum)*eps) goto 100
        enddo

        write(*,*)' GSER : a is too large, ITMAX is too small ',a,itmax
        stop
! note that gamser is now gamma(a,x), not P(a,x) as in the orignal
  100  gamser = tsum*exp(-x + a*log(x))
       return
       end subroutine
        

      
