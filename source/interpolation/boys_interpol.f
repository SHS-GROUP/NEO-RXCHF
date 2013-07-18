        subroutine boys_interpol(m,t,Fmx)
        use table
        implicit none
        integer :: i,m,mult
        double precision :: t,rem,tmp,Fmx
        double precision, parameter ::eps=1.0d-15,tmax=35.0d0,zero=0.0d0
        double precision, parameter ::SQRPI=0.88622692545275801365D0
        double precision, dimension(0:lmax+6) :: F


        
        if(t>Tmax .or. abs(t-Tmax)<eps)then
          F(0) = SQRPI/SQRT(t)
          do  i=1,m
             F(i) = F(i-1)*(1D0*i-0.5D0)/t
          enddo
          Fmx=F(m)

c          Fmx=dfac(2*m-1)/(2.0d0**(m+1)) * sqrt(pi/t**(2*m+1))

          large_t=large_t+1
        elseif(t<Tmax)then
          call tsearch(t,mult,rem)

          if(rem <eps)then
            Fmx=ftable(m,mult)
          else
            Fmx=ftable(m+5,mult)-rem/6.0d0*ftable(m+6,mult)
            do i=4,0,-1
               Fmx=ftable(m+i,mult)-rem/dble(i+1)*Fmx
            enddo
          endif
          small_t=small_t+1
        endif

        end subroutine
