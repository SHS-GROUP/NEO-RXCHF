        function gammln(xx)
        double precision :: gammln,xx
        integer :: j
        double precision :: pi,ser,stp,tmp,x,y,cof(6)
        
        pi=acos(-1.0d0)
        stp=sqrt(2.0d0*pi)
        data cof /76.18009172947146d0,-86.50532032941677d0, 
     $            24.01409824083091d0,-1.231739572450155d0,
     $            0.1208650973866179d-2,-0.5395239384953d-5/

        x=xx
        y=x
        tmp=x+5.5d0
        tmp=(x+0.5d0)*log(tmp)-tmp
        ser=1.000000000190015d0
        do j=1,6
           y=y+1.0d0
           ser=ser+cof(j)/y
        enddo

        gammln=tmp+log(stp*ser/x)
        return
        end function
