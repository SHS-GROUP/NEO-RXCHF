        program large
        use table
        implicit none
        integer :: m,n
        double precision :: t,pi
        real(kind=8) :: ans,tmp1,tmp2
        character(len=32) input
        
        pi=acos(-1.0d0)
        n=iargc()

        if(n==0)then
           print *,'usage: ./prog m t'
        elseif(n==1)then
           t=35.5d0
           call getarg(1,input)
           read(input,*)m
           print *,'m,t: ',m,t
        elseif(n==2)then
           call getarg(1,input)
           read(input,*)m
           call getarg(2,input)
           read(input,*)t
           print *,'m,t: ',m,t
        endif 

        call dble_fac
!        print *,'dfac: ',dfac

        tmp1=2.0d0**(m+1)
        tmp2=sqrt(pi/t**(2*m+1))

        print *,'1: 2.0d0**(m+1) ',tmp1
        print *,'2: sqrt(pi/t**(2*m+1)) ',tmp2
        print *,'3: dfac(2*m-1) ', dfac(2*m-1)
        print *,'Fmt :2*3/1 '

        ans=dfac(2*m-1)/tmp1
        ans=ans*tmp2

        print *,'Fmt: ', ans

        end
