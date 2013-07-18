        subroutine farray(mmax,xj,F)
! returns a whole array of Fm(T)(m=0...lmax+5) for a given T (xj in
! this case)
        implicit none
        integer :: mmax
        double precision :: xj,gln,gamser,gammcf,tmp,mp5,expxj,txj
        double precision, parameter :: SQRPI=0.88622692545275801365D0
        double precision, parameter :: one=1.0d0,zero=0.0d0,eps=1.0d-15
        double precision, dimension(0:mmax) :: F
        integer :: i,j

        
        if (xj > 35.0d0 .or. abs(xj-35.0d0)<eps)then
            if(mmax .lt. zero) then
               write(*,*) 'bad arguments in gammp: mmax',mmax
               stop
            endif 

            F(0) = SQRPI/SQRT(xj)
            do  i=1,mmax
                F(i) = F(i-1)*(1D0*i-0.5D0)/xj
            enddo
        elseif (xj >zero .and. xj < 35.0d0)then
               if(mmax < zero) then
                  write(*,*) 'bad arguments in gammp: mmax',mmax
               stop
               endif 
                mp5=mmax+0.5d0
                expxj=exp(-xj)
                txj=2.0d0*xj
                if(xj .lt. mp5+one)then
                        call gser_m(gamser,mp5,xj)
                        F(mmax)=gamser/(2.0d0*xj**mp5)
                        do i=mmax-1,0,-1
                           F(i)=(F(i+1)*txj + expxj)/(2.0d0*i+one)
                        enddo

                else
                        call gcf(gammcf,mp5,xj,gln)
!       gammcf is eq. 6.2.3 in the numerical recipes (f77)

                        tmp=(one-gammcf)*exp(gln)
                        F(mmax)=tmp/(2.0d0*xj**mp5)
                        do i=mmax-1,0,-1
                           F(i)=(F(i+1)*txj + expxj)/(2.0d0*i+one)
                        enddo
                endif
        elseif(abs(xj-zero)<eps)then
              F(mmax)=one/(2.0d0*dble(mmax)+one)
              expxj=exp(-xj)
              txj=2.0d0*xj

              do i=mmax-1,0,-1
                 F(i)=(F(i+1)*txj + expxj)/(2.0d0*i+one)
              enddo
        endif

        return
        end subroutine
