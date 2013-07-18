        module table
        implicit none
        integer, parameter :: lmax = 32, xmax=360
        double precision, parameter :: del = 0.1d0
!**** Now 6th order expansion (7term Taylor expansion)
! The table size is lmax+6+1 by xmax+1 (up to t=36.d0)
! Currently, table(39,361) : column major order
        double precision, dimension(0:38,0:xmax) :: ftable
c        double precision, dimension(-1:2*lmax-1) :: dfac
        integer(kind=8) :: small_t, large_t
        double precision :: pi
        end module

        subroutine calc_table
        use table
        implicit none
        integer m,jx,i
        double precision xj
        double precision, dimension(0:lmax+6) :: F

        small_t=0
        large_t=0
        do jx=0,xmax
                xj=jx*del
! farray returns a whole array of Fm(T)(m=0...lmax+6) for a given T (xj in this case)
                call farray(lmax+6,xj,F)
!                if(xj ==25.d0 .or. xj==36.0d0)then
!                print *,xj, (F(i), i=0,lmax+5)
!                end if
                ftable(:,jx)=F
        enddo
c        print *,'ftable(:,100) in calc_table',ftable(:,100) 
        end 
        
        
c        subroutine dble_fac
c        use table
c        implicit none
c        integer :: i,j
c        pi=acos(-1.0d0)

c        dfac=1.0d0
c        do i=2,2*lmax-1
c           do j=i,1,-2
c                dfac(i)=dfac(i)*dble(j)  
c           enddo
c        enddo
c        end subroutine
        
