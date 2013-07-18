        subroutine dble_fac
        use table
        implicit none
        integer :: i,j

        dfac=1.0d0
        do i=2,2*lmax-1
           do j=i,1,-2
                dfac(i)=dfac(i)*dble(j)  
           enddo
        enddo
        end subroutine
