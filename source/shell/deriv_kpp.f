        subroutine deriv_kpp(tlim,ulim,vlim,Pmat1,Pmat2,qP1P2,dKpp)
        implicit none
! tlim is max of t1+t2, where 0<=t1<=n1+n1' (shell pair 1/2= elec/nuc)
        integer,intent(in)::tlim,ulim,vlim
        double precision,intent(in) :: Pmat1(3),Pmat2(3),qP1P2
        double precision,intent(out) :: dKpp(0:tlim,0:ulim,0:vlim)
! Local vars
        integer :: i,j,k
        double precision :: dgx,dgy,dgz
        
        do i=0,tlim
             call gdel_x(i,Pmat(1),qP1P2,Pmat2(1),dgx)
          do j=0,ulim
             call gdel_x(j,Pmat(2),qP1P2,Pmat2(2),dgy)
           do k=0,vlim
             call gdel_x(k,Pmat(3),qP1P2,Pmat2(3),dgz)
             dKpp(i,j,k)=dgx*dgy*dgz
           enddo
          enddo
        enddo

        end subroutine
