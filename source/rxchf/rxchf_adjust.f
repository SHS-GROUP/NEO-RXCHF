!=======================================================================
      subroutine RXCHF_adjust_omg2_ints(ng2,GM2_1ICR,GM2_2ICR,GM2sICR)
! Hack for 2 electron (one regular electron) case in a singlet arrangement
!    Sets all GM2sICR = 0
!    Sets all GM2_2ICR = 0
!    Multiplies all GM2_1ICR by 1/2
!=======================================================================
      implicit none

      integer i,ng2
      double precision GM2_1ICR(ng2),GM2_2ICR(ng2),GM2sICR(ng2)
      double precision half,zero

      half=0.50d+00
      zero=0.0d+00

      do i=1,ng2
        GM2_1ICR(i)=half*GM2_1ICR(i)
        GM2_2ICR(i)=zero
        GM2sICR(i)=zero
      end do

      return
      end

