!======================================================================
      subroutine RXCHF_FAE_GAMee(nebf,ngee,GAM_ee,DAE,FAE,E_AE_GAMee)
 
!======================================================================
      implicit none

! Input Variables
      integer nebf
      integer ngee
      double precision GAM_ee(ngee)
      double precision DAE(nebf,nebf)

! Variables Returned
      double precision FAE(nebf,nebf)
      double precision E_AE_GAMee
      
! Local Variables
      double precision zero,two,half
      parameter(zero=0.0d+00,two=2.0d+00,half=5.0d-01)
!     integer nebfLT
      integer ia_12,ia_21
      integer ie1,je1,ie2,je2
!     double precision vee1
!     double precision vee2
      double precision val
      double precision XFAE(nebf,nebf)
      double precision vee1,vee2,valee


      XFAE=zero
      E_AE_GAMee=zero

      do ie2=1,nebf
      do je2=1,nebf

         do ie1=1,nebf
         do je1=1,nebf

          call pack_4D(nebf,nebf,nebf,je2,ie2,je1,ie1,ia_12)
          call pack_4D(nebf,nebf,nebf,je1,ie2,je2,ie1,ia_21)

          vee1=GAM_ee(ia_12)
          vee2=GAM_ee(ia_21)
          valee=(vee1-half*vee2)

          XFAE(ie1,je1)=XFAE(ie1,je1)+DAE(ie2,je2)*valee

          E_AE_GAMee=E_AE_GAMee+half*DAE(ie1,je1)*DAE(ie2,je2)*valee
       

         end do
         end do

      end do
      end do

!     write(*,*)'IN: FAE_GAMee  '
!     write(*,*)'E_AE_GAMee= ',E_AE_GAMee
!     write(*,*)

!     nebfLT=(nebf+nebf*nebf)/2
!     write(*,*)'IN E_from_GAM_ee: EE Contribution to Elec FOCK Matrix:'
!     call print_my_fock(nebf,nebfLT,xfocke)

!  Update the full electronic Fock matrix
      call add2fock(nebf,XFAE,FAE)


      return
      end
