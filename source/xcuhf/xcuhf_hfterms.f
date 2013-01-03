!=======================================================================
      subroutine FAE_ecore(nebf,nebf2,GAM_ecore,DAE,FAE,E_AE_ecore)

! Calculate Alpha or Beta Fock matrices
! and contribution to total energy for core alpha electron terms.
! AO contracted integrals are stored in memory.
!=======================================================================
      implicit none

! Input Variables
      integer nebf,nebf2
      double precision DAE(nebf,nebf)
      double precision GAM_ecore(nebf2)

! Variables Returned
      double precision FAE(nebf,nebf)
      double precision E_AE_ecore

! Local variables
      integer ie1,je1,ia
      double precision val_ecore
      double precision zero
      parameter(zero=0.0d+00)


      E_AE_ecore=zero

      do ie1=1,nebf
         do je1=1,nebf

            call pack_2D(nebf,je1,ie1,ia)
!           val_ecore=GAM_ecore(ie1,je1)
            val_ecore=GAM_ecore(ia)
            FAE(ie1,je1)=val_ecore
            E_AE_ecore = E_AE_ecore + DAE(ie1,je1)*val_ecore

         end do
      end do


      return
      end
!=======================================================================
      subroutine FP_pcore(npbf,npbf2,GAM_pcore,DP,FP,E_pcore)

! Calculate Alpha or Beta Fock matrices
! and contribution to total energy for 2-particle terms.
! AO contracted integrals are stored in memory.
!=======================================================================
      implicit none

! Input Variables
      integer npbf,npbf2
      double precision DP(npbf,npbf)
      double precision GAM_pcore(npbf2)

! Variables Returned
      double precision FP(npbf,npbf)
      double precision E_pcore

! Local variables
      integer ip,jp,ia
      double precision val_pcore
      double precision zero
      parameter(zero=0.0d+00)


      E_pcore = zero

      do ip=1,npbf
         do jp=1,npbf

            call pack_2D(npbf,jp,ip,ia)
!           val_ecore=GAM_ecore(ie1,je1)
            val_pcore=GAM_pcore(ia)
            FP(ip,jp)=val_pcore
            E_pcore=E_pcore+DP(ip,jp)*val_pcore

         end do
      end do


      return
      end
!=======================================================================
      subroutine FAE_GAMep(nebf,npbf,ng1,DAE,DP,GAM_ep,FAE,E_AE_GAMep)

! Calculate Alpha or Beta Fock matrix
! and contribution to total energy for NEO-HF mixed e-p term.
! AO contracted integrals are stored in memory.
!=======================================================================
      implicit none
! Input Variables
      integer nebf
      integer npbf
      integer ng1
      double precision GAM_ep(ng1)
      double precision DAE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FAE(nebf,nebf)
      double precision E_AE_GAMep

! Local variables
!     integer npbfLT
!     integer nebfLT
      integer ia
      integer ie1,je1,ip,jp
      double precision val_ep
      double precision XFAE(nebf,nebf)
      double precision zero
      parameter(zero=0.0d+00)


      XFAE=zero
      E_AE_GAMep=zero

      do ip=1,npbf
         do jp=1,npbf

            do ie1=1,nebf
               do je1=1,nebf

                  call pack_4D(nebf,nebf,npbf,
     x                         je1,ie1,jp,ip,ia)

                  val_ep=GAM_ep(ia)
!                 val_ep=GAM_ep(ip,jp,ie1,je1)
                  XFAE(ie1,je1)=XFAE(ie1,je1)+DP(ip,jp)*val_ep
                  E_AE_GAMep=E_AE_GAMep+DP(ip,jp)*DAE(ie1,je1)*val_ep

               end do
            end do

         end do
      end do

!     npbfLT=(npbf+npbf*npbf)/2
!     nebfLT=(nebf+nebf*nebf)/2
!     write(*,*)'IN E_from_GAM_ep: EP Contribution to Nuc FOCK Matrix:'
!     call print_my_fock(npbf,npbfLT,xfockp)
!     write(*,*)'IN E_from_GAM_ep: EP Contribution to Elec FOCK Matrix:'
!     call print_my_fock(nebf,nebfLT,xfocke)
!  Update the full Fock matrices
      call add2fock(nebf,XFAE,FAE)


      return
      end
!=======================================================================
      subroutine FP_GAMep(nebf,npbf,ng1,DETOT,DP,GAM_ep,FP,E_P_GAMep)

! Calculate QM Particle matrix
! and contribution to total energy for NEO-HF mixed e-p term.
! AO contracted integrals are stored in memory.
!=======================================================================
      implicit none

! Input Variables
      integer nebf
      integer npbf
      integer ng1
      double precision GAM_ep(ng1)
      double precision DETOT(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FP(npbf,npbf)
      double precision E_P_GAMep

! Local variables
!     integer npbfLT
!     integer nebfLT
      integer ia
      integer ie1,je1,ip,jp
      double precision val_ep
      double precision XFP(npbf,npbf)
      double precision zero
      parameter(zero=0.0d+00)


      XFP=zero
      E_P_GAMep=zero

      do ip=1,npbf
         do jp=1,npbf

            do ie1=1,nebf
               do je1=1,nebf

                  call pack_4D(nebf,nebf,npbf,
     x                         je1,ie1,jp,ip,ia)

                  val_ep=GAM_ep(ia)
!                 val_ep=GAM_ep(ip,jp,ie1,je1)
                  XFP(ip,jp)=XFP(ip,jp)+DETOT(ie1,je1)*val_ep
                  E_P_GAMep=E_P_GAMep+DP(ip,jp)*DETOT(ie1,je1)*val_ep

               end do
            end do

         end do
      end do

!     npbfLT=(npbf+npbf*npbf)/2
!     nebfLT=(nebf+nebf*nebf)/2
!     write(*,*)'IN E_from_GAM_ep: EP Contribution to Nuc FOCK Matrix:'
!     call print_my_fock(npbf,npbfLT,xfockp)
!     write(*,*)'IN E_from_GAM_ep: EP Contribution to Elec FOCK Matrix:'
!     call print_my_fock(nebf,nebfLT,xfocke)
!  Update the full Fock matrices
      call add2fock(npbf,XFP,FP)


      return
      end
!======================================================================
      subroutine FAE_GAMee(nebf,ngee,GAM_ee,DAE,DETOT,FAE,E_AE_GAMee)
 
!======================================================================
      implicit none

! Input Variables
      integer nebf
      integer ngee
      double precision GAM_ee(ngee)
      double precision DAE(nebf,nebf)
      double precision DETOT(nebf,nebf)

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
      double precision OMG12,OMG21


      XFAE=zero
      E_AE_GAMee=zero

      do ie2=1,nebf
      do je2=1,nebf

         do ie1=1,nebf
         do je1=1,nebf

          call pack_4D(nebf,nebf,nebf,je2,ie2,je1,ie1,ia_12)
          call pack_4D(nebf,nebf,nebf,je1,ie2,je2,ie1,ia_21)

          OMG12=GAM_ee(ia_12)
          OMG21=GAM_ee(ia_21)

!         write(*,*)'OMG12=',OMG12
!         write(*,*)'OMG21=',OMG21
!         write(*,*)
!----------------Form-Alpha-Fock-Matrix-------------------(

         XFAE(ie1,je1)=XFAE(ie1,je1)
     x                  + DETOT(ie2,je2)*OMG12
     x                  - DAE(ie2,je2)*OMG21

!----------------Form-Alpha-Fock-Matrix-------------------)

!----------------Contribution-to-Total-Energy--------------------(
!           val=OMG12-OMG21
            E_AE_GAMee=E_AE_GAMee+half*DAE(ie1,je1)*(
     x                    DETOT(ie2,je2)*OMG12
     x                  - DAE(ie2,je2)*OMG21   )
       
!           AE_ee=AE_ee+DAE(ie1,je1)*DAE(ie2,je2)*val
!           BE_ee=BE_ee+DBE(ie1,je1)*DBE(ie2,je2)*val
!----------------Contribution-to-Total-Energy--------------------(

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
