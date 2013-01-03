!=======================================================================
      subroutine FAE_OMG1(nebf,npbf,ng1,DAE,DP,FAE,SAE1,
     x                    E_AE_OMG1,S_AE_OMG1)

! Calculate Alpha or Beta Fock matrices
! and contribution to total energy for 2-particle terms.
! AO contracted integrals are stored on disk.
!=======================================================================
      implicit none

! Input Variables
      integer ng1
      integer nebf
      integer npbf
      double precision DAE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FAE(nebf,nebf)
      double precision SAE1(nebf,nebf)
      double precision E_AE_OMG1
      double precision S_AE_OMG1

! Local variables
      integer ia
      integer ie1
      integer je1
      integer ip
      integer jp
      double precision zero
      parameter(zero=0.0d+00)
      double precision val_gam1
      double precision val_gam1s
      double precision XFAE(nebf,nebf)

      XFAE=zero
      SAE1=zero
      E_AE_OMG1=zero
      S_AE_OMG1=zero

      open(801,file='GAM_1.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      open(802,file='GAM_1s.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ip=1,npbf
         do jp=1,npbf

            do ie1=1,nebf
               do je1=1,nebf

                  call pack_4D(nebf,nebf,npbf,
     x                         je1,ie1,jp,ip,ia)

!                 val_gam1=GAM_1(ip,jp,ie1,je1)
!                 val_gam1=GAM_1(ia)
                  read(801,REC=ia) val_gam1
                  XFAE(ie1,je1)=XFAE(ie1,je1)+DP(ip,jp)*val_gam1

                  read(802,REC=ia) val_gam1s
                  SAE1(ie1,je1)=SAE1(ie1,je1)+DP(ip,jp)*val_gam1s

                  E_AE_OMG1=E_AE_OMG1+DP(ip,jp)*DAE(ie1,je1)*val_gam1
                  S_AE_OMG1=S_AE_OMG1+DP(ip,jp)*DAE(ie1,je1)*val_gam1s

               end do
            end do

         end do
      end do

!  Update the full Fock matrices
      call add2fock(nebf,XFAE,FAE)

      close(801)
      close(802)


      return
      end
!=======================================================================
      subroutine UFP_OMG1(nebf,npbf,ng1,DETOT,DP,FP,SP1,
     x                    E_P_OMG1,S_P_OMg1)

! Calculate QM particle Fock matrices
! and contribution to total energy for 2-particle terms.
! AO contracted integrals are stored on disk.
!=======================================================================
      implicit none

! Input Variables
      integer ng1
      integer nebf
      integer npbf
      double precision DETOT(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FP(npbf,npbf)
      double precision SP1(npbf,npbf)
      double precision E_P_OMG1
      double precision S_P_OMG1

! Local variables
      integer ia
      integer ie1
      integer je1
      integer ip
      integer jp
      double precision zero
      parameter(zero=0.0d+00)
      double precision val_gam1
      double precision val_gam1s
      double precision XFP(npbf,npbf)

      XFP=zero
      SP1=zero
      E_P_OMG1=zero
      S_P_OMG1=zero

      open(801,file='GAM_1.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      open(802,file='GAM_1s.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ip=1,npbf
         do jp=1,npbf

            do ie1=1,nebf
               do je1=1,nebf

                  call pack_4D(nebf,nebf,npbf,
     x                         je1,ie1,jp,ip,ia)

!                 val_gam1=GAM_1(ip,jp,ie1,je1)
!                 val_gam1=GAM_1(ia)
                  read(801,REC=ia) val_gam1
                  XFP(ip,jp)=XFP(ip,jp)+DETOT(ie1,je1)*val_gam1

                  read(802,REC=ia) val_gam1s
                  SP1(ip,jp)=SP1(ip,jp)+DETOT(ie1,je1)*val_gam1s

                  E_P_OMG1=E_P_OMG1+DP(ip,jp)*DETOT(ie1,je1)*val_gam1
                  S_P_OMG1=S_P_OMG1+DP(ip,jp)*DETOT(ie1,je1)*val_gam1s

               end do
            end do

         end do
      end do

!  Update the full Fock matrix
      call add2fock(npbf,XFP,FP)

      close(801)
      close(802)


      return
      end

