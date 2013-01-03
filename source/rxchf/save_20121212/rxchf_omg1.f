!=======================================================================
      subroutine RXCHF_FP_OMG1(nebf,npbf,ng1,DBE,DP,FP,SP,
     x                    E_P_OMG1,S_P_OMG1)

! Calculate QM particle Fock matrices
! and contribution to total energy for 2-particle terms.
! AO contracted integrals are stored on disk.
!=======================================================================
      implicit none

! Input Variables
      integer ng1
      integer nebf
      integer npbf
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FP(npbf,npbf)
      double precision SP(npbf,npbf)
      double precision E_P_OMG1
      double precision S_P_OMG1

! Local variables
      integer ia
      integer ie1
      integer je1
      integer ip
      integer jp
      integer ios
      double precision zero
      parameter(zero=0.0d+00)
      double precision val_gam1
      double precision val_gam1s
      double precision XFP(npbf,npbf)

      XFP=zero
      SP=zero
      E_P_OMG1=zero
      S_P_OMG1=zero

      open(801,file='GAM_1.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8,iostat=ios)

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
                  XFP(ip,jp)=XFP(ip,jp)+DBE(ie1,je1)*val_gam1

                  read(802,REC=ia) val_gam1s
                  SP(ip,jp)=SP(ip,jp)+DBE(ie1,je1)*val_gam1s

                  E_P_OMG1=E_P_OMG1+DP(ip,jp)*DBE(ie1,je1)*val_gam1
                  S_P_OMG1=S_P_OMG1+DP(ip,jp)*DBE(ie1,je1)*val_gam1s

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
!=======================================================================
      subroutine RXCHF_FBE_OMG1(nebf,npbf,ng1,DBE,DP,FBE,SBE,
     x                    E_BE_OMG1,S_BE_OMG1)

! Calculate special electron Fock matrices
! and contribution to total energy for 2-particle terms.
! AO contracted integrals are stored on disk.
!=======================================================================
      implicit none

! Input Variables
      integer ng1
      integer nebf
      integer npbf
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)

! Variables Returned
      double precision FBE(nebf,nebf)
      double precision SBE(nebf,nebf)
      double precision E_BE_OMG1
      double precision S_BE_OMG1

! Local variables
      integer ia
      integer ie1
      integer je1
      integer ip
      integer jp
      integer ios
      double precision zero
      parameter(zero=0.0d+00)
      double precision val_gam1
      double precision val_gam1s
      double precision XFBE(nebf,nebf)

      XFBE=zero
      SBE=zero
      E_BE_OMG1=zero
      S_BE_OMG1=zero

      open(801,file='GAM_1.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8,iostat=ios)

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
                  XFBE(ie1,je1)=XFBE(ie1,je1)+DP(ip,jp)*val_gam1

                  read(802,REC=ia) val_gam1s
                  SBE(ie1,je1)=SBE(ie1,je1)+DP(ip,jp)*val_gam1s

                  E_BE_OMG1=E_BE_OMG1+DP(ip,jp)*DBE(ie1,je1)*val_gam1
                  S_BE_OMG1=S_BE_OMG1+DP(ip,jp)*DBE(ie1,je1)*val_gam1s

               end do
            end do

         end do
      end do

!  Update the full Fock matrix
      call add2fock(nebf,XFBE,FBE)

      close(801)
      close(802)


      return
      end

