!======================================================================
      subroutine RXCHF_FOCK(LCMF,LNEOHF,nelec,NAE,NBE,
     x                    nebf,nebf2,npbf,npbf2,ngee,
     x                    ng1,ng2,ng3,ng4,
     x                    SZG2ICR,SZG3ICR,SZG4ICR,
     x                    NG2CHK,NG3CHK,NG4CHK,
     x                    DAE,DBE,DP,
     x                    GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                    GM2_1ICR,GM2_2ICR,GM2sICR,
     x                    GM3_1ICR,GM3_2ICR,
     x                    GM4ICR,
     x                    FP,FAE,FBE, 
     x                    E_OMG1, 
     x                    E_OMG2,
     x                    E_OMG3,
     x                    E_OMG4,
     x                    S_OMG1,
     x                    S_OMG2,
     x                    E_total,
     x                    S_total)

!======================================================================
      implicit none

! Input Variables
      logical LCMF,LNEOHF
      integer nebf,npbf,npbf2,nebf2,ngee,ng1,ng2,ng3,ng4
      integer SZG2ICR,SZG3ICR,SZG4ICR
      integer nelec,NAE,NBE
      integer NG2CHK,NG3CHK,NG4CHK
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GAM_ecore(nebf2)
      double precision GAM_pcore(npbf2)
      double precision GAM_ep(ng1)
      double precision GAM_ee(ngee)
!     double precision GM2ICR(ng2)
!     double precision GM2sICR(ng2)
!     double precision GM3ICR(ng3)
!     double precision GM4ICR(ng4)

      double precision GM2_1ICR(SZG2ICR)
      double precision GM2_2ICR(SZG2ICR)
      double precision GM2sICR(SZG2ICR)
      double precision GM3_1ICR(SZG3ICR)
      double precision GM3_2ICR(SZG3ICR)
      double precision GM4ICR(SZG4ICR)

! Variables Returned
      double precision FAE(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision FP(npbf,npbf)
      double precision E_UHF
      double precision E_total
      double precision S_total
      double precision E_ecore
      double precision E_pcore
      double precision E_GAMep
      double precision E_GAMee
      double precision E_OMG1
      double precision E_OMG2
      double precision E_OMG3
      double precision E_OMG4
      double precision S_OMG1
      double precision S_OMG2

! Local Variables
      double precision XFAE(nebf,nebf)
      double precision XFBE(nebf,nebf)
      double precision XFP(npbf,npbf)

      double precision XSAE(nebf,nebf)
      double precision XSBE(nebf,nebf)
      double precision XSP(npbf,npbf)

      double precision E_AE_OMG2
      double precision E_AE_OMG3
      double precision E_AE_OMG4
      double precision S_AE_OMG2

      double precision E_BE_OMG1
      double precision E_BE_OMG2
      double precision E_BE_OMG3
      double precision E_BE_OMG4
      double precision S_BE_OMG1
      double precision S_BE_OMG2

      double precision E_P_OMG1
      double precision E_P_OMG2
      double precision E_P_OMG3
      double precision E_P_OMG4
      double precision S_P_OMG1
      double precision S_P_OMG2

      double precision Psi_HOMG_Psi

      double precision zero
      parameter(zero=0.0d+00)
      double precision one
      parameter(one=1.0d+00)

      logical rxchfdbg
C ARS( testing
      integer nebflt,npbflt,i,j
      rxchfdbg=.false.
      if(LCMF) rxchfdbg=.true. 
C )
C      write(*,*) LCMF,LNEOHF,nelec,NAE,NBE,
C     x                    nebf,nebf2,npbf,npbf2,ngee,
C     x                    ng1,ng2,ng3,ng4,
C     x                    SZG2ICR,SZG3ICR,SZG4ICR,
C     x                    NG2CHK,NG3CHK,NG4CHK,
C     x                    GM2_1ICR,GM2_2ICR,GM2sICR,
C     x                    GM3_1ICR,GM3_2ICR,GM3_3ICR,
C     x                    GM4ICR

! Initialize
      FAE=zero 
      FBE=zero 
      FP=zero 

      XFAE=zero 
      XSAE=zero 
      XFP=zero 
      XSP=zero 
      XFBE=zero 
      XSBE=zero 

      E_AE_OMG2 = zero 
      E_AE_OMG3 = zero 
      E_AE_OMG4 = zero 
      S_AE_OMG2 = zero 

      E_P_OMG1 = zero 
      E_P_OMG2 = zero 
      E_P_OMG3 = zero 
      E_P_OMG4 = zero 
      S_P_OMG1 = zero 
      S_P_OMG2 = zero 

      E_BE_OMG1 = zero 
      E_BE_OMG2 = zero 
      E_BE_OMG3 = zero 
      E_BE_OMG4 = zero 
      S_BE_OMG1 = zero 
      S_BE_OMG2 = zero 

! Construct Electronic Fock Matrix (regular electrons)
      call RXCHF_IC_construct_FAE(LNEOHF,nelec,NAE,NBE,
     x                              nebf,npbf,nebf2,ngee,
     x                              ng1,ng2,ng3,ng4,
     x                              SZG2ICR,SZG3ICR,SZG4ICR,
     x                              NG2CHK,NG3CHK,NG4CHK,
     x                              DAE,DBE,DP,
     x                              GM2_1ICR,GM2_2ICR,GM2sICR,
     x                              GM3_1ICR,GM3_2ICR,
     x                              GM4ICR,
     x                              XFAE,XSAE,
     x                              E_AE_OMG2,
     x                              E_AE_OMG3,
     x                              E_AE_OMG4,
     x                              S_AE_OMG2)

! Construct QM Particle Fock Matrix
      call RXCHF_IC_construct_FP(LNEOHF,nelec,NAE,NBE,
     x                           nebf,npbf,npbf2,ngee,
     x                           ng1,ng2,ng3,ng4,
     x                           SZG2ICR,SZG3ICR,SZG4ICR,
     x                           NG2CHK,NG3CHK,NG4CHK,
     x                           DAE,DBE,DP,
     x                           GM2_1ICR,GM2_2ICR,GM2sICR,
     x                           GM3_1ICR,GM3_2ICR,
     x                           GM4ICR,
     x                           XFP,XSP,
     x                           E_P_OMG1, 
     x                           E_P_OMG2,
     x                           E_P_OMG3, 
     x                           E_P_OMG4,
     x                           S_P_OMG1,
     x                           S_P_OMG2)

! Construct Electronic Fock Matrix (special electron)
      call RXCHF_IC_construct_FBE(LNEOHF,nelec,NAE,NBE,
     x                            nebf,npbf,nebf2,ngee,
     x                            ng1,ng2,ng3,ng4,
     x                            SZG2ICR,SZG3ICR,SZG4ICR,
     x                            NG2CHK,NG3CHK,NG4CHK,
     x                            DAE,DBE,DP,
     x                            GM2_1ICR,GM2_2ICR,GM2sICR,
     x                            GM3_1ICR,GM3_2ICR,
     x                            GM4ICR,
     x                            XFBE,XSBE,
     x                            E_BE_OMG1, 
     x                            E_BE_OMG2,
     x                            E_BE_OMG3, 
     x                            E_BE_OMG4,
     x                            S_BE_OMG1,
     x                            S_BE_OMG2)

      if (rxchfdbg) then
       write(*,*)
       write(*,*) 'E_P_OMG1  =',E_P_OMG1  
       write(*,*) 'E_BE_OMG1 =',E_BE_OMG1 
       write(*,*)
       write(*,*) 'E_P_OMG2  =',E_P_OMG2 
       write(*,*) 'E_AE_OMG2 =',E_AE_OMG2 
       write(*,*) 'E_BE_OMG2 =',E_BE_OMG2 
       write(*,*)
       write(*,*) 'E_P_OMG3  =',E_P_OMG3 
       write(*,*) 'E_AE_OMG3 =',E_AE_OMG3 
       write(*,*) 'E_BE_OMG3 =',E_BE_OMG3 
       write(*,*)
       write(*,*) 'E_P_OMG4  =',E_P_OMG4 
       write(*,*) 'E_AE_OMG4 =',E_AE_OMG4 
       write(*,*) 'E_BE_OMG4 =',E_BE_OMG4 
       write(*,*)
       write(*,*) 'S_P_OMG1  =',S_P_OMG1 
       write(*,*) 'S_BE_OMG1 =',S_BE_OMG1 
       write(*,*)
       write(*,*) 'S_P_OMG2  =',S_P_OMG2 
       write(*,*) 'S_AE_OMG2 =',S_AE_OMG2 
       write(*,*) 'S_BE_OMG2 =',S_BE_OMG2 
       write(*,*)
      end if

      Psi_HOMG_Psi = E_P_OMG1 + E_P_OMG2 + E_P_OMG3 + E_P_OMG4
      S_total      = S_P_OMG1 + S_P_OMG2

      E_OMG1  = E_P_OMG1 / S_total
      E_OMG2  = E_P_OMG2 / S_total
      E_OMG3  = E_P_OMG3 / S_total
      E_OMG4  = E_P_OMG4 / S_total
      S_OMG1  = S_P_OMG1
      S_OMG2  = S_P_OMG2

      E_total = E_OMG1 + E_OMG2 + E_OMG3 + E_OMG4

! Correct Regular Electron Fock Matrix
      call RXCHF_Fock_correction1(nebf,Psi_HOMG_Psi,S_total,
     x                            XFAE,XSAE,FAE)

! Correct QM Particle Fock Matrix
      call RXCHF_Fock_correction2(npbf,Psi_HOMG_Psi,S_total,
     x                            XFP,XSP,FP)

! Correct Special Electron Fock Matrix
      call RXCHF_Fock_correction2(nebf,Psi_HOMG_Psi,S_total,
     x                           XFBE,XSBE,FBE)

! Fock testing
      if (LCMF) then
!        call RXCHF_Fock_testing1(nebf,FAE,DAE,E_ecore,E_GAMee,E_OMG2)

        call RXCHF_Fock_testing2(npbf,FP,DP)

        call RXCHF_Fock_testing2(nebf,FBE,DBE)

        call UFM_sym_check(nebf,npbf,FAE,FBE,FP)
      end if

C ARS( testing
      if (rxchfdbg) then
        nebflt=nebf*(nebf+1)/2
        npbflt=npbf*(npbf+1)/2
        write(*,*) "DAE:"
        call prt_lower_triangle(nebf,nebflt,DAE)
        write(*,*)
        write(*,*) "DBE:"
        call prt_lower_triangle(nebf,nebflt,DBE)
        write(*,*)
        write(*,*) "DP:"
        call prt_lower_triangle(npbf,npbflt,DP)
        write(*,*)
        write(*,*) "FAE:"
        call prt_lower_triangle(nebf,nebflt,FAE)
        write(*,*)
        write(*,*) "FBE:"
        call prt_lower_triangle(nebf,nebflt,FBE)
        write(*,*)
        write(*,*) "FP:"
        call prt_lower_triangle(npbf,npbflt,FP)
      end if
C )
      return
      end
!======================================================================
      subroutine RXCHF_Fock_correction1(nbf,Psi_HOMG_Psi,S_total,
     x                                  XF,XS,F)
C Corrects FAE (with overlap terms and Gamma_2s terms)
!======================================================================
      implicit none

! Input Variables
      integer nbf
      double precision Psi_HOMG_Psi
      double precision S_total
      double precision XS(nbf,nbf)
      double precision XF(nbf,nbf)

! Variables Returned
      double precision F(nbf,nbf)

! Local Variables
      integer i,j
      double precision COEF1,COEF2
      double precision zero,one
      parameter(zero=0.0d+00,one=1.0d+00)

! Initialize
      F=zero
      COEF1=one/S_total
      COEF2=Psi_HOMG_Psi / (S_Total * S_total)

! Fock Matrix Correction
      do i=1,nbf
         do j=1,nbf

            F(j,i) = COEF1*XF(j,i) + COEF2*XS(j,i)

         end do
      end do


      return
      end
!======================================================================
      subroutine RXCHF_Fock_correction2(nbf,Psi_HOMG_Psi,S_total,
     x                                  XF,XS,F)
C Corrects FP and FEB (with overlap terms and Gamma_1s/Gamma_2s terms)
!======================================================================
      implicit none

! Input Variables
      integer nbf
      double precision Psi_HOMG_Psi
      double precision S_total
      double precision XS(nbf,nbf)
      double precision XF(nbf,nbf)

! Variables Returned
      double precision F(nbf,nbf)

! Local Variables
      integer i,j
      double precision COEF1,COEF2
      double precision zero,one
      parameter(zero=0.0d+00,one=1.0d+00)

! Initialize
      F=zero
      COEF1=one/S_total
      COEF2=Psi_HOMG_Psi / (S_Total * S_total)

! Fock Matrix Correction
      do i=1,nbf
         do j=1,nbf

            F(j,i)= COEF1*XF(j,i) - COEF2*XS(j,i)

         end do
      end do


      return
      end
!======================================================================
      subroutine RXCHF_Fock_testing1(nbf,F,D,E_ecore,E_GAMee,E_OMG2)
C Tests FAE and DAE
!======================================================================
      implicit none

! Input Variables
      integer nbf
      double precision F(nbf,nbf)
      double precision D(nbf,nbf)
      double precision E_ecore,E_GAMee
      double precision E_OMG2 ! On entry already div by S_total

! Local variables
      integer i,j
      double precision fdsum,ans
      double precision zero,two
      parameter(zero=0.0d+00,two=2.0d+00)

      fdsum=zero

      do i=1,nbf
        do j=1,nbf

          fdsum=fdsum+D(j,i)*F(j,i)

        end do
      end do

      ans=E_ecore+two*E_GAMee+E_OMG2

      write(*,*) "Testing Fock matrix routine 1:"
      write(*,*) fdsum,ans

      return
      end
!======================================================================
      subroutine RXCHF_Fock_testing2(nbf,F,D)
C Tests FP and DP or FBE and DBE
!======================================================================
      implicit none

! Input Variables
      integer nbf
      double precision F(nbf,nbf)
      double precision D(nbf,nbf)

! Local variables
      integer i,j
      double precision fdsum,ans
      double precision zero
      parameter(zero=0.0d+00)

      fdsum=zero

      do i=1,nbf
        do j=1,nbf

          fdsum=fdsum+D(j,i)*F(j,i)

        end do
      end do

      ans=zero

      write(*,*) "Testing Fock matrix routine 2:"
      write(*,*) fdsum,ans

      return
      end
!======================================================================
      subroutine RXCHF_IC_construct_FAE(LNEOHF,nelec,NAE,NBE,
     x                              nebf,npbf,nebf2,ngee,
     x                              ng1,ng2,ng3,ng4,
     x                              SZG2ICR,SZG3ICR,SZG4ICR,
     x                              NG2CHK,NG3CHK,NG4CHK,
     x                              DAE,DBE,DP,
     x                              GM2_1ICR,GM2_2ICR,GM2sICR,
     x                              GM3_1ICR,GM3_2ICR,
     x                              GM4ICR,
     x                              FAE,SAE,
     x                              E_AE_OMG2,
     x                              E_AE_OMG3,
     x                              E_AE_OMG4,
     x                              S_AE_OMG2)

! In-core Construct Regular Electronic Fock Matrix
! nelec: total number of electrons
! NAE: number regular electrons
! NBE: number special electrons
!
!======================================================================
      implicit none
! Input Variables
      logical LNEOHF
      integer nebf,npbf,nebf2,ngee,ng1,ng2,ng3,ng4
      integer SZG2ICR,SZG3ICR,SZG4ICR
      integer nelec,NAE,NBE
      integer NG2CHK,NG3CHK,NG4CHK
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DETOT(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM2_1ICR(SZG2ICR),GM2_2ICR(SZG2ICR)
      double precision GM2sICR(SZG2ICR)
      double precision GM3_1ICR(SZG3ICR),GM3_2ICR(SZG3ICR)
      double precision GM4ICR(SZG4ICR)

! Variables Returned
      double precision FAE(nebf,nebf)
      double precision SAE(nebf,nebf)
      double precision E_AE_OMG2
      double precision E_AE_OMG3
      double precision E_AE_OMG4
      double precision S_AE_OMG2

! Local Variables
      double precision zero
      parameter(zero=0.0d+00)

! Initialize
      FAE=zero
      SAE=zero
      E_AE_OMG2  = zero
      E_AE_OMG3  = zero
      E_AE_OMG4  = zero
      S_AE_OMG2  = zero

      call RXCHF_FAE_OMG2(NG2CHK,nebf,npbf,ng2,
     x                    DAE,DBE,DP,GM2_1ICR,GM2_2ICR,GM2sICR,
     x                    FAE,SAE,E_AE_OMG2,S_AE_OMG2)

      if (NAE.gt.1) then

         call RXCHF_FAE_OMG3(NG3CHK,nebf,npbf,ng3,
     x                       DAE,DBE,DP,GM3_1ICR,GM3_2ICR,
     x                       FAE,E_AE_OMG3)

         if (NAE.gt.2) then

            call RXCHF_FAE_OMG4(NG4CHK,nebf,npbf,ng4,
     x                          DAE,DBE,DP,GM4ICR,FAE,E_AE_OMG4)

         end if

      end if

      return
      end
!======================================================================
      subroutine RXCHF_IC_construct_FP(LNEOHF,nelec,NAE,NBE,
     x                             nebf,npbf,npbf2,ngee,
     x                             ng1,ng2,ng3,ng4,
     x                             SZG2ICR,SZG3ICR,SZG4ICR,
     x                             NG2CHK,NG3CHK,NG4CHK,
     x                             DAE,DBE,DP,
     x                             GM2_1ICR,GM2_2ICR,GM2sICR,
     x                             GM3_1ICR,GM3_2ICR,
     x                             GM4ICR,
     x                             FP,SP,
     x                             E_P_OMG1, 
     x                             E_P_OMG2,
     x                             E_P_OMG3, 
     x                             E_P_OMG4,
     x                             S_P_OMG1,
     x                             S_P_OMG2)

! In-core Construct QM-Particle Fock Matrix
! nelec: total number of electrons
! NAE: number regular electrons
! NBE: number special electrons
!
!======================================================================
      implicit none
! Input Variables
      logical LNEOHF
      integer nebf,npbf,npbf2,ngee,ng1,ng2,ng3,ng4
      integer SZG2ICR,SZG3ICR,SZG4ICR
      integer nelec,NAE,NBE
      integer NG2CHK,NG3CHK,NG4CHK
      double precision DAE(nebf,nebf),DBE(nebf,nebf),DP(npbf,npbf)
      double precision GM2_1ICR(SZG2ICR),GM2_2ICR(SZG2ICR) 
      double precision GM2sICR(SZG2ICR)
      double precision GM3_1ICR(SZG2ICR),GM3_2ICR(SZG2ICR)
      double precision GM4ICR(SZG4ICR)

! Variables Returned
      double precision FP(npbf,npbf)
      double precision SP(npbf,npbf)
      double precision E_P_OMG1 
      double precision E_P_OMG2
      double precision E_P_OMG3
      double precision E_P_OMG4
      double precision S_P_OMG1
      double precision S_P_OMG2

! Local Variables
      double precision zero
      parameter(zero=0.0d+00)

! Initialize
      FP=zero    
      SP=zero    
      E_P_OMG1  = zero 
      E_P_OMG2  = zero 
      E_P_OMG3  = zero 
      E_P_OMG4  = zero 
      S_P_OMG1  = zero
      S_P_OMG2  = zero
 
      call RXCHF_FP_OMG1(nebf,npbf,ng1,DBE,DP,FP,SP,
     x                   E_P_OMG1,S_P_OMG1)

      call RXCHF_FP_OMG2(NG2CHK,nebf,npbf,ng2,
     x                   DAE,DBE,DP,GM2_1ICR,GM2_2ICR,GM2sICR,
     x                   FP,SP,E_P_OMG2,S_P_OMG2)

      if (NAE.gt.1) then

         call RXCHF_FP_OMG3(NG3CHK,nebf,npbf,ng3,
     x                      DAE,DBE,DP,GM3_1ICR,GM3_2ICR,
     x                      FP,E_P_OMG3)

         if (NAE.gt.2) then

            call RXCHF_FP_OMG4(NG4CHK,nebf,npbf,ng4,
     x                         DAE,DBE,DP,GM4ICR,FP,E_P_OMG4)

         end if

      end if

      return
      end
!======================================================================
      subroutine RXCHF_IC_construct_FBE(LNEOHF,nelec,NAE,NBE,
     x                             nebf,npbf,nebf2,ngee,
     x                             ng1,ng2,ng3,ng4,
     x                             SZG2ICR,SZG3ICR,SZG4ICR,
     x                             NG2CHK,NG3CHK,NG4CHK,
     x                             DAE,DBE,DP,
     x                             GM2_1ICR,GM2_2ICR,GM2sICR,
     x                             GM3_1ICR,GM3_2ICR,
     x                             GM4ICR,
     x                             FBE,SBE,
     x                             E_BE_OMG1, 
     x                             E_BE_OMG2,
     x                             E_BE_OMG3, 
     x                             E_BE_OMG4,
     x                             S_BE_OMG1,
     x                             S_BE_OMG2)

! In-core Construct special electronic Fock Matrix
! nelec: total number of electrons
! NAE: number regular electrons
! NBE: number special electrons
!
!======================================================================
      implicit none
! Input Variables
      logical LNEOHF
      integer nebf,npbf,nebf2,ngee,ng1,ng2,ng3,ng4
      integer SZG2ICR,SZG3ICR,SZG4ICR
      integer nelec,NAE,NBE
      integer NG2CHK,NG3CHK,NG4CHK
      double precision DAE(nebf,nebf),DBE(nebf,nebf),DP(npbf,npbf)
      double precision GM2_1ICR(SZG2ICR),GM2_2ICR(SZG2ICR)
      double precision GM2sICR(SZG2ICR)
      double precision GM3_1ICR(SZG2ICR),GM3_2ICR(SZG2ICR)
      double precision GM4ICR(SZG4ICR)

! Variables Returned
      double precision FBE(nebf,nebf)
      double precision SBE(nebf,nebf)
      double precision E_BE_OMG1 
      double precision E_BE_OMG2
      double precision E_BE_OMG3 
      double precision E_BE_OMG4
      double precision S_BE_OMG1
      double precision S_BE_OMG2

! Local Variables
      double precision zero
      parameter(zero=0.0d+00)

! Initialize
      FBE=zero    
      SBE=zero    
      E_BE_OMG1  = zero 
      E_BE_OMG2  = zero 
      E_BE_OMG3  = zero 
      E_BE_OMG4  = zero 
      S_BE_OMG1  = zero
      S_BE_OMG2  = zero
 
      call RXCHF_FBE_OMG1(nebf,npbf,ng1,DBE,DP,FBE,SBE,
     x                   E_BE_OMG1,S_BE_OMG1)

      call RXCHF_FBE_OMG2(NG2CHK,nebf,npbf,ng2,
     x                   DAE,DBE,DP,GM2_1ICR,GM2_2ICR,GM2sICR,
     x                   FBE,SBE,E_BE_OMG2,S_BE_OMG2)

      if (NAE.gt.1) then

         call RXCHF_FBE_OMG3(NG3CHK,nebf,npbf,ng3,
     x                      DAE,DBE,DP,GM3_1ICR,GM3_2ICR,
     x                      FBE,E_BE_OMG3)

         if (NAE.gt.2) then

            call RXCHF_FBE_OMG4(NG4CHK,nebf,npbf,ng4,
     x                          DAE,DBE,DP,GM4ICR,FBE,E_BE_OMG4)

         end if

      end if

      return
      end
