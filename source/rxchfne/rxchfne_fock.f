!======================================================================
      subroutine RXCHFne_FOCK(LCMF,LNEOHF,nelec,NAE,NBE,
     x                    nebf,nebf2,npbf,npbf2,ngee,
     x                    ng1,ng2,ng3,ng4,
     x                    SZG2ICR,
     x                    NG2CHK,
     x                    DAE,DBE,DP,
     x                    GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                    GM2ICR,
     x                    FP,FAE,FBE, 
     x                    E_ecore,
     x                    E_GAMee,
     x                    E_OMG1, 
     x                    E_OMG2,
     x                    E_total,
     x                    S_total)

!======================================================================
      implicit none

! Input Variables
      logical LCMF,LNEOHF
      integer nebf,npbf,npbf2,nebf2,ngee,ng1,ng2,ng3,ng4
      integer SZG2ICR  !,SZG3ICR,SZG4ICR
      integer nelec,NAE,NBE
      integer NG2CHK
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

      double precision GM2ICR(SZG2ICR)
C      double precision GM2sICR(SZG2ICR)
C      double precision GM3ICR(SZG3ICR)
C      double precision GM4ICR(SZG4ICR)

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
      double precision XFAE1(nebf,nebf)
      double precision XFAE2(nebf,nebf)
      double precision XFBE(nebf,nebf)
      double precision XFP(npbf,npbf)

      double precision XSBE(nebf,nebf)
      double precision XSP(npbf,npbf)

      double precision E_AE_ecore
      double precision E_AE_GAMee
      double precision E_AE_OMG2

!     double precision E_BE
!     double precision OVLP_BE
      double precision E_BE_ecore
      double precision E_BE_GAMep
      double precision E_BE_GAMee
      double precision E_BE_OMG1
      double precision E_BE_OMG2
      double precision E_BE_OMG3
      double precision E_BE_OMG4
      double precision S_BE_OMG1
      double precision S_BE_OMG2

!     double precision E_P
!     double precision OVLP_P
      double precision E_P_pcore
      double precision E_P_GAMep
      double precision E_P_OMG1
      double precision E_P_OMG2
C      double precision E_P_OMG3
C      double precision E_P_OMG4
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

! Initialize
      FAE=zero 
      FBE=zero 
      FP=zero 

      XFAE1=zero 
      XFAE2=zero 
      XFP=zero 
      XSP=zero 
      XFBE=zero 
      XSBE=zero 

      E_AE_ecore = zero  
      E_AE_GAMee = zero 
      E_AE_OMG2 = zero 

      E_P_OMG1 = zero 
      E_P_OMG2 = zero 
      S_P_OMG1 = zero 

      E_BE_OMG1 = zero 
      E_BE_OMG2 = zero 
      S_BE_OMG1 = zero 

! Construct Electronic Fock Matrix (regular electrons)
      call RXCHFne_IC_construct_FAE(LNEOHF,nelec,NAE,NBE,
     x                              nebf,npbf,nebf2,ngee,
     x                              ng1,ng2,ng3,ng4,
     x                              SZG2ICR,
     x                              NG2CHK,
     x                              DAE,DBE,DP,
     x                              GAM_ecore,GAM_ee,
     x                              GM2ICR,
     x                              XFAE1,XFAE2,
     x                              E_AE_ecore,
     x                              E_AE_GAMee,
     x                              E_AE_OMG2)

! Construct QM Particle Fock Matrix
      call RXCHFne_IC_construct_FP(LNEOHF,nelec,NAE,NBE,
     x                           nebf,npbf,npbf2,ngee,
     x                           ng1,ng2,ng3,ng4,
     x                           SZG2ICR,
     x                           NG2CHK,
     x                           DAE,DBE,DP,
     x                           GM2ICR,
     x                           XFP,XSP,
     x                           E_P_OMG1, 
     x                           E_P_OMG2,
     x                           S_P_OMG1)

! Construct Electronic Fock Matrix (special electron)
      call RXCHFne_IC_construct_FBE(LNEOHF,nelec,NAE,NBE,
     x                            nebf,npbf,nebf2,ngee,
     x                            ng1,ng2,ng3,ng4,
     x                            SZG2ICR,
     x                            NG2CHK,
     x                            DAE,DBE,DP,
     x                            GM2ICR,
     x                            XFBE,XSBE,
     x                            E_BE_OMG1, 
     x                            E_BE_OMG2,
     x                            S_BE_OMG1)

      if (rxchfdbg) then
       write(*,*)
       write(*,*) 'E_AE_ecore=',E_AE_ecore
       write(*,*) 'E_AE_GAMee=',E_AE_GAMee
       write(*,*)
       write(*,*) 'E_P_OMG1  =',E_P_OMG1  
       write(*,*) 'E_BE_OMG1 =',E_BE_OMG1 
       write(*,*)
       write(*,*) 'E_P_OMG2  =',E_P_OMG2 
       write(*,*) 'E_AE_OMG2 =',E_AE_OMG2 
       write(*,*) 'E_BE_OMG2 =',E_BE_OMG2 
       write(*,*)
       write(*,*) 'S_P_OMG1  =',S_P_OMG1 
       write(*,*) 'S_BE_OMG1 =',S_BE_OMG1 
       write(*,*)
      end if

      Psi_HOMG_Psi = E_P_OMG1 + E_P_OMG2  ! (or E_BE_OMG1 + E_BE_OMG2)
      S_total      = S_P_OMG1             ! (or S_BE_OMG1)

      E_ecore = E_AE_ecore
      E_GAMee = E_AE_GAMee
      E_OMG1  = E_P_OMG1 / S_total        ! (or E_BE_OMG1 / S_total)
      E_OMG2  = E_P_OMG2 / S_total        ! (or E_BE_OMG2 / S_total)
                                          ! (or E_AE_OMG2 / S_total)

      E_total = E_ecore + E_GAMee + E_OMG1 + E_OMG2


! Correct Regular Electron Fock Matrix
      call RXCHFne_Fock_correction1(nebf,S_total,
     x                            XFAE1,XFAE2,FAE)

! Correct QM Particle Fock Matrix
      call RXCHFne_Fock_correction2(npbf,Psi_HOMG_Psi,S_total,
     x                            XFP,XSP,FP)

! Correct Special Electron Fock Matrix
      call RXCHFne_Fock_correction2(nebf,Psi_HOMG_Psi,S_total,
     x                           XFBE,XSBE,FBE)

! Fock testing
      if (LCMF) then
        call RXCHF_Fock_testing1(nebf,FAE,DAE,E_ecore,E_GAMee,E_OMG2)

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
      subroutine RXCHFne_Fock_correction1(nbf,S_total,
     x                                  XF1,XF2,F)
C Corrects FEA (with overlap terms)
!======================================================================
      implicit none

! Input Variables
      integer nbf
      double precision S_total
      double precision XF1(nbf,nbf)
      double precision XF2(nbf,nbf)

! Variables Returned
      double precision F(nbf,nbf)

! Local Variables
      integer i,j
      double precision COEF1
      double precision zero,one
      parameter(zero=0.0d+00,one=1.0d+00)

! Initialize
      F=zero
      COEF1=one/S_total

! Fock Matrix Correction
      do i=1,nbf
         do j=1,nbf

            F(j,i) = XF1(j,i) + COEF1*XF2(j,i)

         end do
      end do


      return
      end
!======================================================================
      subroutine RXCHFne_Fock_correction2(nbf,Psi_HOMG_Psi,S_total,
     x                                  XF,XS,F)
C Corrects FP and FEB (with overlap terms and Gamma_1S terms)
!======================================================================
      implicit none

! Input Variables
      integer nbf
      double precision Psi_HOMG_Psi
      double precision S_total
      double precision XS(nbf,nbf)
      double precision XF(nbf,nbf)
      double precision E_ecore,E_GAMee

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
      subroutine RXCHFne_IC_construct_FAE(LNEOHF,nelec,NAE,NBE,
     x                              nebf,npbf,nebf2,ngee,
     x                              ng1,ng2,ng3,ng4,
     x                              SZG2ICR,
     x                              NG2CHK,
     x                              DAE,DBE,DP,
     x                              GAM_ecore,GAM_ee,
     x                              GM2ICR,
     x                              FAE1,FAE2,
     x                              E_AE_ecore,
     x                              E_AE_GAMee,
     x                              E_AE_OMG2)

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
      integer SZG2ICR
      integer nelec,NAE,NBE
      integer NG2CHK
      double precision DAE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DETOT(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GAM_ecore(nebf2)
      double precision GAM_ee(ngee)
      double precision GM2ICR(SZG2ICR)

! Variables Returned
      double precision FAE1(nebf,nebf)
      double precision FAE2(nebf,nebf)
      double precision E_AE_ecore
      double precision E_AE_GAMee
      double precision E_AE_OMG2

! Local Variables
      double precision zero
      parameter(zero=0.0d+00)

! Initialize
      E_AE_ecore = zero
      E_AE_GAMee = zero
      E_AE_OMG2  = zero
      FAE1=zero
      FAE2=zero

      call FAE_ecore(nebf,nebf2,GAM_ecore,DAE,FAE1,E_AE_ecore)

      if (NAE.gt.1) then
       call RXCHFne_FAE_GAMee(nebf,ngee,GAM_ee,DAE,FAE1,E_AE_GAMee)
      end if

      call RXCHFne_FAE_OMG2(NG2CHK,nebf,npbf,ng2,
     x                    DAE,DBE,DP,GM2ICR,FAE2,E_AE_OMG2)

      return
      end
!======================================================================
      subroutine RXCHFne_IC_construct_FP(LNEOHF,nelec,NAE,NBE,
     x                             nebf,npbf,npbf2,ngee,
     x                             ng1,ng2,ng3,ng4,
     x                             SZG2ICR,
     x                             NG2CHK,
     x                             DAE,DBE,DP,
     x                             GM2ICR,
     x                             FP,SP,
     x                             E_P_OMG1, 
     x                             E_P_OMG2,
     x                             S_P_OMG1)

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
      integer SZG2ICR
      integer nelec,NAE,NBE
      integer NG2CHK
      double precision DAE(nebf,nebf),DBE(nebf,nebf),DP(npbf,npbf)
      double precision GM2ICR(SZG2ICR)

! Variables Returned
      double precision FP(npbf,npbf)
      double precision SP(npbf,npbf)
      double precision E_P_OMG1 
      double precision E_P_OMG2
      double precision S_P_OMG1

! Local Variables
      double precision zero
      parameter(zero=0.0d+00)

! Initialize
      FP=zero    
      SP=zero    
      E_P_OMG1  = zero 
      E_P_OMG2  = zero 
      S_P_OMG1  = zero
 
      call RXCHF_FP_OMG1(nebf,npbf,ng1,DBE,DP,FP,SP,
     x                   E_P_OMG1,S_P_OMG1)

      call RXCHFne_FP_OMG2(NG2CHK,nebf,npbf,ng2,NAE,NBE,
     x                   DAE,DBE,DP,GM2ICR,FP,E_P_OMG2)

      return
      end
!======================================================================
      subroutine RXCHFne_IC_construct_FBE(LNEOHF,nelec,NAE,NBE,
     x                             nebf,npbf,nebf2,ngee,
     x                             ng1,ng2,ng3,ng4,
     x                             SZG2ICR,
     x                             NG2CHK,
     x                             DAE,DBE,DP,
     x                             GM2ICR,
     x                             FBE,SBE,
     x                             E_BE_OMG1, 
     x                             E_BE_OMG2,
     x                             S_BE_OMG1)

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
      integer SZG2ICR
      integer nelec,NAE,NBE
      integer NG2CHK
      double precision DAE(nebf,nebf),DBE(nebf,nebf),DP(npbf,npbf)
      double precision GM2ICR(SZG2ICR)

! Variables Returned
      double precision FBE(nebf,nebf)
      double precision SBE(nebf,nebf)
      double precision E_BE_OMG1 
      double precision E_BE_OMG2
      double precision S_BE_OMG1

! Local Variables
      double precision zero
      parameter(zero=0.0d+00)

! Initialize
      FBE=zero    
      SBE=zero    
      E_BE_OMG1  = zero 
      E_BE_OMG2  = zero 
      S_BE_OMG1  = zero
 
      call RXCHF_FBE_OMG1(nebf,npbf,ng1,DBE,DP,FBE,SBE,
     x                   E_BE_OMG1,S_BE_OMG1)

      call RXCHFne_FBE_OMG2(NG2CHK,nebf,npbf,ng2,NAE,NBE,
     x                   DAE,DBE,DP,GM2ICR,FBE,E_BE_OMG2)

      return
      end
