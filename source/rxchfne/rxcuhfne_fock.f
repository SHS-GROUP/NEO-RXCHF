!======================================================================
      subroutine RXCUHF_FOCK(LCMF,LNEOHF,LADDEXCH,nelec,NBE,
     x                       NAalpE,NAbetE,nebf,nebf2,npbf,npbf2,ngee,
     x                       ng1,ng2,ng3,ng4,
     x                       SZG2ICR,
     x                       NG2CHK,
     x                       DAalpE,DAbetE,DBE,DP,
     x                       GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                       GM2ICR,GM2sICR,
     x                       FP,FAalpE,FAbetE,FBE, 
     x                       E_ecore,
     x                       E_GAMee,
     x                       E_OMG1, 
     x                       E_OMG2,
     x                       E_total,
     x                       S_total)

C NAalpE : num alpha regular electrons
C NAbetE : num beta regular electrons
C NBE    : num special electrons (=1 with spin alpha)
!======================================================================
      implicit none

! Input Variables
      logical LCMF,LNEOHF,LADDEXCH
      integer nebf,npbf,npbf2,nebf2,ngee,ng1,ng2,ng3,ng4
      integer SZG2ICR  !,SZG3ICR,SZG4ICR
      integer nelec,NBE,NAalpE,NAbetE
      integer NG2CHK
      double precision DAtotE(nebf,nebf)
      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
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

      double precision GM2ICR(SZG2ICR)    ! Regular OMG2 terms
      double precision GM2sICR(SZG2ICR)   ! Exchange OMG2 terms
C      double precision GM3ICR(SZG3ICR)
C      double precision GM4ICR(SZG4ICR)

! Variables Returned
      double precision FAalpE(nebf,nebf)
      double precision FAbetE(nebf,nebf)
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
      double precision XFAalpE1(nebf,nebf)
      double precision XFAalpE2(nebf,nebf)
      double precision XFAbetE1(nebf,nebf)
      double precision XFAbetE2(nebf,nebf)
      double precision XFBE(nebf,nebf)
      double precision XFP(npbf,npbf)

      double precision XSBE(nebf,nebf)
      double precision XSP(npbf,npbf)

      double precision E_AalpE_ecore
      double precision E_AalpE_GAMee
      double precision E_AbetE_ecore
      double precision E_AbetE_GAMee

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
      FAalpE=zero 
      FAbetE=zero 
      FBE=zero 
      FP=zero 

      XFAalpE1=zero 
      XFAalpE2=zero 
      XFAbetE1=zero 
      XFAbetE2=zero 
      XFP=zero 
      XSP=zero 
      XFBE=zero 
      XSBE=zero 

      E_AalpE_ecore = zero  
      E_AalpE_GAMee = zero 
      E_AbetE_ecore = zero  
      E_AbetE_GAMee = zero 

      E_AE_OMG2 = zero 

      E_P_OMG1 = zero 
      E_P_OMG2 = zero 
      S_P_OMG1 = zero 

      E_BE_OMG1 = zero 
      E_BE_OMG2 = zero 
      S_BE_OMG1 = zero 

! Form total regular electronic density 
      DAtotE=0.0d+00
      call add2fock(nebf,DAalpE,DAtotE)
      call add2fock(nebf,DAbetE,DAtotE)

! Construct alpha/beta Electronic Fock Matrix (regular electrons)
      call RXCUHFne_IC_construct_FE(LNEOHF,LADDEXCH,
     x                            nelec,NAalpE,NAbetE,NBE,
     x                            nebf,npbf,nebf2,ngee,
     x                            ng1,ng2,ng3,ng4,
     x                            SZG2ICR,
     x                            NG2CHK,
     x                            DAalpE,DAbetE,DAtotE,DBE,DP,
     x                            GAM_ecore,GAM_ee,
     x                            GM2ICR,GM2sICR,
     x                            XFAalpE1,XFAalpE2,
     x                            XFAbetE1,XFAbetE2,
     x                            E_AalpE_ecore,
     x                            E_AalpE_GAMee,
     x                            E_AbetE_ecore,
     x                            E_AbetE_GAMee,
     x                            E_AE_OMG2)

! Construct QM Particle Fock Matrix
      call RXCUHFne_IC_construct_FP(LNEOHF,LADDEXCH,
     x                            nelec,NAbetE,NAalpE,NBE,
     x                            nebf,npbf,npbf2,ngee,
     x                            ng1,ng2,ng3,ng4,
     x                            SZG2ICR,
     x                            NG2CHK,
     x                            DAalpE,DAbetE,DAtotE,DBE,DP,
     x                            GM2ICR,GM2sICR,
     x                            XFP,XSP,
     x                            E_P_OMG1, 
     x                            E_P_OMG2,
     x                            S_P_OMG1)

! Construct Electronic Fock Matrix (special electron)
      call RXCUHFne_IC_construct_FBE(LNEOHF,LADDEXCH,
     x                             nelec,NAbetE,NAalpE,NBE,
     x                             nebf,npbf,nebf2,ngee,
     x                             ng1,ng2,ng3,ng4,
     x                             SZG2ICR,
     x                             NG2CHK,
     x                             DAalpE,DAbetE,DAtotE,DBE,DP,
     x                             GM2ICR,GM2sICR,
     x                             XFBE,XSBE,
     x                             E_BE_OMG1, 
     x                             E_BE_OMG2,
     x                             S_BE_OMG1)

      if (rxchfdbg) then
       write(*,*)
       write(*,*) 'E_AalpE_ecore=',E_AalpE_ecore
       write(*,*) 'E_AbetE_ecore=',E_AbetE_ecore
       write(*,*) 'E_AE_ecore=',E_AalpE_ecore+E_AbetE_ecore
       write(*,*)
       write(*,*) 'E_AalpE_GAMee=',E_AalpE_GAMee
       write(*,*) 'E_AbetE_GAMee=',E_AbetE_GAMee
       write(*,*) 'E_AE_GAMee=',E_AalpE_GAMee+E_AbetE_GAMee
       write(*,*)
       write(*,*) 'E_P_OMG1  =',E_P_OMG1  
       write(*,*) 'E_BE_OMG1 =',E_BE_OMG1 
       write(*,*)
       write(*,*) 'E_AE_OMG2 =',E_AE_OMG2
       write(*,*) 'E_P_OMG2  =',E_P_OMG2 
       write(*,*) 'E_BE_OMG2 =',E_BE_OMG2 
       write(*,*)
       write(*,*) 'S_P_OMG1  =',S_P_OMG1 
       write(*,*) 'S_BE_OMG1 =',S_BE_OMG1 
       write(*,*)
      end if

      Psi_HOMG_Psi = E_P_OMG1 + E_P_OMG2  ! (or E_BE_OMG1 + E_BE_OMG2)
      S_total      = S_P_OMG1             ! (or S_BE_OMG1)

      E_ecore = E_AalpE_ecore + E_AbetE_ecore
      E_GAMee = E_AalpE_GAMee + E_AbetE_GAMee
      E_OMG1  = E_P_OMG1 / S_total        ! (or E_BE_OMG1 / S_total)
      E_OMG2  = E_P_OMG2 / S_total        ! (or E_BE_OMG2 / S_total)
                                          ! (or E_AE_OMG2 / S_total)

      E_total = E_ecore + E_GAMee + E_OMG1 + E_OMG2


! Correct Regular Electron Fock Matrix
      call RXCHFne_Fock_correction1(nebf,S_total,
     x                            XFAalpE1,XFAalpE2,FAalpE)
      call RXCHFne_Fock_correction1(nebf,S_total,
     x                            XFAbetE1,XFAbetE2,FAbetE)

! Correct QM Particle Fock Matrix
      call RXCHFne_Fock_correction2(npbf,Psi_HOMG_Psi,S_total,
     x                            XFP,XSP,FP)

! Correct Special Electron Fock Matrix
      call RXCHFne_Fock_correction2(nebf,Psi_HOMG_Psi,S_total,
     x                           XFBE,XSBE,FBE)

! Fock testing
      if (LCMF) then
!        call RXCHF_Fock_testing1(nebf,FAE,DAE,E_ecore,E_GAMee,E_OMG2)

!        call RXCHF_Fock_testing2(npbf,FP,DP)

!        call RXCHF_Fock_testing2(nebf,FBE,DBE)

        call RXCUHF_sym_check(nebf,npbf,FAalpE,FAbetE,FBE,FP)
      end if

C ARS( testing
      if (rxchfdbg) then
        nebflt=nebf*(nebf+1)/2
        npbflt=npbf*(npbf+1)/2
        write(*,*) "DAalpE:"
        call prt_lower_triangle(nebf,nebflt,DAalpE)
        write(*,*)
        write(*,*) "DAbetE:"
        call prt_lower_triangle(nebf,nebflt,DAbetE)
        write(*,*)
        write(*,*) "DBE:"
        call prt_lower_triangle(nebf,nebflt,DBE)
        write(*,*)
        write(*,*) "DP:"
        call prt_lower_triangle(npbf,npbflt,DP)
        write(*,*)
        write(*,*) "FAalpE:"
        call prt_lower_triangle(nebf,nebflt,FAalpE)
        write(*,*)
        write(*,*) "FAbetE:"
        call prt_lower_triangle(nebf,nebflt,FAbetE)
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
      subroutine RXCUHFne_IC_construct_FE(LNEOHF,LADDEXCH,
     x                                  nelec,NAalpE,NAbetE,NBE,
     x                                  nebf,npbf,nebf2,ngee,
     x                                  ng1,ng2,ng3,ng4,
     x                                  SZG2ICR,
     x                                  NG2CHK,
     x                                  DAalpE,DAbetE,DAtotE,DBE,DP,
     x                                  GAM_ecore,GAM_ee,
     x                                  GM2ICR,GM2sICR,
     x                                  FAalpE1,FAalpE2,
     x                                  FAbetE1,FAbetE2,
     x                                  E_AalpE_ecore,
     x                                  E_AalpE_GAMee,
     x                                  E_AbetE_ecore,
     x                                  E_AbetE_GAMee,
     x                                  E_AE_OMG2)

! In-core Construct Regular Electronic Fock Matrix
! nelec: total number of electrons
! NAE: number regular electrons
! NBE: number special electrons
!
!======================================================================
      implicit none
! Input Variables
      logical LNEOHF,LADDEXCH
      integer nebf,npbf,nebf2,ngee,ng1,ng2,ng3,ng4
      integer SZG2ICR
      integer nelec,NAalpE,NAbetE,NBE
      integer NG2CHK
      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
      double precision DAtotE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GAM_ecore(nebf2)
      double precision GAM_ee(ngee)
      double precision GM2ICR(SZG2ICR)    ! Regular OMG2 terms
      double precision GM2sICR(SZG2ICR)   ! Exchange OMG2 terms

! Variables Returned
      double precision FAalpE1(nebf,nebf)
      double precision FAalpE2(nebf,nebf)
      double precision FAbetE1(nebf,nebf)
      double precision FAbetE2(nebf,nebf)
      double precision E_AalpE_ecore
      double precision E_AalpE_GAMee
      double precision E_AbetE_ecore
      double precision E_AbetE_GAMee
      double precision E_AE_OMG2

! Local Variables
      double precision E_AalpE_OMG2
      double precision E_AbetE_OMG2
      double precision zero
      parameter(zero=0.0d+00)

! Initialize
      E_AalpE_ecore = zero
      E_AalpE_GAMee = zero
      E_AbetE_ecore = zero
      E_AbetE_GAMee = zero
      FAalpE1=zero
      FAalpE2=zero
      FAbetE1=zero
      FAbetE2=zero
      E_AE_OMG2  = zero
      E_AalpE_OMG2  = zero
      E_AbetE_OMG2  = zero

      call FAE_ecore(nebf,nebf2,GAM_ecore,DAalpE,FAalpE1,E_AalpE_ecore)
      call FAE_ecore(nebf,nebf2,GAM_ecore,DAbetE,FAbetE1,E_AbetE_ecore)

      if ((NAalpE+NAbetE).gt.1) then
       call FAE_GAMee(nebf,ngee,GAM_ee,DAalpE,DAtotE,
     x                FAalpE1,E_AalpE_GAMee)
       call FAE_GAMee(nebf,ngee,GAM_ee,DAbetE,DAtotE,
     x                FAbetE1,E_AbetE_GAMee)
      end if

      if (LADDEXCH) then
        call RXCUHFne_FE_OMG2(NG2CHK,nebf,npbf,ng2,
     x                      DAalpE,DAbetE,DAtotE,DBE,DP,
     x                      GM2ICR,GM2sICR,
     x                      FAalpE2,FAbetE2,
     x                      E_AE_OMG2)
      else
        call RXCHFne_FAE_OMG2(NG2CHK,nebf,npbf,ng2,
     x                      DAalpE,DBE,DP,GM2ICR,FAalpE2,E_AalpE_OMG2)
        call RXCHFne_FAE_OMG2(NG2CHK,nebf,npbf,ng2,
     x                      DAbetE,DBE,DP,GM2ICR,FAbetE2,E_AbetE_OMG2)
        E_AE_OMG2=E_AalpE_OMG2+E_AbetE_OMG2
      end if

      return
      end
!======================================================================
      subroutine RXCUHFne_IC_construct_FP(LNEOHF,LADDEXCH,
     x                                  nelec,NAalpE,NAbetE,NBE,
     x                                  nebf,npbf,npbf2,ngee,
     x                                  ng1,ng2,ng3,ng4,
     x                                  SZG2ICR,
     x                                  NG2CHK,
     x                                  DAalpE,DAbetE,DAtotE,DBE,DP,
     x                                  GM2ICR,GM2sICR,
     x                                  FP,SP,
     x                                  E_P_OMG1, 
     x                                  E_P_OMG2,
     x                                  S_P_OMG1)

! In-core Construct QM-Particle Fock Matrix
! nelec: total number of electrons
! NAE: number regular electrons
! NBE: number special electrons
!
!======================================================================
      implicit none
! Input Variables
      logical LNEOHF,LADDEXCH
      integer nebf,npbf,npbf2,ngee,ng1,ng2,ng3,ng4
      integer SZG2ICR
      integer nelec,NAalpE,NAbetE,NBE
      integer NG2CHK
      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
      double precision DAtotE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM2ICR(SZG2ICR)    ! Regular OMG2 terms
      double precision GM2sICR(SZG2ICR)   ! Exchange OMG2 terms

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

      if (LADDEXCH) then
        call RXCUHFne_FP_OMG2(NG2CHK,nebf,npbf,ng2,NAalpE,NBE,
     x                      DAalpE,DAtotE,DBE,DP,
     x                      GM2ICR,GM2sICR,FP,E_P_OMG2)
      else
        call RXCHFne_FP_OMG2(NG2CHK,nebf,npbf,ng2,NAalpE+NAbetE,NBE,
     x                     DAtotE,DBE,DP,GM2ICR,FP,E_P_OMG2)
      end if

      return
      end
!======================================================================
      subroutine RXCUHFne_IC_construct_FBE(LNEOHF,LADDEXCH,
     x                                   nelec,NAbetE,NAalpE,NBE,
     x                                   nebf,npbf,nebf2,ngee,
     x                                   ng1,ng2,ng3,ng4,
     x                                   SZG2ICR,
     x                                   NG2CHK,
     x                                   DAalpE,DAbetE,DAtotE,DBE,DP,
     x                                   GM2ICR,GM2sICR,
     x                                   FBE,SBE,
     x                                   E_BE_OMG1, 
     x                                   E_BE_OMG2,
     x                                   S_BE_OMG1)

! In-core Construct special electronic Fock Matrix
! nelec: total number of electrons
! NAE: number regular electrons
! NBE: number special electrons
!
!======================================================================
      implicit none
! Input Variables
      logical LNEOHF,LADDEXCH
      integer nebf,npbf,nebf2,ngee,ng1,ng2,ng3,ng4
      integer SZG2ICR
      integer nelec,NAalpE,NAbetE,NBE
      integer NG2CHK
      double precision DAalpE(nebf,nebf)
      double precision DAbetE(nebf,nebf)
      double precision DAtotE(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GM2ICR(SZG2ICR)    ! Regular OMG2 terms
      double precision GM2sICR(SZG2ICR)   ! Exchange OMG2 terms

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

      if (LADDEXCH) then
        call RXCUHFne_FBE_OMG2(NG2CHK,nebf,npbf,ng2,NAalpE,NBE,
     x                       DAalpE,DAtotE,DBE,DP,
     x                       GM2ICR,GM2sICR,FBE,E_BE_OMG2)
      else
        call RXCHFne_FBE_OMG2(NG2CHK,nebf,npbf,ng2,NAalpE+NAbetE,NBE,
     x                      DAtotE,DBE,DP,GM2ICR,FBE,E_BE_OMG2)
      end if

      return
      end
