C======================================================================
      subroutine RXCHFmult_fock_int(LCMF,nelec,NAE,NBE,
     x                              nebf,nebf2,npbf,npbf2,
     x                              ng1,ng2,ng3,ng4,
     x                              SZG2ICR,SZG3ICR,SZG4ICR,
     x                              NG2CHK,NG3CHK,NG4CHK,
     x                              DAE,DBE,DP,
     x                              GM2ICR,GM3ICR,GM4ICR,
     x                              S_total,XSBE,XSP,
     x                              FP,FAE,FBE, 
     x                              E_OMG2,E_OMG3,E_OMG4,
     x                              E_total)

C Interaction energy and resulting contributions to Fock matrices
C======================================================================
      implicit none

C Input variables
      logical           LCMF
      integer           nebf,npbf,npbf2,nebf2
      integer           ng1,ng2,ng3,ng4
      integer           SZG2ICR,SZG3ICR,SZG4ICR
      integer           nelec,NAE,NBE
      integer           NG2CHK,NG3CHK,NG4CHK
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM2ICR(SZG2ICR)
      double precision  GM3ICR(SZG3ICR)
      double precision  GM4ICR(SZG4ICR)
      double precision  S_total           !  Read in overlap contributions
      double precision  XSBE(nebf,nebf)   !  to Fock matrices obtained
      double precision  XSP(npbf,npbf)    !  from previous XCHF calculation

C Output variables
      double precision  FP(npbf,npbf)
      double precision  FAE(nebf,nebf)
      double precision  FBE(nebf,nebf)
      double precision  E_total
      double precision  E_OMG1
      double precision  E_OMG2
      double precision  E_OMG3
      double precision  E_OMG4

C Local variables
      logical           rxchfdbg
      integer           nebflt,npbflt,i,j
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebf,nebf)
      double precision  XFP(npbf,npbf)

      double precision  E_AE_OMG2
      double precision  E_AE_OMG3
      double precision  E_AE_OMG4

      double precision  E_BE_OMG2
      double precision  E_BE_OMG3
      double precision  E_BE_OMG4

      double precision  E_P_OMG2
      double precision  E_P_OMG3
      double precision  E_P_OMG4

      double precision  Psi_HOMG_Psi

      double precision  zero,one
      parameter(zero=0.0d+00,one=1.0d+00)

      rxchfdbg=.false.
      if(LCMF) rxchfdbg=.true. 

C Initialize
      FAE=zero 
      FBE=zero 
      FP=zero 

      XFAE=zero 
      XFP=zero 
      XFBE=zero 

      E_AE_OMG2 = zero 
      E_AE_OMG3 = zero 
      E_AE_OMG4 = zero 

      E_P_OMG2 = zero 
      E_P_OMG3 = zero 
      E_P_OMG4 = zero 

      E_BE_OMG2 = zero 
      E_BE_OMG3 = zero 
      E_BE_OMG4 = zero 

C Construct Electronic Fock Matrix (regular electrons)
      call RXCHFmult_IC_construct_FAE(nelec,NAE,NBE,
     x                                nebf,npbf,nebf2,
     x                                ng1,ng2,ng3,ng4,
     x                                SZG2ICR,SZG3ICR,SZG4ICR,
     x                                NG2CHK,NG3CHK,NG4CHK,
     x                                DAE,DBE,DP,
     x                                GM2ICR,GM3ICR,GM4ICR,
     x                                XFAE,
     x                                E_AE_OMG2,
     x                                E_AE_OMG3,
     x                                E_AE_OMG4)

C Construct QM Particle Fock Matrix
      call RXCHFmult_IC_construct_FP(nelec,NAE,NBE,
     x                               nebf,npbf,npbf2,
     x                               ng1,ng2,ng3,ng4,
     x                               SZG2ICR,SZG3ICR,SZG4ICR,
     x                               NG2CHK,NG3CHK,NG4CHK,
     x                               DAE,DBE,DP,
     x                               GM2ICR,GM3ICR,GM4ICR,
     x                               XFP,
     x                               E_P_OMG2,
     x                               E_P_OMG3, 
     x                               E_P_OMG4)

C Construct Electronic Fock Matrix (special electrons)
      call RXCHFmult_IC_construct_FBE(nelec,NAE,NBE,
     x                                nebf,npbf,nebf2,
     x                                ng1,ng2,ng3,ng4,
     x                                SZG2ICR,SZG3ICR,SZG4ICR,
     x                                NG2CHK,NG3CHK,NG4CHK,
     x                                DAE,DBE,DP,
     x                                GM2ICR,GM3ICR,GM4ICR,
     x                                XFBE,
     x                                E_BE_OMG2,
     x                                E_BE_OMG3, 
     x                                E_BE_OMG4)

      if (rxchfdbg) then
       write(*,*)
       write(*,*) "Read in S_total:", S_total
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
      end if

      Psi_HOMG_Psi = E_P_OMG2 + E_P_OMG3 + E_P_OMG4

      E_OMG2  = E_P_OMG2 / S_total
      E_OMG3  = E_P_OMG3 / S_total
      E_OMG4  = E_P_OMG4 / S_total

      E_total = E_OMG2 + E_OMG3 + E_OMG4

C Correct Regular Electron Fock Matrix
      call RXCHFmult_FAE_correction(nebf,S_total,XFAE,FAE)

C Correct QM Particle Fock Matrix
      call RXCHF_Fock_correction2(npbf,Psi_HOMG_Psi,S_total,
     x                            XFP,XSP,FP)

C Correct Special Electron Fock Matrix
      call RXCHF_Fock_correction2(nebf,Psi_HOMG_Psi,S_total,
     x                            XFBE,XSBE,FBE)

C Fock testing
      if (LCMF) then
C        call RXCHF_Fock_testing1(nebf,FAE,DAE,E_ecore,E_GAMee,E_OMG2)

        call RXCHF_Fock_testing2(npbf,FP,DP)

        call RXCHF_Fock_testing2(nebf,FBE,DBE)

        call UFM_sym_check(nebf,npbf,FAE,FBE,FP)
      end if

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

      return
      end
C======================================================================
      subroutine RXCHFmult_FAE_correction(nbf,S_total,XF,F)

C Corrects FAE (with overlap terms)
C======================================================================

      implicit none

C Input variables
      integer           nbf
      double precision  S_total
      double precision  XF(nbf,nbf)

C Output variables
      double precision  F(nbf,nbf)

C Local variables
      integer           i,j
      double precision  coeff
      double precision  zero,one
      parameter(zero=0.0d+00,one=1.0d+00)

C Initialize
      F=zero
      coeff=one/S_total

C Fock Matrix Correction
      do i=1,nbf
        do j=1,nbf
          F(j,i) = coeff1*XF(j,i)
        end do
      end do

      return
      end
C======================================================================
      subroutine RXCHFmult_IC_construct_FAE(nelec,NAE,NBE,
     x                                      nebf,npbf,nebf2,
     x                                      ng1,ng2,ng3,ng4,
     x                                      SZG2ICR,SZG3ICR,SZG4ICR,
     x                                      NG2CHK,NG3CHK,NG4CHK,
     x                                      DAE,DBE,DP,
     x                                      GM2ICR,GM3ICR,GM4ICR,
     x                                      FAE,
     x                                      E_AE_OMG2,
     x                                      E_AE_OMG3,
     x                                      E_AE_OMG4)

C In-core Construct Regular Electronic Fock Matrix
C nelec: total number of electrons
C NAE: number regular electrons
C NBE: number special electrons
C======================================================================
      implicit none

C Input variables
      integer           nebf,npbf,nebf2
      integer           ng1,ng2,ng3,ng4
      integer           SZG2ICR,SZG3ICR,SZG4ICR
      integer           nelec,NAE,NBE
      integer           NG2CHK,NG3CHK,NG4CHK
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM2ICR(SZG2ICR)
      double precision  GM3ICR(SZG3ICR)
      double precision  GM4ICR(SZG4ICR)

C Output variables
      double precision  FAE(nebf,nebf)
      double precision  E_AE_OMG2
      double precision  E_AE_OMG3
      double precision  E_AE_OMG4

C Local Variables
      double precision  zero
      parameter(zero=0.0d+00)

C Initialize
      FAE=zero
      E_AE_OMG2  = zero
      E_AE_OMG3  = zero
      E_AE_OMG4  = zero

      call RXCHFmult_FAE_OMG2(NG2CHK,nebf,npbf,ng2,
     x                        DAE,DBE,DP,GM2ICR,
     x                        FAE,E_AE_OMG2)

      if (NBE.gt.1) then

         call RXCHFmult_FAE_OMG3(NG3CHK,nebf,npbf,ng3,
     x                           DAE,DBE,DP,GM3ICR,
     x                           FAE,E_AE_OMG3)

         if (NBE.gt.2) then

            call RXCHFmult_FAE_OMG4(NG4CHK,nebf,npbf,ng4,
     x                              DAE,DBE,DP,GM4ICR,
     x                              FAE,E_AE_OMG4)

         end if

      end if

      return
      end
C======================================================================
      subroutine RXCHFmult_IC_construct_FP(nelec,NAE,NBE,
     x                                     nebf,npbf,npbf2,
     x                                     ng1,ng2,ng3,ng4,
     x                                     SZG2ICR,SZG3ICR,SZG4ICR,
     x                                     NG2CHK,NG3CHK,NG4CHK,
     x                                     DAE,DBE,DP,
     x                                     GM2ICR,GM3ICR,GM4ICR,
     x                                     FP,
     x                                     E_P_OMG2,
     x                                     E_P_OMG3, 
     x                                     E_P_OMG4)

C In-core Construct QM-Particle Fock Matrix
C nelec: total number of electrons
C NAE: number regular electrons
C NBE: number special electrons
C======================================================================
      implicit none

C Input variables
      integer           nebf,npbf,npbf2
      integer           ng1,ng2,ng3,ng4
      integer           SZG2ICR,SZG3ICR,SZG4ICR
      integer           nelec,NAE,NBE
      integer           NG2CHK,NG3CHK,NG4CHK
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM2ICR(SZG2ICR)
      double precision  GM3ICR(SZG3ICR)
      double precision  GM4ICR(SZG4ICR)

C Output variables
      double precision  FP(npbf,npbf)
      double precision  E_P_OMG2
      double precision  E_P_OMG3
      double precision  E_P_OMG4

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)

C Initialize
      FP=zero    
      E_P_OMG2  = zero 
      E_P_OMG3  = zero 
      E_P_OMG4  = zero 
 
      call RXCHFmult_FP_OMG2(NG2CHK,nebf,npbf,ng2,
     x                       DAE,DBE,DP,GM2ICR,
     x                       FP,E_P_OMG2)

      if (NBE.gt.1) then

         call RXCHFmult_FP_OMG3(NG3CHK,nebf,npbf,ng3,
     x                          DAE,DBE,DP,GM3ICR,
     x                          FP,E_P_OMG3)

         if (NBE.gt.2) then

            call RXCHFmult_FP_OMG4(NG4CHK,nebf,npbf,ng4,
     x                             DAE,DBE,DP,GM4ICR,
     x                             FP,E_P_OMG4)

         end if

      end if

      return
      end
C======================================================================
      subroutine RXCHFmult_IC_construct_FBE(nelec,NAE,NBE,
     x                                      nebf,npbf,nebf2,
     x                                      ng1,ng2,ng3,ng4,
     x                                      SZG2ICR,SZG3ICR,SZG4ICR,
     x                                      NG2CHK,NG3CHK,NG4CHK,
     x                                      DAE,DBE,DP,
     x                                      GM2ICR,GM3ICR,GM4ICR,
     x                                      FBE,
     x                                      E_BE_OMG2,
     x                                      E_BE_OMG3, 
     x                                      E_BE_OMG4)

C In-core Construct special electronic Fock Matrix
C nelec: total number of electrons
C NAE: number regular electrons
C NBE: number special electrons
C======================================================================
      implicit none

C Input variables
      integer           nebf,npbf,nebf2
      integer           ng1,ng2,ng3,ng4
      integer           SZG2ICR,SZG3ICR,SZG4ICR
      integer           nelec,NAE,NBE
      integer           NG2CHK,NG3CHK,NG4CHK
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM2ICR(SZG2ICR)
      double precision  GM3ICR(SZG3ICR)
      double precision  GM4ICR(SZG4ICR)

C Output variables
      double precision  FBE(nebf,nebf)
      double precision  E_BE_OMG2
      double precision  E_BE_OMG3 
      double precision  E_BE_OMG4

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)

C Initialize
      FBE=zero    
      E_BE_OMG2  = zero 
      E_BE_OMG3  = zero 
      E_BE_OMG4  = zero 
 
      call RXCHFmult_FBE_OMG2(NG2CHK,nebf,npbf,ng2,
     x                        DAE,DBE,DP,GM2ICR,GM2sICR,
     x                        FBE,E_BE_OMG2)

      if (NBE.gt.1) then

         call RXCHFmult_FBE_OMG3(NG3CHK,nebf,npbf,ng3,
     x                           DAE,DBE,DP,GM3ICR,
     x                           FBE,E_BE_OMG3)

         if (NBE.gt.2) then

            call RXCHFmult_FBE_OMG4(NG4CHK,nebf,npbf,ng4,
     x                              DAE,DBE,DP,GM4ICR,
     x                              FBE,E_BE_OMG4)

         end if

      end if

      return
      end
