C======================================================================
      subroutine RXCHFmult_fock_int(LCMF,nelec,NAE,NBE,
     x                              nebf,nebfBE,npbf,
     x                              SZG2ICR,SZG3ICR,SZG4ICR,
     x                              NG2CHK,NG3CHK,NG4CHK,
     x                              DAE,DBE,DP,
     x                              GM2ICR,GM3ICR,GM4ICR,
     x                              S_total,S_OMG2,XSBE,XSP,
     x                              FP,FAE,FBE, 
     x                              E_OMG2,E_OMG3,E_OMG4,
     x                              E_total)

C Interaction energy and resulting contributions to Fock matrices
C======================================================================
      implicit none

C Input variables
      logical           LCMF
      integer           nebf,nebfBE,npbf
      integer           ng2,ng3,ng4
      integer           SZG2ICR,SZG3ICR,SZG4ICR
      integer           nelec,NAE,NBE
      integer           NG2CHK,NG3CHK,NG4CHK
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebfBE,nebfBE)
      double precision  DP(npbf,npbf)
      double precision  GM2ICR(SZG2ICR)
      double precision  GM3ICR(SZG3ICR)
      double precision  GM4ICR(SZG4ICR)
      double precision  S_total             !  Read in overlap contributions
      double precision  S_OMG2              !  to Fock matrices obtained       
      double precision  XSBE(nebfBE,nebfBE) !  from previous XCHF calculation
      double precision  XSP(npbf,npbf)      !

C Output variables
      double precision  FP(npbf,npbf)
      double precision  FAE(nebf,nebf)
      double precision  FBE(nebfBE,nebfBE)
      double precision  E_total
      double precision  E_OMG1
      double precision  E_OMG2
      double precision  E_OMG3
      double precision  E_OMG4

C Local variables
      logical           rxchfdbg
      integer           nebflt,nebfBElt,npbflt
      integer           i,j
      double precision  XFAE(nebf,nebf)
      double precision  XFBE(nebfBE,nebfBE)
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

      ng2=SZG2ICR
      ng3=SZG3ICR
      ng4=SZG4ICR

C Construct Electronic Fock Matrix (regular electrons)
      call RXCHFmult_IC_construct_FAE(nelec,NAE,NBE,
     x                                nebf,nebfBE,npbf,
     x                                ng2,ng3,ng4,
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
     x                               nebf,nebfBE,npbf,
     x                               ng2,ng3,ng4,
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
     x                                nebf,nebfBE,npbf,
     x                                ng2,ng3,ng4,
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
       write(*,*) "S_OMG2  (read) = ", S_OMG2
       write(*,*) "S_total (read) = ", S_total
       write(*,*)
       write(*,*) 'E_P_OMG2       =',E_P_OMG2 
       write(*,*) 'E_AE_OMG2      =',E_AE_OMG2 
       write(*,*) 'E_BE_OMG2      =',E_BE_OMG2 
       write(*,*)
       write(*,*) 'E_P_OMG3       =',E_P_OMG3 
       write(*,*) 'E_AE_OMG3      =',E_AE_OMG3 
       write(*,*) 'E_BE_OMG3      =',E_BE_OMG3 
       write(*,*)
       write(*,*) 'E_P_OMG4       =',E_P_OMG4 
       write(*,*) 'E_AE_OMG4      =',E_AE_OMG4 
       write(*,*) 'E_BE_OMG4      =',E_BE_OMG4 
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
      call RXCHF_Fock_correction2(nebfBE,Psi_HOMG_Psi,S_total,
     x                            XFBE,XSBE,FBE)

C Fock testing
      if (LCMF) then
        call RXCHFmult_Fock_testing(nebf,nebfBE,npbf,
     x                              FAE,FBE,FP,DAE,DBE,DP,
     x                              E_total,S_total,
     x                              E_OMG3,E_OMG4,S_OMG2)
        call UFM_sym_check2(nebf,FAE)
        call UFM_sym_check2(nebfBE,FBE)
        call UFM_sym_check2(npbf,FP)
      end if

      if (rxchfdbg) then
        nebflt=nebf*(nebf+1)/2
        nebfBElt=nebfBE*(nebfBE+1)/2
        npbflt=npbf*(npbf+1)/2
        write(*,*) "FAE int:"
        call prt_lower_triangle(nebf,nebflt,FAE)
        write(*,*)
        write(*,*) "FBE int:"
        call prt_lower_triangle(nebfBE,nebfBElt,FBE)
        write(*,*)
        write(*,*) "FP int:"
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
          F(j,i) = coeff*XF(j,i)
        end do
      end do

      return
      end
!======================================================================
      subroutine RXCHFmult_Fock_testing(nebf,nebfBE,npbf,
     x                                  FAE,FBE,FP,DAE,DBE,DP,
     x                                  E_total,S_total,
     x                                  E_OMG3,E_OMG4,S_OMG2)
C Tests Fock matrices by summing up with density matrices
!======================================================================
      implicit none

! Input variables
      integer nebf,nebfBE,npbf
      double precision FAE(nebf,nebf)
      double precision FBE(nebfBE,nebfBE)
      double precision FP(npbf,npbf)
      double precision DAE(nebf,nebf)
      double precision DBE(nebfBE,nebfBE)
      double precision DP(npbf,npbf)
      double precision E_total,E_OMG3,E_OMG4 ! On entry already div by S_total
      double precision S_total,S_OMG2

! Local variables
      integer i,j
      double precision fdsum,ans
      double precision zero,two
      parameter(zero=0.0d+00,two=2.0d+00)

      fdsum=zero
      do i=1,npbf
        do j=1,npbf
          fdsum=fdsum+DP(j,i)*FP(j,i)
        end do
      end do
      ans=zero
      write(*,*) "Proton Fock matrix test:"
      write(*,*) fdsum,ans

      fdsum=zero
      do i=1,nebf
        do j=1,nebf
          fdsum=fdsum+DAE(j,i)*FAE(j,i)
        end do
      end do
      ans=E_total
      write(*,*) "Regular electron Fock matrix test:"
      write(*,*) fdsum,ans

      fdsum=zero
      do i=1,nebfBE
        do j=1,nebfBE
          fdsum=fdsum+DBE(j,i)*FBE(j,i)
        end do
      end do
      ans=E_OMG3+two*E_OMG4-S_OMG2*E_total/S_total
      write(*,*) "Special electron Fock matrix test:"
      write(*,*) fdsum,ans

      return
      end
C======================================================================
      subroutine RXCHFmult_IC_construct_FAE(nelec,NAE,NBE,
     x                                      nebf,nebfBE,npbf,
     x                                      ng2,ng3,ng4,
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
      integer           nebf,nebfBE,npbf
      integer           ng2,ng3,ng4
      integer           SZG2ICR,SZG3ICR,SZG4ICR
      integer           nelec,NAE,NBE
      integer           NG2CHK,NG3CHK,NG4CHK
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebfBE,nebfBE)
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

      call RXCHFmult_FAE_OMG2(NG2CHK,nebf,nebfBE,npbf,ng2,
     x                        DAE,DBE,DP,GM2ICR,
     x                        FAE,E_AE_OMG2)

      if (NBE.gt.1) then

         call RXCHFmult_FAE_OMG3(NG3CHK,nebf,nebfBE,npbf,ng3,
     x                           DAE,DBE,DP,GM3ICR,
     x                           FAE,E_AE_OMG3)

         if (NBE.gt.2) then

            call RXCHFmult_FAE_OMG4(NG4CHK,nebf,nebfBE,npbf,ng4,
     x                              DAE,DBE,DP,GM4ICR,
     x                              FAE,E_AE_OMG4)

         end if

      end if

      return
      end
C======================================================================
      subroutine RXCHFmult_IC_construct_FP(nelec,NAE,NBE,
     x                                     nebf,nebfBE,npbf,
     x                                     ng2,ng3,ng4,
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
      integer           nebf,nebfBE,npbf
      integer           ng2,ng3,ng4
      integer           SZG2ICR,SZG3ICR,SZG4ICR
      integer           nelec,NAE,NBE
      integer           NG2CHK,NG3CHK,NG4CHK
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebfBE,nebfBE)
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
 
      call RXCHFmult_FP_OMG2(NG2CHK,nebf,nebfBE,npbf,ng2,
     x                       DAE,DBE,DP,GM2ICR,
     x                       FP,E_P_OMG2)

      if (NBE.gt.1) then

         call RXCHFmult_FP_OMG3(NG3CHK,nebf,nebfBE,npbf,ng3,
     x                          DAE,DBE,DP,GM3ICR,
     x                          FP,E_P_OMG3)

         if (NBE.gt.2) then

            call RXCHFmult_FP_OMG4(NG4CHK,nebf,nebfBE,npbf,ng4,
     x                             DAE,DBE,DP,GM4ICR,
     x                             FP,E_P_OMG4)

         end if

      end if

      return
      end
C======================================================================
      subroutine RXCHFmult_IC_construct_FBE(nelec,NAE,NBE,
     x                                      nebf,nebfBE,npbf,
     x                                      ng2,ng3,ng4,
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
      integer           nebf,nebfBE,npbf
      integer           ng2,ng3,ng4
      integer           SZG2ICR,SZG3ICR,SZG4ICR
      integer           nelec,NAE,NBE
      integer           NG2CHK,NG3CHK,NG4CHK
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebfBE,nebfBE)
      double precision  DP(npbf,npbf)
      double precision  GM2ICR(SZG2ICR)
      double precision  GM3ICR(SZG3ICR)
      double precision  GM4ICR(SZG4ICR)

C Output variables
      double precision  FBE(nebfBE,nebfBE)
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
 
      call RXCHFmult_FBE_OMG2(NG2CHK,nebf,nebfBE,npbf,ng2,
     x                        DAE,DBE,DP,GM2ICR,
     x                        FBE,E_BE_OMG2)

      if (NBE.gt.1) then

         call RXCHFmult_FBE_OMG3(NG3CHK,nebf,nebfBE,npbf,ng3,
     x                           DAE,DBE,DP,GM3ICR,
     x                           FBE,E_BE_OMG3)

         if (NBE.gt.2) then

            call RXCHFmult_FBE_OMG4(NG4CHK,nebf,nebfBE,npbf,ng4,
     x                              DAE,DBE,DP,GM4ICR,
     x                              FBE,E_BE_OMG4)

         end if

      end if

      return
      end
