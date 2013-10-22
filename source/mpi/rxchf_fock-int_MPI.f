C======================================================================
      subroutine RXCHF_fock_int_MPI(nproc,rank,
     x                              LCMF,LADDEXCH,nelec,NAE,NBE,
     x                              nebf,nebf2,npbf,npbf2,
     x                              ng1,ng2,ng3,ng4,
     x                              SZG2ICR,SZG3ICR,SZG4ICR,
     x                              SZG2exICR,SZG3exICR,
     x                              NG2CHK,NG3CHK,NG4CHK,
     x                              DAE,DBE,DP,
     x                              GM2ICR,GM3ICR,GM4ICR,
     x                              GM2exICR,GM3ex1ICR,GM3ex2ICR,
     x                              S_total,S_OMG2,XSBE,XSP,
     x                              FP,FAE,FBE, 
     x                              E_OMG2,E_OMG3,E_OMG4,
     x                              E_total)

C Interaction energy and resulting contributions to Fock matrices
C Loop over integrals, constructing Fock matrices simultaneously
C======================================================================
      implicit none

C Input variables
      integer           nproc,rank
      logical           LCMF
      logical           LADDEXCH
      integer           nebf,npbf,npbf2,nebf2
      integer           ng1,ng2,ng3,ng4
      integer           SZG2ICR,SZG3ICR,SZG4ICR
      integer           SZG2exICR,SZG3exICR
      integer           nelec,NAE,NBE
      integer           NG2CHK,NG3CHK,NG4CHK
      double precision  DAE(nebf,nebf)
      double precision  DBE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GM2ICR(SZG2ICR)
      double precision  GM3ICR(SZG3ICR)
      double precision  GM4ICR(SZG4ICR)
      double precision  GM2exICR(SZG2exICR)
      double precision  GM3ex1ICR(SZG3exICR)
      double precision  GM3ex2ICR(SZG3exICR)
      double precision  S_total           !  Read in overlap contributions
      double precision  S_OMG2            !  to Fock matrices obtained       
      double precision  XSBE(nebf,nebf)   !  from previous XCHF calculation
      double precision  XSP(npbf,npbf)    !

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

      double precision  Psi_HOMG_Psi

      integer*4 ierr

      double precision  zero
      parameter(zero=0.0d+00)

      rxchfdbg=.false.
      if(LCMF) rxchfdbg=.true. 

C Initialize
      FAE=zero 
      FBE=zero 
      FP=zero 

      XFAE=zero 
      XFP=zero 
      XFBE=zero 

      E_OMG2 = zero 
      E_OMG3 = zero 
      E_OMG4 = zero 

C Include GAM2 contributions
      if (LADDEXCH) then
       call RXCHF_OMG2ex_MPI(nproc,rank,
     x                       NG2CHK,nebf,npbf,
     x                       ng2,SZG2ICR,
     x                       DAE,DBE,DP,
     x                       GM2ICR,GM2exICR,
     x                       XFAE,XFBE,XFP,E_OMG2)
      else
       call RXCHF_OMG2_MPI(nproc,rank,
     x                     NG2CHK,nebf,npbf,
     x                     ng2,SZG2ICR,
     x                     DAE,DBE,DP,
     x                     GM2ICR,
     x                     XFAE,XFBE,XFP,E_OMG2)
      end if

C Include GAM3 contributions
      if (NBE.gt.1) then

       if (LADDEXCH) then
        call RXCHF_OMG3ex_MPI(nproc,rank,
     x                        NG3CHK,nebf,npbf,
     x                        ng3,SZG3ICR,
     x                        DAE,DBE,DP,
     x                        GM3ICR,GM3ex1ICR,GM3ex2ICR,
     x                        XFAE,XFBE,XFP,E_OMG3)
       else
        call RXCHF_OMG3_MPI(nproc,rank,
     x                      NG3CHK,nebf,npbf,
     x                      ng3,SZG3ICR,
     x                      DAE,DBE,DP,
     x                      GM3ICR,
     x                      XFAE,XFBE,XFP,E_OMG3)
       end if

C Include GAM4 contributions
       if (NBE.gt.2) then

        call RXCHF_OMG4_MPI(nproc,rank,
     x                      NG4CHK,nebf,npbf,
     x                      ng4,SZG4ICR,
     x                      DAE,DBE,DP,
     x                      GM4ICR,
     x                      XFAE,XFBE,XFP,E_OMG4)

       end if

      end if

      if ((rxchfdbg).and.(rank.eq.0)) then
       write(*,*)
       write(*,*) "S_OMG2  (read) = ", S_OMG2
       write(*,*) "S_total (read) = ", S_total
       write(*,*)
       write(*,*) 'E_OMG2         =',E_OMG2 
       write(*,*)
       write(*,*) 'E_OMG3         =',E_OMG3 
       write(*,*)
       write(*,*) 'E_OMG4         =',E_OMG4 
       write(*,*)
      end if

      Psi_HOMG_Psi = E_OMG2 + E_OMG3 + E_OMG4

      E_OMG2  = E_OMG2 / S_total
      E_OMG3  = E_OMG3 / S_total
      E_OMG4  = E_OMG4 / S_total

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
      if ((LCMF).and.(rank.eq.0)) then
        call RXCHFmult_Fock_testing(nebf,npbf,FAE,FBE,FP,DAE,DBE,DP,
     x                              E_total,S_total,
     x                              E_OMG3,E_OMG4,S_OMG2)
        call UFM_sym_check(nebf,npbf,FAE,FBE,FP)
      end if

      if ((rxchfdbg).and.(rank.eq.0)) then
        nebflt=nebf*(nebf+1)/2
        npbflt=npbf*(npbf+1)/2
        write(*,*) "FAE int:"
        call prt_lower_triangle(nebf,nebflt,FAE)
        write(*,*)
        write(*,*) "FBE int:"
        call prt_lower_triangle(nebf,nebflt,FBE)
        write(*,*)
        write(*,*) "FP int:"
        call prt_lower_triangle(npbf,npbflt,FP)
      end if

      return
      end

