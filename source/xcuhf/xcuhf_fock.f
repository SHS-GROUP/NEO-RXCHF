!======================================================================
      subroutine UHF_FOCK(LCMF,LNEOHF,nelec,NAE,NBE,
     x                    nebf,nebf2,npbf,npbf2,ngee,
     x                    ng1,ng2,ng3,ng4,
     x                    SZG2ICR,SZG3ICR,SZG4ICR,
     x                    NG2CHK,NG3CHK,NG4CHK,
     x                    DAE,DBE,DP,
     x                    GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                    GM2ICR,GM2sICR,GM3ICR,GM4ICR,
     x                    FP,FAE,FBE, 
     x                    E_pcore,
     x                    E_ecore,
     x                    E_GAMep,
     x                    E_GAMee,
     x                    E_OMG1, 
     x                    E_OMG2,
     x                    E_OMG3,
     x                    E_OMG4,
     x                    S_OMG1,
     x                    S_OMG2,
     x                    E_UHF,
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
      double precision DETOT(nebf,nebf)
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
      double precision GM2sICR(SZG2ICR)
      double precision GM3ICR(SZG3ICR)
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

      double precision XSAE1(nebf,nebf)
      double precision XSAE2(nebf,nebf)
      double precision XSBE1(nebf,nebf)
      double precision XSBE2(nebf,nebf)
      double precision XSP1(npbf,npbf)
      double precision XSP2(npbf,npbf)

!     double precision E_AE
!     double precision OVLP_AE
      double precision E_AE_ecore
      double precision E_AE_GAMep
      double precision E_AE_GAMee
      double precision E_AE_OMG1
      double precision E_AE_OMG2
      double precision E_AE_OMG3
      double precision E_AE_OMG4
      double precision S_AE_OMG1
      double precision S_AE_OMG2

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
      double precision E_P_OMG3
      double precision E_P_OMG4
      double precision S_P_OMG1
      double precision S_P_OMG2

      double precision Psi_H_Psi

      double precision zero
      parameter(zero=0.0d+00)
      double precision one
      parameter(one=1.0d+00)

      integer npbfLT,nebfLT
      double precision E_check
!     double precision CHKfockp(npbf,npbf)
      double precision CHKS(npbf,npbf)

! Initialize
       FAE=zero 
       FBE=zero 
       FP=zero 
 
!      SAE=zero 
!      SBE=zero 
!      SP=zero 
 
       XFAE=zero 
       XFBE=zero 
       XFP=zero 
 
       XSAE1=zero 
       XSAE2=zero 
       XSBE1=zero 
       XSBE2=zero 
       XSP1=zero 
       XSP2=zero 
 
!      E_AE = zero  
!      OVLP_AE = zero  
       E_AE_ecore = zero  
       E_AE_GAMep = zero 
       E_AE_GAMee = zero 
       E_AE_OMG1 = zero 
       E_AE_OMG2 = zero 
       E_AE_OMG3 = zero 
       E_AE_OMG4 = zero 
       S_AE_OMG1 = zero 
       S_AE_OMG2 = zero 

!      E_BE = zero 
!      OVLP_BE = zero 
       E_BE_ecore = zero 
       E_BE_GAMep = zero 
       E_BE_GAMee = zero 
       E_BE_OMG1 = zero 
       E_BE_OMG2 = zero 
       E_BE_OMG3 = zero 
       E_BE_OMG4 = zero 
       S_BE_OMG1 = zero 
       S_BE_OMG2 = zero 

!      E_P = zero  
!      OVLP_P  = zero 
       E_P_pcore = zero 
       E_P_GAMep = zero 
       E_P_OMG1 = zero 
       E_P_OMG2 = zero 
       E_P_OMG3 = zero 
       E_P_OMG4 = zero 
       S_P_OMG1 = zero 
       S_P_OMG2 = zero 

! Form total electronic density 
      DETOT=0.0d+00
      call add2fock(nebf,DAE,DETOT)
      call add2fock(nebf,DBE,DETOT)

!     write(*,*)
!     write(*,*)'IN: UHF_FOCK '
!     write(*,*)'SZG2ICR=',SZG2ICR
!     write(*,*)'SZG3ICR=',SZG3ICR
!     write(*,*)'SZG4ICR=',SZG4ICR
!     write(*,*)
!     write(*,*)
!     write(*,*)'IN: UHF_FOCK '
!     write(*,*)'Before call IC_construct_FP, GAM_EE='
!     write(*,*)GAM_ee
!     write(*,*)

! Construct QM Particle Fock Matrix
      call IC_construct_FP(LNEOHF,nelec,NAE,NBE,
     x                           nebf,npbf,npbf2,ngee,
     x                           ng1,ng2,ng3,ng4,
     x                           SZG2ICR,SZG3ICR,SZG4ICR,
     x                           NG2CHK,NG3CHK,NG4CHK,
     x                           DAE,DBE,DETOT,DP,
     x                           GAM_pcore,GAM_ep,
     x                           GM2ICR,GM2sICR,GM3ICR,GM4ICR,
     x                           XFP,XSP1,XSP2, 
     x                           E_P_pcore,
     x                           E_P_GAMep,
     x                           E_P_OMG1, 
     x                           E_P_OMG2,
     x                           E_P_OMG3,
     x                           E_P_OMG4,
     x                           S_P_OMG1,
     x                           S_P_OMG2)

! Construct Alpha Electronic Fock Matrix
      call IC_construct_FAE(LNEOHF,nelec,NAE,NBE,
     x                            nebf,npbf,nebf2,ngee,
     x                            ng1,ng2,ng3,ng4,
     x                            SZG2ICR,SZG3ICR,SZG4ICR,
     x                            NG2CHK,NG3CHK,NG4CHK,
     x                            DAE,DBE,DETOT,DP,
     x                            GAM_ecore,GAM_ep,GAM_ee,
     x                            GM2ICR,GM2sICR,GM3ICR,GM4ICR,
     x                            XFAE,XSAE1,XSAE2, 
     x                            E_AE_ecore,
     x                            E_AE_GAMep,
     x                            E_AE_GAMee,
     x                            E_AE_OMG1, 
     x                            E_AE_OMG2,
     x                            E_AE_OMG3,
     x                            E_AE_OMG4,
     x                            S_AE_OMG1,
     x                            S_AE_OMG2)

! Construct Beta Electronic Fock Matrix
      call IC_construct_FAE(LNEOHF,nelec,NBE,NAE,
     x                            nebf,npbf,nebf2,ngee,
     x                            ng1,ng2,ng3,ng4,
     x                            SZG2ICR,SZG3ICR,SZG4ICR,
     x                            NG2CHK,NG3CHK,NG4CHK,
     x                            DBE,DAE,DETOT,DP,
     x                            GAM_ecore,GAM_ep,GAM_ee,
     x                            GM2ICR,GM2sICR,GM3ICR,GM4ICR,
     x                            XFBE,XSBE1,XSBE2, 
     x                            E_BE_ecore,
     x                            E_BE_GAMep,
     x                            E_BE_GAMee,
     x                            E_BE_OMG1, 
     x                            E_BE_OMG2,
     x                            E_BE_OMG3,
     x                            E_BE_OMG4,
     x                            S_BE_OMG1,
     x                            S_BE_OMG2)

! Total Overlap
       S_OMG1  = S_P_OMG1 
       S_OMG2  = S_P_OMG2 
       S_total = one+S_OMG1+S_OMG2   

!-----------ON-THE-FLY-CODE-TESTING------------------------------------(
      if(LCMF) then

         write(*,*)
         write(*,*)'IN: UHF_FOCK'
         write(*,*)
         write(*,*) 'E_ecore=',(E_AE_ecore+E_BE_ecore)
         write(*,*) 'E_AE_ecore=',E_AE_ecore
         write(*,*) 'E_BE_ecore=',E_BE_ecore
         write(*,*)
         write(*,*) 'E_GAMee=',(E_AE_GAMee+E_BE_GAMee)
         write(*,*) 'E_AE_GAMee=',E_AE_GAMee
         write(*,*) 'E_BE_GAMee=',E_BE_GAMee
         write(*,*)
         write(*,*) 'E_P_GAMep =',E_P_GAMep  
         write(*,*) 'E_E_GAMep =',(E_AE_GAMep+E_BE_GAMep)
         write(*,*) 'E_AE_GAMep=',E_AE_GAMep
         write(*,*) 'E_BE_GAMep=',E_BE_GAMep
         write(*,*)
         write(*,*) 'E_P_OMG1  =',E_P_OMG1  
         write(*,*) 'E_E_OMG1  =',(E_AE_OMG1+E_BE_OMG1)
         write(*,*) 'E_AE_OMG1 =',E_AE_OMG1 
         write(*,*) 'E_BE_OMG1 =',E_BE_OMG1 
         write(*,*)
         write(*,*) 'E_P_OMG2  =',E_P_OMG2 
         write(*,*) 'E_E_OMG2  =',(E_AE_OMG2+E_BE_OMG2)
         write(*,*) 'E_AE_OMG2 =',E_AE_OMG2 
         write(*,*) 'E_BE_OMG2 =',E_BE_OMG2 
         write(*,*)
         write(*,*) 'E_P_OMG3  =',E_P_OMG3 
         write(*,*) 'E_E_OMG3  =',(E_AE_OMG3+E_BE_OMG3)
         write(*,*) 'E_AE_OMG3 =',E_AE_OMG3 
         write(*,*) 'E_BE_OMG3 =',E_BE_OMG3 
         write(*,*)
         write(*,*) 'E_P_OMG4  =',E_P_OMG4 
         write(*,*) 'E_E_OMG4  =',(E_AE_OMG4+E_BE_OMG4)
         write(*,*) 'E_AE_OMG4 =',E_AE_OMG4 
         write(*,*) 'E_BE_OMG4 =',E_BE_OMG4 
         write(*,*)
         write(*,*) 'S_P_OMG1  =',S_P_OMG1 
         write(*,*) 'S_E_OMG1  =',(S_AE_OMG1+S_BE_OMG1)
         write(*,*) 'S_AE_OMG1 =',S_AE_OMG1 
         write(*,*) 'S_BE_OMG1 =',S_BE_OMG1 
         write(*,*)
         write(*,*) 'S_P_OMG2  =',S_P_OMG2 
         write(*,*) 'S_E_OMG2  =',(S_AE_OMG2+S_BE_OMG2)
         write(*,*) 'S_AE_OMG2 =',S_AE_OMG2 
         write(*,*) 'S_BE_OMG2 =',S_BE_OMG2 
         write(*,*)
         write(*,*)'E_E_Total=',( E_AE_ecore+E_BE_ecore
     x                       + E_P_pcore 
     x                       + E_AE_GAMee+E_BE_GAMee
     x                       + E_AE_GAMep+E_BE_GAMep
     x                       + E_AE_OMG1+E_BE_OMG1
     x                       + E_AE_OMG2+E_BE_OMG2
     x                       + E_AE_OMG3+E_BE_OMG3
     x                       + E_AE_OMG4+E_BE_OMG4 ) /
     x                       ( one
     x                       + S_AE_OMG1+S_BE_OMG1
     x                       + S_AE_OMG2+S_BE_OMG2 )
         write(*,*)
      end if ! endif for LCMF
!-----------ON-THE-FLY-CODE-TESTING------------------------------------)

! Energetics
       E_pcore = E_P_pcore / S_total
       E_ecore = (E_AE_ecore+E_BE_ecore) / S_total
       E_GAMep = E_P_GAMep / S_total  
       E_GAMee = (E_AE_GAMee+E_BE_GAMee) / S_total 
       E_OMG1  = E_P_OMG1  / S_total 
       E_OMG2  = E_P_OMG2  / S_total 
       E_OMG3  = E_P_OMG3  / S_total 
       E_OMG4  = E_P_OMG4  / S_total 

       E_UHF = E_pcore + E_ecore + E_GAMep + E_GAMee      

       E_total = E_UHF + E_OMG1 + E_OMG2 + E_OMG3 + E_OMG4

       Psi_H_Psi = E_P_pcore 
     x           + (E_AE_ecore+E_BE_ecore) 
     x           + (E_AE_GAMee+E_BE_GAMee)  
     x           + E_P_GAMep   
     x           + E_P_OMG1   
     x           + E_P_OMG2   
     x           + E_P_OMG3   
     x           + E_P_OMG4   

! Correct QM Particle Fock Matrix
      call UHF_Fock_correction(npbf,Psi_H_Psi,S_total,
     x                         XFP,XSP1,XSP2,FP)

! Correct Alpha Electronic Fock Matrix
      call UHF_Fock_correction(nebf,Psi_H_Psi,S_total,
     x                         XFAE,XSAE1,XSAE2,FAE)

! Correct Beta Electronic Fock Matrix
      call UHF_Fock_correction(nebf,Psi_H_Psi,S_total,
     x                         XFBE,XSBE1,XSBE2,FBE)

!-----------ON-THE-FLY-CODE-TESTING------------------------------------(
      if(LCMF) then

         npbfLT=npbf*(npbf+1)/2
         nebfLT=nebf*(nebf+1)/2

         E_check=E_P_pcore+E_P_GAMep+E_P_OMG1+E_P_OMG2+E_P_OMG3+E_P_OMG4

         CHKS=0.0d+00
         call add2fock(npbf,XSP1,CHKS)
         call add2fock(npbf,XSP2,CHKS)

         call CHK_my_fockp(npbf,npbfLT,E_check,S_total,DP,
     x                     XFP,CHKS)

!        write(*,*)
!        write(*,*)'CHECK ALPHA FOCK MATRIX:'
!        call CHECK_FAE(nebf,nebfLT,
!    x                  DAE,FAE,XSAE1,XSAE2,E_total,S_total,
!    x                  E_P_pcore,E_AE_ecore,E_AE_GAMep,E_AE_GAMee,
!    x                  E_AE_OMG1,E_AE_OMG2,E_AE_OMG3,E_AE_OMG4)

!        write(*,*)
!        write(*,*)'CHECK BETA FOCK MATRIX:'
!        call CHECK_FAE(nebf,nebfLT,
!    x                  DBE,FBE,XSBE1,XSBE2,S_total,
!    x                  E_P_pcore,E_BE_ecore,E_BE_GAMep,E_BE_GAMee,
!    x                  E_BE_OMG1,E_BE_OMG2,E_BE_OMG3,E_BE_OMG4)

!        write(*,*)
!        write(*,*)'CHECK QMP FOCK MATRIX:'
!        call CHECK_FP(npbf,npbfLT,
!    x                 DP,FP,XSP1,XSP2,S_total,
!    x                 E_P_pcore,E_AE_ecore+E_BE_ecore,
!    x                 E_P_GAMep,E_AE_GAMee+E_BE_GAMee,
!    x                 E_P_OMG1,E_P_OMG2,E_P_OMG3,E_P_OMG4)

         call UFM_sym_check(nebf,npbf,FAE,FBE,FP)
!        call NonGem_OVLAP_check(nebf,npbf,nelec,DE,DP)
C ARS( print out Fock matrices
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
C )
      end if
!-----------ON-THE-FLY-CODE-TESTING------------------------------------)


      return
      end
!=======================================================================
      subroutine CHECK_FAE(nebf,nebfLT,
     x                     DE,focke,se1,se2,S_total,
     x                     E_pcore,E_ecore,E_ep,E_ee,
     x                     E_gam1,E_gam2,E_gam3,E_gam4)
!=======================================================================
      implicit none

! Input Variables
      integer nebf
      integer nebfLT
      double precision S_total
      double precision E_pcore
      double precision E_ecore
      double precision E_ep
      double precision E_ee
      double precision E_gam1
      double precision E_gam2
      double precision E_gam3
      double precision E_gam4
      double precision DE(nebf,nebf)
      double precision focke(nebf,nebf)
      double precision se1(nebf,nebf)
      double precision se2(nebf,nebf)

! Local Variables
      integer ie,je,iLT
      double precision se(nebf,nebf)
      double precision SIGMAe_I
      double precision SIGMAe_II
      double precision LTfocke(nebfLT)
      double precision LTDE(nebfLT)
      double precision LTSE(nebfLT)
      double precision E_XCHF
      double precision FCKTRC_e
      double precision Check_e
      double precision one,two,three,four
      parameter(one=1.0d+00,two=2.0d+00,three=3.0d+00,four=4.0d+00)

      double precision TRACEP


      do ie=1,nebf
         do je=1,nebf
            se(ie,je) = se1(ie,je) + two * se2(ie,je)
         end do
      end do

      iLT=0
      do ie=1,nebf
         do je=1,ie

            iLT=iLT+1
            LTfocke(iLT)=focke(ie,je)
            LTDE(iLT)=DE(ie,je)
            LTSE(iLT)=se(ie,je)

         end do
      end do

      E_XCHF = E_ecore
     x       + E_pcore
     x       + E_ep
     x       + E_ee
     x       + E_gam1
     x       + E_gam2
     x       + E_gam3
     x       + E_gam4

      SIGMAe_I =         E_ecore
     x         +         E_ep
     x         + two *   E_ee
     x         +         E_gam1
     x         + two *   E_gam2
     x         + three * E_gam3
     x         + four *  E_gam4

      SIGMAe_II=TRACEP(LTDE,LTSE,nebf)
      FCKTRC_e=TRACEP(LTDE,LTfocke,nebf)

      Check_e= ( SIGMAe_I / S_total )
     x       - ( E_XCHF*SIGMAe_II / S_total**2 )

!     write(*,*)
!     write(*,*)'>>> FOCK CHECK <<<'
!     write(*,*)
!     write(*,*)'S_total   = ',S_total
!     write(*,*)'E_XCHF    = ',E_XCHF
!     write(*,*)'SIGMAe_I  =',SIGMAe_I
!     write(*,*)'SIGMAe_II =',SIGMAe_II
      write(*,*)'===================================='
      write(*,*)'FCKTRC_e  = ',FCKTRC_e
      write(*,*)'Check_e   = ',Check_e
      write(*,*)'===================================='
      write(*,*)


      return
      end
!=======================================================================
      subroutine CHECK_FP(npbf,npbfLT,
     x                    DP,fockp,sp1,sp2,S_total,
     x                    E_pcore,E_ecore,
     x                    E_ep,E_ee,
     x                    E_gam1,E_gam2,E_gam3,E_gam4)
!=======================================================================
      implicit none

! Input Variables
      integer npbf
      integer npbfLT
      double precision S_total
      double precision E_pcore
      double precision E_ecore
      double precision E_ep
      double precision E_ee
      double precision E_gam1
      double precision E_gam2
      double precision E_gam3
      double precision E_gam4
      double precision DP(npbf,npbf)
      double precision fockp(npbf,npbf)
      double precision sp1(npbf,npbf)
      double precision sp2(npbf,npbf)

! Local Variables
      integer ip,jp,iLT
      double precision sp(npbf,npbf)
      double precision SIGMAp_I
      double precision SIGMAp_II
      double precision LTfockp(npbfLT)
      double precision LTDP(npbfLT)
      double precision LTSP(npbfLT)
      double precision E_XCHF
      double precision FCKTRC_p
      double precision Check_p
      double precision one,two,three,four
      parameter(one=1.0d+00,two=2.0d+00,three=3.0d+00,four=4.0d+00)

      double precision TRACEP


      do ip=1,npbf
         do jp=1,npbf
            sp(ip,jp) = sp1(ip,jp) + sp2(ip,jp)
         end do
      end do

      iLT=0
      do ip=1,npbf
         do jp=1,ip

            iLT=iLT+1
            LTfockp(iLT)=fockp(ip,jp)
            LTDP(iLT)=DP(ip,jp)
            LTSP(iLT)=sp(ip,jp)

         end do
      end do

      E_XCHF = E_ecore
     x       + E_pcore
     x       + E_ep
     x       + E_ee
     x       + E_gam1
     x       + E_gam2
     x       + E_gam3
     x       + E_gam4


      FCKTRC_p=TRACEP(LTDP,LTfockp,npbf)

      Check_p= ( E_XCHF / S_total**2 )
     x       - ( E_ecore / S_total )
     x       - ( E_ee / S_total )


!     write(*,*)
!     write(*,*)'>>> FOCK CHECK <<<'
!     write(*,*)
!     write(*,*)'S_total   = ',S_total
!     write(*,*)'E_XCHF    = ',E_XCHF
!     write(*,*)'SIGMAe_I  =',SIGMAe_I
!     write(*,*)'SIGMAe_II =',SIGMAe_II
      write(*,*)'===================================='
      write(*,*)'FCKTRC_p  = ',FCKTRC_p
      write(*,*)'Check_p   = ',Check_p
      write(*,*)'===================================='
      write(*,*)


      return
      end
!=======================================================================
      subroutine UFM_sym_check(nebf,npbf,FAE,FBE,FP)

! Checks that the Fock matrices are symmetric
!=======================================================================
      implicit none

! Input Variables
      integer nebf,npbf
      double precision FAE(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision FP(npbf,npbf)

! Local variables
      logical LASYM
      logical LBSYM
      logical LPSYM
      integer i,j
      double precision val,tolerance
      parameter(tolerance=1.0d-10)

      LASYM=.true.
      LBSYM=.true.
      LPSYM=.true.

      do i=1,nebf
         do j=1,i
            val=FAE(i,j)-FAE(j,i)
            if(abs(val).gt.tolerance) then
              write(*,*)'ALPHA FOCK MATRIX IS NOT SYMMETRIC FOR IJ=',I,J
              write(*,*)'>>>> FM(IJ)=',FAE(i,j)
              write(*,*)'>>>> FM(JI)=',FAE(j,i)
              LASYM=.false.
            end if
         end do
      end do
      do i=1,nebf
         do j=1,i
            val=FBE(i,j)-FBE(j,i)
            if(abs(val).gt.tolerance) then
              write(*,*)'BETA FOCK MATRIX IS NOT SYMMETRIC FOR IJ=',I,J
              write(*,*)'>>>> FM(IJ)=',FBE(i,j)
              write(*,*)'>>>> FM(JI)=',FBE(j,i)
              LBSYM=.false.
            end if
         end do
      end do
      do i=1,npbf
         do j=1,i
            val=FP(i,j)-FP(j,i)
            if(abs(val).gt.tolerance) then
              write(*,*)'QMP FOCK MATRIX IS NOT SYMMETRIC FOR IJ=',I,J
              write(*,*)'>>>> FM(IJ)=',FP(i,j)
              write(*,*)'>>>> FM(JI)=',FP(j,i)
              LPSYM=.false.
            end if
         end do
      end do

      if(LASYM) then
         write(*,*)'FAE IS SYMMETRIC'
      end if
      if(LBSYM) then
         write(*,*)'FBE IS SYMMETRIC'
      end if
      if(LPSYM) then
         write(*,*)'FP  IS SYMMETRIC'
      end if


      return
      end
!=======================================================================
      subroutine UFM_sym_check2(nbf,F)

! Checks that the Fock matrices are symmetric
!=======================================================================
      implicit none

! Input Variables
      integer nbf
      double precision F(nbf,nbf)

! No Variables Returned

! Local variables
      logical LSYM
      logical LOUTPUT
      integer i,j
      double precision val,tolerance,maxdif
      parameter(tolerance=1.0d-10)

      maxdif=0.0d+00
      LSYM=.true.
      LOUTPUT=.false.

      do i=1,nbf
         do j=1,i
            val=F(i,j)-F(j,i)
            if(abs(val).gt.tolerance) then
              maxdif=val
              if(LOUTPUT) then
              write(*,*)' FOCK MATRIX IS NOT SYMMETRIC FOR IJ=',I,J
              write(*,*)'>>>> FM(IJ)=',F(i,j)
              write(*,*)'>>>> FM(JI)=',F(j,i)
              end if
              LSYM=.false.
            end if
         end do
      end do

      if(LSYM) then
         write(*,*)'FOCK MATRIX IS SYMMETRIC'
      else
         write(*,*)' FOCK MATRIX IS NOT SYMMETRIC, MAXDIFF =',maxdif
      end if


      return
      end
!======================================================================
      subroutine UHF_Fock_correction(nbf,Psi_H_Psi,S_total,
     x                               XF,XS1,XS2,F)
!======================================================================
      implicit none

! Input Variables
      integer nbf
      double precision Psi_H_Psi
      double precision S_total
      double precision XS1(nbf,nbf)
      double precision XS2(nbf,nbf)
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
      COEF2=Psi_H_Psi / (S_Total * S_total)

! Fock Matrix Correction
      do i=1,nbf
         do j=1,nbf

            F(i,j) = COEF1*XF(i,j) - COEF2*( XS1(i,j) + XS2(i,j) )

         end do
      end do


      return
      end
!======================================================================
      subroutine IC_construct_FAE(LNEOHF,nelec,NAE,NBE,
     x                            nebf,npbf,nebf2,ngee,
     x                            ng1,ng2,ng3,ng4,
     x                            SZG2ICR,SZG3ICR,SZG4ICR,
     x                            NG2CHK,NG3CHK,NG4CHK,
     x                            DAE,DBE,DETOT,DP,
     x                            GAM_ecore,GAM_ep,GAM_ee,
     x                            GM2ICR,GM2sICR,GM3ICR,GM4ICR,
     x                            FAE,SAE1,SAE2, 
     x                            E_AE_ecore,
     x                            E_AE_GAMep,
     x                            E_AE_GAMee,
     x                            E_AE_OMG1, 
     x                            E_AE_OMG2,
     x                            E_AE_OMG3,
     x                            E_AE_OMG4,
     x                            S_AE_OMG1,
     x                            S_AE_OMG2)

! In-core Construct Alpha electronic Fock Matrix
! nelec: total number of electrons
! NAE: number alpha electrons
! NBE: number beta electrons
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
      double precision GAM_ecore(nebf2)
      double precision GAM_ep(ng1)
      double precision GAM_ee(ngee)
!     double precision GM2ICR(ng2)
!     double precision GM2sICR(ng2)
!     double precision GM3ICR(ng3)
!     double precision GM4ICR(ng4)

      double precision GM2ICR(SZG2ICR)
      double precision GM2sICR(SZG2ICR)
      double precision GM3ICR(SZG3ICR)
      double precision GM4ICR(SZG4ICR)

! Variables Returned
      double precision FAE(nebf,nebf)
      double precision SAE1(nebf,nebf)
      double precision SAE2(nebf,nebf)
      double precision E_AE_ecore
      double precision E_AE_GAMep
      double precision E_AE_GAMee
      double precision E_AE_OMG1
      double precision E_AE_OMG2
      double precision E_AE_OMG3
      double precision E_AE_OMG4
      double precision S_AE_OMG1
      double precision S_AE_OMG2

! Local Variables
      double precision zero
      parameter(zero=0.0d+00)
!     double precision SAE1(nebf,nebf)
!     double precision SAE2(nebf,nebf)


! Initialize
      E_AE_ecore = zero
      E_AE_GAMep = zero
      E_AE_GAMee = zero

      E_AE_OMG1  = zero
      E_AE_OMG2  = zero
      E_AE_OMG3  = zero
      E_AE_OMG4  = zero

      S_AE_OMG1  = zero
      S_AE_OMG2  = zero

      FAE=zero
      SAE1=zero
      SAE2=zero

!     SAE1=zero
!     SAE2=zero

! NEO-UHF Contributions
      call FAE_ecore(nebf,nebf2,GAM_ecore,DAE,FAE,E_AE_ecore)

      call FAE_GAMep(nebf,npbf,ng1,DAE,DP,GAM_ep,FAE,E_AE_GAMep)

      ! More than 1 total electrons (nelec.ge.2)
!     write(*,*)'IN: IC_construct_FAE  '
!     write(*,*)'Before call FAE_GAMee '
!     write(*,*)
!     write(*,*)'DAE: ',DAE
!     write(*,*)'DETOT: ',DETOT
!     write(*,*)'E_AE_GAMee= ',E_AE_GAMee
!     write(*,*)
      call FAE_GAMee(nebf,ngee,GAM_ee,DAE,DETOT,FAE,E_AE_GAMee)
!     write(*,*)'IN: IC_construct_FAE  '
!     write(*,*)'After call FAE_GAMee '
!     write(*,*)'E_AE_GAMee= ',E_AE_GAMee
!     write(*,*)

      if(LNEOHF) RETURN
!     write(*,*)'IN: IC_construct_FAE '
!     write(*,*)'AFTER LNEOHF '

! NEO-XCUHF Contributions

      call FAE_OMG1(nebf,npbf,ng1,DAE,DP,FAE,SAE1,E_AE_OMG1,S_AE_OMG1)

      if(nelec.gt.1) then

         call FAE_OMG2(NG2CHK,NAE,NBE,nebf,npbf,ng2,
     x                DAE,DBE,DETOT,DP,GM2ICR,GM2SICR,
     x                FAE,SAE2,E_AE_OMG2,S_AE_OMG2)

      end if

      if(nelec.gt.2) then

         call FAE_OMG3(NG3CHK,nebf,npbf,ng3,NAE,NBE,
     x                DAE,DBE,DP,GM3ICR,FAE,E_AE_OMG3)

      end if

      if(nelec.gt.3) then

         call FAE_OMG4(NG4CHK,nebf,npbf,ng4,NAE,NBE,
     x                DAE,DBE,DP,GM4ICR,FAE,E_AE_OMG4)

      end if


      return
      end
!======================================================================
      subroutine IC_construct_FP(LNEOHF,nelec,NAE,NBE,
     x                           nebf,npbf,npbf2,ngee,
     x                           ng1,ng2,ng3,ng4,
     x                           SZG2ICR,SZG3ICR,SZG4ICR,
     x                           NG2CHK,NG3CHK,NG4CHK,
     x                           DAE,DBE,DETOT,DP,
     x                           GAM_pcore,GAM_ep,
     x                           GM2ICR,GM2sICR,GM3ICR,GM4ICR,
     x                           FP,SP1,SP2, 
     x                           E_pcore,
     x                           E_P_GAMep,
     x                           E_P_OMG1, 
     x                           E_P_OMG2,
     x                           E_P_OMG3,
     x                           E_P_OMG4,
     x                           S_P_OMG1,
     x                           S_P_OMG2)

! In-core Construct QM-Particle Fock Matrix
! nelec: total number of electrons
! NAE: number alpha electrons
! NBE: number beta electrons
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
      double precision DETOT(nebf,nebf)
      double precision GAM_pcore(npbf2)
      double precision GAM_ep(ng1)
      double precision GAM_ee(ngee)
!     double precision GM2ICR(ng2)
!     double precision GM2sICR(ng2)
!     double precision GM3ICR(ng3)
!     double precision GM4ICR(ng4)

      double precision GM2ICR(SZG2ICR)
      double precision GM2sICR(SZG2ICR)
      double precision GM3ICR(SZG3ICR)
      double precision GM4ICR(SZG4ICR)

! Variables Returned
      double precision FP(npbf,npbf)
      double precision SP1(npbf,npbf)
      double precision SP2(npbf,npbf)
      double precision E_pcore
      double precision E_P_GAMep
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
      SP1=zero    
      SP2=zero    

      E_pcore   = zero 
      E_P_GAMep = zero 
      E_P_OMG1  = zero 
      E_P_OMG2  = zero 
      E_P_OMG3  = zero 
      E_P_OMG4  = zero 
      S_P_OMG1  = zero
      S_P_OMG2  = zero
 

! NEO-UHF Contributions

      call FP_pcore(npbf,npbf2,GAM_pcore,DP,FP,E_pcore)

      call FP_GAMep(nebf,npbf,ng1,DETOT,DP,GAM_ep,FP,E_P_GAMep)

      if(LNEOHF) RETURN

! NEO-XCUHF Contributions

      call UFP_OMG1(nebf,npbf,ng1,DETOT,DP,FP,SP1,
     x              E_P_OMG1,S_P_OMg1)

      if(nelec.gt.1) then

         call UFP_OMG2(NG2CHK,nebf,npbf,ng2,NAE,NBE,
     x                 DAE,DBE,DETOT,DP,GM2ICR,GM2sICR,FP,SP2,
     x                 E_P_OMG2,S_P_OMG2)

      end if

      if(nelec.gt.2) then

         call UFP_OMG3(NG3CHK,nebf,npbf,ng3,NAE,NBE,
     x                 DAE,DBE,DP,GM3ICR,FP,E_P_OMG3)

      end if

      if(nelec.gt.3) then

         call UFP_OMG4(NG4CHK,nebf,npbf,ng4,NAE,NBE,
     x                 DAE,DBE,DP,GM4ICR,FP,E_P_OMG4)

      end if


      return
      end
