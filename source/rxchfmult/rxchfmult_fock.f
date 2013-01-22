!======================================================================
      subroutine rxchfmult_xchf_fock(LNEOHF,LGAM4,LG4DSCF,LG4IC,
     x                      LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                      NG4CHK,NG3CHK,NG2CHK,
     x                      SZG4IC,SZG3IC1,SZG2ICR,
     x                      npebf,nebf,nebf2,npbf,npbf2,nelec,
     x                      ngee,ng1,ng2,ng3,ng4,DE,DP,
     x                      GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                      GM4ICR,GM3IC1,GM2ICR,GM2SICR,
     x                      ng2prm,ngtg1,ng3prm,
     x                      nat,pmass,cat,zan,bcoef1,gamma1,
     x                      KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                      focke,fockp,
     x                      E_total,E_nuc,E_ecore,E_pcore,E_ep,
     x                      E_ee,E_gam1,E_gam2,E_gam3,E_gam4,
     x                      S_total,S_gam1,S_gam2)

C Essentially copied from ../xchf_independent/xchf1_fock.f
!======================================================================
      implicit none

! Input Variables
      logical LCMF     ! Check My Fock: turns on on-the-fly Fock check
      logical LNEOHF   ! Control for doing NEO-HF only
      logical LG4DSCF  ! Control for direct SCF for GAM4 term
      logical LG3DSCF  ! Control for direct SCF for GAM3 term
      logical LG2IC1   ! Control for SCF for in-core GAM2 term
      logical LG3IC1   ! Control for SCF for in-core GAM3 term
      logical LG4IC    ! Control for SCF for in-core GAM4 term
      logical LG2DSCF  ! Control for direct SCF for GAM2 term
      integer NG4CHK   ! Number of chunks for dividing GAM4 evaluation
      integer NG3CHK   ! Number of chunks for dividing GAM3 evaluation
      integer NG2CHK   ! Number of chunks for dividing GAM2 evaluation
      integer SZG2ICR
      integer SZG3IC1
      integer SZG4IC 
      integer nebf
      integer nebf2
      integer npbf
      integer npbf2
      integer nelec
!     integer ngtgmn

      integer ngee
      integer ng1
      integer ng2
      integer ng3
      integer ng4

      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GAM_ecore(nebf2)
      double precision GAM_pcore(npbf2)
      double precision GAM_ep(ng1)
      double precision GAM_ee(ngee)
      double precision GM2ICR(SZG2ICR)
      double precision GM2SICR(SZG2ICR)
      double precision GM3IC1(SZG3IC1)
      double precision GM4ICR(SZG4IC)
!     double precision efct(nebf)
!     double precision pfct(npbf)
!     double precision GAM_1s(ng1)
!     double precision GAM_1(ng1)
!     double precision GAM_2s(ng2)
!     double precision GAM_2(ng2)
!     double precision GAM_3PK(ng3)
!     double precision GAM_4PK(ng4)
!     double precision xxse(nebf,nebf)
!     double precision xxsp(npbf,npbf)
!-----DIRECT-SCF-RELATED-----------------------------------------------(
      integer ngtg1
      integer npebf
      integer ng2prm,ng3prm,nat
!-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
      integer KPESTR(nebf)  ! Map contracted index to primitive start
      integer KPEEND(nebf)  ! Map contracted index to primitive end
!-------Basis Set Info-------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)
!-----DIRECT-SCF-RELATED-----------------------------------------------)

! Variables Returned
      double precision focke(nebf,nebf)
      double precision fockp(npbf,npbf)
      double precision Xfocke(nebf,nebf)
      double precision Xfockp(npbf,npbf)
      double precision E_total
      double precision E_ecore
      double precision E_pcore
      double precision E_ep
      double precision E_ee
      double precision E_gam1
      double precision E_gam2
      double precision E_gam3
      double precision E_gam4
      double precision E_nuc
      double precision S_total
      double precision S_gam1
      double precision S_gam2

! Local variables
      integer npbfLT
      integer nebfLT
      double precision zero,one
      parameter(zero=0.0d+00,one=1.0d+00)
      double precision H_expect 
      double precision E_check 
      double precision E_NEOHF 
      double precision psiHpsi 
      double precision se1(nebf,nebf)
      double precision se2(nebf,nebf)
      double precision sp1(npbf,npbf)
      double precision sp2(npbf,npbf)
      double precision CHKfockp(npbf,npbf)
      double precision CHKS(npbf,npbf)

      logical LGAM4

!----------------(
!     LGAM4=.FALSE.
!     LGAM4=.TRUE.
!----------------)


! Initialize values for energy components
      psiHpsi=zero

      E_total=zero
      E_ecore=zero
      E_pcore=zero
      E_ep=zero
      E_ee=zero
      E_gam1=zero
      E_gam2=zero
      E_gam3=zero
      E_gam4=zero

      S_total=zero
      S_gam1=zero
      S_gam2=zero

! Zero-out the Fock matrices
!     call zero_out(nebf,focke)
!     call zero_out(npbf,fockp)
      focke=zero
      fockp=zero

! Zero-out the XFock matrices
!     call zero_out(nebf,Xfocke)
!     call zero_out(npbf,Xfockp)
      Xfocke=zero
      Xfockp=zero

! Zero-out the Fock-like overlap matrices
!     call zero_out(nebf,se1)
!     call zero_out(npbf,sp1)
      se1=zero
      sp1=zero

!     call zero_out(nebf,se2)
!     call zero_out(npbf,sp2)
      se2=zero
      sp2=zero


! Calculate energy components and build Fock matrices

!  Conventional NEO-HF components:
      call E_from_GAM_ecore(nebf,nebf2,
     * GAM_ecore,DE,focke,E_ecore)

!     nebfLT=(nebf*nebf+nebf)/2
!     write(*,*)'XCHF SIDE... Elec HCORE Matrix:'
!     call print_my_fock(nebf,nebfLT,focke)
!     write(*,*)'GAM_ecore integrals'
!     write(*,*) GAM_ecore

      call E_from_GAM_pcore(npbf,npbf2,
     * GAM_pcore,DP,fockp,E_pcore)

!     npbfLT=(npbf*npbf+npbf)/2
!     write(*,*)'XCHF SIDE... Nuc HCORE Matrix:'
!     call print_my_fock(npbf,npbfLT,fockp)

      call E_from_GAM_ep(nebf,npbf,ng1,
     * GAM_ep,DE,DP,focke,fockp,E_ep)

!     write(*,*)'XCHF SIDE... Nuc Fock Matrix:'
!     call print_my_fock(npbf,npbfLT,fockp)
!     write(*,*)'GAM_ep integrals'
!     write(*,*) GAM_ep

      if(nelec.gt.1) then
         call E_from_GAM_ee(nebf,ngee,
     *    GAM_ee,DE,focke,E_ee)
      else
         E_ee=0.0d+00
      end if
!         nebfLT=(nebf*nebf+nebf)/2
!         write(*,*)'XCHF SIDE... Elec Fock Matrix:'
!         call print_my_fock(nebf,nebfLT,focke)
!         call abrt

! Retain the standard NEOHF values for energy and Fock matrices
!  for use later
      Xfocke=focke
      Xfockp=fockp
      E_NEOHF=E_ecore+E_pcore+E_ep+E_ee+E_nuc
!         write(*,*)'NEOHF Energy is:',E_NEOHF

! IF NEO-HF ONLY GET OUT OF HERE:
      if(LNEOHF) then 
         E_total=E_NEOHF
         S_total=1.0d+00
!-----------ON-THE-FLY-CODE-TESTING------------------------------------(
         if(LCMF) then
            call Fock_sym_check(nebf,npbf,focke,fockp)
            call NonGem_OVLAP_check(nebf,npbf,nelec,DE,DP)
         end if
!-----------ON-THE-FLY-CODE-TESTING------------------------------------)
         RETURN
      end if
!  XC components:
!     if(ngtgmn.gt.0) then

         call S_from_GAM_1s(nebf,npbf,ng1,
     *    DE,DP,se1,sp1,S_gam1)

         call E_from_GAM_1(nebf,npbf,ng1,
     *    DE,DP,focke,fockp,E_gam1)

         if(nelec.gt.1) then

            if(LG2DSCF) then
               call GAM2_DSCF(NG2CHK,nebf,npebf,npbf,
     x                        ng2,ng2prm,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        DE,DP,focke,fockp,se2,sp2,E_gam2,S_gam2)
            else if(LG2IC1) then
!     subroutine EG2IC1(Nchunks,nebf,npbf,ng2,
!    x                  DE,DP,GM2ICR,GM2SICR,
!    x                  FE,SE,FP,SP,E_gam2,S_gam2)
               call EG2IC1(NG2CHK,nebf,npbf,ng2,
     x                     DE,DP,GM2ICR,GM2SICR,
     x                     focke,se2,fockp,sp2,E_gam2,S_gam2)
            else
               call S_from_GAM_2s(nebf,npbf,ng2,
     *         DE,DP,se2,sp2,S_gam2)

               call E_from_GAM_2(nebf,npbf,ng2,
     *         DE,DP,focke,fockp,E_gam2)
            end if

            if(nelec.gt.2) then

!              write(*,*)'about to call E_from_GAM_3'
!              call flshbf(6)
               if(LG3DSCF) then
                  call GAM3_DSCF(NG3CHK,nebf,npebf,npbf,
     x                           ng3,ng3prm,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                           DE,DP,focke,fockp,E_gam3)
               else if(LG3IC1) then
                  call EG3IC1(NG3CHK,nebf,npbf,ng3,
     x                        DE,DP,GM3IC1,focke,fockp,E_gam3)
               else
                  call E_from_GAM_3(nebf,npbf,ng3,
     x             DE,DP,focke,fockp,E_gam3)
               end if
!              write(*,*)'about to call E_from_GAM_4'
!              call flshbf(6)
               if(LGAM4) then
                  if(LG4DSCF) then
                     call GAM4_DSCF(NG4CHK,nebf,npbf,ngee,ng2,ng4,
     x                              DE,DP,focke,fockp,E_gam4)
                  else if(LG4IC) then
                     call EG4ICR(NG4CHK,nebf,npbf,ng4,
     x                           DE,DP,GM4ICR,focke,fockp,E_gam4)
                  else
                     call E_from_GAM_4(nebf,npbf,ng4,
     x                DE,DP,focke,fockp,E_gam4)
                  end if
               end if

            end if

         end if

!     end if

!  Calculate total overlap
! CWS-TEST the non-geminal overlap contribution is included in S_gam1
      S_total=one+S_gam1+S_gam2
!     S_total=S_gam1+S_gam2

!  Calculate total energy
      psiHpsi=(E_ecore+E_pcore+E_ep+E_ee+
     x        E_gam1+E_gam2+E_gam3+E_gam4)/S_total

      E_total=psiHpsi+E_nuc

!     H_expect=E_ecore+E_pcore+E_ep+E_ee+
!    x         E_gam1+E_gam2+E_gam3+E_gam4
      H_expect=psiHpsi/S_total

!  Make XC-related corrections to the Fock matrices
      call XCfock_correction(nebf,npbf,
     *            E_total,S_total,H_expect,
     *            se1,sp1,se2,sp2,focke,fockp,CHKfockp,CHKS)

!-----------ON-THE-FLY-CODE-TESTING------------------------------------(
      if(LCMF) then
         npbfLT=npbf*(npbf+1)/2
         nebfLT=nebf*(nebf+1)/2
!!       E_check=psiHpsi-E_ecore-E_ee
!!       E_check=E_ecore+E_pcore+E_ep+E_ee+E_gam1+E_gam2+E_gam3+E_gam4
         E_check=E_pcore+E_ep+E_gam1+E_gam2+E_gam3+E_gam4
         call CHK_my_fockp(npbf,npbfLT,E_check,S_total,DP,
     x                     CHKfockp,CHKS)

         call Focking_Around(nebf,nebfLT,npbf,npbfLT,
     x                       DE,DP,focke,fockp,se1,se2,sp1,sp2,
     x                       psiHpsi,S_total,E_pcore,E_ecore,
     x                       E_ep,E_ee,E_gam1,E_gam2,E_gam3,E_gam4)

         call Fock_sym_check(nebf,npbf,focke,fockp)
         call NonGem_OVLAP_check(nebf,npbf,nelec,DE,DP)
      end if
!-----------ON-THE-FLY-CODE-TESTING------------------------------------)

      return
      end
