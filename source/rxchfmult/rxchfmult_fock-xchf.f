C======================================================================
      subroutine RXCHFmult_fock_xchf(LGAM4,LG4DSCF,LG4IC,
     x                      LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                      NG4CHK,NG3CHK,NG2CHK,
     x                      SZG4IC,SZG3IC1,SZG2ICR,
     x                      npebf,nebf,nebf2,npbf,npbf2,nelec,
     x                      ngee,ng1,ng2,ng3,ng4,DE,DP,
     x                      GM4ICR,GM3IC1,GM2ICR,GM2SICR,
     x                      ng2prm,ngtg1,ng3prm,
     x                      nat,pmass,cat,zan,bcoef1,gamma1,
     x                      KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                      focke,fockp,SEtot,SPtot,
     x                      E_total,E_gam1,E_gam2,
     x                      E_gam3,E_gam4,
     x                      S_total,S_gam1,S_gam2)

C Adapted from ../xchf_independent/xchf1_fock.f to have G-only terms
C======================================================================
      implicit none

C Input variables
      logical           LCMF             ! Check My Fock: turns on on-the-fly Fock check
      logical           LG4DSCF          ! Control for direct SCF for GAM4 term
      logical           LG3DSCF          ! Control for direct SCF for GAM3 term
      logical           LG2IC1           ! Control for SCF for in-core GAM2 term
      logical           LG3IC1           ! Control for SCF for in-core GAM3 term
      logical           LG4IC            ! Control for SCF for in-core GAM4 term
      logical           LG2DSCF          ! Control for direct SCF for GAM2 term
      integer           NG4CHK           ! Number of chunks for dividing GAM4 evaluation
      integer           NG3CHK           ! Number of chunks for dividing GAM3 evaluation
      integer           NG2CHK           ! Number of chunks for dividing GAM2 evaluation
      integer           SZG2ICR
      integer           SZG3IC1
      integer           SZG4IC 
      integer           nebf
      integer           nebf2
      integer           npbf
      integer           npbf2
      integer           nelec
      integer           ngee
      integer           ng1
      integer           ng2
      integer           ng3
      integer           ng4
      double precision  DE(nebf,nebf)
      double precision  DP(npbf,npbf)
      double precision  GAM_ecore(nebf2)
      double precision  GAM_pcore(npbf2)
      double precision  GAM_ep(ng1)
      double precision  GAM_ee(ngee)
      double precision  GM2ICR(SZG2ICR)
      double precision  GM2SICR(SZG2ICR)
      double precision  GM3IC1(SZG3IC1)
      double precision  GM4ICR(SZG4IC)
C Direct SCF
      integer           ngtg1
      integer           npebf
      integer           ng2prm,ng3prm,nat
      integer           ELCAM(npebf,3)   ! Angular mom for electrons
      integer           NUCAM(npbf,3)    ! Angular mom for quantum nuclei
      integer           AMPEB2C(npebf)   ! Map primitive index to contracted
      integer           KPESTR(nebf)     ! Map contracted index to primitive start
      integer           KPEEND(nebf)     ! Map contracted index to primitive end
      double precision  ELCEX(npebf)     ! Exponents: elec basis
      double precision  NUCEX(npbf)      ! Exponents: nuc basis
      double precision  ELCBFC(npebf,3)  ! Basis centers: elec basis
      double precision  NUCBFC(npbf,3)   ! basis centers: nuc basis
      double precision  AGEBFCC(npebf)   ! Map prim index to contract coef
      double precision  AGNBFCC(npbf)    ! Nuclear contract coef
      double precision  pmass            ! Mass of nonelectron quantum particle 
      double precision  zan(nat)         ! Classical nuclear charges
      double precision  cat(3,nat)       ! XYZ Coordinates of atoms
      double precision  bcoef1(ngtg1)
      double precision  gamma1(ngtg1)

C Output variables
      double precision  focke(nebf,nebf)
      double precision  fockp(npbf,npbf)
      double precision  E_total
      double precision  E_gam1
      double precision  E_gam2
      double precision  E_gam3
      double precision  E_gam4
      double precision  S_total
      double precision  S_gam1
      double precision  S_gam2
      double precision  SEtot(nebf,nebf)  ! Overlap contributions
      double precision  SPtot(npbf,npbf)  !

C Local variables
      logical           LGAM4
      integer           i,j
      integer           npbfLT
      integer           nebfLT
      double precision  H_expect 
      double precision  E_check 
      double precision  E_NEOHF 
      double precision  psiHpsi 
      double precision  se1(nebf,nebf)
      double precision  se2(nebf,nebf)
      double precision  sp1(npbf,npbf)
      double precision  sp2(npbf,npbf)
      double precision  CHKfockp(npbf,npbf)
      double precision  CHKS(npbf,npbf)
      double precision  zero,one
      parameter(zero=0.0d+00,one=1.0d+00)


C Variables for overlap contributions to Fock matrices

! Initialize
      psiHpsi=zero

      E_total=zero
      E_gam1=zero
      E_gam2=zero
      E_gam3=zero
      E_gam4=zero

      S_total=zero
      S_gam1=zero
      S_gam2=zero

      focke=zero
      fockp=zero

      se1=zero
      se2=zero
      sp1=zero
      sp2=zero


! Calculate energy components and build Fock matrices

      call S_from_GAM_1s(nebf,npbf,ng1,DE,DP,se1,sp1,S_gam1)

      call E_from_GAM_1(nebf,npbf,ng1,DE,DP,focke,fockp,E_gam1)

      if(nelec.gt.1) then

         if(LG2DSCF) then
            call GAM2_DSCF(NG2CHK,nebf,npebf,npbf,
     x                     ng2,ng2prm,nat,ngtg1,
     x                     pmass,cat,zan,bcoef1,gamma1,
     x                     KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                     ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                     DE,DP,focke,fockp,se2,sp2,E_gam2,S_gam2)
         else if(LG2IC1) then
            call EG2IC1(NG2CHK,nebf,npbf,ng2,
     x                  DE,DP,GM2ICR,GM2SICR,
     x                  focke,se2,fockp,sp2,E_gam2,S_gam2)
         else
            call S_from_GAM_2s(nebf,npbf,ng2,DE,DP,se2,sp2,S_gam2)

            call E_from_GAM_2(nebf,npbf,ng2,DE,DP,focke,fockp,E_gam2)
         end if

         if(nelec.gt.2) then

            if(LG3DSCF) then
               call GAM3_DSCF(NG3CHK,nebf,npebf,npbf,
     x                        ng3,ng3prm,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        DE,DP,focke,fockp,E_gam3)
            else if(LG3IC1) then
               call EG3IC1(NG3CHK,nebf,npbf,ng3,
     x                     DE,DP,GM3IC1,focke,fockp,E_gam3)
            else
               call E_from_GAM_3(nebf,npbf,ng3,
     x          DE,DP,focke,fockp,E_gam3)
            end if

            if(LGAM4) then
               if(LG4DSCF) then
                  call GAM4_DSCF(NG4CHK,nebf,npbf,ngee,ng2,ng4,
     x                           DE,DP,focke,fockp,E_gam4)
               else if(LG4IC) then
                  call EG4ICR(NG4CHK,nebf,npbf,ng4,
     x                        DE,DP,GM4ICR,focke,fockp,E_gam4)
               else
                  call E_from_GAM_4(nebf,npbf,ng4,
     x             DE,DP,focke,fockp,E_gam4)
               end if
            end if

         end if

      end if

!  Calculate total overlap
      S_total=S_gam1+S_gam2

!  Calculate total energy
      psiHpsi=(E_gam1+E_gam2+E_gam3+E_gam4)/S_total

      E_total=psiHpsi

      H_expect=psiHpsi/S_total

!  Make XC-related corrections to the Fock matrices
      call XCfock_correction(nebf,npbf,
     *            E_total,S_total,H_expect,
     *            se1,sp1,se2,sp2,focke,fockp,CHKfockp,CHKS)

      if(LCMF) then
         call Fock_sym_check(nebf,npbf,focke,fockp)

C Adapt testing below for G ansatz
C         npbfLT=npbf*(npbf+1)/2
C         nebfLT=nebf*(nebf+1)/2
C         E_check=E_pcore+E_ep+E_gam1+E_gam2+E_gam3+E_gam4
C         call CHK_my_fockp(npbf,npbfLT,E_check,S_total,DP,
C     x                     CHKfockp,CHKS)
C
C         call Focking_Around(nebf,nebfLT,npbf,npbfLT,
C     x                       DE,DP,focke,fockp,se1,se2,sp1,sp2,
C     x                       psiHpsi,S_total,E_pcore,E_ecore,
C     x                       E_ep,E_ee,E_gam1,E_gam2,E_gam3,E_gam4)
C
C         call NonGem_OVLAP_check(nebf,npbf,nelec,DE,DP)

        nebflt=nebf*(nebf+1)/2
        npbflt=npbf*(npbf+1)/2
        write(*,*) "FBE XCHF:"
        call prt_lower_triangle(nebf,nebflt,focke)
        write(*,*)
        write(*,*) "FP XCHF:"
        call prt_lower_triangle(npbf,npbflt,fockp)


      end if

C Store overlap contributions to Fock matrices separately for interaction routine
      SEtot=0.0d+00
      do i=1,nebf
        do j=1,nebf
          SEtot(j,i)=se1(j,i)+2.0d+00*se2(j,i)
        end do
      end do

      SPtot=0.0d+00
      do i=1,npbf
        do j=1,npbf
          SPtot(j,i)=sp1(j,i)+sp2(j,i)
        end do
      end do

      return
      end

