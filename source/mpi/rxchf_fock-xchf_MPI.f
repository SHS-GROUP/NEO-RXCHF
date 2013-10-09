C======================================================================
      subroutine RXCHF_fock_xchf_MPI(nproc,rank,
     x                               LCMF,nelec,
     x                               nebf,nebf2,npbf,npbf2,
     x                               ngee,ng1,ng2,ng3,ng4,
     x                               SZG2ICR,SZG3ICR,SZG4ICR,
     x                               NG2CHK,NG3CHK,NG4CHK,
     x                               DE,DP,
     x                               GM2ICR,GM2sICR,GM3ICR,GM4ICR,
     x                               focke,fockp,SEtot,SPtot,
     x                               E_total,E_gam1,E_gam2,
     x                               E_gam3,E_gam4,
     x                               S_total,S_gam1,S_gam2)

C XCHF energy and resulting contributions to Fock matrices
C Loop over integrals, constructing Fock matrices simultaneously
C======================================================================
      implicit none

C Input variables
      integer           nproc,rank
      logical           LCMF             ! Check My Fock: turns on on-the-fly Fock check
      integer           NG2CHK           ! Number of chunks for dividing GAM2 evaluation
      integer           NG3CHK           ! Number of chunks for dividing GAM3 evaluation
      integer           NG4CHK           ! Number of chunks for dividing GAM4 evaluation
      integer           SZG2ICR
      integer           SZG3ICR
      integer           SZG4ICR
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
      double precision  GM2sICR(SZG2ICR)
      double precision  GM3ICR(SZG3ICR)
      double precision  GM4ICR(SZG4ICR)

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


! Include GAM1 contributions

      call S_from_GAM_1s(nebf,npbf,ng1,DE,DP,se1,sp1,S_gam1)
      call E_from_GAM_1(nebf,npbf,ng1,DE,DP,focke,fockp,E_gam1)

! Include GAM2 contributions
      if(nelec.gt.1) then

        call XCHF_OMG2_MPI(nproc,rank,
     x                     NG2CHK,nebf,npbf,
     x                     ng2,SZG2ICR,
     x                     DE,DP,
     x                     GM2ICR,GM2sICR,
     x                     focke,se2,fockp,sp2,E_gam2,S_gam2)

! Include GAM3 contributions
        if(nelec.gt.2) then

          call XCHF_OMG3_MPI(nproc,rank,
     x                       NG3CHK,nebf,npbf,
     x                       ng3,SZG3ICR,
     x                       DE,DP,
     x                       GM3ICR,
     x                       focke,fockp,E_gam3)

! Include GAM4 contributions
          if(nelec.gt.3) then

            call XCHF_OMG4_MPI(nproc,rank,
     x                         NG4CHK,nebf,npbf,
     x                         ng4,SZG4ICR,
     x                         DE,DP,
     x                         GM4ICR,
     x                         focke,fockp,E_gam4)

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
        call RXCHFmult_xchfFock_testing(nebf,npbf,focke,fockp,DE,DP,
     x                                  E_total,S_total,
     x                                  E_gam2,E_gam3,E_gam4,S_gam2)
        call Fock_sym_check(nebf,npbf,focke,fockp)
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

