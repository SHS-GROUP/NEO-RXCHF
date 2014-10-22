!======================================================================
      subroutine RXCHFmult_scf(nelec,NAE,NBE,NPRA,NPRB,NUCST,
     x                         npebf,nebf,nebf2,nebflt,
     x                         npebfBE,nebfBE,nebfBE2,nebfBElt,elindBE,
     x                         npbf,npbf2,npbflt,
     x                         ngtg1,ngee,
     x                         NG2CHK,NG3CHK,NG4CHK,
     x                         read_CE,read_CP,
     x                         LG4DSCF,LG3DSCF,LG2DSCF,
     x                         LSOSCF,LOCBSE,LCMF,LALTBAS,LADDEXCH,
     x                         nat,pmass,cat,zan,
     x                         bcoef1,gamma1,
     x                         KPESTR,KPEEND,
     x                         AMPEB2C,AGEBFCC,
     x                         ELCEX,ELCAM,ELCBFC,
     x                         KPESTR_be,KPEEND_be,
     x                         AMPEB2C_be,AGEBFCC_be,
     x                         ELCEX_be,ELCAM_be,ELCBFC_be,
     x                         AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                         LG2IC1,dimXCHF2,dimINT2,dimINT2ex,
     x                         XCHF_GAM2,INT_GAM2,INT_GAM2ex,XCHF_GAM2s,
     x                         LG3IC1,dimXCHF3,dimINT3,dimINT3ex,
     x                         XCHF_GAM3,INT_GAM3,
     x                         INT_GAM3ex1,INT_GAM3ex2,
     x                         LG4IC,dimXCHF4,dimINT4,
     x                         XCHF_GAM4,INT_GAM4)

!
! PERFORM A NUCLEAR-ELECTRONIC RESTRICTED XC HARTREE-FOCK CALCULATION
! FOR AN NAE-REGULAR ELECTRON NBE-SPECIAL ELECTRON ONE-PROTON SYSTEM
!
!     **DEFINITIONS:
!
!     *FOR REGULAR ELECTRONS:
!     NAE    ::  NUMBER OF REGULAR ELECTRONS
!     DAE    ::  NEW REGULAR ELECTRON DENSITY MATRIX
!     DAE0   ::  OLD REGULAR ELECTRON DENSITY MATRIX
!     VECAE  ::  REGULAR ELECTRON MOS
!     VECAE0 ::  OLD REGULAR ELECTRON MOS
!     AEE    ::  REGULAR ELECTRON ORBITAL EIGENVALUES
!
!     * FOR SPECIAL ELECTRONS:
!     NBE    ::  NUMBER OF SPECIAL ELECTRONS
!     DBE    ::  NEW SPECIAL ELECTRON DENSITY MATRIX
!     DBE0   ::  OLD SPECIAL ELECTRON DENSITY MATRIX
!     VECBE  ::  SPECIAL ELECTRON MOS
!     VECBE0 ::  OLD SPECIAL ELECTRON MOS
!     BEE    ::  SPECIAL ELECTRON ORBITAL EIGENVALUES
!
!     *FOR PROTONS:
!     DP    ::  NEW PROTON DENSITY MATRIX
!     DP0   ::  OLD PROTON DENSITY MATRIX
!     VECP  ::  PROTON MOS
!     EP    ::  PROTON ORBITAL EIGENVALUES
!
!======================================================================
      implicit none
! Input Variables
      logical LNEOHF
      logical LOCBSE   ! Use OCBSE scheme as is (restricted variational freedom for reg/sp elecs)
      logical LOCBSE2  ! Use modified OCBSE scheme (complete variational freedom for reg elecs)
      logical LALTBAS  ! Flag to denote distinct special electron basis
      logical read_CE
      logical read_CP
      logical LGAM4
      logical LCMF
      logical LADDEXCH
      logical LG4DSCF
      logical LG3DSCF
      logical LG2IC1
      logical LG3IC1
      logical LG4IC
      logical LG2DSCF
      integer SZG2ICR
      integer SZG3IC1
      integer SZG4IC
      integer NG2CHK,NG3CHK,NG4CHK
      integer nelec
      integer NAE,NBE
      integer NPRA,NPRB
      integer NUCST
      integer nebf,nebfBE,npbf
      integer nebf2,nebflt
      integer nebfBE2,nebfBElt
      integer npbf2,npbflt
      integer ngee
!-----DIRECT-SCF-RELATED-----------------------------------------------(
      integer ngtg1
      integer npebf,npebfBE
      integer nat
!-------Basis Set Info-------(
      integer ELCAM(npebf,3)                ! Angular mom for electrons
      integer NUCAM(npbf,3)                 ! Angular mom for quantum nuclei
      double precision ELCEX(npebf)         ! Exponents: elec basis
      double precision NUCEX(npbf)          ! Exponents: nuc basis
      double precision ELCBFC(npebf,3)      ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)       ! basis centers: nuc basis
      integer AMPEB2C(npebf)                ! Map primitive index to contracted
      double precision AGEBFCC(npebf)       ! Map prim index to contract coef
      double precision AGNBFCC(npbf)        ! Nuclear contract coef
      integer KPESTR(nebf)                  ! Map contracted index to primitive start
      integer KPEEND(nebf)                  ! Map contracted index to primitive end
! Special electron basis
      integer elindBE(nebfBE)               ! Contracted indices of NBE basis set
      integer ELCAM_be(npebfBE,3)           ! 
      double precision ELCEX_be(npebfBE)    ! 
      double precision ELCBFC_be(npebfBE,3) ! 
      integer AMPEB2C_be(npebfBE)           ! Analogs for special electron basis
      double precision AGEBFCC_be(npebfBE)  ! 
      integer KPESTR_be(nebfBE)             ! 
      integer KPEEND_be(nebfBE)             ! 
!-------Basis Set Info-------)
      double precision pmass                ! Mass of nonelectron quantum particle 
      double precision zan(nat)             ! Classical nuclear charges
      double precision cat(3,nat)           ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)
!-----DIRECT-SCF-RELATED-----------------------------------------------)

      integer          dimXCHF2,dimXCHF3,dimXCHF4
      integer          dimINT2,dimINT3,dimINT4
      integer          dimINT2ex,dimINT3ex
      double precision XCHF_GAM2(dimXCHF2)         ! XCHF GAM2 integrals
      double precision INT_GAM2(dimINT2)           ! Interaction GAM2 integrals
      double precision INT_GAM2ex(dimINT2)         ! Exchange GAM2 integrals
      double precision XCHF_GAM2s(dimXCHF2)        ! XCHF GAM2s integrals
      double precision XCHF_GAM3(dimXCHF3)         ! XCHF GAM3 integrals
      double precision INT_GAM3(dimINT3)           ! Interaction GAM3 integrals
      double precision INT_GAM3ex1(dimINT3)        ! Exchange GAM3 integrals
      double precision INT_GAM3ex2(dimINT3)        ! Exchange GAM3 integrals
      double precision XCHF_GAM4(dimXCHF4)         ! XCHF GAM4 integrals
      double precision INT_GAM4(dimINT4)           ! Interaction GAM4 integrals

! Local variables
      double precision zero,one
      PARAMETER (ZERO=0.0D+00, ONE=1.0D+00) 

      double precision xxse(nebf,nebf)  ! Elec overlap matrix
      double precision xxseBE(nebfBE,nebfBE)  ! Elec overlap matrix
      double precision xxsp(npbf,npbf)  ! Nuc overlap matrix
      double precision GAM_pcore(npbf2)
      double precision GAM_ecore(nebf2)
      double precision GAM_ecoreBE(nebfBE2)
!      double precision GAM_ep(ng1)
      double precision GAM_ee(ngee)

!     integer noccE ! Number of occupied elec orbs
!     integer noccP ! Number of occupied nuc orbs
      integer maxit,maxmicroit

      integer i,ielec
      integer j
      integer k
      integer l

      double precision TOLP
      double precision TOLE
      double precision diffE
      double precision diffAE
      double precision diffBE
      double precision diffP

      double precision E_total

      double precision E_HF
      double precision E_HF_ecore
      double precision E_HF_ee

      double precision E_XCHF
      double precision E_XCHF_gam1
      double precision E_XCHF_gam2
      double precision E_XCHF_gam3
      double precision E_XCHF_gam4

      double precision E_int
      double precision E_int_OMG2
      double precision E_int_OMG3
      double precision E_int_OMG4

      double precision E_nuc

      double precision S_total
      double precision S_gam1
      double precision S_gam2

      double precision DAE(NEBF,NEBF)
      double precision DAE0(NEBF,NEBF)
      double precision DBE(nebfBE,nebfBE)
      double precision DBE0(nebfBE,nebfBE)
      double precision VECAE(NEBF,NEBF)
      double precision VECAE0(NEBF,NEBF)
      double precision VECBE(nebfBE,nebfBE)
      double precision VECBE0(nebfBE,nebfBE)
      double precision AEE(NEBF)
      double precision BEE(nebfBE)
      double precision FAE(nebf,nebf)
      double precision XFAE(nebf,nebf)
      double precision FBE(nebfBE,nebfBE)
      double precision XFBE(nebfBE,nebfBE)

      double precision DP(NPBF,NPBF)
      double precision DP0(NPBF,NPBF)
      double precision VECP(NPBF,NPBF)
      double precision EP(NPBF)
      double precision FP(npbf,npbf)
      double precision XFP(npbf,npbf)

      double precision FAEint(nebf,nebf)
      double precision FBEint(nebfBE,nebfBE)
      double precision FPint(npbf,npbf)

      double precision SBE_XCHF(nebfBE,nebfBE)
      double precision SP_XCHF(npbf,npbf)

      double precision E_total_old
      double precision Delta_E_tot

      logical LDIFFE

!--------SOSCF-RELATED-VARIABLES------------(
      logical LSOSCF,LSOSCFA,LSOSCFB
      logical EIGAVL
      integer NA
      integer NB
      integer ITER
      integer ITSOA ! SOSCF iteration counter
      integer ITSOB ! SOSCF iteration counter
      integer L0,L1
      integer L0b,L1b
      integer L0w,L1w
      integer NFT15
      integer NFT16
      double precision FLT(NEBFLT) !FLT: Lower triangle FAE
      double precision FLTB(NEBFBELT) !FLTB: Lower triangle FBE
      double precision HSTARTA(NPRA)
      double precision GRADA(NPRA)
      double precision PGRADA(NPRA)
      double precision DISPLIA(NPRA)
      double precision DGRADA(NPRA)  ! WRK1
      double precision DISPLA(NPRA)  ! WRK2
      double precision UPDTA(NPRA)   ! WRK3
      double precision DISPLNA(NPRA) ! WRK1+NPR
      double precision DGRADIA(NPRA) ! WRK2+NPR
      double precision UPDTIA(NPRA)  ! WRK3+NPR
      double precision XA(NPRA)
      double precision ORBGRDA
      double precision GA(nebf,nebf) !G(L0,L0)
      double precision WRK(nebf) !WRK(L0)
      double precision, allocatable :: HSTARTB(:)
      double precision, allocatable :: GRADB(:)
      double precision, allocatable :: PGRADB(:)
      double precision, allocatable :: DISPLIB(:)
      double precision, allocatable :: DGRADB(:)  ! WRK1
      double precision, allocatable :: DISPLB(:)  ! WRK2
      double precision, allocatable :: UPDTB(:)   ! WRK3
      double precision, allocatable :: DISPLNB(:) ! WRK1+NPR
      double precision, allocatable :: DGRADIB(:) ! WRK2+NPR
      double precision, allocatable :: UPDTIB(:)  ! WRK3+NPR
      double precision, allocatable :: XB(:)
      double precision GB(nebfBE,nebfBE) !G(L0b,L0b)
      double precision WRKB(nebfBE) !WRK(L0b)
      double precision ORBGRDB
      double precision SMALL
      double precision SOGTOL ! ORBGRAD TOL to activate soscf
!cc   double precision CCC(nebf,nebf) !WRK(L0)
!cc   NPR=(L0-NA)*NA ! Line 2134 RHFCL ?NA is NUM ALPHA E?
!--------SOSCF-RELATED-VARIABLES------------)
C ARS( OCBSE variables
      integer nocca,noccb
      integer nwbf                                    ! # bf in W basis (nebf-nocca)
      integer nwbflt                                  ! one-dimensional size
      double precision, allocatable :: WB(:,:)        ! transformation matrix
      double precision, allocatable :: wFBEw(:,:)     ! Fock matrix in W basis
      double precision, allocatable :: wvecBEw(:,:)   ! eigenvectors in W basis from curr it
      double precision, allocatable :: wBEenw(:)      ! eigenvalues
      double precision, allocatable :: wFLTw(:)       ! flattened wFBEw
      double precision, allocatable :: wGBw(:,:)      ! exponential transformation
      double precision, allocatable :: wWRKw(:)       ! work array (L0w)
C )
C ARS( alt basis variables
      integer ialt,jalt,kalt,lalt
      integer dimint,dimint0
      double precision Ctemp(nebfBE,nebf)
      double precision,allocatable :: Cint(:,:)
      double precision ovlapalt
      double precision, parameter :: tolalt=1.0d-12
C )
C ARS( testing variables
      logical LNOINT
      double precision FBEmo(nebf,nebf)
      double precision work1(nebf,nebf)
      double precision work2(nebf,nebf)
      integer ierr
C )

!--------OUTPUT-FORMATTING---------------------------------------------(
 9000 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS ')

 9050 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS         ',
     * 'ORBGRAD_A ')

 9051 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS         ',
     * 'ORBGRAD_B ')

 9052 FORMAT(/' ITER      TOTAL ENERGY        E CHANGE       ',
     * 'ALPHA DENS       BETA DENS        QMP DENS         ',
     * 'ORBGRAD_A        ORBGRAD_B ')

 9100 FORMAT(1X,I3,F20.10,F17.10,3F17.10)

 9150 FORMAT(1X,I3,F20.10,F17.10,4F17.10)

 9151 FORMAT(1X,I3,F20.10,F17.10,5F17.10)

 9200 FORMAT(/1X,'FINAL NEORXCHF ENERGY IS',F20.10,' AFTER',I4,
     *           ' ITERATIONS')

 9300 FORMAT(/6X,'-----------------------------------------------',/,
     x        6X,'         NEORXCHF ENERGETIC COMPONENTS         ',/,
     x        6X,'-----------------------------------------------',/,
     x       12X,'    E_NUC=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,'  E_ECORE=',1X,F20.10/
     x       12X,'     E_EE=',1X,F20.10/
     x       12X,'     E_HF=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,'   E_GAM1=',1X,F20.10/
     x       12X,'   E_GAM2=',1X,F20.10/
     x       12X,'   E_GAM3=',1X,F20.10/
     x       12X,'   E_GAM4=',1X,F20.10/
     x       12X,'   E_XCHF=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,'   E_OMG2=',1X,F20.10/
     x       12X,'   E_OMG3=',1X,F20.10/
     x       12X,'   E_OMG4=',1X,F20.10/
     x       12X,'    E_int=',1X,F20.10/
     x       12X,'---------------------------------',/,
     x       12X,'  S_TOTAL=',1X,F20.10/
     x       12X,'  E_TOTAL=',1X,F20.10/
     x        6X,'-----------------------------------------------',/)

 9400 FORMAT(/1X,'          INITIAL GUESS ENERGETICS:            ')

 9500 FORMAT(/6X,' ** BEGIN RXCHF SELF-CONSISTENT-FIELD CALCULATION **')

 9610 FORMAT(/1X,' REGULAR ELECTRONIC ORBITALS AND EIGENVALUES:       ')

 9620 FORMAT(/1X,' SPECIAL ELECTRONIC ORBITALS AND EIGENVALUES:       ')

 9700 FORMAT(/1X,'      QM PARTICLE ORBITALS AND EIGENVALUES:         ')

 9800 FORMAT(10X,15(1H-),'START SECOND ORDER SCF',15(1H-))

 2001 FORMAT(/1X,'STARTING MICROITERATIONS FOR ITERATION',1X,I3)

 2000 FORMAT(1X,'CONVERGED ITERATION',1X,I3,1X,'IN',
     x       1X,I3,1X,'MICROITERATIONS',/)
                                           
!--------OUTPUT-FORMATTING---------------------------------------------)

      LOCBSE2=LOCBSE
C      LOCBSE2=.true. 

      if(LOCBSE2) then

       LOCBSE=.false.
       write(*,*) "Using LOCBSE2"

       if (nae.gt.1) then
        nocca=nae/2
       else
        nocca=nae
       end if
       if (nbe.gt.1) then
        noccb=nbe/2
       else
        noccb=nbe
       end if

C       noccvirta=nebf-noccb
C       noccvirtb=nebf-nocca

       if(.not.LALTBAS) then ! distinct sp elec bas => alloc at each it

        nwbf=nebf-nocca
        nwbflt=nwbf*(nwbf+1)/2
        L0w=nwbf
        L1w=nwbf
        if(allocated(WB)) deallocate(WB)
        allocate(WB(nebf,nwbf))
        if(allocated(wFBEw)) deallocate(wFBEw)
        allocate(wFBEw(nwbf,nwbf))
        if(allocated(wvecBEw)) deallocate(wvecBEw)
        allocate(wvecBEw(nwbf,nwbf))
        if(allocated(wBEenw)) deallocate(wBEenw)
        allocate(wBEenw(nwbf))
        if(allocated(wFLTw)) deallocate(wFLTw)
        allocate(wFLTw(nwbflt))
        if(allocated(wGBw)) deallocate(wGBw)
        allocate(wGBw(nwbf,nwbf))
        if(allocated(wWRKw)) deallocate(wWRKw)
        allocate(wWRKw(nwbf))
       else
        dimint0=0    ! Initialize so != dimint for control statement later
       end if

      end if

      if(LOCBSE) write(*,*) "Using LOCBSE"
      LGAM4=.true. ! Always calculate five-particle integrals
C ARS( no interaction
      LNOINT=.false.
C      LNOINT=.true.
      if (LNOINT) then
       write(*,*)
       write(*,*) "******************"
       write(*,*) "  NO INTERACTION  "
       write(*,*) "******************"
       write(*,*)
      end if
C )

!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------(
!      call class_nuc_rep(nat,zan,cat,E_nuc)
      open(800,file='ENUCRP.dat',status='unknown')
      read(800,*) E_nuc
      close(800)
      write(*,*)'READ IN NUC REPULSION'
!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------)

!--------------READ-INTEGRALS-NEEDED-FOR-NEO-HF------------------------(
      call read_nuc_ovlap(npbf,xxsp)
      write(*,*)
      write(*,*)'READ IN NUC OVLAP'
      call read_elec_ovlap(nebf,xxse)
      write(*,*)'READ IN ELEC OVLAP'
      call read_GAM_ecore(nebf,nebf2,GAM_ecore)
      write(*,*)'READ IN GAM_ECORE'
      call read_GAM_pcore(npbf,npbf2,GAM_pcore)
      write(*,*)'READ IN GAM_PCORE'
!      call read_GAM_ep(nebf,npbf,ng1,GAM_ep)
!      write(*,*)'READ IN GAM_EP'
      call read_GAM_ee(nebf,ngee,GAM_ee)
      write(*,*)'READ IN GAM_EE'
      write(*,*)

C store quantities over special electron basis
      call RXCHFmult_contr_mat(nebf,nebfBE,xxse,xxseBE)
      call RXCHFmult_contr_mat(nebf,nebfBE,GAM_ecore,GAM_ecoreBE)

!     write(*,*)
!     write(*,*)'IN: xcuscf '
!     write(*,*)'AFTER READ IN GAM_EE='
!     write(*,*)GAM_ee
!     write(*,*)

!--------------READ-INTEGRALS-NEEDED-FOR-NEO-HF------------------------)

!-------------INITIAL-GUESSES------------------------------------------(
      if(read_CE) then
!        READ IN GUESS FOR E:
!        call read_elec_density(nebf,nelec,DE)
         write(*,*)
         write(*,*) "Reading electronic orbitals"
         if (LALTBAS) then
          write(*,*)
          write(*,*) "  Regular orbitals will be reordered to include"
          write(*,*) "  special electronic basis first"
          write(*,*)
          write(*,*) "  Special orbitals MUST ALREADY BE in order!"
          write(*,*)
         end if
         call RXCHFmult_read_CAE(nebf,nebfBE,LALTBAS,elindBE,
     x                           NAE,DAE,VECAE0)
         call RXCHFmult_read_CBE(nebfBE,NBE,DBE,VECBE0)
         write(*,*) "Done reading electronic orbitals"
         write(*,*)
      else
!       STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
        write(*,*)'ABOUT TO CALL guess_A_elec'
!       call guess_elec(nelec,nebf,xxse,GAM_ecore,DE)
        if ((LOCBSE).or.(LOCBSE2)) then
          call RXCHFmult_guess_elec(LALTBAS,nae,nbe,nebf,nebfBE,
     x                              xxse,xxseBE,
     x                              GAM_ecore,GAM_ecoreBE,
     x                              DAE,DBE,VECAE0,VECBE0)
          write(*,*)'BACK FROM guess_elec for OCBSE'
        else
         call RXCHFmult_guess_A_elec(NAE,nebf,xxse,GAM_ecore,
     x                               DAE,VECAE0)
         call RXCHFmult_guess_A_elec(NBE,nebfBE,xxse,GAM_ecoreBE,
     x                               DBE,VECBE0)
          write(*,*)'BACK FROM guess_elec'
          write(*,*)
        end if
      end if
      if(read_CP) then
!        READ IN GUESS FOR N:
         call read_nuc_density(npbf,1,NUCST,DP)
      else
!        STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
         write(*,*)'ABOUT TO CALL guess_prot'
!        call guess_prot(NUCST,npbf,nebf,xxsp,GAM_pcore,GAM_ep,DE,DP)
         call guess_prot2(NUCST,npbf,xxsp,GAM_pcore,DP)
         write(*,*)'BACK FROM guess_prot'
         write(*,*)
      end if

C ARS( debug: print out initial guess MOs here
      AEe=0.0d+00
      BEe=0.0d+00
      if (LCMF) then
       write(*,*)
       write(*,*) "------------------"
       write(*,*) "INITIAL GUESS MOs:"
       write(*,*) "------------------"
       write(*,*)
       WRITE(*,9610)
       call PREVNU(vecAE0,AEE,nebf,nebf,nebf)
       WRITE(*,9620)
       call PREVNU(vecBE0,BEE,nebfBE,nebfBE,nebfBE)
C vecP not defined yet
C       WRITE(*,9700)
C       call PREVNU(vecp,EP,npbf,npbf,npbf)
       write(*,*)
      end if
C )
!-------------INITIAL-GUESSES------------------------------------------)

!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------(
      if (LSOSCF) then
         SOGTOL=0.40d+00
         SMALL=1.0D-06
         L0=nebf
         L1=nebf
         L0b=nebfBE
         L1b=nebfBE
         LSOSCFA=.true.
         LSOSCFB=.true.
         if((nae.eq.1).or.LOCBSE) LSOSCFA=.FALSE.
         if((nbe.eq.1).or.LOCBSE) LSOSCFB=.FALSE.
      else
         LSOSCFA=.false.
         LSOSCFB=.false.
      end if
      if(LSOSCFA) THEN
         NFT15=15
         OPEN(NFT15, FILE='WORK15', STATUS='UNKNOWN',
     *        ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         NA=nae/2
      end if
      if(LSOSCFB) THEN
         NFT16=16
         OPEN(NFT16, FILE='WORK16', STATUS='UNKNOWN',
     *        ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         NB=nbe/2
         if(.not.LALTBAS) then
! Allocate here (previously on stack)
           if(allocated(XB))      deallocate(XB)
           if(allocated(UPDTIB))  deallocate(UPDTIB)
           if(allocated(DGRADIB)) deallocate(DGRADIB)
           if(allocated(DISPLNB)) deallocate(DISPLNB)
           if(allocated(UPDTB))   deallocate(UPDTB)
           if(allocated(DISPLB))  deallocate(DISPLB)
           if(allocated(DGRADB))  deallocate(DGRADB)
           if(allocated(DISPLIB)) deallocate(DISPLIB)
           if(allocated(PGRADB))  deallocate(PGRADB)
           if(allocated(GRADB))   deallocate(GRADB)
           if(allocated(HSTARTB)) deallocate(HSTARTB)
           allocate(HSTARTB(NPRB))
           allocate(GRADB(NPRB))
           allocate(PGRADB(NPRB))
           allocate(DISPLIB(NPRB))
           allocate(DGRADB(NPRB))
           allocate(DISPLB(NPRB))
           allocate(UPDTB(NPRB))
           allocate(DISPLNB(NPRB))
           allocate(DGRADIB(NPRB))
           allocate(UPDTIB(NPRB))
           allocate(XB(NPRB))
         end if
      end if
!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------)

!
!     SET CONVERGENCE CRITERIA AND MAXIMUM ITERATIONS 
!
      TOLE = 1.0D-06
      TOLP = 1.0D-04
      maxit=100
      maxmicroit=200
      if(LOCBSE) maxit=400
!
!     ZERO OUT 'OLD' DENSITY MATRICES
!
      DAE0=0.0d+00
      DBE0=0.0d+00
      DP0=0.0d+00
!
!     BEGIN XCSCF ITERATIONS
      WRITE(*,9500)
!     WRITE(*,9000)
!
      E_total_old=0.0d+00

      DO I=1,MAXIT

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------(
!     write(*,*)
!     write(*,*)'IN: xcuscf '
!     write(*,*)'Before call to UHF_FOCK, GAM_EE='
!     write(*,*)GAM_ee
!     write(*,*)

C Call HF Fock build for NAE regular electrons
         call RXCHFmult_fock_hf(LCMF,nebf,nebf2,NAE,ngee,
     x                          DAE,GAM_ecore,GAM_ee,
     x                          FAE,E_HF,E_HF_ecore,E_HF_ee)

C Call XCHF Fock build for NBE special electrons and one QM particle
         call RXCHFmult_fock_xchf(LGAM4,LG4DSCF,LG4IC,
     x                   LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                   NG4CHK,NG3CHK,NG2CHK,
     x                   dimXCHF4,dimXCHF3,dimXCHF2,
     x                   npebfBE,nebfBE,nebfBE2,npbf,npbf2,NBE,
     x                   ngee,ngtg1,DBE,DP,
     x                   XCHF_GAM4,XCHF_GAM3,XCHF_GAM2,XCHF_GAM2s,
     x                   nat,pmass,cat,zan,bcoef1,gamma1,
     x                   KPESTR_be,KPEEND_be,
     x                   AMPEB2C_be,AGEBFCC_be,
     x                   ELCEX_be,ELCAM_be,ELCBFC_be,
     x                   AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                   FBE,FP,SBE_XCHF,SP_XCHF,
     x                   E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                   E_XCHF_gam3,E_XCHF_gam4,
     x                   S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
         call RXCHFmult_fock_int(LCMF,LADDEXCH,nelec,NAE,NBE,
     x                           nebf,nebfBE,npbf,
     x                           dimINT2,dimINT3,dimINT4,
     x                           dimINT2ex,dimINT3ex,
     x                           NG2CHK,NG3CHK,NG4CHK,
     x                           DAE,DBE,DP,
     x                           INT_GAM2,INT_GAM3,INT_GAM4,
     x                           INT_GAM2ex,INT_GAM3ex1,INT_GAM3ex2,
     x                           S_total,S_gam2,SBE_XCHF,SP_XCHF,
     x                           FPint,FAEint,FBEint, 
     x                           E_int_OMG2,E_int_OMG3,E_int_OMG4,
     x                           E_int)

C ARS( no interaction
      if(LNOINT) then
          E_int=0.0d+00
          E_HF=0.0d+00
      else
          call add2fock(npbf,FPint,FP)
          call add2fock(nebf,FAEint,FAE)
          call add2fock(nebfBE,FBEint,FBE)
      end if
C )

          IF (LCMF) then
           write(*,*)
           write(*,*) "FAE:"
           call prt_lower_triangle(nebf,nebflt,FAE)
           write(*,*)
           write(*,*) "FBE:"
           call prt_lower_triangle(nebfBE,nebfBElt,FBE)
           write(*,*)
           write(*,*) "FP:"
           call prt_lower_triangle(npbf,npbflt,FP)
           write(*,*)
          END IF  

          E_total=E_HF+E_XCHF+E_int+E_nuc

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------)
         if(I.eq.1) then
            WRITE(*,9400)

      WRITE(*,9300) E_nuc,E_HF_ecore,E_HF_ee,E_HF,
     x  E_XCHF_gam1,E_XCHF_gam2,E_XCHF_gam3,E_XCHF_gam4,E_XCHF,
     x  E_int_OMG2,E_int_OMG3,E_int_OMG4,E_int,
     x  S_total,E_total

         end if

!        Fockp diag
         call UROOTHAN(vecP,EP,xxsp,FP,npbf)
         call construct_DP(nucst,npbf,vecP,DP)

C ARS( reform elec Fock matrices

C ARS( microiterate
!        --> FIND LARGEST CHANGE IN P DENSITY
         CALL DENDIF(DP0,DP,NPBF,DIFFP)
         CALL COPYDEN(DP0,DP,NPBF)
C )

         if(LSOSCFA) ITSOA=0
         if(LSOSCFB) ITSOB=0
         ORBGRDA=0.0d+00
         ORBGRDB=0.0d+00
         PGRADA=0.0d+00
         if(allocated(PGRADB)) PGRADB=0.0D+00

         write(*,2001) I

         if((LSOSCFA).and.(LSOSCFB)) then 
          WRITE(*,9052)
         else if ((LSOSCFA).and.(.not.(LSOSCFB))) then
          WRITE(*,9050)
         else if ((LSOSCFB).and.(.not.(LSOSCFA))) then
          WRITE(*,9051)
         else
          WRITE(*,9000)
         end if

         do ielec=1,maxmicroit

       if((.not.((locbse).or.(locbse2))).or.(ielec.eq.1)) then
C Call HF Fock build for NAE regular electrons
           call RXCHFmult_fock_hf(LCMF,nebf,nebf2,NAE,ngee,
     x                            DAE,GAM_ecore,GAM_ee,
     x                            FAE,E_HF,E_HF_ecore,E_HF_ee)

C Call XCHF Fock build for NBE special electrons and one QM particle
           call RXCHFmult_fock_xchf(LGAM4,LG4DSCF,LG4IC,
     x                   LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                   NG4CHK,NG3CHK,NG2CHK,
     x                   dimXCHF4,dimXCHF3,dimXCHF2,
     x                   npebfBE,nebfBE,nebfBE2,npbf,npbf2,NBE,
     x                   ngee,ngtg1,DBE,DP,
     x                   XCHF_GAM4,XCHF_GAM3,XCHF_GAM2,XCHF_GAM2s,
     x                   nat,pmass,cat,zan,bcoef1,gamma1,
     x                   KPESTR_be,KPEEND_be,
     x                   AMPEB2C_be,AGEBFCC_be,
     x                   ELCEX_be,ELCAM_be,ELCBFC_be,
     x                   AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                   FBE,FP,SBE_XCHF,SP_XCHF,
     x                   E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                   E_XCHF_gam3,E_XCHF_gam4,
     x                   S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
           call RXCHFmult_fock_int(LCMF,LADDEXCH,nelec,NAE,NBE,
     x                             nebf,nebfBE,npbf,
     x                             dimINT2,dimINT3,dimINT4,
     x                             dimINT2ex,dimINT3ex,
     x                             NG2CHK,NG3CHK,NG4CHK,
     x                             DAE,DBE,DP,
     x                             INT_GAM2,INT_GAM3,INT_GAM4,
     x                             INT_GAM2ex,INT_GAM3ex1,INT_GAM3ex2,
     x                             S_total,S_gam2,SBE_XCHF,SP_XCHF,
     x                             FPint,FAEint,FBEint, 
     x                             E_int_OMG2,E_int_OMG3,E_int_OMG4,
     x                             E_int)

C ARS( no interaction
      if(LNOINT) then
            E_int=0.0d+00
            E_HF=0.0d+00
      else
            call add2fock(npbf,FPint,FP)
            call add2fock(nebf,FAEint,FAE)
            call add2fock(nebfBE,FBEint,FBE)
      end if
C )
C ARS( microiteration
          E_total=E_HF+E_XCHF+E_int+E_nuc
C )
      end if
C )

         if (LOCBSE) then  ! hardwired off in place of OCBSE2
! Do OCBSE procedure (restricted solutions for regular and special electrons)

           call RXCHFmult_OCBSE(nebf,nae,nbe,vecAE0,vecBE0,FAE,FBE,xxse,
     x                          vecAE,vecBE,AEe,BEe)

! Form regular electronic density matrix and store stuff for next it
           call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)
           CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
           CALL COPYDEN(DAE0,DAE,NEBF)
           CALL COPYDEN(vecAE0,vecAE,NEBF)

! Form special electronic density matrix and store stuff for next it
           call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)
           CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
           CALL COPYDEN(DBE0,DBE,NEBF)
           CALL COPYDEN(vecBE0,vecBE,NEBF)

! Calculate energy for this it and Fock matrices for next it

C Call HF Fock build for NAE regular electrons
           call RXCHFmult_fock_hf(LCMF,nebf,nebf2,NAE,ngee,
     x                            DAE,GAM_ecore,GAM_ee,
     x                            FAE,E_HF,E_HF_ecore,E_HF_ee)

C Call XCHF Fock build for NBE special electrons and one QM particle
           call RXCHFmult_fock_xchf(LGAM4,LG4DSCF,LG4IC,
     x                   LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                   NG4CHK,NG3CHK,NG2CHK,
     x                   dimXCHF4,dimXCHF3,dimXCHF2,
     x                   npebfBE,nebfBE,nebfBE2,npbf,npbf2,NBE,
     x                   ngee,ngtg1,DBE,DP,
     x                   XCHF_GAM4,XCHF_GAM3,XCHF_GAM2,XCHF_GAM2s,
     x                   nat,pmass,cat,zan,bcoef1,gamma1,
     x                   KPESTR_be,KPEEND_be,
     x                   AMPEB2C_be,AGEBFCC_be,
     x                   ELCEX_be,ELCAM_be,ELCBFC_be,
     x                   AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                   FBE,FP,SBE_XCHF,SP_XCHF,
     x                   E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                   E_XCHF_gam3,E_XCHF_gam4,
     x                   S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
           call RXCHFmult_fock_int(LCMF,LADDEXCH,nelec,NAE,NBE,
     x                             nebf,nebfBE,npbf,
     x                             dimINT2,dimINT3,dimINT4,
     x                             dimINT2ex,dimINT3ex,
     x                             NG2CHK,NG3CHK,NG4CHK,
     x                             DAE,DBE,DP,
     x                             INT_GAM2,INT_GAM3,INT_GAM4,
     x                             INT_GAM2ex,INT_GAM3ex1,INT_GAM3ex2,
     x                             S_total,S_gam2,SBE_XCHF,SP_XCHF,
     x                             FPint,FAEint,FBEint, 
     x                             E_int_OMG2,E_int_OMG3,E_int_OMG4,
     x                             E_int)

C ARS( no interaction
      if(LNOINT) then
            E_int=0.0d+00
            E_HF=0.0d+00
      else
            call add2fock(npbf,FPint,FP)
            call add2fock(nebf,FAEint,FAE)
            call add2fock(nebfBE,FBEint,FBE)
      end if
C )

            IF (LCMF) then
             write(*,*)
             write(*,*) "FAE:"
             call prt_lower_triangle(nebf,nebflt,FAE)
             write(*,*)
             write(*,*) "FBE:"
             call prt_lower_triangle(nebf,nebflt,FBE)
             write(*,*)
             write(*,*) "FP:"
             call prt_lower_triangle(npbf,npbflt,FP)
             write(*,*)
            END IF  

            E_total=E_HF+E_XCHF+E_int+E_nuc

         else if (LOCBSE2) then
! Do OCBSE2 procedure (restricted solutions for special electrons)

! Regular electrons
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------(
         if(LSOSCFA) THEN
          ITER=IELEC
          EIGAVL = ITER.GT.1
         end if
         IF(LSOSCFA .AND.  EIGAVL) THEN                ! first it. skip SOSCF (diag to get EE)
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAE,FLT)
          call SOGRAD(GRADA,FLT,vecAE,WRK,NPRA,NA,L0,L1,NEBFLT,ORBGRDA)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 900  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN   ! only on first SOSCF it. set up approx Hess
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AEE,NPRA,L0,NA,NA)
              END IF
              ITSOA = ITSOA+1
           call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,DISPLA,UPDTA,
     *                 DISPLNA,DGRADIA,UPDTIA,ORBGRDA,NPRA,ITSOA,NFT15)
            call SOTRAN(DISPLIA,vecAE,GA,WRK,NPRA,L0,L1,NA,NA,ORBGRDA)
             CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)
              GO TO 950  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------)

  900 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAE,AEE,xxse,FAE,nebf)
         call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)

  950 CONTINUE
!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
         CALL COPYDEN(DAE0,DAE,NEBF)

C ARS( new stuff

         if(LALTBAS) then

          call RXCHFmult_intersection(nebf,nebf-nocca,
     x                                vecAE(:,nocca+1:nebf),nebfBE,
     x                                xxseBE,dimint,Ctemp)

          if (dimint.le.noccb) then
           write(*,*) "dim of intersection <= # occ special orbitals"
           write(*,*) "dimint:",dimint
           write(*,*) "noccb:",noccb
           return
          end if

          if((i.eq.1).and.(ielec.eq.1)) then
           write(*,*)
           write(*,*) "----------------------------------"
           write(*,*) " Dimension of intersection space:"
           write(*,'(2X,A,1X,I3)') "Max possible (nebfBE)      =",nebfBE
           write(*,'(2X,A,1X,I3)') "Actual (after computation) =",dimint
           write(*,*) "----------------------------------"
           write(*,*)
          end if

          if(dimint.ne.dimint0) then

! Warn about changing dimension if not first iteration
           if(.not.((i.eq.1).and.(ielec.eq.1))) then
            write(*,*)
            write(*,*) "----------------------------------"
            write(*,*) " Dimension of intersection space  "
            write(*,*) " has CHANGED (old,new):           "
            write(*,*) dimint0,dimint
            if(LSOSCFB) write(*,*) " Resetting SOSCF...               "
            write(*,*) "----------------------------------"
            write(*,*)
           end if

! Allocate OCBSE variables
           nwbf=dimint
           nwbflt=nwbf*(nwbf+1)/2
           L0w=nwbf
           L1w=nwbf
           if(allocated(Cint)) deallocate(Cint)
           if(allocated(wWRKw)) deallocate(wWRKw)
           if(allocated(wGBw)) deallocate(wGBw)
           if(allocated(wFLTw)) deallocate(wFLTw)
           if(allocated(wBEenw)) deallocate(wBEenw)
           if(allocated(wvecBEw)) deallocate(wvecBEw)
           if(allocated(wFBEw)) deallocate(wFBEw)
           if(allocated(WB)) deallocate(WB)
           allocate(WB(nebf,nwbf))
           allocate(wFBEw(nwbf,nwbf))
           allocate(wvecBEw(nwbf,nwbf))
           allocate(wBEenw(nwbf))
           allocate(wFLTw(nwbflt))
           allocate(wGBw(nwbf,nwbf))
           allocate(wWRKw(nwbf))
           allocate(Cint(nebfBE,dimint))

! Allocate SOSCF variables
           if(LSOSCFB) then
            NPRB=(dimint-noccb)*noccb
            if(allocated(XB))      deallocate(XB)
            if(allocated(UPDTIB))  deallocate(UPDTIB)
            if(allocated(DGRADIB)) deallocate(DGRADIB)
            if(allocated(DISPLNB)) deallocate(DISPLNB)
            if(allocated(UPDTB))   deallocate(UPDTB)
            if(allocated(DISPLB))  deallocate(DISPLB)
            if(allocated(DGRADB))  deallocate(DGRADB)
            if(allocated(DISPLIB)) deallocate(DISPLIB)
            if(allocated(PGRADB))  deallocate(PGRADB)
            if(allocated(GRADB))   deallocate(GRADB)
            if(allocated(HSTARTB)) deallocate(HSTARTB)
            allocate(HSTARTB(NPRB))
            allocate(GRADB(NPRB))
            allocate(PGRADB(NPRB))
            allocate(DISPLIB(NPRB))
            allocate(DGRADB(NPRB))
            allocate(DISPLB(NPRB))
            allocate(UPDTB(NPRB))
            allocate(DISPLNB(NPRB))
            allocate(DGRADIB(NPRB))
            allocate(UPDTIB(NPRB))
            allocate(XB(NPRB))
            ITSOB=0
            ORBGRDB=0.0d+00
            PGRADB=0.0d+00
           end if

          end if

          do ialt=1,dimint
            do jalt=1,nebfBE
              Cint(jalt,ialt)=Ctemp(jalt,ialt)
            end do
          end do

! Debug: All CBE vectors should be orthogonal to occ CAE vectors
          do ialt=1,nocca
          do jalt=1,dimint
            ovlapalt=zero
            do kalt=1,nebf
            do lalt=1,nebfBE
              ovlapalt=ovlapalt+vecAE(kalt,ialt)*Cint(lalt,jalt)
     x                         *xxse(kalt,lalt)
            end do
            end do
            if (abs(ovlapalt).gt.tolalt) then
             write(*,*)
             write(*,*) "******* ERROR *******"
             write(*,*) "Calculated intersection basis not orthogonal"
             write(*,*) "to occupied regular vectors"
             write(*,*) "reg occ index, int index, ovlap:",
     x                  ialt,jalt,ovlapalt
             write(*,*)
            end if
          end do
          end do

          call RXCHFmult_OCBSE_transF(LCMF,nebfBE,nwbf,Cint,
     x                                FBE,WB,wFBEw)

         else

! Transform FBE (calculated at end of previous iteration) to new W basis
!  - W updated with new vecA from this iteration
!  - vecBE in AO basis from previous iteration is transformed to new W basis
!    (relevant for SOSCF only)
          call RXCHFmult_OCBSE_transF(LCMF,nebf,nwbf,
     x                               vecAE(:,nocca+1:nebf),FBE,WB,wFBEw)

         end if

!-----------------------POSSIBLE-SOSCF-BETA----------------------------(
         if(LSOSCFB) THEN
           ITER=IELEC
           EIGAVL = ITER.GT.1
! Turn off SOSCF if dimint changes (until stable again)
           if((LALTBAS).and.(dimint.ne.dimint0)) EIGAVL=.false.
         end if
         IF(LSOSCFB .AND.  EIGAVL) THEN                ! first it. skip SOSCF (diag to get EE)
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nwbf,nwbfLT,wFBEw,wFLTw)
! Transform vecBE that was used to build FBE into new W basis
           call RXCHFmult_OCBSE_transVt(LCMF,nebfBE,nwbf,WB,
     x                                  xxseBE,vecBE,wvecBEw)
           call SOGRAD(GRADB,wFLTw,wvecBEw,wWRKw,NPRB,NB,
     x                 L0w,L1w,NwBFLT,ORBGRDB)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
              IF(ITSOB.EQ.0) THEN   ! only on first SOSCF it. set up approx Hess
                 WRITE(*,9800)
                 call SOHESS(HSTARTB,wBEenw,NPRB,L0w,NB,NB)
              END IF
              ITSOB = ITSOB+1
           call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,DISPLB,UPDTB,
     *                 DISPLNB,DGRADIB,UPDTIB,ORBGRDB,NPRB,ITSOB,NFT16)
              call SOTRAN(DISPLIB,wvecBEw,wGBw,wWRKw,NPRB,
     x                    L0w,L1w,NB,NB,ORBGRDB)
              CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
              call RXCHFmult_OCBSE_transV(LCMF,nebfBE,nwbf,WB,wvecBEw,! eigenvalues useless
     x                                    wBEenw,vecBE,BEe)
              call RXCHFmult_construct_DE(NBE,nebfBE,vecBE,DBE)
              GO TO 450  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-BETA----------------------------)

! No SOSCF
!  - Diagonalize Fock matrix in W basis of this iteration
!  - Obtain updated vecBE in W basis of this iteration
         call RXCHFmult_OCBSE_diag(LCMF,nebfBE,nwbf,WB,wFBEw,
     x                             wvecBEw,wBEenw,vecBE,BEe)
         call RXCHFmult_construct_DE(NBE,nebfBE,vecBE,DBE)

  450 CONTINUE

         CALL DENDIF(DBE0,DBE,NEBFBE,DIFFBE)
         CALL COPYDEN(DBE0,DBE,NEBFBE)

C )

! Calculate energy for this it and Fock matrices for next it

C Call HF Fock build for NAE regular electrons
           call RXCHFmult_fock_hf(LCMF,nebf,nebf2,NAE,ngee,
     x                            DAE,GAM_ecore,GAM_ee,
     x                            FAE,E_HF,E_HF_ecore,E_HF_ee)

C Call XCHF Fock build for NBE special electrons and one QM particle
           call RXCHFmult_fock_xchf(LGAM4,LG4DSCF,LG4IC,
     x                   LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                   NG4CHK,NG3CHK,NG2CHK,
     x                   dimXCHF4,dimXCHF3,dimXCHF2,
     x                   npebfBE,nebfBE,nebfBE2,npbf,npbf2,NBE,
     x                   ngee,ngtg1,DBE,DP,
     x                   XCHF_GAM4,XCHF_GAM3,XCHF_GAM2,XCHF_GAM2s,
     x                   nat,pmass,cat,zan,bcoef1,gamma1,
     x                   KPESTR_be,KPEEND_be,
     x                   AMPEB2C_be,AGEBFCC_be,
     x                   ELCEX_be,ELCAM_be,ELCBFC_be,
     x                   AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                   FBE,FP,SBE_XCHF,SP_XCHF,
     x                   E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                   E_XCHF_gam3,E_XCHF_gam4,
     x                   S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
           call RXCHFmult_fock_int(LCMF,LADDEXCH,nelec,NAE,NBE,
     x                             nebf,nebfBE,npbf,
     x                             dimINT2,dimINT3,dimINT4,
     x                             dimINT2ex,dimINT3ex,
     x                             NG2CHK,NG3CHK,NG4CHK,
     x                             DAE,DBE,DP,
     x                             INT_GAM2,INT_GAM3,INT_GAM4,
     x                             INT_GAM2ex,INT_GAM3ex1,INT_GAM3ex2,
     x                             S_total,S_gam2,SBE_XCHF,SP_XCHF,
     x                             FPint,FAEint,FBEint, 
     x                             E_int_OMG2,E_int_OMG3,E_int_OMG4,
     x                             E_int)

C ARS( no interaction
      if(LNOINT) then
            E_int=0.0d+00
            E_HF=0.0d+00
      else
            call add2fock(npbf,FPint,FP)
            call add2fock(nebf,FAEint,FAE)
            call add2fock(nebfBE,FBEint,FBE)
      end if
C )

            IF (LCMF) then
             write(*,*)
             write(*,*) "FAE:"
             call prt_lower_triangle(nebf,nebflt,FAE)
             write(*,*)
             write(*,*) "FBE:"
             call prt_lower_triangle(nebfBE,nebfBElt,FBE)
             write(*,*)
             write(*,*) "FP:"
             call prt_lower_triangle(npbf,npbflt,FP)
             write(*,*)
            END IF  

            E_total=E_HF+E_XCHF+E_int+E_nuc


            if(LALTBAS) dimint0=dimint

         else

!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------(
         if(LSOSCFA) THEN
          ITER=IELEC
          EIGAVL = ITER.GT.1
         end if
         IF(LSOSCFA .AND.  EIGAVL) THEN                ! first it. skip SOSCF (diag to get EE)
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebf,nebfLT,FAE,FLT)
          call SOGRAD(GRADA,FLT,vecAE,WRK,NPRA,NA,L0,L1,NEBFLT,ORBGRDA)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDA.LT.SOGTOL  .OR.  ITSOA.GT.0) THEN
              IF(ITSOA.EQ.0) THEN   ! only on first SOSCF it. set up approx Hess
              WRITE(*,9800)
                 call SOHESS(HSTARTA,AEE,NPRA,L0,NA,NA)
              END IF
              ITSOA = ITSOA+1
           call SONEWT(HSTARTA,GRADA,PGRADA,DISPLIA,DGRADA,DISPLA,UPDTA,
     *                 DISPLNA,DGRADIA,UPDTIA,ORBGRDA,NPRA,ITSOA,NFT15)
            call SOTRAN(DISPLIA,vecAE,GA,WRK,NPRA,L0,L1,NA,NA,ORBGRDA)
             CALL DCOPY(NPRA,GRADA,1,PGRADA,1)
              call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)
              GO TO 750  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-ALPHA---------------------------)

  700 CONTINUE
!        Diagonalize Electronic Fock Matrices
!        call ROOTHAN(DAE,vecAE,AEE,xxse,FAE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecAE,AEE,xxse,FAE,nebf)
         call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)

  750 CONTINUE
!        --> FIND LARGEST CHANGE IN Alpha E DENSITY
         CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
         CALL COPYDEN(DAE0,DAE,NEBF)

!-----------------------POSSIBLE-SOSCF-BETA----------------------------(
        if(LSOSCFB) THEN
         ITER=IELEC
         EIGAVL = ITER.GT.1
        end if
         IF(LSOSCFB .AND.  EIGAVL) THEN                ! first it. skip SOSCF (diag to get EE)
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nebfBE,nebfBELT,FBE,FLTB)
          call SOGRAD(GRADB,FLTB,vecBE,WRKB,NPRB,NB,
     x                L0b,L1b,NEBFBELT,ORBGRDB)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
              IF(ITSOB.EQ.0) THEN   ! only on first SOSCF it. set up approx Hess
             WRITE(*,9800)
                 call SOHESS(HSTARTB,BEE,NPRB,L0b,NB,NB)
              END IF
              ITSOB = ITSOB+1
           call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,DISPLB,UPDTB,
     *                 DISPLNB,DGRADIB,UPDTIB,ORBGRDB,NPRB,ITSOB,NFT16)
            call SOTRAN(DISPLIB,vecBE,GB,WRKB,NPRB,
     x                  L0b,L1b,NB,NB,ORBGRDB)
             CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
              call RXCHFmult_construct_DE(NBE,nebfBE,vecBE,DBE)
              GO TO 850  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-BETA----------------------------)

  800 CONTINUE
!        call ROOTHAN(DBE,vecBE,BEE,xxse,FBE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecBE,BEE,xxseBE,FBE,nebfBE)
         call RXCHFmult_construct_DE(NBE,nebfBE,vecBE,DBE)

  850 CONTINUE
!        --> FIND LARGEST CHANGE IN Beta E DENSITY
         CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
         CALL COPYDEN(DBE0,DBE,NEBF)

         end if ! end if for not ocbse or ocbse2

C ARS( microiterate
C!        --> FIND LARGEST CHANGE IN P DENSITY
C         CALL DENDIF(DP0,DP,NPBF,DIFFP)
C         CALL COPYDEN(DP0,DP,NPBF)
C )

!        --> CALCULATE CHANGE IN TOTAL ENERGY
         Delta_E_tot=E_total-E_total_old
         E_total_old=E_total

!        --> PRINT SUMMARY OF THIS ITERATION
         if((LSOSCFA).and.(LSOSCFB)) then
            WRITE(*,9151) IELEC,E_total,Delta_E_tot,
     x                    DIFFAE,DIFFBE,DIFFP,ORBGRDA,ORBGRDB
         else if ((LSOSCFA).and.(.not.(LSOSCFB))) then
            WRITE(*,9150) IELEC,E_total,Delta_E_tot,
     x                    DIFFAE,DIFFBE,DIFFP,ORBGRDA
         else if ((LSOSCFB).and.(.not.(LSOSCFA))) then
            WRITE(*,9150) IELEC,E_total,Delta_E_tot,
     x                    DIFFAE,DIFFBE,DIFFP,ORBGRDB
         else
            WRITE(*,9100) IELEC,E_total,Delta_E_tot,
     x                    DIFFAE,DIFFBE,DIFFP
         end if
C ARS( debug: print out MOs here
      if (LCMF) then
       WRITE(*,9610)
       call PREVNU(vecAE,AEE,nebf,nebf,nebf)
       WRITE(*,9620)
       call PREVNU(vecBE,BEE,nebfBE,nebfBE,nebfBE)
       WRITE(*,9700)
       call PREVNU(vecp,EP,npbf,npbf,npbf)
      end if
C )
! Output the vectors for this iteration for restart if necessary:
         call write_MOs(860,nebf,VECAE)
         call write_MOs(861,nebfBE,VECBE)
         call write_MOs(853,npbf,VECP)

         LDIFFE=( (DIFFAE.LT.TOLE).and.(DIFFBE.LT.TOLE) )
         IF(LDIFFE) GOTO 200
         IF(IELEC.EQ.MAXIT) GOTO 10

       END DO  ! microiterations

  200 CONTINUE
!      IF WE GET HERE - MICROITERATION CONVERGENCE ACHIEVED
       write(*,2000) i,ielec

       IF(DIFFP.LT.TOLP) GOTO 100
       IF(I.EQ.MAXIT) GOTO 10

      END DO   ! iterations
 
  10  CONTINUE
!     IF WE GET HERE SOMETHING WENT WRONG

      if(LSOSCFA) close(NFT15)
      if(LSOSCFB) close(NFT16)

      WRITE(*,*)
      WRITE(*,*)'WARNING:  ITERATION LIMIT EXCEEDED'
      E_total=zero
      WRITE(*,*)
!     STOP
!
  100 CONTINUE
!     IF WE GET HERE WE ARE DONE - CONVERGENCE ACHIEVED

      if(LSOSCFA) close(NFT15)
      if(LSOSCFB) close(NFT16)

C ARS( test convergence
C      if(LSOSCFB) then
C
C! Rebuild Fock matrices
CC Call HF Fock build for NAE regular electrons
C         call RXCHFmult_fock_hf(LCMF,nebf,nebf2,NAE,ngee,
C     x                          DAE,GAM_ecore,GAM_ee,
C     x                          FAE,E_HF,E_HF_ecore,E_HF_ee)
C
CC Call XCHF Fock build for NBE special electrons and one QM particle
C         call RXCHFmult_fock_xchf(LGAM4,LG4DSCF,LG4IC,
C     x                   LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
C     x                   NG4CHK,NG3CHK,NG2CHK,
C     x                   dimXCHF4,dimXCHF3,dimXCHF2,
C     x                   npebfBE,nebfBE,nebfBE2,npbf,npbf2,NBE,
C     x                   ngee,ngtg1,DBE,DP,
C     x                   XCHF_GAM4,XCHF_GAM3,XCHF_GAM2,XCHF_GAM2s,
C     x                   nat,pmass,cat,zan,bcoef1,gamma1,
C     x                   KPESTR_be,KPEEND_be,
C     x                   AMPEB2C_be,AGEBFCC_be,
C     x                   ELCEX_be,ELCAM_be,ELCBFC_be,
C     x                   AGNBFCC,NUCEX,NUCAM,NUCBFC,
C     x                   FBE,FP,SBE_XCHF,SP_XCHF,
C     x                   E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
C     x                   E_XCHF_gam3,E_XCHF_gam4,
C     x                   S_total,S_gam1,S_gam2)
CC Call interaction Fock build for all particles
C         call RXCHFmult_fock_int(LCMF,LADDEXCH,nelec,NAE,NBE,
C     x                           nebf,nebfBE,npbf,
C     x                           dimINT2,dimINT3,dimINT4,
C     x                           dimINT2ex,dimINT3ex,
C     x                           NG2CHK,NG3CHK,NG4CHK,
C     x                           DAE,DBE,DP,
C     x                           INT_GAM2,INT_GAM3,INT_GAM4,
C     x                           INT_GAM2ex,INT_GAM3ex1,INT_GAM3ex2,
C     x                           S_total,S_gam2,SBE_XCHF,SP_XCHF,
C     x                           FPint,FAEint,FBEint, 
C     x                           E_int_OMG2,E_int_OMG3,E_int_OMG4,
C     x                           E_int)
C
CC ARS( no interaction
C      if(LNOINT) then
C          E_int=0.0d+00
C          E_HF=0.0d+00
C      else
C          call add2fock(npbf,FPint,FP)
C          call add2fock(nebf,FAEint,FAE)
C          call add2fock(nebfBE,FBEint,FBE)
C      end if
CC )
C
C          IF (LCMF) then
C           write(*,*)
C           write(*,*) "FAE:"
C           call prt_lower_triangle(nebf,nebflt,FAE)
C           write(*,*)
C           write(*,*) "FBE:"
C           call prt_lower_triangle(nebf,nebflt,FBE)
C           write(*,*)
C           write(*,*) "FP:"
C           call prt_lower_triangle(npbf,npbflt,FP)
C           write(*,*)
C          END IF  
C
C          E_total=E_HF+E_XCHF+E_int+E_nuc
C
C! Calculate updated Fock matrix in MO basis
C      call fock2mobasis(nebf,FBE,vecBE,FBEmo)
C      write(*,*)
C      write(*,*) "FBE in MO basis:"
C      write(*,*)
C      call PREVNU(FBEmo,BEE,nebf,nebf,nebf)
C      CALL COPYDEN(vecBE0,vecBE,NEBF)
C
C! Diagonalize in MO basis
C      CALL RS(nebf,nebf,FBEmo,BEE,2,vecBE,work1,work2,IERR)
C      call fock2mobasis(nebf,FBEmo,vecBE,FBEmo)
C      write(*,*)
C      write(*,*) "FBE after diagonalization in MO basis:"
C      write(*,*)
C      call PREVNU(FBEmo,BEE,nebf,nebf,nebf)
C
C! Output density matrix
C      write(*,*)
C      write(*,*) "DBE before diagonalization:"
C      write(*,*)
C      call PREVNU(DBE,BEE,nebf,nebf,nebf)
C
CC! Diagonalize only nuclear part without reforming
C         call UROOTHAN(vecP,EP,xxsp,FP,npbf)
C         call construct_DP(nucst,npbf,vecP,DP)
C
CC! Reform Fock matrices before diagonalizing electronic parts
CCC Call HF Fock build for NAE regular electrons
CC         call RXCHFmult_fock_hf(LCMF,nebf,nebf2,NAE,ngee,
CC     x                          DAE,GAM_ecore,GAM_ee,
CC     x                          FAE,E_HF,E_HF_ecore,E_HF_ee)
CC
CCC Call XCHF Fock build for NBE special electrons and one QM particle
CC         call RXCHFmult_fock_xchf(LGAM4,LG4DSCF,LG4IC,
CC     x                   LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
CC     x                   NG4CHK,NG3CHK,NG2CHK,
CC     x                   dimXCHF4,dimXCHF3,dimXCHF2,
CC     x                   npebfBE,nebfBE,nebfBE2,npbf,npbf2,NBE,
CC     x                   ngee,ngtg1,DBE,DP,
CC     x                   XCHF_GAM4,XCHF_GAM3,XCHF_GAM2,XCHF_GAM2s,
CC     x                   nat,pmass,cat,zan,bcoef1,gamma1,
CC     x                   KPESTR_be,KPEEND_be,
CC     x                   AMPEB2C_be,AGEBFCC_be,
CC     x                   ELCEX_be,ELCAM_be,ELCBFC_be,
CC     x                   AGNBFCC,NUCEX,NUCAM,NUCBFC,
CC     x                   FBE,FP,SBE_XCHF,SP_XCHF,
CC     x                   E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
CC     x                   E_XCHF_gam3,E_XCHF_gam4,
CC     x                   S_total,S_gam1,S_gam2)
CCC Call interaction Fock build for all particles
CC         call RXCHFmult_fock_int(LCMF,LADDEXCH,nelec,NAE,NBE,
CC     x                           nebf,nebfBE,npbf,
CC     x                           dimINT2,dimINT3,dimINT4,
CC     x                           dimINT2ex,dimINT3ex,
CC     x                           NG2CHK,NG3CHK,NG4CHK,
CC     x                           DAE,DBE,DP,
CC     x                           INT_GAM2,INT_GAM3,INT_GAM4,
CC     x                           INT_GAM2ex,INT_GAM3ex1,INT_GAM3ex2,
CC     x                           S_total,S_gam2,SBE_XCHF,SP_XCHF,
CC     x                           FPint,FAEint,FBEint, 
CC     x                           E_int_OMG2,E_int_OMG3,E_int_OMG4,
CC     x                           E_int)
CC
CCC ARS( no interaction
CC      if (LNOINT) then
CC          E_int=0.0d+00
CC          E_HF=0.0d+00
CC      else
CC          call add2fock(npbf,FPint,FP)
CC          call add2fock(nebf,FAEint,FAE)
CC          call add2fock(nebfBE,FBEint,FBE)
CC      end if
CCC )
CC
CC          IF (LCMF) then
CC           write(*,*)
CC           write(*,*) "FAE:"
CC           call prt_lower_triangle(nebf,nebflt,FAE)
CC           write(*,*)
CC           write(*,*) "FBE:"
CC           call prt_lower_triangle(nebf,nebflt,FBE)
CC           write(*,*)
CC           write(*,*) "FP:"
CC           call prt_lower_triangle(npbf,npbflt,FP)
CC           write(*,*)
CC          END IF  
CC
CC          E_total=E_HF+E_XCHF+E_int+E_nuc
C
CC! Diagonalize post-FP-diagonalized electronic Fock matrices
C         call UROOTHAN(vecAE,AEE,xxse,FAE,nebf)
C         call RXCHFmult_construct_DE(NAE,nebf,vecAE,DAE)
C
C         call UROOTHAN(vecBE,BEE,xxse,FBE,nebf)
C         call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)
C
C! Calculate updated Fock matrix in MO basis
C      call fock2mobasis(nebf,FBE,vecBE,FBEmo)
C      write(*,*)
C      write(*,*) "FBE in MO basis before rebuild:"
C      write(*,*)
C      call PREVNU(FBEmo,BEE,nebf,nebf,nebf)
C      CALL COPYDEN(vecBE0,vecBE,NEBF)
C
C! Output density matrix
C      write(*,*)
C      write(*,*) "DBE after diagonalization:"
C      write(*,*)
C      call PREVNU(DBE,BEE,nebf,nebf,nebf)
C
C
C! Output energy and density change
C         Delta_E_tot=E_total-E_total_old
C         E_total_old=E_total
C
C         CALL DENDIF(DP0,DP,NPBF,DIFFP)
C         CALL COPYDEN(DP0,DP,NPBF)
C         CALL DENDIF(DAE0,DAE,NEBF,DIFFAE)
C         CALL COPYDEN(DAE0,DAE,NEBF)
C         CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
C         CALL COPYDEN(DBE0,DBE,NEBF)
C
C         write(*,*) "After diag without reforming after FP diag:"
C         write(*,*) "E_total:        ",E_total
C         write(*,*) "Delta_E_total:  ",Delta_E_tot
C         write(*,*) "Max DAE change: ",DIFFAE
C         write(*,*) "Max DBE change: ",DIFFBE
C         write(*,*) "Max DP change: ",DIFFP
C
C! One more build and diagonalization
CC Call XCHF Fock build for NBE special electrons and one QM particle
C         call RXCHFmult_fock_xchf(LGAM4,LG4DSCF,LG4IC,
C     x                   LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
C     x                   NG4CHK,NG3CHK,NG2CHK,
C     x                   dimXCHF4,dimXCHF3,dimXCHF2,
C     x                   npebfBE,nebfBE,nebfBE2,npbf,npbf2,NBE,
C     x                   ngee,ngtg1,DBE,DP,
C     x                   XCHF_GAM4,XCHF_GAM3,XCHF_GAM2,XCHF_GAM2s,
C     x                   nat,pmass,cat,zan,bcoef1,gamma1,
C     x                   KPESTR_be,KPEEND_be,
C     x                   AMPEB2C_be,AGEBFCC_be,
C     x                   ELCEX_be,ELCAM_be,ELCBFC_be,
C     x                   AGNBFCC,NUCEX,NUCAM,NUCBFC,
C     x                   FBE,FP,SBE_XCHF,SP_XCHF,
C     x                   E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
C     x                   E_XCHF_gam3,E_XCHF_gam4,
C     x                   S_total,S_gam1,S_gam2)
C          E_total=E_HF+E_XCHF+E_int+E_nuc
C
C! Calculate updated Fock matrix in (new) MO basis
C      call fock2mobasis(nebf,FBE,vecBE,FBEmo)
C      write(*,*)
C      write(*,*) "FBE in MO basis:"
C      write(*,*)
C      call PREVNU(FBEmo,BEE,nebf,nebf,nebf)
C
C! Calculate updated Fock matrix in (old) MO basis
C      call fock2mobasis(nebf,FBE,vecBE0,FBEmo)
C      write(*,*)
C      write(*,*) "FBE in old MO basis:"
C      write(*,*)
C      call PREVNU(FBEmo,BEE,nebf,nebf,nebf)
C
C! Output density matrix
C      write(*,*)
C      write(*,*) "DBE before diagonalization:"
C      write(*,*)
C      call PREVNU(DBE,BEE,nebf,nebf,nebf)
C
C! Diagonalize
C      call UROOTHAN(vecBE,BEE,xxse,FBE,nebf)
C      call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)
C
C! Output density matrix
C      write(*,*)
C      write(*,*) "DBE after diagonalization:"
C      write(*,*)
C      call PREVNU(DBE,BEE,nebf,nebf,nebf)
C
C! Output energy and density change
C      Delta_E_tot=E_total-E_total_old
C      E_total_old=E_total
C
C      CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
C      CALL COPYDEN(DBE0,DBE,NEBF)
C
C      write(*,*) "After last FBE diag:"
C      write(*,*) "E_total:        ",E_total
C      write(*,*) "Delta_E_total:  ",Delta_E_tot
C      write(*,*) "Max DBE change: ",DIFFBE
C
C      end if
C )

!     PRINT FINAL ENERGY AND PUNCH THE ORBITALS
      WRITE(*,9200) E_total,I

      WRITE(*,9300) E_nuc,E_HF_ecore,E_HF_ee,E_HF,
     x  E_XCHF_gam1,E_XCHF_gam2,E_XCHF_gam3,E_XCHF_gam4,E_XCHF,
     x  E_int_OMG2,E_int_OMG3,E_int_OMG4,E_int,
     x  S_total,E_total

!  OUTPUT ELEC AND NUC EIGENVALUES AND EIGENVECTORS
      WRITE(*,9610)
      call PREVNU(vecAE,AEE,nebf,nebf,nebf)
      WRITE(*,9620)
      call PREVNU(vecBE,BEE,nebfBE,nebfBE,nebfBE)
      WRITE(*,9700)
      call PREVNU(vecp,EP,npbf,npbf,npbf)

! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------(
!     subroutine write_MOs(IFIL,nbf,VEC)
!     IFIL=852 :: FinalCE.dat
!     IFIL=853 :: FinalCP.dat
!     IFIL=860 :: FinalCAE.dat
!     IFIL=861 :: FinalCBE.dat
      call write_MOs(860,nebf,VECAE)
      call write_MOs(861,nebfBE,VECBE)
      call write_MOs(853,npbf,VECP)
! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------)
!
      if(LSOSCFB) then
       if(allocated(XB))      deallocate(XB)
       if(allocated(UPDTIB))  deallocate(UPDTIB)
       if(allocated(DGRADIB)) deallocate(DGRADIB)
       if(allocated(DISPLNB)) deallocate(DISPLNB)
       if(allocated(UPDTB))   deallocate(UPDTB)
       if(allocated(DISPLB))  deallocate(DISPLB)
       if(allocated(DGRADB))  deallocate(DGRADB)
       if(allocated(DISPLIB)) deallocate(DISPLIB)
       if(allocated(PGRADB))  deallocate(PGRADB)
       if(allocated(GRADB))   deallocate(GRADB)
       if(allocated(HSTARTB)) deallocate(HSTARTB)
      end if

      if(LALTBAS) then
       if(allocated(Cint))    deallocate(Cint)
      end if

      if(LOCBSE2) then
       if(allocated(wWRKw))   deallocate(wWRKw)
       if(allocated(wGBw))    deallocate(wGBw)
       if(allocated(wFLTw))   deallocate(wFLTw)
       if(allocated(wBEenw))  deallocate(wBEenw)
       if(allocated(wvecBEw)) deallocate(wvecBEw)
       if(allocated(wFBEw))   deallocate(wFBEw)
       if(allocated(WB))      deallocate(WB)
      end if

      RETURN
      END
!======================================================================
      SUBROUTINE RXCHFmult_guess_A_elec(NAE,nebf,xxse,GAM_ecore,DAE,CAE)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial regular electron guess density
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer NAE
      double precision xxse(nebf,nebf)
      double precision GAM_ecore(nebf,nebf)
! Variables Returned
      double precision DAE(nebf,nebf)
! Local variables
      double precision CAE(nebf,nebf)
      double precision EVF(nebf)


      call UROOTHAN(CAE,EVF,xxse,GAM_ecore,nebf)

      call RXCHFmult_construct_DE(NAE,nebf,CAE,DAE)


      RETURN
      END
!======================================================================
      subroutine RXCHFmult_guess_elec(LALTBAS,nae,nbe,nebf,nebfBE,
     x                                xxse,xxseBE,
     x                                GAM_ecore,GAM_ecoreBE,
     x                                DAE,DBE,CAE,CBE)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial regular and special electronic guess density
!======================================================================
      implicit none
! Input Variables
      logical LALTBAS
      integer nebf,nebfBE
      integer nae
      integer nbe
      double precision xxse(nebf,nebf)
      double precision xxseBE(nebfBE,nebfBE)
      double precision GAM_ecore(nebf,nebf)
      double precision GAM_ecoreBE(nebfBE,nebfBE)
! Variables Returned
      double precision DAE(nebf,nebf)
      double precision DBE(nebfBE,nebfBE)
      double precision CAE(nebf,nebf)
      double precision CBE(nebfBE,nebfBE)
! Local variables
      integer i,j,k,l
      integer dimint
      integer nocca,noccb
      integer maxind,currind
      logical debug
      double precision ovlap,maxovlap
      double precision C(nebf,nebf)
      double precision C2(nebfBE,nebfBE)
      double precision Ctemp(nebfBE,nebf)
      double precision EVF(nebf)
      double precision EVF2(nebfBE)

      double precision,allocatable :: Cvirt(:,:),tempvec(:)
      double precision,allocatable :: Cint(:,:),Cint_tr(:,:)

      double precision, parameter   :: zero=0.0d+00
      double precision, parameter   :: tol=1.0d-12

      debug=.false.

      if (nae.gt.1) then
       nocca=nae/2
      else
       nocca=nae
      end if
      if (nbe.gt.1) then
       noccb=nbe/2
      else
       noccb=nbe
      end if

      if(allocated(Cvirt)) deallocate(Cvirt)
      allocate(Cvirt(nebf,nebf-nocca))

      DAE=zero
      DBE=zero
      CAE=zero
      CBE=zero

      call UROOTHAN(C,EVF,xxse,GAM_ecore,nebf)

! Store first nocca evectors as occ reg elec vectors
      do i=1,nocca
        do j=1,nebf
          CAE(j,i)=C(j,i)
        end do
      end do

      if(.not.LALTBAS) then

! Store next noccb evectors as occ spec elec vectors
       do i=1,noccb
         do j=1,nebf
           CBE(j,i)=C(j,i+nocca)
         end do
       end do

! Store nebf-nocca-noccb remaining evectors as virt reg and spec elec vectors
       do i=nocca+noccb+1,nebf
         do j=1,nebf
           CAE(j,i-noccb)=C(j,i)
           CBE(j,i-nocca)=C(j,i)
         end do
       end do

      else

! Store nebf-nocca remaining evectors as virt elec vectors
       do i=nocca+1,nebf
         do j=1,nebf
           Cvirt(j,i-nocca)=C(j,i)
         end do
       end do

! Find intersection of restricted basis set and virt elec vectors
       call RXCHFmult_intersection(nebf,nebf-nocca,Cvirt,nebfBE,xxseBE,
     x                             dimint,Ctemp)
       if (dimint.le.noccb) then
        write(*,*) "dim of intersection <= # occ special orbitals"
        write(*,*) "dimint:",dimint
        write(*,*) "noccb:",noccb
        return
       end if

       write(*,*)
       write(*,*) "----------------------------------"
       write(*,*) " Dimension of intersection space:"
       write(*,'(2X,A,1X,I3)') "Max possible (nebfBE)      =",nebfBE
       write(*,'(2X,A,1X,I3)') "Actual (after computation) =",dimint
       write(*,*) "----------------------------------"
       write(*,*)

       if(allocated(Cint)) deallocate(Cint)
       allocate(Cint(nebfBE,dimint))
       do i=1,dimint
         do j=1,nebfBE
           Cint(j,i)=Ctemp(j,i)
         end do
       end do

! Debug: All CBE vectors should be orthogonal to occ CAE vectors
       do i=1,nocca
       do j=1,dimint
         ovlap=zero
         do k=1,nebf
         do l=1,nebfBE
           ovlap=ovlap+CAE(k,i)*Cint(l,j)*xxse(k,l)
         end do
         end do
         if (abs(ovlap).gt.tol) then
          write(*,*)
          write(*,*) "******* ERROR *******"
          write(*,*) "Calculated intersection basis is not orthogonal"
          write(*,*) "to occupied regular vectors"
          write(*,*) "reg occ index, int index, ovlap:",i,j,ovlap
          write(*,*)
         end if
       end do
       end do

! Check orthonormality of new basis
       if(debug) then
        write(*,*) "Overlaps amongst Cint bfs (should be onormal):"
        do i=1,dimint
        do j=1,dimint
          ovlap=zero
          do k=1,nebfBE
          do l=1,nebfBE
            ovlap=ovlap+Cint(k,i)*Cint(l,j)*xxseBE(k,l)
          end do
          end do
          write(*,*) "i,j,ovlap:",i,j,ovlap
        end do
        end do
       end if

! Fill in regular virtuals with new intersection basis
       do i=1,min(nebf-nocca,dimint)
         do j=1,nebfBE
           CAE(j,i+nocca)=Cint(j,i)
         end do
       end do

! Find how "similar" the original Cvirt are as compared to the new basis
! and fill special electronic vectors with those with greatest overlap
       currind=0
       do i=1,nebf-nocca
         maxovlap=zero
         maxind=0
         do j=1,dimint
           ovlap=zero
           do k=1,nebf
           do l=1,nebfBE
             ovlap=ovlap+Cvirt(k,i)*Cint(l,j)*xxse(k,l)
           end do
           end do
           if(abs(ovlap).gt.maxovlap) then
            maxovlap=abs(ovlap)
            maxind=j
           end if
         end do
         if ((i.le.nebfBE).and.(maxind.ne.0)) then
          currind=currind+1
          do j=1,nebfBE
            CBE(j,currind)=Cint(j,maxind)
          end do
         end if
       end do

C! Special occ guess is obtained by diagonalizing truncated overlap
C! matrix and projecting onto intersection space
C       call UROOTHAN(C2,EVF2,xxseBE,GAM_ecoreBE,nebfBE)
C
C       if(allocated(Cint_tr)) deallocate(Cint_tr)
C       allocate(Cint_tr(dimint,nebfBE))
C       do i=1,dimint
C         do j=1,nebfBE
C           Cint_tr(i,j)=Cint(j,i)
C         end do
C       end do
C
C       if(allocated(tempvec)) deallocate(tempvec)
C       allocate(tempvec(dimint))
C
C       do i=1,nebfBE
C         call RXCHF_matmult(dimint,nebfBE,nebfBE,1,
C     x                      Cint_tr,C2(:,i),tempvec)
C         do j=1,nebfBE
C           do k=1,dimint
C             CBE(j,i)=CBE(j,i)+tempvec(k)*Cint(j,k)
C           end do
C         end do
C       end do
C
C       if(allocated(tempvec)) deallocate(tempvec)

       if(allocated(Cint)) deallocate(Cint)

      end if

      call RXCHFmult_construct_DE(NAE,nebf,CAE,DAE)
      call RXCHFmult_construct_DE(NBE,nebfBE,CBE,DBE)

C      if(allocated(Cint_tr)) deallocate(Cint_tr)
      if(allocated(Cvirt)) deallocate(Cvirt)

      return
      end
!======================================================================
      subroutine RXCHFmult_OCBSE(nebf,nae,nbe,vecAE0,vecBE0,FAE,FBE,
     x                       Selec,vecAE,vecBE,AEen,BEen)
! 
!     Perform OCBSE procedure
!       vecAE0: Regular electron coefficients from previous iteration
!       vecBE0: Special electron coefficients from previous iteration
!       vecAE:  Regular electron coefficients from current iteration
!       vecBE:  Special electron coefficients from current iteration
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nae,nbe
      double precision vecAE0(nebf,nebf),vecBE0(nebf,nebf)
      double precision FAE(nebf,nebf),FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
! Variables Returned
      double precision vecAE(nebf,nebf),vecBE(nebf,nebf)
      double precision AEen(nebf),BEen(nebf)
! Local variables
      integer nocca,noccb,nvirt 
      integer noccvirta,noccvirtb
      double precision zero
      parameter(zero=0.0d+00)

      if (nae.gt.1) then
       nocca=nae/2
      else
       nocca=nae
      end if
      if (nbe.gt.1) then
       noccb=nbe/2
      else
       noccb=nbe
      end if
      nvirt=nebf-nocca-noccb
      noccvirta=nebf-noccb
      noccvirtb=nebf-nocca

      vecAE=zero
      vecBE=zero
      AEen=zero
      BEen=zero

      call RXCHF_OCBSE_driver(nebf,nae,nbe,nocca,noccb,
     x                        nvirt,noccvirta,noccvirtb,
     x                        vecAE0,vecBE0,FAE,FBE,Selec,
     x                        vecAE,vecBE,AEen,BEen)

      return
      end
!======================================================================
      subroutine RXCHFmult_OCBSE2(nebf,npebf,nae,nbe,vecAE,vecBE0,FBE,
     x                        Selec,elcam,elcbfc,elcex,
     x                        ampeb2c,agebfcc,kpestr,kpeend,
     x                        vecBE,BEen)
! 
!     Perform OCBSE procedure for special electrons only
!       vecBE0: Special electron coefficients from previous iteration
!       vecAE:  Regular electron coefficients from current iteration
!       vecBE:  Special electron coefficients from current iteration
!
!======================================================================
      implicit none
! Input Variables
      integer nebf,npebf
      integer nae,nbe
      double precision vecAE(nebf,nebf),vecBE0(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
      integer ampeb2c(npebf)               ! Map prim index to contr index
      integer kpestr(nebf)                 ! Map contr index to prim start
      integer kpeend(nebf)                 ! Map contr index to prim end
      integer elcam(npebf,3)               ! Angular mom for electrons
      double precision agebfcc(npebf)      ! Map prim index to contr coeff
      double precision elcex(npebf)        ! Exponents: elec basis
      double precision elcbfc(npebf,3)     ! Basis centers: elec basis
! Variables Returned
      double precision vecBE(nebf,nebf)
      double precision BEen(nebf)
! Local variables
      integer nocca,noccb,nvirt 
      integer noccvirta,noccvirtb
      double precision zero
      parameter(zero=0.0d+00)

      if (nae.gt.1) then
       nocca=nae/2
      else
       nocca=nae
      end if
      if (nbe.gt.1) then
       noccb=nbe/2
      else
       noccb=nbe
      end if
      nvirt=nebf-nocca-noccb
      noccvirta=nebf-noccb
      noccvirtb=nebf-nocca

      vecBE=zero
      BEen=zero

      call RXCHFmult_OCBSE_driver2(nebf,npebf,nae,nbe,nocca,noccb,
     x                         nvirt,noccvirta,noccvirtb,
     x                         vecAE,vecBE0,FBE,Selec,
     x                         elcam,elcbfc,elcex,
     x                         ampeb2c,agebfcc,kpestr,kpeend,
     x                         vecBE,BEen)

      return
      end
!======================================================================
      subroutine RXCHFmult_OCBSE_driver2(nebf,npebf,nae,nbe,nocca,noccb,
     x                               nvirt,noccvirta,noccvirtb,
     x                               vecAE,vecBE0,FBE,Selec,
     x                               elcam,elcbfc,elcex,
     x                               ampeb2c,agebfcc,kpestr,kpeend,
     x                               vecBE,BEen)
!
!     OCBSE Procedure:
!       - Construct transformation matrix for special electrons
!       - Transform FBE and diagonalize to obtain vecBE
!
!  Switching to virtual space (no remnant of occupied in proj basis)
!
!======================================================================
      implicit none
! Input Variables
      integer nebf,npebf
      integer nae,nbe
      integer nocca,noccb,nvirt
      integer noccvirta,noccvirtb
      double precision vecAE(nebf,nebf),vecBE0(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision Selec(nebf,nebf)
      integer ampeb2c(npebf)               ! Map prim index to contr index
      integer kpestr(nebf)                 ! Map contr index to prim start
      integer kpeend(nebf)                 ! Map contr index to prim end
      integer elcam(npebf,3)               ! Angular mom for electrons
      double precision agebfcc(npebf)      ! Map prim index to contr coeff
      double precision elcex(npebf)        ! Exponents: elec basis
      double precision elcbfc(npebf,3)     ! Basis centers: elec basis
! Variables Returned
      double precision vecBE(nebf,nebf)
      double precision BEen(nebf)
! Local variables
      integer i,j
      integer ie1,je1
      integer iec1,jec1
      integer ie1_start,je1_start
      integer ie1_end,je1_end
      integer mo1,mo2,motodiscard
      integer i1,j1,k1
      integer l1,m1,n1
      double precision a1,b1
      double precision cof_ie1,cof_je1
      double precision ans,ovlap,ovlapcheck,maxovlap
      double precision Amat1(3), Bmat1(3)
      double precision ovlaparr(nebf),projorb(nebf) ! fix2
      double precision WB(nebf,noccvirtb)
      double precision WBtrans(noccvirtb,nebf)
      double precision wFBEw(noccvirtb,noccvirtb)
      double precision wSBw(noccvirtb,noccvirtb)
      double precision AUXB(nebf,noccvirtb)
      double precision xAEen(noccvirta)
      double precision xBEen(noccvirtb)
      double precision xvecBE(noccvirtb,noccvirtb)
      double precision blockvecBE(nebf,noccvirtb)
      double precision zero
      parameter(zero=0.0d+00)

      logical debug
      debug=.false.

! Initialize
      WB=zero
      WBtrans=zero
      wFBEw=zero
      wSBw=zero
      AUXB=zero
      xBEen=zero
      xvecBE=zero
      blockvecBE=zero

      if (debug) then
      write(*,*) "nae,nbe:",nae,nbe
      write(*,*) "nocca,noccb:",nocca,noccb
      write(*,*) "MATRIX vecAE:"
      call PREVNU(vecAE,BEen,nebf,nebf,nebf)
      write(*,*) "MATRIX Previous vecBE:"
      call PREVNU(vecBE0,BEen,nebf,nebf,nebf)
      write(*,*) "MATRIX Previous S:"
      call PREVNU(Selec,BEen,nebf,nebf,nebf)
      end if

!!!!!!!!!!!!!!!! Solve for special electronic solution !!!!!!!!!!!!!!!!

! Form special electronic transformation matrix
      do i=1,noccvirtb
        do j=1,nebf
          WB(j,i)=vecAE(j,nocca+i)
        end do
      end do

      if (debug) then
      write(*,*) "MATRIX WB:"
      call PREVNU(WB,xBEen,noccvirtb,nebf,nebf)
      end if

! Form transpose
      do i=1,noccvirtb
        do j=1,nebf
          WBtrans(i,j)=WB(j,i)
        end do
      end do

! Transform FBE as Wtrans * FBE * W
      call RXCHF_matmult(nebf,nebf,nebf,noccvirtb,
     x                   FBE,WB,AUXB)
      call RXCHF_matmult(noccvirtb,nebf,nebf,noccvirtb,
     x                   WBtrans,AUXB,wFBEw)

! Transform AO overlap matrix as Wtrans * S * W
      call RXCHF_matmult(nebf,nebf,nebf,noccvirtb,
     x                   Selec,WB,AUXB)
      call RXCHF_matmult(noccvirtb,nebf,nebf,noccvirtb,
     x                   WBtrans,AUXB,wSBw)

      if (debug) then
      write(*,*) "MATRIX wSBw:"
      call PREVNU(wSBw,xBEen,noccvirtb,noccvirtb,noccvirtb)
      end if

! Diagonalize transformed FBE to obtain solutions in new basis
      call UROOTHAN(xvecBE,xBEen,wSBw,wFBEw,noccvirtb)

      if (debug) then
      WRITE(*,*) "MATRIX xvecBE:"
      call PREVNU(xvecBE,xBEen,noccvirtb,noccvirtb,noccvirtb)
      end if

! Transform evectors to AO basis as vecBE = W * xvecBE
      call RXCHF_matmult(nebf,noccvirtb,noccvirtb,noccvirtb,
     x                   WB,xvecBE,blockvecBE)

! Pass evectors to output variables (zeros for unfilled noccb part)
      do i=1,noccvirtb
        BEen(i)=xBEen(i)
        do j=1,nebf
          vecBE(j,i)=blockvecBE(j,i)
        end do
      end do

      if (debug) then
      WRITE(*,*) "MATRIX New vecBE:"
      call PREVNU(vecBE,BEen,nebf,nebf,nebf)
      end if

      return
      end


      subroutine invertEXPmat(nbf,A,Ainv)
! Works only for exponential matrix exp(A) where A has zero diagonal
      implicit none
! Input variables
      integer nbf
      double precision A(nbf,nbf)
! Output variables
      double precision Ainv(nbf,nbf)
! Local variables
      integer i,j

      do i=1,nbf
        do j=1,nbf
          if (i.ne.j) then
           Ainv(j,i)=-A(j,i)
          else
           Ainv(j,i)=A(j,i)
          end if
        end do
      end do

      return
      end


      subroutine fock2mobasis(nbf,F,C,Fmo)
      implicit none
! Input variables
      integer nbf
      double precision F(nbf,nbf)
      double precision C(nbf,nbf)
! Output variables
      double precision Fmo(nbf,nbf)
! Local variables
      integer i,j
      double precision Ctrans(nbf,nbf)
      double precision aux(nbf,nbf)

      do i=1,nbf
        do j=1,nbf
          Ctrans(j,i)=C(i,j)
        end do
      end do

      call RXCHF_matmult(nbf,nbf,nbf,nbf,F,C,aux)
      call RXCHF_matmult(nbf,nbf,nbf,nbf,Ctrans,aux,Fmo)

      return
      end

!======================================================================
      subroutine RXCHFmult_OCBSE_phase(nebf,nwbf,
     x                                 Selec,vecAE0,vecAE)
!
!     OCBSE Procedure:
!       - Phase vectors of current it against those of previous it
!
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nwbf
      double precision vecAE0(nebf,nwbf)
      double precision Selec(nebf,nebf)
! Input/Outuput Variables
      double precision vecAE(nebf,nwbf)
! Local variables
      integer i,j
      integer maxind
      double precision vecAEp(nebf,nwbf)
      double precision val
      double precision maxovlap
      double precision zero1(nebf),zero2(nwbf)
      double precision zero
      parameter(zero=0.0d+00)

      logical debug
      debug=.false.

! Initialize
      vecAEp=zero
      zero2=zero

      if (debug) then
       write(*,*) "PHASE vecAE0 to phase against:"
       call PREVNU(vecAE0,zero2,nwbf,nebf,nebf)
      end if

      if (debug) then
       write(*,*) "PHASE vecAE before phasing:"
       call PREVNU(vecAE,zero2,nwbf,nebf,nebf)
      end if

! Iterate over each of the [nwbf] vectors
      do i=1,nwbf

! For vec_i of the old set, find vec_j of the new set that matches most
        maxind=0
        maxovlap=zero
        do j=1,nwbf
          call RXCHFmult_vecovlap(nebf,vecAE0(:,i),vecAE(:,j),Selec,val)
          if(abs(val).gt.abs(maxovlap)) then
           maxind=j
           maxovlap=val
          end if
        end do

! Phase found vector if needed
        if(maxovlap.ge.0) then
         vecAEp(:,i)=vecAE(:,maxind)
        else
         vecAEp(:,i)=-vecAE(:,maxind)
        end if

! Destroy vector to ensure it is not found again for another vec_i
        vecAE(:,maxind)=zero

      end do

      vecAE(:,:)=vecAEp(:,:)

      if (debug) then
       write(*,*) "PHASE vecAE after phasing:"
       call PREVNU(vecAE,zero2,nwbf,nebf,nebf)
      end if

      return
      end

!======================================================================
      subroutine RXCHFmult_vecovlap(n,x,y,Selec,ans)
!
!     Calculate overlap of two MOs <x | y>
!
!======================================================================
      implicit none
! Input Variables
      integer n
      double precision x(n)
      double precision y(n)
      double precision Selec(n,n)
! Output Variables
      double precision ans
! Local Variables
      integer i,j

      ans=0.0d+00

      do i=1,n
      do j=1,n
        ans=ans+x(i)*y(j)*Selec(i,j)
      end do
      end do

      return
      end

!======================================================================
      subroutine RXCHFmult_OCBSE_transF(debug,nebf,nwbf,
     x                                  vecAE,FBE,WB,wFBEw)
!
!     OCBSE Procedure:
!       - Construct transformation matrix for special electrons (WB)
!       - Transform FBE (wFBEw)
!
!  Transformation to virtual space of regular electrons (in vecAE)
!  vecAE on input are the virtuals only (no occupied) to facilitate
!  portability to restricted basis set formulation
!
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nwbf
      double precision vecAE(nebf,nwbf)
      double precision FBE(nebf,nebf)
      logical debug
! Variables Returned
      double precision WB(nebf,nwbf)
      double precision wFBEw(nwbf,nwbf)
! Local variables
      integer i,j
      double precision AUXB(nebf,nwbf)
      double precision zeroarr(nebf)
      double precision zero, one
      parameter(zero=0.0d+00, one=1.0d+00)

! Initialize
      WB=zero
      wFBEw=zero
      AUXB=zero
      zeroarr=zero

      if (debug) then
       write(*,*) "nwbf:",nwbf
       write(*,*) "MATRIX vecAE:"
       call PREVNU(vecAE,zeroarr,nwbf,nebf,nebf)
      end if

! Form special electronic transformation matrix
      do i=1,nwbf
        do j=1,nebf
          WB(j,i)=vecAE(j,i)
        end do
      end do

      if (debug) then
       write(*,*) "MATRIX WB:"
       call PREVNU(WB,zeroarr,nwbf,nebf,nebf)
      end if

! Transform FBE as Wtrans * FBE * W
      call dgemm('n','n',nebf,nwbf,nebf,one,FBE,nebf,WB,nebf,
     x           zero,AUXB,nebf)
      call dgemm('t','n',nwbf,nwbf,nebf,one,WB,nebf,AUXB,nebf,
     x           zero,wFBEw,nwbf)

      if (debug) then
       write(*,*) "MATRIX wFBEw:"
       call PREVNU(wFBEw,zeroarr,nwbf,nwbf,nwbf)
      end if

      return
      end

!======================================================================
      subroutine RXCHFmult_OCBSE_diag(debug,nebf,nwbf,WB,wFBEw,
     x                                wvecBEw,wBEenw,vecBE,BEen)
!
!     OCBSE Procedure:
!       - Read in transformation (WB) 
!       - Read in transformed Fock matrix (wFBEw)
!       - Diagonalize Fock matrix (wvecBEw,wBEenw)
!       - Transform eigenvectors (vecBE,BEen)
!
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nwbf
      double precision WB(nebf,nwbf)
      double precision wFBEw(nwbf,nwbf)
      logical debug
! Variables Returned
      double precision BEen(nebf)
      double precision wBEenw(nwbf)
      double precision vecBE(nebf,nebf)
      double precision wvecBEw(nwbf,nwbf)
! Local variables
      integer i,j
      integer istat
      double precision workq(1)
      double precision, allocatable :: work(:)
      double precision zero
      parameter(zero=0.0d+00)

! Initialize
      wBEenw=zero
      wvecBEw=zero

! Diagonalize transformed FBE to obtain solutions in new basis

! Query work array size for diagonalization and allocate work array
      istat=0
      wvecBEw(:,:)=wFBEw(:,:)
      call dsyev("v","l",nwbf,wvecBEw,nwbf,wBEenw,
     x           workq,-1,istat)
      if (istat.ne.0) then
       write(*,*) "Error in dsyev query"
      end if
      if(allocated(work)) deallocate(work)
      allocate(work(int(workq(1))))

! Diagonalize
      istat=0
      wvecBEw(:,:)=wFBEw(:,:)
      call dsyev("v","l",nwbf,wvecBEw,nwbf,wBEenw,
     x           work,int(workq(1)),istat)
      if (istat.ne.0) then
       write(*,*) "Error in dsyev"
      end if
      if(allocated(work)) deallocate(work)

      if (debug) then
       WRITE(*,*) "MATRIX wvecBEw:"
       call PREVNU(wvecBEw,wBEenw,nwbf,nwbf,nwbf)
      end if

      call RXCHFmult_OCBSE_transV(debug,nebf,nwbf,WB,wvecBEw,wBEenw,
     x                            vecBE,BEen)

      return
      end

!======================================================================
      subroutine RXCHFmult_OCBSE_transV(debug,nebf,nwbf,WB,wvecBEw,
     x                                  wBEenw,vecBE,BEen)
!
!     OCBSE Procedure:
!       - Transform eigenvectors (wvecBEw,wBEenw) -> (vecBE,BEen)
!       - Fill in extra parts of transformed solutions with zeros
!
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nwbf
      double precision WB(nebf,nwbf)
      double precision wBEenw(nwbf)
      double precision wvecBEw(nwbf,nwbf)
      logical debug
! Variables Returned
      double precision vecBE(nebf,nebf)
      double precision BEen(nebf)
! Local variables
      integer i,j
      double precision blockvecBE(nebf,nwbf)
      double precision zero, one
      parameter(zero=0.0d+00, one=1.0d+00)

! Initialize
      blockvecBE=zero
      vecBE=zero
      BEen=zero

! Transform evectors to AO basis as vecBE = W * xvecBE
      call dgemm('n','n',nebf,nwbf,nwbf,one,WB,nebf,wvecBEw,nwbf,
     x           zero,blockvecBE,nebf)

! Pass evectors to output variables (zeros for unfilled part)
      do i=1,nwbf
        BEen(i)=wBEenw(i)
        do j=1,nebf
          vecBE(j,i)=blockvecBE(j,i)
        end do
      end do

      if (debug) then
       WRITE(*,*) "MATRIX New vecBE:"
       call PREVNU(vecBE,BEen,nebf,nebf,nebf)
      end if

      return
      end

!======================================================================
      subroutine RXCHFmult_OCBSE_transVt(debug,nebf,nwbf,WB,
     x                                   Selec,vecBE,wvecBEw)
!
!     OCBSE Procedure:
!       - Transform vectors from AO basis to W basis vecBE -> wvecBEw
!
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nwbf
      double precision WB(nebf,nwbf)
      double precision vecBE(nebf,nebf)
      double precision Selec(nebf,nebf)
      logical debug
! Variables Returned
      double precision wvecBEw(nwbf,nwbf)
! Local variables
      integer i,j
      double precision testvec(nebf,nebf)
      double precision WBinv(nwbf,nebf)
      double precision blockvecBE(nwbf,nebf)
      double precision wSBw(nwbf,nwbf)
      double precision zeroarr(nebf)
      double precision zero, one
      parameter(zero=0.0d+00, one=1.0d+00)

! Initialize
      WBinv=zero
      blockvecBE=zero
      wSBw=zero
      wvecBEw=zero
      zeroarr=zero

      if (debug) then
       WRITE(*,*) "MATRIX Pretransformed vecBE:"
       call PREVNU(vecBE,zeroarr,nebf,nebf,nebf)
      end if

      if (debug) then
       WRITE(*,*) "MATRIX Pretransformed WB:"
       call PREVNU(WB,zeroarr,nwbf,nebf,nebf)
      end if

! W^(-1) = W^t * S
      call dgemm('t','n',nwbf,nebf,nebf,one,WB,nebf,Selec,nebf,
     x           zero,WBinv,nwbf)

! Test: W * W^(-1) = S(W basis) = I
      if(debug) then
       call dgemm('n','n',nwbf,nwbf,nebf,one,WBinv,nwbf,WB,nebf,
     x            zero,wSBw,nwbf)
       write(*,*) "MATRIX wSBw:"
       call PREVNU(wSBw,zeroarr,nwbf,nwbf,nwbf)
      end if

! Transform evectors to W basis as xvecBE = (W^t * S_AO) * vecBE
      call dgemm('n','n',nwbf,nebf,nebf,one,WBinv,nwbf,vecBE,nebf,
     x           zero,blockvecBE,nwbf)

! Pass evectors to output variables (unpassed part should be zero)
      do i=1,nwbf
        do j=1,nwbf
          wvecBEw(j,i)=blockvecBE(j,i)
        end do
      end do

      if (debug) then
       WRITE(*,*) "MATRIX Transformed wvecBEw:"
       call PREVNU(wvecBEw,zeroarr,nwbf,nwbf,nwbf)
      end if

! Lowdin orthogonalize
      call RXCHF_loworth(nwbf,nwbf,wvecBEw)

      if (debug) then
       WRITE(*,*) "MATRIX Lowdin-orth wvecBEw:"
       call PREVNU(wvecBEw,zeroarr,nwbf,nwbf,nwbf)
      end if

      if (debug) then
       call dgemm('n','n',nebf,nebf,nwbf,one,WB,nebf,blockvecBE,nwbf,
     x            zero,testvec,nebf)
       WRITE(*,*) "MATRIX Pretransformed TEST vecBE:"
       call PREVNU(testvec,zeroarr,nebf,nebf,nebf)
      end if

      return
      end

      subroutine checkovlap(m,n,c,S)
      implicit none

      integer m,n
      double precision c(n,m),S(n,n)

      integer i,j
      double precision ct(m,n),aux(m,n),cSc(m,m)
      double precision zero1(m)

      zero1=0.0d+00
      ct=0.0d+00
      aux=0.0d+00

      do i=1,n
      do j=1,m
        ct(j,i)=c(i,j)
      end do
      end do

      call RXCHF_matmult(m,n,n,n,ct,S,aux)
      call RXCHF_matmult(m,n,n,m,aux,c,cSc)

      write(*,*) "MATRIX cSc:"
      call PREVNU(cSc,zero1,m,m,m)

      return
      end

!======================================================================
      subroutine RXCHF_sochgbas(npr,nocc,nwbf,nebf,W0,W,Sao,G,H,D)

! Overwrites gradient (G), Hessian (H) and displacement vector (D)
! which have been calculated in the W0 basis with their corresponding
! values in the W basis provided W0 and W share a common AO basis
!
!    npr : Number of rotation parameters
!   nocc : Number of occupied MOs
!   nwbf : Number of W basis functions
!   nebf : Number of AO basis functions
!     W0 : W basis of previous iteration
!      W : W basis of current iteration
!    Sao : Overlap matrix in AO basis
!      G : Gradient of previous iteration to be transformed
!          stored as (i-1)*nocc + j + 1
!          for i=occ (outer index) and j=virt (inner index)
!      H : Hessian of previous iteration to be transformed
!          stored as 2-dim version of G
!      D : Displacement of previous iteration to be transformed
!          stored as G
!======================================================================
      implicit none

! Input variables
      integer          :: npr
      integer          :: nocc
      integer          :: nwbf
      integer          :: nebf
      double precision :: W0(nebf,nwbf)
      double precision :: W(nebf,nwbf)
      double precision :: Sao(nebf,nebf)

! Input/Output variables
      double precision :: G(npr)
      double precision :: H(npr,npr)
      double precision :: D(npr)

! Local variables
      integer          :: i,j,k,l
      integer          :: ia,ib,ic,id
      integer          :: ind1,ind2,ind3,ind4
      integer          :: nvirt
      double precision :: A(nwbf,nwbf)
      double precision :: At(nwbf,nwbf)
      double precision :: Wt(nwbf,nebf)
      double precision :: aux(nebf,nwbf)
      double precision :: aux2(nwbf,nwbf)
      double precision :: G0(npr)
      double precision :: H0(npr,npr)
      double precision :: D0(npr)
      double precision :: zero1(npr)
      double precision :: zero2(nwbf)
C ARS(
      double precision :: W0t(nwbf,nebf)
      double precision :: B(nwbf,nwbf)
      double precision :: Bt(nwbf,nwbf)
      double precision :: Gtest(npr)
      double precision :: Htest(npr,npr)
      double precision :: Dtest(npr)
C )

      logical          :: debug

      double precision, parameter :: zero=0.0d+00
      double precision, parameter :: one=1.0d+00

      debug=.false.

      A=zero
      Wt=zero
      aux=zero
      aux2=zero
      zero1=zero
      zero2=zero
      nvirt=nwbf-nocc
      G0(:)=G(:)
      H0(:,:)=H(:,:)
      D0(:)=D(:)
      G=zero
      H=zero
      D=zero
C ARS(
      W0t=zero
      B=zero
      Bt=zero
      Gtest=zero
      Htest=zero
      Dtest=zero
C )

      if((nocc*nvirt).ne.npr) then
       write(*,*) "nvirt != npr in RXCHF_sochbas"
       write(*,*) "Exiting..."
       call abrt
      end if

      if(debug) then
       write(*,*) "Initial gradient:"
       write(*,*) G0
       write(*,*)
       write(*,*) "Initial displacement:"
       write(*,*) D0
       write(*,*)
       write(*,*) "Initial Hessian:"
       call PREVNU(H0,zero1,npr,npr,npr)
      end if

! Form transformation matrix A from old W basis to new W basis
      do i=1,nebf
      do j=1,nwbf
        Wt(j,i)=W(i,j)
      end do
      end do
      call RXCHF_matmult(nebf,nebf,nebf,nwbf,Sao,W0,aux)
      call RXCHF_matmult(nwbf,nebf,nebf,nwbf,Wt,aux,A)

      if(debug) then
       write(*,*) "W0 -> W transformation:"
       call PREVNU(A,zero2,nwbf,nwbf,nwbf)
      end if

! Lowdin orthogonalize
      call RXCHF_loworth(nwbf,nwbf,A)

      if (debug) then
       WRITE(*,*) "Lowdin-orth W0 -> W transformation:"
       call PREVNU(A,zero2,nwbf,nwbf,nwbf)
      end if

! Check overlaps
      if(debug) then

       do i=1,nebf
       do j=1,nwbf
         W0t(j,i)=W0(i,j)
       end do
       end do

       aux=zero
       call RXCHF_matmult(nebf,nebf,nebf,nwbf,Sao,W0,aux)
       call RXCHF_matmult(nwbf,nebf,nebf,nwbf,W0t,aux,aux2)
       write(*,*) "W0^t S W0:"
       call PREVNU(aux2,zero2,nwbf,nwbf,nwbf)

       aux=zero
       aux2=zero
       call RXCHF_matmult(nebf,nebf,nebf,nwbf,Sao,W,aux)
       call RXCHF_matmult(nwbf,nebf,nebf,nwbf,Wt,aux,aux2)
       write(*,*) "W^t S W:"
       call PREVNU(aux2,zero2,nwbf,nwbf,nwbf)

       do i=1,nwbf
       do j=1,nwbf
         At(j,i)=A(i,j)
       end do
       end do

       aux2=zero
       call RXCHF_matmult(nwbf,nwbf,nwbf,nwbf,At,A,aux2)
       write(*,*) "A^t A:"
       call PREVNU(aux2,zero2,nwbf,nwbf,nwbf)

       aux2=zero
       call RXCHF_matmult(nwbf,nwbf,nwbf,nwbf,A,At,aux2)
       write(*,*) "A A^t:"
       call PREVNU(aux2,zero2,nwbf,nwbf,nwbf)

      end if

! Transform quantities
      ind1=0
      do i=1,nocc
      do ia=nocc+1,nwbf
        ind1=ind1+1                    ! ind1: (i,a)

        ind2=0
        do j=1,nocc
        do ib=nocc+1,nwbf
          ind2=ind2+1                  ! ind2: (j,b)

          G(ind1)=G(ind1)+G0(ind2)*
     x            (A(i,j)*A(ia,ib)+A(i,ib)*A(ia,j))
!     x            (A(j,i)*A(ib,ia)+A(ib,i)*A(j,ia))
          D(ind1)=D(ind1)+D0(ind2)*
     x            (A(i,j)*A(ia,ib)+A(i,ib)*A(ia,j))
!     x            (A(j,i)*A(ib,ia)+A(ib,i)*A(j,ia))

          ind3=0
          do k=1,nocc
          do ic=nocc+1,nwbf
            ind3=ind3+1                ! ind3: (k,c)

            ind4=0
            do l=1,nocc
            do id=nocc+1,nwbf
              ind4=ind4+1              ! ind4: (l,d)

              H(ind1,ind2)=H(ind1,ind2)+H0(ind3,ind4)*
     x                     (A(i,k)*A(ia,ic)*A(j,l)*A(ib,id)
     x                     +A(i,ic)*A(ia,k)*A(j,id)*A(ib,l)
     x                     +A(i,k)*A(ia,ic)*A(j,id)*A(ib,l)
     x                     +A(i,ic)*A(ia,k)*A(j,l)*A(ib,id))
!     x                     (A(k,i)*A(ic,ia)*A(l,j)*A(id,ib)
!     x                     +A(ic,i)*A(k,ia)*A(id,j)*A(l,ib)
!     x                     +A(k,i)*A(ic,ia)*A(id,j)*A(l,ib)
!     x                     +A(ic,i)*A(k,ia)*A(l,j)*A(id,ib))

            end do
            end do

          end do
          end do

        end do
        end do

      end do
      end do

      if(debug) then
       write(*,*) "Final gradient:"
       write(*,*) G
       write(*,*)
       write(*,*) "Final displacement:"
       write(*,*) D
       write(*,*)
       write(*,*) "Final Hessian:"
       call PREVNU(H,zero1,npr,npr,npr)
      end if

C ARS( check back transformation
      if(debug) then
       aux=zero
       call RXCHF_matmult(nebf,nebf,nebf,nwbf,Sao,W,aux)
       call RXCHF_matmult(nwbf,nebf,nebf,nwbf,W0t,aux,B)
       write(*,*) "W -> W0 transformation:"
       call PREVNU(B,zero2,nwbf,nwbf,nwbf)

       call RXCHF_loworth(nwbf,nwbf,B)

       if (debug) then
        WRITE(*,*) "Lowdin-orth W -> W0 transformation:"
        call PREVNU(B,zero2,nwbf,nwbf,nwbf)
       end if

       do i=1,nwbf
       do j=1,nwbf
         Bt(j,i)=B(i,j)
       end do
       end do

       aux2=zero
       call RXCHF_matmult(nwbf,nwbf,nwbf,nwbf,Bt,B,aux2)
       write(*,*) "B^t B:"
       call PREVNU(aux2,zero2,nwbf,nwbf,nwbf)

       aux2=zero
       call RXCHF_matmult(nwbf,nwbf,nwbf,nwbf,B,Bt,aux2)
       write(*,*) "B B^t:"
       call PREVNU(aux2,zero2,nwbf,nwbf,nwbf)

       ind1=0
       do i=1,nocc
       do ia=nocc+1,nwbf
         ind1=ind1+1                    ! ind1: (i,a)
       
         ind2=0
         do j=1,nocc
         do ib=nocc+1,nwbf
           ind2=ind2+1                  ! ind2: (j,b)
       
           Gtest(ind1)=Gtest(ind1)+G(ind2)*
     x             (A(i,j)*A(ia,ib)+A(i,ib)*A(ia,j))
!     x            (A(j,i)*A(ib,ia)+A(ib,i)*A(j,ia))
           Dtest(ind1)=Dtest(ind1)+D(ind2)*
     x             (A(i,j)*A(ia,ib)+A(i,ib)*A(ia,j))
!     x            (A(j,i)*A(ib,ia)+A(ib,i)*A(j,ia))
       
           ind3=0
           do k=1,nocc
           do ic=nocc+1,nwbf
             ind3=ind3+1                ! ind3: (k,c)
       
             ind4=0
             do l=1,nocc
             do id=nocc+1,nwbf
               ind4=ind4+1              ! ind4: (l,d)
       
               Htest(ind1,ind2)=Htest(ind1,ind2)+H(ind3,ind4)*
     x                      (A(i,k)*A(ia,ic)*A(j,l)*A(ib,id)
     x                      +A(i,ic)*A(ia,k)*A(j,id)*A(ib,l)
     x                      +A(i,k)*A(ia,ic)*A(j,id)*A(ib,l)
     x                      +A(i,ic)*A(ia,k)*A(j,l)*A(ib,id))
!     x                     (A(k,i)*A(ic,ia)*A(l,j)*A(id,ib)
!     x                     +A(ic,i)*A(k,ia)*A(id,j)*A(l,ib)
!     x                     +A(k,i)*A(ic,ia)*A(id,j)*A(l,ib)
!     x                     +A(ic,i)*A(k,ia)*A(l,j)*A(id,ib))
       
             end do
             end do
       
           end do
           end do
       
         end do
         end do
       
       end do
       end do
       
       if(debug) then
        write(*,*) "Back-transformed gradient:"
        write(*,*) Gtest
        write(*,*)
        write(*,*) "Back-transformed displacement:"
        write(*,*) Dtest
        write(*,*)
        write(*,*) "Back-transformed Hessian:"
        call PREVNU(Htest,zero1,npr,npr,npr)
       end if
      end if
C )

      return
      end

!======================================================================
      subroutine RXCHF_sonewt(n,it,H0,H,G0,G,D0,D)

! Determines Hessian and displacement vector at current iteration
! from current gradient and previous iteration Hessian, displacement
! and gradient using Fischer & Almlof 1992 JPC Eqs 1-3
! Also resets SOSCF or scales displacement if too large and skips
! Hessian update procedure if close to convergence using protocols
! from GAMESS SOSCF routines
!
!      n : Number of rotation parameters
!     it : Current iteration
!     H0 : Hessian of previous iteration
!      H : Hessian of current iteration
!     G0 : Gradient of previous iteration
!      G : Gradient of current iteration
!     D0 : Displacement of previous iteration
!      D : Displacement of current iteration
!======================================================================
      implicit none

! Input variables
      integer          :: n
      integer          :: it
      double precision :: H0(n,n)
      double precision :: G0(n),G(n)
      double precision :: D0(n)

! Output variables
      double precision :: H(n,n)
      double precision :: D(n)

! Local variables
      integer          :: i,j
      double precision :: maxG
      double precision :: alp
      double precision :: del(n),v(n)
      double precision :: E
      double precision :: c1,c2
      double precision :: aux(n,n)
      double precision :: sqcdf,scal
      double precision :: ddot
      logical          :: debug

      double precision, parameter :: zero=0.0d+00
      double precision, parameter :: one=1.0d+00
      double precision, parameter :: toobig=1.0d+00
      double precision, parameter :: bigrot=0.1d+00
      double precision, parameter :: small=1.0d-08

      debug=.false.

      H=zero
      D=zero
      aux=zero
      del=zero
      v=zero

      maxG=zero
      do i=1,n
        if(abs(G(i)).gt.maxG) maxG=abs(G(i))
      end do

      if (it.eq.1) then

       H(:,:)=H0(:,:)
       D(:)=D0(:)

      else if (maxG.lt.small) then

       H(:,:)=H0(:,:)
C       call RXCHF_matmult(n,n,n,1,H0,G,D)
       call dgemv('n',n,n,one,H0,n,G,1,zero,D,1)
       do i=1,n
         D(i)=-D(i)
       end do

       if(debug) then
        write(*,*) "small displacement:",D
       end if

      else

! Apply Eqs 1-3 from Fischer & Almlof
       do i=1,n
         del(i)=G(i)-G0(i)
       end do

       alp=ddot(n,D0,1,del,1)
       alp=one/alp

       if(debug) then
        write(*,*) "alp:",alp
       end if

C       if (abs(alp).gt.(1.0d+00/small)) then
C
C        H(:,:)=H0(:,:)
C        call dgemv('n',n,n,one,H0,n,G,1,zero,D,1)
C        do i=1,n
C          D(i)=-D(i)
C        end do
C        write(*,*) "small displacement:",D
C
C       else

C       call RXCHF_matmult(n,n,n,1,H0,del,v)
        call dgemv('n',n,n,one,H0,n,del,1,zero,v,1)

        if(debug) then
         write(*,*) "v:",v
        end if

        c1=ddot(n,del,1,v,1)
        c1=alp*c1
        c1=c1+one
        c1=alp*c1
        c2=-alp
        if(debug) then
         write(*,*) "c1,c2:",c1,c2
        end if
        do i=1,n
        do j=1,n
          E=c1*D0(j)*D0(i)
          E=E+c2*(D0(j)*v(i)+v(j)*D0(i))
          H(j,i)=H0(j,i)+E
        end do
        end do

C       call RXCHF_matmult(n,n,n,1,H,G,D)
        call dgemv('n',n,n,one,H,n,G,1,zero,D,1)
        do i=1,n
          D(i)=-D(i)
        end do

        if(debug) then
         write(*,*) "large displacement:",D
        end if

C       end if

      end if

! Scale displacement in case value is too large (poached from GAMESS)
      SQCDF = dSQRT(DDOT(n,D,1,D,1)/dble(n))
      IF(SQCDF.GT.TOOBIG  .AND.  it.GT.5) THEN
         WRITE(*,9011) SQCDF
         it = 0
      END IF
      IF(SQCDF.GT.BIGROT) THEN
         IF(it.GT.0) WRITE(*,9020) SQCDF
         SCAL=dSQRT(BIGROT/SQCDF)
C ARS( turn off
         write(*,*) "scal:",scal
         write(*,*) "not doing anything"
         SCAL=one
C )
         CALL DSCAL(n,SCAL,D,1)
      END IF

 9011 FORMAT(1X,'*** RESETTTING SOSCF, UPON ENCOUNTERING A HUGE',
     *          ' TOTAL ROTATION=',1P,E12.3)
 9020 FORMAT(1X,'SOSCF IS SCALING ROTATION ANGLE MATRIX, SQCDF=',F12.6)

      return
      end


      subroutine moovlap(nbf,coeffs1,coeffs2,Sao,ovlap)
      implicit none

      integer, intent(in)  :: nbf
      real*8,  intent(in)  :: coeffs1(nbf),coeffs2(nbf)
      real*8,  intent(in)  :: Sao(nbf,nbf)

      real*8,  intent(out) :: ovlap

      integer              :: i,j
      real*8, parameter    :: zero=0.0d+00

      ovlap=zero

      do i=1,nbf
      do j=1,nbf
        ovlap=ovlap+coeffs1(i)*coeffs2(j)*Sao(i,j)
      end do
      end do

      return
      end


      subroutine printmat(mat,dim1,dim2)
      implicit none

      integer, intent(in)    :: dim1,dim2
      real*8,  intent(in)    :: mat(dim1,dim2)

      integer :: i,j

      do i=1,dim1
      write(*,9000) (mat(i,j),j=1,dim2)
      end do

 9000 format(20(2X,G15.6))

      return
      end subroutine printmat

