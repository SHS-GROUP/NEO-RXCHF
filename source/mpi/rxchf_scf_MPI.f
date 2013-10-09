!======================================================================
      subroutine RXCHF_scf_MPI(nproc,rank,
     x                         nelec,NAE,NBE,NPRA,NPRB,NEBFLT,NUCST,
     x                         npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                         ngtg1,ng1,ng2,ng3,ng4,
     x                         NG2CHK,NG3CHK,NG4CHK,
     x                         read_CE,read_CP,
     x                         LG4DSCF,LG3DSCF,LG2DSCF,
     x                         LSOSCF,LOCBSE,LCMF,LADDEXCH,
     x                         ng2prm,ng3prm,nat,pmass,cat,zan,
     x                         bcoef1,gamma1,
     x                         KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                         ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
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
      include 'mpif.h'

! Input Variables
      integer nproc,rank
      logical LNEOHF
      logical LOCBSE   ! Use OCBSE scheme as is (restricted variational freedom for reg/sp elecs)
      logical LOCBSE2  ! Use modified OCBSE scheme (complete variational freedom for reg elecs)
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
      integer NEBFLT
      integer NUCST
      integer nebf
      integer npbf
      integer nebf2
      integer npbf2
      integer ngee
      integer ng1             ! Number of GAM1 integrals
      integer ng2             ! Number of GAM2 integrals
      integer ng3             ! Number of GAM3 integrals
      integer ng4             ! Number of GAM3 integrals
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
      double precision xxsp(npbf,npbf)  ! Nuc overlap matrix
      double precision GAM_pcore(npbf2)
      double precision GAM_ecore(nebf2)
      double precision GAM_ep(ng1)
      double precision GAM_ee(ngee)

!     integer noccE ! Number of occupied elec orbs
!     integer noccP ! Number of occupied nuc orbs
      integer maxit,maxmicroit

      integer i,ielec
      integer j
      integer k
      integer l

      integer NPBFLT

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
      double precision DBE(NEBF,NEBF)
      double precision DBE0(NEBF,NEBF)
      double precision VECAE(NEBF,NEBF)
      double precision VECAE0(NEBF,NEBF)
      double precision VECBE(NEBF,NEBF)
      double precision VECBE0(NEBF,NEBF)
      double precision AEE(NEBF)
      double precision BEE(NEBF)
      double precision FAE(nebf,nebf)
      double precision XFAE(nebf,nebf)
      double precision FBE(nebf,nebf)
      double precision XFBE(nebf,nebf)

      double precision DP(NPBF,NPBF)
      double precision DP0(NPBF,NPBF)
      double precision VECP(NPBF,NPBF)
      double precision EP(NPBF)
      double precision FP(npbf,npbf)
      double precision XFP(npbf,npbf)

      double precision FAEint(nebf,nebf)
      double precision FBEint(nebf,nebf)
      double precision FPint(npbf,npbf)

      double precision SBE_XCHF(nebf,nebf)
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
      integer L0w,L1w
      integer NFT15
      integer NFT16
      double precision FLT(NEBFLT) !FLT: Lower triangle focke
      double precision HSTARTA(NPRA)
      double precision HSTARTB(NPRB)
      double precision GRADA(NPRA)
      double precision GRADB(NPRB)
      double precision PGRADA(NPRA)
      double precision PGRADB(NPRB)
      double precision DISPLIA(NPRA)
      double precision DISPLIB(NPRB)
      double precision DGRADA(NPRA)  ! WRK1
      double precision DGRADB(NPRB)  ! WRK1
      double precision DISPLA(NPRA)  ! WRK2
      double precision DISPLB(NPRB)  ! WRK2
      double precision UPDTA(NPRA)   ! WRK3
      double precision UPDTB(NPRB)   ! WRK3
      double precision DISPLNA(NPRA) ! WRK1+NPR
      double precision DISPLNB(NPRB) ! WRK1+NPR
      double precision DGRADIA(NPRA) ! WRK2+NPR
      double precision DGRADIB(NPRB) ! WRK2+NPR
      double precision UPDTIA(NPRA)  ! WRK3+NPR
      double precision UPDTIB(NPRB)  ! WRK3+NPR
      double precision ORBGRDA
      double precision ORBGRDB
      double precision SMALL
      double precision SOGTOL ! ORBGRAD TOL to activate soscf
      double precision XA(NPRA)
      double precision XB(NPRB)
      double precision GA(nebf,nebf) !G(L0,L0)
      double precision GB(nebf,nebf) !G(L0,L0)
      double precision WRK(nebf) !WRK(L0)
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
C ARS( testing variables
      logical LNOINT
      double precision FBEmo(nebf,nebf)
      double precision work1(nebf,nebf)
      double precision work2(nebf,nebf)
C )
      integer*4 ierr

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
       if(rank.eq.0) write(*,*) "Using LOCBSE2"

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
      end if

      if((LOCBSE).and.(rank.eq.0)) write(*,*) "Using LOCBSE"
      LGAM4=.true. ! Always calculate five-particle integrals
C ARS( no interaction
      LNOINT=.false.
C      LNOINT=.true.
      if ((LNOINT).and.(rank.eq.0)) then
       write(*,*)
       write(*,*) "******************"
       write(*,*) "  NO INTERACTION  "
       write(*,*) "******************"
       write(*,*)
      end if
C )

!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------(
!      call class_nuc_rep(nat,zan,cat,E_nuc)
      if (rank.eq.0) then
       open(800,file='ENUCRP.dat',status='unknown')
       read(800,*) E_nuc
       close(800)
       write(*,*)'READ IN NUC REPULSION'
      end if
      call MPI_BCAST(E_nuc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!----------CALCULATE-CLASSICAL-NUCLEAR-REPULSION-ENERGY----------------)

!--------------READ-INTEGRALS-NEEDED-FOR-NEO-HF------------------------(
      nebf2=nebf*nebf
      npbf2=npbf*npbf
      if (rank.eq.0) then
       call read_nuc_ovlap(npbf,xxsp)
       write(*,*)
       write(*,*)'READ IN NUC OVLAP'
       call read_elec_ovlap(nebf,xxse)
       write(*,*)'READ IN ELEC OVLAP'
       call read_GAM_ecore(nebf,nebf2,GAM_ecore)
       write(*,*)'READ IN GAM_ECORE'
       call read_GAM_pcore(npbf,npbf2,GAM_pcore)
       write(*,*)'READ IN GAM_PCORE'
       call read_GAM_ep(nebf,npbf,ng1,GAM_ep)
       write(*,*)'READ IN GAM_EP'
       call read_GAM_ee(nebf,ngee,GAM_ee)
       write(*,*)'READ IN GAM_EE'
       write(*,*)
      end if
      call MPI_BCAST(xxsp,npbf2,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xxse,nebf2,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(GAM_ecore,nebf2,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(GAM_pcore,npbf2,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(GAM_ep,ng1,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(GAM_ee,ngee,MPI_DOUBLE_PRECISION,
     x               0,MPI_COMM_WORLD,ierr)

! Testing
      write(*,*) "xxsp:"
      call prt_lower_triangle(npbf,(npbf+npbf*npbf)/2,xxsp)
      write(*,*) "xxse:"
      call prt_lower_triangle(nebf,(nebf+nebf*nebf)/2,xxse)

!--------------READ-INTEGRALS-NEEDED-FOR-NEO-HF------------------------)

!-------------INITIAL-GUESSES------------------------------------------(
      if(read_CE) then

       if (rank.eq.0) then
!       READ IN GUESS FOR E:
!       call read_elec_density(nebf,nelec,DE)
        call RXCHFmult_read_CAE(nebf,NAE,DAE,VECAE0)
        call RXCHFmult_read_CBE(nebf,NBE,DBE,VECBE0)
       end if

       call MPI_BCAST(DAE,nebf2,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(VECAE0,nebf2,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DBE,nebf2,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(VECBE0,nebf2,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)

! Testing
       write(*,*) "DAE:"
       call prt_lower_triangle(nebf,(nebf+nebf*nebf)/2,DAE)
       write(*,*) "VECAE0:"
       call prt_lower_triangle(nebf,(nebf+nebf*nebf)/2,VECAE0)
       write(*,*) "DBE:"
       call prt_lower_triangle(nebf,(nebf+nebf*nebf)/2,DBE)
       write(*,*) "VECBE0:"
       call prt_lower_triangle(nebf,(nebf+nebf*nebf)/2,VECBE0)

      else
!       STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
        if(rank.eq.0) write(*,*)'ABOUT TO CALL guess_A_elec'
!       call guess_elec(nelec,nebf,xxse,GAM_ecore,DE)
        if ((LOCBSE).or.(LOCBSE2)) then
          call RXCHFmult_guess_elec(nae,nbe,nebf,xxse,GAM_ecore,
     x                          DAE,DBE,VECAE0,VECBE0)
          if(rank.eq.0) write(*,*)'BACK FROM guess_elec for OCBSE'
        else
         call RXCHFmult_guess_A_elec(NAE,nebf,xxse,GAM_ecore,DAE,VECAE0)
         call RXCHFmult_guess_A_elec(NBE,nebf,xxse,GAM_ecore,DBE,VECBE0)
          if(rank.eq.0) write(*,*)'BACK FROM guess_elec'
        end if
      end if
      if(read_CP) then
!        READ IN GUESS FOR N:

       if (rank.eq.0) then
        call read_nuc_density(npbf,1,NUCST,DP)
       end if

       call MPI_BCAST(DP,npbf2,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)

! Testing
       write(*,*) "DP:"
       call prt_lower_triangle(npbf,(npbf+npbf*npbf)/2,DP)

      else
!        STANDARD GUESS:  HCORE FOR NUC AND ELEC DENSITIES:
         if(rank.eq.0) write(*,*)'ABOUT TO CALL guess_prot'
!        call guess_prot(NUCST,npbf,nebf,xxsp,GAM_pcore,GAM_ep,DE,DP)
         call guess_prot2(NUCST,npbf,xxsp,GAM_pcore,DP)
         if(rank.eq.0) write(*,*)'BACK FROM guess_prot'
      end if

C ARS( debug: print out initial guess MOs here
      if ((LCMF).and.(rank.eq.0)) then
       write(*,*)
       write(*,*) "------------------"
       write(*,*) "INITIAL GUESS MOs:"
       write(*,*) "------------------"
       write(*,*)
       WRITE(*,9610)
       call PREVNU(vecAE0,AEE,nebf,nebf,nebf)
       WRITE(*,9620)
       call PREVNU(vecBE0,BEE,nebf,nebf,nebf)
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
         LSOSCFA=.true.
         LSOSCFB=.true.
         if((nae.eq.1).or.LOCBSE) LSOSCFA=.FALSE.
C         if((nbe.eq.1).or.(LOCBSE).or.(LOCBSE2)) LSOSCFB=.FALSE.
         if((nbe.eq.1).or.(LOCBSE)) LSOSCFB=.FALSE.
      else
         LSOSCFA=.false.
         LSOSCFB=.false.
      end if
      if(LSOSCFA) THEN
         NFT15=15
         if (rank.eq.0) then
          OPEN(NFT15, FILE='WORK15', STATUS='UNKNOWN',
     *         ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         end if
         NA=nae/2
      end if
      if(LSOSCFB) THEN
         NFT16=16
         if (rank.eq.0) then
          OPEN(NFT16, FILE='WORK16', STATUS='UNKNOWN',
     *         ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
         end if
         NB=nbe/2
      end if
!-------------SETUP-FOR-POSSIBLE-SOSCF---------------------------------)

!     SET CONVERGENCE CRITERIA AND MAXIMUM ITERATIONS 
      TOLE = 1.0D-06
      TOLP = 1.0D-04
      maxit=100
      maxmicroit=200
      if(LOCBSE) maxit=400

!     ZERO OUT 'OLD' DENSITY MATRICES
      DAE0=0.0d+00
      DBE0=0.0d+00
      DP0=0.0d+00

!     BEGIN XCSCF ITERATIONS
      if(rank.eq.0) WRITE(*,9500)
      E_total_old=0.0d+00

      DO I=1,MAXIT

C Call HF Fock build for NAE regular electrons
         call RXCHFmult_fock_hf(LCMF,nebf,nebf2,NAE,ngee,
     x                          DAE,GAM_ecore,GAM_ee,
     x                          FAE,E_HF,E_HF_ecore,E_HF_ee)

C Call XCHF Fock build for NBE special electrons and one QM particle
         call RXCHFmult_fock_xchf(LGAM4,LG4DSCF,LG4IC,
     x                   LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                   NG4CHK,NG3CHK,NG2CHK,
     x                   dimXCHF4,dimXCHF3,dimXCHF2,
     x                   npebf,nebf,nebf2,npbf,npbf2,NBE,
     x                   ngee,ng1,ng2,ng3,ng4,DBE,DP,
     x                   XCHF_GAM4,XCHF_GAM3,XCHF_GAM2,XCHF_GAM2s,
     x                   ng2prm,ngtg1,ng3prm,
     x                   nat,pmass,cat,zan,bcoef1,gamma1,
     x                   KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                   ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                   FBE,FP,SBE_XCHF,SP_XCHF,
     x                   E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                   E_XCHF_gam3,E_XCHF_gam4,
     x                   S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
         call RXCHF_fock_int_MPI(nproc,rank,
     x                           LCMF,LADDEXCH,nelec,NAE,NBE,
     x                           nebf,nebf2,npbf,npbf2,
     x                           ng1,ng2,ng3,ng4,
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
          call add2fock(nebf,FBEint,FBE)
      end if
C )

          IF ((LCMF).and.(rank.eq.0)) then
           npbflt=npbf*(npbf+1)/2
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

!--------------FORM-FOCK-MATRICES-AND-CALC-ENERGY-COMPONENTS-----------)
         if((I.eq.1).and.(rank.eq.0)) then
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
         PGRADB=0.0d+00

         if (rank.eq.0) then
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
         end if

         do ielec=1,maxmicroit

       if((.not.((locbse).or.(locbse2))).or.(ielec.eq.1)) then
C Call HF Fock build for NAE regular electrons
           call RXCHFmult_fock_hf(LCMF,nebf,nebf2,NAE,ngee,
     x                            DAE,GAM_ecore,GAM_ee,
     x                            FAE,E_HF,E_HF_ecore,E_HF_ee)

C Call XCHF Fock build for NBE special electrons and one QM particle
           call RXCHFmult_fock_xchf(LGAM4,LG4DSCF,LG4IC,
     x                     LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                     NG4CHK,NG3CHK,NG2CHK,
     x                     dimXCHF4,dimXCHF3,dimXCHF2,
     x                     npebf,nebf,nebf2,npbf,npbf2,NBE,
     x                     ngee,ng1,ng2,ng3,ng4,DBE,DP,
     x                     XCHF_GAM4,XCHF_GAM3,XCHF_GAM2,XCHF_GAM2s,
     x                     ng2prm,ngtg1,ng3prm,
     x                     nat,pmass,cat,zan,bcoef1,gamma1,
     x                     KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                     ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                     FBE,FP,SBE_XCHF,SP_XCHF,
     x                     E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                     E_XCHF_gam3,E_XCHF_gam4,
     x                     S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
           call RXCHFmult_fock_int(LCMF,LADDEXCH,nelec,NAE,NBE,
     x                             nebf,nebf2,npbf,npbf2,
     x                             ng1,ng2,ng3,ng4,
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
            call add2fock(nebf,FBEint,FBE)
      end if
C )
C ARS( microiteration
          E_total=E_HF+E_XCHF+E_int+E_nuc
C )
      end if
C )

         if (LOCBSE) then
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
     x                     LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                     NG4CHK,NG3CHK,NG2CHK,
     x                     dimXCHF4,dimXCHF3,dimXCHF2,
     x                     npebf,nebf,nebf2,npbf,npbf2,NBE,
     x                     ngee,ng1,ng2,ng3,ng4,DBE,DP,
     x                     XCHF_GAM4,XCHF_GAM3,XCHF_GAM2,XCHF_GAM2s,
     x                     ng2prm,ngtg1,ng3prm,
     x                     nat,pmass,cat,zan,bcoef1,gamma1,
     x                     KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                     ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                     FBE,FP,SBE_XCHF,SP_XCHF,
     x                     E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                     E_XCHF_gam3,E_XCHF_gam4,
     x                     S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
           call RXCHFmult_fock_int(LCMF,LADDEXCH,nelec,NAE,NBE,
     x                             nebf,nebf2,npbf,npbf2,
     x                             ng1,ng2,ng3,ng4,
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
            call add2fock(nebf,FBEint,FBE)
      end if
C )

            IF ((LCMF).and.(rank.eq.0)) then
             npbflt=npbf*(npbf+1)/2
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

C ARS( OCBSE/SOSCF

! Transform FBE (calculated at end of previous iteration) to new W basis
!  - W updated with new vecA from this iteration
!  - vecBE in old W basis still from previous iteration
         call RXCHFmult_OCBSE_transF(nebf,nocca,nwbf,
     x                               vecAE,FBE,WB,wFBEw)

!-----------------------POSSIBLE-SOSCF-BETA----------------------------(
         if(LSOSCFB) THEN
           ITER=IELEC
           EIGAVL = ITER.GT.1
         end if
         IF(LSOSCFB .AND.  EIGAVL) THEN                ! first it. skip SOSCF (diag to get EE)
!!!!!!      --> SETUP LOWER TRIANGLE FOCKE FOR SOSCF
           call pack_LT(nwbf,nwbfLT,wFBEw,wFLTw)
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
              call RXCHFmult_OCBSE_transV(nebf,nwbf,WB,wvecBEw,wBEenw, ! eigenvalues useless
     x                                    vecBE,BEe)
              call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)
              GO TO 450  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-BETA----------------------------)

! No SOSCF
!  - Diagonalize Fock matrix in W basis of this iteration
!  - Obtain updated vecBE in W basis of this iteration
         call RXCHFmult_OCBSE_diag(nebf,nwbf,WB,wFBEw,
     x                             wvecBEw,wBEenw,vecBE,BEe)
         call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)

  450 CONTINUE

         CALL DENDIF(DBE0,DBE,NEBF,DIFFBE)
         CALL COPYDEN(DBE0,DBE,NEBF)

C )

! Calculate energy for this it and Fock matrices for next it

C Call HF Fock build for NAE regular electrons
           call RXCHFmult_fock_hf(LCMF,nebf,nebf2,NAE,ngee,
     x                            DAE,GAM_ecore,GAM_ee,
     x                            FAE,E_HF,E_HF_ecore,E_HF_ee)

C Call XCHF Fock build for NBE special electrons and one QM particle
           call RXCHFmult_fock_xchf(LGAM4,LG4DSCF,LG4IC,
     x                     LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                     NG4CHK,NG3CHK,NG2CHK,
     x                     dimXCHF4,dimXCHF3,dimXCHF2,
     x                     npebf,nebf,nebf2,npbf,npbf2,NBE,
     x                     ngee,ng1,ng2,ng3,ng4,DBE,DP,
     x                     XCHF_GAM4,XCHF_GAM3,XCHF_GAM2,XCHF_GAM2s,
     x                     ng2prm,ngtg1,ng3prm,
     x                     nat,pmass,cat,zan,bcoef1,gamma1,
     x                     KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                     ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                     FBE,FP,SBE_XCHF,SP_XCHF,
     x                     E_XCHF,E_XCHF_gam1,E_XCHF_gam2,
     x                     E_XCHF_gam3,E_XCHF_gam4,
     x                     S_total,S_gam1,S_gam2)

C Call interaction Fock build for all particles
           call RXCHFmult_fock_int(LCMF,LADDEXCH,nelec,NAE,NBE,
     x                             nebf,nebf2,npbf,npbf2,
     x                             ng1,ng2,ng3,ng4,
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
            call add2fock(nebf,FBEint,FBE)
      end if
C )

            IF ((LCMF).and.(rank.eq.0)) then
             npbflt=npbf*(npbf+1)/2
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
           call pack_LT(nebf,nebfLT,FBE,FLT)
          call SOGRAD(GRADB,FLT,vecBE,WRK,NPRB,NB,L0,L1,NEBFLT,ORBGRDB)
!!!!!!      IF(ORBGRD.LT.SMALL) THEN
!!!!!!         DIFF = ZERO
!!!!!!         CVGING=.TRUE.
!!!!!!         GO TO 700  ! Check on convergence behavior
!!!!!!      END IF
            IF(ORBGRDB.LT.SOGTOL  .OR.  ITSOB.GT.0) THEN
              IF(ITSOB.EQ.0) THEN   ! only on first SOSCF it. set up approx Hess
             WRITE(*,9800)
                 call SOHESS(HSTARTB,BEE,NPRB,L0,NB,NB)
              END IF
              ITSOB = ITSOB+1
           call SONEWT(HSTARTB,GRADB,PGRADB,DISPLIB,DGRADB,DISPLB,UPDTB,
     *                 DISPLNB,DGRADIB,UPDTIB,ORBGRDB,NPRB,ITSOB,NFT16)
            call SOTRAN(DISPLIB,vecBE,GB,WRK,NPRB,L0,L1,NB,NB,ORBGRDB)
             CALL DCOPY(NPRB,GRADB,1,PGRADB,1)
              call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)
              GO TO 850  ! Use the new C's to form new density (change)
            END IF
         END IF
!-----------------------POSSIBLE-SOSCF-BETA----------------------------)

  800 CONTINUE
!        call ROOTHAN(DBE,vecBE,BEE,xxse,FBE,nebf,nelec,1,NUCST)
         call UROOTHAN(vecBE,BEE,xxse,FBE,nebf)
         call RXCHFmult_construct_DE(NBE,nebf,vecBE,DBE)

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
         if (rank.eq.0) then
          if((LSOSCFA).and.(LSOSCFB)) then
             WRITE(*,9151) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP,ORBGRDA,ORBGRDB
          else if ((LSOSCFA).and.(.not.(LSOSCFB))) then
             WRITE(*,9150) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP,ORBGRDA
          else if ((LSOSCFB).and.(.not.(LSOSCFA))) then
             WRITE(*,9150) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP,ORBGRDB
          else
             WRITE(*,9100) IELEC,E_total,Delta_E_tot,
     x                     DIFFAE,DIFFBE,DIFFP
          end if
C ARS( debug: print out MOs here
          if (LCMF) then
           WRITE(*,9610)
           call PREVNU(vecAE,AEE,nebf,nebf,nebf)
           WRITE(*,9620)
           call PREVNU(vecBE,BEE,nebf,nebf,nebf)
           WRITE(*,9700)
           call PREVNU(vecp,EP,npbf,npbf,npbf)
          end if
C )
! Output the vectors for this iteration for restart if necessary:
          call write_MOs(860,nebf,VECAE)
          call write_MOs(861,nebf,VECBE)
          call write_MOs(853,npbf,VECP)
         end if

         LDIFFE=( (DIFFAE.LT.TOLE).and.(DIFFBE.LT.TOLE) )
         IF(LDIFFE) GOTO 200
         IF(IELEC.EQ.MAXIT) GOTO 10

       END DO  ! microiterations

  200 CONTINUE
!      IF WE GET HERE - MICROITERATION CONVERGENCE ACHIEVED
       if(rank.eq.0) write(*,2000) i,ielec

       IF(DIFFP.LT.TOLP) GOTO 100
       IF(I.EQ.MAXIT) GOTO 10

      END DO   ! iterations
 
  10  CONTINUE
!     IF WE GET HERE SOMETHING WENT WRONG

      if (rank.eq.0) then
       if(LSOSCFA) close(NFT15)
       if(LSOSCFB) close(NFT16)

       WRITE(*,*)
       WRITE(*,*)'WARNING:  ITERATION LIMIT EXCEEDED'
       WRITE(*,*)
      end if

      E_total=zero

  100 CONTINUE
!     IF WE GET HERE WE ARE DONE - CONVERGENCE ACHIEVED

      if (rank.eq.0) then

       if(LSOSCFA) close(NFT15)
       if(LSOSCFB) close(NFT16)

!     PRINT FINAL ENERGY AND PUNCH THE ORBITALS
       WRITE(*,9200) E_total,I

       WRITE(*,9300) E_nuc,E_HF_ecore,E_HF_ee,E_HF,
     x   E_XCHF_gam1,E_XCHF_gam2,E_XCHF_gam3,E_XCHF_gam4,E_XCHF,
     x   E_int_OMG2,E_int_OMG3,E_int_OMG4,E_int,
     x   S_total,E_total

!  OUTPUT ELEC AND NUC EIGENVALUES AND EIGENVECTORS
       WRITE(*,9610)
       call PREVNU(vecAE,AEE,nebf,nebf,nebf)
       WRITE(*,9620)
       call PREVNU(vecBE,BEE,nebf,nebf,nebf)
       WRITE(*,9700)
       call PREVNU(vecp,EP,npbf,npbf,npbf)

! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------(
!     subroutine write_MOs(IFIL,nbf,VEC)
!     IFIL=852 :: FinalCE.dat
!     IFIL=853 :: FinalCP.dat
!     IFIL=860 :: FinalCAE.dat
!     IFIL=861 :: FinalCBE.dat
       call write_MOs(860,nebf,VECAE)
       call write_MOs(861,nebf,VECBE)
       call write_MOs(853,npbf,VECP)
! PUNCH-OUT-THE-FINAL-VECTORS-FOR-E-AND-NUC----------------------------)

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

      return
      end

