C=======================================================================
      subroutine RXCHF_driver_MPI(nproc,rank,
     x                            nblocks,blockrank,
     x                            nelec,nae,nbe,nae_alp,nae_beta,nucst,
     x                            nebf,npebf,npbf,nat,ngtg1,ngee,
     x                            nebfBE,elindBE,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            ng2chk,ng3chk,ng4chk,
     x                            read_CE,read_CP,
     x                            read_GAM2,read_GAM3,read_GAM4,
     x                            LG2IC,LG3IC,LG4IC,
     x                            LG2DSCF,LG3DSCF,LG4DSCF,
     x                            LRXCUHF,
     x                            LSOSCF,LDIIS,LSORXCHF,OCBSETYP,
     x                            LALTBAS,LREORDGS,
     x                            LDBG,EXCHLEV)

C Driver to calculate RXCHF integrals for nbe > 1
C   XCHF_GAM* : integrals needed for XCHF contribution
C    INT_GAM* : integrals needed for interaction contribution
C  INT_GAM*ex : integrals needed for exchange contribution
C
C Maximum dimension of integral needed is nbe + 1
C In this routine, both XCHF_GAM? and INT_GAM? are calculated up to
C this dimension and XCHF_GAM[max] is deallocated at the end
C  - this doesn't incur any additional cost (aside from temp memory)
C  - the case of XCHF_GAM4 (nbe > 3) is dealt with separately at the end
C=======================================================================
      implicit none
      include 'mpif.h'

C Input variables
      integer nproc,rank               ! MPI variables
      integer nblocks
      integer blockrank
      integer ng3block
      integer blockstart,blockend
      integer nelec                    ! Total number of electrons
      integer nae                      ! Number of regular electrons
      integer nae_alp                  ! Number of regular alpha electrons
      integer nae_beta                 ! Number of regular beta electrons
      integer nbe                      ! Number of special electrons
      integer nucst                    ! Nuclear state for which to form density
      integer nebf                     ! Number of contracted electronic basis functions
      integer npebf                    ! Number of primitive electronic basis functions
      integer npbf                     ! Number of proton basis functions
      integer nebfBE                   ! Size of elec basis set to use for NBE elecs
      integer npebfBE                  ! Analog for NBE basis set
      integer nat                      ! Number of classical nuclei
      integer ngtg1                    ! Number of geminal functions
      integer ng1,ng2,ng3,ng4,ngee     ! Numbers of integrals (contracted)
      integer ng2_int,ng2_xchf
      integer ng3_int,ng3_xchf
      integer ng4_int,ng4_xchf
      integer ng1prm,ng2prm,ng3prm     ! Numbers of integrals (primitive)
      integer ng2chk,ng3chk,ng4chk     ! Number of chunks to split integral calculations
      integer ELCAM(npebf,3)           ! Angular mom for electrons
      integer NUCAM(npbf,3)            ! Angular mom for quantum nuclei
      integer AMPEB2C(npebf)           ! Map primitive index to contracted
      integer KPESTR(nebf)             ! Map contracted index to primitive start
      integer KPEEND(nebf)             ! Map contracted index to primitive end
      integer KPESTR_be(nebfBE)        ! Analog for NBE basis set
      integer KPEEND_be(nebfBE)        ! Analog for NBE basis set
      integer elindBE(nebfBE)          ! Contracted indices of NBE basis set
      logical LG2IC                    ! In-core storage of 3-particle integrals
      logical LG3IC                    ! In-core storage of 4-particle integrals
      logical LG4IC                    ! In-core storage of 5-particle integrals
      logical LG2DSCF                  ! Direct SCF for 3-particle integrals
      logical LG3DSCF                  ! Direct SCF for 4-particle integrals
      logical LG4DSCF                  ! Direct SCF for 5-particle integrals
      logical LRXCUHF                  ! Flag for RXCUHF (separate SCF routine)
      logical LSOSCF                   ! SOSCF where applicable
      logical LDIIS                    ! DIIS where applicable
      logical LSORXCHF                 ! SORXCHF where applicable
      logical LALTBAS                  ! Flag to denote distinct special electron basis
      logical LREORDGS                 ! Flag to reorder electronic guess if LALTBAS
      logical LDBG                     ! Debugging flag
      logical read_CE,read_CP          ! Read in orbitals
      logical read_GAM2                ! 
      logical read_GAM3                ! Read in integrals
      logical read_GAM4                ! 
      logical LADDEXCH
      integer EXCHLEV                  ! 0=ne; 1=ae; 2=fe
      integer OCBSETYP                 ! 1=OCBSE; 2=OCBSE2
      double precision pmass           ! Mass of nonelectron quantum particle 
      double precision zan(nat)        ! Classical nuclear charges
      double precision cat(3,nat)      ! Classical nuclear coordinates (au)
      double precision bcoef1(ngtg1)   ! Geminal b_k
      double precision gamma1(ngtg1)   ! Geminal gamma_k
      double precision ELCEX(npebf)    ! Electronic basis function exponents
      double precision NUCEX(npbf)     ! Nuclear basis function exponents
      double precision ELCBFC(npebf,3) ! Electronic basis function centers
      double precision NUCBFC(npbf,3)  ! Nuclear basis function centers
      double precision AGEBFCC(npebf)  ! Map primitive index to contract coeff
      double precision AGNBFCC(npbf)   ! Nuclear contraction coeff

C Local variables
      character*4  istring              ! File I/O variables
      integer      nebf2,npbf2          !
      integer      nebflt,npbflt        ! Convenient quantities
      integer      nebfBE2,nebfBElt     !
      integer      npra,nprb            ! Number of distinct electron pairs
      integer      dimXCHF2             !
      integer      dimXCHF3             ! Dimensions of XCHF integral arrays
      integer      dimXCHF4             !
      integer      dimINT2              !
      integer      dimINT3              ! Dimensions of interaction integral arrays
      integer      dimINT4              !
      integer      dimINT2ex            ! Dimensions of exchange integral arrays
      integer      dimINT3ex            !
      logical      LXCHFSYMM            ! Flag to symmetrize XCHF3 integrals (nbe > 2)
      integer      unitno
      integer      i,j,istat
      integer      contrind,primind
      integer      numprims,primstart,primend
      integer      currcontrind,currprimind
      logical      laddbasis                     ! Working variable
      integer      XELCAM(npebf,3)               !
      integer      XAMPEB2C(npebf)               ! Dummy basis set variables
      integer      XKPESTR(nebf)                 !
      integer      XKPEEND(nebf)                 !
      double precision              :: XELCEX(npebf)        !
      double precision              :: XELCBFC(npebf,3)     ! Dummy basis set variables
      double precision              :: XAGEBFCC(npebf)      !
      double precision              :: wtime,wtime1,wtime2  ! Timing variables
      double precision, allocatable :: XCHF_GAM2(:)         ! 3-particle XCHF integrals
      double precision, allocatable :: XCHF_GAM2s(:)        ! 3-particle XCHF overlap integrals
      double precision, allocatable :: XCHF_GAM3(:)         ! 4-particle XCHF integrals
      double precision, allocatable :: XCHF_GAM4(:)         ! 4-particle XCHF integrals
      double precision, allocatable :: INT_GAM2(:)          ! 3-particle interaction integrals
      double precision, allocatable :: INT_GAM3(:)          ! 4-particle interaction integrals
      double precision, allocatable :: INT_GAM4(:)          ! 4-particle interaction integrals
      double precision, allocatable :: INT_GAM2ex(:)        ! 3-particle exchange integrals
      double precision, allocatable :: INT_GAM3ex1(:)       ! 4-particle exchange integrals
      double precision, allocatable :: INT_GAM3ex2(:)       ! 4-particle exchange integrals
      integer,          allocatable :: ELCAM_be(:,:)        !
      integer,          allocatable :: AMPEB2C_be(:)        ! 
      double precision, allocatable :: ELCEX_be(:)          ! Analogs for NBE basis set 
      double precision, allocatable :: ELCBFC_be(:,:)       ! 
      double precision, allocatable :: AGEBFCC_be(:)        ! 
      integer(kind=4),  allocatable :: int_pinds(:)         ! Packed prot indices of nontrival INT_GAM4
      integer(kind=8),  allocatable :: int_einds(:)         ! Packed elec indices of nontrival INT_GAM4
      integer(kind=4),  allocatable :: xchf_pinds(:)        ! Packed prot indices of nontrival XCHF_GAM4
      integer(kind=8),  allocatable :: xchf_einds(:)        ! Packed elec indices of nontrival XCHF_GAM4
#if MPI32
      integer*4    ierr
      integer*4    npebf_,nebf_
      integer*4    npebfBE_,nebfBE_
#else
      integer      ierr
      integer      npebf_,nebf_
      integer      npebfBE_,nebfBE_
#endif

C Reordering of basis set and guess if restricted basis set calculation
      if (LALTBAS) then

C Reorder electronic basis set such that special electron subset is first
C Do on master process then pass necessary variables to slave processes
       if(rank.eq.0) then

        currcontrind=1
        currprimind=1
        npebfBE=0

C First add bfs that are also in special electronic set
        do i=1,nebfBE
          contrind=elindBE(i)
          primstart=KPESTR(contrind)
          primend=KPEEND(contrind)
          numprims=(primend-primstart)+1
          XKPESTR(currcontrind)=currprimind
          XKPEEND(currcontrind)=currprimind+numprims-1
          KPESTR_be(currcontrind)=currprimind
          KPEEND_be(currcontrind)=currprimind+numprims-1
          do primind=primstart,primend
            XELCAM(currprimind,1)=ELCAM(primind,1)
            XELCAM(currprimind,2)=ELCAM(primind,2)
            XELCAM(currprimind,3)=ELCAM(primind,3)
            XAMPEB2C(currprimind)=currcontrind
            XELCEX(currprimind)=ELCEX(primind)
            XELCBFC(currprimind,1)=ELCBFC(primind,1)
            XELCBFC(currprimind,2)=ELCBFC(primind,2)
            XELCBFC(currprimind,3)=ELCBFC(primind,3)
            XAGEBFCC(currprimind)=AGEBFCC(primind)
            currprimind=currprimind+1
            npebfBE=npebfBE+1
          end do
          currcontrind=currcontrind+1
        end do

C Add remaining bfs
        do i=1,nebf
          laddbasis=.true.
          do j=1,nebfBE
            contrind=elindBE(j)
            if (i.eq.contrind) laddbasis=.false.
          end do
          if (laddbasis) then
           primstart=KPESTR(i)
           primend=KPEEND(i)
           numprims=(primend-primstart)+1
           XKPESTR(currcontrind)=currprimind
           XKPEEND(currcontrind)=currprimind+numprims-1
           do primind=primstart,primend
             XELCAM(currprimind,1)=ELCAM(primind,1)
             XELCAM(currprimind,2)=ELCAM(primind,2)
             XELCAM(currprimind,3)=ELCAM(primind,3)
             XAMPEB2C(currprimind)=currcontrind
             XELCEX(currprimind)=ELCEX(primind)
             XELCBFC(currprimind,1)=ELCBFC(primind,1)
             XELCBFC(currprimind,2)=ELCBFC(primind,2)
             XELCBFC(currprimind,3)=ELCBFC(primind,3)
             XAGEBFCC(currprimind)=AGEBFCC(primind)
             currprimind=currprimind+1
           end do
           currcontrind=currcontrind+1
          end if
        end do

        if(allocated(AMPEB2C_be)) deallocate(AMPEB2C_be)
        allocate( AMPEB2C_be(npebfBE),stat=istat )
        if(allocated(ELCEX_be)) deallocate(ELCEX_be)
        allocate( ELCEX_be(npebfBE),stat=istat )
        if(allocated(AGEBFCC_be)) deallocate(AGEBFCC_be)
        allocate( AGEBFCC_be(npebfBE),stat=istat )
        if(allocated(ELCAM_be)) deallocate(ELCAM_be)
        allocate( ELCAM_be(npebfBE,3),stat=istat )
        if(allocated(ELCBFC_be)) deallocate(ELCBFC_be)
        allocate( ELCBFC_be(npebfBE,3),stat=istat )

C Transfer basis set variables
        do i=1,nebf
          KPESTR(i)=XKPESTR(i)
          KPEEND(i)=XKPEEND(i)
        end do
        do i=1,npebf
          ELCAM(i,1)=XELCAM(i,1)
          ELCAM(i,2)=XELCAM(i,2)
          ELCAM(i,3)=XELCAM(i,3)
          AMPEB2C(i)=XAMPEB2C(i)
          ELCEX(i)=XELCEX(i)
          ELCBFC(i,1)=XELCBFC(i,1)
          ELCBFC(i,2)=XELCBFC(i,2)
          ELCBFC(i,3)=XELCBFC(i,3)
          AGEBFCC(i)=XAGEBFCC(i)
          if (i.le.npebfBE) then
           ELCAM_be(i,1)=XELCAM(i,1)
           ELCAM_be(i,2)=XELCAM(i,2)
           ELCAM_be(i,3)=XELCAM(i,3)
           AMPEB2C_be(i)=XAMPEB2C(i)
           ELCEX_be(i)=XELCEX(i)
           ELCBFC_be(i,1)=XELCBFC(i,1)
           ELCBFC_be(i,2)=XELCBFC(i,2)
           ELCBFC_be(i,3)=XELCBFC(i,3)
           AGEBFCC_be(i)=XAGEBFCC(i)
          end if
        end do

C Output reordered basis
        write(*,*)
        write(*,*) "Reordered basis:"
        write(*,*)
        write(*,*)'PRIM  CONT    ANG     EXPONENT CONTRACT  -X- -Y- -Z-'
        write(*,*)'INDEX INDEX   MOM                 COEF'
        do i=1,npebf
          write(*,1000) i,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),
     x                  ELCAM(i,3),ELCEX(i),AGEBFCC(i),
     x                  ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
        end do
        write(*,*)
        write(*,*)' CHECK CONTRACTED ELECTRONIC BASIS FUNCTIONS '
        write(*,*)'CONT INDEX    KPESTR     KPEEND'
        do i=1,nebf
           write(*,2000) i,KPESTR(i),KPEEND(i)
        end do

C Output special electron basis
        write(*,*)
        write(*,*) "Special electron basis:"
        write(*,*)
        write(*,*)'PRIM  CONT    ANG     EXPONENT CONTRACT  -X- -Y- -Z-'
        write(*,*)'INDEX INDEX   MOM                 COEF'
        do i=1,npebfBE
          write(*,1000) i,AMPEB2C_be(i),ELCAM_be(i,1),ELCAM_be(i,2),
     x                  ELCAM_be(i,3),ELCEX_be(i),AGEBFCC_be(i),
     x                  ELCBFC_be(i,1),ELCBFC_be(i,2),ELCBFC_be(i,3)
        end do
        write(*,*)
        write(*,*)' CHECK CONTRACTED ELECTRONIC BASIS FUNCTIONS '
        write(*,*)'CONT INDEX    KPESTR     KPEEND'
        do i=1,nebfBE
           write(*,2000) i,KPESTR_be(i),KPEEND_be(i)
        end do

       end if ! rank 0

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Broadcast static variables including array dimensions to slave processes
       call MPI_BCAST(npebfBE,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)

! Allocate arrays on slave processes
       if (rank.gt.0) then
        if(allocated(AMPEB2C_be)) deallocate(AMPEB2C_be)
        allocate( AMPEB2C_be(npebfBE),stat=istat )
        if(allocated(ELCEX_be)) deallocate(ELCEX_be)
        allocate( ELCEX_be(npebfBE),stat=istat )
        if(allocated(AGEBFCC_be)) deallocate(AGEBFCC_be)
        allocate( AGEBFCC_be(npebfBE),stat=istat )
        if(allocated(ELCAM_be)) deallocate(ELCAM_be)
        allocate( ELCAM_be(npebfBE,3),stat=istat )
        if(allocated(ELCBFC_be)) deallocate(ELCBFC_be)
        allocate( ELCBFC_be(npebfBE,3),stat=istat )
       end if

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C Broadcast array data
#if MPI32
       nebf_=int(nebf,kind=4)
       nebfBE_=int(nebfBE,kind=4)
       npebf_=int(npebf,kind=4)
       npebfBE_=int(npebfBE,kind=4)
#else
       nebf_=nebf
       nebfBE_=nebfBE
       npebf_=npebf
       npebfBE_=npebfBE
#endif
       call MPI_BCAST(AMPEB2C_be,npebfBE_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ELCEX_be,npebfBE_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(AGEBFCC_be,npebfBE_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ELCAM_be,3*npebfBE_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ELCBFC_be,3*npebfBE_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(AMPEB2C,npebf_,MPI_INTEGER8,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ELCEX,npebf_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(AGEBFCC,npebf_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ELCAM,3*npebf_,MPI_INTEGER8,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ELCBFC,3*npebf_,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(KPESTR,nebf_,MPI_INTEGER8,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(KPEEND,nebf_,MPI_INTEGER8,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(KPESTR_be,nebfBE_,MPI_INTEGER8,
     x                0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(KPEEND_be,nebfBE_,MPI_INTEGER8,
     x                0,MPI_COMM_WORLD,ierr)

      else ! not restricted basis set - duplicate orig basis set info

       npebfBE=npebf
       KPESTR_be=KPESTR
       KPEEND_be=KPEEND
       if(allocated(AMPEB2C_be)) deallocate(AMPEB2C_be)
       allocate( AMPEB2C_be(npebfBE),stat=istat )
       if(allocated(ELCEX_be)) deallocate(ELCEX_be)
       allocate( ELCEX_be(npebfBE),stat=istat )
       if(allocated(AGEBFCC_be)) deallocate(AGEBFCC_be)
       allocate( AGEBFCC_be(npebfBE),stat=istat )
       if(allocated(ELCAM_be)) deallocate(ELCAM_be)
       allocate( ELCAM_be(npebfBE,3),stat=istat )
       if(allocated(ELCBFC_be)) deallocate(ELCBFC_be)
       ELCBFC_be=ELCBFC
       AMPEB2C_be=AMPEB2C
       ELCEX_be=ELCEX
       AGEBFCC_be=AGEBFCC
       ELCAM_be=ELCAM
       ELCBFC_be=ELCBFC

      end if

      if(rank.eq.0) then

C Calculate inexpensive integrals usually calculated in main driver
C since the all-electron basis set may have been reordered
       call class_nuc_rep(nat,zan,cat)

       call elec_ovlap(npebf,nebf,nebf*nebf,
     x                 AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

       call check_elec_ovlap(nebf)

       call nuc_ovlap(npbf,npbf*npbf,
     x                 AGNBFCC,NUCEX,NUCAM,NUCBFC)

       call check_nuc_ovlap(npbf)

       call calc_GAM_epcore(nebf,npebf,npbf,nebf*nebf,npbf*npbf,
     x                      nat,pmass,zan,cat,
     x                      AMPEB2C,AGEBFCC,AGNBFCC,
     x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

       call calc_GAM_ee(nebf,npebf,ngee,
     x                  AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

! Process electronic guess (must be given) if LREORDGS (defaults to True)
       if ((LALTBAS).and.(LREORDGS)) then

        call process_elec_guess(LDBG,nebf,nebfBE,elindBE,nbe/2) ! assumes closed shell of sp elecs

! Otherwise one must ensure that the guess is already consistent with the
! reordered basis set of this run (i.e. generated from a previous calc)
       else

        write(*,*) "NOT REORDERING ELECTRONIC GUESS"
        write(*,*) " => assumes MO coeffs are wrt reordered basis set"

       end if

      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C Calculate integrals
      wtime = MPI_WTIME()

C Initialize dimensions
      dimXCHF2  = 1
      dimXCHF3  = 1
      dimXCHF4  = 1
      dimINT2   = 1
      dimINT3   = 1
      dimINT4   = 1
      dimINT2ex = 1
      dimINT3ex = 1

C Set logical flag for -ae
      if(EXCHLEV.gt.0) then
       LADDEXCH=.true.
      else
       LADDEXCH=.false.
      end if

C Variable for file I/O
      write(istring,'(I4.4)') rank
      unitno=rank+200

C nbe >= 1
C  - calculate two-particle XCHF integrals and write to disk
C     => XCHF_GAM1
C  - calculate three-particle integrals and store in memory
C     => INT_GAM2
C     => INT_GAM2ex only needed if RXCHF-ae
C     => XCHF_GAM2 only needed if nbe >= 2

! Calculate two-particle integrals and write to disk with master process
      ng1=nebfBE*nebfBE*npbf*npbf
      ng1prm=npebfBE*npebfBE*npbf*npbf

      if (rank.eq.0) then
       write(*,*)
       write(*,*) "---------------------------"
       write(*,*) " Calculating:   XCHF_GAM1  "
       write(*,*) "---------------------------"
       write(*,*)

       call RXCHFmult_GAM1_MD(nebfBE,npebfBE,npbf,ng1,ng1prm,nat,
     x                        ngtg1,pmass,cat,zan,bcoef1,gamma1,
     x                        AMPEB2C_be,AGEBFCC_be,AGNBFCC,ELCEX_be,
     x                        NUCEX,ELCAM_be,NUCAM,ELCBFC_be,NUCBFC)
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Calculate total number of three-particle integrals
      ng2_int=nebf*nebf*nebfBE*nebfBE*npbf*npbf
      ng2_xchf=nebfBE*nebfBE*nebfBE*nebfBE*npbf*npbf

! Have each process calculate and store its ng2/proc integrals
! Assign residual integrals to last process
      if (rank.eq.(nproc-1)) then
       dimINT2=ng2_int/nproc+mod(ng2_int,nproc)
       dimXCHF2=ng2_xchf/nproc+mod(ng2_xchf,nproc)
       if (LADDEXCH) dimINT2ex=ng2_int/nproc+mod(ng2_int,nproc)
      else
       dimINT2=ng2_int/nproc
       dimXCHF2=ng2_xchf/nproc
       if (LADDEXCH) dimINT2ex=ng2_int/nproc
      end if

! These arrays now local to each process and correspond to only its chunk of integrals
      if(allocated(INT_GAM2)) deallocate(INT_GAM2)
      allocate(INT_GAM2(dimINT2))
      if(allocated(INT_GAM2ex)) deallocate(INT_GAM2ex)
      allocate(INT_GAM2ex(dimINT2ex))
      if(allocated(XCHF_GAM2)) deallocate(XCHF_GAM2)
      allocate(XCHF_GAM2(dimXCHF2))
      if(allocated(XCHF_GAM2s)) deallocate(XCHF_GAM2s)
      allocate(XCHF_GAM2s(dimXCHF2))

      if (read_GAM2) then

! Read chunk with process (restart job requires same number of processes)
       if (rank.eq.0) then
        write(*,*)
        write(*,*) "---------------------------"
        write(*,*) " Reading:         INT_GAM2 "
       end if
       call RXCHF_readint_mpi(nproc,rank,dimINT2,17,
     x                        "INT_GAM2-"//istring//".ufm",
     x                        INT_GAM2)

       if (LADDEXCH) then
        if(rank.eq.0) write(*,*) "                INT_GAM2ex "
        call RXCHF_readint_mpi(nproc,rank,dimINT2ex,19,
     x                         "INT_GAM2ex-"//istring//".ufm",
     x                         INT_GAM2ex)
       end if

       if (nbe.gt.1) then
        if (rank.eq.0) then
         write(*,*) "                 XCHF_GAM2 "
         write(*,*) "                XCHF_GAM2s "
        end if
        call RXCHF_readint_mpi(nproc,rank,dimXCHF2,18,
     x                         "XCHF_GAM2-"//istring//".ufm",
     x                         XCHF_GAM2)
        call RXCHF_readint_mpi(nproc,rank,dimXCHF2,19,
     x                         "XCHF_GAM2s-"//istring//".ufm",
     x                         XCHF_GAM2s)
       end if

       if (rank.eq.0) then
        write(*,*) "---------------------------"
        write(*,*)
       end if

      else

       if (rank.eq.0) then
        write(*,*)
        write(*,*) "---------------------------"
        write(*,*) " Calculating:     INT_GAM2 "
        if (LADDEXCH) then
         write(*,*) "                INT_GAM2ex "
        end if
        write(*,*) "                 XCHF_GAM2 "
        write(*,*) "                XCHF_GAM2s "
        write(*,*) "---------------------------"
        write(*,*)
       end if

       if (LADDEXCH) then
        call RXCHF_GAM2ex_MPI(nproc,rank,
     x                        ng2chk,nebf,npebf,npbf,
     x                        ng2_int,dimINT2,ng2prm,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        XCHF_GAM2,INT_GAM2,INT_GAM2ex,XCHF_GAM2s)
       else
        call RXCHF_GAM2_MPI(nproc,rank,
     x                      ng2chk,nebf,npebf,nebfBE,npebfBE,npbf,
     x                      ng2_int,ng2_xchf,dimINT2,dimXCHF2,
     x                      nat,ngtg1,pmass,cat,zan,bcoef1,gamma1,
     x                      KPESTR,KPEEND,
     x                      AMPEB2C,AGEBFCC,
     x                      ELCEX,ELCAM,ELCBFC,
     x                      KPESTR_be,KPEEND_be,
     x                      AMPEB2C_be,AGEBFCC_be,
     x                      ELCEX_be,ELCAM_be,ELCBFC_be,
     x                      AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                      XCHF_GAM2,INT_GAM2,XCHF_GAM2s)
       end if

! Write integrals to disk (chunk calculated by each process)
       call RXCHF_writeint_mpi(nproc,rank,dimINT2,17,
     x                         "INT_GAM2-"//istring//".ufm",
     x                         INT_GAM2)
       if(rank.eq.0) write(*,*) "INT_GAM2 written to disk"

       if (LADDEXCH) then
        call RXCHF_writeint_mpi(nproc,rank,dimINT2ex,19,
     x                          "INT_GAM2ex-"//istring//".ufm",
     x                          INT_GAM2ex)
        if(rank.eq.0) write(*,*) "INT_GAM2ex written to disk"
       end if

      end if

      if (nbe.le.1) then

        if (rank.eq.0) then
         write(*,*)
         write(*,*) "NBE = ",NBE," <= 1"
         write(*,*) "so deallocating XCHF_GAM2"
         write(*,*)
        end if

        dimXCHF2=1
        if(allocated(XCHF_GAM2s)) deallocate(XCHF_GAM2s)
        if(allocated(XCHF_GAM2)) deallocate(XCHF_GAM2)
        allocate(XCHF_GAM2(dimXCHF2))
        allocate(XCHF_GAM2s(dimXCHF2))

        if(allocated(INT_GAM3)) deallocate(INT_GAM3)
        if(allocated(INT_GAM3ex1)) deallocate(INT_GAM3ex1)
        if(allocated(INT_GAM3ex2)) deallocate(INT_GAM3ex2)
        if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
        if(allocated(int_pinds)) deallocate(int_pinds)
        if(allocated(int_einds)) deallocate(int_einds)
        if(allocated(INT_GAM4)) deallocate(INT_GAM4)
        if(allocated(xchf_pinds)) deallocate(xchf_pinds)
        if(allocated(xchf_einds)) deallocate(xchf_einds)
        if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
        allocate(INT_GAM3(dimINT3))
        allocate(INT_GAM3ex1(dimINT3ex))
        allocate(INT_GAM3ex2(dimINT3ex))
        allocate(XCHF_GAM3(dimXCHF3))
        allocate(int_pinds(dimINT4))
        allocate(int_einds(dimINT4))
        allocate(INT_GAM4(dimINT4))
        allocate(xchf_pinds(dimXCHF4))
        allocate(xchf_einds(dimXCHF4))
        allocate(XCHF_GAM4(dimXCHF4))

      else

! Write integrals to disk (chunk calculated by each process)
        if (.not.(read_GAM2)) then
         call RXCHF_writeint_mpi(nproc,rank,dimXCHF2,18,
     x                           "XCHF_GAM2-"//istring//".ufm",
     x                           XCHF_GAM2)
         if(rank.eq.0) write(*,*) "XCHF_GAM2 written to disk"

         call RXCHF_writeint_mpi(nproc,rank,dimXCHF2,19,
     x                           "XCHF_GAM2s-"//istring//".ufm",
     x                           XCHF_GAM2s)
         if(rank.eq.0) write(*,*) "XCHF_GAM2s written to disk"
        end if

C nbe >= 2
C  - calculate four-particle integrals and store in memory
C     => INT_GAM3
C     => INT_GAM3ex only needed if RXCHF-ae
C     => XCHF_GAM3 only needed if nbe >= 3

! Calculate total number of four-particle integrals
        ng3_int=nebf*nebf*nebfBE*nebfBE*nebfBE*nebfBE*npbf*npbf
        ng3_xchf=nebfBE*nebfBE*nebfBE*nebfBE*nebfBE*nebfBE*npbf*npbf

! Have each process calculate and store its ng3/proc integrals
! Assign residual integrals to last process
        if (rank.eq.(nproc-1)) then
         dimINT3=ng3_int/nproc+mod(ng3_int,nproc)
         dimXCHF3=ng3_xchf/nproc+mod(ng3_xchf,nproc)
         if (LADDEXCH) dimINT3ex=ng3_int/nproc+mod(ng3_int,nproc)
        else
         dimINT3=ng3_int/nproc
         dimXCHF3=ng3_xchf/nproc
         if (LADDEXCH) dimINT3ex=ng3_int/nproc
        end if

        if (read_GAM3) then

! Read chunk with process (restart job requires same number of processes)

         if(allocated(INT_GAM3)) deallocate(INT_GAM3)
         allocate(INT_GAM3(dimINT3))

         if (rank.eq.0) then
          write(*,*)
          write(*,*) "---------------------------"
          write(*,*) " Reading:         INT_GAM3 "
         end if
         call RXCHF_readint_mpi(nproc,rank,dimINT3,17,
     x                          "INT_GAM3-"//istring//".ufm",
     x                          INT_GAM3)

         if (LADDEXCH) then

          if(allocated(INT_GAM3ex1)) deallocate(INT_GAM3ex1)
          allocate(INT_GAM3ex1(dimINT3ex))
          if(allocated(INT_GAM3ex2)) deallocate(INT_GAM3ex2)
          allocate(INT_GAM3ex2(dimINT3ex))

          if (rank.eq.0) then
           write(*,*) "               INT_GAM3ex1 "
           write(*,*) "               INT_GAM3ex2 "
          end if
          call RXCHF_readint_mpi(nproc,rank,dimINT3ex,20,
     x                           "INT_GAM3ex1-"//istring//".ufm",
     x                           INT_GAM3ex1)
          call RXCHF_readint_mpi(nproc,rank,dimINT3ex,20,
     x                           "INT_GAM3ex2-"//istring//".ufm",
     x                           INT_GAM3ex2)
         end if

         if (nbe.gt.2) then

          if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
          allocate(XCHF_GAM3(dimXCHF3))

          if(rank.eq.0) write(*,*) "                 XCHF_GAM3 "
          call RXCHF_readint_mpi(nproc,rank,dimXCHF3,18,
     x                           "XCHF_GAM3-"//istring//".ufm",
     x                           XCHF_GAM3)
         end if

         if (rank.eq.0) then
          write(*,*) "---------------------------"
          write(*,*)
         end if

        else ! not read_gam3, so calculate

         if (rank.eq.0) then
          write(*,*)
          write(*,*) "---------------------------"
          write(*,*) " Calculating:     INT_GAM3 "
          if (LADDEXCH) then
           write(*,*) "               INT_GAM3ex1 "
           write(*,*) "               INT_GAM3ex2 "
          end if
          if(nbe.gt.2) then
           write(*,*) "                 XCHF_GAM3 "
          end if
          write(*,*) "---------------------------"
          write(*,*)
         end if

         call get_mpi_range(ng3_int,nblocks,
     x                      blockrank,blockstart,blockend)
         if (blockrank.eq.(nblocks-1)) then
          ng3block=ng3_int/nblocks+mod(ng3_int,nblocks)
          blockend=ng3_int
         else
          ng3block=ng3_int/nblocks
         end if
         if(rank.eq.0) then
          write(*,*) "Computing block ",blockrank," of ",nblocks,
     x               " total blocks"
          write(*,*) "ng3,blockstart,blockend:",
     x                ng3_int,blockstart,blockend
         end if
         if (rank.eq.(nproc-1)) then
          dimINT3=ng3block/nproc+mod(ng3block,nproc)

          dimXCHF3=ng3block/nproc+mod(ng3block,nproc) ! use ng3_int range for XCHF if -ae, override below if -ne
          if (LADDEXCH) dimINT3ex=ng3block/nproc+mod(ng3block,nproc)
         else
          dimINT3=ng3block/nproc
          dimXCHF3=ng3block/nproc
          if (LADDEXCH) dimINT3ex=ng3block/nproc
         end if

         if(allocated(INT_GAM3)) deallocate(INT_GAM3)
         allocate(INT_GAM3(dimINT3))

         if(LADDEXCH) then
          if(allocated(INT_GAM3ex1)) deallocate(INT_GAM3ex1)
          allocate(INT_GAM3ex1(dimINT3ex))
          if(allocated(INT_GAM3ex2)) deallocate(INT_GAM3ex2)
          allocate(INT_GAM3ex2(dimINT3ex))
          if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
          allocate(XCHF_GAM3(dimXCHF3))
         end if

! Only symmetrize XCHF_GAM3 if nbe > 2 (saves significant MPI comm, only applicable to -ae)
         if (nbe.le.2) then
          LXCHFSYMM=.false.
         else
          LXCHFSYMM=.true.
         end if

         if (LADDEXCH) then
          call RXCHF_GAM3ex_MPI(nproc,rank,
     x                          nblocks,blockrank,
     x                          blockstart,blockend,
     x                          ng3chk,nebf,npebf,npbf,
     x                          ng3block,dimINT3,ng3prm,nat,ngtg1,
     x                          pmass,cat,zan,bcoef1,gamma1,
     x                          KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                          ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                          LXCHFSYMM,XCHF_GAM3,INT_GAM3,
     x                          INT_GAM3ex1,INT_GAM3ex2)
         else
          call RXCHF_GAM3_INT_MPI(nproc,rank,
     x                            nblocks,blockrank,
     x                            blockstart,blockend,
     x                            ng3chk,
     x                            nebf,npebf,nebfBE,npebfBE,npbf,
     x                            ng3block,dimINT3,ng3prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,
     x                            AMPEB2C,AGEBFCC,
     x                            ELCEX,ELCAM,ELCBFC,
     x                            KPESTR_be,KPEEND_be,
     x                            AMPEB2C_be,AGEBFCC_be,
     x                            ELCEX_be,ELCAM_be,ELCBFC_be,
     x                            AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                            INT_GAM3)
         end if

! Write integrals to disk (chunk calculated by each process)
         call RXCHF_writeint_mpi(nproc,rank,dimINT3,17,
     x                           "INT_GAM3-"//istring//".ufm",
     x                           INT_GAM3)
         if(rank.eq.0) write(*,*) "INT_GAM3 written to disk"

         if (LADDEXCH) then
          call RXCHF_writeint_mpi(nproc,rank,dimINT3ex,20,
     x                            "INT_GAM3ex1-"//istring//".ufm",
     x                            INT_GAM3ex1)
          if(rank.eq.0) write(*,*) "INT_GAM3ex1 written to disk"

          call RXCHF_writeint_mpi(nproc,rank,dimINT3ex,20,
     x                            "INT_GAM3ex2-"//istring//".ufm",
     x                            INT_GAM3ex2)
          if(rank.eq.0) write(*,*) "INT_GAM3ex2 written to disk"
         end if

        end if

        if (nbe.le.2) then

          if(LADDEXCH) then
           if (rank.eq.0) then
            write(*,*)
            write(*,*) "NBE = ",NBE," <= 2"
            write(*,*) "so deallocating XCHF_GAM3"
            write(*,*)
           end if
          end if

          dimXCHF3=1
          if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
          allocate(XCHF_GAM3(dimXCHF3))

          if(allocated(int_pinds)) deallocate(int_pinds)
          if(allocated(int_einds)) deallocate(int_einds)
          if(allocated(INT_GAM4)) deallocate(INT_GAM4)
          if(allocated(xchf_pinds)) deallocate(xchf_pinds)
          if(allocated(xchf_einds)) deallocate(xchf_einds)
          if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
          allocate(int_pinds(dimINT4))
          allocate(int_einds(dimINT4))
          allocate(INT_GAM4(dimINT4))
          allocate(xchf_pinds(dimXCHF4))
          allocate(xchf_einds(dimXCHF4))
          allocate(XCHF_GAM4(dimXCHF4))

        else

          if(.not.(read_GAM3.or.LADDEXCH)) then
           if (rank.eq.0) then
            write(*,*)
            write(*,*) "---------------------------"
            write(*,*) " Calculating:    XCHF_GAM3 "
            write(*,*) "---------------------------"
            write(*,*)
           end if
          end if

! If -ne, calculate XCHF_GAM3 separately using XCHF-specific block info
          if(.not.(LADDEXCH)) then

           call get_mpi_range(ng3_xchf,nblocks,
     x                        blockrank,blockstart,blockend)
           if (blockrank.eq.(nblocks-1)) then
            ng3block=ng3_xchf/nblocks+mod(ng3_xchf,nblocks)
            blockend=ng3_xchf
           else
            ng3block=ng3_xchf/nblocks
           end if
           if(rank.eq.0) then
            write(*,*) "Computing block ",blockrank," of ",nblocks,
     x                 " total blocks"
            write(*,*) "ng3,blockstart,blockend:",
     x                  ng3_xchf,blockstart,blockend
           end if
           if (rank.eq.(nproc-1)) then
            dimXCHF3=ng3block/nproc+mod(ng3block,nproc)
           else
            dimXCHF3=ng3block/nproc
           end if

           if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
           allocate(XCHF_GAM3(dimXCHF3))

           call RXCHF_GAM3_XCHF_MPI(nproc,rank,
     x                              nblocks,blockrank,
     x                              blockstart,blockend,
     x                              ng3chk,nebfBE,npebfBE,npbf,
     x                              ng3block,dimXCHF3,ng3prm,nat,ngtg1,
     x                              pmass,cat,zan,bcoef1,gamma1,
     x                              KPESTR_be,KPEEND_be,
     x                              AMPEB2C_be,AGEBFCC_be,
     x                              ELCEX_be,ELCAM_be,ELCBFC_be,
     x                              AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                              XCHF_GAM3)

          end if

! Write integrals to disk (chunk calculated by each process)
! Writes XCHF_GAM3 whether it was calculated previously (-ae) or just before (-ne)
! with appropriate dimXCHF
          if (.not.(read_GAM3)) then
           call RXCHF_writeint_mpi(nproc,rank,dimXCHF3,18,
     x                             "XCHF_GAM3-"//istring//".ufm",
     x                             XCHF_GAM3)
           if(rank.eq.0) write(*,*) "XCHF_GAM3 written to disk"
          end if

C Only process five-particle integrals if not just an integral calculation run
          if(nblocks.eq.1) then

C nbe >= 3
C  - calculate interaction five-particle integrals and store in memory
C     => INT_GAM4
C  - assume there is enough memory for each process to store GAM2s

! Calculate total number of four-particle integrals
            ng4_int=nebf*nebf*nebfBE*nebfBE*
     x              nebfBE*nebfBE*nebfBE*nebfBE*npbf*npbf
            ng4_xchf=nebfBE*nebfBE*nebfBE*nebfBE*
     x               nebfBE*nebfBE*nebfBE*nebfBE*npbf*npbf

            if (.not.(read_GAM4)) then

              if (rank.eq.0) then
               write(*,*)
               write(*,*) "---------------------------"
               write(*,*) " Calculating:     INT_GAM4 "
               write(*,*) "---------------------------"
               write(*,*)
              end if

! Dimension of all INT_GAM4 integrals
              if (rank.eq.(nproc-1)) then
               dimINT4=ng4_int/nproc+mod(ng4_int,nproc)
              else
               dimINT4=ng4_int/nproc
              end if
 
! Calculate INT_GAM4 (routine writes nontrivial INT_GAM4 to disk)
              call RXCHF_GAM4_INT_MPI(nproc,rank,
     x                                ng4chk,nebf,nebfBE,npbf,ngee,
     x                                ng2_xchf,ng4_int,
     x                                dimXCHF2,dimINT4,
     x                                XCHF_GAM2s)

              if(nbe.le.3) then

                if(allocated(xchf_pinds)) deallocate(xchf_pinds)
                if(allocated(xchf_einds)) deallocate(xchf_einds)
                if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
                allocate(xchf_pinds(dimXCHF4))
                allocate(xchf_einds(dimXCHF4))
                allocate(XCHF_GAM4(dimXCHF4))

              else

C nbe >= 4
C  - calculate XCHF five-particle integrals and store in memory
C     => XCHF_GAM4
C  - assume there is enough memory for each process to store GAM2s

               if (rank.eq.0) then
                write(*,*)
                write(*,*) "---------------------------"
                write(*,*) " Calculating:    XCHF_GAM4 "
                write(*,*) "---------------------------"
                write(*,*)
               end if

! Dimension of all XCHF_GAM4 integrals
               if (rank.eq.(nproc-1)) then
                dimXCHF4=ng4_xchf/nproc+mod(ng4_xchf,nproc)
               else
                dimXCHF4=ng4_xchf/nproc
               end if

! Calculate XCHF_GAM4 (routine writes nontrivial GAM4s to disk)
               call RXCHF_GAM4_XCHF_MPI(nproc,rank,
     x                                  ng4chk,nebf,nebfBE,npbf,ngee,
     x                                  ng2_xchf,ng4_xchf,
     x                                  dimXCHF2,dimXCHF4,
     x                                  XCHF_GAM2s)

              end if ! nbe >= 3

            end if ! not read_GAM4

! Always read nontrivial INT_GAM4

            if (rank.eq.0) then
             write(*,*)
             write(*,*) "---------------------------"
             write(*,*) " Reading:         INT_GAM4 "
             write(*,*) "---------------------------"
             write(*,*)
            end if

!            call RXCHF_readint_mpi(nproc,rank,dimINT4,17,
!     x                             "INT_GAM4-"//istring//".ufm",
!     x                             INT_GAM4)

            open(unit=unitno,file="INT_GAM4-"//istring//".ufm",
     x           form="unformatted",action="read")
            read(unitno) dimINT4 ! # of nontrivial INT_GAM4 integrals

            if(allocated(int_pinds)) deallocate(int_pinds)
            allocate(int_pinds(dimINT4))
            if(allocated(int_einds)) deallocate(int_einds)
            allocate(int_einds(dimINT4))
            if(allocated(INT_GAM4)) deallocate(INT_GAM4)
            allocate(INT_GAM4(dimINT4))
            
            call RXCHF_readbit_mpi(nproc,rank,dimINT4,unitno,
     x                             int_pinds,int_einds,INT_GAM4)

            close(unitno)

            if (nbe.le.3) then

! XCHF_GAM4 not needed for nbe <= 3
              dimXCHF4=1
              if(allocated(xchf_pinds)) deallocate(xchf_pinds)
              allocate(xchf_pinds(dimXCHF4))
              if(allocated(xchf_einds)) deallocate(xchf_einds)
              allocate(xchf_einds(dimXCHF4))
              if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
              allocate(XCHF_GAM4(dimXCHF4))

            else

! Read nontrivial XCHF_GAM4 for nbe > 3

              if (rank.eq.0) then
               write(*,*)
               write(*,*) "---------------------------"
               write(*,*) " Reading:        XCHF_GAM4 "
               write(*,*) "---------------------------"
               write(*,*)
              end if

!              call RXCHF_readint_mpi(nproc,rank,dimXCHF4,18,
!     x                               "XCHF_GAM4-"//istring//".ufm",
!     x                               XCHF_GAM4)

              open(unit=unitno,file="XCHF_GAM4-"//istring//".ufm",
     x             form="unformatted",action="read")
              read(unitno) dimXCHF4 ! # of nontrivial XCHF_GAM4 integrals

              if(allocated(xchf_pinds)) deallocate(xchf_pinds)
              allocate(xchf_pinds(dimXCHF4))
              if(allocated(xchf_einds)) deallocate(xchf_einds)
              allocate(xchf_einds(dimXCHF4))
              if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
              allocate(XCHF_GAM4(dimXCHF4))

              call RXCHF_readbit_mpi(nproc,rank,dimXCHF4,unitno,
     x                               xchf_pinds,xchf_einds,XCHF_GAM4)

              close(unitno)

            end if ! nbe >= 4

          end if ! nblocks == 1

        end if ! nbe >= 3

      end if ! nbe >= 2

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      wtime1 = MPI_WTIME() - wtime

C Kick-off SCF
      wtime  = MPI_WTIME()

      nebf2=nebf*nebf
      nebflt=nebf*(nebf+1)/2

      nebfBE2=nebfBE*nebfBE
      nebfBElt=nebfBE*(nebfBE+1)/2

      npbf2=npbf*npbf
      npbflt=npbf*(npbf+1)/2

      if ((read_GAM3).or.(nblocks.eq.1)) then
       if(LRXCUHF) then
        call RXCHF_uscf_MPI(nproc,rank,
     x                      nelec,nae,nae_alp,nae_beta,
     x                      nbe,nebflt,npbflt,nucst,
     x                      npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                      ngtg1,ng1,ng2,ng3,ng4,
     x                      NG2CHK,NG3CHK,NG4CHK,
     x                      read_CE,read_CP,
     x                      LG4DSCF,LG3DSCF,LG2DSCF,
     x                      LSOSCF,LDIIS,LSORXCHF,
     x                      OCBSETYP,LDBG,LADDEXCH,
     x                      ng2prm,ng3prm,nat,pmass,cat,zan,
     x                      bcoef1,gamma1,
     x                      KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                      LG2IC,dimXCHF2,dimINT2,dimINT2ex,
     x                      XCHF_GAM2,INT_GAM2,
     x                      INT_GAM2ex,XCHF_GAM2s,
     x                      LG3IC,dimXCHF3,dimINT3,dimINT3ex,
     x                      XCHF_GAM3,INT_GAM3,
     x                      INT_GAM3ex1,INT_GAM3ex2,
     x                      LG4IC,dimXCHF4,dimINT4,
     x                      int_pinds,int_einds,INT_GAM4,
     x                      xchf_pinds,xchf_einds,XCHF_GAM4)
       else
        call RXCHF_scf_MPI(nproc,rank,
     x                     nelec,nae,nbe,nucst,
     x                     npebf,nebf,nebf2,nebflt,
     x                     npebfBE,nebfBE,nebfBE2,nebfBElt,elindBE,
     x                     npbf,npbf2,npbflt,
     x                     ngtg1,ngee,ng1,
     x                     ng2_int,ng2_xchf,
     x                     ng3_int,ng3_xchf,
     x                     ng4_int,ng4_xchf,
     x                     NG2CHK,NG3CHK,NG4CHK,
     x                     read_CE,read_CP,
     x                     LG4DSCF,LG3DSCF,LG2DSCF,
     x                     LSOSCF,LDIIS,LSORXCHF,
     x                     OCBSETYP,LDBG,LADDEXCH,
     x                     LALTBAS,LREORDGS,
     x                     nat,pmass,cat,zan,
     x                     bcoef1,gamma1,
     x                     KPESTR,KPEEND,
     x                     AMPEB2C,AGEBFCC,
     x                     ELCEX,ELCAM,ELCBFC,
     x                     KPESTR_be,KPEEND_be,
     x                     AMPEB2C_be,AGEBFCC_be,
     x                     ELCEX_be,ELCAM_be,ELCBFC_be,
     x                     AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                     LG2IC,dimXCHF2,dimINT2,dimINT2ex,
     x                     XCHF_GAM2,INT_GAM2,INT_GAM2ex,XCHF_GAM2s,
     x                     LG3IC,dimXCHF3,dimINT3,dimINT3ex,
     x                     XCHF_GAM3,INT_GAM3,
     x                     INT_GAM3ex1,INT_GAM3ex2,
     x                     LG4IC,dimXCHF4,dimINT4,
     x                     int_pinds,int_einds,INT_GAM4,
     x                     xchf_pinds,xchf_einds,XCHF_GAM4)
       end if
      end if

      wtime2 = MPI_WTIME() - wtime

C Cleanup
      if(allocated(XCHF_GAM4))   deallocate(XCHF_GAM4)
      if(allocated(xchf_einds))  deallocate(xchf_einds)
      if(allocated(xchf_pinds))  deallocate(xchf_pinds)
      if(allocated(INT_GAM4))    deallocate(INT_GAM4)
      if(allocated(int_einds))   deallocate(int_einds)
      if(allocated(int_pinds))   deallocate(int_pinds)
      if(allocated(XCHF_GAM3))   deallocate(XCHF_GAM3)
      if(allocated(INT_GAM3ex2)) deallocate(INT_GAM3ex2)
      if(allocated(INT_GAM3ex1)) deallocate(INT_GAM3ex1)
      if(allocated(INT_GAM3))    deallocate(INT_GAM3)
      if(allocated(XCHF_GAM2s))  deallocate(XCHF_GAM2s)
      if(allocated(XCHF_GAM2))   deallocate(XCHF_GAM2)
      if(allocated(INT_GAM2ex))  deallocate(INT_GAM2ex)
      if(allocated(INT_GAM2))    deallocate(INT_GAM2)

      if(allocated(ELCBFC_be))   deallocate(ELCBFC_be)
      if(allocated(ELCAM_be))    deallocate(ELCAM_be)
      if(allocated(AGEBFCC_be))  deallocate(AGEBFCC_be)
      if(allocated(ELCEX_be))    deallocate(ELCEX_be)
      if(allocated(AMPEB2C_be))  deallocate(AMPEB2C_be)

C Print timing summary
      if (rank.eq.0) then
       write(*,*)
       write(*,*) "FINISHED RXCHFMULT CALCULATION"
       write(*,*)
       write(*,3000) wtime1,wtime2
       write(*,*)
      end if


 1000 format(1X,I3,I6,I5,I3,I3,F12.6,F10.6,F10.6,F10.6,F10.6)
 2000 format(1X,I3,I6,I5)
 3000 FORMAT(/8X,'  +--------------------------------------+',/,
     X        8X,'  |    TIMING SUMMARY FOR CALCULATION    |',/,
     x        8X,'  +--------------------------------------+',/,
     x        8X,'    TIME TO EVALUATE INTEGRALS:',1X,F12.4/
     x        8X,'                  TIME FOR SCF:',1X,F12.4/)

      return
      end
C=======================================================================
      subroutine RXCHF_readint_mpi(nproc,rank,n,namelen,fname,arr)

C Reads integrals from file "fname" into an n-dimensional array [arr]
C Integrals must be written to an unformatted file in chunks of size
C nchunk 
C=======================================================================
      implicit none

C Input variables
      integer nproc,rank
      integer namelen,n
      character(len=namelen) fname
C Output variables
      double precision arr(n)
C Local variables
      integer i
      integer unitno
      integer nchunk
      parameter(nchunk=100000000)

      unitno=200+rank

C Initialize
      arr=0.0d+00

C Read in from file
      open(unit=unitno,file=fname,form="unformatted")
      do i=1,n/nchunk
        read(unitno) arr((i-1)*nchunk+1:i*nchunk)
      end do
      if(mod(n,nchunk).ne.0) then
        read(unitno) arr((n/nchunk)*nchunk+1:n)
      end if
      close(unitno)

      return
      end

C=======================================================================
      subroutine RXCHF_readbit_mpi(nproc,rank,n,unitno,pinds,einds,arr)

C Reads integrals from already-opened file with unit number [unitno]
C into an n-dimensional array [arr]
C Also stored are n-dimensional proton indices [ip jp] and
C n-dimensional electron indices [i1 j1 i2 j2 i3 j3 i4 j4]
C Indices and integrals must be written to an unformatted file
C in chunks of size [nchunk]
C=======================================================================
      implicit none

C Input variables
      integer nproc,rank
      integer n,unitno
C Output variables
      integer(kind=4)  pinds(n)
      integer(kind=8)  einds(n)
      double precision arr(n)
C Local variables
      integer i
      integer nchunk
      parameter(nchunk=100000000)

C Initialize
      pinds=0
      einds=0
      arr=0.0d+00

C Read packed proton indices
      do i=1,n/nchunk
        read(unitno) pinds((i-1)*nchunk+1:i*nchunk)
      end do
      if(mod(n,nchunk).ne.0) then
        read(unitno) pinds((n/nchunk)*nchunk+1:n)
      end if

C Read packed electron indices
      do i=1,n/nchunk
        read(unitno) einds((i-1)*nchunk+1:i*nchunk)
      end do
      if(mod(n,nchunk).ne.0) then
        read(unitno) einds((n/nchunk)*nchunk+1:n)
      end if

C Read integrals
      do i=1,n/nchunk
        read(unitno) arr((i-1)*nchunk+1:i*nchunk)
      end do
      if(mod(n,nchunk).ne.0) then
        read(unitno) arr((n/nchunk)*nchunk+1:n)
      end if

      return
      end

C=======================================================================
      subroutine RXCHF_writeint_mpi(nproc,rank,n,namelen,fname,arr)

C Writes integrals from n-dimensional array [arr] to file "fname"
C Integrals must be read from an unformatted file in chunks of size
C nchunk 
C=======================================================================
      implicit none

C Input variables
      integer nproc,rank
      integer namelen,n
      character(len=namelen) fname
      double precision arr(n)
C Local variables
      integer i
      integer unitno
      integer nchunk
      parameter(nchunk=100000000)

      unitno=200+rank

C Write to file
      open(unit=unitno,file=fname,form="unformatted")
      do i=1,n/nchunk
        write(unitno) arr((i-1)*nchunk+1:i*nchunk)
      end do
      if(mod(n,nchunk).ne.0) then
        write(unitno) arr((n/nchunk)*nchunk+1:n)
      end if
      close(unitno)

      return
      end

