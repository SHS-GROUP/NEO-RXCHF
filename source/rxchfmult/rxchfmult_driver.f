C=======================================================================
      subroutine RXCHFmult_driver(nelec,nae,nbe,nucst,
     x                            nebf,npebf,npbf,nat,ngtg1,ngee,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            nebfBE,elindBE,
     x                            ng2chk,ng3chk,ng4chk,
     x                            read_CE,read_CP,
     x                            read_GAM2,read_GAM3,read_GAM4,
     x                            LG2IC,LG3IC,LG4IC,
     x                            LG2DSCF,LG3DSCF,LG4DSCF,
     x                            LSOSCF,LOCBSE,LCMF,LADDEXCH)

C Driver to calculate RXCHF integrals for nbe > 1
C   XCHF_GAM* : integrals needed for XCHF contribution
C    INT_GAM* : integrals needed for interaction contribution
C  INT_GAM*ex : integrals needed for exchange contribution
C
C Maximum dimension of integral needed is nbe + 1
C=======================================================================
      implicit none
      include "omp_lib.h"

C Input variables
      integer nelec                    ! Total number of electrons
      integer nae                      ! Number of regular electrons
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
      integer ng1prm,ng2prm,ng3prm     ! Numbers of integrals (primitive)
      integer ng2chk,ng3chk,ng4chk     ! Number of chunks to split integral calculations

      integer ELCAM(npebf,3)           ! Angular mom for electrons
      integer NUCAM(npbf,3)            ! Angular mom for quantum nuclei
      integer AMPEB2C(npebf)           ! Map primitive index to contracted
      integer KPESTR(nebf)             ! Map contracted index to primitive start
      integer KPEEND(nebf)             ! Map contracted index to primitive end
      integer KPESTR_be(nebf)          ! Analog for NBE basis set
      integer KPEEND_be(nebf)          ! Analog for NBE basis set
      integer elindBE(nebfBE)          ! Contracted indices of NBE basis set

      logical LG2IC                    ! In-core storage of 3-particle integrals
      logical LG3IC                    ! In-core storage of 4-particle integrals
      logical LG4IC                    ! In-core storage of 5-particle integrals
      logical LG2DSCF                  ! Direct SCF for 3-particle integrals
      logical LG3DSCF                  ! Direct SCF for 4-particle integrals
      logical LG4DSCF                  ! Direct SCF for 5-particle integrals
      logical LSOSCF                   ! SOSCF where applicable
      logical LOCBSE                   ! OCBSE2 procedure
      logical LCMF                     ! On-the-fly Fock matrix check and debugging
      logical LADDEXCH                 ! Flag to activate approximate exchange
      logical read_CE,read_CP          ! Read in orbitals
      logical read_GAM2                ! 
      logical read_GAM3                ! Read in integrals
      logical read_GAM4                ! 
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
      integer i,j,istat
      integer contrind,primind
      integer numprims,primstart,primend
      integer currcontrind,currprimind
      integer nebf2,nebflt                  !
      integer nebfBE2,nebfBElt              ! Convenient quantities
      integer npbf2,npbflt                  !
      integer npra,nprb                     ! Number of distinct electron pairs
      logical LALTBAS                       ! Flag to denote distinct special electron basis
      logical laddbasis                     ! Working variable

      integer XELCAM(npebf,3)               ! 
      integer XAMPEB2C(npebf)               !
      integer XKPESTR(nebf)                 ! 
      integer XKPEEND(nebf)                 ! Dummy basis set variables
      double precision XELCEX(npebf)        ! 
      double precision XELCBFC(npebf,3)     ! 
      double precision XAGEBFCC(npebf)      ! 

      integer dimXCHF2                      !
      integer dimXCHF3                      ! Dimensions of XCHF integral arrays
      integer dimXCHF4                      !
      integer dimINT2                       !
      integer dimINT3                       ! Dimensions of interaction integral arrays
      integer dimINT4                       !
      integer dimINT2ex                     ! Dimensions of exchange integral arrays
      integer dimINT3ex                     !

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

      integer, allocatable ::          ELCAM_be(:,:)        !
      integer, allocatable ::          AMPEB2C_be(:)        ! 
      double precision, allocatable :: ELCEX_be(:)          ! Analogs for NBE basis set 
      double precision, allocatable :: ELCBFC_be(:,:)       ! 
      double precision, allocatable :: AGEBFCC_be(:)        ! 


C Reorder electronic basis set such that special electron subset is first
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
      write(*,*)'PRIM  CONT    ANG       EXPONENT CONTRACT  -X- -Y- -Z-'
      write(*,*)'INDEX INDEX   MOM                  COEF'
      do i=1,npebf
        write(*,1000) i,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),
     x                ELCAM(i,3),ELCEX(i),AGEBFCC(i),
     x                ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
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
      write(*,*)'PRIM  CONT    ANG       EXPONENT CONTRACT  -X- -Y- -Z-'
      write(*,*)'INDEX INDEX   MOM                  COEF'
      do i=1,npebfBE
        write(*,1000) i,AMPEB2C_be(i),ELCAM_be(i,1),ELCAM_be(i,2),
     x                ELCAM_be(i,3),ELCEX_be(i),AGEBFCC_be(i),
     x                ELCBFC_be(i,1),ELCBFC_be(i,2),ELCBFC_be(i,3)
      end do
      write(*,*)
      write(*,*)' CHECK CONTRACTED ELECTRONIC BASIS FUNCTIONS '
      write(*,*)'CONT INDEX    KPESTR     KPEEND'
      do i=1,nebfBE
         write(*,2000) i,KPESTR_be(i),KPEEND_be(i)
      end do

C If RXCHF-ae cannot have distinct regular/special electronic basis sets
      if (nebf.ne.nebfBE) then
       LALTBAS=.true.
      else
       LALTBAS=.false.
      end if
      if ((LALTBAS).and.(LADDEXCH)) then
       write(*,*) "Cannot have restricted special electronic basis set"
       write(*,*) "for RXCHF-ae. Exiting..."
       return
      end if

C Calculate inexpensive integrals usually calculated in main driver
C since the all-electron basis set was reordered
      call class_nuc_rep(nat,zan,cat)

      call elec_ovlap(npebf,nebf,nebf*nebf,
     x                AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

      call check_elec_ovlap(nebf)

      call nuc_ovlap(npbf,npbf*npbf,
     x               AGNBFCC,NUCEX,NUCAM,NUCBFC)

      call check_nuc_ovlap(npbf)

      call calc_GAM_epcore(nebf,npebf,npbf,nebf*nebf,npbf*npbf,
     x                     nat,pmass,zan,cat,
     x                     AMPEB2C,AGEBFCC,AGNBFCC,
     x                     ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      call calc_GAM_ee(nebf,npebf,ngee,
     x                 AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

C Calculate geminal integrals
      wtime = omp_get_wtime()

C Initialize dimensions
      dimXCHF2  = 1
      dimXCHF3  = 1
      dimXCHF4  = 1
      dimINT2   = 1
      dimINT3   = 1
      dimINT4   = 1
      dimINT2ex = 1
      dimINT3ex = 1

C nbe >= 2
C  - calculate two-particle integrals and write to disk
C     => XCHF_GAM1
C     => XCHF_GAM1s
C  - calculate three-particle integrals and store in memory
C     => INT_GAM2
C     => INT_GAM2ex only needed if RXCHF-ae
C     => XCHF_GAM2
C     => XCHF_GAM2s
C  - calculate four-particle integrals and store in memory
C     => INT_GAM3
C     => INT_GAM3ex1 only needed if RXCHF-ae
C     => INT_GAM3ex2 only needed if RXCHF-ae

      write(*,*)
      write(*,*) "---------------------------"
      write(*,*) " Calculating:   XCHF_GAM1  "
      write(*,*) "---------------------------"
      write(*,*)

      ng1=nebfBE*nebfBE*npbf*npbf
      ng1prm=npebfBE*npebfBE*npbf*npbf
      call RXCHFmult_GAM1_MD(nebfBE,npebfBE,npbf,ng1,ng1prm,nat,
     x                       ngtg1,pmass,cat,zan,bcoef1,gamma1,
     x                       AMPEB2C_be,AGEBFCC_be,AGNBFCC,ELCEX_be,
     x                       NUCEX,ELCAM_be,NUCAM,ELCBFC_be,NUCBFC)

      dimINT2=nebf*nebf*nebfBE*nebfBE*npbf*npbf
      if (LADDEXCH) dimINT2ex=dimINT2
      dimXCHF2=nebfBE*nebfBE*nebfBE*nebfBE*npbf*npbf

      if(allocated(INT_GAM2)) deallocate(INT_GAM2)
      allocate(INT_GAM2(dimINT2))
      if(allocated(INT_GAM2ex)) deallocate(INT_GAM2ex)
      allocate(INT_GAM2ex(dimINT2ex))
      if(allocated(XCHF_GAM2)) deallocate(XCHF_GAM2)
      allocate(XCHF_GAM2(dimXCHF2))
      if(allocated(XCHF_GAM2s)) deallocate(XCHF_GAM2s)
      allocate(XCHF_GAM2s(dimXCHF2))

      if (read_GAM2) then

       write(*,*)
       write(*,*) "---------------------------"

       write(*,*) " Reading:         INT_GAM2 "
       call RXCHFmult_readint(12,dimINT2,"INT_GAM2.ufm",INT_GAM2)

       if (LADDEXCH) then
        write(*,*) "                INT_GAM2ex "
        call RXCHFmult_readint(14,dimINT2ex,"INT_GAM2ex.ufm",INT_GAM2ex)
       end if

       if (nbe.gt.1) then
        write(*,*) "                 XCHF_GAM2 "
        write(*,*) "                XCHF_GAM2s "
        call RXCHFmult_readint(13,dimXCHF2,"XCHF_GAM2.ufm",XCHF_GAM2)
        call RXCHFmult_readint(14,dimXCHF2,"XCHF_GAM2s.ufm",XCHF_GAM2s)
       end if

       write(*,*) "---------------------------"
       write(*,*)

      else

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

       if (LADDEXCH) then
        call RXCHFmult_GAM2_IC1ex(ng2chk,nebf,npebf,npbf,
     x                            dimINT2,ng2prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            XCHF_GAM2,INT_GAM2,
     x                            INT_GAM2ex,XCHF_GAM2s)
       else
        call RXCHFmult_GAM2(ng2chk,
     x                      nebf,npebf,nebfBE,npebfBE,npbf,
     x                      dimINT2,dimXCHF2,nat,ngtg1,
     x                      pmass,cat,zan,bcoef1,gamma1,
     x                      KPESTR,KPEEND,
     x                      AMPEB2C,AGEBFCC,
     x                      ELCEX,ELCAM,ELCBFC,
     x                      KPESTR_be,KPEEND_be,
     x                      AMPEB2C_be,AGEBFCC_be,
     x                      ELCEX_be,ELCAM_be,ELCBFC_be,
     x                      AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                      INT_GAM2,XCHF_GAM2,XCHF_GAM2s)
       end if

       open(unit=20,file="INT_GAM2.ufm",form="unformatted")
       write(20) INT_GAM2
       close(20)
       write(*,*) "INT_GAM2 written to disk"
       if (LADDEXCH) then
        open(unit=21,file="INT_GAM2ex.ufm",form="unformatted")
        write(21) INT_GAM2ex
        close(21)
        write(*,*) "INT_GAM2ex written to disk"
       end if
       open(unit=22,file="XCHF_GAM2.ufm",form="unformatted")
       write(22) XCHF_GAM2
       close(22)
       write(*,*) "XCHF_GAM2 written to disk"
       open(unit=23,file="XCHF_GAM2s.ufm",form="unformatted")
       write(23) XCHF_GAM2s
       close(23)
       write(*,*) "XCHF_GAM2s written to disk"

      end if

      dimINT3=nebf*nebf*nebfBE*nebfBE*nebfBE*nebfBE*npbf*npbf
      if (LADDEXCH) dimINT3ex=dimINT3
      if(allocated(INT_GAM3)) deallocate(INT_GAM3)
      allocate(INT_GAM3(dimINT3))
      if(allocated(INT_GAM3ex1)) deallocate(INT_GAM3ex1)
      allocate(INT_GAM3ex1(dimINT3ex))
      if(allocated(INT_GAM3ex2)) deallocate(INT_GAM3ex2)
      allocate(INT_GAM3ex2(dimINT3ex))

      if (read_GAM3) then

       write(*,*)
       write(*,*) "---------------------------"

       write(*,*) " Reading:         INT_GAM3 "
       call RXCHFmult_readint(12,dimINT3,"INT_GAM3.ufm",INT_GAM3)

       if (LADDEXCH) then
        write(*,*) "               INT_GAM3ex1 "
        write(*,*) "               INT_GAM3ex2 "
        call RXCHFmult_readint(15,dimINT3ex,
     x                         "INT_GAM3ex1.ufm",INT_GAM3ex1)
        call RXCHFmult_readint(15,dimINT3ex,
     x                         "INT_GAM3ex2.ufm",INT_GAM3ex2)
       end if

       write(*,*) "---------------------------"
       write(*,*)

      else

       write(*,*)
       write(*,*) "---------------------------"
       write(*,*) " Calculating:     INT_GAM3 "
       if (LADDEXCH) then
        write(*,*) "               INT_GAM3ex1 "
        write(*,*) "               INT_GAM3ex2 "
       end if
       write(*,*) "---------------------------"
       write(*,*)

       if (LADDEXCH) then

        call RXCHFmult_GAM3_IC1ex(ng3chk,nebf,npebf,npbf,
     x                            dimINT3,ng3prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            INT_GAM3,INT_GAM3ex1,INT_GAM3ex2)
       else
        call RXCHFmult_GAM3_INT(ng3chk,
     x                          nebf,npebf,nebfBE,npebfBE,npbf,
     x                          dimINT3,nat,ngtg1,
     x                          pmass,cat,zan,bcoef1,gamma1,
     x                          KPESTR,KPEEND,
     x                          AMPEB2C,AGEBFCC,
     x                          ELCEX,ELCAM,ELCBFC,
     x                          KPESTR_be,KPEEND_be,
     x                          AMPEB2C_be,AGEBFCC_be,
     x                          ELCEX_be,ELCAM_be,ELCBFC_be,
     x                          AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                          INT_GAM3)
       end if

       open(unit=24,file="INT_GAM3.ufm",form="unformatted")
       write(24) INT_GAM3
       close(24)
       write(*,*) "INT_GAM3 written to disk"
       if (LADDEXCH) then
        open(unit=25,file="INT_GAM3ex1.ufm",form="unformatted")
        write(25) INT_GAM3ex1
        close(25)
        write(*,*) "INT_GAM3ex1 written to disk"
        open(unit=26,file="INT_GAM3ex2.ufm",form="unformatted")
        write(26) INT_GAM3ex2
        close(26)
        write(*,*) "INT_GAM3ex2 written to disk"
       end if

      end if

      if (nbe.le.2) then

       if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
       if(allocated(INT_GAM4)) deallocate(INT_GAM4)
       if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
       allocate(XCHF_GAM3(dimXCHF3))
       allocate(INT_GAM4(dimINT4))
       allocate(XCHF_GAM4(dimXCHF4))

      else

C nbe >= 3
C  - calculate four-particle integrals and store in memory
C     => XCHF_GAM3
C  - calculate five-particle integrals and store in memory
C     => INT_GAM4

       dimXCHF3=nebfBE*nebfBE*nebfBE*nebfBE*nebfBE*nebfBE*npbf*npbf
       if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
       allocate(XCHF_GAM3(dimXCHF3))

       if (read_GAM3) then

        write(*,*)
        write(*,*) "---------------------------"
        write(*,*) " Reading:        XCHF_GAM3 "
        write(*,*) "---------------------------"
        write(*,*)

        call RXCHFmult_readint(13,dimXCHF3,"XCHF_GAM3.ufm",XCHF_GAM3)

       else

        write(*,*)
        write(*,*) "---------------------------"
        write(*,*) " Calculating:    XCHF_GAM3 "
        write(*,*) "---------------------------"
        write(*,*)
      
        
        call RXCHFmult_GAM3_XCHF(ng3chk,
     x                           nebfBE,npebfBE,npbf,
     x                           dimXCHF3,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           KPESTR_be,KPEEND_be,
     x                           AMPEB2C_be,AGEBFCC_be,
     x                           ELCEX_be,ELCAM_be,ELCBFC_be,
     x                           AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                           XCHF_GAM3)

        open(unit=27,file="XCHF_GAM3.ufm",form="unformatted")
        write(27) XCHF_GAM3
        close(27)
        write(*,*) "XCHF_GAM3 written to disk"

       end if

       dimINT4=npbf*npbf*nebf*nebf*
     x         nebfBE*nebfBE*nebfBE*nebfBE*nebfBE*nebfBE
       if(allocated(INT_GAM4)) deallocate(INT_GAM4)
       allocate(INT_GAM4(dimINT4))

       if (read_GAM4) then

        write(*,*)
        write(*,*) "---------------------------"
        write(*,*) " Reading:         INT_GAM4 "
        write(*,*) "---------------------------"
        write(*,*)

        call RXCHFmult_readint(12,dimINT4,"INT_GAM4.ufm",INT_GAM4)

       else

        write(*,*)
        write(*,*) "---------------------------"
        write(*,*) " Calculating:     INT_GAM4 "
        write(*,*) "---------------------------"
        write(*,*)

        call RXCHFmult_GAM4_INT(ng4chk,nebf,nebfBE,npbf,
     x                          ngee,dimXCHF2,dimINT4,
     x                          XCHF_GAM2s,INT_GAM4)

        open(unit=28,file="INT_GAM4.ufm",form="unformatted")
        write(28) INT_GAM4
        close(28)
        write(*,*) "INT_GAM4 written to disk"

       end if

       if (nbe.le.3) then

        if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
        allocate(XCHF_GAM4(dimXCHF4))

       else

C nbe >= 4
C  - calculate five-particle integrals and store in memory
C     => XCHF_GAM4

        dimXCHF4=npbf*npbf*nebfBE*nebfBE*
     x           nebfBE*nebfBE*nebfBE*nebfBE*nebfBE*nebfBE
        if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
        allocate(XCHF_GAM4(dimXCHF4))

        if (read_GAM4) then

         write(*,*)
         write(*,*) "---------------------------"
         write(*,*) " Reading:        XCHF_GAM4 "
         write(*,*) "---------------------------"
         write(*,*)

         call RXCHFmult_readint(13,dimXCHF4,"XCHF_GAM4.ufm",XCHF_GAM4)

        else

         write(*,*)
         write(*,*) "---------------------------"
         write(*,*) " Calculating:    XCHF_GAM4 "
         write(*,*) "---------------------------"
         write(*,*)

         call RXCHFmult_GAM4_XCHF(ng4chk,nebf,nebfBE,npbf,
     x                            ngee,dimXCHF2,dimXCHF4,
     x                            XCHF_GAM2s,XCHF_GAM4)

         open(unit=29,file="XCHF_GAM4.ufm",form="unformatted")
         write(29) XCHF_GAM4
         close(29)
         write(*,*) "XCHF_GAM4 written to disk"

        end if ! read GAM4

       end if ! nbe >= 4

      end if ! nbe >= 3

      wtime1 = omp_get_wtime() - wtime

C Integral testing(
C     write(*,*)
C     write(*,*) "Writing integrals to disk"
C     write(*,*)
C
C     write(*,*) "Writing 3-particle integrals"
C     open(unit=31,file="INT_GAM2.out")
C     open(unit=32,file="XCHF_GAM2.out")
C     open(unit=33,file="XCHF_GAM2s.out")
C     do i=1,dimINT2
C       write(31,'(G20.12)') INT_GAM2(i)
C     end do
C     do i=1,dimXCHF2
C       write(32,'(G20.12)') XCHF_GAM2(i)
C       write(33,'(G20.12)') XCHF_GAM2s(i)
C     end do
C     close(33)
C     close(32)
C     close(31)
C     write(*,*) "Done"
C
C     write(*,*) "Writing 4-particle integrals"
C     open(unit=41,file="INT_GAM3.out")
C     open(unit=42,file="XCHF_GAM3.out")
C     do i=1,dimINT3
C       write(41,'(G20.12)') INT_GAM3(i)
C     end do
C     do i=1,dimXCHF3
C       write(42,'(G20.12)') XCHF_GAM3(i)
C     end do
C     close(42)
C     close(41)
C     write(*,*) "Done"
C
C     write(*,*) "Writing 5-particle integrals"
C     open(unit=51,file="INT_GAM4.out")
C     open(unit=52,file="XCHF_GAM4.out")
C     do i=1,dimINT4
C       write(51,'(G20.12)') INT_GAM4(i)
C     end do
C     do i=1,dimXCHF4
C       write(52,'(G20.12)') XCHF_GAM4(i)
C     end do
C     close(52)
C     close(51)
C     write(*,*) "Done"
C )

C Kick-off SCF
      wtime  = omp_get_wtime()

      nebf2=nebf*nebf
      nebflt=nebf*(nebf+1)/2

      nebfBE2=nebfBE*nebfBE
      nebfBElt=nebfBE*(nebfBE+1)/2

      npbf2=npbf*npbf
      npbflt=npbf*(npbf+1)/2

      if(nae.gt.1) then
       npra=(nebf-(nae/2))*(nae/2) ! occ-vir pairs for regular elecs
      else
       npra=nae*(nebf-nae)
      end if
      if (LOCBSE) then ! account for nocca less virtual orbitals
       if (nae.gt.1) then
        nprb=((nebfBE-nae/2)-(nbe/2))*(nbe/2)
       else
        nprb=((nebfBE-nae)-(nbe/2))*(nbe/2)
       end if
      else
       nprb=(nebfBE-(nbe/2))*(nbe/2) ! occ-vir pairs for special elecs
      end if

      call RXCHFmult_scf(nelec,nae,nbe,npra,nprb,nucst,
     x                   npebf,nebf,nebf2,nebflt,
     x                   npebfBE,nebfBE,nebfBE2,nebfBElt,elindBE,
     x                   npbf,npbf2,npbflt,
     x                   ngtg1,ngee,
     x                   NG2CHK,NG3CHK,NG4CHK,
     x                   read_CE,read_CP,
     x                   LG4DSCF,LG3DSCF,LG2DSCF,
     x                   LSOSCF,LOCBSE,LCMF,LALTBAS,LADDEXCH,
     x                   nat,pmass,cat,zan,
     x                   bcoef1,gamma1,
     x                   KPESTR,KPEEND,
     x                   AMPEB2C,AGEBFCC,
     x                   ELCEX,ELCAM,ELCBFC,
     x                   KPESTR_be,KPEEND_be,
     x                   AMPEB2C_be,AGEBFCC_be,
     x                   ELCEX_be,ELCAM_be,ELCBFC_be,
     x                   AGNBFCC,NUCEX,NUCAM,NUCBFC,
     x                   LG2IC,dimXCHF2,dimINT2,dimINT2ex,
     x                   XCHF_GAM2,INT_GAM2,INT_GAM2ex,XCHF_GAM2s,
     x                   LG3IC,dimXCHF3,dimINT3,dimINT3ex,
     x                   XCHF_GAM3,INT_GAM3,
     x                   INT_GAM3ex1,INT_GAM3ex2,
     x                   LG4IC,dimXCHF4,dimINT4,
     x                   XCHF_GAM4,INT_GAM4)

      wtime2 = omp_get_wtime() - wtime

C Cleanup
      if(allocated(XCHF_GAM4))   deallocate(XCHF_GAM4)
      if(allocated(INT_GAM4))    deallocate(INT_GAM4)
      if(allocated(XCHF_GAM3))   deallocate(XCHF_GAM3)
      if(allocated(INT_GAM3ex2)) deallocate(INT_GAM3ex2)
      if(allocated(INT_GAM3ex1)) deallocate(INT_GAM3ex1)
      if(allocated(INT_GAM3))    deallocate(INT_GAM3)
      if(allocated(XCHF_GAM2s))  deallocate(XCHF_GAM2s)
      if(allocated(XCHF_GAM2))   deallocate(XCHF_GAM2)
      if(allocated(INT_GAM2ex))  deallocate(INT_GAM2ex)
      if(allocated(INT_GAM2))    deallocate(INT_GAM2)

      if(allocated(ELCBFC_be))  deallocate(ELCBFC_be)
      if(allocated(ELCAM_be))   deallocate(ELCAM_be)
      if(allocated(AGEBFCC_be)) deallocate(AGEBFCC_be)
      if(allocated(ELCEX_be))   deallocate(ELCEX_be)
      if(allocated(AMPEB2C_be)) deallocate(AMPEB2C_be)

C Print timing summary
      write(*,*)
      write(*,*) "FINISHED RXCHFMULT CALCULATION"
      write(*,*)
      write(*,3000) wtime1,wtime2
      write(*,*)


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
      subroutine RXCHFmult_readint(namelen,n,fname,arr)

C Reads integrals from file "fname" into an n-dimensional array [arr]
C Integrals must be written to an unformatted file as a one-array-write
C=======================================================================
      implicit none

C Input variables
      integer namelen,n
      character(len=namelen) fname
C Output variables
      double precision arr(n)
C Local variables
      integer i

C Initialize
      arr=0.0d+00

C Read in from file
      open(unit=20,file=fname,form="unformatted")
      read(20) arr
      close(20)

      return
      end

