C=======================================================================
      subroutine RXCHF_driver_MPI(nproc,rank,
     x                            nblocks,blockrank,
     x                            nelec,nae,nbe,nucst,
     x                            nebf,npebf,npbf,nat,ngtg1,
     x                            ng1,ng2,ng3,ng4,ngee,
     x                            ng1prm,ng2prm,ng3prm,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            ng2chk,ng3chk,ng4chk,
     x                            read_CE,read_CP,
     x                            read_GAM2,read_GAM3,read_GAM4,
     x                            LG2IC,LG3IC,LG4IC,
     x                            LG2DSCF,LG3DSCF,LG4DSCF,
     x                            LSOSCF,LDIIS,LSORXCHF,OCBSETYP,
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
      integer nbe                      ! Number of special electrons
      integer nucst                    ! Nuclear state for which to form density
      integer nebf                     ! Number of contracted electronic basis functions
      integer npebf                    ! Number of primitive electronic basis functions
      integer npbf                     ! Number of proton basis functions
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
      logical LG2IC                    ! In-core storage of 3-particle integrals
      logical LG3IC                    ! In-core storage of 4-particle integrals
      logical LG4IC                    ! In-core storage of 5-particle integrals
      logical LG2DSCF                  ! Direct SCF for 3-particle integrals
      logical LG3DSCF                  ! Direct SCF for 4-particle integrals
      logical LG4DSCF                  ! Direct SCF for 5-particle integrals
      logical LSOSCF                   ! SOSCF where applicable
      logical LDIIS                    ! DIIS where applicable
      logical LSORXCHF                 ! SORXCHF where applicable
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
#if MPI32
      integer*4    ierr
#else
      integer      ierr
#endif
      character*4  istring              ! File I/O variables
      integer      nebf2,npbf2          ! Convenient quantities
      integer      nebflt,npbflt        !
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

C nbe >= 1
C  - calculate two-particle XCHF integrals and write to disk
C     => XCHF_GAM1
C  - calculate three-particle integrals and store in memory
C     => INT_GAM2
C     => INT_GAM2ex only needed if RXCHF-ae
C     => XCHF_GAM2 only needed if nbe >= 2

! Calculate two-particle integrals and write to disk with master process
      if (rank.eq.0) then
       write(*,*)
       write(*,*) "---------------------------"
       write(*,*) " Calculating:   XCHF_GAM1  "
       write(*,*) "---------------------------"
       write(*,*)

       call RXCHFmult_GAM1_OMP_MD(nebf,npebf,npbf,ng1,ng1prm,nat,
     x                            ngtg1,pmass,cat,zan,bcoef1,gamma1,
     x                            AMPEB2C,AGEBFCC,AGNBFCC,ELCEX,
     x                            NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Have each process calculate and store its ng2/proc integrals
! Assign residual integrals to last process
      if (rank.eq.(nproc-1)) then
       dimINT2=ng2/nproc+mod(ng2,nproc)
       dimXCHF2=ng2/nproc+mod(ng2,nproc)
       if (LADDEXCH) dimINT2ex=ng2/nproc+mod(ng2,nproc)
      else
       dimINT2=ng2/nproc
       dimXCHF2=ng2/nproc
       if (LADDEXCH) dimINT2ex=ng2/nproc
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
     x                        ng2,dimINT2,ng2prm,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        XCHF_GAM2,INT_GAM2,INT_GAM2ex,XCHF_GAM2s)
       else
        call RXCHF_GAM2_MPI(nproc,rank,
     x                      ng2chk,nebf,npebf,npbf,
     x                      ng2,dimINT2,ng2prm,nat,ngtg1,
     x                      pmass,cat,zan,bcoef1,gamma1,
     x                      KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
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
        if(allocated(INT_GAM4)) deallocate(INT_GAM4)
        if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
        allocate(INT_GAM3(dimINT3))
        allocate(INT_GAM3ex1(dimINT3ex))
        allocate(INT_GAM3ex2(dimINT3ex))
        allocate(XCHF_GAM3(dimXCHF3))
        allocate(INT_GAM4(dimINT4))
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

! Have each process calculate and store its ng3/proc integrals
! Assign residual integrals to last process
        if (rank.eq.(nproc-1)) then
         dimINT3=ng3/nproc+mod(ng3,nproc)
         dimXCHF3=ng3/nproc+mod(ng3,nproc)
         if (LADDEXCH) dimINT3ex=ng3/nproc+mod(ng3,nproc)
        else
         dimINT3=ng3/nproc
         dimXCHF3=ng3/nproc
         if (LADDEXCH) dimINT3ex=ng3/nproc
        end if

! These arrays now local to each process and correspond to only its chunk of integrals
        if(allocated(INT_GAM3)) deallocate(INT_GAM3)
        allocate(INT_GAM3(dimINT3))
        if(allocated(INT_GAM3ex1)) deallocate(INT_GAM3ex1)
        allocate(INT_GAM3ex1(dimINT3ex))
        if(allocated(INT_GAM3ex2)) deallocate(INT_GAM3ex2)
        allocate(INT_GAM3ex2(dimINT3ex))
        if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
        allocate(XCHF_GAM3(dimXCHF3))

        if (read_GAM3) then

! Read chunk with process (restart job requires same number of processes)
         if (rank.eq.0) then
          write(*,*)
          write(*,*) "---------------------------"
          write(*,*) " Reading:         INT_GAM3 "
         end if
         call RXCHF_readint_mpi(nproc,rank,dimINT3,17,
     x                          "INT_GAM3-"//istring//".ufm",
     x                          INT_GAM3)

         if (LADDEXCH) then
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
          if(rank.eq.0) write(*,*) "                 XCHF_GAM3 "
          call RXCHF_readint_mpi(nproc,rank,dimXCHF3,18,
     x                           "XCHF_GAM3-"//istring//".ufm",
     x                           XCHF_GAM3)
         end if

         if (rank.eq.0) then
          write(*,*) "---------------------------"
          write(*,*)
         end if

        else

         if (rank.eq.0) then
          write(*,*)
          write(*,*) "---------------------------"
          write(*,*) " Calculating:     INT_GAM3 "
          if (LADDEXCH) then
           write(*,*) "               INT_GAM3ex1 "
           write(*,*) "               INT_GAM3ex2 "
          end if
          write(*,*) "                 XCHF_GAM3 "
          write(*,*) "---------------------------"
          write(*,*)
         end if

        call get_mpi_range(ng3,nblocks,blockrank,blockstart,blockend)
        if (blockrank.eq.(nblocks-1)) then
         ng3block=ng3/nblocks+mod(ng3,nblocks)
         blockend=ng3
        else
         ng3block=ng3/nblocks
        end if
        if(rank.eq.0) then
         write(*,*) "Computing block ",blockrank," of ",nblocks,
     x              " total blocks"
         write(*,*) "ng3,blockstart,blockend:",ng3,blockstart,blockend
        end if
        if (rank.eq.(nproc-1)) then
         dimINT3=ng3block/nproc+mod(ng3block,nproc)
         dimXCHF3=ng3block/nproc+mod(ng3block,nproc)
         if (LADDEXCH) dimINT3ex=ng3block/nproc+mod(ng3block,nproc)
        else
         dimINT3=ng3block/nproc
         dimXCHF3=ng3block/nproc
         if (LADDEXCH) dimINT3ex=ng3block/nproc
        end if

! Only symmetrize XCHF_GAM3 if nbe > 2 (saves significant MPI comm)
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
          call RXCHF_GAM3_MPI(nproc,rank,
     x                        nblocks,blockrank,
     x                        blockstart,blockend,
     x                        ng3chk,nebf,npebf,npbf,
     x                        ng3block,dimINT3,ng3prm,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        LXCHFSYMM,XCHF_GAM3,INT_GAM3)
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

          if (rank.eq.0) then
           write(*,*)
           write(*,*) "NBE = ",NBE," <= 2"
           write(*,*) "so deallocating XCHF_GAM3"
           write(*,*)
          end if

          dimXCHF3=1
          if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
          allocate(XCHF_GAM3(dimXCHF3))

          if(allocated(INT_GAM4)) deallocate(INT_GAM4)
          if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
          allocate(INT_GAM4(dimINT4))
          allocate(XCHF_GAM4(dimXCHF4))

        else

! Write integrals to disk (chunk calculated by each process)
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

! Have each process calculate and store its ng4/proc integrals
! Assign residual integrals to last process
            if (rank.eq.(nproc-1)) then
             dimINT4=ng4/nproc+mod(ng4,nproc)
            else
             dimINT4=ng4/nproc
            end if

! This array now local to each process and correspond to only its chunk of integrals
            if(allocated(INT_GAM4)) deallocate(INT_GAM4)
            allocate(INT_GAM4(dimINT4))

            if (read_GAM4) then

! Read chunk with process (restart job requires same number of processes)
             if (rank.eq.0) then
              write(*,*)
              write(*,*) "---------------------------"
              write(*,*) " Reading:         INT_GAM4 "
              write(*,*) "---------------------------"
              write(*,*)
             end if

             call RXCHF_readint_mpi(nproc,rank,dimINT4,17,
     x                              "INT_GAM4-"//istring//".ufm",
     x                              INT_GAM4)

            else

             if(nbe.le.3) then     ! Otherwise calc INT_GAM4 and XCHF_GAM4 together

              if (rank.eq.0) then
               write(*,*)
               write(*,*) "---------------------------"
               write(*,*) " Calculating:     INT_GAM4 "
               write(*,*) "---------------------------"
               write(*,*)
              end if

              call RXCHF_GAM4_MPI(nproc,rank,
     x                            ng4chk,nebf,npbf,ngee,ng2,ng4,
     x                            dimXCHF2,dimINT4,
     x                            XCHF_GAM2s,INT_GAM4)

C! Write integrals to disk (chunk calculated by each process)
C            call RXCHF_writeint_mpi(nproc,rank,dimINT4,17,
C     x                              "INT_GAM4-"//istring//".ufm",
C     x                              INT_GAM4)
C            if(rank.eq.0) write(*,*) "INT_GAM4 written to disk"

             end if

            end if

            if (nbe.le.3) then

              dimXCHF4=1
              if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
              allocate(XCHF_GAM4(dimXCHF4))

            else

C nbe >= 4
C  - calculate XCHF five-particle integrals and store in memory
C     => XCHF_GAM4
C  - assume there is enough memory for each process to store GAM2s

! Have each process calculate and store its ng4/proc integrals
! Assign residual integrals to last process
              if (rank.eq.(nproc-1)) then
               dimXCHF4=ng4/nproc+mod(ng4,nproc)
              else
               dimXCHF4=ng4/nproc
              end if

! This array now local to each process and correspond to only its chunk of integrals
              if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
              allocate(XCHF_GAM4(dimXCHF4))

              if (read_GAM4) then

! Read chunk with process (restart job requires same number of processes)
               if (rank.eq.0) then
                write(*,*)
                write(*,*) "---------------------------"
                write(*,*) " Reading:        XCHF_GAM4 "
                write(*,*) "---------------------------"
                write(*,*)
               end if

               call RXCHF_readint_mpi(nproc,rank,dimXCHF4,18,
     x                                "XCHF_GAM4-"//istring//".ufm",
     x                                XCHF_GAM4)

              else

               if (rank.eq.0) then
                write(*,*)
                write(*,*) "---------------------------"
                write(*,*) " Calculating:     INT_GAM4 "
                write(*,*) " Calculating:    XCHF_GAM4 "
                write(*,*) "---------------------------"
                write(*,*)
               end if

               call XCHF_GAM4_MPI(nproc,rank,
     x                            ng4chk,nebf,npbf,ngee,ng2,ng4,
     x                            dimXCHF2,dimXCHF4,
     x                            XCHF_GAM2s,INT_GAM4,XCHF_GAM4)

! Write integrals to disk (chunk calculated by each process)
             call RXCHF_writeint_mpi(nproc,rank,dimINT4,17,
     x                               "INT_GAM4-"//istring//".ufm",
     x                               INT_GAM4)
             if(rank.eq.0) write(*,*) "INT_GAM4 written to disk"
             call RXCHF_writeint_mpi(nproc,rank,dimXCHF4,18,
     x                               "XCHF_GAM4-"//istring//".ufm",
     x                               XCHF_GAM4)
             if(rank.eq.0) write(*,*) "XCHF_GAM4 written to disk"

              end if ! read GAM4

            end if ! nbe >= 4

          end if ! nblocks == 1

        end if ! nbe >= 3

      end if ! nbe >= 2

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      wtime1 = MPI_WTIME() - wtime

C Kick-off SCF
      wtime  = MPI_WTIME()

      nebf2=nebf*nebf
      npbf2=npbf*npbf
      nebflt=nebf*(nebf+1)/2
      npbflt=npbf*(npbf+1)/2

      if ((read_GAM3).or.(nblocks.eq.1)) then
       call RXCHF_scf_MPI(nproc,rank,
     x                    nelec,nae,nbe,nebflt,npbflt,nucst,
     x                    npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                    ngtg1,ng1,ng2,ng3,ng4,
     x                    NG2CHK,NG3CHK,NG4CHK,
     x                    read_CE,read_CP,
     x                    LG4DSCF,LG3DSCF,LG2DSCF,
     x                    LSOSCF,LDIIS,LSORXCHF,
     x                    OCBSETYP,LDBG,LADDEXCH,
     x                    ng2prm,ng3prm,nat,pmass,cat,zan,
     x                    bcoef1,gamma1,
     x                    KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                    ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                    LG2IC,dimXCHF2,dimINT2,dimINT2ex,
     x                    XCHF_GAM2,INT_GAM2,INT_GAM2ex,XCHF_GAM2s,
     x                    LG3IC,dimXCHF3,dimINT3,dimINT3ex,
     x                    XCHF_GAM3,INT_GAM3,
     x                    INT_GAM3ex1,INT_GAM3ex2,
     x                    LG4IC,dimXCHF4,dimINT4,
     x                    XCHF_GAM4,INT_GAM4)
      end if

      wtime2 = MPI_WTIME() - wtime

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

C Print timing summary
      if (rank.eq.0) then
       write(*,*)
       write(*,*) "FINISHED RXCHFMULT CALCULATION"
       write(*,*)
       write(*,3000) wtime1,wtime2
       write(*,*)
      end if


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

