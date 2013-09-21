C=======================================================================
      subroutine RXCHF_driver_MPI(nproc,rank,
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
     x                            LSOSCF,LOCBSE,LCMF,LADDEXCH)

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
      integer ierr
      integer nebf2,npbf2,nebflt        ! Convenient quantities
      integer npra,nprb                 ! Number of distinct electron pairs
      integer dimXCHF2                  !
      integer dimXCHF3                  ! Dimensions of XCHF integral arrays
      integer dimXCHF4                  !
      integer dimINT2                   !
      integer dimINT3                   ! Dimensions of interaction integral arrays
      integer dimINT4                   !
      integer dimINT2ex                 ! Dimensions of exchange integral arrays
      integer dimINT3ex                 !
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

      dimINT2=ng2
      dimXCHF2=ng2
      if (LADDEXCH) dimINT2ex=ng2

      if(allocated(INT_GAM2)) deallocate(INT_GAM2)
      allocate(INT_GAM2(dimINT2))
      if(allocated(INT_GAM2ex)) deallocate(INT_GAM2ex)
      allocate(INT_GAM2ex(dimINT2ex))
      if(allocated(XCHF_GAM2)) deallocate(XCHF_GAM2)
      allocate(XCHF_GAM2(dimXCHF2))
      if(allocated(XCHF_GAM2s)) deallocate(XCHF_GAM2s)
      allocate(XCHF_GAM2s(dimXCHF2))

      if (read_GAM2) then

! Read with master process
       if (rank.eq.0) then
        write(*,*)
        write(*,*) "---------------------------"
 
        write(*,*) " Reading:         INT_GAM2 "
        call RXCHFmult_readint(ng2,"INT_GAM2.ufm",INT_GAM2)

        if (LADDEXCH) then
         write(*,*) "                INT_GAM2ex "
         call RXCHFmult_readint(ng2,"INT_GAM2ex.ufm",INT_GAM2ex)
        end if

        if (nbe.gt.1) then
         write(*,*) "                 XCHF_GAM2 "
         write(*,*) "                XCHF_GAM2s "
         call RXCHFmult_readint(ng2,"XCHF_GAM2.ufm",XCHF_GAM2)
         call RXCHFmult_readint(ng2,"XCHF_GAM2s.ufm",XCHF_GAM2s)
        end if

        write(*,*) "---------------------------"
        write(*,*)
       end if

! Broadcast to all processes
       call MPI_BCAST(INT_GAM2,dimINT2,MPI_DOUBLE_PRECISION,
     x                0,MPI_COMM_WORLD,ierr)
       if (LADDEXCH) then
        call MPI_BCAST(INT_GAM2ex,dimINT2ex,MPI_DOUBLE_PRECISION,
     x                 0,MPI_COMM_WORLD,ierr)
       end if
       if (nbe.gt.1) then
        call MPI_BCAST(XCHF_GAM2,dimXCHF2,MPI_DOUBLE_PRECISION,
     x                 0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(XCHF_GAM2s,dimXCHF2,MPI_DOUBLE_PRECISION,
     x                 0,MPI_COMM_WORLD,ierr)
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
     x                        ng2,ng2prm,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        XCHF_GAM2,INT_GAM2,INT_GAM2ex,XCHF_GAM2s)
       else
        call RXCHF_GAM2_MPI(nproc,rank,
     x                      ng2chk,nebf,npebf,npbf,
     x                      ng2,ng2prm,nat,ngtg1,
     x                      pmass,cat,zan,bcoef1,gamma1,
     x                      KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                      XCHF_GAM2,INT_GAM2,XCHF_GAM2s)
       end if

! Write integrals to disk with master process
       if (rank.eq.0) then
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

! Write integrals to disk with master process
        if (rank.eq.0) then
         if (.not.(read_GAM2)) then
          open(unit=22,file="XCHF_GAM2.ufm",form="unformatted")
          write(22) XCHF_GAM2
          close(22)
          write(*,*) "XCHF_GAM2 written to disk"
          open(unit=23,file="XCHF_GAM2s.ufm",form="unformatted")
          write(23) XCHF_GAM2s
          close(23)
          write(*,*) "XCHF_GAM2s written to disk"
         end if
        end if

!!!!!!! debug
      if(allocated(XCHF_GAM2s))  deallocate(XCHF_GAM2s)
      if(allocated(XCHF_GAM2))   deallocate(XCHF_GAM2)
      if(allocated(INT_GAM2ex))  deallocate(INT_GAM2ex)
      if(allocated(INT_GAM2))    deallocate(INT_GAM2)
      return
!!!!!!!

C nbe >= 2
C  - calculate four-particle integrals and store in memory
C     => INT_GAM3
C     => INT_GAM3ex only needed if RXCHF-ae
C     => XCHF_GAM3 only needed if nbe >= 3

        dimINT3=ng3
        dimXCHF3=ng3
        if (LADDEXCH) dimINT3ex=ng3

        if(allocated(INT_GAM3)) deallocate(INT_GAM3)
        allocate(INT_GAM3(dimINT3))
        if(allocated(INT_GAM3ex1)) deallocate(INT_GAM3ex1)
        allocate(INT_GAM3ex1(dimINT3ex))
        if(allocated(INT_GAM3ex2)) deallocate(INT_GAM3ex2)
        allocate(INT_GAM3ex2(dimINT3ex))
        if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
        allocate(XCHF_GAM3(dimXCHF3))

        if (read_GAM3) then

! Read with master process
         if (rank.eq.0) then
          write(*,*)
          write(*,*) "---------------------------"

          write(*,*) " Reading:         INT_GAM3 "
          call RXCHFmult_readint(ng3,"INT_GAM3.ufm",INT_GAM3)

          if (LADDEXCH) then
           write(*,*) "               INT_GAM3ex1 "
           write(*,*) "               INT_GAM3ex2 "
           call RXCHFmult_readint(ng3,"INT_GAM3ex1.ufm",INT_GAM3ex1)
           call RXCHFmult_readint(ng3,"INT_GAM3ex2.ufm",INT_GAM3ex2)
          end if

          if (nbe.gt.2) then
           write(*,*) "                 XCHF_GAM3 "
           call RXCHFmult_readint(ng3,"XCHF_GAM3.ufm",XCHF_GAM3)
          end if

          write(*,*) "---------------------------"
          write(*,*)
         end if

! Broadcast to all processes
         call MPI_BCAST(INT_GAM3,dimINT3,MPI_DOUBLE_PRECISION,
     x                  0,MPI_COMM_WORLD,ierr)
         if (LADDEXCH) then
          call MPI_BCAST(INT_GAM3ex1,dimINT3ex,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(INT_GAM3ex2,dimINT3ex,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
         end if
         if (nbe.gt.2) then
          call MPI_BCAST(XCHF_GAM3,dimXCHF3,MPI_DOUBLE_PRECISION,
     x                   0,MPI_COMM_WORLD,ierr)
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

         if (LADDEXCH) then
          call RXCHFmult_GAM3_IC1ex(ng3chk,nebf,npebf,npbf,
     x                            ng3,ng3prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            XCHF_GAM3,INT_GAM3,
     x                            INT_GAM3ex1,INT_GAM3ex2)
         else
          call RXCHFmult_GAM3_IC1(ng3chk,nebf,npebf,npbf,
     x                            ng3,ng3prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            XCHF_GAM3,INT_GAM3)
         end if

! Write integrals to disk with master process
         if (rank.eq.0) then
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

! Write integrals to disk with master process
          if (rank.eq.0) then
           if (.not.(read_GAM3)) then
            open(unit=27,file="XCHF_GAM3.ufm",form="unformatted")
            write(27) XCHF_GAM3
            close(27)
            write(*,*) "XCHF_GAM3 written to disk"
           end if
          end if

C nbe >= 3
C  - calculate interaction five-particle integrals and store in memory
C     => INT_GAM4

          dimINT4=ng4
          if(allocated(INT_GAM4)) deallocate(INT_GAM4)
          allocate(INT_GAM4(dimINT4))

          if (read_GAM4) then

! Read with master process
           if (rank.eq.0) then
            write(*,*)
            write(*,*) "---------------------------"
            write(*,*) " Reading:         INT_GAM4 "
            write(*,*) "---------------------------"
            write(*,*)

            call RXCHFmult_readint(ng4,"INT_GAM4.ufm",INT_GAM4)
           end if

! Broadcast to all processes
           call MPI_BCAST(INT_GAM4,dimINT4,MPI_DOUBLE_PRECISION,
     x                    0,MPI_COMM_WORLD,ierr)

          else

           if (rank.eq.0) then
            write(*,*)
            write(*,*) "---------------------------"
            write(*,*) " Calculating:     INT_GAM4 "
            write(*,*) "---------------------------"
            write(*,*)
           end if

           call RXCHFmult_GAM4_ICR(ng4chk,nebf,npbf,ngee,ng2,ng4,
     x                             XCHF_GAM2s,INT_GAM4)

! Write integrals to disk with master process
           if (rank.eq.0) then
            open(unit=28,file="INT_GAM4.ufm",form="unformatted")
            write(28) INT_GAM4
            close(28)
            write(*,*) "INT_GAM4 written to disk"
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

            dimXCHF4=ng4
            if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
            allocate(XCHF_GAM4(dimXCHF4))

            if (read_GAM4) then

! Read with master process
             if (rank.eq.0) then
              write(*,*)
              write(*,*) "---------------------------"
              write(*,*) " Reading:        XCHF_GAM4 "
              write(*,*) "---------------------------"
              write(*,*)

              call RXCHFmult_readint(ng4,"XCHF_GAM4.ufm",XCHF_GAM4)
             end if

! Broadcast to all processes
             call MPI_BCAST(XCHF_GAM4,dimXCHF4,MPI_DOUBLE_PRECISION,
     x                      0,MPI_COMM_WORLD,ierr)

            else

             if (rank.eq.0) then
              write(*,*)
              write(*,*) "---------------------------"
              write(*,*) " Calculating:    XCHF_GAM4 "
              write(*,*) "---------------------------"
              write(*,*)
             end if

             call GAM4_ICR(ng4chk,nebf,npbf,ngee,ng2,ng4,
     x                     XCHF_GAM2s,XCHF_GAM4)

! Write integrals to disk with master process
             if (rank.eq.0) then
              open(unit=29,file="XCHF_GAM4.ufm",form="unformatted")
              write(29) XCHF_GAM4
              close(29)
              write(*,*) "XCHF_GAM4 written to disk"
             end if

            end if ! read GAM4

          end if ! nbe >= 4

        end if ! nbe >= 3

      end if ! nbe >= 2

      wtime1 = MPI_WTIME() - wtime

C Kick-off SCF
      wtime  = MPI_WTIME()

      nebf2=nebf*nebf
      npbf2=npbf*npbf
      if(nae.gt.1) then
       npra=(nebf-(nae/2))*(nae/2) ! occ-vir pairs for regular elecs
      else
         npra=nae*(nebf-nae)
      end if
      if(nbe.gt.1) then
       if (LOCBSE) then ! account for nocca less virtual orbitals
        if (nae.gt.1) then
         nprb=((nebf-nae/2)-(nbe/2))*(nbe/2)
        else
         nprb=((nebf-nae)-(nbe/2))*(nbe/2)
        end if
       else
        nprb=(nebf-(nbe/2))*(nbe/2) ! occ-vir pairs for special elecs
       end if
      else
       nprb=nbe*(nebf-nbe)
      end if
      nebflt=nebf*(nebf+1)/2

      call RXCHFmult_scf(nelec,nae,nbe,npra,nprb,nebflt,nucst,
     x                   npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                   ngtg1,ng1,ng2,ng3,ng4,
     x                   NG2CHK,NG3CHK,NG4CHK,
     x                   read_CE,read_CP,
     x                   LG4DSCF,LG3DSCF,LG2DSCF,
     x                   LSOSCF,LOCBSE,LCMF,LADDEXCH,
     x                   ng2prm,ng3prm,nat,pmass,cat,zan,
     x                   bcoef1,gamma1,
     x                   KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                   ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                   LG2IC,dimXCHF2,dimINT2,dimINT2ex,
     x                   XCHF_GAM2,INT_GAM2,INT_GAM2ex,XCHF_GAM2s,
     x                   LG3IC,dimXCHF3,dimINT3,dimINT3ex,
     x                   XCHF_GAM3,INT_GAM3,
     x                   INT_GAM3ex1,INT_GAM3ex2,
     x                   LG4IC,dimXCHF4,dimINT4,
     x                   XCHF_GAM4,INT_GAM4)

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

