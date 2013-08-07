C=======================================================================
      subroutine RXCHFmult_driver(nelec,nae,nbe,nucst,
     x                            nebf,npebf,npbf,nat,ngtg1,
     x                            ng1,ng2,ng3,ng4,ngee,
     x                            ng1prm,ng2prm,ng3prm,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            ng2chk,ng3chk,ng4chk,
     x                            read_CE,read_CP,
     x                            LG2IC,LG3IC,LG4IC,
     x                            LG2DSCF,LG3DSCF,LG4DSCF,
     x                            LSOSCF,LOCBSE,LCMF)

C Driver to calculate RXCHF integrals for nbe > 1
C   XCHF_GAM* : integrals needed for XCHF contribution
C    INT_GAM* : integrals needed for interaction contribution
C
C Maximum dimension of integral needed is nbe + 1
C In this routine, both XCHF_GAM? and INT_GAM? are calculated up to
C this dimension and XCHF_GAM[max] is deallocated at the end
C  - this doesn't incur any additional cost (aside from temp memory)
C  - the case of XCHF_GAM4 (nbe > 3) is dealt with separately at the end
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
      logical read_CE,read_CP          ! Read in orbitals
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
      integer nebf2,npbf2,nebflt        ! Convenient quantities
      integer npra,nprb                 ! Number of distinct electron pairs
      integer dimXCHF2                  !
      integer dimXCHF3                  ! Dimensions of XCHF integral arrays
      integer dimXCHF4                  !
      integer dimINT2                   !
      integer dimINT3                   ! Dimensions of interaction integral arrays
      integer dimINT4                   !
      double precision              :: wtime,wtime1,wtime2  ! Timing variables
      double precision, allocatable :: XCHF_GAM2(:)         ! 3-particle XCHF integrals
      double precision, allocatable :: XCHF_GAM2s(:)        ! 3-particle XCHF overlap integrals
      double precision, allocatable :: XCHF_GAM3(:)         ! 4-particle XCHF integrals
      double precision, allocatable :: XCHF_GAM4(:)         ! 4-particle XCHF integrals
      double precision, allocatable :: INT_GAM2(:)          ! 3-particle interaction integrals
      double precision, allocatable :: INT_GAM3(:)          ! 4-particle interaction integrals
      double precision, allocatable :: INT_GAM4(:)          ! 4-particle interaction integrals

C Calculate integrals
      wtime = omp_get_wtime()

C Initialize dimensions
      dimXCHF2 = 1
      dimXCHF3 = 1
      dimXCHF4 = 1
      dimINT2  = 1
      dimINT3  = 1
      dimINT4  = 1

C nbe >= 1
C  - calculate two-particle XCHF integrals and write to disk
C     => XCHF_GAM1
C  - calculate three-particle integrals and store in memory
C     => INT_GAM2
C     => XCHF_GAM2 only needed if nbe >= 2

      write(*,*)
      write(*,*) "---------------------------"
      write(*,*) " Calculating:   XCHF_GAM1  "
      write(*,*) "---------------------------"
      write(*,*)

      call RXCHFmult_GAM1_OMP_MD(nebf,npebf,npbf,ng1,ng1prm,nat,
     x                           ngtg1,pmass,cat,zan,bcoef1,gamma1,
     x                           AMPEB2C,AGEBFCC,AGNBFCC,ELCEX,
     x                           NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      write(*,*)
      write(*,*) "---------------------------"
      write(*,*) " Calculating:     INT_GAM2 "
      write(*,*) "                 XCHF_GAM2 "
      write(*,*) "---------------------------"
      write(*,*)

      dimINT2=ng2
      dimXCHF2=ng2
      if(allocated(INT_GAM2)) deallocate(INT_GAM2)
      allocate(INT_GAM2(dimINT2))
      if(allocated(XCHF_GAM2)) deallocate(XCHF_GAM2)
      allocate(XCHF_GAM2(dimXCHF2))
      if(allocated(XCHF_GAM2s)) deallocate(XCHF_GAM2s)
      allocate(XCHF_GAM2s(dimXCHF2))

      call RXCHFmult_GAM2_IC1(ng2chk,nebf,npebf,npbf,
     x                        ng2,ng2prm,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        XCHF_GAM2,INT_GAM2,XCHF_GAM2s)

      if (nbe.le.1) then

        write(*,*)
        write(*,*) "NBE = ",NBE," <= 1"
        write(*,*) "so deallocating XCHF_GAM2"
        write(*,*)

        dimXCHF2=1
        if(allocated(XCHF_GAM2s)) deallocate(XCHF_GAM2s)
        if(allocated(XCHF_GAM2)) deallocate(XCHF_GAM2)
        allocate(XCHF_GAM2(dimXCHF2))
        allocate(XCHF_GAM2s(dimXCHF2))

        if(allocated(INT_GAM3)) deallocate(INT_GAM3)
        if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
        if(allocated(INT_GAM4)) deallocate(INT_GAM4)
        if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
        allocate(INT_GAM3(dimINT3))
        allocate(XCHF_GAM3(dimXCHF3))
        allocate(INT_GAM4(dimINT4))
        allocate(XCHF_GAM4(dimXCHF4))

      else

C nbe >= 2
C  - calculate four-particle integrals and store in memory
C     => INT_GAM3
C     => XCHF_GAM3 only needed if nbe >= 3

        write(*,*) "---------------------------"
        write(*,*) " Calculating:     INT_GAM3 "
        write(*,*) "                 XCHF_GAM3 "
        write(*,*) "---------------------------"
        write(*,*)

        dimINT3=ng3
        dimXCHF3=ng3
        if(allocated(INT_GAM3)) deallocate(INT_GAM3)
        allocate(INT_GAM3(dimINT3))
        if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
        allocate(XCHF_GAM3(dimXCHF3))

        call RXCHFmult_GAM3_IC1(ng3chk,nebf,npebf,npbf,
     x                          ng3,ng3prm,nat,ngtg1,
     x                          pmass,cat,zan,bcoef1,gamma1,
     x                          KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                          ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                          XCHF_GAM3,INT_GAM3)

        if (nbe.le.2) then

          write(*,*)
          write(*,*) "NBE = ",NBE," <= 2"
          write(*,*) "so deallocating XCHF_GAM3"
          write(*,*)

          dimXCHF3=1
          if(allocated(XCHF_GAM3)) deallocate(XCHF_GAM3)
          allocate(XCHF_GAM3(dimXCHF3))

          if(allocated(INT_GAM4)) deallocate(INT_GAM4)
          if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
          allocate(INT_GAM4(dimINT4))
          allocate(XCHF_GAM4(dimXCHF4))

        else

C nbe >= 3
C  - calculate interaction five-particle integrals and store in memory
C     => INT_GAM4

          write(*,*)
          write(*,*) "---------------------------"
          write(*,*) " Calculating:     INT_GAM4 "
          write(*,*) "---------------------------"
          write(*,*)

          dimINT4=ng4
          if(allocated(INT_GAM4)) deallocate(INT_GAM4)
          allocate(INT_GAM4(dimINT4))

          call RXCHFmult_GAM4_ICR(ng4chk,nebf,npbf,ngee,ng2,ng4,
     x                            XCHF_GAM2s,INT_GAM4)

          if (nbe.le.3) then

            dimXCHF4=1
            if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
            allocate(XCHF_GAM4(dimXCHF4))

          else

C nbe >= 4
C  - calculate XCHF five-particle integrals and store in memory
C     => XCHF_GAM4

            write(*,*)
            write(*,*) "---------------------------"
            write(*,*) " Calculating:    XCHF_GAM4 "
            write(*,*) "---------------------------"
            write(*,*)

            dimXCHF4=ng4
            if(allocated(XCHF_GAM4)) deallocate(XCHF_GAM4)
            allocate(XCHF_GAM4(dimXCHF4))

            call GAM4_ICR(ng4chk,nebf,npbf,ngee,ng2,ng4,
     x                    XCHF_GAM2s,XCHF_GAM4)

          end if ! nbe >= 4

        end if ! nbe >= 3

      end if ! nbe >= 2

      wtime1 = omp_get_wtime() - wtime

C      write(*,*) "dimINT2:",dimINT2
C      write(*,*) "dimINT3:",dimINT3
C      write(*,*) "dimINT4:",dimINT4
C      write(*,*) "dimXCHF2:",dimXCHF2
C      write(*,*) "dimXCHF3:",dimXCHF3
C      write(*,*) "dimXCHF4:",dimXCHF4

C ARS(
      write(*,*) "Writing INT_GAM2 to disk..."
      open(unit=20,file="INT_GAM2.ufm",form="unformatted")
      write(20) INT_GAM2
      close(20)
      write(*,*) "Done."
      write(*,*)

      write(*,*) "Writing XCHF_GAM2 to disk..."
      open(unit=21,file="XCHF_GAM2.ufm",form="unformatted")
      write(21) XCHF_GAM2
      close(21)
      write(*,*) "Done."

      write(*,*) "Writing XCHF_GAM2s to disk..."
      open(unit=22,file="XCHF_GAM2s.ufm",form="unformatted")
      write(22) XCHF_GAM2s
      close(22)
      write(*,*) "Done."

      write(*,*) "Writing INT_GAM3 to disk..."
      open(unit=23,file="INT_GAM3.ufm",form="unformatted")
      write(23) INT_GAM3
      close(23)
      write(*,*) "Done."
      write(*,*)
C )

C Kick-off SCF
      wtime  = omp_get_wtime()

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
     x                   LSOSCF,LOCBSE,LCMF,
     x                   ng2prm,ng3prm,nat,pmass,cat,zan,
     x                   bcoef1,gamma1,
     x                   KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                   ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                   LG2IC,dimXCHF2,dimINT2,
     x                   XCHF_GAM2,INT_GAM2,XCHF_GAM2s,
     x                   LG3IC,dimXCHF3,dimINT3,
     x                   XCHF_GAM3,INT_GAM3,
     x                   LG4IC,dimXCHF4,dimINT4,
     x                   XCHF_GAM4,INT_GAM4)

      wtime2 = omp_get_wtime() - wtime

C Cleanup
      if(allocated(XCHF_GAM4))  deallocate(XCHF_GAM4)
      if(allocated(INT_GAM4))   deallocate(INT_GAM4)
      if(allocated(XCHF_GAM3))  deallocate(XCHF_GAM3)
      if(allocated(INT_GAM3))   deallocate(INT_GAM3)
      if(allocated(XCHF_GAM2s)) deallocate(XCHF_GAM2s)
      if(allocated(XCHF_GAM2))  deallocate(XCHF_GAM2)
      if(allocated(INT_GAM2))   deallocate(INT_GAM2)

C Print timing summary
      write(*,*)
      write(*,*) "FINISHED RXCHFMULT CALCULATION"
      write(*,*)
      write(*,3000) wtime1,wtime2
      write(*,*)


 3000 FORMAT(/8X,'  +--------------------------------------+',/,
     X        8X,'  |    TIMING SUMMARY FOR CALCULATION    |',/,
     x        8X,'  +--------------------------------------+',/,
     x        8X,'    TIME TO EVALUATE INTEGRALS:',1X,F12.4/
     x        8X,'                  TIME FOR SCF:',1X,F12.4/)

      return
      end

