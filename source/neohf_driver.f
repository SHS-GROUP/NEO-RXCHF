C=======================================================================
      subroutine NEOHF_driver(nelec,nae,nbe,nucst,
     x                        nebf,npebf,npbf,nat,ngtg,ngee,
     x                        pmass,cat,zan,bgem,ggem,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        read_CE,read_CP,
     x                        LXCUHF,LSOSCF,LDBG)

C Driver to interface to serial NEO-HF code
C  - should be easily extendable to doing XCHF calculations with the serial code
C
C=======================================================================
      implicit none
      include 'mpif.h'

C Input variables
      integer nelec                    ! Total number of electrons
      integer nae                      ! Number of alpha electrons
      integer nbe                      ! Number of beta electrons
      integer nucst                    ! Nuclear state for which to form density
      integer nebf                     ! Number of contracted electronic basis functions
      integer npebf                    ! Number of primitive electronic basis functions
      integer npbf                     ! Number of proton basis functions
      integer nat                      ! Number of classical nuclei
      integer ngtg                     ! Number of geminal functions (dummy)
      integer ngee                     ! Numbers of 2-elec integrals (contracted)
      integer ELCAM(npebf,3)           ! Angular mom for electrons
      integer NUCAM(npbf,3)            ! Angular mom for quantum nuclei
      integer AMPEB2C(npebf)           ! Map primitive index to contracted
      integer KPESTR(nebf)             ! Map contracted index to primitive start
      integer KPEEND(nebf)             ! Map contracted index to primitive end
      logical LXCUHF                   ! Flag for NEO-UHF (separate SCF routine)
      logical LSOSCF                   ! Flag for SOSCF
      logical LDBG                     ! Debugging flag
      logical read_CE,read_CP          ! Read in orbitals
      double precision pmass           ! Mass of nonelectron quantum particle 
      double precision zan(nat)        ! Classical nuclear charges
      double precision cat(3,nat)      ! Classical nuclear coordinates (au)
      double precision bgem(ngtg)      ! Geminal b_k (dummy)
      double precision ggem(ngtg)      ! Geminal gamma_k (dummy)
      double precision ELCEX(npebf)    ! Electronic basis function exponents
      double precision NUCEX(npbf)     ! Nuclear basis function exponents
      double precision ELCBFC(npebf,3) ! Electronic basis function centers
      double precision NUCBFC(npbf,3)  ! Nuclear basis function centers
      double precision AGEBFCC(npebf)  ! Map primitive index to contract coeff
      double precision AGNBFCC(npbf)   ! Nuclear contraction coeff

C Local variables
      integer          nebf2,npbf2          ! Convenient quantities
      integer          nebflt,npbflt        !
      integer          npr,npra,nprb        ! Number of distinct electron pairs
      integer          ng1                  ! Number of ep integrals
      integer          intdum               !
      logical          logdum               ! Dummy variables
      double precision arrdum(1)            !
      double precision wtime,wtime1,wtime2  ! Timing variables


C Calculate integrals and store on disk
      wtime = MPI_WTIME()

      call class_nuc_rep(nat,zan,cat)

      call elec_ovlap(npebf,nebf,nebf*nebf,
     x                AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

      call check_elec_ovlap(nebf)

      call nuc_ovlap(npbf,npbf*npbf,
     x                AGNBFCC,NUCEX,NUCAM,NUCBFC)

      call check_nuc_ovlap(npbf)

      call calc_GAM_epcore(nebf,npebf,npbf,nebf*nebf,npbf*npbf,
     x                     nat,pmass,zan,cat,
     x                     AMPEB2C,AGEBFCC,AGNBFCC,
     x                     ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      ng1=nebf*nebf*npbf*npbf
      call calc_GAM_ep(nebf,npebf,npbf,ng1,
     x                 AMPEB2C,AGEBFCC,AGNBFCC,
     x                 ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

      call calc_GAM_ee(nebf,npebf,ngee,
     x                 AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

      wtime1 = MPI_WTIME() - wtime

C Set dummy variables for serial XCHF code call
      intdum = 1
      arrdum = 0.0d+00
      logdum = .true.

C Set SCF variables
      nebf2=nebf*nebf
      npbf2=npbf*npbf
      nebflt=nebf*(nebf+1)/2
      npbflt=npbf*(npbf+1)/2
      if(nelec.gt.1) then
       npr=(nebf-(nelec/2))*(nelec/2) ! Number occ-vir pairs
      else
       npr=nelec*(nebf-nelec)
      end if

C Set UHF variables
      if(LXCUHF) then
       npra=nae*(nebf-nae)
       nprb=nbe*(nebf-nbe)
      end if

C SCF call
      wtime  = MPI_WTIME()

      if(LXCUHF) then
       call xcuscf(nelec,nae,nbe,npra,nprb,nebflt,nucst,
     x             npebf,nebf,nebf2,npbf,npbf2,ngee,
     x             ngtg,ng1,intdum,intdum,intdum,intdum,intdum,intdum,
     x             read_CE,read_CP,
     x             logdum,logdum,LDBG,LSOSCF,
     x             intdum,intdum,nat,pmass,cat,zan,bgem,ggem,
     x             KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x             ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x             intdum,arrdum,arrdum,
     x             intdum,arrdum,
     x             intdum,arrdum)
      else
       call xcscf(nelec,npr,nebflt,nucst,
     x            npebf,nebf,nebf2,npbf,npbf2,ngee,
     x            ngtg,ng1,intdum,intdum,intdum,intdum,intdum,intdum,
     x            read_CE,read_CP,
     x            logdum,logdum,logdum,logdum,logdum,LSOSCF,LDBG,
     x            intdum,intdum,nat,pmass,cat,zan,bgem,ggem,
     x            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x            logdum,intdum,arrdum,arrdum,
     x            logdum,intdum,arrdum,
     x            logdum,intdum,arrdum)
      end if

      wtime2 = MPI_WTIME() - wtime

C Print timing summary
      write(*,*)
      write(*,*) "FINISHED NEO-HF CALCULATION"
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

