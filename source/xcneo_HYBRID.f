!=======================================================================
      program xcneo_hybrid

!=======================================================================
      implicit none
!     include 'mpif.h'
      include 'omp_lib.h'

!     integer nproc,myid,ierr


!     call mpi_set(nproc,myid,ierr)

!     call xcneo_hybrid_driver(nproc,myid)
      call xcneo_driver

!     call mpi_final(ierr)

      end
!=======================================================================
      subroutine xcneo_driver

!=======================================================================
      implicit none
!     include 'mpif.h'
      include 'omp_lib.h'

!     integer nproc,myid
!-------Basis Set Info-------(
      integer,allocatable :: ELCAM(:,:)  ! Angular mom for electrons
      integer,allocatable :: NUCAM(:,:)   ! Angular mom for quantum nuclei
      double precision,allocatable :: ELCEX(:) ! Exponents: elec basis
      double precision,allocatable :: NUCEX(:)  ! Exponents: nuc basis
      double precision,allocatable :: ELCBFC(:,:) ! Basis centers: elec basis
      double precision,allocatable :: NUCBFC(:,:)  ! basis centers: nuc basis
      integer,allocatable :: AMPEB2C(:) ! Map primitive index to contracted
      double precision,allocatable :: AGEBFCC(:) ! Map prim index to contract coef
      double precision,allocatable :: AGNBFCC(:)  ! Nuclear contract coef
      integer,allocatable :: KPESTR(:)  ! Map contracted index to primitive start
      integer,allocatable :: KPEEND(:)  ! Map contracted index to primitive end
      double precision,allocatable :: zan(:) ! Classical nuclear charges
      double precision,allocatable :: cat(:,:) ! XYZ Coordinates of classical atoms
      double precision,allocatable :: bcoef1(:)
      double precision,allocatable :: gamma1(:)
      integer nat
      integer ngtg1
      integer nelec
      integer NAE
      integer NAalpE,NAbetE
      integer NBE
      double precision pmass    ! Mass of nonelectron quantum particle 
!-------Basis Set Info-------)
      double precision,allocatable :: GM2ICR(:)
      double precision,allocatable :: GM2sICR(:)
      double precision,allocatable :: GM3IC1(:)
      double precision,allocatable :: GM4ICR(:)
      double precision,allocatable :: GM2_1ICR(:)
      double precision,allocatable :: GM2_2ICR(:)
      double precision,allocatable :: GM3_1IC1(:)
      double precision,allocatable :: GM3_2IC1(:)
      integer SZG2ICR
      integer SZG3IC1
      integer SZG4IC
      integer i,j,idum,istat
      integer NUCST
      integer NG2CHK
      integer NG3CHK
      integer NG4CHK
      logical LNEOHF
      logical LXCUHF
      logical LXCROHF
      logical LG2DSCF
      logical LG3DSCF
      logical LG2IC1 ! GAM2 Integrals In-Core algo-1: contracted
      logical LG3IC1 ! GAM3 Integrals In-Core algo-1
      logical LG3IC2 ! GAM3 Integrals In-Core algo-2
      logical LG4DSCF
      logical LG4IC  ! GAM4 Integrals In-Core 
      logical read_CE
      logical read_CP
      logical read_GAM2
      logical read_GAM3
      logical read_GAM4
      logical LGAM4
      logical LCMF
      logical LSOSCF
      logical LRXCHF
      logical LRXCUHF
      logical LOCBSE

      double precision a2bohr,bohr2a
      parameter(bohr2a=0.529177249d+00)
      parameter(a2bohr=1.0d+00/0.529177249d+00)

      integer npebf,nebf,npbf
      integer nebf2,npbf2,NPR,NEBFLT
      integer NPRA,NPRB

      integer ngee
      integer ng1,ng2,ng3,ng4
      integer ng1prm,ng2prm,ng3prm,ng4prm

      double precision wtime,wtime1,wtime2


       wtime = omp_get_wtime()


!-------READ-INPUT-FILE-AND-ALLOCATE-MEMORY-FOR-BASIS-SET--------------(
!     if(myid.eq.0) then
         open(unit=9,file='basis_definition.inp')

         read(9,*)ngtg1
         if(allocated(bcoef1)) deallocate(bcoef1)
         allocate( bcoef1(ngtg1),stat=istat )
         if(allocated(gamma1)) deallocate(gamma1)
         allocate( gamma1(ngtg1),stat=istat )
         do i=1,ngtg1
            read(9,*)bcoef1(i),gamma1(i)
         end do

         read(9,*)nat
         if(allocated(zan)) deallocate(zan)
         allocate( zan(nat),stat=istat )
         if(allocated(cat)) deallocate(cat)
         allocate( cat(3,nat),stat=istat )
         do i=1,nat
            read(9,*)zan(i),cat(1,i),cat(2,i),cat(3,i)
            cat(1,i)=a2bohr*cat(1,i)
            cat(2,i)=a2bohr*cat(2,i)
            cat(3,i)=a2bohr*cat(3,i)
         end do

         read(9,*)nebf
         read(9,*)npebf
         if(allocated(AMPEB2C)) deallocate(AMPEB2C)
         allocate( AMPEB2C(npebf),stat=istat )
         if(allocated(ELCEX)) deallocate(ELCEX)
         allocate( ELCEX(npebf),stat=istat )
         if(allocated(AGEBFCC)) deallocate(AGEBFCC)
         allocate( AGEBFCC(npebf),stat=istat )
         if(allocated(ELCAM)) deallocate(ELCAM)
         allocate( ELCAM(npebf,3),stat=istat )
         if(allocated(ELCBFC)) deallocate(ELCBFC)
         allocate( ELCBFC(npebf,3),stat=istat )
         do i=1,npebf
            read(9,*)idum,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),ELCAM(i,3),
     x    ELCEX(i),AGEBFCC(i),ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
            do j=1,3
               ELCBFC(i,j)=a2bohr*ELCBFC(i,j)
            end do
         end do

         read(9,*)npbf
         if(allocated(NUCEX)) deallocate(NUCEX)
         allocate( NUCEX(npbf),stat=istat )
         if(allocated(AGNBFCC)) deallocate(AGNBFCC)
         allocate( AGNBFCC(npbf),stat=istat )
         if(allocated(NUCAM)) deallocate(NUCAM)
         allocate( NUCAM(npbf,3),stat=istat )
         if(allocated(NUCBFC)) deallocate(NUCBFC)
         allocate( NUCBFC(npbf,3),stat=istat )
         do i=1,npbf
            read(9,*)idum,idum,NUCAM(i,1),NUCAM(i,2),NUCAM(i,3),
     x    NUCEX(i),AGNBFCC(i),NUCBFC(i,1),NUCBFC(i,2),NUCBFC(i,3)
            do j=1,3
               NUCBFC(i,j)=a2bohr*NUCBFC(i,j)
            end do
         end do

         read(9,*)pmass
         read(9,*)nelec
         read(9,*)NAE
         read(9,*)NBE
         read(9,*) NUCST
         read(9,*) LNEOHF
         read(9,*) LXCUHF
         read(9,*) LXCROHF
         read(9,*) LRXCHF
         read(9,*) LRXCUHF
         read(9,*) read_CE
         read(9,*) read_CP
         read(9,*) read_GAM2
         read(9,*) read_GAM3
         read(9,*) read_GAM4
         read(9,*) NG4CHK
         read(9,*) LGAM4
         read(9,*) LG4DSCF
         read(9,*) LG4IC
         read(9,*) NG3CHK
         read(9,*) LG3DSCF
         read(9,*) LG3IC1
         read(9,*) LG3IC2
         read(9,*) NG2CHK
         read(9,*) LG2DSCF
         read(9,*) LG2IC1
         read(9,*) LCMF
         read(9,*) LSOSCF
         read(9,*) LOCBSE

         close(9)

         if(allocated(KPESTR)) deallocate(KPESTR)
         allocate( KPESTR(nebf),stat=istat )
         if(allocated(KPEEND)) deallocate(KPEEND)
         allocate( KPEEND(nebf),stat=istat )
         call make_KPE(nebf,npebf,AMPEB2C,KPESTR,KPEEND)

      ngee=nebf*nebf*nebf*nebf

      ng1=nebf*nebf*npbf*npbf
      ng2=nebf*nebf*nebf*nebf*npbf*npbf
      ng3=nebf*nebf*nebf*nebf*nebf*nebf*npbf*npbf
      ng4=nebf*nebf*nebf*nebf*nebf*nebf*nebf*nebf*npbf*npbf

      ng1prm=npebf*npebf*npbf*npbf
      ng2prm=npebf*npebf*npebf*npebf*npbf*npbf
      ng3prm=npebf*npebf*npebf*npebf*npebf*npebf*npbf*npbf
c     ng4prm=npebf*npebf*npebf*npebf*npebf*npebf*npebf*npebf*npbf*npbf
      ng4prm=1

      if (((LRXCHF).or.(LRXCUHF)).and.(nelec.ge.4)) then
       write(*,*) "Currently only nelec<4 is supported"
       write(*,*) "Exiting..."
       return
      end if

      if((LRXCHF).and.(nae.eq.1)) LRXCUHF=.false.  ! So that PsH hacks work
      if(LRXCHF.and.LRXCUHF) then
       LRXCHF=.false.
       write(*,*) "Overriding LRXCHF since LRXCUHF=.TRUE."
      end if
      if(LNEOHF.and.LRXCHF) then
       LRXCHF=.false.
       write(*,*) "Overriding LRXCHF since LNEOHF=.TRUE."
      end if
      if(LNEOHF.and.LRXCUHF) then
       LRXCUHF=.false.
       LXCUHF=.false.
       write(*,*) "Overriding LRXCUHF since LNEOHF=.TRUE."
      end if

      if(LRXCUHF) then
       write(*,*) "Replacing NAE/NBE accordingly with LRXCUHF request"

C Ensure num beta reg elecs > num alpha reg elecs since special electron is assigned spin alpha
       if (NAE.le.NBE) then  
        NAalpE=NAE
        NAbetE=NBE
       else
        NAalpE=NBE
        NAbetE=NAE
       end if

       NAE=NAalpE+NAbetE     ! From now, NAE = num regular electrons
       NBE=1                 ! From now, NBE = num special electrons = 1

      end if

      write(*,*)
      write(*,*)'nat   =',nat
      write(*,*)'npebf =',npebf
      write(*,*)'nebf  =',nebf
      write(*,*)'npbf  =',npbf
      write(*,*)'ngtg1 =',ngtg1
      write(*,*)'ngee  =',ngee
      write(*,*)'ng1   =',ng1
      write(*,*)'ng2   =',ng2
      write(*,*)'ng3   =',ng3
      write(*,*)'ng4   =',ng4
      write(*,*)'ng1prm=',ng1prm
      write(*,*)'ng2prm=',ng2prm
      write(*,*)'ng3prm=',ng3prm
c     write(*,*)'ng4prm=',ng4prm
      write(*,*)
      write(*,*)'PMASS   =',PMASS
      write(*,*)'nelec   =',nelec

      if ((LRXCHF).or.(LRXCUHF)) then
       write(*,*)'NAE     =',NAE,'= total number of regular electrons'
       write(*,*)'NBE     =',NBE,'= number of special electrons'
      else
       write(*,*)'NAE     =',NAE
       write(*,*)'NBE     =',NBE
      end if

      write(*,*)'NUCST   =',NUCST
      write(*,*)'LNEOHF  =',LNEOHF
      write(*,*)'LXCUHF  =',LXCUHF
      write(*,*)'LXCROHF =',LXCROHF
      write(*,*)'LRXCHF  =',LRXCHF
      write(*,*)'LRXCUHF  =',LRXCUHF
      write(*,*)'read_CE =',read_CE
      write(*,*)'read_CP =',read_CP
      write(*,*)'READ_GAM2=',read_GAM2
      write(*,*)'READ_GAM3=',read_GAM3
      write(*,*)'READ_GAM4=',read_GAM4
      write(*,*)'NG4CHK  =',NG4CHK
      write(*,*)'LGAM4   =',LGAM4
      write(*,*)'LG4DSCF =',LG4DSCF
      write(*,*)'LG4IC   =',LG4IC 
      write(*,*)'NG3CHK  =',NG3CHK
      write(*,*)'LG3DSCF =',LG3DSCF
      write(*,*)'LG3IC1  =',LG3IC1
      write(*,*)'LG3IC2  =',LG3IC2
      write(*,*)'NG2CHK  =',NG2CHK
      write(*,*)'LG2DSCF =',LG2DSCF
      write(*,*)'LG2IC1  =',LG2IC1
      write(*,*)'LCMF    =',LCMF
      write(*,*)'LSOSCF  =',LSOSCF
      write(*,*)'LOCBSE  =',LOCBSE
      if(LRXCUHF) then
       write(*,*)'NAalpE =',NAalpE,'= number of alpha regular electrons'
       write(*,*)'NAbetE =',NAbetE,'= number of beta regular electrons'
      end if
      write(*,*) "Geminal parameters: k, b_k, gamm_k"
         do i=1,ngtg1
            write(*,*) i,bcoef1(i),gamma1(i)
         end do
      write(*,*)


         write(*,*)' CHECK CONTRACTED ELECTRONIC BASIS FUNCTIONS '
         write(*,*)'CONT INDEX    KPESTR     KPEEND'
         do i=1,nebf
            write(*,8000) i,KPESTR(i),KPEEND(i)
         end do

         WRITE(*,*)
         WRITE(*,*)'ELECTRONIC BASIS FUNCTIONS:'
         WRITE(*,*)
      WRITE(*,*)'PRIM  CONT    ANG       EXPONENT CONTRACT  -X- -Y- -Z-'
         WRITE(*,*)'INDEX INDEX   MOM                  COEF'
         DO i=1,npebf
           WRITE(*,9000) i,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),ELCAM(i,3),
     x    ELCEX(i),AGEBFCC(i),ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
         END DO
! Normalize the contraction coefficients for elec basis functions
!        call ELCNORM(npebf,AGEBFCC,ELCEX,ELCAM,ELCBFC)
!        call ELCNORM2(npebf,nebf,
!    x                 AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)
!        call ELCNORM3(npebf,nebf,KPESTR,KPEEND,
!    x                 AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)
         call ELCNORM3(npebf,nebf,
     x                 AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

         WRITE(*,*)'ELECTRONIC BASIS FUNCTIONS:'
         WRITE(*,*)'CONTRACT COEFF HAVE BEEN NORMALIZED'
         WRITE(*,*)
      WRITE(*,*)'PRIM  CONT    ANG       EXPONENT CONTRACT  -X- -Y- -Z-'
         WRITE(*,*)'INDEX INDEX   MOM                  COEF'
         DO i=1,npebf
           WRITE(*,9000) i,AMPEB2C(i),ELCAM(i,1),ELCAM(i,2),ELCAM(i,3),
     x    ELCEX(i),AGEBFCC(i),ELCBFC(i,1),ELCBFC(i,2),ELCBFC(i,3)
         END DO

         WRITE(*,*)
         WRITE(*,*)'NUCLEAR BASIS FUNCTIONS:'
         WRITE(*,*)
      WRITE(*,*)'PRIM  CONT    ANG       EXPONENT CONTRACT  -X- -Y- -Z-'
         WRITE(*,*)'INDEX INDEX   MOM                  COEF'
         DO i=1,npbf
           WRITE(*,9000) i,i,NUCAM(i,1),NUCAM(i,2),NUCAM(i,3),
     x    NUCEX(i),AGNBFCC(i),NUCBFC(i,1),NUCBFC(i,2),NUCBFC(i,3)
         END DO

! Normalize the contraction coefficients for nuc basis functions
         call NUCNORM(npbf,AGNBFCC,NUCEX,NUCAM,NUCBFC)

         WRITE(*,*)
         WRITE(*,*)'NUCLEAR BASIS FUNCTIONS:'
         WRITE(*,*)'CONTRACT COEFF HAVE BEEN NORMALIZED'
         WRITE(*,*)
      WRITE(*,*)'PRIM  CONT    ANG       EXPONENT CONTRACT  -X- -Y- -Z-'
         WRITE(*,*)'INDEX INDEX   MOM                  COEF'
         DO i=1,npbf
           WRITE(*,9000) i,i,NUCAM(i,1),NUCAM(i,2),NUCAM(i,3),
     x    NUCEX(i),AGNBFCC(i),NUCBFC(i,1),NUCBFC(i,2),NUCBFC(i,3)
         END DO
!-------READ-INPUT-FILE-AND-ALLOCATE-MEMORY-FOR-BASIS-SET--------------)

         ng1prm=(npebf**2)*npbf**2
         ng2prm=(npebf**4)*npbf**2
         ng3prm=(npebf**6)*npbf**2
!        ng4prm=(npebf**8)*npbf**2

!----CALCULATE-INEXPENSIVE-INTEGRALS-ON-MASTER-NODE--------------------(
! NOTE:  GEMINAL INTEGRALS WILL STILL BE OVER AVAILABLE OMP THREADS
! Standard NEO-HF integrals:
         call class_nuc_rep(nat,zan,cat)

         call elec_ovlap(npebf,nebf,nebf*nebf,
     x                   AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

         call check_elec_ovlap(nebf)

         call nuc_ovlap(npbf,npbf*npbf,
     x                  AGNBFCC,NUCEX,NUCAM,NUCBFC)

         call check_nuc_ovlap(npbf)

         call calc_GAM_epcore(nebf,npebf,npbf,nebf*nebf,npbf*npbf,
     x                        nat,pmass,zan,cat,
     x                        AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

         write(*,*)'all done in calc_GAM_epcore'

         call calc_GAM_ep(nebf,npebf,npbf,ng1,
     x                    AMPEB2C,AGEBFCC,AGNBFCC,
     x                    ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

         write(*,*)'all done in calc_GAM_ep'

         call calc_GAM_ee(nebf,npebf,ngee,
     x                    AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

         write(*,*)'all done in calc_GAM_ee'
!----CALCULATE-INEXPENSIVE-INTEGRALS-ON-MASTER-NODE--------------------)

         if(LNEOHF.or.(nelec.lt.2)) then
            SZG2ICR=1
            if(allocated(GM2ICR)) deallocate(GM2ICR)
            allocate( GM2ICR(SZG2ICR),stat=istat )
            if(allocated(GM2_1ICR)) deallocate(GM2_1ICR)
            allocate( GM2_1ICR(SZG2ICR),stat=istat )
            if(allocated(GM2_2ICR)) deallocate(GM2_2ICR)
            allocate( GM2_2ICR(SZG2ICR),stat=istat )
            if(allocated(GM2sICR)) deallocate(GM2sICR)
            allocate( GM2sICR(SZG2ICR),stat=istat )
         end if

         if(LNEOHF.or.(nelec.lt.3)) then
            SZG3IC1=1
            if(allocated(GM3IC1)) deallocate(GM3IC1)
            allocate( GM3IC1(SZG3IC1),stat=istat )
            if(allocated(GM3_1IC1)) deallocate(GM3_1IC1)
            allocate( GM3_1IC1(SZG3IC1),stat=istat )
            if(allocated(GM3_2IC1)) deallocate(GM3_2IC1)
            allocate( GM3_2IC1(SZG3IC1),stat=istat )
         end if

         if(LNEOHF.or.(nelec.le.3)) then
            SZG4IC=1
            if(allocated(GM4ICR)) deallocate(GM4ICR)
            allocate( GM4ICR(SZG4IC),stat=istat )
         end if

         if(.NOT.LNEOHF) then

          if ((LRXCHF).or.(LRXCUHF)) then
            call RXCHF_GAM1_OMP_MD(nebf,npebf,npbf,ng1,ng1prm,nat,ngtg1,
     x                       pmass,cat,zan,bcoef1,gamma1,
     x                       AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
          else
            call GAM1_OMP_MD(nebf,npebf,npbf,ng1,ng1prm,nat,ngtg1,
     x                       pmass,cat,zan,bcoef1,gamma1,
     x                       AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
          end if

            if(nelec.gt.1) then

              if(.not.(read_GAM2.or.LG2DSCF)) then

                if(LG2IC1) then

                   SZG2ICR=ng2

                   if(allocated(GM2sICR)) deallocate(GM2sICR)
                   allocate( GM2sICR(SZG2ICR),stat=istat )

                   if ((LRXCHF).or.(LRXCUHF)) then

                       if(allocated(GM2_1ICR)) deallocate(GM2_1ICR)
                       allocate( GM2_1ICR(SZG2ICR),stat=istat )
                       if(allocated(GM2_2ICR)) deallocate(GM2_2ICR)
                       allocate( GM2_2ICR(SZG2ICR),stat=istat )

                       call RXCHF_GAM2_IC1(NG2CHK,nebf,npebf,
     x                           npbf,ng2,ng2prm,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                           GM2_1ICR,GM2_2ICR,GM2sICR)

C Hacks for 2 electron (one regular electron) case in a singlet arrangement:
C   - Set all GM2sICR = 0
C   - Set all GM2_2ICR = 0
C   - Multiply all GM2_1ICR by 1/2
                       if (nae.eq.1) then

                         write(*,*) "Performing hack for 1 reg elec"
                         call RXCHF_adjust_omg2_ints(ng2,GM2_1ICR,
     x                                                   GM2_2ICR,
     x                                                   GM2sICR)

                       end if


                   else

                       if(allocated(GM2ICR)) deallocate(GM2ICR)
                       allocate( GM2ICR(SZG2ICR),stat=istat )

                       call GAM2_IC1(NG2CHK,nebf,npebf,npbf,
     x                           ng2,ng2prm,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                           GM2ICR,GM2sICR)
                   end if

                else

!                 call GAM2_OMP_MD(nebf,npebf,npbf,ng2,ng2prm,nat,ngtg1,
!    x                             pmass,cat,zan,bcoef1,gamma1,
!    x                             AMPEB2C,AGEBFCC,AGNBFCC,
!    x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

                   if ((LRXCHF).or.(LRXCUHF)) then
                      write(*,*) "************************************"
                      write(*,*) "Only in-core GAM2 coded for RXC(U)HF"
                      write(*,*) "************************************"
                      return
                   end if

                   call GAM2_CONV(NG2CHK,nebf,npebf,npbf,
     x                       ng2,ng2prm,nat,ngtg1,
     x                       pmass,cat,zan,bcoef1,gamma1,
     x                       KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

                   SZG2ICR=1
                   if(allocated(GM2ICR)) deallocate(GM2ICR)
                   allocate( GM2ICR(SZG2ICR),stat=istat )
                   if(allocated(GM2SICR)) deallocate(GM2SICR)
                   allocate( GM2SICR(SZG2ICR),stat=istat )

                end if !end if for LG2IC1

              end if ! end if for read and LG2DSCF
            end if ! end if for nelec.gt.1

              if(nelec.gt.2) then
                if(.not.(read_GAM3.or.LG3DSCF)) then
      
                  if(LG3IC1) then

                     SZG3IC1=ng3

                   if ((LRXCHF).or.(LRXCUHF)) then

                       if(allocated(GM3_1IC1)) deallocate(GM3_1IC1)
                       allocate( GM3_1IC1(SZG3IC1),stat=istat )
                       if(allocated(GM3_2IC1)) deallocate(GM3_2IC1)
                       allocate( GM3_2IC1(SZG3IC1),stat=istat )

                       call RXCHF_GAM3_IC1(NG3CHK,nebf,npebf,npbf,
     x                            ng3,ng3prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            GM3_1IC1,GM3_2IC1)

                   else

                       if(allocated(GM3IC1)) deallocate(GM3IC1)
                       allocate( GM3IC1(SZG3IC1),stat=istat )

                       call GAM3_IC1(NG3CHK,nebf,npebf,npbf,
     x                            ng3,ng3prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            GM3IC1)
                   end if

                  else if(LG3IC2) then

                   if ((LRXCHF).or.(LRXCUHF)) then
                      write(*,*) "************************************"
                      write(*,*) "Only in-core GAM3 coded for RXC(U)HF"
                      write(*,*) "************************************"
                      return
                   end if

                     SZG3IC1=ng3

                     if(allocated(GM3IC1)) deallocate(GM3IC1)
                     allocate( GM3IC1(SZG3IC1),stat=istat )

                     call GAM3_IC2(NG3CHK,nebf,npebf,npbf,
     x                            ng3,ng3prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                            GM3IC1)

                  else 

!                  call GAM3_OMP_MD(nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
!    x                              pmass,cat,zan,bcoef1,gamma1,
!    x                              AMPEB2C,AGEBFCC,AGNBFCC,
!    x                              ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
      
                   if ((LRXCHF).or.(LRXCUHF)) then
                      write(*,*) "************************************"
                      write(*,*) "Only in-core GAM3 coded for RXC(U)HF"
                      write(*,*) "************************************"
                      return
                   end if

                     call GAM3_CONV(NG3CHK,nebf,npebf,npbf,
     x                            ng3,ng3prm,nat,ngtg1,
     x                            pmass,cat,zan,bcoef1,gamma1,
     x                            KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                            ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

                     SZG3IC1=1
                     if(allocated(GM3IC1)) deallocate(GM3IC1)
                     allocate( GM3IC1(SZG3IC1),stat=istat )

                   end if

                end if

                if(nelec.gt.3) then

                   if(LGAM4.and.(.not.LG4DSCF)) then

                      if(LG4IC) then

                         if(.NOT.LG2IC1) then
                             write(*,*)'LG2IC1 MUST BE TRUE IF LG4IC'
                             return
                         end if

                         SZG4IC=ng4

                         if(allocated(GM4ICR)) deallocate(GM4ICR)
                         allocate( GM4ICR(SZG4IC),stat=istat )

                    ! UHF/RXCHF calculation requires different OMG4 integrals:
                         if ((LXCUHF).or.(LRXCHF).or.(LRXCUHF)) then
                          call OMG4_ICR(NG4CHK,nebf,npbf,ngee,ng2,ng4,
     x                                  GM2SICR,GM4ICR)
                         else
!                   call GAM4_ICR(NG4CHK,nebf,npbf,ngee,ng2,ng4,GM4ICR)
                            call GAM4_ICR(NG4CHK,nebf,npbf,ngee,ng2,ng4,
     x                                    GM2SICR,GM4ICR)
                         end if ! end if for LXCUHF or RXCHF or RXCUHF

                      end if ! end if for LG4IC

                   else

                      SZG4IC=1
                      if(allocated(GM4ICR)) deallocate(GM4ICR)
                      allocate( GM4ICR(SZG4IC),stat=istat )

                   end if

                end if !endif for nelec.gt.3

              end if  ! end if for nelec gt 2

         end if  !end if for NOT.neohf


         nebf2=nebf*nebf
         npbf2=npbf*npbf
         if (LRXCHF) then
          if(nae.gt.1) then
             NPRA=(nebf-(nae/2))*(nae/2) ! OCC-VIR PAIRS for reg elecs
          else
             NPRA=nae*(nebf-nae)
          end if
         else
          if(nelec.gt.1) then
             NPR=(nebf-(nelec/2))*(nelec/2) !Number OCC-VIR PAIRS
          else
             NPR=nelec*(nebf-nelec)
          end if
         end if
         NEBFLT=nebf*(nebf+1)/2

      if(LRXCUHF) then
         NPRA=NAalpE*(nebf-NAalpE)
         NPRB=NAbetE*(nebf-NAbetE)
      end if

      if(LXCUHF) then
! nelec = 1 not allowed for XCUHF:
         if(nelec.eq.1) then
            write(*,*)'XCUHF and NELEC=1 NOT PERMITTED, ABORTING'
            RETURN
         end if
         nebf2=nebf*nebf
         npbf2=npbf*npbf
         NPRA=NAE*(nebf-NAE)
         NPRB=NBE*(nebf-NBE)
         NEBFLT=nebf*(nebf+1)/2
      end if

!     end if ! end if for myid=0


         wtime1 = omp_get_wtime() - wtime
         wtime  = omp_get_wtime()

!-------------------NEO-SCF--------------------------------------------(
!     call xcscf(myid,nproc,nelec,NPR,NEBFLT,NUCST,
!    x           nebf,nebf2,npbf,npbf2,ngee,
!    x           ngtg1,ng1,ng2,ng3,ng4,NG3CHK,NG4CHK,
!    x           read_CE,read_CP,LNEOHF,
!    x           pmass,cat,zan,bcoef1,gamma1,
!    x           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
!    x           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)


! Hack(
         if(LG3IC2) LG3IC1=.TRUE.
! Hack)

         if(LXCUHF) then

            call xcuscf(nelec,NAE,NBE,NPRA,NPRB,NEBFLT,NUCST,
     x                  npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                  ngtg1,ng1,ng2,ng3,ng4,NG4CHK,NG3CHK,NG2CHK,
     x                  read_CE,read_CP,
     x                  LNEOHF,LGAM4,LCMF,LSOSCF,
     x                  ng2prm,ng3prm,nat,pmass,cat,zan,bcoef1,gamma1,
     x                  KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                  ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                  SZG2ICR,GM2ICR,GM2SICR,
     x                  SZG3IC1,GM3IC1,
     x                  SZG4IC,GM4ICR)


         elseif(LXCROHF) then

            write(*,*)'NO ROHF YET!'
            RETURN

         elseif(LRXCHF) then

            call xcrxchf(nelec,NAE,NBE,NPRA,NEBFLT,NUCST,
     x                   npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                   ngtg1,ng1,ng2,ng3,ng4,NG2CHK,NG3CHK,NG4CHK,
     x                   read_CE,read_CP,
     x                   LNEOHF,LGAM4,LCMF,LSOSCF,LOCBSE,
     x                   ng2prm,ng3prm,nat,pmass,cat,zan,bcoef1,gamma1,
     x                   KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                   ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                   SZG2ICR,GM2_1ICR,GM2_2ICR,GM2sICR,
     x                   LG3IC1,SZG3IC1,GM3_1IC1,GM3_2IC1,
     x                   LG4IC,SZG4IC,GM4ICR)

         elseif(LRXCUHF) then

            call xcrxcuhf(nelec,NAalpE,NAbetE,NBE,NPRA,NPRB,NEBFLT,
     x                    NUCST,npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                    ngtg1,ng1,ng2,ng3,ng4,NG2CHK,NG3CHK,NG4CHK,
     x                    read_CE,read_CP,
     x                    LNEOHF,LGAM4,LCMF,LSOSCF,LOCBSE,
     x                    ng2prm,ng3prm,nat,pmass,cat,zan,bcoef1,gamma1,
     x                    KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                    ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                    SZG2ICR,GM2_1ICR,GM2_2ICR,GM2sICR,
     x                    LG3IC1,SZG3IC1,GM3_1IC1,GM3_2IC1,
     x                    LG4IC,SZG4IC,GM4ICR)

         else

            call xcscf(nelec,NPR,NEBFLT,NUCST,
     x                 npebf,nebf,nebf2,npbf,npbf2,ngee,
     x                 ngtg1,ng1,ng2,ng3,ng4,NG4CHK,NG3CHK,NG2CHK,
     x                 read_CE,read_CP,
     x                 LNEOHF,LGAM4,LG4DSCF,LG3DSCF,LG2DSCF,LCMF,
     x                 ng2prm,ng3prm,nat,pmass,cat,zan,bcoef1,gamma1,
     x                 KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                 ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                 LG2IC1,SZG2ICR,GM2ICR,GM2SICR,
     x                 LG3IC1,SZG3IC1,GM3IC1,
     x                 LG4IC,SZG4IC,GM4ICR)
         end if
!-------------------NEO-SCF--------------------------------------------)

         if(allocated(GM2ICR)) deallocate(GM2ICR)
         if(allocated(GM2_1ICR)) deallocate(GM2_1ICR)
         if(allocated(GM2_2ICR)) deallocate(GM2_2ICR)
         if(allocated(GM2sICR)) deallocate(GM2sICR)
         if(allocated(GM3IC1)) deallocate(GM3IC1)
         if(allocated(GM3_1IC1)) deallocate(GM3_1IC1)
         if(allocated(GM3_2IC1)) deallocate(GM3_2IC1)
         if(allocated(GM4ICR)) deallocate(GM4ICR)

         wtime2 = omp_get_wtime() - wtime
         write(*,3000) wtime1,wtime2



 3000 FORMAT(/8X,'  +--------------------------------------+',/,
     X        8X,'  |    TIMING SUMMARY FOR CALCULATION    |',/,
     x        8X,'  +--------------------------------------+',/,
     x        8X,'    TIME TO EVALUATE INTEGRALS:',1X,F12.4/
     x        8X,'                  TIME FOR SCF:',1X,F12.4/)

 8000 format(/1X,I3,I6,I5)
 9000 format(/1X,I3,I6,I5,I3,I3,F12.6,F10.6,F10.6,F10.6,F10.6)


      return
      end
!=======================================================================
      subroutine mpi_set(nproc,myid,ierr)
! Setup for MPI
!=======================================================================
      implicit none
!     include 'mpif.h'
! Variables Returned
      integer nproc,myid,ierr

!     call MPI_INIT(ierr)
!     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
!     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

      return
      end
!=======================================================================
      subroutine mpi_final(ierr)

!=======================================================================
      implicit none
!     include 'mpif.h'
! Variables Returned
      integer ierr

!     call MPI_FINALIZE(ierr)

      return
      end

!======================================================================
      subroutine make_KPE(nebf,npebf,AMPEB2C,KPESTR,KPEEND)
 
!     Create a map between contracted EBF index and:
!     1) Beginning primitive index of the contracted shell:  KPESTR()
!     2) Ending primitive index of the contracted shell: KPEEND()
!======================================================================
      implicit none
! Input Variables
      integer nebf,npebf
      integer AMPEB2C(npebf)
! Variables Returned
      integer KPESTR(nebf)
      integer KPEEND(nebf)
! Local Variables
      integer ichk,ist,iend
      integer iep,iec


      ichk=1
      ist=1
      do iep=1,npebf
         iec=AMPEB2C(iep)
         if(iec.ne.ichk) then
            ichk=ichk+1
            iend=iep-1
            KPESTR(iec-1)=ist
            KPEEND(iec-1)=iend
            ist=iep
         end if
         if(iep.eq.npebf) then
            iend=iep
            KPESTR(iec)=ist
            KPEEND(iec)=iend
         end if
      end do


      return
      end
!=======================================================================
      subroutine NUCNORM(npbf,AGNBFCC,NUCEX,NUCAM,NUCBFC)
! Normalize the contraction coefficients of nuclear basis set
!=======================================================================
      implicit none
! Input Variables
      integer npbf
      integer NUCAM(npbf,3)  ! Angular mom for quantum nuclei
      double precision NUCEX(npbf) ! Exponents: nuc basis
      double precision NUCBFC(npbf,3) ! Basis centers: nuc basis
! Variables Returned
      double precision AGNBFCC(npbf) ! Map prim index to contract coef
! Local Variables
      integer ip
      integer I1,J1,K1
      double precision ans,A1,Amat1(3)

      do ip=1,npbf

         A1=NUCEX(ip)
         I1=NUCAM(ip,1)
         J1=NUCAM(ip,2)
         K1=NUCAM(ip,3)
!        Amat1(1)=NUCBFC(ip,1)
!        Amat1(2)=NUCBFC(ip,2)
!        Amat1(3)=NUCBFC(ip,3)
         Amat1(1)=0.0d+00
         Amat1(2)=0.0d+00
         Amat1(3)=0.0d+00

         call gfovlap(I1,J1,K1,A1,Amat1,
     2                I1,J1,K1,A1,Amat1,
     3                ans)

         ans=1.0d+00/sqrt(ans)
         AGNBFCC(ip)=AGNBFCC(ip)*ans
      end do

      return
      end

