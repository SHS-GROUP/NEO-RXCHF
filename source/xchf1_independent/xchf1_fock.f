!======================================================================
      subroutine xchf1_fock(LNEOHF,LGAM4,LG4DSCF,LG4IC,
     x                      LG3DSCF,LG3IC1,LG2DSCF,LG2IC1,LCMF,
     x                      NG4CHK,NG3CHK,NG2CHK,
     x                      SZG4IC,SZG3IC1,SZG2ICR,
     x                      npebf,nebf,nebf2,npbf,npbf2,nelec,
     x                      ngee,ng1,ng2,ng3,ng4,DE,DP,
     x                      GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
     x                      GM4ICR,GM3IC1,GM2ICR,GM2SICR,
     x                      ng2prm,ngtg1,ng3prm,
     x                      nat,pmass,cat,zan,bcoef1,gamma1,
     x                      KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                      ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                      focke,fockp,Xfocke,Xfockp,
     x                      E_total,E_nuc,E_ecore,E_pcore,E_ep,
     x                      E_ee,E_gam1,E_gam2,E_gam3,E_gam4,
     x                      S_total,S_gam1,S_gam2)

!======================================================================
      implicit none

! Input Variables
      logical LCMF     ! Check My Fock: turns on on-the-fly Fock check
      logical LNEOHF   ! Control for doing NEO-HF only
      logical LG4DSCF  ! Control for direct SCF for GAM4 term
      logical LG3DSCF  ! Control for direct SCF for GAM3 term
      logical LG2IC1   ! Control for SCF for in-core GAM2 term
      logical LG3IC1   ! Control for SCF for in-core GAM3 term
      logical LG4IC    ! Control for SCF for in-core GAM4 term
      logical LG2DSCF  ! Control for direct SCF for GAM2 term
      integer NG4CHK   ! Number of chunks for dividing GAM4 evaluation
      integer NG3CHK   ! Number of chunks for dividing GAM3 evaluation
      integer NG2CHK   ! Number of chunks for dividing GAM2 evaluation
      integer SZG2ICR
      integer SZG3IC1
      integer SZG4IC 
      integer nebf
      integer nebf2
      integer npbf
      integer npbf2
      integer nelec
!     integer ngtgmn

      integer ngee
      integer ng1
      integer ng2
      integer ng3
      integer ng4

      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision GAM_ecore(nebf2)
      double precision GAM_pcore(npbf2)
      double precision GAM_ep(ng1)
      double precision GAM_ee(ngee)
      double precision GM2ICR(SZG2ICR)
      double precision GM2SICR(SZG2ICR)
      double precision GM3IC1(SZG3IC1)
      double precision GM4ICR(SZG4IC)
!     double precision efct(nebf)
!     double precision pfct(npbf)
!     double precision GAM_1s(ng1)
!     double precision GAM_1(ng1)
!     double precision GAM_2s(ng2)
!     double precision GAM_2(ng2)
!     double precision GAM_3PK(ng3)
!     double precision GAM_4PK(ng4)
!     double precision xxse(nebf,nebf)
!     double precision xxsp(npbf,npbf)
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

! Variables Returned
      double precision focke(nebf,nebf)
      double precision fockp(npbf,npbf)
      double precision Xfocke(nebf,nebf)
      double precision Xfockp(npbf,npbf)
      double precision E_total
      double precision E_ecore
      double precision E_pcore
      double precision E_ep
      double precision E_ee
      double precision E_gam1
      double precision E_gam2
      double precision E_gam3
      double precision E_gam4
      double precision E_nuc
      double precision S_total
      double precision S_gam1
      double precision S_gam2

! Local variables
      integer npbfLT
      integer nebfLT
      double precision zero,one
      parameter(zero=0.0d+00,one=1.0d+00)
      double precision H_expect 
      double precision E_check 
      double precision E_NEOHF 
      double precision psiHpsi 
      double precision se1(nebf,nebf)
      double precision se2(nebf,nebf)
      double precision sp1(npbf,npbf)
      double precision sp2(npbf,npbf)
      double precision CHKfockp(npbf,npbf)
      double precision CHKS(npbf,npbf)

      logical LGAM4

!----------------(
!     LGAM4=.FALSE.
!     LGAM4=.TRUE.
!----------------)


! Initialize values for energy components
      psiHpsi=zero

      E_total=zero
      E_ecore=zero
      E_pcore=zero
      E_ep=zero
      E_ee=zero
      E_gam1=zero
      E_gam2=zero
      E_gam3=zero
      E_gam4=zero

      S_total=zero
      S_gam1=zero
      S_gam2=zero

! Zero-out the Fock matrices
!     call zero_out(nebf,focke)
!     call zero_out(npbf,fockp)
      focke=zero
      fockp=zero

! Zero-out the XFock matrices
!     call zero_out(nebf,Xfocke)
!     call zero_out(npbf,Xfockp)
      Xfocke=zero
      Xfockp=zero

! Zero-out the Fock-like overlap matrices
!     call zero_out(nebf,se1)
!     call zero_out(npbf,sp1)
      se1=zero
      sp1=zero

!     call zero_out(nebf,se2)
!     call zero_out(npbf,sp2)
      se2=zero
      sp2=zero


! Calculate energy components and build Fock matrices

!  Conventional NEO-HF components:
      call E_from_GAM_ecore(nebf,nebf2,
     * GAM_ecore,DE,focke,E_ecore)

!     nebfLT=(nebf*nebf+nebf)/2
!     write(*,*)'XCHF SIDE... Elec HCORE Matrix:'
!     call print_my_fock(nebf,nebfLT,focke)
!     write(*,*)'GAM_ecore integrals'
!     write(*,*) GAM_ecore

      call E_from_GAM_pcore(npbf,npbf2,
     * GAM_pcore,DP,fockp,E_pcore)

!     npbfLT=(npbf*npbf+npbf)/2
!     write(*,*)'XCHF SIDE... Nuc HCORE Matrix:'
!     call print_my_fock(npbf,npbfLT,fockp)

      call E_from_GAM_ep(nebf,npbf,ng1,
     * GAM_ep,DE,DP,focke,fockp,E_ep)

!     write(*,*)'XCHF SIDE... Nuc Fock Matrix:'
!     call print_my_fock(npbf,npbfLT,fockp)
!     write(*,*)'GAM_ep integrals'
!     write(*,*) GAM_ep

      if(nelec.gt.1) then
         call E_from_GAM_ee(nebf,ngee,
     *    GAM_ee,DE,focke,E_ee)
      else
         E_ee=0.0d+00
      end if
!         nebfLT=(nebf*nebf+nebf)/2
!         write(*,*)'XCHF SIDE... Elec Fock Matrix:'
!         call print_my_fock(nebf,nebfLT,focke)
!         call abrt

! Retain the standard NEOHF values for energy and Fock matrices
!  for use later
      Xfocke=focke
      Xfockp=fockp
      E_NEOHF=E_ecore+E_pcore+E_ep+E_ee+E_nuc
!         write(*,*)'NEOHF Energy is:',E_NEOHF

! IF NEO-HF ONLY GET OUT OF HERE:
      if(LNEOHF) then 
         E_total=E_NEOHF
         S_total=1.0d+00
!-----------ON-THE-FLY-CODE-TESTING------------------------------------(
         if(LCMF) then
            call Fock_sym_check(nebf,npbf,focke,fockp)
            call NonGem_OVLAP_check(nebf,npbf,nelec,DE,DP)
         end if
!-----------ON-THE-FLY-CODE-TESTING------------------------------------)
         RETURN
      end if
!  XC components:
!     if(ngtgmn.gt.0) then

         call S_from_GAM_1s(nebf,npbf,ng1,
     *    DE,DP,se1,sp1,S_gam1)

         call E_from_GAM_1(nebf,npbf,ng1,
     *    DE,DP,focke,fockp,E_gam1)

         if(nelec.gt.1) then

            if(LG2DSCF) then
               call GAM2_DSCF(NG2CHK,nebf,npebf,npbf,
     x                        ng2,ng2prm,nat,ngtg1,
     x                        pmass,cat,zan,bcoef1,gamma1,
     x                        KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                        ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                        DE,DP,focke,fockp,se2,sp2,E_gam2,S_gam2)
            else if(LG2IC1) then
!     subroutine EG2IC1(Nchunks,nebf,npbf,ng2,
!    x                  DE,DP,GM2ICR,GM2SICR,
!    x                  FE,SE,FP,SP,E_gam2,S_gam2)
               call EG2IC1(NG2CHK,nebf,npbf,ng2,
     x                     DE,DP,GM2ICR,GM2SICR,
     x                     focke,se2,fockp,sp2,E_gam2,S_gam2)
            else
               call S_from_GAM_2s(nebf,npbf,ng2,
     *         DE,DP,se2,sp2,S_gam2)

               call E_from_GAM_2(nebf,npbf,ng2,
     *         DE,DP,focke,fockp,E_gam2)
            end if

            if(nelec.gt.2) then

!              write(*,*)'about to call E_from_GAM_3'
!              call flshbf(6)
               if(LG3DSCF) then
                  call GAM3_DSCF(NG3CHK,nebf,npebf,npbf,
     x                           ng3,ng3prm,nat,ngtg1,
     x                           pmass,cat,zan,bcoef1,gamma1,
     x                           KPESTR,KPEEND,AMPEB2C,AGEBFCC,AGNBFCC,
     x                           ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     x                           DE,DP,focke,fockp,E_gam3)
               else if(LG3IC1) then
                  call EG3IC1(NG3CHK,nebf,npbf,ng3,
     x                        DE,DP,GM3IC1,focke,fockp,E_gam3)
               else
                  call E_from_GAM_3(nebf,npbf,ng3,
     x             DE,DP,focke,fockp,E_gam3)
               end if
!              write(*,*)'about to call E_from_GAM_4'
!              call flshbf(6)
               if(LGAM4) then
                  if(LG4DSCF) then
                     call GAM4_DSCF(NG4CHK,nebf,npbf,ngee,ng2,ng4,
     x                              DE,DP,focke,fockp,E_gam4)
                  else if(LG4IC) then
                     call EG4ICR(NG4CHK,nebf,npbf,ng4,
     x                           DE,DP,GM4ICR,focke,fockp,E_gam4)
                  else
                     call E_from_GAM_4(nebf,npbf,ng4,
     x                DE,DP,focke,fockp,E_gam4)
                  end if
               end if

            end if

         end if

!     end if

!  Calculate total overlap
! CWS-TEST the non-geminal overlap contribution is included in S_gam1
      S_total=one+S_gam1+S_gam2
!     S_total=S_gam1+S_gam2

!  Calculate total energy
      psiHpsi=(E_ecore+E_pcore+E_ep+E_ee+
     x        E_gam1+E_gam2+E_gam3+E_gam4)/S_total

      E_total=psiHpsi+E_nuc

!     H_expect=E_ecore+E_pcore+E_ep+E_ee+
!    x         E_gam1+E_gam2+E_gam3+E_gam4
      H_expect=psiHpsi/S_total

!  Make XC-related corrections to the Fock matrices
      call XCfock_correction(nebf,npbf,
     *            E_total,S_total,H_expect,
     *            se1,sp1,se2,sp2,focke,fockp,CHKfockp,CHKS)

!-----------ON-THE-FLY-CODE-TESTING------------------------------------(
      if(LCMF) then
         npbfLT=npbf*(npbf+1)/2
         nebfLT=nebf*(nebf+1)/2
!!       E_check=psiHpsi-E_ecore-E_ee
!!       E_check=E_ecore+E_pcore+E_ep+E_ee+E_gam1+E_gam2+E_gam3+E_gam4
         E_check=E_pcore+E_ep+E_gam1+E_gam2+E_gam3+E_gam4
         call CHK_my_fockp(npbf,npbfLT,E_check,S_total,DP,
     x                     CHKfockp,CHKS)

         call Focking_Around(nebf,nebfLT,npbf,npbfLT,
     x                       DE,DP,focke,fockp,se1,se2,sp1,sp2,
     x                       psiHpsi,S_total,E_pcore,E_ecore,
     x                       E_ep,E_ee,E_gam1,E_gam2,E_gam3,E_gam4)

         call Fock_sym_check(nebf,npbf,focke,fockp)
         call NonGem_OVLAP_check(nebf,npbf,nelec,DE,DP)

         write(*,*) "DE:"
         call prt_lower_triangle(nebf,nebflt,de)
         write(*,*)
         write(*,*) "DP:"
         call prt_lower_triangle(npbf,npbflt,dp)
         write(*,*)
         write(*,*) "FE:"
         call prt_lower_triangle(nebf,nebflt,focke)
         write(*,*)
         write(*,*) "FP:"
         call prt_lower_triangle(npbf,npbflt,fockp)

      end if
!-----------ON-THE-FLY-CODE-TESTING------------------------------------)

      return
      end
!=======================================================================
      subroutine Fock_sym_check(nebf,npbf,focke,fockp)

! Checks that the elec and nuc Fock matrices are symmetric
!=======================================================================
      implicit none
! Input Variables
      integer nebf,npbf
      double precision focke(nebf,nebf),fockp(npbf,npbf)
! Local variables
      logical symE,symP
      integer i,j
      double precision val,tolerance
      parameter(tolerance=1.0d-10)

      symE=.true.
      symP=.true.

      do i=1,nebf
         do j=1,i
            val=focke(i,j)-focke(j,i)
            if(abs(val).gt.tolerance) then
              write(*,*)'ELEC FOCK MATRIX IS NOT SYMMETRIC FOR IJ=',I,J
              write(*,*)'>>>> FE(IJ)=',focke(i,j)
              write(*,*)'>>>> FE(JI)=',focke(j,i)
              symE=.false.
            end if
         end do
      end do

      do i=1,npbf
         do j=1,i
            val=fockp(i,j)-fockp(j,i)
            if(abs(val).gt.tolerance) then
              write(*,*)'NUC FOCK MATRIX IS NOT SYMMETRIC FOR IJ=',I,J
              write(*,*)'>>>> FP(IJ)=',fockp(i,j)
              write(*,*)'>>>> FP(JI)=',fockp(j,i)
              symP=.false.
            end if
         end do
      end do

      if(symE.and.symP) then
         write(*,*)'NUC and ELEC Fock matrices are symmetric'
      end if


      return
      end
!======================================================================
      subroutine NonGem_OVLAP_check(nebf,npbf,nelec,DE,DP)

! Checks that the total NEO-HF overlap is unity
!======================================================================
      implicit none

! Input Variables
      integer nebf
      integer npbf
      integer nelec
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)

! Local variables
      logical EUNITY,PUNITY
      integer ia 
      integer ie1
      integer je1
      integer ip
      integer jp
      double precision val
      double precision SE_TOTAL,SP_TOTAL
      double precision zero
      parameter(zero=0.0d+00)
      double precision tolerance
      parameter(tolerance=1.0d-10)

      EUNITY=.true.
      PUNITY=.true.

      open(814,file='eovlap.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      SE_TOTAL=zero
      do ie1=1,nebf
         do je1=1,nebf
!           Map the 2-index contracted integral to 1-D:
            call pack_2D(nebf,je1,ie1,ia)
            read(814,REC=ia) val
            SE_TOTAL=SE_TOTAL+DE(ie1,je1)*val
         end do
      end do
      close(814)
! Divide by the number of electrons:
      SE_TOTAL=SE_TOTAL/dble(nelec)

      open(815,file='novlap.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      SP_TOTAL=zero
      do ip=1,npbf
         do jp=1,npbf
!           Map the 2-index contracted integral to 1-D:
            call pack_2D(npbf,jp,ip,ia)
            read(815,REC=ia) val
            SP_TOTAL=SP_TOTAL+DP(ip,jp)*val
         end do
      end do
      close(815)

      val=SE_TOTAL-1.0d+00
      if(abs(val).gt.tolerance) then
          write(*,*)'NEOHF ELEC OVLP PROBLEM'
          write(*,*)'SE_TOTAL=',SE_TOTAL
          write(*,*)'SP_TOTAL=',SP_TOTAL
          EUNITY=.false.
      end if

      val=SP_TOTAL-1.0d+00
      if(abs(val).gt.tolerance) then
          write(*,*)'NEOHF NUC OVLP PROBLEM'
          write(*,*)'SE_TOTAL=',SE_TOTAL
          write(*,*)'SP_TOTAL=',SP_TOTAL
          PUNITY=.false.
      end if

      if(EUNITY.and.PUNITY) then
         write(*,*)'NEOHF NORMALIZED TO UNITY'
      end if


      return
      end

!======================================================================
      subroutine CHK_my_fockp(npbf,npbfLT,E_check,S_total,DP,
     x           CHKfockp,CHKS)
!======================================================================
      implicit none

! Input Variables
      integer npbf
      integer npbfLT
      double precision E_check
      double precision S_total
      double precision DP(npbf,npbf)
      double precision CHKfockp(npbf,npbf)
      double precision CHKS(npbf,npbf)

! Local Variables
      integer i,j,iLT
      double precision E_TRP
      double precision S_TRP
      double precision LTfockp(npbfLT)
      double precision LTDP(npbfLT)
      double precision LTS(npbfLT)

      double precision TRACEP


      iLT=0
      do i=1,npbf
         do j=1,i

            iLT=iLT+1
            LTfockp(iLT)=CHKfockp(i,j)
            LTDP(iLT)=DP(i,j)
            LTS(iLT)=CHKS(i,j)

         end do
      end do

!c    E_TRP=TRACEP(DAE,SCR,L1)
!-This line commented out for independent version
      E_TRP=TRACEP(LTDP,LTfockp,npbf)
!-This line commented out for independent version
      S_TRP=TRACEP(LTDP,LTS,npbf)

      write(*,*)
      write(*,*)'>>>IN CHECK_MY_FOCK<<<'
      write(*,*)
      write(*,*)'E_check=',E_check
      write(*,*)'E_TRP=',E_TRP
      write(*,*)
      write(*,*)'S_total=',S_total
      write(*,*)'S_TRP=',S_TRP
      write(*,*)


      return
      end
!=======================================================================
      subroutine Focking_Around(nebf,nebfLT,npbf,npbfLT,
     x                          DE,DP,focke,fockp,se1,se2,sp1,sp2,
     x                          E_total,S_total,E_pcore,E_ecore,
     x                          E_ep,E_ee,E_gam1,E_gam2,E_gam3,E_gam4)
!=======================================================================
      implicit none

! Input Variables
      integer nebf
      integer npbf
      integer nebfLT
      integer npbfLT
      double precision S_total
      double precision E_total
      double precision E_pcore
      double precision E_ecore
      double precision E_ep
      double precision E_ee
      double precision E_gam1
      double precision E_gam2
      double precision E_gam3
      double precision E_gam4
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision focke(nebf,nebf)
      double precision fockp(npbf,npbf)
      double precision se1(nebf,nebf)
      double precision se2(nebf,nebf)
      double precision sp1(npbf,npbf)
      double precision sp2(npbf,npbf)

! Local Variables
      integer ie,je,iLT
      integer ip,jp
      double precision se(nebf,nebf)
      double precision sp(npbf,npbf)
      double precision SIGMAe_I
      double precision SIGMAe_II
      double precision SIGMAp_I
      double precision SIGMAp_II
      double precision LTfocke(nebfLT)
      double precision LTfockp(npbfLT)
      double precision LTDE(nebfLT)
      double precision LTDP(npbfLT)
      double precision LTSE(nebfLT)
      double precision LTSP(npbfLT)
      double precision E_XCHF
      double precision FCKTRC_e
      double precision FCKTRC_p
      double precision Check_e
      double precision Check_p
      double precision one,two,three,four
      parameter(one=1.0d+00,two=2.0d+00,three=3.0d+00,four=4.0d+00)

      double precision TRACEP


      do ie=1,nebf
         do je=1,nebf
            se(ie,je) = se1(ie,je) + two * se2(ie,je)
         end do
      end do


      do ip=1,npbf
         do jp=1,npbf
            sp(ip,jp) = sp1(ip,jp) + sp2(ip,jp)
         end do
      end do


      iLT=0
      do ie=1,nebf
         do je=1,ie

            iLT=iLT+1
            LTfocke(iLT)=focke(ie,je)
            LTDE(iLT)=DE(ie,je)
            LTSE(iLT)=se(ie,je)

         end do
      end do

      iLT=0
      do ip=1,npbf
         do jp=1,ip

            iLT=iLT+1
            LTfockp(iLT)=fockp(ip,jp)
            LTDP(iLT)=DP(ip,jp)
            LTSP(iLT)=sp(ip,jp)

         end do
      end do

      E_XCHF = E_ecore
     x       + E_pcore
     x       + E_ep
     x       + E_ee
     x       + E_gam1
     x       + E_gam2
     x       + E_gam3
     x       + E_gam4

      SIGMAe_I =         E_ecore
     x         +         E_ep
     x         + two *   E_ee
     x         +         E_gam1
     x         + two *   E_gam2
     x         + three * E_gam3
     x         + four *  E_gam4

      SIGMAe_II=TRACEP(LTDE,LTSE,nebf) 
!     SIGMAp_II=TRACEP(LTDP,LTSP,npbf) 

      FCKTRC_e=TRACEP(LTDE,LTfocke,nebf)
      FCKTRC_p=TRACEP(LTDP,LTfockp,npbf)

      Check_e= ( SIGMAe_I / S_total ) 
     x       - ( E_XCHF*SIGMAe_II / S_total**2 ) 

      Check_p= ( E_XCHF / S_total**2 ) 
     x       - ( E_ecore / S_total ) 
     x       - ( E_ee / S_total ) 


      write(*,*)
      write(*,*)'>>> FOCK CHECK <<<'
      write(*,*)
!     write(*,*)'S_total   = ',S_total
!     write(*,*)'E_XCHF    = ',E_XCHF
!     write(*,*)'SIGMAe_I  =',SIGMAe_I
!     write(*,*)'SIGMAe_II =',SIGMAe_II
      write(*,*)'===================================='
      write(*,*)'FCKTRC_e  = ',FCKTRC_e
      write(*,*)'Check_e   = ',Check_e
      write(*,*)
      write(*,*)'FCKTRC_p  = ',FCKTRC_p
      write(*,*)'Check_p   = ',Check_p
      write(*,*)'===================================='
      write(*,*)


      return
      end

!*MODULE MTHLIB  *DECK TRACEP
      DOUBLE PRECISION FUNCTION TRACEP(A,B,N)
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
      DIMENSION A(*),B(*)
!
!     ----- COMPUTE 1 ELECTRON EXPECTATION VALUE ----
!     THE DENSITY MATRIX IS STORED IN A IN TRIANGULAR ORDER.
!     NOTE THAT THE OFF DIAGONAL ELEMENTS OF THE DENSITY ARE
!     TO BE DOUBLED.
!     THE PROPERTY INTEGRALS ARE IN B IN TRIANGULAR ORDER.
!     THE EXPECTATION VALUE IS THE DOT PRODUCT OF THESE MATRICES.
!
      N2 = (N*N+N)/2
      TRACE = DDOT(N2,A,1,B,1)
!                 DOUBLE EVERYTHING, THEN SUBTRACT THE DIAGONAL
      TRACE = TRACE+TRACE
      K = 0
      DO 100 I = 1,N
         K = K+I
         TRACE = TRACE - A(K)*B(K)
  100 CONTINUE
      TRACEP = TRACE
      RETURN
      END

!======================================================================
      subroutine XCfock_correction(nebf,npbf,
     *            E_total,S_total,H_expect,
     *            se1,sp1,se2,sp2,focke,fockp,CHKfockp,CHKS)

!======================================================================
      implicit none

      integer nebf
      integer npbf

      double precision E_total
      double precision S_total
      double precision H_expect
      double precision se1(nebf,nebf)
      double precision sp1(npbf,npbf)
      double precision se2(nebf,nebf)
      double precision sp2(npbf,npbf)
      double precision focke(nebf,nebf)
      double precision fockp(npbf,npbf)
      double precision CHKfockp(npbf,npbf)
      double precision CHKS(npbf,npbf)

! Local variables
      integer ip
      integer jp
      integer ie1
      integer je1
      double precision one,two
      parameter(one=1.0d+00,two=2.0d+00)


      do ie1=1,nebf
         do je1=1,nebf

            focke(ie1,je1)=(one/S_total)*focke(ie1,je1)-
     x                    H_expect*(se1(ie1,je1)+two*se2(ie1,je1))

         end do
      end do

      do ip=1,npbf
         do jp=1,npbf

!           CHKfockp(ip,jp)=(one/S_total)*fockp(ip,jp)
            CHKfockp(ip,jp)=fockp(ip,jp)
            CHKS(ip,jp)=sp1(ip,jp)+sp2(ip,jp)

            fockp(ip,jp)=(one/S_total)*fockp(ip,jp)-
     x                    H_expect*(sp1(ip,jp)+sp2(ip,jp))

         end do
      end do


      return
      end

!======================================================================
      subroutine E_from_GAM_ecore(nebf,nebf2,
     * GAM_ecore,DE,focke,E_ecore)

!======================================================================
      implicit none

! Input Variables
      integer nebf,nebf2

      double precision focke(nebf,nebf)
      double precision DE(nebf,nebf)
!     double precision GAM_ecore(nebf,nebf)
      double precision GAM_ecore(nebf2)
      double precision E_ecore

! Local variables
      integer ie1,je1,ia
      double precision val_ecore
      double precision zero
      parameter(zero=0.0d+00)


      E_ecore=zero

      do ie1=1,nebf
         do je1=1,nebf

            call pack_2D(nebf,je1,ie1,ia)
!           val_ecore=GAM_ecore(ie1,je1)
            val_ecore=GAM_ecore(ia)
            focke(ie1,je1)=val_ecore
            E_ecore=E_ecore+DE(ie1,je1)*val_ecore

         end do
      end do


      return
      end

!======================================================================
      subroutine E_from_GAM_pcore(npbf,npbf2,
     * GAM_pcore,DP,fockp,E_pcore)

!======================================================================
      implicit none

! Input Variables
      integer npbf,npbf2
      double precision fockp(npbf,npbf)
      double precision DP(npbf,npbf)
!     double precision GAM_pcore(npbf,npbf)
      double precision GAM_pcore(npbf2)
      double precision E_pcore

! Local variables
      integer ip,jp,ia
      double precision val_pcore
      double precision zero
      parameter(zero=0.0d+00)


      E_pcore=zero

      do ip=1,npbf
         do jp=1,npbf

            call pack_2D(npbf,jp,ip,ia)
            val_pcore=GAM_pcore(ia)
!           val_pcore=GAM_pcore(ip,jp)
            fockp(ip,jp)=val_pcore
            E_pcore=E_pcore+DP(ip,jp)*val_pcore
!       write(*,*)'IN E_from_GAM_pcore'
!       write(*,*)'DP=',DP(ip,jp)
!       write(*,*)'val_pcore=',val_pcore
!       write(*,*)'E_pcore=',E_pcore

         end do
      end do


      return
      end

!======================================================================
      subroutine E_from_GAM_ee(nebf,ngee,
     * GAM_ee,DE,focke,E_ee)
!
!======================================================================
      implicit none

      double precision zero,two,half
      parameter(zero=0.0d+00,two=2.0d+00,half=5.0d-01)
      integer nebf
      integer ngee
      double precision GAM_ee(ngee)
      double precision DE(nebf,nebf)
      double precision focke(nebf,nebf)
      double precision E_ee

! Local variables
      integer nebfLT
      integer ia1,ia2
      integer ie1,je1,ie2,je2
      double precision vee1
      double precision vee2
      double precision val
      double precision xfocke(nebf,nebf)

      E_ee=zero

      call zero_out(nebf,xfocke)

      do ie2=1,nebf
         do je2=1,nebf

            do ie1=1,nebf
               do je1=1,nebf

                  call pack_4D(nebf,nebf,nebf,
     x                         je2,ie2,je1,ie1,ia1)

                  call pack_4D(nebf,nebf,nebf,
     x                         je1,ie2,je2,ie1,ia2)
   
                  vee1=GAM_ee(ia1)
                  vee2=GAM_ee(ia2)
!                 vee1=GAM_ee(ie1,je1,ie2,je2)
!                 vee2=GAM_ee(ie1,je2,ie2,je1)
                  val=(vee1-half*vee2)*half
                  xfocke(ie1,je1)=xfocke(ie1,je1)+two*DE(ie2,je2)*val
                  E_ee=E_ee+DE(ie1,je1)*DE(ie2,je2)*val

               end do
            end do

         end do
      end do

!     nebfLT=(nebf+nebf*nebf)/2
!     write(*,*)'IN E_from_GAM_ee: EE Contribution to Elec FOCK Matrix:'
!     call print_my_fock(nebf,nebfLT,xfocke)

!  Update the full electronic Fock matrix
      call add2fock(nebf,xfocke,focke)

      return
      end

!======================================================================
      subroutine E_from_GAM_ep(nebf,npbf,ng1,
     * GAM_ep,DE,DP,focke,fockp,E_ep)
 
!======================================================================
      implicit none

      integer nebf
      integer npbf
      integer ng1

      double precision GAM_ep(ng1)
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision focke(nebf,nebf)
      double precision fockp(npbf,npbf)
      double precision E_ep

! Local variables
      integer npbfLT
      integer nebfLT
      integer ia
      integer ie1,je1,ip,jp
      double precision val_ep
      double precision xfocke(nebf,nebf)
      double precision xfockp(npbf,npbf)
      double precision zero
      parameter(zero=0.0d+00)


      E_ep=zero
      
      call zero_out(nebf,xfocke)
      call zero_out(npbf,xfockp)

      do ip=1,npbf
         do jp=1,npbf

            do ie1=1,nebf
               do je1=1,nebf

                  call pack_4D(nebf,nebf,npbf,
     x                         je1,ie1,jp,ip,ia)

                  val_ep=GAM_ep(ia)
!                 val_ep=GAM_ep(ip,jp,ie1,je1)

                  xfocke(ie1,je1)=xfocke(ie1,je1)+DP(ip,jp)*val_ep
                  xfockp(ip,jp)=xfockp(ip,jp)+DE(ie1,je1)*val_ep
!                 write(*,*)'MAKING FOCK P: EP CONTRIBUTION'
!                 write(*,*)'DE(ie1,je1)=',DE(ie1,je1)
!                 write(*,*)'val_ep=',val_ep
!                 write(*,*)'xfockp(ip,jp)=',xfockp(ip,jp)
!                 write(*,*)

                  E_ep=E_ep+DE(ie1,je1)*DP(ip,jp)*val_ep

               end do
            end do

         end do
      end do

!     npbfLT=(npbf+npbf*npbf)/2
!     nebfLT=(nebf+nebf*nebf)/2
!     write(*,*)'IN E_from_GAM_ep: EP Contribution to Nuc FOCK Matrix:'
!     call print_my_fock(npbf,npbfLT,xfockp)
!     write(*,*)'IN E_from_GAM_ep: EP Contribution to Elec FOCK Matrix:'
!     call print_my_fock(nebf,nebfLT,xfocke)
!  Update the full Fock matrices
      call add2fock(nebf,xfocke,focke)
      call add2fock(npbf,xfockp,fockp)


      return
      end

!======================================================================
      subroutine E_from_GAM_1(nebf,npbf,ng1,
     *    DE,DP,focke,fockp,E_gam1)
!
!======================================================================
      implicit none

      double precision zero
      parameter(zero=0.0d+00)

      integer ng1
      integer nebf
      integer npbf

!     double precision GAM_1(ng1)
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision focke(nebf,nebf)
      double precision fockp(npbf,npbf)
      double precision E_gam1

! Local variables
      integer ia
      integer ie1
      integer je1
      integer ip
      integer jp
      double precision val_gam1
      double precision xfocke(nebf,nebf)
      double precision xfockp(npbf,npbf)

      E_gam1=zero
      
      call zero_out(nebf,xfocke)
      call zero_out(npbf,xfockp)

      open(801,file='GAM_1.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
 

      do ip=1,npbf
         do jp=1,npbf

            do ie1=1,nebf
               do je1=1,nebf

                  call pack_4D(nebf,nebf,npbf,
     x                         je1,ie1,jp,ip,ia)

!                 val_gam1=GAM_1(ip,jp,ie1,je1)
!                 val_gam1=GAM_1(ia)
                  read(801,REC=ia) val_gam1

                  xfocke(ie1,je1)=xfocke(ie1,je1)+DP(ip,jp)*val_gam1
                  xfockp(ip,jp)=xfockp(ip,jp)+DE(ie1,je1)*val_gam1

                  E_gam1=E_gam1+DE(ie1,je1)*DP(ip,jp)*val_gam1

               end do
            end do

         end do
      end do

!  Update the full Fock matrices
      call add2fock(nebf,xfocke,focke)
      call add2fock(npbf,xfockp,fockp)

      close(801)

      return
      end

!======================================================================
      subroutine S_from_GAM_1s(nebf,npbf,ng1,
     *    DE,DP,se1,sp1,S_gam1)
!
!======================================================================
      implicit none

      double precision zero
      parameter(zero=0.0d+00)

      integer ng1
      integer nebf
      integer npbf

!     double precision GAM_1s(ng1)
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision se1(nebf,nebf)
      double precision sp1(npbf,npbf)
      double precision S_gam1

! Local variables
      integer ia 
      integer ie1
      integer je1
      integer ip
      integer jp
      double precision val_gam1s

      S_gam1=zero
      
      call zero_out(nebf,se1)
      call zero_out(npbf,sp1)

      open(802,file='GAM_1s.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
 
      do ip=1,npbf
         do jp=1,npbf

            do ie1=1,nebf
               do je1=1,nebf

                  call pack_4D(nebf,nebf,npbf,
     x                         je1,ie1,jp,ip,ia)

!                 val_gam1s=GAM_1s(ip,jp,ie1,je1)
!                 val_gam1s=GAM_1s(ia)
                  read(802,REC=ia) val_gam1s

                  se1(ie1,je1)=se1(ie1,je1)+DP(ip,jp)*val_gam1s
                  sp1(ip,jp)=sp1(ip,jp)+DE(ie1,je1)*val_gam1s

                  S_gam1=S_gam1+DE(ie1,je1)*DP(ip,jp)*val_gam1s

               end do
            end do

         end do
      end do

      close(802)

      return
      end

!======================================================================
      subroutine E_from_GAM_2(nebf,npbf,ng2,
     *    DE,DP,focke,fockp,E_gam2)
!
!======================================================================
      implicit none

      double precision zero,half,two
      parameter(zero=0.0d+00,half=0.5d+00,two=2.0d+00)

      integer ng2
      integer nebf
      integer npbf

!CWS-IO
!     double precision GAM_2(ng2)
!     double precision GAM_2(1)
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision focke(nebf,nebf)
      double precision fockp(npbf,npbf)
      double precision E_gam2

! Local variables
      integer ia1
      integer ia2

      integer ie1
      integer je1
      integer ie2
      integer je2
      integer ip
      integer jp

      double precision val_gam2
      double precision xia1
      double precision xia2
      double precision xfocke(nebf,nebf)
      double precision xfockp(npbf,npbf)

      open(803,file='GAM_2.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
 

      E_gam2=zero
      
      call zero_out(nebf,xfocke)
      call zero_out(npbf,xfockp)

      do ip=1,npbf
      do jp=1,npbf

         do ie1=1,nebf
         do je1=1,nebf

            do ie2=1,nebf
            do je2=1,nebf

                call index_GAM_2PK(nebf,npbf,
     x                             ip,jp,
     x                             ie1,je1,
     x                             ie2,je2,ia1)

                call index_GAM_2PK(nebf,npbf,
     x                             ip,jp,
     x                             ie1,je2,
     x                             ie2,je1,ia2)

!cc         val_gam2=GAM_2(ip,jp,ie1,je1,ie2,je2)-
!cc  *          half*GAM_2(ip,jp,ie1,je2,ie2,je1)
!           val_gam2=GAM_2(ia1)-
!    *          half*GAM_2(ia2)

            read(803,REC=ia1) xia1
            read(803,REC=ia2) xia2
            val_gam2=xia1-half*xia2

            xfocke(ie1,je1)=xfocke(ie1,je1)+
     *             DP(ip,jp)*DE(ie2,je2)*two*val_gam2

            xfockp(ip,jp)=xfockp(ip,jp)+
     *             DE(ie1,je1)*DE(ie2,je2)*val_gam2

            E_gam2=E_gam2+DE(ie1,je1)*DE(ie2,je2)*DP(ip,jp)*val_gam2

            end do
            end do

         end do
         end do

      end do
      end do
!     do ie1=1,nebf
!     do je1=1,nebf

!        do ie2=1,nebf
!        do je2=1,nebf

!           do ip=1,npbf
!           do jp=1,npbf

!           val_gam2=GAM_2(ip,jp,ie1,je1,ie2,je2)-
!    *          half*GAM_2(ip,jp,ie1,je2,ie2,je1)

!           xfocke(ie1,je1)=xfocke(ie1,je1)+
!    *             DP(ip,jp)*DE(ie2,je2)*two*val_gam2

!           xfockp(ip,jp)=xfockp(ip,jp)+
!    *             DE(ie1,je1)*DE(ie2,je2)*val_gam2

!           E_gam2=E_gam2+DE(ie1,je1)*DE(ie2,je2)*DP(ip,jp)*val_gam2

!           end do
!           end do

!        end do
!        end do

!     end do
!     end do

!  Update the full Fock matrices
      call add2fock(nebf,xfocke,focke)
      call add2fock(npbf,xfockp,fockp)

      close(803)

      return
      end

!======================================================================
      subroutine S_from_GAM_2s(nebf,npbf,ng2,
     *    DE,DP,se2,sp2,S_gam2)
!
!======================================================================
      implicit none

      double precision zero,half
      parameter(zero=0.0d+00,half=0.5d+00)

      integer ng2
      integer nebf
      integer npbf

!CWS-IO
!     double precision GAM_2s(ng2)
!     double precision GAM_2s(1)
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision se2(nebf,nebf)
      double precision sp2(npbf,npbf)
      double precision S_gam2

! Local variables
      integer ia1
      integer ia2
      integer ie1
      integer je1
      integer ie2
      integer je2
      integer ip
      integer jp
      double precision val_gam2s
      double precision xia1
      double precision xia2

      open(804,file='GAM_2s.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)


      S_gam2=zero
      
      call zero_out(nebf,se2)
      call zero_out(npbf,sp2)

      do ip=1,npbf
      do jp=1,npbf

         do ie1=1,nebf
         do je1=1,nebf

            do ie2=1,nebf
            do je2=1,nebf

!           val_gam2s=GAM_2s(ip,jp,ie1,je1,ie2,je2)

                call index_GAM_2PK(nebf,npbf,
     x                             ip,jp,
     x                             ie1,je1,
     x                             ie2,je2,ia1)

                call index_GAM_2PK(nebf,npbf,
     x                             ip,jp,
     x                             ie1,je2,
     x                             ie2,je1,ia2)

!cc-new organization here
!cc-originally in gam_2.f
!  calc <ie1 ie2 ip|g(1)g(2)|je1 je2 jp> - 1/2<ie1 ie2 ip|g(1)g(2)|je2 je1 jp>
!                       val_gam2s=GAM_2s(ia1)-
!    *                       half*GAM_2s(ia2)

                        read(804,REC=ia1) xia1
                        read(804,REC=ia2) xia2
                        val_gam2s=xia1-half*xia2

!c                      val_gam2s=GAM_2s(ip,jp,ie1,je1,ie2,je2)-
!c   *                       half*GAM_2s(ip,jp,ie1,je2,ie2,je1)
!cc                     GAM_2s(ip,jp,ie1,je1,ie2,je2)=val_gam2s
!cc

!  gam_2s integrals are already this way (gam_2.f):
!           val_gam2s=GAM_2s(ip,jp,ie1,je1,ie2,je2)-
!    *           half*GAM_2s(ip,jp,ie1,je2,ie2,je1)

            se2(ie1,je1)=se2(ie1,je1)+
     *                   DP(ip,jp)*DE(ie2,je2)*val_gam2s
            sp2(ip,jp)=sp2(ip,jp)+
     *                   DE(ie1,je1)*DE(ie2,je2)*val_gam2s

            S_gam2=S_gam2+DE(ie1,je1)*DE(ie2,je2)*DP(ip,jp)*val_gam2s

            end do
            end do

         end do
         end do

      end do
      end do

      close(804)


      return
      end

!======================================================================
      subroutine E_from_GAM_3(nebf,npbf,ng3,
     *    DE,DP,focke,fockp,E_gam3)
!
!======================================================================
      implicit none

      double precision zero,half,two,three,four,eight
      parameter(zero=0.0d+00,half=0.5d+00,two=2.0d+00)
      parameter(three=3.0d+00,four=4.0d+00,eight=8.0d+00)

      integer nebf
      integer npbf
      integer ng3

!     double precision GAM_3PK(ng3)
!     double precision GAM_3PK(1)
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision focke(nebf,nebf)
      double precision fockp(npbf,npbf)
      double precision E_gam3

! Local variables
      integer ie1
      integer je1
      integer ie2
      integer je2
      integer ie3
      integer je3
      integer ip
      integer jp

      integer ia_123 
      integer ia_213
      integer ia_312
      integer ia_132
      integer ia_231
      integer ia_321

      double precision val_gam3
      double precision xfocke(nebf,nebf)
      double precision xfockp(npbf,npbf)

!CWS-IO
      double precision x123
      double precision x213
      double precision x312
      double precision x132
      double precision x231
      double precision x321
 
!CWS-IO
      open(805,file='GAM_3.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      E_gam3=zero
      
      call zero_out(nebf,xfocke)
      call zero_out(npbf,xfockp)

      do ip=1,npbf
      do jp=1,npbf

            do ie1=1,nebf
            do je1=1,nebf

                  do ie2=1,nebf
                  do je2=1,nebf

                     do ie3=1,nebf
                     do je3=1,nebf
!CCC
!    AS Packed--> GAM_3(ip,jp,je1,ie1,je2,ie2,je3,ie3)
!                   = GAM_3PK(je3,ie3,je2,ie2,je1,ie1,jp,ip)

      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je1,ie2,je2,ie3,je3,ia_123)
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je2,ie2,je1,ie3,je3,ia_213)
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je3,ie2,je1,ie3,je2,ia_312)
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je1,ie2,je3,ie3,je2,ia_132)
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je2,ie2,je3,ie3,je1,ia_231)
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je3,ie2,je2,ie3,je1,ia_321)

!CWS-IO  start
                        read(805,REC=ia_123) x123
                        read(805,REC=ia_213) x213
                        read(805,REC=ia_312) x312
                        read(805,REC=ia_132) x132
                        read(805,REC=ia_231) x231
                        read(805,REC=ia_321) x321

                  val_gam3=eight*x123 -
     x                      four*x213 +
     x                       two*x312 -
     x                      four*x132 +
     x                       two*x231 -
     x                      four*x321  
!CWS-IO  end

!                 val_gam3=eight*GAM_3PK(ia_123)-
!    x                      four*GAM_3PK(ia_213)+
!    x                       two*GAM_3PK(ia_312)-
!    x                      four*GAM_3PK(ia_132)+
!    x                       two*GAM_3PK(ia_231)-
!    x                      four*GAM_3PK(ia_321) 

!cc               val_gam3=eight*GAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3)-
!cc  x                      four*GAM_3(ip,jp,ie1,je2,ie2,je1,ie3,je3)+
!cc  x                       two*GAM_3(ip,jp,ie1,je3,ie2,je1,ie3,je2)-
!cc  x                      four*GAM_3(ip,jp,ie1,je1,ie2,je3,ie3,je2)+
!cc  x                       two*GAM_3(ip,jp,ie1,je2,ie2,je3,ie3,je1)-
!cc  x                      four*GAM_3(ip,jp,ie1,je3,ie2,je2,ie3,je1) 


!CCC
                  val_gam3=val_gam3/eight
!CWS-debug
!                    write(*,*)'ie1=',ie1
!                    write(*,*)'je1=',je1
!                    write(*,*)'ie2=',ie2
!                    write(*,*)'je2=',je2
!                    write(*,*)'ie3=',ie3
!                    write(*,*)'je3=',je3
!                    write(*,*)'vx =',val_gam3
!CWS-debug


            xfocke(ie1,je1)=xfocke(ie1,je1)+
     *             DP(ip,jp)*DE(ie2,je2)*DE(ie3,je3)*three*val_gam3

            xfockp(ip,jp)=xfockp(ip,jp)+
     *             DE(ie1,je1)*DE(ie2,je2)*DE(ie3,je3)*val_gam3

            E_gam3=E_gam3+
     x            DE(ie1,je1)*DE(ie2,je2)*DE(ie3,je3)*DP(ip,jp)*val_gam3

                     end do
                     end do

                  end do
                  end do

            end do
            end do

      end do
      end do

!CWS-IO
      close(805)

!CWS-debug
!     stop
!CWS-debug
!  Update the full Fock matrices
      call add2fock(nebf,xfocke,focke)
      call add2fock(npbf,xfockp,fockp)


      return
      end

!======================================================================
      subroutine E_from_GAM_4(nebf,npbf,ng4,
     *    DE,DP,focke,fockp,E_gam4)
!
!======================================================================
      implicit none
        
      double precision zero
      double precision half
      double precision two
      double precision three
      double precision four
      double precision eight 
      parameter(zero=0.0d+00,half=0.5d+00,two=2.0d+00)
      parameter(three=3.0d+00,four=4.0d+00,eight=8.0d+00)

      integer nebf
      integer npbf
      integer ng4

!CWS-IO
!     double precision GAM_4PK(ng4)
!     double precision GAM_4PK(1)
!     double precision efct(nebf)
!     double precision pfct(npbf)
      double precision DE(nebf,nebf)
      double precision DP(npbf,npbf)
      double precision focke(nebf,nebf)
      double precision fockp(npbf,npbf)
      double precision E_gam4

! Local variables

!CWS_screenG5
      logical Lread
      double precision tolerance
      parameter(tolerance=1.0d-02) 
      double precision ans_scr
!CWS_screenG5

      integer ie1
      integer je1
      integer ie2
      integer je2
      integer ie3
      integer je3
      integer ie4
      integer je4
      integer ip
      integer jp
      integer ia
      double precision val_gam4
      double precision xfocke(nebf,nebf)
      double precision xfockp(npbf,npbf)

      E_gam4=zero

 
!           write(*,*)'>>>i>>E_from_GAM_4<<<<<'
!           write(*,*)
!CWS-IO
      open(806,file='GAM_4.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      
      call zero_out(nebf,xfocke)
      call zero_out(npbf,xfockp)

      do ip=1,npbf
      do jp=1,npbf

        do ie1=1,nebf
        do je1=1,nebf

!CWS_screenG5
!         Lread=.true.
!         call screenG5(ip,jp,ie1,je1,
!    x                  nebf,npbf,efct,pfct,ans_scr)
!         ans_scr=abs(ans_scr)
!         if(ans_scr.lt.tolerance) Lread=.false.
!CWS_screenG5

          do ie2=1,nebf
          do je2=1,nebf

            do ie3=1,nebf
            do je3=1,nebf

              do ie4=1,nebf
              do je4=1,nebf

                 call index_GAM_4PK2(nebf,npbf,ip,jp,
     *                                        ie1,je1,
     *                                        ie2,je2,
     *                                        ie3,je3,
     *                                        ie4,je4,ia)
                  
!CWS-IO
!                 read(404,*)ia,val_gam4
!                 read(404,REC=ia) val_gam4
!CWS_screenG5
                  val_gam4=zero
!                 if(Lread) read(806,REC=ia) val_gam4
                  read(806,REC=ia) val_gam4
!CWS_screenG5
!CWS-EXPERIMENT
!                 val_gam4=0.0d+00
!                 write(*,*)'ia= ',ia,' val_gam4= ',val_gam4
!                 val_gam4=GAM_4PK(ia)
                  val_gam4=val_gam4/(four*four)
!           write(*,*)'ip =',ip
!           write(*,*)'jp =',jp
!           write(*,*)'ie1=',ie1
!           write(*,*)'je2=',je2
!           write(*,*)'ie2=',ie2
!           write(*,*)'je2=',je2
!           write(*,*)'ie3=',ie3
!           write(*,*)'je3',je3
!           write(*,*)'ie4=',ie4
!           write(*,*)'je4=',je4
!           write(*,*)'ia =',ia
!           write(*,*)'gam4(',ia,')=',val_gam4
!           write(*,*)

            xfocke(ie1,je1)=xfocke(ie1,je1)+
     *       DP(ip,jp)*DE(ie2,je2)*DE(ie3,je3)*DE(ie4,je4)*four*val_gam4

            xfockp(ip,jp)=xfockp(ip,jp)+
     *       DE(ie1,je1)*DE(ie2,je2)*DE(ie3,je3)*DE(ie4,je4)*val_gam4

            E_gam4=E_gam4+
     x       DE(ie1,je1)*DE(ie2,je2)*
     x       DE(ie3,je3)*DE(ie4,je4)*DP(ip,jp)*val_gam4

              end do
              end do

            end do
            end do

          end do
          end do

        end do
        end do

      end do
      end do

!CWS-IO
      close(806)

!  Update the full Fock matrices
      call add2fock(nebf,xfocke,focke)
      call add2fock(npbf,xfockp,fockp)


      return
      end

!======================================================================
      subroutine zero_out(nbf,amat)
!
!======================================================================

      parameter(zero=0.0d+00)

      integer nbf
      double precision amat(nbf,nbf)

! Local variables
      integer i
      integer j

      do i=1,nbf
         do j=1,nbf

            amat(i,j)=zero

         end do
      end do


      return
      end

!======================================================================
      subroutine add2fock(nbf,f_part,f_full)
!
!======================================================================

      integer nbf
      double precision f_full(nbf,nbf)
      double precision f_part(nbf,nbf)

! Local variables
      integer i
      integer j

      do i=1,nbf
         do j=1,nbf

            f_full(i,j)=f_full(i,j)+f_part(i,j)

         end do
      end do


      return
      end

!=======================================================================
      subroutine get_nuc_density(N1,N2,XDP,DP)

!=======================================================================
      implicit none

      integer IR,IW,IP,IS,IPK,IDAF,NAV,IODA
      integer N1,N2
      double precision XDP(N2),DP(N1,N1)

      integer ii,jj,kk

!-This line commented out for independent version
!cc   COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)


! Read lower-triangle-packed Nuclear density from disk
!-This line commented out for independent version
!cc   CALL daread(IDAF,IODA,XDP,N2,445,0)

! Expand lower triangle packed array to full symmetric matrix

      kk=0            
      do ii=1,N1
         do jj=1,ii
            kk=kk+1
            DP(ii,jj)=XDP(kk)
            DP(jj,ii)=XDP(kk)
         end do
      end do


      return
      end

!=======================================================================
      subroutine get_elec_density(L1,L2,XDE,DE)

!=======================================================================
      implicit none

      integer L1,L2
      double precision XDE(L2),DE(L1,L1)

      integer ii,jj,kk


! Expand lower triangle packed array to full symmetric matrix

      kk=0            
      do ii=1,L1
         do jj=1,ii
            kk=kk+1
            DE(ii,jj)=XDE(kk)
            DE(jj,ii)=XDE(kk)
         end do
      end do


      return
      end

!======================================================================
      subroutine Fock_Difference(nbf,fock,Xfock,Dfock)
!======================================================================
      implicit none
! Input Variables
      integer nbf
      double precision fock(nbf,nbf)
      double precision Xfock(nbf,nbf)
! Variables Returned
      double precision Dfock(nbf,nbf)
! Local Variables
      integer i
      integer j

      do i=1,nbf
         do j=1,nbf
            Dfock(i,j)=fock(i,j)-Xfock(i,j)
         end do
      end do

      return
      end

!======================================================================
      subroutine Pack_Fock(nbf,nbfLT,Dfock,Pfock)
!======================================================================
      implicit none
! Input Variables
      integer nbf
      integer nbfLT
      double precision Dfock(nbf,nbf)
! Variables Returned
      double precision Pfock(nbfLT)
! Local Variables
      integer i
      integer j
      integer ia

      ia=0
      do i=1,nbf
         do j=1,i
            ia=ia+1
            Pfock(ia)=Dfock(i,j)
         end do
      end do

      return
      end

!======================================================================
      subroutine Add_Fock(nLT,Pfock,GMSFCK)
!======================================================================
      implicit none
! Input Variables
      integer nLT
      double precision Pfock(nLT)
! Variables Returned
      double precision GMSFCK(nLT)
! Local Variables
      integer ia

      ia=0
      do ia=1,nLT
         GMSFCK(ia)=GMSFCK(ia)+Pfock(ia)
      end do

      return
      end
!
!     subroutine get_nuc_density(N1,N2,XDP,DP)
!     subroutine CHK_my_fockp(npbf,npbfLT,E_check,S_total,DP,
!    x           CHKfockp,CHKS)
!
!
!=======================================================================
!     subroutine MEM1_xchf1_fock(NFCK,nebfLT,nLT,XDE,GMSFCK,EXCHF1)

!=======================================================================
!     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     PARAMETER (MXATM=2000)
!     parameter (mxrt=100)

!     dimension XDE(nebfLT)
!     dimension GMSFCK(nLT)

!     COMMON /ENRGYS/ ENUC,EELEC,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
!    *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2)
!     COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
!    *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
!     COMMON /NUCMOI/ NUNIQN,IUNIQN(20),IUNIQT(20),NQMNUC,IQMNUC(20),
!    *                IQNTYP(20),NUMNB,NUCST,NAUXNB,IAUXNB(20),NUMULT,
!    *                NNA,NNB,NTAUXB
!     COMMON /XXXXXX/ NPEBFX

! NFCK :: Controls which Fock matrix is to be updated
!         1=> Electronic
!         2=> Nuclear
!
! GMSFCK :: The Fock matrix to be updated---either nuclear or electronic
! nLT    :: Lower triangle dimension of Fock matrix to be updated
! XDE    :: Lower triangle electronic density matrix
! nebfLT :: Lower triangle dimension of electronic matrices
!
! nelec  :: number of electrons
! npebf  :: number of primitive electronic basis functions
!           currently calculated in apbf.f and cludgely passed
!           through a common block (XXXX) to here...
! nebf   :: number of contracted electronic basis functions
! nebfLT :: dimension for lower-triangle packing for contracted 
!           electronic matrix elements
! ngam1  :: number of gamma_1 integrals (2-particle e-p integrals - G2)
! ngam2  :: number of gamma_2 integrals (3-particle e-p integrals - G3)
! ngam3  :: number of gamma_3 integrals (4-particle e-p integrals - G4)
! ngam4  :: number of gamma_4 integrals (5-particle e-p integrals - G5)
 
!     write(*,*)'In MEM1_xchf1_fock'
!     write(*,*)'npebf=',npebfx
!     write(*,*)'nebfLT=',nebfLT
!     call flshbf(6)
!     nelec=ne
!     nebf=num
!     npebf=NPEBFX
!     npbf=numnb

!     nebf2=nebf*nebf
!     npbf2=npbf*npbf
!     npbfLT=(npbf*npbf+npbf)/2
!     ngee=(nebf**4)
!     ngam1=(nebf**2)*npbf*npbf
!     ngam2=(nebf**4)*npbf*npbf
!     ngam3=(nebf**6)*npbf*npbf
!     ngam4=(nebf**8)*npbf*npbf

!     call MEM2_xchf1_fock(NFCK,nelec,npebf,nebf,npbf,nebf2,npbf2,
!    x                     nebfLT,npbfLT,nLT,
!    x                     ngee,ngam1,ngam2,ngam3,ngam4,
!    x                     enuc,XDE,GMSFCK,EXCHF1)


!     return
!     end

!=======================================================================
!     subroutine MEM2_xchf1_fock(NFCK,nelec,npebf,nebf,npbf,nebf2,npbf2,
!    x                           nebfLT,npbfLT,nLT,
!    x                           ngee,ngam1,ngam2,ngam3,ngam4,
!    x                           e_nuc,XDE,GMSFCK,EXCHF1)

!=======================================================================
!     implicit none

! Input Variables
!     integer NFCK
!     integer nelec
!     integer npebf
!     integer nebf
!     integer npbf
!     integer nebf2
!     integer npbf2
!     integer nebfLT
!     integer npbfLT
!     integer nLT
!     integer ngee
!     integer ngam1
!     integer ngam2
!     integer ngam3
!     integer ngam4
!     double precision e_nuc
!     double precision XDE(nebfLT)
!     double precision GMSFCK(nLT)
! Variables Returned
!     double precision EXCHF1

! Local Variables
!     double precision DE(nebf,nebf)
!     double precision DP(npbf,npbf)
!     double precision XDP(npbfLT)
!     double precision GAM_ecore(nebf2)
!     double precision GAM_pcore(npbf2)
!     double precision GAM_ep(ngam1)
!     double precision GAM_ee(ngee)
!     double precision focke(nebf,nebf)  ! Total XCHF E Fock Matrix
!     double precision fockp(npbf,npbf)  ! Total XCHF P Fock Matrix
!     double precision Xfocke(nebf,nebf) ! NEOHF E Fock Matrix
!     double precision Xfockp(npbf,npbf) ! NEOHF P Fock Matrix
!     double precision Dfocke(nebf,nebf) ! XCHF_FOCK - NEOHF_Fock for E
!     double precision Dfockp(npbf,npbf) ! XCHF_FOCK - NEOHF_Fock for P
!     double precision Pfocke(nebfLT) ! XCHF_FOCK - NEOHF_Fock for E LT
!     double precision Pfockp(npbfLT) ! XCHF_FOCK - NEOHF_Fock for P LT
      
!     double precision E_NEOHF
!c    double precision E_check1
!c    double precision E_check2
!     double precision E_total,E_ecore,E_pcore,E_ep
!     double precision E_ee,E_gam1,E_gam2,E_gam3,E_gam4
!     double precision S_total,S_gam1,S_gam2

!c    write(*,*)'In MEM2_xchf1_fock'
!c    write(*,*)'npebf=',npebf
!c    call flshbf(6)


!     write(*,*)'about to call get_nuc_density'
!     call get_nuc_density(npbf,npbfLT,XDP,DP)
!c    write(*,*)'back from get_nuc_density'
!c    write(*,*)'Back from get_nuc_density'
!c    write(*,*)'In MEM2_xchf1_fock Nuc Density is:'
!c    call print_my_fock(npbf,npbfLT,DP)
!c    write(*,*) 
!c    call get_elec_density(nebf,nebfLT,XDE,DE)
!c    write(*,*)'Back from get_elec_density'
!c    write(*,*)'In MEM2_xchf1_fock elec Density is:'
!c    call print_my_fock(nebf,nebfLT,DE)
!c    write(*,*) 
!c    write(*,*)'Calling calc_GAM_epcore'
!     call calc_GAM_epcore(nebf,npebf,npbf,nebf2,npbf2,
!    x                     GAM_ecore,GAM_pcore)
!c    write(*,*)'Back from calc_GAM_epcore'
!c    write(*,*)'Calling calc_GAM_ee'
!     call calc_GAM_ee(nebf,npebf,ngee,GAM_ee)
!c    write(*,*)'Back from calc_GAM_ee'
!c    write(*,*)'Calling calc_GAM_ep'
!c    call flshbf(6)
!     call calc_GAM_ep(nebf,npebf,npbf,ngam1,GAM_ep)
!c    write(*,*)'Back from calc_GAM_ep'
!c    write(*,*)'Calling xchf1_fock'
!c    call flshbf(6)

!     call xchf1_fock(nebf,nebf2,npbf,npbf2,nelec,
!    x                ngee,ngam1,ngam2,ngam3,ngam4,DE,DP,
!    x                GAM_ecore,GAM_pcore,GAM_ep,GAM_ee,
!    x                focke,fockp,Xfocke,Xfockp,
!    x                E_total,E_nuc,E_ecore,E_pcore,E_ep,
!    x                E_ee,E_gam1,E_gam2,E_gam3,E_gam4,
!    x                S_total,S_gam1,S_gam2)

!     E_NEOHF=E_ecore+E_pcore+E_ep+E_ee+E_nuc
!c    E_NEOHF=E_NEOHF/S_total
!c    E_check1=E_total-(E_ecore+E_ee)
!CCCC
!     if(NFCK.eq.1) then
!c    if(NFCK.eq.2) then
!        write(*,*)
!        write(*,*)'S_total=',S_total
!        write(*,*)'E_ecore=',E_ecore/S_total
!        write(*,*)'E_ee=',E_ee/S_total
!        write(*,*)'E_pcore=',E_pcore/S_total
!        write(*,*)'E_ep=',E_ep/S_total
!        write(*,*)'E_gam1=',E_gam1/S_total
!        write(*,*)'E_gam2=',E_gam2/S_total
!        write(*,*)'E_gam3=',E_gam3/S_total
!        write(*,*)'E_gam4=',E_gam4/S_total
!        write(*,*)'E_nuc=',E_nuc
!        write(*,*)'E_NEOHF=',E_NEOHF
!        write(*,*)'E_total=',E_total
!        write(*,*)
!     end if
!     EXCHF1=E_total
!CCCC

! Update either electronic or nuclear Fock matrix
!     if(NFCK.eq.1) then
!c       write(*,*)'UPDATE:  In xchf1_fock elec XC1 Fock is:'
!c       call print_my_fock(nebf,nebfLT,Xfocke)
!c       write(*,*) 
!        call Fock_Difference(nebf,focke,Xfocke,Dfocke)
!        call Pack_Fock(nebf,nebfLT,Dfocke,Pfocke)
!        call Add_Fock(nLT,Pfocke,GMSFCK)
!     else if(NFCK.eq.2) then
!c       write(*,*)'UPDATE:  In xchf1_fock Nuc XC1 Fock is:'
!c       call print_my_fock(npbf,npbfLT,Xfockp)
!c       write(*,*) 
!        call Fock_Difference(npbf,fockp,Xfockp,Dfockp)
!        call Pack_Fock(npbf,npbfLT,Dfockp,Pfockp)
!        call Add_Fock(nLT,Pfockp,GMSFCK)
!     end if


!     return
!     end

