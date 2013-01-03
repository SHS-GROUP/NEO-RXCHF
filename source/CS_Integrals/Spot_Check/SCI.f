C Spot Check Integral
C *    *     *
C======================================================================
      program spot_check
C======================================================================
      implicit none

      call UT_single_int_test

      end

C======================================================================
      subroutine UT_single_int_test

C This is a wrapper for unit testing neo-xc integrals.
C
C The type of integral tested (e.g. G2 G3 G4)
C is determined by the number of basis functions 
C given in the user's input file.  Flags are then 
C set accordingly:
C
C do_G2 :: 1e/1p
C do_G3 :: 2e/1p
C do_G4 :: 3e/1p
C
C  Only one integral class can be tested at a time.
C
C INPUT FILE FORMAT:
C
C gamma1
C gamma2
C Nbra
C I J K A x y z
C I J K A x y z
C Nket
C L M N B x y z
C L M N B x y z
C cmatx cmaty cmatz ==> coordinates of classical nucleus
C ZNUC ==> charge on classical nucleus
C nquad_coulomb (should always be >0)
C USE_RI
C USE_MD
C USE_AE
C USE_AC
C======================================================================
      implicit none

      integer ia,ii
      integer NQUAD_coul
      integer NBra
      integer NKet

      integer I(5),J(5),K(5)
      integer L(5),M(5),N(5)

      double precision A(5)
      double precision B(5)

      double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2
      double precision zero
      parameter (zero=0.0d+00)

      double precision Cmat(3)
      double precision Amat(5,3)
      double precision Bmat(5,3)

      logical use_RI
      logical orth_abs
      logical use_MD
      logical use_AE
      logical use_AC
      logical do_G2
      logical do_G3
      logical do_G4

C Initialize variables:
      GAMMA1=zero
      GAMMA2=zero
      ZNUC=zero

      Cmat(1)=zero
      Cmat(2)=zero
      Cmat(3)=zero

      do ia=1,5
         I(ia)=0
         J(ia)=0
         K(ia)=0
         A(ia)=zero
         L(ia)=0
         M(ia)=0
         N(ia)=0
         B(ia)=zero

         do ii=1,3
            Amat(ia,ii)=zero
            Bmat(ia,ii)=zero
         end do

      end do

      use_RI=.false.
      orth_abs=.false.
      use_MD=.false.
      use_AE=.false.
      use_AC=.false.
      do_G2=.false.
      do_G3=.false.
      do_G4=.false.

C Read input file
      open(unit=9,file='basis_definition.inp')
      read(9,*)gamma1
      read(9,*)gamma2

C Read Bra:
      read(9,*)Nbra
      do ia=1,Nbra
      read(9,*)I(ia),J(ia),K(ia),A(ia),Amat(ia,1),Amat(ia,2),Amat(ia,3)
      end do

C Read Ket:
      read(9,*)Nket

      if(Nbra.ne.Nket) then
         write(*,*)'Error:  Nbra must equal Nket'
         stop
      end if

      do ia=1,Nket
      read(9,*)L(ia),M(ia),N(ia),B(ia),Bmat(ia,1),Bmat(ia,2),Bmat(ia,3)
      end do

C Read classical nucleus info (for VPC and VEC type integrals)
      read(9,*) Cmat(1),Cmat(2),Cmat(3)
      read(9,*) ZNUC

C Control for desired integrals

      if(Nbra.eq.2) do_G2=.true.
      if(Nbra.eq.3) do_G3=.true.
      if(Nbra.eq.4) do_G4=.true.
c     if(Nbra.eq.5) do_G5=.true.

      read(9,*) NQUAD_coul
      read(9,*) use_RI
      read(9,*) orth_abs
      read(9,*) use_MD
      read(9,*) use_AE
      read(9,*) use_AC

C Execute integral calculation

      call integral_driver(use_RI,orth_abs,
     x                     use_MD,use_AC,use_AE,
     x                     do_G2,do_G3,do_G4,
     x                     NQUAD_coul,GAMMA1,GAMMA2,ZNUC,
     x                     Amat,Bmat,Cmat,
     x                     I,J,K,A,L,M,N,B)



      return
      end

C======================================================================
      subroutine integral_driver(use_RI,orth_abs,
     x                           use_MD,use_AC,use_AE,
     x                           do_G2,do_G3,do_G4,
     x                           NQUAD_coul,GAMMA1,GAMMA2,ZNUC,
     x                           Amat,Bmat,Cmat,
     x                           I,J,K,A,L,M,N,B)


C======================================================================
      implicit none 

      integer NQUAD_coul
      integer I(5),J(5),K(5)
      integer L(5),M(5),N(5)
      double precision A(5)
      double precision B(5)

      double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2

      double precision Cmat(3)
      double precision Amat(5,3)
      double precision Bmat(5,3)

      logical orth_abs
      logical use_AC
      logical use_AE
      logical use_RI
      logical use_MD
      logical do_G2
      logical do_G3
      logical do_G4

C Local variables
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Amat3(3)
      double precision Amat4(3)
      double precision Amat5(3)

      double precision Bmat1(3)
      double precision Bmat2(3)
      double precision Bmat3(3)
      double precision Bmat4(3)
      double precision Bmat5(3)

      double precision G2_xggs_AC
      double precision G2_xggs_MD
      double precision G2_xggVEC_AC
      double precision G2_xggVEC_MD
      double precision G2_xggVEC_RI
      double precision G2_xggVEP_AC
      double precision G2_xggVEP_MD
      double precision G2_xggVEP_RI

      double precision G3_xgsg_MD
      double precision G3_xgTE_MD
      double precision G3_xgTEg1_MD
      double precision G3_xgTEg2_MD
      double precision G3_xgTEg3_MD
      double precision G3_xgTPg_MD
      double precision G3_xgVEC_MD
      double precision G3_xgVPCg_MD
      double precision G3_xgVECg1_MD
      double precision G3_xgVECg2_MD
      double precision G3_xgVEP_MD
      double precision G3_xgVEE_MD
      double precision G3_xgVEPg1_MD
      double precision G3_xgVEPg2_MD
      double precision G3_xgVEEg1_MD
      double precision G3_xgVEEg2_MD

      double precision G3_xgVEC_RI
      double precision G3_xgVPCg_RI
      double precision G3_xgVECg1_RI
      double precision G3_xgVECg2_RI
      double precision G3_xgVEP_RI
      double precision G3_xgVEE_RI
      double precision G3_xgVEPg1_RI
      double precision G3_xgVEPg2_RI
      double precision G3_xgVEEg1_RI
      double precision G3_xgVEEg2_RI

      double precision G3_xgVEC_AE
      double precision G3_xgVPCg_AE
      double precision G3_xgVECg1_AE
      double precision G3_xgVECg2_AE
      double precision G3_xgVEP_AE
      double precision G3_xgVEE_AE
      double precision G3_xgVEPg1_AE
      double precision G3_xgVEPg2_AE
      double precision G3_xgVEEg1_AE

      double precision G3_xgsg_AC
      double precision G3_xgTE_AC
      double precision G3_xgTEg1_AC
      double precision G3_xgTEg2_AC
      double precision G3_xgTEg3_AC
      double precision G3_xgTPg_AC
      double precision G3_xgVEC_AC
      double precision G3_xgVPCg_AC
      double precision G3_xgVECg1_AC
      double precision G3_xgVECg2_AC
      double precision G3_xgVEP_AC
      double precision G3_xgVEE_AC
      double precision G3_xgVEPg1_AC
      double precision G3_xgVEPg2_AC
      double precision G3_xgVEEg1_AC
      double precision G3_xgVEEg2_AC

      double precision G4_xgVEPg_MD
      double precision G4_xgVEPg_RI
      double precision G4_xgVEPg_AC
      double precision G4_xgVEPg_AE
      double precision G4_xgVEEg2_MD
      double precision G4_xgVEEg2_RI
      double precision G4_xgVEEg2_AC
      double precision G4_xgVEEg2_AE
      double precision xgVepg_SONLY1
      double precision xgVepg_SONLY2
      double precision xgVeeg_SONLY1
      double precision xgVeeg_SONLY2
      
      write(*,*)
      write(*,*)'**************************************'
      write(*,*)'    NEO-XC Integral Test Driver       '
      write(*,*)'**************************************'
      write(*,*)

C If RI has been requested, orthonormalize
C the auxilliary basis set (ABS)
      if(use_RI) call OrthoABS

      Amat1(1)=Amat(1,1)
      Amat1(2)=Amat(1,2)
      Amat1(3)=Amat(1,3)

      Amat2(1)=Amat(2,1)
      Amat2(2)=Amat(2,2)
      Amat2(3)=Amat(2,3)

      Amat3(1)=Amat(3,1)
      Amat3(2)=Amat(3,2)
      Amat3(3)=Amat(3,3)

      Amat4(1)=Amat(4,1)
      Amat4(2)=Amat(4,2)
      Amat4(3)=Amat(4,3)

      Amat5(1)=Amat(5,1)
      Amat5(2)=Amat(5,2)
      Amat5(3)=Amat(5,3)

      Bmat1(1)=Bmat(1,1)
      Bmat1(2)=Bmat(1,2)
      Bmat1(3)=Bmat(1,3)

      Bmat2(1)=Bmat(2,1)
      Bmat2(2)=Bmat(2,2)
      Bmat2(3)=Bmat(2,3)

      Bmat3(1)=Bmat(3,1)
      Bmat3(2)=Bmat(3,2)
      Bmat3(3)=Bmat(3,3)

      Bmat4(1)=Bmat(4,1)
      Bmat4(2)=Bmat(4,2)
      Bmat4(3)=Bmat(4,3)

      Bmat5(1)=Bmat(5,1)
      Bmat5(2)=Bmat(5,2)
      Bmat5(3)=Bmat(5,3)

     
      if(do_G2) then
         write(*,*)
         write(*,*)'**************************************'
         write(*,*)'           G2 Integrals             '
         write(*,*)'**************************************'
         write(*,*)
         if(use_AC) then
            call G2_AC_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                        I(2),J(2),K(2),A(2),Amat2,
     x                        L(1),M(1),N(1),B(1),Bmat1,
     x                        L(2),M(2),N(2),B(2),Bmat2,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        G2_xggs_AC,G2_xggVEC_AC,G2_xggVEP_AC)
            write(*,*)
            write(*,*)'---------- AC VALUES ------------'
            write(*,*)'XGGS  =',G2_xggs_AC
            write(*,*)'XGGVEC=',G2_xggVEC_AC
            write(*,*)'XGGVEP=',G2_xggVEP_AC
         end if
         if(use_MD) then
            call G2_MD_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                        I(2),J(2),K(2),A(2),Amat2,
     x                        L(1),M(1),N(1),B(1),Bmat1,
     x                        L(2),M(2),N(2),B(2),Bmat2,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        G2_xggs_MD,G2_xggVEC_MD,G2_xggVEP_MD)
            write(*,*)
            write(*,*)'---------- MD VALUES ------------'
            write(*,*)'XGGS  =',G2_xggs_MD
            write(*,*)'XGGVEC=',G2_xggVEC_MD
            write(*,*)'XGGVEP=',G2_xggVEP_MD
         end if
         if(use_AE) then
            write(*,*)
            write(*,*)'---------- AE VALUES ------------'
            write(*,*)'Sorry, Aux Expansion approach is not'
            write(*,*)'implemented for G2 integrals.'
            write(*,*)
         end if
         if(use_RI) then
            call G2_RI_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                        I(2),J(2),K(2),A(2),Amat2,
     x                        L(1),M(1),N(1),B(1),Bmat1,
     x                        L(2),M(2),N(2),B(2),Bmat2,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        G2_xggVEC_RI,G2_xggVEP_RI)
            write(*,*)
            write(*,*)'---------- RI VALUES ------------'
            write(*,*)'XGGVEC=',G2_xggVEC_RI
            write(*,*)'XGGVEP=',G2_xggVEP_RI
         end if
      end if

      if(do_G3) then
         write(*,*)
         write(*,*)'**************************************'
         write(*,*)'           G3 Integrals             '
         write(*,*)'**************************************'
         write(*,*)
         if(use_AC) then
            call G3_AC_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                        I(2),J(2),K(2),A(2),Amat2,
     x                        I(3),J(3),K(3),A(3),Amat3,
     x                        L(1),M(1),N(1),B(1),Bmat1,
     x                        L(2),M(2),N(2),B(2),Bmat2,
     x                        L(3),M(3),N(3),B(3),Bmat3,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        G3_xgsg_AC,G3_xgTE_AC,G3_xgTEg1_AC,
     x                        G3_xgTEg3_AC,
     x                        G3_xgTEg2_AC,G3_xgTPg_AC,
     x                        G3_xgVEC_AC,G3_xgVPCg_AC,
     x                        G3_xgVECg1_AC,G3_xgVECg2_AC,G3_xgVEP_AC,
     x                        G3_xgVEE_AC,G3_xgVEPg1_AC,G3_xgVEPg2_AC,
     x                        G3_xgVEEg1_AC,G3_xgVEEg2_AC)

            write(*,*)
            write(*,*)'---------- AC VALUES ------------'
            write(*,*)'xgSg   =',G3_xgsg_AC
            write(*,*)'xgTE   =',G3_xgTE_AC
            write(*,*)'xgTEg1 =',G3_xgTEg1_AC
            write(*,*)'xgTEg2 =',G3_xgTEg2_AC
            write(*,*)'xgTEg3 =',G3_xgTEg3_AC
            write(*,*)'xgTPg  =',G3_xgTPg_AC
            write(*,*)'xgVEC  =',G3_xgVEC_AC
            write(*,*)'xgVPCg =',G3_xgVPCg_AC
            write(*,*)'xgVECg1=',G3_xgVECg1_AC
            write(*,*)'xgVECg2=',G3_xgVECg2_AC
            write(*,*)'xgVEP  =',G3_xgVEP_AC
            write(*,*)'xgVEE  =',G3_xgVEE_AC
            write(*,*)'xgVEPg1=',G3_xgVEPg1_AC
            write(*,*)'xgVEPg2=',G3_xgVEPg2_AC
            write(*,*)'xgVEEg1=',G3_xgVEEg1_AC
            write(*,*)'xgVEEg2=',G3_xgVEEg2_AC
         end if
         if(use_MD) then
            call G3_MD_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                        I(2),J(2),K(2),A(2),Amat2,
     x                        I(3),J(3),K(3),A(3),Amat3,
     x                        L(1),M(1),N(1),B(1),Bmat1,
     x                        L(2),M(2),N(2),B(2),Bmat2,
     x                        L(3),M(3),N(3),B(3),Bmat3,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        G3_xgsg_MD,G3_xgTE_MD,G3_xgTEg1_MD,
     x                        G3_xgTEg3_MD,
     x                        G3_xgTEg2_MD,G3_xgTPg_MD,
     x                        G3_xgVEC_MD,G3_xgVPCg_MD,
     x                        G3_xgVECg1_MD,G3_xgVECg2_MD,G3_xgVEP_MD,
     x                        G3_xgVEE_MD,G3_xgVEPg1_MD,G3_xgVEPg2_MD,
     x                        G3_xgVEEg1_MD,G3_xgVEEg2_MD)


            write(*,*)
            write(*,*)'---------- MD VALUES ------------'
            write(*,*)'xgSg   =',G3_xgsg_MD
            write(*,*)'xgTE   =',G3_xgTE_MD
            write(*,*)'xgTEg1 =',G3_xgTEg1_MD
            write(*,*)'xgTEg2 =',G3_xgTEg2_MD
            write(*,*)'xgTEg3 =',G3_xgTEg3_MD
            write(*,*)'xgTPg  =',G3_xgTPg_MD
            write(*,*)'xgVEC  =',G3_xgVEC_MD
            write(*,*)'xgVPCg =',G3_xgVPCg_MD
            write(*,*)'xgVECg1=',G3_xgVECg1_MD
            write(*,*)'xgVECg2=',G3_xgVECg2_MD
            write(*,*)'xgVEP  =',G3_xgVEP_MD
            write(*,*)'xgVEE  =',G3_xgVEE_MD
            write(*,*)'xgVEPg1=',G3_xgVEPg1_MD
            write(*,*)'xgVEPg2=',G3_xgVEPg2_MD
            write(*,*)'xgVEEg1=',G3_xgVEEg1_MD
            write(*,*)'xgVEEg2=',G3_xgVEEg2_MD
         end if

         if(use_AE) then
            call G3_AE_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                        I(2),J(2),K(2),A(2),Amat2,
     x                        I(3),J(3),K(3),A(3),Amat3,
     x                        L(1),M(1),N(1),B(1),Bmat1,
     x                        L(2),M(2),N(2),B(2),Bmat2,
     x                        L(3),M(3),N(3),B(3),Bmat3,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        G3_xgVEC_AE,G3_xgVPCg_AE,
     x                        G3_xgVECg1_AE,G3_xgVECg2_AE,G3_xgVEP_AE,
     x                        G3_xgVEE_AE,G3_xgVEPg1_AE,G3_xgVEPg2_AE,
     x                        G3_xgVEEg1_AE)
            write(*,*)
            write(*,*)'---------- AE VALUES ------------'
            write(*,*)'xgVEC  =',G3_xgVEC_AE
            write(*,*)'xgVPCg =',G3_xgVPCg_AE
            write(*,*)'xgVECg1=',G3_xgVECg1_AE
            write(*,*)'xgVECg2=',G3_xgVECg2_AE
            write(*,*)'xgVEP  =',G3_xgVEP_AE
            write(*,*)'xgVEE  =',G3_xgVEE_AE
            write(*,*)'xgVEPg1=',G3_xgVEPg1_AE
            write(*,*)'xgVEPg2=',G3_xgVEPg2_AE
            write(*,*)'xgVEEg1=',G3_xgVEEg1_AE

         end if
         if(use_RI) then
            call G3_RI_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                        I(2),J(2),K(2),A(2),Amat2,
     x                        I(3),J(3),K(3),A(3),Amat3,
     x                        L(1),M(1),N(1),B(1),Bmat1,
     x                        L(2),M(2),N(2),B(2),Bmat2,
     x                        L(3),M(3),N(3),B(3),Bmat3,
     x                        gamma1,gamma2,
     x                        ZNUC,Cmat,
     x                        G3_xgVEC_RI,G3_xgVPCg_RI,
     x                        G3_xgVECg1_RI,G3_xgVECg2_RI,G3_xgVEP_RI,
     x                        G3_xgVEE_RI,G3_xgVEPg1_RI,G3_xgVEPg2_RI,
     x                        G3_xgVEEg1_RI,G3_xgVEEg2_RI)

            write(*,*)
            write(*,*)'---------- RI VALUES ------------'
            write(*,*)'xgVEC  =',G3_xgVEC_RI
            write(*,*)'xgVPCg =',G3_xgVPCg_RI
            write(*,*)'xgVECg1=',G3_xgVECg1_RI
            write(*,*)'xgVECg2=',G3_xgVECg2_RI
            write(*,*)'xgVEP  =',G3_xgVEP_RI
            write(*,*)'xgVEE  =',G3_xgVEE_RI
            write(*,*)'xgVEPg1=',G3_xgVEPg1_RI
            write(*,*)'xgVEPg2=',G3_xgVEPg2_RI
            write(*,*)'xgVEEg1=',G3_xgVEEg1_RI
            write(*,*)'xgVEEg2=',G3_xgVEEg2_RI
         end if
      end if

      if(do_G4) then
         write(*,*)
         write(*,*)'**************************************'
         write(*,*)'           G4 Integrals             '
         write(*,*)'**************************************'
         write(*,*)
         if(use_AC) then
            call G4_AC_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                        I(2),J(2),K(2),A(2),Amat2,
     x                        I(3),J(3),K(3),A(3),Amat3,
     x                        I(4),J(4),K(4),A(4),Amat4,
     x                        L(1),M(1),N(1),B(1),Bmat1,
     x                        L(2),M(2),N(2),B(2),Bmat2,
     x                        L(3),M(3),N(3),B(3),Bmat3,
     x                        L(4),M(4),N(4),B(4),Bmat4,
     x                        gamma1,gamma2,
     x                        G4_xgVEPg_AC,G4_xgVEEg2_AC)

            write(*,*)
            write(*,*)'---------- AC VALUES ------------'
            write(*,*)'xgVEPg =',G4_xgVEPg_AC
            write(*,*)'xgVEEg2=',G4_xgVEEg2_AC
         end if
         if(use_MD) then
            call G4_MD_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                        I(2),J(2),K(2),A(2),Amat2,
     x                        I(3),J(3),K(3),A(3),Amat3,
     x                        I(4),J(4),K(4),A(4),Amat4,
     x                        L(1),M(1),N(1),B(1),Bmat1,
     x                        L(2),M(2),N(2),B(2),Bmat2,
     x                        L(3),M(3),N(3),B(3),Bmat3,
     x                        L(4),M(4),N(4),B(4),Bmat4,
     x                        gamma1,gamma2,
     x                        G4_xgVEPg_MD,G4_xgVEEg2_MD)

            call G4_sonly_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                           I(2),J(2),K(2),A(2),Amat2,
     x                           I(3),J(3),K(3),A(3),Amat3,
     x                           I(4),J(4),K(4),A(4),Amat4,
     x                           L(1),M(1),N(1),B(1),Bmat1,
     x                           L(2),M(2),N(2),B(2),Bmat2,
     x                           L(3),M(3),N(3),B(3),Bmat3,
     x                           L(4),M(4),N(4),B(4),Bmat4,
     x                           gamma1,gamma2,
     x                           xgVeeg_SONLY1,
     x                           xgVeeg_SONLY2,
     x                           xgVepg_SONLY1,
     x                           xgVepg_SONLY2)

            write(*,*)
            write(*,*)'---------- MD VALUES ------------'
            write(*,*)'xgVEPg =',G4_xgVEPg_MD
            write(*,*)'xgVEEg2=',G4_xgVEEg2_MD
            write(*,*)'xgVEPg  SO1 =',xgVepg_SONLY1
            write(*,*)'xgVEPg  SO2 =',xgVepg_SONLY2
            write(*,*)'xgVEEg2 SO1 =',xgVeeg_SONLY1
            write(*,*)'xgVEEg2 SO2 =',xgVeeg_SONLY2

         end if

         if(use_AE) then
            call G4_AE_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                        I(2),J(2),K(2),A(2),Amat2,
     x                        I(3),J(3),K(3),A(3),Amat3,
     x                        I(4),J(4),K(4),A(4),Amat4,
     x                        L(1),M(1),N(1),B(1),Bmat1,
     x                        L(2),M(2),N(2),B(2),Bmat2,
     x                        L(3),M(3),N(3),B(3),Bmat3,
     x                        L(4),M(4),N(4),B(4),Bmat4,
     x                        gamma1,gamma2,
     x                        G4_xgVEEg2_AE,G4_xgVEPg_AE)
            write(*,*)
            write(*,*)'---------- AE VALUES ------------'
            write(*,*)'xgVEPg =',G4_xgVEPg_AE
            write(*,*)'xgVEEg2=',G4_xgVEEg2_AE

         end if
         if(use_RI) then
            call G4_RI_DRIVER(I(1),J(1),K(1),A(1),Amat1,
     x                        I(2),J(2),K(2),A(2),Amat2,
     x                        I(3),J(3),K(3),A(3),Amat3,
     x                        I(4),J(4),K(4),A(4),Amat4,
     x                        L(1),M(1),N(1),B(1),Bmat1,
     x                        L(2),M(2),N(2),B(2),Bmat2,
     x                        L(3),M(3),N(3),B(3),Bmat3,
     x                        L(4),M(4),N(4),B(4),Bmat4,
     x                        gamma1,gamma2,
     x                        G4_xgVEPg_RI,G4_xgVEEg2_RI)

            write(*,*)
            write(*,*)'---------- RI VALUES ------------'
            write(*,*)'xgVEPg =',G4_xgVEPg_RI
            write(*,*)'xgVEEg2=',G4_xgVEEg2_RI
         end if
      end if



      RETURN
      END

