C RI_G3_xgVEE.f  
C RI_G3_xgVEEg1.f
C RI_G3_xgVEEg2.f
C RI_G3_xgVEP.f  
C RI_G3_xgVEPg1.f
C RI_G3_xgVEPg2.f
C RI_G3_xgVECg1.f
C RI_G3_xgVPCg.f

C======================================================================
      subroutine RI_G3_xgVPCg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        Cmat,ZNUC,
     *                        xgVPCg)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3

      double precision  A1
      double precision  A2
      double precision  A3
      double precision  B1
      double precision  B2
      double precision  B3

      double precision  ZNUC
      double precision  Cmat(3)

      double precision  gamA12
      double precision  gamA13
      double precision  gamA23
      double precision  gamB12
      double precision  gamB13
      double precision  gamB23

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Amat3(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)
      double precision  Bmat3(3)

C Variables returned
      double precision  xgVPCg

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Mab
      double precision  xvpc
      double precision  xggs
      double precision  abs_norm_a
      double precision  abs_norm_b
      double precision  ans
c     double precision  gamma1
c     double precision  gamma2

      integer  NauxBF
      integer  i
      integer  a
      integer  b
      integer  ia
      integer  LX(1000),MX(1000),NX(1000)

      double precision BX(1000)
      double precision BmatX(3,1000)
      double precision BmatX1(3)
      double precision BmatX2(3)

cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'USING RI'
      do i=1,1000
         LX(i)=0
         MX(i)=0
         NX(i)=0
         BX(i)=zero
         BmatX(1,i)=zero
         BmatX(2,i)=zero
         BmatX(3,i)=zero
      end do

c Read RI auxilliary basis input file
      open(unit=912,file='aux_p_basis.inp',status='old')
      read(912,*) NauxBF

      do i=1,NauxBF
      read(912,*)LX(i),MX(i),NX(i),
     x BX(i),BmatX(1,i),BmatX(2,i),BmatX(3,i)
      end do

C  Open ABS normalization constant file for p ABS
      open(914,file='ABSNORMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C  Open transformation matrix product file for p ABS
      open(920,file='ABSMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)


      xgVPCg=zero
      ans=zero

      do b=1,NauxBF

         read(914,REC=b) abs_norm_b
         BmatX1(1)=Bmatx(1,b)
         BmatX1(2)=Bmatx(2,b)
         BmatX1(3)=Bmatx(3,b)
 
c  Calculate Overlap 
         call G3ovlap(I1,J1,K1,A1,Amat1,
     *                I2,J2,K2,A2,Amat2,
     *                I3,J3,K3,A3,Amat3,
     *                L1,M1,N1,B1,Bmat1,
     *                L2,M2,N2,B2,Bmat2,
     x                LX(b),MX(b),NX(b),BX(b),BmatX1,
     *                gamA12,gamA13,gamA23,
     *                gamB12,gamB13,gamB23,
     *                xggs)

         do a=1,NauxBF

            read(914,REC=a) abs_norm_a
            BmatX2(1)=Bmatx(1,a)
            BmatX2(2)=Bmatx(2,a)
            BmatX2(3)=Bmatx(3,a)

C  Read in Mab
            call pack_2D(NauxBF,a,b,ia)
            read(920,REC=ia) Mab 

c  Calculate VeC 
            call gfvec(L3,M3,N3,B3,Bmat3,
     x                 LX(a),MX(a),NX(a),BX(a),BmatX2,
     x                 Cmat,xvpc)

            ans=xvpc*xggs*Mab*abs_norm_a*abs_norm_b

            xgVPCg=ans+xgVPCg

         end do

      end do

         xgVPCg=ZNUC*xgVPCg

      close(920)
      close(914)
      close(912)
      

      return
      end

c        call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
c    x                      L3,M3,N3,B3,Bmat3,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      LX(k3B),MX(k3B),NX(k3B),BX(k3B),BmatX1,
c    x                      gamma1,gamma2,xggsB)
C
c        call pgiovlap(I2,J2,K2,A2,Amat2,
c    x                 L3,M3,N3,B3,Bmat3,
c    x                 L2,M2,N2,B2,Bmat2,
c    x                 LX(a),MX(a),NX(a),BX(a),BmatX1,
c    x                 gamma1,gamma2,xggsB)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'CHECK'
c     write(*,*)'K1=',K1
c     write(*,*)'K2=',K2
c     write(*,*)'K3=',K3
c     write(*,*)'N1=',N1
c     write(*,*)'N2=',N2
c     write(*,*)'N3=',N3
c     write(*,*)'xvec=',xvec
c     write(*,*)'xggs=',xggs
c     xgVEC=ZNUC*xvec*xggs


C======================================================================
      subroutine RI_G3_xgVECg1(I1,J1,K1,A1,Amat1,
     *                         I2,J2,K2,A2,Amat2,
     *                         I3,J3,K3,A3,Amat3,
     *                         L1,M1,N1,B1,Bmat1,
     *                         L2,M2,N2,B2,Bmat2,
     *                         L3,M3,N3,B3,Bmat3,
     *                         gamA12,gamA13,gamA23,
     *                         gamB12,gamB13,gamB23,
     *                         Cmat,ZNUC,
     *                         xgVECg1)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3

      double precision  A1
      double precision  A2
      double precision  A3
      double precision  B1
      double precision  B2
      double precision  B3

      double precision  ZNUC
      double precision  Cmat(3)

      double precision  gamA12
      double precision  gamA13
      double precision  gamA23
      double precision  gamB12
      double precision  gamB13
      double precision  gamB23

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Amat3(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)
      double precision  Bmat3(3)

C Variables returned
      double precision  xgVECg1

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Mab
      double precision  xvec
      double precision  xggs
      double precision  abs_norm_a
      double precision  abs_norm_b
      double precision  ans
c     double precision  gamma1
c     double precision  gamma2

      integer  NauxBF
      integer  i
      integer  a
      integer  b
      integer  ia
      integer  LX(1000),MX(1000),NX(1000)

      double precision BX(1000)
      double precision BmatX(3,1000)
      double precision BmatX1(3)
      double precision BmatX2(3)

cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'USING RI'
      do i=1,1000
         LX(i)=0
         MX(i)=0
         NX(i)=0
         BX(i)=zero
         BmatX(1,i)=zero
         BmatX(2,i)=zero
         BmatX(3,i)=zero
      end do

c Read RI auxilliary basis input file
      open(unit=911,file='aux_e_basis.inp',status='old')
      read(911,*) NauxBF

      do i=1,NauxBF
      read(911,*)LX(i),MX(i),NX(i),
     x BX(i),BmatX(1,i),BmatX(2,i),BmatX(3,i)
      end do

C  Open ABS normalization constant file for p ABS
      open(913,file='ABSNORME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C  Open transformation matrix product file for p ABS
      open(919,file='ABSME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)


      xgVECg1=zero
      ans=zero

      do b=1,NauxBF

         read(913,REC=b) abs_norm_b
         BmatX1(1)=Bmatx(1,b)
         BmatX1(2)=Bmatx(2,b)
         BmatX1(3)=Bmatx(3,b)
 
c  Calculate Overlap 
         call G3ovlap(I1,J1,K1,A1,Amat1,
     *                I2,J2,K2,A2,Amat2,
     *                I3,J3,K3,A3,Amat3,
     x                LX(b),MX(b),NX(b),BX(b),BmatX1,
     *                L2,M2,N2,B2,Bmat2,
     *                L3,M3,N3,B3,Bmat3,
     *                gamA12,gamA13,gamA23,
     *                gamB12,gamB13,gamB23,
     *                xggs)

         do a=1,NauxBF

            read(913,REC=a) abs_norm_a
            BmatX2(1)=Bmatx(1,a)
            BmatX2(2)=Bmatx(2,a)
            BmatX2(3)=Bmatx(3,a)

C  Read in Mab
            call pack_2D(NauxBF,a,b,ia)
            read(919,REC=ia) Mab 

c  Calculate VeC 
            call gfvec(L1,M1,N1,B1,Bmat1,
     x                 LX(a),MX(a),NX(a),BX(a),BmatX2,
     x                 Cmat,xvec)

            ans=xvec*xggs*Mab*abs_norm_a*abs_norm_b

            xgVECg1=ans+xgVECg1

         end do

      end do

         xgVECg1=-1.0d+00*ZNUC*xgVECg1

      close(919)
      close(913)
      close(911)
      

      return
      end

c        call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
c    x                      L3,M3,N3,B3,Bmat3,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      LX(k3B),MX(k3B),NX(k3B),BX(k3B),BmatX1,
c    x                      gamma1,gamma2,xggsB)
C
c        call pgiovlap(I2,J2,K2,A2,Amat2,
c    x                 L3,M3,N3,B3,Bmat3,
c    x                 L2,M2,N2,B2,Bmat2,
c    x                 LX(a),MX(a),NX(a),BX(a),BmatX1,
c    x                 gamma1,gamma2,xggsB)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'CHECK'
c     write(*,*)'K1=',K1
c     write(*,*)'K2=',K2
c     write(*,*)'K3=',K3
c     write(*,*)'N1=',N1
c     write(*,*)'N2=',N2
c     write(*,*)'N3=',N3
c     write(*,*)'xvec=',xvec
c     write(*,*)'xggs=',xggs
c     xgVEC=ZNUC*xvec*xggs


C======================================================================
      subroutine RI_G3_xgVEPg2(I1,J1,K1,A1,Amat1,
     *                         I2,J2,K2,A2,Amat2,
     *                         I3,J3,K3,A3,Amat3,
     *                         L1,M1,N1,B1,Bmat1,
     *                         L2,M2,N2,B2,Bmat2,
     *                         L3,M3,N3,B3,Bmat3,
     *                         gamA12,gamA13,gamA23,
     *                         gamB12,gamB13,gamB23,
     *                         xgVEPg2)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3

      double precision  A1
      double precision  A2
      double precision  A3
      double precision  B1
      double precision  B2
      double precision  B3

c     double precision  ZNUC
c     double precision  Cmat(3)

      double precision  gamA12
      double precision  gamA13
      double precision  gamA23
      double precision  gamB12
      double precision  gamB13
      double precision  gamB23

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Amat3(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)
      double precision  Bmat3(3)

C Variables returned
      double precision  xgVEPg2

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Mab
      double precision  xvep
      double precision  xggs
      double precision  abs_norm_a
      double precision  abs_norm_b
      double precision  ans
      double precision  gamma1
      double precision  gamma2

      integer  NauxBF
      integer  i
      integer  a
      integer  b
      integer  ia
      integer  LX(1000),MX(1000),NX(1000)

      double precision BX(1000)
      double precision BmatX(3,1000)
      double precision BmatX1(3)
      double precision BmatX2(3)

cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'USING RI'
      do i=1,1000
         LX(i)=0
         MX(i)=0
         NX(i)=0
         BX(i)=zero
         BmatX(1,i)=zero
         BmatX(2,i)=zero
         BmatX(3,i)=zero
      end do

c Read RI auxilliary basis input file
      open(unit=912,file='aux_p_basis.inp',status='old')
      read(912,*) NauxBF

      do i=1,NauxBF
      read(912,*)LX(i),MX(i),NX(i),
     x BX(i),BmatX(1,i),BmatX(2,i),BmatX(3,i)
      end do

C  Open ABS normalization constant file for p ABS
      open(914,file='ABSNORMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C  Open transformation matrix product file for p ABS
      open(920,file='ABSMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)


      xgVEPg2=zero
      ans=zero

      do b=1,NauxBF

         read(914,REC=b) abs_norm_b
         BmatX1(1)=Bmatx(1,b)
         BmatX1(2)=Bmatx(2,b)
         BmatX1(3)=Bmatx(3,b)
 
         gamma1=gamA23
         gamma2=gamB23

c  Calculate Overlap 
c        call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
c    x                      I3,J3,K3,A3,Amat3,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      LX(b),MX(b),NX(b),BX(b),BmatX1,
c    x                      gamma1,gamma2,xggs)
C
         call pgiovlap(I2,J2,K2,A2,Amat2,
     x                 I3,J3,K3,A3,Amat3,
     x                 L2,M2,N2,B2,Bmat2,
     x                 LX(b),MX(b),NX(b),BX(b),BmatX1,
     x                 gamma1,gamma2,xggs)


         do a=1,NauxBF

            read(914,REC=a) abs_norm_a
            BmatX2(1)=Bmatx(1,a)
            BmatX2(2)=Bmatx(2,a)
            BmatX2(3)=Bmatx(3,a)

C  Read in Mab
            call pack_2D(NauxBF,a,b,ia)
            read(920,REC=ia) Mab 

c  Calculate Vep 
            call gfvee(I1,J1,K1,A1,Amat1,
     x                 LX(a),MX(a),NX(a),BX(a),BmatX2,
     x                 L1,M1,N1,B1,Bmat1,
     x                 L3,M3,N3,B3,Bmat3,
     x                 xvep)

            ans=xvep*xggs*Mab*abs_norm_a*abs_norm_b

            xgVEPg2=ans+xgVEPg2

         end do

      end do

      xgVEPg2=-1.0d+00*xgVEPg2

      close(920)
      close(914)
      close(912)
      

      return
      end

c        call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
c    x                      L3,M3,N3,B3,Bmat3,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      LX(k3B),MX(k3B),NX(k3B),BX(k3B),BmatX1,
c    x                      gamma1,gamma2,xggsB)
C
c        call pgiovlap(I2,J2,K2,A2,Amat2,
c    x                 L3,M3,N3,B3,Bmat3,
c    x                 L2,M2,N2,B2,Bmat2,
c    x                 LX(a),MX(a),NX(a),BX(a),BmatX1,
c    x                 gamma1,gamma2,xggsB)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'CHECK'
c     write(*,*)'K1=',K1
c     write(*,*)'K2=',K2
c     write(*,*)'K3=',K3
c     write(*,*)'N1=',N1
c     write(*,*)'N2=',N2
c     write(*,*)'N3=',N3
c     write(*,*)'xvec=',xvec
c     write(*,*)'xggs=',xggs
c     xgVEC=ZNUC*xvec*xggs


C======================================================================
      subroutine RI_G3_xgVEPg1(I1,J1,K1,A1,Amat1,
     *                         I2,J2,K2,A2,Amat2,
     *                         I3,J3,K3,A3,Amat3,
     *                         L1,M1,N1,B1,Bmat1,
     *                         L2,M2,N2,B2,Bmat2,
     *                         L3,M3,N3,B3,Bmat3,
     *                         gamA12,gamA13,gamA23,
     *                         gamB12,gamB13,gamB23,
     *                         xgVEPg1)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3

      double precision  A1
      double precision  A2
      double precision  A3
      double precision  B1
      double precision  B2
      double precision  B3

      double precision  gamA12
      double precision  gamA13
      double precision  gamA23
      double precision  gamB12
      double precision  gamB13
      double precision  gamB23

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Amat3(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)
      double precision  Bmat3(3)

C Variables returned
      double precision  xgVEPg1

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Ma1b1
      double precision  Ma3b3
      double precision  xvep
      double precision  xggs
      double precision  abs_norm_a1
      double precision  abs_norm_a3
      double precision  abs_norm_b1
      double precision  abs_norm_b3
      double precision  ans

      integer  i,ia
      integer  ia1,ib1
      integer  ia3,ib3

      integer  NauxBFE
      integer  NauxBFP

      integer  LXE(1000),MXE(1000),NXE(1000)
      integer  LXP(1000),MXP(1000),NXP(1000)

      double precision BXE(1000)
      double precision BXP(1000)

      double precision BmatXE(3,1000)
      double precision BmatXP(3,1000)

      double precision BmatXa1(3)
      double precision BmatXb1(3)

      double precision BmatXa3(3)
      double precision BmatXb3(3)

cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'USING RI'
      do i=1,1000
         LXE(i)=0
         MXE(i)=0
         NXE(i)=0
         BXE(i)=zero
         BmatXE(1,i)=zero
         BmatXE(2,i)=zero
         BmatXE(3,i)=zero

         LXP(i)=0
         MXP(i)=0
         NXP(i)=0
         BXP(i)=zero
         BmatXP(1,i)=zero
         BmatXP(2,i)=zero
         BmatXP(3,i)=zero
      end do

C>>>>>>>>>> READ ABS for E <<<<<<<<<<
      open(unit=911,file='aux_e_basis.inp',status='old')
      read(911,*) NauxBFE

      do i=1,NauxBFE
      read(911,*)LXE(i),MXE(i),NXE(i),
     x BXE(i),BmatXE(1,i),BmatXE(2,i),BmatXE(3,i)
      end do
C  Open ABS normalization constant file for e ABS
      open(913,file='ABSNORME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
C  Open transformation matrix product file for e ABS
      open(919,file='ABSME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C>>>>>>>>>> READ ABS for P <<<<<<<<<<
      open(unit=912,file='aux_p_basis.inp',status='old')
      read(912,*) NauxBFP

      do i=1,NauxBFP
      read(912,*)LXP(i),MXP(i),NXP(i),
     x BXP(i),BmatXP(1,i),BmatXP(2,i),BmatXP(3,i)
      end do
C  Open ABS normalization constant file for p ABS
      open(914,file='ABSNORMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
C  Open transformation matrix product file for p ABS
      open(920,file='ABSMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)


      xgVEPg1=zero
      ans=zero

      do ib1=1,NauxBFE

                read(913,REC=ib1) abs_norm_b1
         BmatXb1(1)=BmatXE(1,ib1)
         BmatXb1(2)=BmatXE(2,ib1)
         BmatXb1(3)=BmatXE(3,ib1)
 
         do ib3=1,NauxBFP

                   read(914,REC=ib3) abs_norm_b3
            BmatXb3(1)=BmatXP(1,ib3)
            BmatXb3(2)=BmatXP(2,ib3)
            BmatXb3(3)=BmatXP(3,ib3)
 

c  Calculate Overlap 
            call G3ovlap(LXE(ib1),MXE(ib1),NXE(ib1),BXE(ib1),BmatXb1,
     x                   I2,J2,K2,A2,Amat2,
     x                   LXP(ib3),MXP(ib3),NXP(ib3),BXP(ib3),BmatXb3,
     x                   L1,M1,N1,B1,Bmat1,
     x                   L2,M2,N2,B2,Bmat2,
     x                   L3,M3,N3,B3,Bmat3,
     x                   gamA12,gamA13,gamA23,
     x                   gamB12,gamB13,gamB23,
     x                   xggs)

            do ia1=1,NauxBFE

                      read(913,REC=ia1) abs_norm_a1
               BmatXa1(1)=BmatXE(1,ia1)
               BmatXa1(2)=BmatXE(2,ia1)
               BmatXa1(3)=BmatXE(3,ia1)
 
               do ia3=1,NauxBFP

                         read(914,REC=ia3) abs_norm_a3
                  BmatXa3(1)=BmatXP(1,ia3)
                  BmatXa3(2)=BmatXP(2,ia3)
                  BmatXa3(3)=BmatXP(3,ia3)

C  Read in Mab
                  call pack_2D(NauxBFE,ia1,ib1,ia)
                  read(919,REC=ia) Ma1b1 
                  call pack_2D(NauxBFP,ia3,ib3,ia)
                  read(920,REC=ia) Ma3b3 

c  Calculate Vep 
                  call gfvee(I1,J1,K1,A1,Amat1,
     x                       I3,J3,K3,A3,Amat3,
     x                      LXE(ia1),MXE(ia1),NXE(ia1),BXE(ia1),BmatXa1,
     x                      LXP(ia3),MXP(ia3),NXP(ia3),BXP(ia3),BmatXa3,
     x                       xvep)

                  ans=abs_norm_a1*abs_norm_b1*Ma1b1*
     x                abs_norm_a3*abs_norm_b3*Ma3b3*
     x                xvep*xggs

                  xgVEPg1=ans+xgVEPg1

               end do
            end do
         end do
      end do

      xgVEPg1=-1.0d+00*xgVEPg1

      close(919)
      close(913)
      close(911)
      
      close(920)
      close(914)
      close(912)

      return
      end

c        call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
c    x                      L3,M3,N3,B3,Bmat3,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      LX(k3B),MX(k3B),NX(k3B),BX(k3B),BmatX1,
c    x                      gamma1,gamma2,xggsB)
C
c        call pgiovlap(I2,J2,K2,A2,Amat2,
c    x                 L3,M3,N3,B3,Bmat3,
c    x                 L2,M2,N2,B2,Bmat2,
c    x                 LX(a),MX(a),NX(a),BX(a),BmatX1,
c    x                 gamma1,gamma2,xggsB)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'CHECK'
c     write(*,*)'K1=',K1
c     write(*,*)'K2=',K2
c     write(*,*)'K3=',K3
c     write(*,*)'N1=',N1
c     write(*,*)'N2=',N2
c     write(*,*)'N3=',N3
c     write(*,*)'xvec=',xvec
c     write(*,*)'xggs=',xggs
c     xgVEC=ZNUC*xvec*xggs


C======================================================================
      subroutine RI_G3_xgVEP(I1,J1,K1,A1,Amat1,
     *                       I2,J2,K2,A2,Amat2,
     *                       I3,J3,K3,A3,Amat3,
     *                       L1,M1,N1,B1,Bmat1,
     *                       L2,M2,N2,B2,Bmat2,
     *                       L3,M3,N3,B3,Bmat3,
     *                       gamA12,gamA13,gamA23,
     *                       gamB12,gamB13,gamB23,
     *                       xgVEP)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3

      double precision  A1
      double precision  A2
      double precision  A3
      double precision  B1
      double precision  B2
      double precision  B3

c     double precision  ZNUC
c     double precision  Cmat(3)

      double precision  gamA12
      double precision  gamA13
      double precision  gamA23
      double precision  gamB12
      double precision  gamB13
      double precision  gamB23

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Amat3(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)
      double precision  Bmat3(3)

C Variables returned
      double precision  xgVEP

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Mab
      double precision  xvep
      double precision  xggs
      double precision  abs_norm_a
      double precision  abs_norm_b
      double precision  ans
      double precision  gamma1
      double precision  gamma2

      integer  NauxBF
      integer  i
      integer  a
      integer  b
      integer  ia
      integer  LX(1000),MX(1000),NX(1000)

      double precision BX(1000)
      double precision BmatX(3,1000)
      double precision BmatX1(3)
      double precision BmatX2(3)

cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'USING RI'
      do i=1,1000
         LX(i)=0
         MX(i)=0
         NX(i)=0
         BX(i)=zero
         BmatX(1,i)=zero
         BmatX(2,i)=zero
         BmatX(3,i)=zero
      end do

c Read RI auxilliary basis input file
      open(unit=912,file='aux_p_basis.inp',status='old')
      read(912,*) NauxBF

      do i=1,NauxBF
      read(912,*)LX(i),MX(i),NX(i),
     x BX(i),BmatX(1,i),BmatX(2,i),BmatX(3,i)
      end do

C  Open ABS normalization constant file for p ABS
      open(914,file='ABSNORMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C  Open transformation matrix product file for p ABS
      open(920,file='ABSMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)


      xgVEP=zero
      ans=zero

      do b=1,NauxBF

         read(914,REC=b) abs_norm_b
         BmatX1(1)=Bmatx(1,b)
         BmatX1(2)=Bmatx(2,b)
         BmatX1(3)=Bmatx(3,b)
 
         gamma1=gamA23
         gamma2=zero

c  Calculate Overlap 
c        call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
c    x                      I3,J3,K3,A3,Amat3,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      LX(b),MX(b),NX(b),BX(b),BmatX1,
c    x                      gamma1,gamma2,xggs)
C
         call pgiovlap(I2,J2,K2,A2,Amat2,
     x                 I3,J3,K3,A3,Amat3,
     x                 L2,M2,N2,B2,Bmat2,
     x                 LX(b),MX(b),NX(b),BX(b),BmatX1,
     x                 gamma1,gamma2,xggs)


         do a=1,NauxBF

            read(914,REC=a) abs_norm_a
            BmatX2(1)=Bmatx(1,a)
            BmatX2(2)=Bmatx(2,a)
            BmatX2(3)=Bmatx(3,a)

C  Read in Mab
            call pack_2D(NauxBF,a,b,ia)
            read(920,REC=ia) Mab 

c  Calculate Vep 
            call gfvee(I1,J1,K1,A1,Amat1,
     x                 LX(a),MX(a),NX(a),BX(a),BmatX2,
     x                 L1,M1,N1,B1,Bmat1,
     x                 L3,M3,N3,B3,Bmat3,
     x                 xvep)

            ans=xvep*xggs*Mab*abs_norm_a*abs_norm_b

            xgVEP=ans+xgVEP

         end do

      end do

      xgVEP=-1.0d+00*xgVEP

      close(920)
      close(914)
      close(912)
      

      return
      end

c        call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
c    x                      L3,M3,N3,B3,Bmat3,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      LX(k3B),MX(k3B),NX(k3B),BX(k3B),BmatX1,
c    x                      gamma1,gamma2,xggsB)
C
c        call pgiovlap(I2,J2,K2,A2,Amat2,
c    x                 L3,M3,N3,B3,Bmat3,
c    x                 L2,M2,N2,B2,Bmat2,
c    x                 LX(a),MX(a),NX(a),BX(a),BmatX1,
c    x                 gamma1,gamma2,xggsB)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'CHECK'
c     write(*,*)'K1=',K1
c     write(*,*)'K2=',K2
c     write(*,*)'K3=',K3
c     write(*,*)'N1=',N1
c     write(*,*)'N2=',N2
c     write(*,*)'N3=',N3
c     write(*,*)'xvec=',xvec
c     write(*,*)'xggs=',xggs
c     xgVEC=ZNUC*xvec*xggs


C======================================================================
      subroutine RI_G3_xgVEEg2(I1,J1,K1,A1,Amat1,
     *                         I2,J2,K2,A2,Amat2,
     *                         I3,J3,K3,A3,Amat3,
     *                         L1,M1,N1,B1,Bmat1,
     *                         L2,M2,N2,B2,Bmat2,
     *                         L3,M3,N3,B3,Bmat3,
     *                         gamA12,gamA13,gamA23,
     *                         gamB12,gamB13,gamB23,
     *                         xgVEEg2)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3

      double precision  A1
      double precision  A2
      double precision  A3
      double precision  B1
      double precision  B2
      double precision  B3

c     double precision  ZNUC
c     double precision  Cmat(3)

      double precision  gamA12
      double precision  gamA13
      double precision  gamA23
      double precision  gamB12
      double precision  gamB13
      double precision  gamB23

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Amat3(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)
      double precision  Bmat3(3)

C Variables returned
      double precision  xgVEEg2

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Ma1b1
      double precision  Ma2b2
      double precision  xvee
      double precision  xggs
      double precision  abs_norm_a1
      double precision  abs_norm_a2
      double precision  abs_norm_b1
      double precision  abs_norm_b2
      double precision  ans
c     double precision  gamma1
c     double precision  gamma2

      integer  NauxBF
      integer  i
      integer  ia1
      integer  ib1
      integer  ia2
      integer  ib2
      integer  ia
      integer  LX(1000),MX(1000),NX(1000)

      double precision BX(1000)
      double precision BmatX(3,1000)
      double precision BmatXa1(3)
      double precision BmatXb1(3)
      double precision BmatXa2(3)
      double precision BmatXb2(3)

cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'USING RI'
      do i=1,1000
         LX(i)=0
         MX(i)=0
         NX(i)=0
         BX(i)=zero
         BmatX(1,i)=zero
         BmatX(2,i)=zero
         BmatX(3,i)=zero
      end do

c Read RI auxilliary basis input file
      open(unit=911,file='aux_e_basis.inp',status='old')
      read(911,*) NauxBF

      do i=1,NauxBF
      read(911,*)LX(i),MX(i),NX(i),
     x BX(i),BmatX(1,i),BmatX(2,i),BmatX(3,i)
      end do

C  Open ABS normalization constant file for e ABS
      open(913,file='ABSNORME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C  Open transformation matrix product file for e ABS
      open(919,file='ABSME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)


      xgVEEg2=zero
      ans=zero

      do ib1=1,NauxBF

         read(913,REC=ib1) abs_norm_b1
         BmatXb1(1)=Bmatx(1,ib1)
         BmatXb1(2)=Bmatx(2,ib1)
         BmatXb1(3)=Bmatx(3,ib1)
 
         do ib2=1,NauxBF

            read(913,REC=ib2) abs_norm_b2
            BmatXb2(1)=Bmatx(1,ib2)
            BmatXb2(2)=Bmatx(2,ib2)
            BmatXb2(3)=Bmatx(3,ib2)
 

c  Calculate Overlap 
            call G3ovlap(I1,J1,K1,A1,Amat1,
     x                   L2,M2,N2,B2,Bmat2,
     x                   I3,J3,K3,A3,Amat3,
     x                   LX(ib1),MX(ib1),NX(ib1),BX(ib1),BmatXb1,
     x                   LX(ib2),MX(ib2),NX(ib2),BX(ib2),BmatXb2,
     x                   L3,M3,N3,B3,Bmat3,
     x                   gamA12,gamA13,gamA23,
     x                   gamB12,gamB13,gamB23,
     x                   xggs)

            do ia1=1,NauxBF

                     read(913,REC=ia1) abs_norm_a1
               BmatXa1(1)=Bmatx(1,ia1)
               BmatXa1(2)=Bmatx(2,ia1)
               BmatXa1(3)=Bmatx(3,ia1)
 
               do ia2=1,NauxBF

                        read(913,REC=ia2) abs_norm_a2
                  BmatXa2(1)=Bmatx(1,ia2)
                  BmatXa2(2)=Bmatx(2,ia2)
                  BmatXa2(3)=Bmatx(3,ia2)

C  Read in Mab
                  call pack_2D(NauxBF,ia1,ib1,ia)
                  read(919,REC=ia) Ma1b1 
                  call pack_2D(NauxBF,ia2,ib2,ia)
                  read(919,REC=ia) Ma2b2 

c  Calculate Vep 
                  call gfvee(LX(ia1),MX(ia1),NX(ia1),BX(ia1),BmatXa1,
     x                       I2,J2,K2,A2,Amat2,
     x                       L1,M1,N1,B1,Bmat1,
     x                       LX(ia2),MX(ia2),NX(ia2),BX(ia2),BmatXa2,
     x                       xvee)

                  ans=abs_norm_a1*abs_norm_b1*Ma1b1*
     x                abs_norm_a2*abs_norm_b2*Ma2b2*
     x                xvee*xggs

                  xgVEEg2=ans+xgVEEg2

               end do
            end do
         end do
      end do

      close(919)
      close(913)
      close(911)
      

      return
      end

c        call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
c    x                      L3,M3,N3,B3,Bmat3,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      LX(k3B),MX(k3B),NX(k3B),BX(k3B),BmatX1,
c    x                      gamma1,gamma2,xggsB)
C
c        call pgiovlap(I2,J2,K2,A2,Amat2,
c    x                 L3,M3,N3,B3,Bmat3,
c    x                 L2,M2,N2,B2,Bmat2,
c    x                 LX(a),MX(a),NX(a),BX(a),BmatX1,
c    x                 gamma1,gamma2,xggsB)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'CHECK'
c     write(*,*)'K1=',K1
c     write(*,*)'K2=',K2
c     write(*,*)'K3=',K3
c     write(*,*)'N1=',N1
c     write(*,*)'N2=',N2
c     write(*,*)'N3=',N3
c     write(*,*)'xvec=',xvec
c     write(*,*)'xggs=',xggs
c     xgVEC=ZNUC*xvec*xggs


C======================================================================
      subroutine RI_G3_xgVEEg1(I1,J1,K1,A1,Amat1,
     *                         I2,J2,K2,A2,Amat2,
     *                         I3,J3,K3,A3,Amat3,
     *                         L1,M1,N1,B1,Bmat1,
     *                         L2,M2,N2,B2,Bmat2,
     *                         L3,M3,N3,B3,Bmat3,
     *                         gamA12,gamA13,gamA23,
     *                         gamB12,gamB13,gamB23,
     *                         xgVEEg1)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3

      double precision  A1
      double precision  A2
      double precision  A3
      double precision  B1
      double precision  B2
      double precision  B3

c     double precision  ZNUC
c     double precision  Cmat(3)

      double precision  gamA12
      double precision  gamA13
      double precision  gamA23
      double precision  gamB12
      double precision  gamB13
      double precision  gamB23

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Amat3(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)
      double precision  Bmat3(3)

C Variables returned
      double precision  xgVEEg1

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Mab
      double precision  xvee
      double precision  xggs
      double precision  abs_norm_a
      double precision  abs_norm_b
      double precision  ans
      double precision  gamma1
      double precision  gamma2

      integer  NauxBF
      integer  i
      integer  a
      integer  b
      integer  ia
      integer  LX(1000),MX(1000),NX(1000)

      double precision BX(1000)
      double precision BmatX(3,1000)
      double precision BmatX1(3)
      double precision BmatX2(3)

cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'USING RI'
      do i=1,1000
         LX(i)=0
         MX(i)=0
         NX(i)=0
         BX(i)=zero
         BmatX(1,i)=zero
         BmatX(2,i)=zero
         BmatX(3,i)=zero
      end do

c Read RI auxilliary basis input file
      open(unit=911,file='aux_e_basis.inp',status='old')
      read(911,*) NauxBF

      do i=1,NauxBF
      read(911,*)LX(i),MX(i),NX(i),
     x BX(i),BmatX(1,i),BmatX(2,i),BmatX(3,i)
      end do

C  Open ABS normalization constant file for e ABS
      open(913,file='ABSNORME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C  Open transformation matrix product file for e ABS
      open(919,file='ABSME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)


      xgVEEg1=zero
      ans=zero

      do b=1,NauxBF

         read(913,REC=b) abs_norm_b
         BmatX1(1)=Bmatx(1,b)
         BmatX1(2)=Bmatx(2,b)
         BmatX1(3)=Bmatx(3,b)
 
         gamma1=gamA13
         gamma2=gamB13

c  Calculate Overlap 
c        call cws_gam1_xggs(I1,J1,K1,A1,Amat1,
c    x                      I3,J3,K3,A3,Amat3,
c    x                      LX(b),MX(b),NX(b),BX(b),BmatX1,
c    x                      L3,M3,N3,B3,Bmat3,
c    x                      gamma1,gamma2,xggs)
C
         call pgiovlap(I1,J1,K1,A1,Amat1,
     x                 I3,J3,K3,A3,Amat3,
     x                 LX(b),MX(b),NX(b),BX(b),BmatX1,
     x                 L3,M3,N3,B3,Bmat3,
     x                 gamma1,gamma2,xggs)


         do a=1,NauxBF

            read(913,REC=a) abs_norm_a
            BmatX2(1)=Bmatx(1,a)
            BmatX2(2)=Bmatx(2,a)
            BmatX2(3)=Bmatx(3,a)

C  Read in Mab
            call pack_2D(NauxBF,a,b,ia)
            read(919,REC=ia) Mab 

c  Calculate Vep 
            call gfvee(LX(a),MX(a),NX(a),BX(a),BmatX2,
     x                 I2,J2,K2,A2,Amat2,
     x                 L1,M1,N1,B1,Bmat1,
     x                 L2,M2,N2,B2,Bmat2,
     x                 xvee)

            ans=xvee*xggs*Mab*abs_norm_a*abs_norm_b

            xgVEEg1=ans+xgVEEg1

         end do

      end do

      close(919)
      close(913)
      close(911)
      

      return
      end

c        call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
c    x                      L3,M3,N3,B3,Bmat3,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      LX(k3B),MX(k3B),NX(k3B),BX(k3B),BmatX1,
c    x                      gamma1,gamma2,xggsB)
C
c        call pgiovlap(I2,J2,K2,A2,Amat2,
c    x                 L3,M3,N3,B3,Bmat3,
c    x                 L2,M2,N2,B2,Bmat2,
c    x                 LX(a),MX(a),NX(a),BX(a),BmatX1,
c    x                 gamma1,gamma2,xggsB)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'CHECK'
c     write(*,*)'K1=',K1
c     write(*,*)'K2=',K2
c     write(*,*)'K3=',K3
c     write(*,*)'N1=',N1
c     write(*,*)'N2=',N2
c     write(*,*)'N3=',N3
c     write(*,*)'xvec=',xvec
c     write(*,*)'xggs=',xggs
c     xgVEC=ZNUC*xvec*xggs



C======================================================================
      subroutine RI_G3_xgVEE(I1,J1,K1,A1,Amat1,
     *                       I2,J2,K2,A2,Amat2,
     *                       I3,J3,K3,A3,Amat3,
     *                       L1,M1,N1,B1,Bmat1,
     *                       L2,M2,N2,B2,Bmat2,
     *                       L3,M3,N3,B3,Bmat3,
     *                       gamA12,gamA13,gamA23,
     *                       gamB12,gamB13,gamB23,
     *                       xgVEE)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3

      double precision  A1
      double precision  A2
      double precision  A3
      double precision  B1
      double precision  B2
      double precision  B3

c     double precision  ZNUC
c     double precision  Cmat(3)

      double precision  gamA12
      double precision  gamA13
      double precision  gamA23
      double precision  gamB12
      double precision  gamB13
      double precision  gamB23

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Amat3(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)
      double precision  Bmat3(3)

C Variables returned
      double precision  xgVEE

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Mab
      double precision  xvee
      double precision  xggs
      double precision  abs_norm_a
      double precision  abs_norm_b
      double precision  ans
      double precision  gamma1
      double precision  gamma2

      integer  NauxBF
      integer  i
      integer  a
      integer  b
      integer  ia
      integer  LX(1000),MX(1000),NX(1000)

      double precision BX(1000)
      double precision BmatX(3,1000)
      double precision BmatX1(3)
      double precision BmatX2(3)

cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'USING RI'
      do i=1,1000
         LX(i)=0
         MX(i)=0
         NX(i)=0
         BX(i)=zero
         BmatX(1,i)=zero
         BmatX(2,i)=zero
         BmatX(3,i)=zero
      end do

c Read RI auxilliary basis input file
      open(unit=911,file='aux_e_basis.inp',status='old')
      read(911,*) NauxBF

      do i=1,NauxBF
      read(911,*)LX(i),MX(i),NX(i),
     x BX(i),BmatX(1,i),BmatX(2,i),BmatX(3,i)
      end do

C  Open ABS normalization constant file for e ABS
      open(913,file='ABSNORME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C  Open transformation matrix product file for e ABS
      open(919,file='ABSME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)


      xgVEE=zero
      ans=zero

      do b=1,NauxBF

         read(913,REC=b) abs_norm_b
         BmatX1(1)=Bmatx(1,b)
         BmatX1(2)=Bmatx(2,b)
         BmatX1(3)=Bmatx(3,b)
 
         gamma1=gamA23
         gamma2=zero

c  Calculate Overlap 
c        call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
c    x                      I3,J3,K3,A3,Amat3,
c    x                      LX(b),MX(b),NX(b),BX(b),BmatX1,
c    x                      L3,M3,N3,B3,Bmat3,
c    x                      gamma1,gamma2,xggs)
C
         call pgiovlap(I2,J2,K2,A2,Amat2,
     x                 I3,J3,K3,A3,Amat3,
     x                 LX(b),MX(b),NX(b),BX(b),BmatX1,
     x                 L3,M3,N3,B3,Bmat3,
     x                 gamma1,gamma2,xggs)


         do a=1,NauxBF

            read(913,REC=a) abs_norm_a
            BmatX2(1)=Bmatx(1,a)
            BmatX2(2)=Bmatx(2,a)
            BmatX2(3)=Bmatx(3,a)

C  Read in Mab
            call pack_2D(NauxBF,a,b,ia)
            read(919,REC=ia) Mab 

c  Calculate Vep 
            call gfvee(I1,J1,K1,A1,Amat1,
     x                 LX(a),MX(a),NX(a),BX(a),BmatX2,
     x                 L1,M1,N1,B1,Bmat1,
     x                 L2,M2,N2,B2,Bmat2,
     x                 xvee)

            ans=xvee*xggs*Mab*abs_norm_a*abs_norm_b

            xgVEE=ans+xgVEE

         end do

      end do

      close(919)
      close(913)
      close(911)
      

      return
      end

c        call cws_gam1_xggs(I2,J2,K2,A2,Amat2,
c    x                      L3,M3,N3,B3,Bmat3,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      LX(k3B),MX(k3B),NX(k3B),BX(k3B),BmatX1,
c    x                      gamma1,gamma2,xggsB)
C
c        call pgiovlap(I2,J2,K2,A2,Amat2,
c    x                 L3,M3,N3,B3,Bmat3,
c    x                 L2,M2,N2,B2,Bmat2,
c    x                 LX(a),MX(a),NX(a),BX(a),BmatX1,
c    x                 gamma1,gamma2,xggsB)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c     write(*,*)'CHECK'
c     write(*,*)'K1=',K1
c     write(*,*)'K2=',K2
c     write(*,*)'K3=',K3
c     write(*,*)'N1=',N1
c     write(*,*)'N2=',N2
c     write(*,*)'N3=',N3
c     write(*,*)'xvec=',xvec
c     write(*,*)'xggs=',xggs
c     xgVEC=ZNUC*xvec*xggs

