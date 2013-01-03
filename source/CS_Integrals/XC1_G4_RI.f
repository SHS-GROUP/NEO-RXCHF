C RI_G4_xgVEEg2.f
C RI_G4_xgVEPg.f

C======================================================================
      subroutine RI_G4_xgVEPg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        I4,J4,K4,A4,Amat4,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L4,M4,N4,B4,Bmat4,
     x                        gamA14,gamA24,gamA34,
     x                        gamB14,gamB24,gamB34,
     *                        xgVEPg)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  I4,J4,K4
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3
      integer  L4,M4,N4

      double precision  A1
      double precision  A2
      double precision  A3
      double precision  A4
      double precision  B1
      double precision  B2
      double precision  B3
      double precision  B4

      double precision  gamA14
      double precision  gamA24
      double precision  gamA34
      double precision  gamB14
      double precision  gamB24
      double precision  gamB34

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Amat3(3)
      double precision  Amat4(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)
      double precision  Bmat3(3)
      double precision  Bmat4(3)

C Variables returned
      double precision  xgVEPg

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Mab
      double precision  xvee
      double precision  xggs
      double precision  abs_norm_a
      double precision  abs_norm_b
      double precision  ans

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

C  Open ABS normalization constant file for e ABS
      open(914,file='ABSNORMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C  Open transformation matrix product file for e ABS
      open(920,file='ABSMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)


      xgVEPg=zero
      ans=zero

      do b=1,NauxBF

         read(914,REC=b) abs_norm_b
         BmatX1(1)=Bmatx(1,b)
         BmatX1(2)=Bmatx(2,b)
         BmatX1(3)=Bmatx(3,b)
 
c  Calculate Overlap 
         call G3ovlap(I1,J1,K1,A1,Amat1,
     *                I2,J2,K2,A2,Amat2,
     *                LX(b),MX(b),NX(b),BX(b),BmatX1,
     *                L1,M1,N1,B1,Bmat1,
     *                L2,M2,N2,B2,Bmat2,
     *                L4,M4,N4,B4,Bmat4,
     *                zero,gamA14,zero,
     *                zero,zero,gamB24,
     *                xggs)

c    *                gamA12,gamA13,gamA23,
c    *                gamB12,gamB13,gamB23,


         do a=1,NauxBF

            read(914,REC=a) abs_norm_a
            BmatX2(1)=Bmatx(1,a)
            BmatX2(2)=Bmatx(2,a)
            BmatX2(3)=Bmatx(3,a)

C  Read in Mab
            call pack_2D(NauxBF,a,b,ia)
            read(920,REC=ia) Mab 

c  Calculate Vep 
            call gfvee(I3,J3,K3,A3,Amat3,
     x                 I4,J4,K4,A4,Amat4,
     x                 L3,M3,N3,B3,Bmat3,
     x                 LX(a),MX(a),NX(a),BX(a),BmatX2,
     x                 xvee)

            ans=xvee*xggs*Mab*abs_norm_a*abs_norm_b

            xgVEPg=ans+xgVEPg

         end do

      end do

      close(920)
      close(914)
      close(912)
      
      xgVEPg=-1.0d+00*xgVEPg


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
      subroutine RI_G4_xgVEEg2(I1,J1,K1,A1,Amat1,
     *                         I2,J2,K2,A2,Amat2,
     *                         I3,J3,K3,A3,Amat3,
     *                         I4,J4,K4,A4,Amat4,
     *                         L1,M1,N1,B1,Bmat1,
     *                         L2,M2,N2,B2,Bmat2,
     *                         L3,M3,N3,B3,Bmat3,
     *                         L4,M4,N4,B4,Bmat4,
     x                         gamA14,gamA24,gamA34,
     x                         gamB14,gamB24,gamB34,
     *                         xgVEEg2)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  I3,J3,K3
      integer  I4,J4,K4
      integer  L1,M1,N1
      integer  L2,M2,N2
      integer  L3,M3,N3
      integer  L4,M4,N4

      double precision  A1
      double precision  A2
      double precision  A3
      double precision  A4
      double precision  B1
      double precision  B2
      double precision  B3
      double precision  B4

      double precision  gamA14
      double precision  gamA24
      double precision  gamA34
      double precision  gamB14
      double precision  gamB24
      double precision  gamB34

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Amat3(3)
      double precision  Amat4(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)
      double precision  Bmat3(3)
      double precision  Bmat4(3)

C Variables returned
      double precision  xgVEEg2

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Mab
      double precision  xvee
      double precision  xggs
      double precision  abs_norm_a
      double precision  abs_norm_b
      double precision  ans

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


      xgVEEg2=zero
      ans=zero

      do b=1,NauxBF

         read(913,REC=b) abs_norm_b
         BmatX1(1)=Bmatx(1,b)
         BmatX1(2)=Bmatx(2,b)
         BmatX1(3)=Bmatx(3,b)
 
c  Calculate Overlap 
         call G3ovlap(I1,J1,K1,A1,Amat1,
     *                I3,J3,K3,A3,Amat3,
     *                I4,J4,K4,A4,Amat4,
     *                LX(b),MX(b),NX(b),BX(b),BmatX1,
     *                L3,M3,N3,B3,Bmat3,
     *                L4,M4,N4,B4,Bmat4,
     *                zero,gamA14,zero,
     *                zero,zero,gamB34,
     *                xggs)

c    x                   gamA12,gamA13,gamA23,
c    x                   gamB12,gamB13,gamB23,


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

            xgVEEg2=ans+xgVEEg2

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

