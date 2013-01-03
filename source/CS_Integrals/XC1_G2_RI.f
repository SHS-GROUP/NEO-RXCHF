C======================================================================
      subroutine RI_G2_xgVEC(I1,J1,K1,A1,Amat1,
     *                       I2,J2,K2,A2,Amat2,
     *                       L1,M1,N1,B1,Bmat1,
     *                       L2,M2,N2,B2,Bmat2,
     *                       gamma1,gamma2,
     *                       Cmat,ZNUC,
     *                       xgVEC)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  L1,M1,N1
      integer  L2,M2,N2

      double precision  A1
      double precision  A2
      double precision  B1
      double precision  B2

      double precision  ZNUC
      double precision  Cmat(3)

      double precision  gamma1
      double precision  gamma2

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)

C Variables returned
      double precision  xgVEC

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Mab
      double precision  xvec
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

C  Open ABS normalization constant file for p ABS
      open(913,file='ABSNORME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C  Open transformation matrix product file for p ABS
      open(919,file='ABSME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)


      xgVEC=zero
      ans=zero

      do b=1,NauxBF

         read(913,REC=b) abs_norm_b
         BmatX1(1)=Bmatx(1,b)
         BmatX1(2)=Bmatx(2,b)
         BmatX1(3)=Bmatx(3,b)
 
c  Calculate Overlap 
c        call cws_gam1_xggs(I1,J1,K1,A1,Amat1,
c    x                      I2,J2,K2,A2,Amat2,
c    x                      LX(b),MX(b),NX(b),BX(b),BmatX1,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      gamma1,gamma2,xggs)
C
         call pgiovlap(I1,J1,K1,A1,Amat1,
     x                 I2,J2,K2,A2,Amat2,
     x                 LX(b),MX(b),NX(b),BX(b),BmatX1,
     x                 L2,M2,N2,B2,Bmat2,
     x                 gamma1,gamma2,xggs)


         do a=1,NauxBF

            read(913,REC=a) abs_norm_a
            BmatX2(1)=Bmatx(1,a)
            BmatX2(2)=Bmatx(2,a)
            BmatX2(3)=Bmatx(3,a)

C  Read in Mab
            call pack_2D(NauxBF,a,b,ia)
            read(919,REC=ia) Mab 

c  Calculate VeC 
            call gfvec(LX(a),MX(a),NX(a),BX(a),BmatX2,
     x                 L1,M1,N1,B1,Bmat1,
     x                 Cmat,xvec)

            ans=xvec*xggs*Mab*abs_norm_a*abs_norm_b

            xgVEC=ans+xgVEC

         end do

      end do

         xgVEC=-1.0d+00*ZNUC*xgVEC

      close(919)
      close(913)
      close(911)
      

      return
      end

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
      subroutine RI_G2_xgVPC(I1,J1,K1,A1,Amat1,
     *                       I2,J2,K2,A2,Amat2,
     *                       L1,M1,N1,B1,Bmat1,
     *                       L2,M2,N2,B2,Bmat2,
     *                       gamma1,gamma2,
     *                       Cmat,ZNUC,
     *                       xgVPC)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  L1,M1,N1
      integer  L2,M2,N2

      double precision  A1
      double precision  A2
      double precision  B1
      double precision  B2

      double precision  ZNUC
      double precision  Cmat(3)

      double precision  gamma1
      double precision  gamma2

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)

C Variables returned
      double precision  xgVPC

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Mab
      double precision  xvpc
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


      xgVPC=zero
      ans=zero

      do b=1,NauxBF

         read(914,REC=b) abs_norm_b
         BmatX1(1)=Bmatx(1,b)
         BmatX1(2)=Bmatx(2,b)
         BmatX1(3)=Bmatx(3,b)
 
c  Calculate Overlap 
c        call cws_gam1_xggs(I1,J1,K1,A1,Amat1,
c    x                      LX(b),MX(b),NX(b),BX(b),BmatX1,
c    x                      L1,M1,N1,B1,Bmat1,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      gamma1,gamma2,xggs)
C
         call pgiovlap(I1,J1,K1,A1,Amat1,
     x                 LX(b),MX(b),NX(b),BX(b),BmatX1,
     x                 L1,M1,N1,B1,Bmat1,
     x                 L2,M2,N2,B2,Bmat2,
     x                 gamma1,gamma2,xggs)


         do a=1,NauxBF

            read(914,REC=a) abs_norm_a
            BmatX2(1)=Bmatx(1,a)
            BmatX2(2)=Bmatx(2,a)
            BmatX2(3)=Bmatx(3,a)

C  Read in Mab
            call pack_2D(NauxBF,a,b,ia)
            read(920,REC=ia) Mab 

c  Calculate VeC 
            call gfvec(I2,J2,K2,A2,Amat2,
     x                 LX(a),MX(a),NX(a),BX(a),BmatX2,
     x                 Cmat,xvpc)

            ans=xvpc*xggs*Mab*abs_norm_a*abs_norm_b

            xgVPC=ans+xgVPC

         end do

      end do

         xgVPC=ZNUC*xgVPC

      close(920)
      close(914)
      close(912)
      

      return
      end

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
      subroutine RI_G2_xgVEP(I1,J1,K1,A1,Amat1,
     *                       I2,J2,K2,A2,Amat2,
     *                       L1,M1,N1,B1,Bmat1,
     *                       L2,M2,N2,B2,Bmat2,
     *                       gamma1,gamma2,
     *                       xgVEP)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  L1,M1,N1
      integer  L2,M2,N2

      double precision  A1
      double precision  A2
      double precision  B1
      double precision  B2

      double precision  gamma1
      double precision  gamma2

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)

C Variables returned
      double precision  xgVEP

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Ma1b1
      double precision  Ma2b2
      double precision  xvep
      double precision  xggs
      double precision  abs_norm_a1
      double precision  abs_norm_a2
      double precision  abs_norm_b1
      double precision  abs_norm_b2
      double precision  ans

      integer  i,ia
      integer  ia1,ib1
      integer  ia2,ib2

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

      double precision BmatXa2(3)
      double precision BmatXb2(3)

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


      xgVEP=zero

      do ib1=1,NauxBFE

                read(913,REC=ib1) abs_norm_b1
         BmatXb1(1)=BmatXE(1,ib1)
         BmatXb1(2)=BmatXE(2,ib1)
         BmatXb1(3)=BmatXE(3,ib1)
 
         do ib2=1,NauxBFP

                   read(914,REC=ib2) abs_norm_b2
            BmatXb2(1)=BmatXP(1,ib2)
            BmatXb2(2)=BmatXP(2,ib2)
            BmatXb2(3)=BmatXP(3,ib2)
 

c  Calculate Overlap 
c        call cws_gam1_xggs(I1,J1,K1,A1,Amat1,
c    x                      LXP(ib2),MXP(ib2),NXP(ib2),BXP(ib2),BmatXb2,
c    x                      LXE(ib1),MXE(ib1),NXE(ib1),BXE(ib1),BmatXb1,
c    x                      L2,M2,N2,B2,Bmat2,
c    x                      gamma1,gamma2,xggs)
C
         call pgiovlap(I1,J1,K1,A1,Amat1,
     x                 LXP(ib2),MXP(ib2),NXP(ib2),BXP(ib2),BmatXb2,
     x                 LXE(ib1),MXE(ib1),NXE(ib1),BXE(ib1),BmatXb1,
     x                 L2,M2,N2,B2,Bmat2,
     x                 gamma1,gamma2,xggs)

            do ia1=1,NauxBFE

                      read(913,REC=ia1) abs_norm_a1
               BmatXa1(1)=BmatXE(1,ia1)
               BmatXa1(2)=BmatXE(2,ia1)
               BmatXa1(3)=BmatXE(3,ia1)
 
               ans=zero

               do ia2=1,NauxBFP

                         read(914,REC=ia2) abs_norm_a2
                  BmatXa2(1)=BmatXP(1,ia2)
                  BmatXa2(2)=BmatXP(2,ia2)
                  BmatXa2(3)=BmatXP(3,ia2)

C  Read in Mab
c                 call pack_2D(NauxBFE,ia1,ib1,ia)
                  call pack_2D(NauxBFE,ib1,ia1,ia)
                  read(919,REC=ia) Ma1b1 
c                 call pack_2D(NauxBFP,ia2,ib2,ia)
                  call pack_2D(NauxBFP,ib2,ia2,ia)
                  read(920,REC=ia) Ma2b2 

c  Calculate Vep 
c                 call gfvee(I1,J1,K1,A1,Amat1,
c    x                       I2,J2,K2,A2,Amat2,
c    x                      LXE(ia1),MXE(ia1),NXE(ia1),BXE(ia1),BmatXa1,
c    x                      LXP(ia2),MXP(ia2),NXP(ia2),BXP(ia2),BmatXa2,
c    x                       xvep)
                 call gfvee(LXE(ia1),MXE(ia1),NXE(ia1),BXE(ia1),BmatXa1,
     x                      I2,J2,K2,A2,Amat2,
     x                      L1,M1,N1,B1,Bmat1,
     x                      LXP(ia2),MXP(ia2),NXP(ia2),BXP(ia2),BmatXa2,
     x                      xvep)

                  ans=abs_norm_a1*abs_norm_b1*Ma1b1*
     x                abs_norm_a2*abs_norm_b2*Ma2b2*
     x                xvep*xggs

                  xgVEP=ans+xgVEP

               end do
            end do
         end do
      end do

      xgVEP=-1.0d+00*xgVEP

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
      subroutine RI_G2_xgVEPg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        gamma1,gamma2,
     *                        xgVEPg)

C======================================================================
      implicit none

C Input variables
      integer  I1,J1,K1
      integer  I2,J2,K2
      integer  L1,M1,N1
      integer  L2,M2,N2

      double precision  A1
      double precision  A2
      double precision  B1
      double precision  B2

      double precision  gamma1
      double precision  gamma2

      double precision  Amat1(3)
      double precision  Amat2(3)
      double precision  Bmat1(3)
      double precision  Bmat2(3)

C Variables returned
      double precision  xgVEPg

C Local variables
      double precision  zero
      parameter(zero=0.0d+00)
      double precision  Ma1b1
      double precision  Ma2b2
      double precision  xvep
      double precision  xggs
      double precision  abs_norm_a1
      double precision  abs_norm_a2
      double precision  abs_norm_b1
      double precision  abs_norm_b2
      double precision  ans

      integer  i,ia
      integer  ia1,ib1
      integer  ia2,ib2

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

      double precision BmatXa2(3)
      double precision BmatXb2(3)

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


      xgVEPg=zero
      ans=zero

      do ib1=1,NauxBFE

                read(913,REC=ib1) abs_norm_b1
         BmatXb1(1)=BmatXE(1,ib1)
         BmatXb1(2)=BmatXE(2,ib1)
         BmatXb1(3)=BmatXE(3,ib1)
 
         do ib2=1,NauxBFP

                   read(914,REC=ib2) abs_norm_b2
            BmatXb2(1)=BmatXP(1,ib2)
            BmatXb2(2)=BmatXP(2,ib2)
            BmatXb2(3)=BmatXP(3,ib2)
 

c  Calculate Overlap 
c        call cws_gam1_xggs(LXE(ib1),MXE(ib1),NXE(ib1),BXE(ib1),BmatXb1,
c    x                      I2,J2,K2,A2,Amat2,
c    x                      L1,M1,N2,B1,Bmat1,
c    x                      LXP(ib2),MXP(ib2),NXP(ib2),BXP(ib2),BmatXb2,
c    x                      gamma1,gamma2,xggs)
C
         call pgiovlap(LXE(ib1),MXE(ib1),NXE(ib1),BXE(ib1),BmatXb1,
     x                 I2,J2,K2,A2,Amat2,
     x                 L1,M1,N2,B1,Bmat1,
     x                 LXP(ib2),MXP(ib2),NXP(ib2),BXP(ib2),BmatXb2,
     x                 gamma1,gamma2,xggs)

            do ia1=1,NauxBFE

                      read(913,REC=ia1) abs_norm_a1
               BmatXa1(1)=BmatXE(1,ia1)
               BmatXa1(2)=BmatXE(2,ia1)
               BmatXa1(3)=BmatXE(3,ia1)
 
               do ia2=1,NauxBFP

                         read(914,REC=ia2) abs_norm_a2
                  BmatXa2(1)=BmatXP(1,ia2)
                  BmatXa2(2)=BmatXP(2,ia2)
                  BmatXa2(3)=BmatXP(3,ia2)

C  Read in Mab
                  call pack_2D(NauxBFE,ia1,ib1,ia)
                  read(919,REC=ia) Ma1b1 
                  call pack_2D(NauxBFP,ia2,ib2,ia)
                  read(920,REC=ia) Ma2b2 

c  Calculate Vep 
                  call gfvee(I1,J1,K1,A1,Amat1,
     x                      LXP(ia2),MXP(ia2),NXP(ia2),BXP(ia2),BmatXa2,
     x                      LXE(ia1),MXE(ia1),NXE(ia1),BXE(ia1),BmatXa1,
     x                       L2,M2,N2,B2,Bmat2,
     x                       xvep)

                  ans=abs_norm_a1*abs_norm_b1*Ma1b1*
     x                abs_norm_a2*abs_norm_b2*Ma2b2*
     x                xvep*xggs

                  xgVEPg=ans+xgVEPg

               end do
            end do
         end do
      end do

      xgVEPg=-1.0d+00*xgVEPg

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

