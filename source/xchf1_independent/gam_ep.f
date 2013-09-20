C======================================================================
      subroutine calc_GAM_ep(nebf,npebf,npbf,ng1,
     x                       AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
c    x                       GAM_ep)

C======================================================================
      implicit none
C Input Variables
      integer nebf,npebf,npbf,ng1
C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
C-------Basis Set Info-------)
C Local variables
C--------------------------------(
C Basis set-related local variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2

      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
C--------------------------------)
      integer ia
      integer ip,jp,ie1,je1
      integer iec1,jec1
      double precision GAM_ep(ng1)
      double precision Cof_ip,Cof_jp
      double precision Cof_ie1,Cof_je1
      double precision ans
      double precision zero
      parameter(zero=0.0d+00)

C  zero out the 1-D arrays to hold the contracted integrals
      do ia=1,ng1
         GAM_ep(ia)=zero
      end do

      do ip=1,npbf
         do jp=1,npbf
            do ie1=1,npebf
               do je1=1,npebf

C  Map from primitive BF indices to contracted indices
c                 call MPEB2C(ie1,iec1)
c                 call MPEB2C(je1,jec1)
                  iec1=AMPEB2C(ie1)
                  jec1=AMPEB2C(je1)
C  Get primitive Electron Basis Function Contraction Coefficients
c                 call GEBFCC(ie1,Cof_ie1)
c                 call GEBFCC(je1,Cof_je1)
                  Cof_ie1=AGEBFCC(ie1)
                  Cof_je1=AGEBFCC(je1)
C  Get Nuclear Basis Function Contraction Coefficients
c                 call GNBFCC(ip,Cof_ip)
c                 call GNBFCC(jp,Cof_jp)
                  Cof_ip=AGNBFCC(ip)
                  Cof_jp=AGNBFCC(jp)
C  Map the 4-index contracted integral to 1-D:
                  call pack_4D(nebf,nebf,npbf,
     x                         jec1,iec1,jp,ip,ia)

C Get Basis set info:
                  A1=ELCEX(ie1)
                  I1=ELCAM(ie1,1)
                  J1=ELCAM(ie1,2)
                  K1=ELCAM(ie1,3)
                  Amat1(1)=ELCBFC(ie1,1)
                  Amat1(2)=ELCBFC(ie1,2)
                  Amat1(3)=ELCBFC(ie1,3)

                  A2=NUCEX(ip)
                  I2=NUCAM(ip,1)
                  J2=NUCAM(ip,2)
                  K2=NUCAM(ip,3)
                  Amat2(1)=NUCBFC(ip,1)
                  Amat2(2)=NUCBFC(ip,2)
                  Amat2(3)=NUCBFC(ip,3)

                  B1=ELCEX(je1)
                  L1=ELCAM(je1,1)
                  M1=ELCAM(je1,2)
                  N1=ELCAM(je1,3)
                  Bmat1(1)=ELCBFC(je1,1)
                  Bmat1(2)=ELCBFC(je1,2)
                  Bmat1(3)=ELCBFC(je1,3)

                  B2=NUCEX(jp)
                  L2=NUCAM(jp,1)
                  M2=NUCAM(jp,2)
                  N2=NUCAM(jp,3)
                  Bmat2(1)=NUCBFC(jp,1)
                  Bmat2(2)=NUCBFC(jp,2)
                  Bmat2(3)=NUCBFC(jp,3)

C  Calculate primitive integrals
                  call xcalc_GAM_ep(I1,J1,K1,A1,Amat1,
     x                              I2,J2,K2,A2,Amat2,
     x                              L1,M1,N1,B1,Bmat1,
     x                              L2,M2,N2,B2,Bmat2,ans)


c     write(*,*)'ia=',ia
c     write(*,*)'Cof_ie1=',Cof_ie1
c     write(*,*)'Cof_je1=',Cof_je1
c     write(*,*)'Cof_ip=',Cof_ip
c     write(*,*)'Cof_jp=',Cof_jp
c     write(*,*)'ans=',ans
                  GAM_ep(ia)=GAM_ep(ia)+ans*
     x                      Cof_ip*Cof_jp*Cof_ie1*Cof_je1
c     write(*,*)'GAM_ep(ia)=',GAM_ep(ia)
c     write(*,*)
c                 GAM_ep(ip,jp,ie1,je1)=ans*
c    x                                  pfct(ip)*pfct(jp)*
c    x                                  efct(ie1)*efct(je1)

               end do
            end do
         end do
      end do

C Write electron-proton integrals to file
      open(811,file='GAM_ep.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ip=1,npbf
         do jp=1,npbf
            do iec1=1,nebf
               do jec1=1,nebf

C  Map the 4-index contracted integral to 1-D:
                  call pack_4D(nebf,nebf,npbf,
     x                         jec1,iec1,jp,ip,ia)

                  write(811,REC=ia) GAM_ep(ia)

               end do
            end do
         end do
      end do

      close(811)


      return
      end

C======================================================================
      subroutine xcalc_GAM_ep(I1,J1,K1,A1,Amat1,
     x                        I2,J2,K2,A2,Amat2,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L2,M2,N2,B2,Bmat2,ans)

C======================================================================
      implicit none
C Input Variables
      integer I1,J1,K1 
      integer I2,J2,K2         
      integer L1,M1,N1
      integer L2,M2,N2
      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
C Variables Returned
      double precision ans
C Local Variables
      double precision coul_sign
      LOGICAL DEBUG


c     debug=.true.
      debug=.false.

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      call gfvee(I1,J1,K1,A1,Amat1,
     1           I2,J2,K2,A2,Amat2,
     2           L1,M1,N1,B1,Bmat1,
     3           L2,M2,N2,B2,Bmat2,
     4           ans)

      coul_sign=-1.d+00
      ans=ans*coul_sign

C>>>>>>>>>> debug
      if(debug) then
         write(*,*)'>>> GAM_ep <<<'
c        write(*,*)'ijkl=',ie1,je1,ip,jp
         write(*,*)'I1 J1 K1 Alpha1 Amat1',I1,J1,K1,A1,Amat1
         write(*,*)'I2 J2 K2 Alpha2 Amat2',I2,J2,K2,A2,Amat2
         write(*,*)'L1 M1 N1 Beta1 Bmat1 ',L1,M1,N1,B1,Bmat1
         write(*,*)'L2 M2 N2 Beta2 Bmat2 ',L2,M2,N2,B2,Bmat2
         write(*,*)'ans=',ans
         write(*,*)
      end if

      return
      end


