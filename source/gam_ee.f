C======================================================================
      subroutine calc_GAM_ee(nebf,npebf,ngee,
     x                       AMPEB2C,AGEBFCC,ELCEX,ELCAM,ELCBFC)

C======================================================================
      implicit none

      integer nebf
      integer npebf
      integer ngee
c     double precision efct(nebf)
      double precision GAM_ee(ngee)

C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
c     integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
c     double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
c     double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
c     double precision AGNBFCC(npbf)  ! Nuclear contract coef
C-------Basis Set Info-------)

C Local variables
      integer ie1,je1,ie2,je2
      integer iec1,jec1,iec2,jec2
      integer ia
      double precision Cof_ie1,Cof_je1
      double precision Cof_ie2,Cof_je2
      double precision ans
      double precision zero
      parameter(zero=0.0d+00)

      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2
      double precision A1,Amat1(3) 
      double precision A2,Amat2(3) 
      double precision B1,Bmat1(3) 
      double precision B2,Bmat2(3) 


      write(*,*)
      write(*,*)'**************************************'
      write(*,*)'    Computing GAM_ee Integrals    '
      write(*,*)
      write(*,*)'nebf =',nebf
      write(*,*)'npebf=',npebf
      write(*,*)'ngee =',ngee
      write(*,*)'**************************************'
      write(*,*)

C  zero out the 1-D array to hold the contracted integrals
      do ia=1,ngee
         GAM_ee(ia)=zero
      end do

      do ie1=1,npebf
         do je1=1,npebf
            do ie2=1,npebf
               do je2=1,npebf

C  Map from primitive BF indices to contracted indices
c                 call MPEB2C(ie1,iec1)
c                 call MPEB2C(ie2,iec2)
c                 call MPEB2C(je1,jec1)
c                 call MPEB2C(je2,jec2)
                  iec1=AMPEB2C(ie1)
                  iec2=AMPEB2C(ie2)
                  jec1=AMPEB2C(je1)
                  jec2=AMPEB2C(je2)
C  Get primitive Electron Basis Function Contraction Coefficients
c                 call GEBFCC(ie1,Cof_ie1)
c                 call GEBFCC(ie2,Cof_ie2)
c                 call GEBFCC(je1,Cof_je1)
c                 call GEBFCC(je2,Cof_je2)
                  Cof_ie1=AGEBFCC(ie1)
                  Cof_ie2=AGEBFCC(ie2)
                  Cof_je1=AGEBFCC(je1)
                  Cof_je2=AGEBFCC(je2)
C  Map the 4-index contracted integral to 1-D:
                  call pack_4D(nebf,nebf,nebf,
     x                         jec2,iec2,jec1,iec1,ia)

C Get Basis set info:
                  A1=ELCEX(ie1)
                  I1=ELCAM(ie1,1)
                  J1=ELCAM(ie1,2)
                  K1=ELCAM(ie1,3)
                  Amat1(1)=ELCBFC(ie1,1)
                  Amat1(2)=ELCBFC(ie1,2)
                  Amat1(3)=ELCBFC(ie1,3)

                  A2=ELCEX(ie2)
                  I2=ELCAM(ie2,1)
                  J2=ELCAM(ie2,2)
                  K2=ELCAM(ie2,3)
                  Amat2(1)=ELCBFC(ie2,1)
                  Amat2(2)=ELCBFC(ie2,2)
                  Amat2(3)=ELCBFC(ie2,3)

                  B1=ELCEX(je1)
                  L1=ELCAM(je1,1)
                  M1=ELCAM(je1,2)
                  N1=ELCAM(je1,3)
                  Bmat1(1)=ELCBFC(je1,1)
                  Bmat1(2)=ELCBFC(je1,2)
                  Bmat1(3)=ELCBFC(je1,3)

                  B2=ELCEX(je2)
                  L2=ELCAM(je2,1)
                  M2=ELCAM(je2,2)
                  N2=ELCAM(je2,3)
                  Bmat2(1)=ELCBFC(je2,1)
                  Bmat2(2)=ELCBFC(je2,2)
                  Bmat2(3)=ELCBFC(je2,3)


C Calculate the primitive integral
C  calc <ie1 ie2|V^ee|je1 je2> 
c                 call xcalc_GAM_ee(ie1,je1,ie2,je2,ans)
                  call xcalc_GAM_ee(I1,J1,K1,A1,Amat1,
     x                              I2,J2,K2,A2,Amat2,
     x                              L1,M1,N1,B1,Bmat1,
     x                              L2,M2,N2,B2,Bmat2,
     x                              ans)

                  GAM_ee(ia)=GAM_ee(ia)+ans*
     x                      Cof_ie1*Cof_je1*
     x                      Cof_ie2*Cof_je2

c                 write(*,*)'ia=',ia,' GAM_ee=',GAM_ee(ia)

               end do
            end do
         end do
      end do

C Write 2 electron integrals to file
      open(810,file='GAM_ee.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do iec1=1,nebf
         do jec1=1,nebf
            do iec2=1,nebf
               do jec2=1,nebf

C  Map the 4-index contracted integral to 1-D:
                  call pack_4D(nebf,nebf,nebf,
     x                         jec2,iec2,jec1,iec1,ia)

                  write(810,REC=ia) GAM_ee(ia)

               end do
            end do
         end do
      end do

      close(810)
c     do ia=1,ngee
c        write(*,*)'ia=',ia,' GAM_ee=',GAM_ee(ia)
c     end do

      return
      end

C======================================================================
c     subroutine xcalc_GAM_ee(ie1,je1,ie2,je2,ans)
      subroutine xcalc_GAM_ee(I1,J1,K1,A1,Amat1,
     x                        I2,J2,K2,A2,Amat2,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L2,M2,N2,B2,Bmat2,
     x                        ans)

C======================================================================
      implicit none

      LOGICAL DEBUG

C Input Variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2
      double precision A1,Amat1(3) 
      double precision A2,Amat2(3) 
      double precision B1,Bmat1(3) 
      double precision B2,Bmat2(3) 

c     integer ie1,je1
c     integer ie2,je2

C Variables Returned
      double precision ans

C Local Variables

C     debug=.true.
      debug=.false.

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c     call get_BF(1,ie1,I1,J1,K1,A1,Amat1)
c     call get_BF(1,je1,L1,M1,N1,B1,Bmat1)
c     call get_BF(1,ie2,I2,J2,K2,A2,Amat2)
c     call get_BF(1,je2,L2,M2,N2,B2,Bmat2)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


      call gfvee(I1,J1,K1,A1,Amat1,
     1           I2,J2,K2,A2,Amat2,
     2           L1,M1,N1,B1,Bmat1,
     3           L2,M2,N2,B2,Bmat2,
     4           ans)

C>>>>>>>>>> debug
      if(debug) then
         write(*,*)'>>> GAM_ee <<<'
c        write(*,*)'ijkl=',ie1,je1,ie2,je2
         write(*,*)'I1 J1 K1 Alpha1 Amat1',I1,J1,K1,A1,Amat1
         write(*,*)'I2 J2 K2 Alpha2 Amat2',I2,J2,K2,A2,Amat2
         write(*,*)'L1 M1 N1 Beta1 Bmat1 ',L1,M1,N1,B1,Bmat1
         write(*,*)'L2 M2 N2 Beta2 Bmat2 ',L2,M2,N2,B2,Bmat2
         write(*,*)'ans=',ans
         write(*,*)
      end if


      return
      end

