C----------DIRECT-ACCESS-FILES-FOR-XCHF--------------------------------(
c     ----------XC-NEO------------
c          801  ::   GAM_1.ufm  
c          802  ::   GAM_1s.ufm 
c          803  ::   GAM_2.ufm  
c          804  ::   GAM_2s.ufm 
c          805  ::   GAM_3.ufm  
c          806  ::   GAM_4.ufm  
c     -------STANDARD-NEO---------
c          810  ::   GAM_ee.ufm 
c          811  ::   GAM_ep.ufm 
c          812  ::   GAM_ecr.ufm
c          813  ::   GAM_pcr.ufm
c          814  ::   eovlap.ufm 
c          815  ::   novlap.ufm 
C----------DIRECT-ACCESS-FILES-FOR-XCHF--------------------------------)
C
C          800  ::   ENUCRP.dat  classical nuclear repulsion energy
C
C---------FORMATTED-FILES-FOR-MO-COEFFICIENTS--------------------------(
c          850  ::   guessCE.inp               
c          851  ::   guessCP.inp               
c          852  ::   FinalCE.dat
c          853  ::   FinalCP.dat
C---------FORMATTED-FILES-FOR-MO-COEFFICIENTS--------------------------)
c
C=======================================================================
      subroutine read_GAM_ee(ne,ngee,GAM_ee)
C=======================================================================
      implicit none
C Input Variables
      integer ne,ngee
C Variables Returned
      double precision GAM_ee(ngee)
C Local Variables
      integer ia
      integer ip,jp
      integer ie1,je1
      integer ie2,je2
      double precision ans

      open(810,file='GAM_ee.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ie1=1,ne
       do je1=1,ne
        do ie2=1,ne
         do je2=1,ne

            call pack_4D(ne,ne,ne,
     x                   je2,ie2,je1,ie1,ia)

            read(810,REC=ia) ans

            GAM_ee(ia)=ans

         end do
        end do
       end do
      end do

      close(810)


      return
      end

C=======================================================================
      subroutine read_GAM_ep(ne,np,ng1,GAM_ep)
C=======================================================================
      implicit none
C Input Variables
      integer ne,np,ng1
C Variables Returned
      double precision GAM_ep(ng1)
C Local Variables
      integer ia
      integer ip,jp
      integer ie1,je1
      double precision ans

      open(811,file='GAM_ep.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ip=1,np
       do jp=1,np
        do ie1=1,ne
         do je1=1,ne

C  Map the 4-index contracted integral to 1-D:
            call pack_4D(ne,ne,np,
     x                   je1,ie1,jp,ip,ia)

            read(811,REC=ia) ans 

            GAM_ep(ia)=ans

         end do
        end do
       end do
      end do

      close(811)


      return
      end

C=======================================================================
      subroutine read_GAM_ecore(ne,ne2,GAM_ecore)
C=======================================================================
      implicit none
C Input Variables
      integer ne,ne2
C Variables Returned
      double precision GAM_ecore(ne2)
C Local Variables
      integer ia
      integer ie,je
      double precision ans

      open(812,file='GAM_ecr.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ie=1,ne
         do je=1,ne

C  Map the 2-index contracted integral to 1-D:
            call pack_2D(ne,je,ie,ia)

            read(812,REC=ia) GAM_ecore(ia)

         end do
      end do

      close(812)



      return
      end

C=======================================================================
      subroutine read_GAM_pcore(np,np2,GAM_pcore)
C=======================================================================
      implicit none
C Input Variables
      integer np,np2
C Variables Returned
      double precision GAM_pcore(np2)
C Local Variables
      integer ia
      integer ip,jp
      double precision ans

      open(813,file='GAM_pcr.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ip=1,np
         do jp=1,np

C  Map the 2-index contracted integral to 1-D:
            call pack_2D(np,jp,ip,ia)

            read(813,REC=ia) GAM_pcore(ia)

         end do
      end do

      close(813)



      return
      end

C=======================================================================
      subroutine read_nuc_ovlap(np,GS)
C=======================================================================
      implicit none
C Input Variables
      integer np
C Variables Returned
      double precision GS(np,np)
C Local Variables
      integer ia
      integer ip,jp
      double precision ans

      open(815,file='novlap.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ip=1,np
         do jp=1,np

C  Map the 2-index contracted integral to 1-D:
            call pack_2D(np,jp,ip,ia)

            read(815,REC=ia) ans
            GS(ip,jp)=ans

         end do
      end do

      close(815)



      return
      end

C=======================================================================
      subroutine read_elec_ovlap(ne,GS)
C=======================================================================
      implicit none
C Input Variables
      integer ne
C Variables Returned
      double precision GS(ne,ne)
C Local Variables
      integer ia
      integer ie1,je1
      double precision ans

      open(814,file='eovlap.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ie1=1,ne
         do je1=1,ne

C  Map the 2-index contracted integral to 1-D:
            call pack_2D(ne,je1,ie1,ia)

            read(814,REC=ia) ans 
            GS(ie1,je1)=ans

         end do
      end do

      close(814)



      return
      end
C=======================================================================
      subroutine read_elec_density(nebf,n_particle,DE)
C=======================================================================
      implicit none
C Input Variables
      integer nebf, n_particle
C Variables Returned
      double precision DE(nebf,nebf)
C Local Variables
      integer ia
      integer nocc,kstart,klast
      integer ie1,je1
      integer I,J,NMOS,IC,IMIN,IMAX,NUM1,JJ,ICC,NSTM,MODJ,MODIC
      integer K
      double precision ans
      double precision prefac
      double precision C(nebf,nebf)
      double precision VEC(nebf,nebf)


C--------FORMATTING----------------------------------------------------(
 9040 FORMAT(I2,I3,5E15.8)
 9060 FORMAT(1X,'*** ERROR IN READMO:   PROBLEM READING ORBITALS!'/
     *       1X,'POSSIBLY A DAMAGED OR MANGLED ORBITAL INPUT GROUP?'/
     *       1X,'ERROR OCCURED AT ORBITAL=',I6,' (MODULUS 100=',I4,'),'/
     *       1X,'         ITS LINE NUMBER=',I6/
     *       1X,'DATA READ FROM INPUT WAS ORBITAL=',I6,' LINE=',I6)
C--------FORMATTING----------------------------------------------------)

      open(850,file='guessCE.inp',status='unknown')

      NSTM=0

      NMOS=nebf
      NUM1=nebf

C-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------(
      DO 280 J = 1,NMOS
         IMAX = 0
         IC = 0
  240    CONTINUE
            IMIN = IMAX+1
            IMAX = IMAX+5
            IC = IC+1
            IF(IMAX .GT. NUM1) IMAX = NUM1
            READ(850,9040) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
     *                                        IMAX+NSTM)
c           READ(IR,9040,END=300,ERR=300) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
c    *                                                       IMAX+NSTM)
            MODJ  = MOD(J ,100 )
            MODIC = MOD(IC,1000)
            IF(JJ.EQ.MODJ . AND.  ICC.EQ.MODIC) GO TO 260
               WRITE(*,9060) J,MODJ,IC,JJ,ICC
               STOP
  260       CONTINUE
         IF(IMAX .LT. NUM1) GO TO 240
  280 CONTINUE
C-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------)

      close(850)

      C=VEC

C-----FORM-DENSITY-MATRIX----------------------------------------------(
c     if(i_particle.eq.1) then
      if(n_particle.gt.1) then
         nocc=n_particle/2
         prefac=2.0d+00
      else
         nocc=1
         prefac=1.0d+00
      end if
      kstart=1
      klast=nocc
c     end if

c     if(i_particle.eq.2) then
c        nocc=1
c        prefac=1.0d+00
c        kstart=NUCST
c        klast=NUCST
c     end if

      DO ie1=1,nebf
         DO je1=1,nebf

            DE(ie1,je1)=0.0d+00

            do k=kstart,klast
               DE(ie1,je1) = DE(ie1,je1) + prefac*C(ie1,k)*C(je1,k)
            end do

         END DO
      END DO
C-----FORM-DENSITY-MATRIX----------------------------------------------)




      return
      end

C=======================================================================
      subroutine read_nuc_density(npbf,n_particle,NUCST,DP)
C=======================================================================
      implicit none
C Input Variables
      integer npbf, n_particle,NUCST
C Variables Returned
      double precision DP(npbf,npbf)
C Local Variables
      integer ia
      integer nocc,kstart,klast
      integer ip,jp
      integer I,J,NMOS,IC,IMIN,IMAX,NUM1,JJ,ICC,NSTM,MODJ,MODIC
      integer K
      double precision ans
      double precision prefac
      double precision C(npbf,npbf)
      double precision VEC(npbf,npbf)


C--------FORMATTING----------------------------------------------------(
 9040 FORMAT(I2,I3,5E15.8)
 9060 FORMAT(1X,'*** ERROR IN NUC READMO: PROBLEM READING ORBITALS!'/
     *       1X,'POSSIBLY A DAMAGED OR MANGLED ORBITAL INPUT GROUP?'/
     *       1X,'ERROR OCCURED AT ORBITAL=',I6,' (MODULUS 100=',I4,'),'/
     *       1X,'         ITS LINE NUMBER=',I6/
     *       1X,'DATA READ FROM INPUT WAS ORBITAL=',I6,' LINE=',I6)
C--------FORMATTING----------------------------------------------------)

      open(851,file='guessCP.inp',status='unknown')

      NMOS=npbf
      NUM1=npbf
      NSTM=0

C-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------(
      DO 280 J = 1,NMOS
         IMAX = 0
         IC = 0
  240    CONTINUE
            IMIN = IMAX+1
            IMAX = IMAX+5
            IC = IC+1
            IF(IMAX .GT. NUM1) IMAX = NUM1
            READ(851,9040) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
     *                                        IMAX+NSTM)
c           READ(IR,9040,END=300,ERR=300) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
c    *                                                       IMAX+NSTM)
            MODJ  = MOD(J ,100 )
            MODIC = MOD(IC,1000)
            IF(JJ.EQ.MODJ . AND.  ICC.EQ.MODIC) GO TO 260
               WRITE(*,9060) J,MODJ,IC,JJ,ICC
               STOP
  260       CONTINUE
         IF(IMAX .LT. NUM1) GO TO 240
  280 CONTINUE
C-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------)

      close(851)

      C=VEC

C-----FORM-DENSITY-MATRIX----------------------------------------------(
c     if(i_particle.eq.1) then
c     if(n_particle.gt.1) then
c        nocc=n_particle/2
c        prefac=2.0d+00
c     else
c        nocc=1
c        prefac=1.0d+00
c     end if
c     kstart=1
c     klast=nocc
c     end if

c     if(i_particle.eq.2) then
         nocc=1
         prefac=1.0d+00
         kstart=NUCST
         klast=NUCST
c     end if

      DO ip=1,npbf
         DO jp=1,npbf

            DP(ip,jp)=0.0d+00

            do k=kstart,klast
               DP(ip,jp) = DP(ip,jp) + prefac*C(ip,k)*C(jp,k)
            end do

         END DO
      END DO
C-----FORM-DENSITY-MATRIX----------------------------------------------)




      return
      end

C=======================================================================
      subroutine write_MOs(IFIL,nbf,VEC)
C     IFIL=852 :: FinalCE.dat
C     IFIL=853 :: FinalCP.dat
C     IFIL=860 :: FinalCAE.dat
C     IFIL=861 :: FinalCBE.dat
C     IFIL=870 :: FinalCAalpE.dat
C     IFIL=871 :: FinalCAbetE.dat
C=======================================================================
      implicit none
C Input Variables
      integer IFIL,nbf
      double precision VEC(nbf,nbf)
C Variables Returned
C Local Variables

      call PUSQLF(IFIL,VEC,nbf,nbf,nbf)


      return
      end
C*MODULE MTHLIB  *DECK PUSQLF
      SUBROUTINE PUSQLF(LUFILE,V,M,N,LDV)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
c     LOGICAL GOPARR,DSKWRK,MASWRK
C
      DIMENSION V(LDV,M)
C
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ----- PUNCH A RECTANGULAR MATRIX WITH ORDERING LABELS -----
C     -V- IS -N- ROWS BY -M- COLUMNS, WITH TRUE LEAD DIMENSION -LDV-
C     NOTE THAT -PUSQLF- IS AN EXACT CLONE OF -PUSQL- EXCEPT A UNIT
C     NUMBER IS TO BE GIVEN AS AN ARGUMENT.
C
c     IF (.NOT.MASWRK) RETURN

      if(LUFILE.eq.852) then
         open(LUFILE,file='FinalCE.dat',status='unknown')
      else if(LUFILE.eq.853) then
         open(LUFILE,file='FinalCP.dat',status='unknown')
      else if(LUFILE.eq.860) then
         open(LUFILE,file='FinalCAE.dat',status='unknown')
      else if(LUFILE.eq.861) then
         open(LUFILE,file='FinalCBE.dat',status='unknown')
      else if(LUFILE.eq.870) then
         open(LUFILE,file='FinalCAalpE.dat',status='unknown')
      else if(LUFILE.eq.871) then
         open(LUFILE,file='FinalCAbetE.dat',status='unknown')
      end if
C
      DO J = 1,M
         IC = 0
         MAX = 0
  100    CONTINUE
            MIN = MAX+1
            MAX = MAX+5
            IC = IC+1
            IF (MAX .GT. N) MAX = N
            MODJ  = MOD(J ,100 )
            MODIC = MOD(IC,1000)
            WRITE (LUFILE,9008) MODJ,MODIC,(V(I,J),I = MIN,MAX)
         IF (MAX.LT.N) GO TO 100
      ENDDO
      close(LUFILE)
      RETURN
C
 9008 FORMAT(I2,I3,1P,5E15.8)
      END

C======================================================================
      subroutine prt_lower_triangle(nbf,nbfLT,mat)
C======================================================================
      implicit none
      integer nbf
      integer nbfLT
      double precision mat(nbf,nbf)
      double precision mat_LT(nbfLT)
      integer i
      integer j
      integer ia

C Pack Fock into lower triangle
      ia=0
      do i=1,nbf
         do j=1,i
            ia=ia+1
            mat_LT(ia)=mat(i,j)
         end do
      end do

c     CALL PRTRIL(fockLT,nbf)
      CALL PRTRI(mat_LT,nbf)


      return
      end

C*MODULE MTHLIB  *DECK PRTRI
      SUBROUTINE PRTRI(D,N)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      DIMENSION D(*)
C
c     COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
c     COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ----- PRINT SYMMETRIC MATRIX -D- OF DIMENSION -N- -----
C
c     IF (MASWRK) THEN
      MAX = 5
c     IF (NPRINT .EQ. 6) MAX = 10
      MM1 = MAX-1
      DO 120 I0=1,N,MAX
         IL = MIN(N,I0+MM1)
         WRITE(*,9008)
         WRITE(*,9028) (I,I=I0,IL)
         WRITE(*,9008)
         IL = -1
         DO 100 I=I0,N
            IL=IL+1
            J0=I0+(I*I-I)/2
            JL=J0+MIN(IL,MM1)
            WRITE(*,9048) I,(D(J),J=J0,JL)
  100    CONTINUE
  120 CONTINUE
c     END IF
      RETURN
 9008 FORMAT(1X)
 9028 FORMAT(6X,10(4X,I4,4X))
 9048 FORMAT(I5,1X,10F12.7)
      END

c     SUBROUTINE PREVNU(V,E,M,LDV,ISTART,IEND)
      SUBROUTINE PREVNU(V,E,M,LDV,NBF)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
c     LOGICAL GOPARR,DSKWRK,MASWRK
c
c     parameter (MXAO=8192)
C
c     CHARACTER*8 QNUN,QNN
c     CHARACTER*10 PBFLAB
C
      DIMENSION V(LDV,M),E(M)
C
c     COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
c     COMMON /NUCMON/ QNUN(20),QNN(20),PBFLAB(MXAO)
c     COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ----- PRINT OUT EIGENDATA (VECTORS AND VALUES) -----
C     THE ROWS WILL BE LABELED WITH THE BASIS FUNCTION TAGS.
C     -V- IS N X M, WITH TRUE LEADING DIMENSION -LDV-
C
c     IF (MASWRK) THEN
      MAX = 5
c     IF (NPRINT .EQ. 6) MAX = 10
      IMAX = 0
C
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      WRITE (*,9008)
      WRITE (*,9028) (I,I = IMIN,IMAX)
      WRITE (*,9008)
      WRITE (*,9068) (E(I),I = IMIN,IMAX)
      WRITE (*,9008)
      J = 0
c     DO 120 IDX = ISTART,IEND
      DO 120 IDX = 1,NBF
         J = J + 1
c        WRITE (IW,9048) J,PBFLAB(IDX),(V(J,I),I = IMIN,IMAX)
         WRITE (*,9048) J,(V(J,I),I = IMIN,IMAX)
  120 CONTINUE
      IF (IMAX .LT. M) GO TO 100
c     ENDIF
      RETURN
C
 9008 FORMAT(1X)
 9028 FORMAT(17X,10(4X,I4,3X))
c9048 FORMAT(I5,2X,A10,10F11.6)
 9048 FORMAT(I5,2X,'LABEL     ',10F11.6)
 9068 FORMAT(17X,10F11.4)
      END

