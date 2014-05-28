!=======================================================================
      subroutine RXCHFmult_construct_DE(nelec,nebf,CE,DE)

! Form electronic density matrix 
!=======================================================================
      implicit none

! Input Variables
      integer nelec,nebf
      double precision CE(nebf,nebf)

! Variables Returned
      double precision DE(nebf,nebf)

! Local Variables
      integer i,j,k
      integer nocc
      double precision coeff
      double precision one,two
      parameter(one=1.0d+00,two=2.0d+00)

      if(nelec.gt.1) then
       coeff=two
       nocc=nelec/2
      else
       coeff=one
       nocc=1
      end if

      DO i=1,nebf
         DO j=1,nebf

          DE(j,i)=0.0d+00
          do k=1,nocc
            DE(j,i) = DE(j,i) + coeff*CE(j,k)*CE(i,k)
          end do

         END DO
      END DO


      return
      end
!=======================================================================
      subroutine RXCHFmult_read_CAE(nebf,nebfBE,NAE,LALTBAS,elindBE,
     x                              LTCSCF,regmos_1,regmos_2,
     x                              DAE,C)
!=======================================================================
      implicit none
! Input Variables
      integer nebf,NAE
      integer nebfBE              ! num bfs in special elec basis
      logical LALTBAS             ! flag for distinct special elec basis
      integer elindBE(nebfBE)     ! Contracted indices of NBE basis set
      logical LTCSCF
      integer regmos_1(nae)
      integer regmos_2(nae)
! Variables Returned
      double precision DAE(nebf,nebf)
! Local Variables
      integer ia
      integer nocc
      integer ie1,je1
      integer I,J,NMOS,IC,IMIN,IMAX,NUM1,JJ,ICC,NSTM,MODJ,MODIC
      integer k
      integer contrind,currcontrind
      logical laddbasis
      double precision ans
      double precision coeff
      double precision C(nebf,nebf)
      double precision VEC(nebf,nebf)
      double precision one,two
      parameter(one=1.0d+00,two=2.0d+00)


!--------FORMATTING----------------------------------------------------(
 9040 FORMAT(I2,I3,5E15.8)
 9060 FORMAT(1X,'*** ERROR IN READ_CAE:   PROBLEM READING ORBITALS!'/
     *       1X,'POSSIBLY A DAMAGED OR MANGLED ORBITAL INPUT GROUP?'/
     *       1X,'ERROR OCCURED AT ORBITAL=',I6,' (MODULUS 100=',I4,'),'/
     *       1X,'         ITS LINE NUMBER=',I6/
     *       1X,'DATA READ FROM INPUT WAS ORBITAL=',I6,' LINE=',I6)
!--------FORMATTING----------------------------------------------------)

      open(862,file='guessCAE.inp',status='unknown')

      NSTM=0

      NMOS=nebf
      NUM1=nebf

!-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------(
      DO 280 J = 1,NMOS
         IMAX = 0
         IC = 0
  240    CONTINUE
            IMIN = IMAX+1
            IMAX = IMAX+5
            IC = IC+1
            IF(IMAX .GT. NUM1) IMAX = NUM1
            READ(862,9040) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
     *                                        IMAX+NSTM)
!           READ(IR,9040,END=300,ERR=300) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
!    *                                                       IMAX+NSTM)
            MODJ  = MOD(J ,100 )
            MODIC = MOD(IC,1000)
            IF(JJ.EQ.MODJ . AND.  ICC.EQ.MODIC) GO TO 260
               WRITE(*,9060) J,MODJ,IC,JJ,ICC
               STOP
  260       CONTINUE
         IF(IMAX .LT. NUM1) GO TO 240
  280 CONTINUE
!-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------)

      close(862)

      if (LALTBAS) then

C Reorder electronic basis set such that special electron subset is first
       currcontrind=1

C First add bfs that are also in special electronic set
       do i=1,nebfBE
         contrind=elindBE(i)
         do j=1,nebf
           C(currcontrind,j)=VEC(contrind,j)
         end do
         currcontrind=currcontrind+1
       end do

C Add remaining bfs
       do i=1,nebf
         laddbasis=.true.
         do j=1,nebfBE
           contrind=elindBE(j)
           if (i.eq.contrind) laddbasis=.false.
         end do
         if (laddbasis) then
          do j=1,nebf
            C(currcontrind,j)=VEC(i,j)
          end do
          currcontrind=currcontrind+1
         end if
       end do

      else

       C=VEC

      end if

!-----FORM-DENSITY-MATRIX----------------------------------------------(
      if(LTCSCF) then
       call TCSCF_construct_DAE(NAE,nebf,regmos_1,regmos_2,C,DAE)
      else
       call RXCHFmult_construct_DE(NAE,nebf,C,DAE)
      end if
!-----FORM-DENSITY-MATRIX----------------------------------------------)


      return
      end
!=======================================================================
      subroutine RXCHFmult_read_CBE(nebf,NBE,
     x                              LTCSCF,spemos_1,spemos_2,
     x                              DBE,C)
!=======================================================================
      implicit none
! Input Variables
      integer nebf,NBE
      logical LTCSCF
      integer spemos_1(nbe)
      integer spemos_2(nbe)
! Variables Returned
      double precision DBE(nebf,nebf)
! Local Variables
      integer ia
      integer nocc
      integer ie1,je1
      integer I,J,NMOS,IC,IMIN,IMAX,NUM1,JJ,ICC,NSTM,MODJ,MODIC
      integer K
      double precision ans
      double precision coeff
      double precision C(nebf,nebf)
      double precision VEC(nebf,nebf)
      double precision one,two
      parameter(one=1.0d+00,two=2.0d+00)


!--------FORMATTING----------------------------------------------------(
 9040 FORMAT(I2,I3,5E15.8)
 9060 FORMAT(1X,'*** ERROR IN READ_CBE:   PROBLEM READING ORBITALS!'/
     *       1X,'POSSIBLY A DAMAGED OR MANGLED ORBITAL INPUT GROUP?'/
     *       1X,'ERROR OCCURED AT ORBITAL=',I6,' (MODULUS 100=',I4,'),'/
     *       1X,'         ITS LINE NUMBER=',I6/
     *       1X,'DATA READ FROM INPUT WAS ORBITAL=',I6,' LINE=',I6)
!--------FORMATTING----------------------------------------------------)

      open(863,file='guessCBE.inp',status='unknown')

      NSTM=0

      NMOS=nebf
      NUM1=nebf

!-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------(
      DO 280 J = 1,NMOS
         IMAX = 0
         IC = 0
  240    CONTINUE
            IMIN = IMAX+1
            IMAX = IMAX+5
            IC = IC+1
            IF(IMAX .GT. NUM1) IMAX = NUM1
            READ(863,9040) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
     *                                        IMAX+NSTM)
!           READ(IR,9040,END=300,ERR=300) JJ,ICC,(VEC(I,J),I=IMIN+NSTM,
!    *                                                       IMAX+NSTM)
            MODJ  = MOD(J ,100 )
            MODIC = MOD(IC,1000)
            IF(JJ.EQ.MODJ . AND.  ICC.EQ.MODIC) GO TO 260
               WRITE(*,9060) J,MODJ,IC,JJ,ICC
               STOP
  260       CONTINUE
         IF(IMAX .LT. NUM1) GO TO 240
  280 CONTINUE
!-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------)

      close(863)

      C=VEC

!-----FORM-DENSITY-MATRIX----------------------------------------------(
      if(LTCSCF) then
       call TCSCF_construct_DBE(NBE,nebf,spemos_1,spemos_2,C,DBE)
      else
       call RXCHFmult_construct_DE(NBE,nebf,C,DBE)
      end if
!-----FORM-DENSITY-MATRIX----------------------------------------------)


      return
      end
C=======================================================================
      subroutine TCSCF_read_nuc_density(npbf,n_particle,NUCST,DP)
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
      call TCSCF_construct_DP(npbf,C,DP)
C-----FORM-DENSITY-MATRIX----------------------------------------------)

      return
      end
!=======================================================================
      subroutine TCSCF_construct_DP(npbf,CP,DP)

! Form TC-SCF nuclear density matrix 
!=======================================================================
      implicit none

! Input Variables
      integer npbf
      double precision CP(npbf,npbf)

! Variables Returned
      double precision DP(npbf,npbf)

! Local Variables
      integer ip,jp
      double precision half
      parameter(half=0.50d+00)

C Convention:
C   configuration 1 has MO 1 occupied and MO 2 unoccupied
C   configuration 2 has MO 2 occupied and MO 1 unoccupied
C
C *********** PROBABLY HAS TO BE CHANGED ************
C
      DP=0.0d+00

      DO ip=1,npbf
         DO jp=1,npbf

           DP(ip,jp) = half*(CP(ip,1)*CP(jp,1)+CP(ip,2)*CP(jp,2))

         END DO
      END DO

      return
      end
!=======================================================================
      subroutine TCSCF_construct_DAE(nelec,nebf,mos_1,mos_2,CE,DE)

! Form TCSCF regular electronic density matrix 
!=======================================================================
      implicit none

! Input Variables
      integer nelec,nebf
      double precision CE(nebf,nebf)
      integer mos_1(nelec)
      integer mos_2(nelec)

! Variables Returned
      double precision DE(nebf,nebf)

! Local Variables
      integer i,j,k
      double precision half
      parameter(half=0.50d+00)

      DO i=1,nebf
         DO j=1,nebf

          DE(j,i)=0.0d+00
          do k=1,nelec
            DE(j,i) = DE(j,i) + half*
     x              ( CE(j,mos_1(k))*CE(i,mos_1(k))
     x              + CE(j,mos_2(k))*CE(i,mos_2(k)) )
          end do

         END DO
      END DO


      return
      end
!=======================================================================
      subroutine TCSCF_construct_DBE(nelec,nebf,mos_1,mos_2,CE,DE)

! Form TCSCF special electronic density matrix 
!=======================================================================
      implicit none

! Input Variables
      integer nelec,nebf
      double precision CE(nebf,nebf)
      integer mos_1(nelec)
      integer mos_2(nelec)

! Variables Returned
      double precision DE(nebf,nebf)

! Local Variables
      integer i,j,k
      double precision half
      parameter(half=0.50d+00)

      DO i=1,nebf
         DO j=1,nebf

          DE(j,i)=0.0d+00
          do k=1,nelec
            DE(j,i) = DE(j,i) + half*
     x              ( CE(j,mos_1(k))*CE(i,mos_1(k))
     x              + CE(j,mos_2(k))*CE(i,mos_2(k)) )
          end do

         END DO
      END DO


      return
      end

