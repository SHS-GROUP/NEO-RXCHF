!=======================================================================
      subroutine RXCHF_construct_DAE(NAE,nebf,CAE,DAE)

! Form regular electronic density matrix 
!=======================================================================
      implicit none

! Input Variables
      integer NAE,nebf
      double precision CAE(nebf,nebf)

! Variables Returned
      double precision DAE(nebf,nebf)

! Local Variables
      integer i,j,k
      integer nocc
      double precision coeff
      double precision one,two
      parameter(one=1.0d+00,two=2.0d+00)

      if(NAE.gt.1) then
       coeff=two
       nocc=nae/2
      else
       coeff=one
       nocc=1
      end if

      DO i=1,nebf
         DO j=1,nebf

          DAE(j,i)=0.0d+00
          do k=1,nocc
            DAE(j,i) = DAE(j,i) + coeff*CAE(j,k)*CAE(i,k)
          end do

         END DO
      END DO


      return
      end
!=======================================================================
      subroutine RXCHF_construct_DBE(NBE,nebf,CBE,DBE)

! Form special electronic density matrix 
!=======================================================================
      implicit none

! Input Variables
      integer NBE,nebf
      double precision CBE(nebf,nebf)

! Variables Returned
      double precision DBE(nebf,nebf)

! Local Variables
      integer i,j,k
      integer nocc

      nocc=1

      DO i=1,nebf
         DO j=1,nebf

          DBE(j,i)=0.0d+00
          do k=1,nocc
            DBE(j,i) = DBE(j,i) + CBE(j,k)*CBE(i,k)
          end do

         END DO
      END DO


      return
      end
!=======================================================================
      subroutine RXCHF_read_CAE(nebf,NAE,DAE,C)
!=======================================================================
      implicit none
! Input Variables
      integer nebf,NAE
! Variables Returned
      double precision DAE(nebf,nebf)
! Local Variables
      integer ia
      integer nocc
      integer ie1,je1
      integer I,J,NMOS,IC,IMIN,IMAX,NUM1,JJ,ICC,NSTM,MODJ,MODIC
      integer k
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

      C=VEC

!-----FORM-DENSITY-MATRIX----------------------------------------------(
      if(NAE.gt.1) then
       coeff=two
       nocc=nae/2
      else
       coeff=one
       nocc=1
      end if

      DO ie1=1,nebf
         DO je1=1,nebf

            DAE(je1,ie1)=0.0d+00

            do k=1,nocc
               DAE(je1,ie1) = DAE(je1,ie1) + coeff*C(je1,k)*C(ie1,k)
            end do

         END DO
      END DO
!-----FORM-DENSITY-MATRIX----------------------------------------------)


      return
      end
!=======================================================================
      subroutine RXCHF_read_CBE(nebf,NBE,DBE,C)
!=======================================================================
      implicit none
! Input Variables
      integer nebf,NBE
! Variables Returned
      double precision DBE(nebf,nebf)
! Local Variables
      integer ia
      integer nocc
      integer ie1,je1
      integer I,J,NMOS,IC,IMIN,IMAX,NUM1,JJ,ICC,NSTM,MODJ,MODIC
      integer K
      double precision ans
      double precision prefac
      double precision C(nebf,nebf)
      double precision VEC(nebf,nebf)


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
      nocc=1 

      DO ie1=1,nebf
         DO je1=1,nebf

            DBE(je1,ie1)=0.0d+00

            do k=1,nocc
               DBE(je1,ie1) = DBE(je1,ie1) + C(je1,k)*C(ie1,k)
            end do

         END DO
      END DO
!-----FORM-DENSITY-MATRIX----------------------------------------------)


      return
      end
