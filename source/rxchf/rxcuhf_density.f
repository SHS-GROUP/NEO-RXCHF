!=======================================================================
      subroutine RXCUHF_construct_DAE(NAE,nebf,CAE,DAE)

! Form alpha (or beta) electronic density matrix 
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

      DO i=1,nebf
         DO j=1,nebf

          DAE(j,i)=0.0d+00
          do k=1,nae
            DAE(j,i) = DAE(j,i) + CAE(j,k)*CAE(i,k)
          end do

         END DO
      END DO


      return
      end
!=======================================================================
      subroutine RXCUHF_read_CAE(fname,nebf,NAE,DAE,C)
!=======================================================================
      implicit none
! Input Variables
      integer nebf,NAE
      character*15 fname
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

      write(*,*) "Reading orbitals from ",fname
      open(862,file=fname,status='unknown')

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

      DO ie1=1,nebf
         DO je1=1,nebf

            DAE(je1,ie1)=0.0d+00

            do k=1,nae
               DAE(je1,ie1) = DAE(je1,ie1) + C(je1,k)*C(ie1,k)
            end do

         END DO
      END DO
!-----FORM-DENSITY-MATRIX----------------------------------------------)


      return
      end
