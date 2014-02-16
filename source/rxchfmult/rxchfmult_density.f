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
      subroutine RXCHFmult_read_CAE(nebf,NAE,DAE,C)
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
      subroutine RXCHFmult_read_CBE(nebf,NBE,DBE,C)
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
      if(NBE.gt.1) then
       coeff=two
       nocc=nbe/2
      else
       coeff=one
       nocc=1
      end if

      DO ie1=1,nebf
         DO je1=1,nebf

            DBE(je1,ie1)=0.0d+00

            do k=1,nocc
               DBE(je1,ie1) = DBE(je1,ie1) + coeff*C(je1,k)*C(ie1,k)
            end do

         END DO
      END DO
!-----FORM-DENSITY-MATRIX----------------------------------------------)


      return
      end
!=======================================================================
      subroutine RXCHFmult_read_CE(nebf,nebfBE,elindBE,xxse,xxseBE,
     x                             NAE,DAE,CAE,NBE,DBE,CBE)
!=======================================================================
      implicit none
! Input Variables
      integer nebf,NAE
      integer nebfBE,NBE          ! num bfs in special elec basis
      integer elindBE(nebfBE)     ! Contracted indices of NBE basis set
      double precision xxse(nebf,nebf)
      double precision xxseBE(nebfBE,nebfBE)
! Variables Returned
      double precision DAE(nebf,nebf)
      double precision DBE(nebfBE,nebfBE)
      double precision CAE(nebf,nebf)
      double precision CBE(nebfBE,nebfBE)
! Local Variables
      integer ia
      integer nocca,noccb
      integer ie1,je1
      integer I,J,NMOS,IC,IMIN,IMAX,NUM1,JJ,ICC,NSTM,MODJ,MODIC
      integer k,l
      integer contrind,currcontrind
      logical laddbasis
      integer dimint
      integer currind,maxind
      double precision ans
      double precision ovlap,maxovlap
      double precision coeffa,coeffb
      double precision VECA(nebf,nebf)
      double precision VECB(nebf,nebf)
      double precision Ctemp(nebfBE,nebf)

      double precision,allocatable :: Cint(:,:)

      double precision zero,one,two,tol
      parameter(zero=0.0d+00,one=1.0d+00,two=2.0d+00,tol=1.0d-12)


!--------FORMATTING----------------------------------------------------(
 9040 FORMAT(I2,I3,5E15.8)
 9060 FORMAT(1X,'*** ERROR IN READ_CE:    PROBLEM READING ORBITALS!'/
     *       1X,'POSSIBLY A DAMAGED OR MANGLED ORBITAL INPUT GROUP?'/
     *       1X,'ERROR OCCURED AT ORBITAL=',I6,' (MODULUS 100=',I4,'),'/
     *       1X,'         ITS LINE NUMBER=',I6/
     *       1X,'DATA READ FROM INPUT WAS ORBITAL=',I6,' LINE=',I6)
!--------FORMATTING----------------------------------------------------)

      if(NAE.gt.1) then
       coeffa=two
       nocca=nae/2
      else
       coeffa=one
       nocca=1
      end if
      CAE=zero

      if(NBE.gt.1) then
       coeffb=two
       noccb=nbe/2
      else
       coeffb=one
       noccb=1
      end if
      CBE=zero

C Read in regular electron guess
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
            READ(862,9040) JJ,ICC,(VECA(I,J),I=IMIN+NSTM,
     *                                         IMAX+NSTM)
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

C Reorder electronic basis set such that special electron subset is first
      currcontrind=1

C First add bfs that are also in special electronic set
      do i=1,nebfBE
        contrind=elindBE(i)
        do j=1,nebf
          CAE(currcontrind,j)=VECA(contrind,j)
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
           CAE(currcontrind,j)=VECA(i,j)
         end do
         currcontrind=currcontrind+1
        end if
      end do

C Read in special electron guess
      open(863,file='guessCBE.inp',status='unknown')

      NSTM=0

      NMOS=nebfBE
      NUM1=nebfBE

!-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------(
      DO 380 J = 1,NMOS
         IMAX = 0
         IC = 0
  340    CONTINUE
            IMIN = IMAX+1
            IMAX = IMAX+5
            IC = IC+1
            IF(IMAX .GT. NUM1) IMAX = NUM1
            READ(863,9040) JJ,ICC,(VECB(I,J),I=IMIN+NSTM,
     *                                         IMAX+NSTM)
            MODJ  = MOD(J ,100 )
            MODIC = MOD(IC,1000)
            IF(JJ.EQ.MODJ . AND.  ICC.EQ.MODIC) GO TO 360
               WRITE(*,9060) J,MODJ,IC,JJ,ICC
               STOP
  360       CONTINUE
         IF(IMAX .LT. NUM1) GO TO 340
  380 CONTINUE
!-----READ-GAMESS-STYLE-VEC-GROUP--------------------------------------)

      close(863)

C Find intersection of occupied special space and virtual regular space
      call RXCHFmult_intersection(nebf,nebf-nocca,CAE(:,nocca+1:nebf),
     x                            nebfBE,xxseBE,dimint,Ctemp)
      if (dimint.le.noccb) then
       write(*,*) "dim of intersection <= # occ special orbitals"
       write(*,*) "dimint:",dimint
       write(*,*) "noccb:",noccb
       return
      end if

      write(*,*)
      write(*,*) "----------------------------------"
      write(*,*) " Dimension of intersection space:"
      write(*,'(2X,A,1X,I3)') "Max possible (nebfBE)      =",nebfBE
      write(*,'(2X,A,1X,I3)') "Actual (after computation) =",dimint
      write(*,*) "----------------------------------"
      write(*,*)

      if(allocated(Cint)) deallocate(Cint)
      allocate(Cint(nebfBE,dimint))
      do i=1,dimint
        do j=1,nebfBE
          Cint(j,i)=Ctemp(j,i)
        end do
      end do

! Debug: All CBE vectors should be orthogonal to occ CAE vectors
       do i=1,nocca
       do j=1,dimint
         ovlap=zero
         do k=1,nebf
         do l=1,nebfBE
           ovlap=ovlap+CAE(k,i)*Cint(l,j)*xxse(k,l)
         end do
         end do
         if (abs(ovlap).gt.tol) then
          write(*,*)
          write(*,*) "******* ERROR *******"
          write(*,*) "Calculated intersection basis is not orthogonal"
          write(*,*) "to occupied regular vectors"
          write(*,*) "reg occ index, int index, ovlap:",i,j,ovlap
          write(*,*)
         end if
       end do
       end do

! Check orthonormality of new basis
C       if(debug) then
        write(*,*) "Overlaps amongst Cint bfs (should be onormal):"
        do i=1,dimint
        do j=1,dimint
          ovlap=zero
          do k=1,nebfBE
          do l=1,nebfBE
            ovlap=ovlap+Cint(k,i)*Cint(l,j)*xxseBE(k,l)
          end do
          end do
          write(*,*) "i,j,ovlap:",i,j,ovlap
        end do
        end do
C       end if

! Fill in regular virtuals with new intersection basis
       do i=1,min(nebf-nocca,dimint)
         do j=1,nebfBE
           CAE(j,i+nocca)=Cint(j,i)
         end do
       end do

! Find how "similar" the original CBE are as compared to the new basis
! and fill special electronic vectors with those with greatest overlap
       currind=0
       do i=1,nebfBE
         maxovlap=zero
         maxind=0
         do j=1,dimint
           ovlap=zero
           do k=1,nebfBE
           do l=1,nebfBE
             ovlap=ovlap+VECB(k,i)*Cint(l,j)*xxseBE(k,l)
           end do
           end do
           write(*,*) "overlap for orig MO",i,"with new MO",j,
     x                "is:",ovlap
           if(abs(ovlap).gt.maxovlap) then
            maxovlap=abs(ovlap)
            maxind=j
           end if
         end do
         write(*,*) "Max ovlap for orig MO",i,"is with new MO",maxind,
     x              "with overlap:",maxovlap
         if (maxind.ne.0) then
          currind=currind+1
          do j=1,nebfBE
            CBE(j,currind)=Cint(j,maxind)
          end do
         end if
       end do
       
! Debug: All CBE vectors should be orthogonal to occ CAE vectors
       do i=1,nocca
       do j=1,dimint
         ovlap=zero
         do k=1,nebf
         do l=1,nebfBE
           ovlap=ovlap+CAE(k,i)*CBE(l,j)*xxse(k,l)
         end do
         end do
         if (abs(ovlap).gt.tol) then
          write(*,*)
          write(*,*) "******* ERROR *******"
          write(*,*) "Calculated intersection basis is not orthogonal"
          write(*,*) "to occupied regular vectors"
          write(*,*) "reg occ index, int index, ovlap:",i,j,ovlap
          write(*,*)
         end if
       end do
       end do

      if(allocated(Cint)) deallocate(Cint)

C Form density matrices
      DO ie1=1,nebf
      DO je1=1,nebf

         DAE(je1,ie1)=zero

         do k=1,nocca
            DAE(je1,ie1) = DAE(je1,ie1) + coeffa*CAE(je1,k)*CAE(ie1,k)
         end do

      END DO
      END DO

      DO ie1=1,nebfBE
      DO je1=1,nebfBE

         DBE(je1,ie1)=zero

         do k=1,noccb
            DBE(je1,ie1) = DBE(je1,ie1) + coeffb*CBE(je1,k)*CBE(ie1,k)
         end do

      END DO
      END DO

      return
      end
!=======================================================================
      subroutine RXCHFmult_read_CAE_ocbse3(nebf,nebfBE,elindBE,xxse,
     x                                     NAE,DAE,CAE,NBE,CBE)
!=======================================================================
      implicit none
! Input Variables
      integer nebf,NAE
      integer nebfBE,NBE          ! num bfs in special elec basis
      integer elindBE(nebfBE)     ! Contracted indices of NBE basis set
      double precision xxse(nebf,nebf)
      double precision CBE(nebfBE,nebfBE)
! Variables Returned
      double precision DAE(nebf,nebf)
      double precision CAE(nebf,nebf)
! Local Variables
      integer ia
      integer nocca,noccb
      integer ie1,je1
      integer I,J,NMOS,IC,IMIN,IMAX,NUM1,JJ,ICC,NSTM,MODJ,MODIC
      integer k,l
      integer contrind,currcontrind
      logical laddbasis
      integer dimint
      integer currind,maxind
      double precision ans
      double precision ovlap,maxovlap
      double precision coeffa
      double precision VECA(nebf,nebf)
      double precision Ctemp(nebf,nebf)

      double precision zero,one,two,tol
      parameter(zero=0.0d+00,one=1.0d+00,two=2.0d+00,tol=1.0d-10)


!--------FORMATTING----------------------------------------------------(
 9040 FORMAT(I2,I3,5E15.8)
 9060 FORMAT(1X,'*** ERROR IN READ_CE:    PROBLEM READING ORBITALS!'/
     *       1X,'POSSIBLY A DAMAGED OR MANGLED ORBITAL INPUT GROUP?'/
     *       1X,'ERROR OCCURED AT ORBITAL=',I6,' (MODULUS 100=',I4,'),'/
     *       1X,'         ITS LINE NUMBER=',I6/
     *       1X,'DATA READ FROM INPUT WAS ORBITAL=',I6,' LINE=',I6)
!--------FORMATTING----------------------------------------------------)

      if(NBE.gt.1) then
       noccb=nbe/2
      else
       noccb=1
      end if

      if(NAE.gt.1) then
       coeffa=two
       nocca=nae/2
      else
       coeffa=one
       nocca=1
      end if
      CAE=zero

C Read in regular electron guess
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
            READ(862,9040) JJ,ICC,(VECA(I,J),I=IMIN+NSTM,
     *                                         IMAX+NSTM)
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

C Reorder electronic basis set such that special electron subset is first
      currcontrind=1

C First add bfs that are also in special electronic set
      do i=1,nebfBE
        contrind=elindBE(i)
        do j=1,nebf
          CAE(currcontrind,j)=VECA(contrind,j)
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
           CAE(currcontrind,j)=VECA(i,j)
         end do
         currcontrind=currcontrind+1
        end if
      end do

      Ctemp=zero
      do i=1,nebfBE
        do k=1,nebfBE
          Ctemp(k,i)=CBE(k,i)
        end do
      end do

C Replace nebfBE guess regular virtual MOs with closest overlapping
C special electron MOs to ensure linear independence
      do i=1,nebfBE
        maxovlap=zero
        do j=nocca+1,nebf
          call moovlap(nebf,CAE(:,j),Ctemp(:,i),xxse,ovlap)
          if (abs(ovlap).ge.maxovlap) then
           maxind=j
           maxovlap=abs(ovlap)
          end if
        end do
        CAE(:,maxind)=zero
        do k=1,nebfBE
          CAE(k,maxind)=Ctemp(k,i)
        end do
      end do

C G-S orthogonalize
      do i=2,nebf
        do j=i,1,-1
          call moovlap(nebf,CAE(:,i),CAE(:,j),xxse,ovlap)
          do k=1,nebf
            CAE(k,i)=CAE(k,i)-ovlap*CAE(k,j)
          end do
        end do
        call moovlap(nebf,CAE(:,i),CAE(:,i),xxse,ovlap)
        do k=1,nebf
          CAE(k,i)=CAE(k,i)/dsqrt(ovlap)
        end do
      end do

! Check orthogonality
      do i=1,nebf
        do j=1,nebf
          call moovlap(nebf,CAE(:,i),CAE(:,j),xxse,ovlap)
          if ((abs(ovlap).gt.tol).and.(i.ne.j)) then
           write(*,*) "G-S procedure yielded MOs:",i,j,
     x                "that are not orthogonal:",ovlap
          end if
        end do
      end do

C Form density matrices
      DO ie1=1,nebf
      DO je1=1,nebf

         DAE(je1,ie1)=zero

         do k=1,nocca
            DAE(je1,ie1) = DAE(je1,ie1) + coeffa*CAE(je1,k)*CAE(ie1,k)
         end do

      END DO
      END DO

      return
      end


