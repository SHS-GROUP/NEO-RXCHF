C=======================================================================
      subroutine OrthoABS
C  
C  Main driver for orthonormalization of electron and proton ABS
C Orthonormalize the ABS
C => Normalize the ABS
C => Calculate Overlap Matrix
C => Calculate transformation matrix X
C=======================================================================
      implicit none
C Input Variables
c     integer
c     double precision
C Variables Returned
C Local Variables
      integer NauxBFe
      integer NauxBFp

      write(*,*)
      write(*,*) 'Orthogonalization of ABS...'
      write(*,*)

      open(unit=911,file='aux_e_basis.inp',status='unknown')
      read(911,*) NauxBFe
c     write(*,*) 'In OrthoABS:  NauxBFe=',NauxBFe
      close(911)
      open(unit=912,file='aux_p_basis.inp',status='unknown')
      read(912,*) NauxBFp
      close(912)

      call norm_abs_e
      call Smat_abs_e
      call ortabs_e(NauxBFe)
c     call test_ort_e

      call norm_abs_p
      call Smat_abs_p
      call ortabs_p(NauxBFp)
c     call test_ort_p
c     stop

      call drive_contract_X(NauxBFe,NauxBFp)

      write(*,*)
      write(*,*) '**Orthogonalization of ABS completed.**'
      write(*,*)

      return
      end

C======================================================================
      subroutine drive_contract_X(NauxBFe,NauxBFp)
C======================================================================
      implicit none
C Input Variables
      integer NauxBFe
      integer NauxBFp
C Local Variables
      integer a,b,ia
      double precision Mab


C  Open transformation matrix product file for e ABS
      open(919,file='ABSME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do a=1,NauxBFe
         do b=1,NauxBFe
            call contract_XE(a,b,NauxBFe,Mab)
            call pack_2D(NauxBFE,a,b,ia)
            write(919,REC=ia) Mab
         end do
      end do
      close(919)


C  Open transformation matrix product file for p ABS
      open(920,file='ABSMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do a=1,NauxBFp
         do b=1,NauxBFp
            call contract_XP(a,b,NauxBFP,Mab)
            call pack_2D(NauxBFP,a,b,ia)
            write(920,REC=ia) Mab
         end do
      end do
      close(920)


      return
      end

C======================================================================
      subroutine contract_XE(a,b,NauxBF,XX)

C======================================================================
      implicit none

C Input Variables
      integer a
      integer b
      integer NauxBF
C Variables Returned
      double precision XX
C Local Variables
      integer i
      integer ia
      double precision X_ai
      double precision X_bi
      double precision ans

C  Open transformation matrix file for e ABS
      open(917,file='ABSXE.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      XX=0.0d+00

      do i=1,NauxBF

          call pack_2D(NauxBF,a,i,ia)
          read(917,REC=ia) X_ai
          call pack_2D(NauxBF,b,i,ia)
          read(917,REC=ia) X_bi

          ans=X_ai*X_bi

          XX=XX+ans

      end do


      return
      end

C======================================================================
      subroutine contract_XP(a,b,NauxBF,XX)

C======================================================================
      implicit none

C Input Variables
      integer a
      integer b
      integer NauxBF
C Variables Returned
      double precision XX
C Local Variables
      integer i
      integer ia
      double precision X_ai
      double precision X_bi
      double precision ans

C  Open transformation matrix file for e ABS
      open(918,file='ABSXP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      XX=0.0d+00

      do i=1,NauxBF

          call pack_2D(NauxBF,a,i,ia)
          read(918,REC=ia) X_ai
          call pack_2D(NauxBF,b,i,ia)
          read(918,REC=ia) X_bi

          ans=X_ai*X_bi

          XX=XX+ans

      end do


      return
      end

C=======================================================================
      subroutine norm_abs_e
C  Normalize the ABS for electrons
C=======================================================================
      implicit none
C Input Variables
c     integer 
c     double precision
C Variables retuned
c     integer 
c     double precision
C Local Variables
      integer k 
      integer NauxBF 
      integer LX,MX,NX
      double precision ans
      double precision abs_norm_e
      double precision BX
      double precision BmatX(3)

C  Open file to store normalization factors for e ABS
      open(913,file='ABSNORME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
c  Read RI auxilliary basis input file
      open(unit=911,file='aux_e_basis.inp',status='unknown')
      read(911,*) NauxBF
 
      do k=1,NauxBF

         read(911,*) LX,MX,NX,BX,BmatX(1),BmatX(2),BmatX(3)

         call gfovlap(LX,MX,NX,BX,BmatX,
     2                LX,MX,NX,BX,BmatX,ans)

         abs_norm_e=1.0d+00/sqrt(ans)

c        write(*,*)'ans=',ans
c        write(*,*)'abs_norm_e=',abs_norm_e
c        write(*,*)
         write(913,REC=k) abs_norm_e

      end do

      close(911)
      close(913)

      return
      end

C=======================================================================
      subroutine norm_abs_p
C  Normalize the ABS for protons
C=======================================================================
      implicit none
C Input Variables
c     integer 
c     double precision
C Variables retuned
c     integer 
c     double precision
C Local Variables
      integer k 
      integer NauxBF 
      integer LX,MX,NX
      double precision ans
      double precision abs_norm_p
      double precision BX
      double precision BmatX(3)

C  Open file to store normalization factors for p ABS
      open(914,file='ABSNORMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
c  Read RI auxilliary basis input file
      open(unit=912,file='aux_p_basis.inp',status='unknown')
      read(912,*) NauxBF
 
      do k=1,NauxBF

         read(912,*) LX,MX,NX,BX,BmatX(1),BmatX(2),BmatX(3)

         call gfovlap(LX,MX,NX,BX,BmatX,
     2                LX,MX,NX,BX,BmatX,ans)

         abs_norm_p=1.0d+00/sqrt(ans)
c        write(*,*)'ans=',ans
c        write(*,*)'abs_norm_p=',abs_norm_p
c        write(*,*)

         write(914,REC=k) abs_norm_p

      end do

      close(912)
      close(914)

      return
      end

C=======================================================================
      subroutine Smat_abs_e
C  Calculate overlap matrix for e ABS
C=======================================================================
      implicit none
C Input Variables
c     integer 
c     double precision
C Variables retuned
c     integer 
c     double precision
C Local Variables
      integer NauxBF
      integer i
      integer j
      integer ia
      integer IX,JX,KX
      integer LX,MX,NX
      double precision AX,AmatX(3)
      double precision BX,BmatX(3)
      double precision ans
      double precision smat
      double precision abs_norm_i
      double precision abs_norm_j
      integer  LXE(1000),MXE(1000),NXE(1000)
      double precision BXE(1000)
      double precision XmatE(3,1000)
      double precision zero
      parameter(zero=0.0d+00)


C  Open file to store overlap matrix for e ABS
      open(915,file='ABSSE.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
C  Open file with normalization factors for e ABS
      open(913,file='ABSNORME.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
c  Read RI auxilliary basis input file
      open(unit=911,file='aux_e_basis.inp',status='unknown')
      read(911,*) NauxBF
C Read in aux basis:
      do i=1,1000
         LXE(i)=0
         MXE(i)=0
         NXE(i)=0
         BXE(i)=zero
         XmatE(1,i)=zero
         XmatE(2,i)=zero
         XmatE(3,i)=zero
      end do

      do i=1,NauxBF
      read(911,*)LXE(i),MXE(i),NXE(i),BXE(i),
     x XmatE(1,i),XmatE(2,i),XmatE(3,i)
      end do


      do i=1,NauxBF

c        read(911,*) IX,JX,KX,AX,AmatX(1),AmatX(2),AmatX(3)
         read(913,REC=i) abs_norm_i
         IX=LXE(i)
         JX=MXE(i)
         KX=NXE(i)
         AX=BXE(i)
         AmatX(1)=XmatE(1,i)
         AmatX(2)=XmatE(2,i)
         AmatX(3)=XmatE(3,i)


         do j=1,NauxBF

c           read(911,*) LX,MX,NX,BX,BmatX(1),BmatX(2),BmatX(3)
            read(913,REC=j) abs_norm_j
            LX=LXE(j)
            MX=MXE(j)
            NX=NXE(j)
            BX=BXE(j)
            BmatX(1)=XmatE(1,j)
            BmatX(2)=XmatE(2,j)
            BmatX(3)=XmatE(3,j)


            call pack_2D(NauxBF,i,j,ia)

            call gfovlap(IX,JX,KX,AX,AmatX,
     2                   LX,MX,NX,BX,BmatX,ans)

            smat=abs_norm_i*abs_norm_j*ans
c           write(*,*)'abs_norm_i=',abs_norm_i
c           write(*,*)'abs_norm_j=',abs_norm_j
c           write(*,*)'smat=',smat

            write(915,REC=ia) smat

         end do
      end do

c           write(*,*)
      close(911)
      close(913)
      close(915)

      return
      end

C=======================================================================
      subroutine Smat_abs_p
C  Calculate overlap matrix for p ABS
C=======================================================================
      implicit none
C Input Variables
c     integer 
c     double precision
C Variables retuned
c     integer 
c     double precision
C Local Variables
      integer NauxBF
      integer i
      integer j
      integer ia
      integer IX,JX,KX
      integer LX,MX,NX
      double precision AX,AmatX(3)
      double precision BX,BmatX(3)
      double precision ans
      double precision smat
      double precision abs_norm_i
      double precision abs_norm_j
      integer  LXP(1000),MXP(1000),NXP(1000)
      double precision BXP(1000)
      double precision XmatP(3,1000)
      double precision zero
      parameter(zero=0.0d+00)

C  Open file to store overlap matrix for p ABS
      open(916,file='ABSSP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
C  Open file with normalization factors for e ABS
      open(914,file='ABSNORMP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
c  Read RI auxilliary basis input file
      open(unit=912,file='aux_p_basis.inp',status='unknown')
      read(912,*) NauxBF

C Read in aux basis:
      do i=1,1000
         LXP(i)=0
         MXP(i)=0
         NXP(i)=0
         BXP(i)=zero
         XmatP(1,i)=zero
         XmatP(2,i)=zero
         XmatP(3,i)=zero
      end do

      do i=1,NauxBF
      read(912,*)LXP(i),MXP(i),NXP(i),BXP(i),
     x XmatP(1,i),XmatP(2,i),XmatP(3,i)
      end do

      do i=1,NauxBF

c        read(412,*) IX,JX,KX,AX,AmatX(1),AmatX(2),AmatX(3)
         read(914,REC=i) abs_norm_i
         IX=LXP(i)
         JX=MXP(i)
         KX=NXP(i)
         AX=BXP(i)
         AmatX(1)=XmatP(1,i)
         AmatX(2)=XmatP(2,i)
         AmatX(3)=XmatP(3,i)

         do j=1,NauxBF

c           read(412,*) LX,MX,NX,BX,BmatX(1),BmatX(2),BmatX(3)
            read(914,REC=j) abs_norm_j
            LX=LXP(j)
            MX=MXP(j)
            NX=NXP(j)
            BX=BXP(j)
            BmatX(1)=XmatP(1,j)
            BmatX(2)=XmatP(2,j)
            BmatX(3)=XmatP(3,j)

            call pack_2D(NauxBF,i,j,ia)

            call gfovlap(IX,JX,KX,AX,AmatX,
     2                   LX,MX,NX,BX,BmatX,ans)

            smat=abs_norm_i*abs_norm_j*ans

            write(916,REC=ia) smat

         end do
      end do

      close(912)
      close(914)
      close(916)

      return
      end

C======================================================================
      subroutine ortabs_e(NauxBF)
C  Orthogonalize the e ABS
C======================================================================
      implicit none
C Input variables
      integer NauxBF
C Variables Returned

C Local Variables
      integer i
      integer j
      integer ia
      double precision S(NauxBF,NauxBF)
      double precision X(NauxBF,NauxBF)

c     write(*,*)'In ortabs_e'

C  Open overlap matrix file for e ABS
      open(915,file='ABSSE.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
C  Open file to store transformation matrix for e ABS
      open(917,file='ABSXE.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C  Read ABS overlap matrix into memory
      do i=1,NauxBF
         do j=1,NauxBF

            call pack_2D(NauxBF,i,j,ia)
            read(915,REC=ia) S(i,j)

         end do
      end do

      call ortabs(NauxBF,S,X)

C Store transformation matrix for use later
      do i=1,NauxBF
         do j=1,NauxBF

            call pack_2D(NauxBF,i,j,ia)
            write(917,REC=ia) X(i,j)
c           write(*,*)'i=',i,'j=',j,'X=',X(i,j)

         end do
      end do


      close(915)
      close(917)


      return
      end

C======================================================================
      subroutine ortabs_p(NauxBF)
C  Orthogonalize the p ABS
C======================================================================
      implicit none
C Input variables
      integer NauxBF
C Variables Returned

C Local Variables
      integer i
      integer j
      integer ia
      double precision S(NauxBF,NauxBF)
      double precision X(NauxBF,NauxBF)

C  Open overlap matrix file for p ABS
      open(916,file='ABSSP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)
C  Open file to store transformation matrix for e ABS
      open(918,file='ABSXP.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

C  Read ABS overlap matrix into memory
      do i=1,NauxBF
         do j=1,NauxBF

            call pack_2D(NauxBF,i,j,ia)
            read(916,REC=ia) S(i,j)

         end do
      end do

      call ortabs(NauxBF,S,X)

C Store transformation matrix for use later
      do i=1,NauxBF
         do j=1,NauxBF

            call pack_2D(NauxBF,i,j,ia)
            write(918,REC=ia) X(i,j)

         end do
      end do


      close(916)
      close(918)


      return
      end

C======================================================================
      subroutine ortabs(NB,S,X)
C CONSTRUCT ORTHONORMALIZING TRANSFORMATION MATRIX
C     ==============================
C     ON INPUT:
C     S         ::  OVERLAP MATRIX
C     NB        ::  NUMBER OF BASIS FUNCTIONS
C     ==============================
C     WORKING:
C     EV        ::  EIGENVALUES OF S MATRIX
C     EVECS     ::  EIGENVECTORS OF S MATRIX
C     EVECST    ::  TRANSPOSE OF EIGENVECTORS OF S MATRIX
C     X         ::  TRANSFORMATION MATRIX
C     XP        ::  TRANSPOSE OF TRANFORMATION MATRIX
C     FV1       ::  WORK SPACE
C     FV2       ::  WORK SPACE
C     FV3       ::  WORK SPACE
C     ==============================
C     OUTPUT:
C     X         ::  TRANSFORMATION MATRIX
C     ============================== 
C======================================================================
      implicit none
C Input variables
      integer NB
      double precision S(NB,NB)
C Variables Returned
      double precision X(NB,NB)
C Local Variables
      logical DEBUG
      integer i
      integer j
      integer IERR
      double precision ZERO
      parameter(ZERO=0.0d+00)
      double precision EVS(NB)
      double precision DEVS(NB,NB)
      double precision EVECS(NB,NB)
      double precision EVECST(NB,NB)
      double precision FV1(NB)
      double precision FV2(NB)
      double precision FV3(NB,NB)
C
C
C     ---> DIAGONALIZE OVERLAP MATRIX
C
c     DEBUG=.TRUE.
      DEBUG=.FALSE.
      CALL RS(NB,NB,S,EVS,2,EVECS,FV1,FV2,IERR)
      IF(DEBUG) THEN
         WRITE(*,*)
         WRITE(*,*)'---- OVERLAP MATRIX --> EIGENVALUES ----'
         WRITE(*,*)
         DO I=1,NB
               WRITE(*,*)'EVS(',I,')=',EVS(I)
         END DO
         WRITE(*,*)
         WRITE(*,*)'---- OVERLAP MATRIX --> EIGENVECTORS ----'
         WRITE(*,*)
         DO I=1,NB
            DO J=1,NB
               WRITE(*,*)'EVECS(',I,J,')=',EVECS(I,J)
            END DO
         END DO
      END IF
C           
C     ---- SYMMETRIC ORTHOGONALIZATION SCHEME ----
C
C     ---> FORM DIAGONAL EIGENVALUE MATRIX
C     ---> GET TRANSPOSE OF VECTORS FROM DIAGONALIZATION OF S
C
      DO I=1,NB
         DO J=1,NB
            DEVS(I,J)=ZERO
            EVECST(I,J)=EVECS(J,I)
         END DO
         DEVS(I,I)=1.0d+00/SQRT(EVS(I))
      END DO
      IF(DEBUG) THEN
         WRITE(*,*)
         WRITE(*,*)'---- S^-(1/2) EIGENVALUE MATRIX ----'
         WRITE(*,*)
         DO I=1,NB
            DO J=1,NB
               WRITE(*,*)'DEVS(',I,J,')=',DEVS(I,J)
            END DO
         END DO
      END IF
C
C     ---> FORM TRANSFORMATION MATRIX 
C
C     MULTIPLY EVECS.DEVS.EVECST = X
C     
      CALL MATMULT(NB,NB,NB,NB,DEVS,EVECST,FV3)
      CALL MATMULT(NB,NB,NB,NB,EVECS,FV3,X)
      IF(DEBUG) THEN
         WRITE(*,*)
         WRITE(*,*)'---- X MATRIX ----'
         WRITE(*,*)
         DO I=1,NB
            DO J=1,NB
               WRITE(*,*)'X(',I,J,')=',X(I,J)
            END DO
         END DO
      END IF
C
C     FORM TRANSPOSE OF TRANSFORMATION MATRIX
C
c     DO I=1,NB
c        DO J=1,NB
c           XP(I,J)=X(J,I)
c        END DO
c     END DO
C
C     TRANSFORM FOCK MATRIX 
C 
C     MULTIPLY XP.F.X = FP
C
c     CALL MATMULT(NB,NB,NB,NB,F,X,FV4)
c     CALL MATMULT(NB,NB,NB,NB,XP,FV4,FP)



      return
      end
