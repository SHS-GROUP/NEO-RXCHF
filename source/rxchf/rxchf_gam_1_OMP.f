C======================================================================
      subroutine RXCHF_GAM1_OMP_MD(nebf,npebf,npbf,ng1,ng1prm,nat,ngtg1,
     x                       pmass,cat,zan,bcoef1,gamma1,
     x                       AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)

C======================================================================
      implicit none
      include 'omp_lib.h'
c     include 'mpif.h'
C Input Variables
      integer nebf   ! Number of contracted elec basis functions
      integer npbf   ! Number of nuclear basis functions
      integer npebf  ! Number of primitive elec basis functions
      integer ng1    ! Number of contracted GAM_1 integrals
      integer ng1prm ! Number primitive gamma1 integrals
      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas
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
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)
C Local variables
      integer ia   ! Packing index
      integer ip   ! Nuc basis index
      integer jp   ! Nuc basis index
      integer ie1  ! Primitive elec basis func index 
      integer je1  ! Primitive elec basis func index 
      integer iec1 ! Contracted elec basis func index 
      integer jec1 ! Contracted elec basis func index 
      double precision GAM_1(ng1)
c     double precision XGAM_1(ng1)
      double precision GAM_1S(ng1)
c     double precision XGAM_1S(ng1)
      double precision Cof_ie1
      double precision Cof_je1
      double precision Cof_ip
      double precision Cof_jp
      double precision ans
      double precision ansE
      double precision ansS
      double precision zero
      parameter(zero=0.0d+00)
C--------------------------------(
C Basis set-related local variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3

      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision A3,Amat3(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
      double precision B3,Bmat3(3)
C--------------------------------)
C--------------------------------(
C MPI-Related Local variables
c     integer myid,ierr,nprocs
c     integer loopi,iLP
c     integer istart,iend
c     integer loop_map(ng1prm,4)
C--------------------------------)
C---OPENMP-RELATED-VARIABLES-----(
      integer id
      integer loopi,iLP
      double precision wtime
      integer loop_map(ng1prm,4)
C---OPENMP-RELATED-VARIABLES-----)

c     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
c     call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime()
C---OPENMP-TIMING------------------------------------------------------)

C-----------INITIALIZE-DATA-STRUCTURES-ON-MASTER-----------------------(
c     IF(myid .eq. 0) THEN
         write(*,*)
         write(*,*)'**************************************'
         write(*,*)'       GAM_1 Integrals    '
         write(*,*)
         write(*,*)'METHOD = MD'
         write(*,*)'npebf  =',npebf
         write(*,*)'nebf   =',nebf
         write(*,*)'npbf   =',npbf
         write(*,*)'ng1    =',ng1
         write(*,*)' Available processors: ',omp_get_num_procs()
         write(*,*)' Available threads     ',omp_get_max_threads()
         write(*,*)' Threads in use        ',omp_get_num_threads()
         write(*,*)'**************************************'
         write(*,*)


C  zero out the 1-D arrays to hold the contracted integrals
         do ia=1,ng1
            GAM_1(ia)=zero
            GAM_1s(ia)=zero
c           XGAM_1(ia)=zero
c           XGAM_1s(ia)=zero
         end do
C Compress nested loops
         Loopi=0
         do ip=1,npbf
            do jp=1,npbf
               do ie1=1,npebf
                  do je1=1,npebf

                     Loopi=Loopi+1
                     loop_map(Loopi,1)=je1
                     loop_map(Loopi,2)=ie1
                     loop_map(Loopi,3)=jp
                     loop_map(Loopi,4)=ip

                  end do
               end do
            end do
         end do
c     END IF  ! end if for myid
C-----------INITIALIZE-DATA-STRUCTURES-ON-MASTER-----------------------)
c     write(*,*)'1 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(GAM_2,ng2,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(GAM_1,ng1,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(XGAM_1,ng1,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(GAM_2s,ng2,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(GAM_1s,ng1,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(XGAM_1s,ng1,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     write(*,*)'1 AFTER MPI BARRIER MYID=',myid
c     call MPI_BCAST(loop_map,ng1prm*4,mpi_integer,0,
c    x     MPI_COMM_WORLD,ierr)
c     write(*,*)'BCAST LoopMap IERR=',ierr
c     write(*,*)'2 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'2 AFTER MPI BARRIER MYID=',myid

c     write(*,*)'16 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'16 AFTER MPI BARRIER MYID=',myid
c     call loop_size(1,ng1prm,nprocs,myid,istart,iend)
c     write(*,*)'17 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'17 AFTER MPI BARRIER MYID=',myid
c     write(*,*)'MYID  =',myid
c     write(*,*)'ISTART=',istart
c     write(*,*)'IEND  =',iend

c     if(myid.eq.1) then
c        write(*,*)'PROC1: AGEBFCC is:'
c        write(*,*) AGEBFCC
c        write(*,*)'PROC1: AMPEB2C is:'
c        write(*,*) AMPEB2C
c     end if
c     if(myid.eq.0) then
c        write(*,*)'MYID=',myid
c        write(*,*)'zan is:',zan
c        write(*,*)'cat is:',cat
c        write(*,*)'pmass=',pmass
c        write(*,*)
c        write(*,*)'CHECK AMPEB2C, MYID=',myid
c        write(*,*) AMPEB2C
c     end if
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     if(myid.eq.1) then
c        write(*,*)'MYID=',myid
c        write(*,*)'zan is:',zan
c        write(*,*)'cat is:',cat
c        write(*,*)'pmass=',pmass
c        write(*,*)
c        write(*,*)'CHECK AMPEB2C, MYID=',myid
c        write(*,*) AMPEB2C
c     end if
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf)
!$ompx shared(ng1)
!$ompx shared(ng1prm)
!$ompx shared(GAM_1,GAM_1s)
!$ompx private(iLp) 
!$ompx private(ia)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(iec1,jec1)
!$ompx private(A1,I1,J1,K1,Amat1)
!$ompx private(A2,I2,J2,K2,Amat2)
!$ompx private(B1,L1,M1,N1,Bmat1)
!$ompx private(B2,L2,M2,N2,Bmat2)
!$ompx private(Cof_ie1,Cof_je1)
!$ompx private(Cof_ip,Cof_jp)
!$ompx private(ansE,ansS)
!$ompx private(id)
CCCCCCCCCCCCCCCCCC!$ompx reduction(+:GAM_1,GAM_1s)

      id= omp_get_thread_num()
      write(*,*)' Hello from process ',id
      if(id.eq.0) then
         write(*,*)'Threads in use', omp_get_num_threads()
      end if

!$omp do
      do iLP=1,ng1prm
c     do ip=1,npbf
c        do jp=1,npbf
c           do ie1=1,npebf
c              do je1=1,npebf
C  Map loop indices
         je1=loop_map(iLP,1)
         ie1=loop_map(iLP,2)
         jp =loop_map(iLP,3)
         ip =loop_map(iLP,4)

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
c                    call xcalc_GAM1_MD(ip,jp,ie1,je1,ansE,ansS)
c                    ansS=0.0d+00
c                    call xcalc_GAM1s_MD(ip,jp,ie1,je1,ansS)
         call RXCHF_xcalc_GAM1_MD(I1,J1,K1,A1,Amat1,
     x                      I2,J2,K2,A2,Amat2,
     x                      L1,M1,N1,B1,Bmat1,
     x                      L2,M2,N2,B2,Bmat2,
     x                      nat,ngtg1,
     x                      pmass,cat,zan,
     x                      bcoef1,gamma1,
     x                      ansE,ansS)


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
     x                jec1,iec1,jp,ip,ia)

         GAM_1(ia)=GAM_1(ia)+ansE*
     x             Cof_ip*Cof_jp*Cof_ie1*Cof_je1

         GAM_1s(ia)=GAM_1s(ia)+ansS*
     x              Cof_ip*Cof_jp*Cof_ie1*Cof_je1


c              end do
c           end do
c        end do
c     end do
      end do

!$omp end do
!$omp end parallel
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

C---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime() - wtime
      write(*,*)'TIME TO CALCULATE GAM_1 INTEGRALS: ',wtime
C---OPENMP-TIMING------------------------------------------------------)

c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c     call MPI_REDUCE(GAM_1,XGAM_1,ng1,
c    x                mpi_double_precision,mpi_sum,0,
c    x                MPI_COMM_WORLD,ierr)

c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c     call MPI_REDUCE(GAM_1s,XGAM_1s,ng1,
c    x                mpi_double_precision,mpi_sum,0,
c    x                MPI_COMM_WORLD,ierr)

c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)



c     IF(myid .eq. 0) THEN
C  Open files to store integrals
         open(801,file='GAM_1.ufm',form='unformatted',
     x    status='unknown',access='direct',RECL=8)

         open(802,file='GAM_1s.ufm',form='unformatted',
     x    status='unknown',access='direct',RECL=8)

c        GAM_1=XGAM_1
c        GAM_1s=XGAM_1s

         do ip=1,npbf
            do jp=1,npbf
               do iec1=1,nebf
                  do jec1=1,nebf

                     call pack_4D(nebf,nebf,npbf,
     x                            jec1,iec1,jp,ip,ia)

C  Write to file
                     write(801,REC=ia) GAM_1(ia)
                     write(802,REC=ia) GAM_1s(ia)


                  end do
               end do
            end do
         end do

         close(801)
         close(802)

C ENDIF OF MYID=0: 
c     END IF

c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)


      return
      end

C======================================================================
      subroutine RXCHF_xcalc_GAM1_MD(I1,J1,K1,A1,Amat1,
     x                         I2,J2,K2,A2,Amat2,
     x                         L1,M1,N1,B1,Bmat1,
     x                         L2,M2,N2,B2,Bmat2,
     x                         nat,ngtg1,
     x                         pmass,cat,zan,
     x                         bcoef1,gamma1,
     x                         ansE,ansS)

C======================================================================
      implicit none
C Input Variables
      integer nat,ngtg1
C--------------------------------(
C Basis set-related variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3

      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision A3,Amat3(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
      double precision B3,Bmat3(3)
C--------------------------------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)
C Variables Returned
      double precision ansE
      double precision ansS
C Local Variables
      integer iii
      integer K,L
      logical debug

      double precision ans
      double precision ZNUC
      double precision Cmat(3)
      double precision OXX
      double precision OGX
      double precision OGG
      double precision TEXG  
      double precision TPXG  
      double precision TEGX  
      double precision TPGX  
      double precision TEGG  
      double precision TPGG  
      double precision VCEGX 
      double precision VCEGG 
      double precision VCPGX 
      double precision VCPGG 
      double precision VMGX  
      double precision VMGG  
      double precision GV 
      double precision BCFV 
      double precision GW 
      double precision BCFW 
      double precision ogxix,ogxiy,ogxiz 
      double precision tegxi,tpgxi 
      double precision VCEGXI 
      double precision VCPGXI 
      double precision VCEGXIA 
      double precision VCPGXIA 
      double precision VMGXI 
      double precision texgi,tpxgi 
      double precision oggix,oggiy,oggiz 
      double precision teggi,tpggi
      double precision VCEGGI
      double precision VCPGGI
      double precision VCEGGIA
      double precision VCPGGIA
      double precision VMGGI
      double precision dum1,dum2,dum3,sval

      double precision S   
      double precision TE  
      double precision TP  
      double precision VCE 
      double precision VCP 
      double precision VMX 

      double precision ZERO,one 
      parameter(ZERO=0.0d+00,one=1.0d+00)


C     debug=.true.
      debug=.false.

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

c     call get_BF(1,ie1,I1,J1,K1,A1,Amat1)
c     call get_BF(1,je1,L1,M1,N1,B1,Bmat1)
c     call get_BF(2,ip,I2,J2,K2,A2,Amat2)
c     call get_BF(2,jp,L2,M2,N2,B2,Bmat2)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C
      OGX=ZERO
      OGG=ZERO
      TEXG=ZERO
      TPXG=ZERO
      TEGX=ZERO
      TPGX=ZERO
      TEGG=ZERO
      TPGG=ZERO
      VCEGX=ZERO
      VCEGG=ZERO
      VCPGX=ZERO
      VCPGG=ZERO
      VMGX=ZERO
      VMGG=ZERO
C
      DO K=1,NGTG1
           GV=GAMMA1(K)
         BCFV=BCOEF1(K)
C
c     SUBROUTINE pgike1(I1,J1,K1,alp1,Amat1,
c    1                  I2,J2,K2,alp2,Amat2,
c    2                  L1,M1,N1,beta1,Bmat1,
c    3                  L2,M2,N2,beta2,Bmat2,
c    4                  gamA,gamB,xmass1,ans,sval)

C         CALL PGIKE1(I1,J1,K1,A1,Amat1,
C     1               I2,J2,K2,A2,Amat2,
C     2               L1,M1,N1,B1,Bmat1,
C     3               L2,M2,N2,B2,Bmat2,
C     4               gv,zero,one,tegxi,sval)
c    5               ogxix,ogxiy,ogxiz)
c        write(*,*)'electronic kinetic gx integral = ',tegxi
C
C         CALL PGIKE1(I2,J2,K2,A2,Amat2,
C     1               I1,J1,K1,A1,Amat1,
C     2               L2,M2,N2,B2,Bmat2,
C     3               L1,M1,N1,B1,Bmat1,
C     4               gv,zero,pmass,tpgxi,sval)
c    5               ogxix,ogxiy,ogxiz)
C
c        OGX=OGX+BCFV*OGXIX*OGXIY*OGXIZ
C         tegx=tegx+bcfv*tegxi
C         tpgx=tpgx+bcfv*tpgxi
C
C         VCEGXI = ZERO
C         VCPGXI = ZERO
C         DO III=1,NAT
C            cmat(1)=cat(1,III)
C            cmat(2)=cat(2,III)
C            cmat(3)=cat(3,III)
C            ZNUC=ZAN(III)
C
c           CALL PGIVEC(I1,J1,K1,A1,Amat1,
c    1                  I2,J2,K2,A2,Amat2,
c    2                  L1,M1,N1,B1,Bmat1,
c    3                  L2,M2,N2,B2,Bmat2,
c    4                  GV,zero,CMAT,VCEGXIA)
c           VCEGXI=VCEGXIA*ZNUC+VCEGXI
c           write(*,*)'AC:  VCEGXIA=',VCEGXIA
C            call G2_MD_xggvec(I1,J1,K1,A1,Amat1,
C     1                        I2,J2,K2,A2,Amat2,
C     2                        L1,M1,N1,B1,Bmat1,
C     3                        L2,M2,N2,B2,Bmat2,
C     4                        GV,zero,Cmat,ZNUC,
C     5                        VCEGXIA)
C            VCEGXI=VCEGXIA*ZNUC+VCEGXI
c           write(*,*)'MD:  VCEGXIA=',VCEGXIA
c           write(*,*)

C
c           CALL PGIVEC(I2,J2,K2,A2,Amat2,
c    1                  I1,J1,K1,A1,Amat1,
c    2                  L2,M2,N2,B2,Bmat2,
c    3                  L1,M1,N1,B1,Bmat1,
c    4                  GV,zero,CMAT,VCPGXIA)
c           VCPGXI=VCPGXIA*ZNUC+VCPGXI
C            CALL G2_MD_xggvec(I2,J2,K2,A2,Amat2,
C     x                        I1,J1,K1,A1,Amat1,
C     x                        L2,M2,N2,B2,Bmat2,
C     x                        L1,M1,N1,B1,Bmat1,
C     x                        GV,zero,CMAT,ZNUC,
C     x                        VCPGXIA)
C            VCPGXI=VCPGXIA*ZNUC+VCPGXI
C
C         END DO
C
c        CALL PGIVEE(I1,J1,K1,A1,Amat1,
c    1               I2,J2,K2,A2,Amat2,
c    2               L1,M1,N1,B1,Bmat1,
c    3               L2,M2,N2,B2,Bmat2,
c    4               GV,ZERO,VMGXI)
C         call G2_MD_xggvee(I1,J1,K1,A1,Amat1,
C     1                     I2,J2,K2,A2,Amat2,
C     2                     L1,M1,N1,B1,Bmat1,
C     3                     L2,M2,N2,B2,Bmat2,
C     4                     GV,zero,VMGXI)

C
C         VCEGX=VCEGX+BCFV*VCEGXI
C         VCPGX=VCPGX+BCFV*VCPGXI
C         VMGX=VMGX+BCFV*VMGXI
C
C         CALL PGIKE1(I1,J1,K1,A1,Amat1,
C     1               I2,J2,K2,A2,Amat2,
C     2               L1,M1,N1,B1,Bmat1,
C     3               L2,M2,N2,B2,Bmat2,
C     4               zero,gv,one,texgi,sval)
c    5               dum1,dum2,dum3)
c        write(*,*)'electronic kinetic xg integral = ',texgi
C
C         CALL PGIKE1(I2,J2,K2,A2,Amat2,
C     1               I1,J1,K1,A1,Amat1,
C     2               L2,M2,N2,B2,Bmat2,
C     3               L1,M1,N1,B1,Bmat1,
C     4               zero,gv,pmass,tpxgi,sval)
c    5               dum1,dum2,dum3)
C
C         texg=texg+bcfv*texgi
C         tpxg=tpxg+bcfv*tpxgi
C
         DO L=1,NGTG1
            GW=GAMMA1(L)
            BCFW=BCOEF1(L)
C
            CALL PGIKE1(I1,J1,K1,A1,Amat1,
     1                  I2,J2,K2,A2,Amat2,
     2                  L1,M1,N1,B1,Bmat1,
     3                  L2,M2,N2,B2,Bmat2,
     4                  gv,gw,one,teggi,sval)
c    5                  oggix,oggiy,oggiz)
c        write(*,*)'electronic kinetic gg integral = ',teggi
c        write(*,*)
C
            CALL PGIKE1(I2,J2,K2,A2,Amat2,
     1                  I1,J1,K1,A1,Amat1,
     2                  L2,M2,N2,B2,Bmat2,
     3                  L1,M1,N1,B1,Bmat1,
     4                  gv,gw,pmass,tpggi,sval)
c    5                  oggix,oggiy,oggiz)
C
c           OGG=OGG+BCFV*BCFW*OGGIX*OGGIY*OGGIZ
            tegg=tegg+bcfv*bcfw*teggi
            tpgg=tpgg+bcfv*bcfw*tpggi
C
            VCEGGI = ZERO
            VCPGGI = ZERO
            DO III=1,NAT
               Cmat(1)=cat(1,III)
               Cmat(2)=cat(2,III)
               Cmat(3)=cat(3,III)
               ZNUC=ZAN(III)
C
c              CALL PGIVEC(I1,J1,K1,A1,Amat1,
c    1                     I2,J2,K2,A2,Amat2,
c    2                     L1,M1,N1,B1,Bmat1,
c    3                     L2,M2,N2,B2,Bmat2,
c    4                     GV,GW,CMAT,VCEGGIA)
c              VCEGGI=VCEGGIA*ZNUC+VCEGGI
               call G2_MD_xggvec(I1,J1,K1,A1,Amat1,
     1                           I2,J2,K2,A2,Amat2,
     2                           L1,M1,N1,B1,Bmat1,
     3                           L2,M2,N2,B2,Bmat2,
     4                           GV,GW,Cmat,ZNUC,
     5                           VCEGGIA)
               VCEGGI=VCEGGIA*ZNUC+VCEGGI
C
c              CALL PGIVEC(I2,J2,K2,A2,Amat2,
c    1                     I1,J1,K1,A1,Amat1,
c    2                     L2,M2,N2,B2,Bmat2,
c    3                     L1,M1,N1,B1,Bmat1,
c    4                     GV,GW,CMAT,VCPGGIA)
c              VCPGGI=VCPGGIA*ZNUC+VCPGGI
               CALL G2_MD_xggvec(I2,J2,K2,A2,Amat2,
     x                           I1,J1,K1,A1,Amat1,
     x                           L2,M2,N2,B2,Bmat2,
     x                           L1,M1,N1,B1,Bmat1,
     x                           GV,GW,CMAT,ZNUC,
     x                           VCPGGIA)
               VCPGGI=VCPGGIA*ZNUC+VCPGGI
C
            END DO
C
c           CALL PGIVEE(I1,J1,K1,A1,Amat1,
c    1                  I2,J2,K2,A2,Amat2,
c    2                  L1,M1,N1,B1,Bmat1,
c    3                  L2,M2,N2,B2,Bmat2,
c    4                  GV,GW,VMGGI)
            call G2_MD_xggvee(I1,J1,K1,A1,Amat1,
     1                        I2,J2,K2,A2,Amat2,
     2                        L1,M1,N1,B1,Bmat1,
     3                        L2,M2,N2,B2,Bmat2,
     4                        GV,GW,VMGGI)
C
            VCEGG=VCEGG+BCFV*BCFW*VCEGGI
            VCPGG=VCPGG+BCFV*BCFW*VCPGGI
            VMGG=VMGG+BCFV*BCFW*VMGGI
         END DO
      END DO
C
c     S=OXX+2*OGX+OGG
c     TE=TEXX+TEGX+TEXG+TEGG
c     TP=TPXX+TPGX+TPXG+TPGG
c     VCE=VCEXX+2*VCEGX+VCEGG
c     VCP=VCPXX+2*VCPGX+VCPGG
c     VMX=VMXX+2*VMGX+VMGG

c     S=2.0d+00*OGX+OGG
      TE=TEGX+TEXG+TEGG
      TP=TPGX+TPXG+TPGG
      VCE=2.0d+00*VCEGX+VCEGG
      VCP=2.0d+00*VCPGX+VCPGG
      VMX=2.0d+00*VMGX+VMGG

      ansE=TE+TP-VCE+VCP-VMX
c     ansS=S
C
c     HE =TE-VCE
c     HP =TP+VCP
c     HEP=VMX
c     IF(DEBUG) THEN
c        WRITE(*,*)
c        WRITE(*,*)'INDICES => IJKL =',II,JJ,KK,LL
c        WRITE(*,*)'TE  = ',TE
c        WRITE(*,*)'TP  = ',TP
c        WRITE(*,*)'VCE = ',VCE
c        WRITE(*,*)'VCP = ',VCP
c        WRITE(*,*)'VMX = ',VMX
c        WRITE(*,*)'HE  = ',HE
c        WRITE(*,*)'HP  = ',HP
c        WRITE(*,*)'HEP = ',HEP
c        WRITE(*,*)
c     END IF
C


      OXX=ZERO
      OGX=ZERO
      OGG=ZERO
C   
c     call pgiovlap_1D(I1,A1,amat1(1),
c    &                 I2,A2,amat2(1),
c    &                 L1,B1,bmat1(1),
c    &                 L2,B2,bmat2(1),
c    &                 zero,zero,OXXX)
C
c     call pgiovlap_1D(J1,A1,amat1(2),
c    &                 J2,A2,amat2(2),
c    &                 M1,B1,bmat1(2),
c    &                 M2,B2,bmat2(2),
c    &                 zero,zero,OXXY)
C
c     call pgiovlap_1D(K1,A1,amat1(3),
c    &                 K2,A2,amat2(3),
c    &                 N1,B1,bmat1(3),
c    &                 N2,B2,bmat2(3),
c    &                 zero,zero,OXXZ)
C
c    oxx=oxxx+oxxy+oxxz

C
      DO K=1,NGTG1
           GV=GAMMA1(K)
         BCFV=BCOEF1(K)
C
c        call pgiovlap_1D(I1,A1,amat1(1),
c    &                    I2,A2,amat2(1),
c    &                    L1,B1,bmat1(1),
c    &                    L2,B2,bmat2(1),
c    &                    GV,zero,OGXIX)

c        call pgiovlap_1D(J1,A1,amat1(2),
c    &                    J2,A2,amat2(2),
c    &                    M1,B1,bmat1(2),
c    &                    M2,B2,bmat2(2),
c    &                    GV,zero,OGXIY)

c        call pgiovlap_1D(K1,A1,amat1(3),
c    &                    K2,A2,amat2(3),
c    &                    N1,B1,bmat1(3),
c    &                    N2,B2,bmat2(3),
c    &                    GV,zero,OGXIZ)
C
c        OGX=OGX+BCFV*OGXIX*OGXIY*OGXIZ
C
C         call G2_MD_xggs(I1,J1,K1,A1,Amat1,
C     1                   I2,J2,K2,A2,Amat2,
C     2                   L1,M1,N1,B1,Bmat1,
C     3                   L2,M2,N2,B2,Bmat2,
C     4                   GV,zero,OGXIX)


C         OGX=OGX+BCFV*OGXIX
C
         DO L=1,NGTG1
            GW=GAMMA1(L)
            BCFW=BCOEF1(L)
C
c           call pgiovlap_1D(I1,A1,amat1(1),
c    &                       I2,A2,amat2(1),
c    &                       L1,B1,bmat1(1),
c    &                       L2,B2,bmat2(1),
c    &                       GV,GW,OGGIX)
C
c           call pgiovlap_1D(J1,A1,amat1(2),
c    &                       J2,A2,amat2(2),
c    &                       M1,B1,bmat1(2),
c    &                       M2,B2,bmat2(2),
c    &                       GV,GW,OGGIY)
C
c           call pgiovlap_1D(K1,A1,amat1(3),
c    &                       K2,A2,amat2(3),
c    &                       N1,B1,bmat1(3),
c    &                       N2,B2,bmat2(3),
c    &                       GV,GW,OGGIZ)
C
c           OGG=OGG+BCFV*BCFW*OGGIX*OGGIY*OGGIZ

            call G2_MD_xggs(I1,J1,K1,A1,Amat1,
     1                      I2,J2,K2,A2,Amat2,
     2                      L1,M1,N1,B1,Bmat1,
     3                      L2,M2,N2,B2,Bmat2,
     4                      GV,GW,OGGIX)
C
            OGG=OGG+BCFV*BCFW*OGGIX


         end do
      end do

      S=OXX+2.0d+00*OGX+OGG
      ansS=s

C ARS( testing
C      if (ansE.gt.1.0d0) then
C       write(*,*)
C       write(*,*) "large ansE:",ansE
C       write(*,*) "TE,TP,VCE,VCP,VMX:"
C       write(*,*) TE,TP,VCE,VCP,VMX
C       write(*,*) "I1,J1,K1,I2,J2,K2,L1,M1,N1,L2,M2,N2:"
C       write(*,*) I1,J1,K1,I2,J2,K2,L1,M1,N1,L2,M2,N2
C       write(*,*) "A1,A2,B1,B2,GV,GW:"
C       write(*,*) A1,A2,B1,B2,GV,GW
C       write(*,*)
C      end if
C      if (ansS.gt.1.0d0) then
C       write(*,*)
C       write(*,*) "large ansS:",ansS
C       write(*,*) "OXX,OGX,OGG:"
C       write(*,*) OXX,OGX,OGG
C       write(*,*) "I1,J1,K1,I2,J2,K2,L1,M1,N1,L2,M2,N2:"
C       write(*,*) I1,J1,K1,I2,J2,K2,L1,M1,N1,L2,M2,N2
C       write(*,*) "A1,A2,B1,B2,GV,GW:"
C       write(*,*) A1,A2,B1,B2,GV,GW
C       write(*,*)
C      end if
C )
      return
      end



