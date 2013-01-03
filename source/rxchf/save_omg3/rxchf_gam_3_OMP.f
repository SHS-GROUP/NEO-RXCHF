C======================================================================
      subroutine RXCHF_GAM3_OMP_MD(nebf,npebf,npbf,ng3,ng3prm,nat,ngtg1,
     x                       pmass,cat,zan,bcoef1,gamma1,
     x                       AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
c    x                       XGAM_2,XGAM_2s)

C======================================================================
      implicit none
      include 'omp_lib.h'
c     include 'mpif.h'

C Input Variables
      integer ng3    ! Number contracted gamma3 integrals
      integer ng3prm ! Number primitive gamma3 integrals
      integer npebf  ! Number primitive electronic basis functions
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions

      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas

C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
C-------Basis Set Info-------)

      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1) 
      double precision gamma1(ngtg1)

      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
C Local variables
      integer ia
      integer ia_123
      integer ia_132
      integer ia_213
      integer ia_231
      integer ia_312
      integer ia_321
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3
      integer iec1,jec1,iec2,jec2,iec3,jec3
      double precision Cof_ie1,Cof_je1
      double precision Cof_ie2,Cof_je2
      double precision Cof_ie3,Cof_je3
      double precision Cof_ip,Cof_jp

      double precision zero,half,six
      parameter(zero=0.0d+00,half=0.5d+00,six=6.0d+00)
      double precision ans

c     double precision ansE
c     double precision ansS

      double precision x123
      double precision x132
      double precision x213
      double precision x231
      double precision x312
      double precision x321
      double precision xxxx

c     double precision val_gam2s
C-------------------------------------
      double precision GAM_3(ng3)
c     double precision XGAM_3(ng3)
C--------------------------------(
C Basis set-related local variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1,Amat1(3) 
      double precision A2,Amat2(3) 
      double precision A3,Amat3(3) 
      double precision A4,Amat4(3) 
      double precision B1,Bmat1(3) 
      double precision B2,Bmat2(3) 
      double precision B3,Bmat3(3) 
      double precision B4,Bmat4(3) 
C--------------------------------)
C--------------------------------(
C MPI-Related Local variables
c     integer myid,ierr,nprocs
c     integer loopi,iLP
c     integer istart,iend
c     integer loop_map(ng3prm,8)
C--------------------------------)
C---OPENMP-RELATED-VARIABLES-----(
      integer id
      integer loopi,iLP
      double precision wtime
      integer loop_map(ng3prm,8)
C---OPENMP-RELATED-VARIABLES-----)
CCWS_int_stat
c     integer iii
c     double precision tol
CCWS_int_stat

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
         write(*,*)'         GAM_3 Integrals    '
         write(*,*)
         write(*,*)'METHOD = MD'
         write(*,*)'npebf=',npebf
         write(*,*)'nebf=',nebf
         write(*,*)'npbf=',npbf
         write(*,*)'nprim_g2 =',ng3prm
         write(*,*)'ng2 =',ng3
         write(*,*)' Available processors: ',omp_get_num_procs()
         write(*,*)' Available threads     ',omp_get_max_threads()
         write(*,*)' Threads in use        ',omp_get_num_threads()
c        write(*,*)'NPROC =',nprocs
         write(*,*)'**************************************'
         write(*,*)
c     end if
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     IF(myid .eq. 0) THEN
C Initialize 1-D arrays to hold the contracted integrals
         do ia=1,ng3
            GAM_3(ia)=zero
c           XGAM_3(ia)=zero
         end do
c     end if

c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c     IF(myid .eq. 0) THEN
c     write(*,*)'GOT HERE 1 MYID=',myid
C Compress nested loops to 1-dimension
         Loopi=0
         do ip=1,npbf
         do jp=1,npbf
            do ie1=1,npebf
            do je1=1,npebf
               do ie2=1,npebf
               do je2=1,npebf
                  do ie3=1,npebf
                  do je3=1,npebf

                           Loopi=Loopi+1
                           loop_map(Loopi,1)=je3
                           loop_map(Loopi,2)=ie3
                           loop_map(Loopi,3)=je2
                           loop_map(Loopi,4)=ie2
                           loop_map(Loopi,5)=je1
                           loop_map(Loopi,6)=ie1
                           loop_map(Loopi,7)=jp
                           loop_map(Loopi,8)=ip

                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do
c     write(*,*)'GOT HERE 2 MYID=',myid
C ENDIF OF MYID=0: 
c     else
c       loop_map=0.0d+00
c     END IF
C-----------INITIALIZE-DATA-STRUCTURES-ON-MASTER-----------------------)

c     write(*,*)'1 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     call MPI_BCAST(GAM_3,ng3,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)
c     call MPI_BCAST(XGAM_3,ng3,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)
c     write(*,*)'1 AFTER MPI BARRIER MYID=',myid
c     call MPI_BCAST(loop_map,ng3prm*8,mpi_integer,0,
c    x     MPI_COMM_WORLD,ierr)
c     write(*,*)'BCAST LoopMap IERR=',ierr
c     write(*,*)'2 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'2 AFTER MPI BARRIER MYID=',myid


c     write(*,*)'16 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'16 AFTER MPI BARRIER MYID=',myid
c     call loop_size(1,ng3prm,nprocs,myid,istart,iend)
c     write(*,*)'17 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'17 AFTER MPI BARRIER MYID=',myid
c     if(myid.eq.0) then
c        write(*,*)'MYID  =',myid
c        write(*,*)'ISTART=',istart
c        write(*,*)'IEND  =',iend
c        write(*,*)'ngtg11 is:',ngtg1
c        write(*,*)'bcoef1 is:',bcoef1
c        write(*,*)'gamma1 is:',gamma1
c        write(*,*)'ELCEX is:',ELCEX
c        write(*,*)'NUCEX is:',NUCEX
c        write(*,*)'ELCAM is:',ELCAM
c        write(*,*)'AMPEB2C is:',AMPEB2C
c        write(*,*)'AGEBFCC is:',AGEBFCC
c        write(*,*)'AGNBFCC is:',AGNBFCC
c        write(*,*)'ELCBFC is:',ELCBFC
c     end if
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     if(myid.eq.1) then
c        write(*,*)'MYID  =',myid
c        write(*,*)'ISTART=',istart
c        write(*,*)'IEND  =',iend
c        write(*,*)'ngtg11 is:',ngtg1
c        write(*,*)'bcoef1 is:',bcoef1
c        write(*,*)'gamma1 is:',gamma1
c        write(*,*)'ELCEX is:',ELCEX
c        write(*,*)'NUCEX is:',NUCEX
c        write(*,*)'ELCAM is:',ELCAM
c        write(*,*)'AMPEB2C is:',AMPEB2C
c        write(*,*)'AGEBFCC is:',AGEBFCC
c        write(*,*)'AGNBFCC is:',AGNBFCC
c        write(*,*)'ELCBFC is:',ELCBFC
c     end if
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c     if(myid.eq.1) then
c        write(*,*)'PROC1: AGEBFCC is:'
c        write(*,*) AGEBFCC
c        write(*,*)'PROC1: AMPEB2C is:'
c        write(*,*) AMPEB2C
c     end if
c     if(myid.eq.0) then 
c        write(*,*)'CHECK AMPEB2C, MYID=',myid
c        write(*,*) AMPEB2C
c     end if
c     if(myid.eq.1) then 
c        write(*,*)'CHECK AMPEB2C, MYID=',myid
c        write(*,*) AMPEB2C
c     end if
c     call mpi_final(ierr)
c     stop

C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf)
!$ompx shared(ng3)
!$ompx shared(ng3prm)
!$ompx shared(GAM_3)
!$ompx private(iLp) 
!$ompx private(ia)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(ie3,je3) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(iec3,jec3)
!$ompx private(A1,I1,J1,K1,Amat1)
!$ompx private(A2,I2,J2,K2,Amat2)
!$ompx private(A3,I3,J3,K3,Amat3)
!$ompx private(A4,I4,J4,K4,Amat4)
!$ompx private(B1,L1,M1,N1,Bmat1)
!$ompx private(B2,L2,M2,N2,Bmat2)
!$ompx private(B3,L3,M3,N3,Bmat3)
!$ompx private(B4,L4,M4,N4,Bmat4)
!$ompx private(Cof_ie1,Cof_je1)
!$ompx private(Cof_ie2,Cof_je2)
!$ompx private(Cof_ie3,Cof_je3)
!$ompx private(Cof_ip,Cof_jp)
!$ompx private(ans)
!$ompx private(id)
CCCCCCCCCCCCCCCCCmpx reduction(+:GAM_3)

      id= omp_get_thread_num()
      write(*,*)' Hello from process ',id
      if(id.eq.0) then
         write(*,*)'Threads in use', omp_get_num_threads()
      end if

!$omp do
c     istart=1
c     iend=2
c     do iLP=istart,iend
      do iLP=1,ng3prm
c        do ip=1,npbf
c        do jp=1,npbf
c           do ie1=1,npebf
c           do je1=1,npebf
c              do ie2=1,npebf
c              do je2=1,npebf
c                 do ie3=1,npebf
c                 do je3=1,npebf

c        write(*,*)'DBG MYID=',myid
C  Map loop indices
         je3=loop_map(iLP,1)
         ie3=loop_map(iLP,2)
         je2=loop_map(iLP,3)
         ie2=loop_map(iLP,4)
         je1=loop_map(iLP,5)
         ie1=loop_map(iLP,6)
         jp =loop_map(iLP,7)
         ip =loop_map(iLP,8)
c                 write(*,*)'**DBG iLP= MYID=',iLP,myid
c                 write(*,*)'**DBG IP = JP =  MYID=',ip,jp,myid
c                 write(*,*)'**DBG IE1= JE1=  MYID=',ie1,je1,myid
c                 write(*,*)'**DBG IE2= JE2=  MYID=',ie2,je2,myid
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

         A3=ELCEX(ie3)
         I3=ELCAM(ie3,1)
         J3=ELCAM(ie3,2)
         K3=ELCAM(ie3,3)
         Amat3(1)=ELCBFC(ie3,1)
         Amat3(2)=ELCBFC(ie3,2)
         Amat3(3)=ELCBFC(ie3,3)

         A4=NUCEX(ip)
         I4=NUCAM(ip,1)
         J4=NUCAM(ip,2)
         K4=NUCAM(ip,3)
         Amat4(1)=NUCBFC(ip,1)
         Amat4(2)=NUCBFC(ip,2)
         Amat4(3)=NUCBFC(ip,3)

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

         B3=ELCEX(je3)
         L3=ELCAM(je3,1)
         M3=ELCAM(je3,2)
         N3=ELCAM(je3,3)
         Bmat3(1)=ELCBFC(je3,1)
         Bmat3(2)=ELCBFC(je3,2)
         Bmat3(3)=ELCBFC(je3,3)

         B4=NUCEX(jp)
         L4=NUCAM(jp,1)
         M4=NUCAM(jp,2)
         N4=NUCAM(jp,3)
         Bmat4(1)=NUCBFC(jp,1)
         Bmat4(2)=NUCBFC(jp,2)
         Bmat4(3)=NUCBFC(jp,3)

         call RXCHF_xcalc_GAM3_MD(I1,J1,K1,A1,Amat1,
     x                      I2,J2,K2,A2,Amat2,
     x                      I3,J3,K3,A3,Amat3,
     x                      I4,J4,K4,A4,Amat4,
     x                      L1,M1,N1,B1,Bmat1,
     x                      L2,M2,N2,B2,Bmat2,
     x                      L3,M3,N3,B3,Bmat3,
     x                      L4,M4,N4,B4,Bmat4,
     x                      nat,ngtg1,
     x                      pmass,cat,zan,
     x                      bcoef1,gamma1,
     x                      ans)

c                       call underflow(ans)

C  Map from primitive BF indices to contracted indices
                        iec1=AMPEB2C(ie1)
                        iec2=AMPEB2C(ie2)
                        iec3=AMPEB2C(ie3)
                        jec1=AMPEB2C(je1)
                        jec2=AMPEB2C(je2)
                        jec3=AMPEB2C(je3)
C  Get primitive Electron Basis Function Contraction Coefficients 
                        Cof_ie1=AGEBFCC(ie1)
                        Cof_ie2=AGEBFCC(ie2)
                        Cof_ie3=AGEBFCC(ie3)
                        Cof_je1=AGEBFCC(je1)
                        Cof_je2=AGEBFCC(je2)
                        Cof_je3=AGEBFCC(je3)
C  Get Nuclear Basis Function Contraction Coefficients
                        Cof_ip=AGNBFCC(ip)
                        Cof_jp=AGNBFCC(jp)
C  Map the 8-index contracted integral to 1-D:
                        call index_GAM_3PK(nebf,npbf,
     x                                     ip,jp,
     x                                     iec1,jec1,
     x                                     iec2,jec2,
     x                                     iec3,jec3,ia)
c                 write(*,*)'DBG NEBF= NPBF=  MYID=',nebf,npbf,myid
c                 write(*,*)'DBG IP  = JP  =  MYID=',ip,jp,myid
c                 write(*,*)'DBG IEC1= JEC1=  MYID=',iec1,jec1,myid
c                 write(*,*)'DBG IEC2= JEC2=  MYID=',iec2,jec2,myid
c                 write(*,*)'DBG IEC2= JEC2=  MYID=',iec3,jec3,myid
c                 write(*,*)'DBG MYID=  IA=  ANS=',myid,ia,ans

                        GAM_3(ia)=GAM_3(ia)+ans
     x                            *Cof_ip*Cof_jp
     x                            *Cof_ie1*Cof_je1
     x                            *Cof_ie2*Cof_je2
     x                            *Cof_ie3*Cof_je3
  
c     end do
c                 end do
c                 end do
c              end do
c              end do
c           end do
c           end do
c        end do
c        end do
         end do
!$omp end do
!$omp end parallel      
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

C---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime() - wtime
      write(*,*)'TIME TO CALCULATE GAM_3 INTEGRALS: ',wtime
C---OPENMP-TIMING------------------------------------------------------)

c     write(*,*)'AFTER LOOP, BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'AFTER LOOP, AFTER MPI BARRIER MYID=',myid

c     call mpi_final(ierr)
c     stop
c     if(myid.eq.0) then
c        write(*,*)'BEFORE MPIREDUCE XGAM_2 is:'
c        write(*,*) XXGAM_2
c     end if

c     if(myid.eq.0) then
c        write(*,*)'MYID=',myid
c        write(*,*)'XGAM3 is:',XGAM_3
c     end if
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     if(myid.eq.1) then
c        write(*,*)'MYID=',myid
c        write(*,*)'XGAM3 is:',XGAM_3
c     end if

c     call MPI_REDUCE(XGAM_3,GAM_3,ng3,
c    x                mpi_double_precision,mpi_sum,0,
c    x                MPI_COMM_WORLD,ierr)
c     write(*,*)'after mpi_reduce ierr=',ierr
c     XGAM_3=XXGAM_3
c     if(myid.eq.0) then
c        write(*,*)'AFTER MPIREDUCE XGAM_3 is:'
c        write(*,*) XGAM_3
c     end if

c     write(*,*)'16 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'16 AFTER MPI BARRIER MYID=',myid

c     IF(myid .eq. 0) THEN
CCWS-IO
         write(*,*)
         write(*,*)'FINISHED CALCULATING CONTRACTED GAM_3 INTEGRALS'
         write(*,*)'BEGINNING SYMMETRIZATION OF GAM_3 INTEGRALS'
         write(*,*)
c        write(*,*)'GAM3 is now:',GAM_3
c        write(*,*)

         open(805,file='GAM_3.ufm',form='unformatted',
     x    status='unknown',access='direct',RECL=8)

      do ip=1,npbf
      do jp=1,npbf
         do ie1=1,nebf
         do je1=1,nebf
            do ie2=1,nebf
            do je2=1,nebf
               do ie3=1,nebf
               do je3=1,nebf

C  GAM_2 Symmetrization:
C  Determine packing indices for XGAM_3 integral matrices
C
C              As Packed-->       XGAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je1,ie2,je2,ie3,je3,ia)
                             ia_123=ia
                               x123=GAM_3(ia_123)

C              As Packed-->       XGAM_3(je2,ie2,je3,ie3,je1,ie1,jp,ip)
c                           XGAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie1,je1,ie3,je3,ie2,je2,ia)
                             ia_132=ia
                               x132=GAM_3(ia_132)

C              As Packed-->       XGAM_3(je3,ie3,je1,ie1,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie2,je2,ie1,je1,ie3,je3,ia)
                             ia_213=ia
                               x213=GAM_3(ia_213)

C              As Packed-->       XGAM_3(je1,ie1,je3,ie3,je2,ie2,jp,ip)
c                           XGAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie2,je2,ie3,je3,ie1,je1,ia)
                             ia_231=ia
                               x231=GAM_3(ia_231)

C              As Packed-->       XGAM_3(je2,ie2,je1,ie1,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie3,je3,ie1,je1,ie2,je2,ia)
                             ia_312=ia
                               x312=GAM_3(ia_312)

C              As Packed-->       XGAM_3(je1,ie1,je2,ie2,je3,ie3,jp,ip)
c                           XGAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1) 
      call index_GAM_3PK(nebf,npbf,ip,jp,ie3,je3,ie2,je2,ie1,je1,ia)
                             ia_321=ia
                               x321=GAM_3(ia_321)

                       xxxx=(x123+x132+x213+x231+x312+x321)/six

c                      write(*,*)'xxxx=',xxxx
cc                     GAM_3(ip,jp,ie1,je1,ie2,je2,ie3,je3)=xxxx 
c                      GAM_3PK(ia_123)=xxxx 
c                      call put_GAM3(ia_123,ng3,xxxx)
                       write(805,REC=ia_123) xxxx

cc                     GAM_3(ip,jp,ie1,je1,ie3,je3,ie2,je2)=xxxx
c                      GAM_3PK(ia_132)=xxxx 
c                      call put_GAM3(ia_132,ng3,xxxx)
                       write(805,REC=ia_132) xxxx

cc                     GAM_3(ip,jp,ie2,je2,ie1,je1,ie3,je3)=xxxx
c                      GAM_3PK(ia_213)=xxxx 
c                      call put_GAM3(ia_213,ng3,xxxx)
                       write(805,REC=ia_213) xxxx

cc                     GAM_3(ip,jp,ie2,je2,ie3,je3,ie1,je1)=xxxx
c                      GAM_3PK(ia_231)=xxxx 
c                      call put_GAM3(ia_231,ng3,xxxx)
                       write(805,REC=ia_231) xxxx

cc                     GAM_3(ip,jp,ie3,je3,ie1,je1,ie2,je2)=xxxx
c                      GAM_3PK(ia_312)=xxxx 
c                      call put_GAM3(ia_312,ng3,xxxx)
                       write(805,REC=ia_312) xxxx

cc                     GAM_3(ip,jp,ie3,je3,ie2,je2,ie1,je1)=xxxx
c                      GAM_3PK(ia_321)=xxxx 
c                      call put_GAM3(ia_321,ng3,xxxx)
                       write(805,REC=ia_321) xxxx


               end do
               end do
            end do
            end do
         end do
         end do
      end do
      end do

      close(805)
      write(*,*)'ALL DONE WITH SYMMETRIZATION'
C ENDIF OF MYID=0: 
c     END IF

c     write(*,*)'16 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'16 AFTER MPI BARRIER MYID=',myid

      return
      end

C======================================================================
      subroutine RXCHF_xcalc_GAM3_MD(I1,J1,K1,A1,Amat1,
     x                         I2,J2,K2,A2,Amat2,
     x                         I3,J3,K3,A3,Amat3,
     x                         I4,J4,K4,A4,Amat4,
     x                         L1,M1,N1,B1,Bmat1,
     x                         L2,M2,N2,B2,Bmat2,
     x                         L3,M3,N3,B3,Bmat3,
     x                         L4,M4,N4,B4,Bmat4,
     x                         nat,ngtg1,
     x                         pmass,cat,zan,
     x                         bcoef1,gamma1,
     x                         ans1,ans2,ans3)

C======================================================================
      implicit none

C Input Variables
      integer nat
      integer ngtg1
      double precision pmass
      double precision zan(nat)
      double precision cat(3,nat)
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)

      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1,Amat1(3) 
      double precision A2,Amat2(3) 
      double precision A3,Amat3(3) 
      double precision A4,Amat4(3) 
      double precision B1,Bmat1(3) 
      double precision B2,Bmat2(3) 
      double precision B3,Bmat3(3) 
      double precision B4,Bmat4(3) 

C Variables Returned
      double precision ans1,ans2,ans3

C Local Variables
      integer iii    ! Index for looping over natoms
      integer ik,il  ! Indices for geminal loops
      integer iat

      double precision gamA14
      double precision gamA24
      double precision gamA34
      double precision gamB14
      double precision gamB24
      double precision gamB34
      double precision gamA
      double precision gamB

      double precision cmat(3)
      double precision znuc

C      double precision gVEE
C      double precision xgVEE
      double precision xx,yy,zz
      double precision zero,half,one,two,four
      parameter(zero=0.0d+00,one=1.0d+00,two=2.0d+00,four=4.0d+00)
      parameter(half=0.5d+00)

      double precision gHEg
      double precision gVEEg1
      double precision gVEEg2
      double precision gVEEg3
      double precision xgHEg
      double precision xgVEEg1 
      double precision xgVEEg2 
      double precision xgVEEg3 

      double precision xmass,coulomb_sign
      double precision xke,Vc
      double precision val_vec 
    
      integer NQUAD_coul
      integer NQUAD_ovlap


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     BASIS FUNCTIONS: ASSIGN CENTERS, EXPONENTS, ANGULAR MOM.
C

C     *****ie1 ::  electron 1 bra *****
C     *****ie2 ::  electron 2 bra *****      
C     *****ie3 ::  electron 3 bra *****      
C     *****ip  ::  proton bra     *****

C     *****je1 ::  electron 1 ket *****
C     *****je2 ::  electron 2 ket *****
C     *****je3 ::  electron 3 ket *****
C     *****jp  ::  proton ket     *****

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C      gVEE=zero

C      DO IK=1,NGTG1

C         gamA=gamma1(ik)
C         gamB=zero

C  xgVEE
C  --- g(1,p)VEE(2,3)---
C         call G2_MD_xggs(I1,J1,K1,A1,Amat1,
C     1                   I4,J4,K4,A4,Amat4,
C     2                   L1,M1,N1,B1,Bmat1,
C     3                   L4,M4,N4,B4,Bmat4,
C     4                   gamA,gamB,xx)

c        call pgiovlap(I1,J1,K1,A1,Amat1,
c    x                 I4,J4,K4,A4,Amat4,
c    x                 L1,M1,N1,B1,Bmat1,
c    x                 L4,M4,N4,B4,Bmat4,
c    x                 gamA,gamB,xx)


C         call gfvee(I2,J2,K2,A2,Amat2,
C     x              I3,J3,K3,A3,Amat3,
C     x              L2,M2,N2,B2,Bmat2,
C     x              L3,M3,N3,B3,Bmat3,
C     x              yy)

C         call underflow(xx)
C         call underflow(yy)

C         xgVEE=xx*yy
C         gVEE=gVEE+(bcoef1(ik)*xgVEE)

C End 1 gamma loop
C      end do
C Begin 2 gamma loop

      gHEg=zero
      gVEEg1=zero
      gVEEg2=zero
      gVEEg3=zero

      DO IK=1,NGTG1
         DO IL=1,NGTG1

C>>>>>>>>>>>>>>>>>>>>  xgHEg <<<<<<<<<<<<<<<<<<<<
C Factors as <g(1)g(2)><h^e(3)>
c           ndim=3
c           natom=nat
c           xmass=one
c           coulomb_sign=-one
c           call o3_Hcore_val(NDIM,NATOM,
c    x                        xmass,zan,cmat,
c    x                        coulomb_sign,
c    x                        I3,J3,K3,A3,Amat3,
c    x                        L3,M3,N3,B3,Bmat3,
c    x                        xx)
CCCCCC-rather than call o3_Hcore_val evaluate 
C everything right here...
C
C-RXCHFexch--- T^e(e3) --
C           -- T^e(3)  --RXCHFexch-

            xmass=one
            coulomb_sign=-one

            call gfke(I3,J3,K3,A3,Amat3,
     x                L3,M3,N3,B3,Bmat3,
     x                xmass,xke)
C
C-RXCHFexch--- V^{eC}(e3) --
C           -- V^{eC}(3)  --RXCHFexch-
C
            Vc = ZERO
            do iat=1,nat
                 Cmat(1)=cat(1,iat)
                 Cmat(2)=cat(2,iat)
                 Cmat(3)=cat(3,iat)
                 call gfvec(I3,J3,K3,A3,Amat3,
     x                      L3,M3,N3,B3,Bmat3,
     x                      Cmat,val_vec)
                 Vc = Vc + (zan(iat)*val_vec)
             end do
c            hcore = xke + (Vc*coulomb_sign)
             xx = xke + (Vc*coulomb_sign)
CCCCCC

c           gamA14 = gamma1(ik)
c           gamA24 = ZERO
c           gamA34 = ZERO
c           gamB14 = ZERO
c           gamB24 = gamma1(il)
c           gamB34 = ZERO

            gamA = gamma1(ik)
            gamB = gamma1(il)

C-RXCHFexch---g(e1,p1) g(e2,p1)--
C           --gA(1,4) gB(2,4)--RXCHFexch-

c           call G4ovlap_typ1(I1,J1,K1,A1,Amat1,
c           call G3ovlap(I1,J1,K1,A1,Amat1,
c    x                   I2,J2,K2,A2,Amat2,
c    x                   I4,J4,K4,A4,Amat4,
c    x                   L1,M1,N1,B1,Bmat1,
c    x                   L2,M2,N2,B2,Bmat2,
c    x                   L4,M4,N4,B4,Bmat4,
c    x                   ZERO,gamA,ZERO,
c    x                   ZERO,ZERO,gamB,yy)
cc   4                gamA12,gamA13,gamA23,
cc   4                gamB12,gamB13,gamB23,sval)
            call G3_MD_xggs(I1,J1,K1,A1,Amat1,
     x                      I2,J2,K2,A2,Amat2,
     x                      I4,J4,K4,A4,Amat4,
     x                      L1,M1,N1,B1,Bmat1,
     x                      L2,M2,N2,B2,Bmat2,
     x                      L4,M4,N4,B4,Bmat4,
     x                      ZERO,gamA,ZERO,     ! gam13 = g(1,4)
     x                      ZERO,ZERO,gamB,yy)  ! gam23 = g(2,4)

C-RXCHFexch------g(e1,p1) V^{ep}(e3,p1) g(e2,p1)--
C              --gA(1,4) V^{ep}(3,4) gB(2,4)--RXCHFexch-

            coulomb_sign  = -ONE
            gamA14 = gamma1(ik)
            gamA24 = ZERO
            gamA34 = ZERO
            gamB14 = ZERO
            gamB24 = gamma1(il)
            gamB34 = ZERO
!           call rys_G4vee_r34(NQUAD_coul,NQUAD_ovlap,
c           call interface_G4vee_r34(NQUAD_coul,NQUAD_ovlap,
c    x                               I1,J1,K1,A1,Amat1,
c    x                               I2,J2,K2,A2,Amat2,
c    x                               I3,J3,K3,A3,Amat3,
c    x                               I4,J4,K4,A4,Amat4,
c    x                               L1,M1,N1,B1,Bmat1,
c    x                               L2,M2,N2,B2,Bmat2,
c    x                               L3,M3,N3,B3,Bmat3,
c    x                               L4,M4,N4,B4,Bmat4,
c    x                               gamA14,gamA24,gamA34,
c    x                               gamB14,gamB24,gamB34,zz)
            call G4_MD_xgVepg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        I4,J4,K4,A4,Amat4,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gamA14,gamB24,
     *                        zz)
c           call G4Vep_AUX_g14g24V34(I1,J1,K1,A1,Amat1,
c    *                               I2,J2,K2,A2,Amat2,
c    *                               I3,J3,K3,A3,Amat3,
c    *                               I4,J4,K4,A4,Amat4,
c    *                               L1,M1,N1,B1,Bmat1,
c    *                               L2,M2,N2,B2,Bmat2,
c    *                               L3,M3,N3,B3,Bmat3,
c    *                               L4,M4,N4,B4,Bmat4,
c    *                               gamA14,zero,
c    *                               zero,gamB24,
c    *                               zz)

             call underflow(xx)
             call underflow(yy)
             call underflow(zz)

             xgHEg = ((xx*yy) + (zz*coulomb_sign) )


C>>>>>>>>>>>>>>>>>>>>  xgVEEg1 <<<<<<<<<<<<<<<<<<<<
C  --- g(1,p)VEE(2,3)g(1,p)---
C Factors as <g(1)g(1)><V^{ee}(2,3)>
            gamA=gamma1(ik)
            gamB=gamma1(il)
            call G2_MD_xggs(I1,J1,K1,A1,Amat1,
     1                      I4,J4,K4,A4,Amat4,
     2                      L1,M1,N1,B1,Bmat1,
     3                      L4,M4,N4,B4,Bmat4,
     4                      gamA,gamB,xx)

c           call pgiovlap(I1,J1,K1,A1,Amat1,
c    x                    I4,J4,K4,A4,Amat4,
c    x                    L1,M1,N1,B1,Bmat1,
c    x                    L4,M4,N4,B4,Bmat4,
c    x                    gamA,gamB,xx)


            call gfvee(I2,J2,K2,A2,Amat2,
     x                 I3,J3,K3,A3,Amat3,
     x                 L2,M2,N2,B2,Bmat2,
     x                 L3,M3,N3,B3,Bmat3,
     x                 yy)


            call underflow(xx)
            call underflow(yy)

            xgVEEg1=xx*yy

C>>>>>>>>>>>>>>>>>>>>  xgVEEg2 <<<<<<<<<<<<<<<<<<<<
C ans = <GA(1)GA(2)GA(3)GB(4)|g(1,4)g(2,4)g(3,4)/(r1-r2)|GB(1)GB(2)GB(3)GB(4)>
C ans = <ie1 ie2 ie3 ip|g(1,p)g(2,p)g(3,p)/(r1-r2)|je1 je2 je3 jp>
            gamA14=gamma1(ik)
            gamA24=zero
            gamA34=zero
            gamB14=zero
            gamB24=zero
            gamB34=gamma1(il)
CCWS- 11-08-2010(:  There is a problem with the rys_G4Vee_r12
C                   AC routine, so call a CS routine to 
C                   evaluate the xgVEEg2 integral:
c           NQUAD_coul=5
c           NQUAD_ovlap=5
c           call rys_G4Vee_r12(NQUAD_coul,NQUAD_ovlap,
c    x                              I1,J1,K1,A1,Amat1,
c    x                              I2,J2,K2,A2,Amat2,
c    x                              I3,J3,K3,A3,Amat3,
c    x                              I4,J4,K4,A4,Amat4,
c    x                              L1,M1,N1,B1,Bmat1,
c    x                              L2,M2,N2,B2,Bmat2,
c    x                              L3,M3,N3,B3,Bmat3,
c    x                              L4,M4,N4,B4,Bmat4,
c    x                              gamA14,gamA24,gamA34,
c    x                              gamB14,gamB24,gamB34,
c    x                              xgVEEg2)
C
C-RXCHFexch------g(e1,p1) V^{ee}(e2,e3) g(e3,p1)--
C              --gA(1,4) V^{ee}(2,3) gB(3,4)--
C              --gA(3,4) V^{ee}(2,1) gB(1,4)--
C              --gB(3,4) V^{ee}(2,1) gA(1,4)--
C              --gA(1,4) V^{ee}(2,1) gB(3,4)--RXCHFexch-
C
            call G4_MD_xgVeeg(I3,J3,K3,A3,Amat3,
     *                        I2,J2,K2,A2,Amat2,
     *                        I1,J1,K1,A1,Amat1,
     *                        I4,J4,K4,A4,Amat4,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gamA14,gamB34,
     *                        xgVEEg2)
CCWS- 11-08-2010 )
c           call G4Vee_AUX_g14g34V12(I1,J1,K1,A1,Amat1,
c    *                               I2,J2,K2,A2,Amat2,
c    *                               I3,J3,K3,A3,Amat3,
c    *                               I4,J4,K4,A4,Amat4,
c    *                               L1,M1,N1,B1,Bmat1,
c    *                               L2,M2,N2,B2,Bmat2,
c    *                               L3,M3,N3,B3,Bmat3,
c    *                               L4,M4,N4,B4,Bmat4,
c    *                               gamA14,gamB34,
c    *                               xgVEEg2)

            gamA14=gamma1(ik)
            gamA24=zero
            gamA34=zero
            gamB14=zero
            gamB24=zero
            gamB34=gamma1(il)

C-RXCHFexch------g(e1,p1) V^{ee}(e1,e3) g(e2,p1)--
C              --gA(1,4) V^{ee}(1,3) gB(2,4)--
C              --gA(1,4) V^{ee}(1,2) gB(3,4)--RXCHFexch-

            call G4_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                        I3,J3,K3,A3,Amat3,
     *                        I2,J2,K2,A2,Amat2,
     *                        I4,J4,K4,A4,Amat4,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L4,M4,N4,B4,Bmat4,
     *                        gamA14,gamB34,
     *                        xgVEEg3)


C  Sum terms against bcoeff
            gHEg  =gHEg  +(bcoef1(ik)*bcoef1(il)*xgHEg)
            gVEEg1=gVEEg1+(bcoef1(ik)*bcoef1(il)*xgVEEg1)
            gVEEg2=gVEEg2+(bcoef1(ik)*bcoef1(il)*xgVEEg2)
            gVEEg3=gVEEg3+(bcoef1(ik)*bcoef1(il)*xgVEEg3)

C End 2 gamma loop
         end do
      end do


C Total integral build
      ans1 = gVEEg1*half
      ans2 = gHEg + gVEEg3*half
      ans3 = gVEEg2*half

CCWS-debug
c                    write(*,*)'========='
c                    write(*,*)'ie1=',ie1
c                    write(*,*)'je1=',je1
c                    write(*,*)'ie2=',ie2
c                    write(*,*)'je2=',je2
c                    write(*,*)'ie3=',ie3
c                    write(*,*)'je3=',je3
c                    write(*,*)'gHEg   =',gHEg
c                    write(*,*)'gVEEg1 =',gVEEg1
c                    write(*,*)'gVEEg2 =',gVEEg2
c                    write(*,*)'ans    =',ans
c                    write(*,*)'========='
CCWS-debug



      return
      end

