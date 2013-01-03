C======================================================================
      subroutine GAM2_OMP_MD(nebf,npebf,npbf,ng2,ng2prm,nat,ngtg1,
     x                       pmass,cat,zan,bcoef1,gamma1,
     x                       AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC)
c    x                       XGAM_2,XGAM_2s)

C======================================================================
      implicit none
      include 'omp_lib.h'
c     include 'mpif.h'

C Input Variables
      integer ng2    ! Number contracted gamma2 integrals
      integer ng2prm ! Number primitive gamma2 integrals
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
      integer ia,ia1,ia2
      integer ip,jp,ie1,je1,ie2,je2
      integer iec1,jec1,iec2,jec2
      double precision Cof_ie1,Cof_je1
      double precision Cof_ie2,Cof_je2
      double precision Cof_ip,Cof_jp

      double precision zero,half
      parameter(zero=0.0d+00,half=0.5d+00)
      double precision ans
      double precision ansE
      double precision ansS
      double precision x12
      double precision x21
c     double precision val_gam2s
C-------------------------------------
      double precision,allocatable :: XGAM_2(:)
      double precision,allocatable :: XGAM_2s(:)
c     double precision XXGAM_2(ng2)
c     double precision XXGAM_2s(ng2)
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
      integer istat
C--------------------------------(
C MPI-Related Local variables
c     integer myid,ierr,nprocs
c     integer loopi,iLP
c     integer istart,iend
c     integer loop_map(ng2prm,6)
C--------------------------------)
C---OPENMP-RELATED-VARIABLES-----(
      integer id
      integer loopi,iLP
      double precision wtime
      integer,allocatable :: loop_map(:,:)
C---OPENMP-RELATED-VARIABLES-----)

CCWS_int_stat
c     integer iii
c     double precision tol
CCWS_int_stat

c     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
c     call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(allocated(XGAM_2)) deallocate(XGAM_2)
      allocate( XGAM_2(ng2),stat=istat )

      if(allocated(XGAM_2s)) deallocate(XGAM_2s)
      allocate( XGAM_2s(ng2),stat=istat )

      if(allocated(loop_map)) deallocate(loop_map)
      allocate( loop_map(ng2prm,6),stat=istat )

C---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime()
C---OPENMP-TIMING------------------------------------------------------)

C-----------INITIALIZE-DATA-STRUCTURES-ON-MASTER-----------------------(
c     IF(myid .eq. 0) THEN
         write(*,*)
         write(*,*)'**************************************'
         write(*,*)'         GAM_2 Integrals    '
         write(*,*)
         write(*,*)'METHOD   = MD'
         write(*,*)'npebf    =',npebf
         write(*,*)'nebf     =',nebf
         write(*,*)'npbf     =',npbf
         write(*,*)'nprim_g2 =',ng2prm
         write(*,*)'ng2      =',ng2
c        write(*,*)'NPROC    =',nprocs
         write(*,*)' Available processors: ',omp_get_num_procs()
         write(*,*)' Available threads     ',omp_get_max_threads()
         write(*,*)' Threads in use        ',omp_get_num_threads()
         write(*,*)'**************************************'
         write(*,*)
c     end if
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     IF(myid .eq. 0) THEN
C Initialize 1-D arrays to hold the contracted integrals
         do ia=1,ng2
c           GAM_2(ia)=zero
c           GAM_2s(ia)=zero
            XGAM_2(ia)=zero
            XGAM_2s(ia)=zero
c           XXGAM_2(ia)=zero
c           XXGAM_2s(ia)=zero
         end do
c     end if

c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c     IF(myid .eq. 0) THEN
c     write(*,*)'GOT HERE 1 MYID=',myid
C Compress nested loops 
         Loopi=0
         do ip=1,npbf
            do jp=1,npbf
               do ie1=1,npebf
                  do je1=1,npebf
                     do ie2=1,npebf
                        do je2=1,npebf

                           Loopi=Loopi+1
                           loop_map(Loopi,1)=je2
                           loop_map(Loopi,2)=ie2
                           loop_map(Loopi,3)=je1
                           loop_map(Loopi,4)=ie1
                           loop_map(Loopi,5)=jp
                           loop_map(Loopi,6)=ip

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
c     END IF  ! end if for myid
C-----------INITIALIZE-DATA-STRUCTURES-ON-MASTER-----------------------)

c     write(*,*)'1 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(GAM_2,ng2,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(XGAM_2,ng2,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(XXGAM_2,ng2,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(GAM_2s,ng2,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(XGAM_2s,ng2,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(XXGAM_2s,ng2,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     write(*,*)'1 AFTER MPI BARRIER MYID=',myid
c     call MPI_BCAST(loop_map,ng2prm*6,mpi_integer,0,
c    x     MPI_COMM_WORLD,ierr)
c     write(*,*)'BCAST LoopMap IERR=',ierr
c     write(*,*)'2 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'2 AFTER MPI BARRIER MYID=',myid


c     write(*,*)'16 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'16 AFTER MPI BARRIER MYID=',myid
c     call loop_size(1,ng2prm,nprocs,myid,istart,iend)
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
c     call mpi_final(ierr)
c     stop

C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf)
!$ompx shared(ng2)
!$ompx shared(ng2prm)
!$ompx shared(XGAM_2,XGAM_2s)
!$ompx private(iLp) 
!$ompx private(ia)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(iec1,jec1)
!$ompx private(iec2,jec2)
!$ompx private(A1,I1,J1,K1,Amat1)
!$ompx private(A2,I2,J2,K2,Amat2)
!$ompx private(A3,I3,J3,K3,Amat3)
!$ompx private(B1,L1,M1,N1,Bmat1)
!$ompx private(B2,L2,M2,N2,Bmat2)
!$ompx private(B3,L3,M3,N3,Bmat3)
!$ompx private(Cof_ie1,Cof_je1)
!$ompx private(Cof_ie2,Cof_je2)
!$ompx private(Cof_ip,Cof_jp)
!$ompx private(ansE,ansS)
!$ompx private(id)
CCCCCCCCCCCCCCCCCCCCCCmpx reduction(+:XGAM_2,XGAM_2s)

      id= omp_get_thread_num()
      write(*,*)' Hello from process ',id
      if(id.eq.0) then
         write(*,*)'Threads in use', omp_get_num_threads()
      end if

!$omp do
c     istart=1
c     iend=2
c     do iLP=istart,iend
      do iLP=1,ng2prm
c     do ip=1,npbf
c        do jp=1,npbf
c           do ie1=1,npebf
c              do je1=1,npebf
c                 do ie2=1,npebf
c                    do je2=1,npebf

c        write(*,*)'DBG MYID=',myid
C  Map loop indices
         je2=loop_map(iLP,1)
         ie2=loop_map(iLP,2)
         je1=loop_map(iLP,3)
         ie1=loop_map(iLP,4)
         jp =loop_map(iLP,5)
         ip =loop_map(iLP,6)
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

         A3=NUCEX(ip)
         I3=NUCAM(ip,1)
         J3=NUCAM(ip,2)
         K3=NUCAM(ip,3)
         Amat3(1)=NUCBFC(ip,1)
         Amat3(2)=NUCBFC(ip,2)
         Amat3(3)=NUCBFC(ip,3)

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

         B3=NUCEX(jp)
         L3=NUCAM(jp,1)
         M3=NUCAM(jp,2)
         N3=NUCAM(jp,3)
         Bmat3(1)=NUCBFC(jp,1)
         Bmat3(2)=NUCBFC(jp,2)
         Bmat3(3)=NUCBFC(jp,3)

C  Calculate primitive integral
c                       call xcalc_GAM2_MD(ip,jp,
c    x                                     ie1,je1,
c    x                                     ie2,je2,
c    x                                     nat,ngtg1,
c    x                                     pmass,cat,zan,
c    x                                     bcoef1,gamma1,
c    x                                     ansE,ansS)
         call xcalc_GAM2_MD(I1,J1,K1,A1,Amat1,
     x                      I2,J2,K2,A2,Amat2,
     x                      I3,J3,K3,A3,Amat3,
     x                      L1,M1,N1,B1,Bmat1,
     x                      L2,M2,N2,B2,Bmat2,
     x                      L3,M3,N3,B3,Bmat3,
     x                      nat,ngtg1,
     x                      pmass,cat,zan,
     x                      bcoef1,gamma1,
     x                      ansE,ansS)
c                       ansE=1.0d-08
c                       ansS=1.0d-08

c                       call underflow(ansE)
c                       call underflow(ansS)

C  Map from primitive BF indices to contracted indices
c                       call MPEB2C(ie1,iec1)
c                       call MPEB2C(ie2,iec2)
c                       call MPEB2C(je1,jec1)
c                       call MPEB2C(je2,jec2)
                        iec1=AMPEB2C(ie1)
                        iec2=AMPEB2C(ie2)
                        jec1=AMPEB2C(je1)
                        jec2=AMPEB2C(je2)
C  Get primitive Electron Basis Function Contraction Coefficients 
c                       call GEBFCC(ie1,Cof_ie1)
c                       call GEBFCC(ie2,Cof_ie2)
c                       call GEBFCC(je1,Cof_je1)
c                       call GEBFCC(je2,Cof_je2)
                        Cof_ie1=AGEBFCC(ie1)
                        Cof_ie2=AGEBFCC(ie2)
                        Cof_je1=AGEBFCC(je1)
                        Cof_je2=AGEBFCC(je2)
C  Get Nuclear Basis Function Contraction Coefficients
c                       call GNBFCC(ip,Cof_ip)
c                       call GNBFCC(jp,Cof_jp)
                        Cof_ip=AGNBFCC(ip)
                        Cof_jp=AGNBFCC(jp)
C  Map the 6-index contracted integral to 1-D:
                        call index_GAM_2PK(nebf,npbf,
     x                                     ip,jp,
     x                                     iec1,jec1,
     x                                     iec2,jec2,ia)
c                 write(*,*)'DBG NEBF= NPBF=  MYID=',nebf,npbf,myid
c                 write(*,*)'DBG IP  = JP  =  MYID=',ip,jp,myid
c                 write(*,*)'DBG IEC1= JEC1=  MYID=',iec1,jec1,myid
c                 write(*,*)'DBG IEC2= JEC2=  MYID=',iec2,jec2,myid
c                 write(*,*)'DBG IA  =        MYID=',ia,myid

                        XGAM_2(ia)=XGAM_2(ia)+ansE
     x                            *Cof_ip*Cof_jp
     x                            *Cof_ie1*Cof_je1
     x                            *Cof_ie2*Cof_je2
  
                        XGAM_2s(ia)=XGAM_2s(ia)+ansS
     x                             *Cof_ip*Cof_jp
     x                             *Cof_ie1*Cof_je1
     x                             *Cof_ie2*Cof_je2


c                    end do
c                 end do
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
      write(*,*)'TIME TO CALCULATE GAM_2 INTEGRALS: ',wtime
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

c     call MPI_REDUCE(XXGAM_2,XGAM_2,ng2,
c    x                mpi_double_precision,mpi_sum,0,
c    x                MPI_COMM_WORLD,ierr)
c     write(*,*)'after mpi_reduce ierr=',ierr
c     XGAM_2=XXGAM_2
c     if(myid.eq.0) then
c        write(*,*)'AFTER MPIREDUCE XGAM_2 is:'
c        write(*,*) XGAM_2
c     end if

c     write(*,*)'16 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'16 AFTER MPI BARRIER MYID=',myid
c     call MPI_REDUCE(XXGAM_2s,XGAM_2s,ng2,
c    x                mpi_double_precision,mpi_sum,0,
c    x                MPI_COMM_WORLD,ierr)

c     write(*,*)'16 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'16 AFTER MPI BARRIER MYID=',myid

c     IF(myid .eq. 0) THEN
CCWS-IO
         open(803,file='GAM_2.ufm',form='unformatted',
     x    status='unknown',access='direct',RECL=8)

         open(804,file='GAM_2s.ufm',form='unformatted',
     x    status='unknown',access='direct',RECL=8)


      do ip=1,npbf
         do jp=1,npbf
            do iec1=1,nebf
               do jec1=1,nebf
                  do iec2=1,nebf
                     do jec2=1,nebf

C  GAM_2 Symmetrization:
                        call index_GAM_2PK(nebf,npbf,
     x                                     ip,jp,
     x                                     iec1,jec1,
     x                                     iec2,jec2,ia1)

                        call index_GAM_2PK(nebf,npbf,
     x                                     ip,jp,
     x                                     iec2,jec2,
     x                                     iec1,jec1,ia2)

cc                      x12=XGAM_2(ip,jp,iec1,jec1,iec2,jec2)
cc                      x21=XGAM_2(ip,jp,iec2,jec2,iec1,jec1)
                        x12=XGAM_2(ia1)
                        x21=XGAM_2(ia2)
                        ans=half*(x12+x21)
                        write(803,REC=ia1) ans
c                       GAM_2(ia1)=half*(x12+x21)

cc                      x12=XGAM_2s(ip,jp,iec1,jec1,iec2,jec2)
cc                      x21=XGAM_2s(ip,jp,iec2,jec2,iec1,jec1)
                        x12=XGAM_2s(ia1)
                        x21=XGAM_2s(ia2)
                        ans=half*(x12+x21)
                        write(804,REC=ia1) ans
c                       GAM_2s(ia1)=half*(x12+x21)


                     end do
                  end do
               end do
            end do
         end do
      end do

      close(803)
      close(804)
C ENDIF OF MYID=0: 
c     END IF

c     write(*,*)'16 BEFORE MPI BARRIER MYID=',myid
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(*,*)'16 AFTER MPI BARRIER MYID=',myid

      if(allocated(loop_map)) deallocate(loop_map)
      if(allocated(XGAM_2s)) deallocate(XGAM_2s)
      if(allocated(XGAM_2)) deallocate(XGAM_2)


      return
      end

C======================================================================
      subroutine xcalc_GAM2_MD(I1,J1,K1,A1,Amat1,
     x                         I2,J2,K2,A2,Amat2,
     x                         I3,J3,K3,A3,Amat3,
     x                         L1,M1,N1,B1,Bmat1,
     x                         L2,M2,N2,B2,Bmat2,
     x                         L3,M3,N3,B3,Bmat3,
     x                         nat,ngtg1,
     x                         pmass,cat,zan,
     x                         bcoef1,gamma1,
     x                         ansE,ansS)

c     subroutine xcalc_GAM2_MD(ip,jp,
c    x                         ie1,je1,
c    x                         ie2,je2,
c    x                         nat,ngtg1,
c    x                         pmass,cat,zan,
c    x                         bcoef1,gamma1,
c    x                         ansE,ansS)
C======================================================================
      implicit none

C Input Variables
c     integer ie1,ie2,ip
c     integer je1,je2,jp
      integer nat
      integer ngtg1
      double precision pmass
      double precision zan(nat)
      double precision cat(3,nat)
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)

C Variables Returned
      double precision ansE
      double precision ansS

C Local Variables
      integer iii    ! Index for looping over natoms
      integer ik,il  ! Indices for geminal loops

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

      double precision cmat(3)
      double precision znuc
      double precision vec_ans
      double precision xmass 

      double precision gamA12,gamA13,gamA23 
      double precision gamB12,gamB13,gamB23 

      double precision zero,half,two,four
      parameter(zero=0.0d+00,half=0.5d+00,two=2.0d+00,four=4.0d+00)

      double precision xgTE      
      double precision xgVEC     
      double precision xgVEE     
      double precision xgVEP     

      double precision gTE      
      double precision gVEC     
      double precision gVEE     
      double precision gVEP     

      double precision xgTPg      
      double precision xgTEg1     
      double precision xgTEg2     
      double precision xgTEg3     
      double precision xgVPCg      
      double precision xgVECg1     
      double precision xgVECg2     
      double precision xgVEPg1     
      double precision xgVEPg2     
      double precision xgVEEg1     
      double precision xgVEEg2     
      double precision xgsg     

      double precision gTPg      
      double precision gTEg1     
      double precision gTEg2     
      double precision gTEg3     
      double precision gVPCg      
      double precision gVECg1     
      double precision gVECg2     
      double precision gVEPg1     
      double precision gVEPg2     
      double precision gVEEg1     
      double precision gVEEg2     
      double precision gsg     

C******
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C     BASIS FUNCTIONS: ASSIGN CENTERS, EXPONENTS, ANGULAR MOM.
C

C     *****ie1 ::  electron 1 bra *****
C     *****ie2 ::  electron 2 bra *****      
C     *****ip  ::  proton bra     *****

C     *****je1 ::  electron 1 ket *****
C     *****je2 ::  electron 2 ket *****
C     *****jp  ::  proton ket     *****

c     call get_BF(1,ie1,I1,J1,K1,A1,Amat1)
c     call get_BF(1,ie2,I2,J2,K2,A2,Amat2)
c     call get_BF(2,ip,I3,J3,K3,A3,Amat3)

c     call get_BF(1,je1,L1,M1,N1,B1,Bmat1)
c     call get_BF(1,je2,L2,M2,N2,B2,Bmat2)
c     call get_BF(2,jp,L3,M3,N3,B3,Bmat3)


CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      gTE =zero
      gVEC=zero
      gVEE=zero
      gVEP=zero

      DO IK=1,NGTG1

C        xgTE  (2,3)
C        --g(e2,p1) T^e(e1)--
         call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
         gamA23=GAMMA1(IK)
         xmass=1.0d+00
         call G3_MD_KE(I1,J1,K1,A1,Amat1,
     *                 I2,J2,K2,A2,Amat2,
     *                 I3,J3,K3,A3,Amat3,
     *                 L1,M1,N1,B1,Bmat1,
     *                 L2,M2,N2,B2,Bmat2,
     *                 L3,M3,N3,B3,Bmat3,
     *                 gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23,
     *                 xmass,xgTE)

C        xgVEC  (2,3)
C        --g(e2,p1) V^{eC}(e1)--
         xgVEC=zero
         DO III=1,NAT
            Cmat(1)=cat(1,III)
            Cmat(2)=cat(2,III)
            Cmat(3)=cat(3,III)
            ZNUC=ZAN(III)
 
            call G3_MD_xgVxCg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        Cmat,ZNUC,
     *                        vec_ans)

            call underflow(vec_ans)
            xgVEC=xgVEC-ZNUC*vec_ans
 
         END DO

C        xgVEP  (2,3)
C        --g(e2,p1) V^{ep}(e1,p1)--
         call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                     I3,J3,K3,A3,Amat3,
     *                     I2,J2,K2,A2,Amat2,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L3,M3,N3,B3,Bmat3,
     *                     L2,M2,N2,B2,Bmat2,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     xgVEP)

         call underflow(xgVEP)
         xgVEP = -xgVEP

C        xgVEE  (2,3)
C        --g(e2,p1) V^{ee}(e1,e2)--
         call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                     I2,J2,K2,A2,Amat2,
     *                     I3,J3,K3,A3,Amat3,
     *                     L1,M1,N1,B1,Bmat1,
     *                     L2,M2,N2,B2,Bmat2,
     *                     L3,M3,N3,B3,Bmat3,
     *                     gamA12,gamA13,gamA23,
     *                     gamB12,gamB13,gamB23,
     *                     xgVEE)
         call underflow(xgVEE)


C  Sum up 1 gamma loop terms
         gTE= gTE +bcoef1(ik)*xgTE
         gVEC=gVEC+bcoef1(ik)*xgVEC
         gVEP=gVEP+bcoef1(ik)*xgVEP
         gVEE=gVEE+bcoef1(ik)*xgVEE

      end do
C --- end of 1 gamma loop ---
C
C --- 2 gamma loop ---
C
      gTPg =zero
      gTEg1=zero
      gTEg2=zero
      gTEg3=zero

      gVPCg =zero
      gVECg1=zero
      gVECg2=zero

      gVEPg1=zero
      gVEPg2=zero
      gVEEg1=zero
      gVEEg2=zero

      gsg = zero

      DO IK=1,NGTG1
         DO IL=1,NGTG1

C%%%%%%%%%%--2 Gamma Kinetic Energy Integrals--%%%%%%%%%%

C           xgTPg  
C           --gA(e1,p1) T^p(p) gB(e2,p1)--
C           --gA(1,3) T^p(3) gB(2,3)--
C           --gA(3,1) T^p(1) gB(2,1)--
C           --gA(2,1) T^p(1) gB(3,1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                 gamB12,gamB13,gamB23)
            gamA12=gamma1(IK)
            gamB13=gamma1(IL)
            xmass=pmass
            call G3_MD_KE(I3,J3,K3,A3,Amat3,
     x                    I1,J1,K1,A1,Amat1,
     *                    I2,J2,K2,A2,Amat2,
     *                    L3,M3,N3,B3,Bmat3,
     *                    L1,M1,N1,B1,Bmat1,
     *                    L2,M2,N2,B2,Bmat2,
     *                    gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23,
     *                    xmass,xgTPg)

C             xgTEg1 (1,3) (2,3)
C-------------g(e1,p1) T^e(e1) g(e2,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     1                    gamB12,gamB13,gamB23)
            gamA13=GAMMA1(IK)
            gamB23=GAMMA1(IL)
            xmass=1.0d+00
            call G3_MD_KE(I1,J1,K1,A1,Amat1,
     *                    I2,J2,K2,A2,Amat2,
     *                    I3,J3,K3,A3,Amat3,
     *                    L1,M1,N1,B1,Bmat1,
     *                    L2,M2,N2,B2,Bmat2,
     *                    L3,M3,N3,B3,Bmat3,
     *                    gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23,
     *                    xmass,xgTEg1)


C             xgTEg2 (2,3) (1,3)
C-------------g(e2,p1) T^e(e1) g(e1,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23)
            gamA23=GAMMA1(IK)
            gamB13=GAMMA1(IL)
            xmass=1.0d+00
            call G3_MD_KE(I1,J1,K1,A1,Amat1,
     *                    I2,J2,K2,A2,Amat2,
     *                    I3,J3,K3,A3,Amat3,
     *                    L1,M1,N1,B1,Bmat1,
     *                    L2,M2,N2,B2,Bmat2,
     *                    L3,M3,N3,B3,Bmat3,
     *                    gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23,
     *                    xmass,xgTEg2)

C             xgTEg3 (2,3) (2,3)
C-------------g(e2,p1)T^e(e1)g(e2,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23)
            gamA23=GAMMA1(IK)
            gamB23=GAMMA1(IL)
            xmass=1.0d+00
            call G3_MD_KE(I1,J1,K1,A1,Amat1,
     *                    I2,J2,K2,A2,Amat2,
     *                    I3,J3,K3,A3,Amat3,
     *                    L1,M1,N1,B1,Bmat1,
     *                    L2,M2,N2,B2,Bmat2,
     *                    L3,M3,N3,B3,Bmat3,
     *                    gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23,
     *                    xmass,xgTEg3)


C%%%%%%%%%%--End of 2 Gamma Kinetic Energy Integrals--%%%%%%%%%%
C
C%%%%%%%%%%--2 Gamma VEC Integrals--%%%%%%%%%%

            xgVPCg =zero
            xgVECg1=zero
            xgVECg2=zero

            DO III=1,NAT
               Cmat(1)=cat(1,III)
               Cmat(2)=cat(2,III)
               Cmat(3)=cat(3,III)
               ZNUC=ZAN(III)
C
C              xgVPCg (1,2)  (1,3)
C              --g(e1,p1) V^{PC}(p1) g(e2,p1)--
C
               call zero_gam(gamA12,gamA13,gamA23,
     *                       gamB12,gamB13,gamB23)
               gamA12=GAMMA1(IK)
               gamB13=GAMMA1(IL)

               call G3_MD_xgVxCg(I3,J3,K3,A3,Amat3,
     x                           I1,J1,K1,A1,Amat1,
     x                           I2,J2,K2,A2,Amat2,
     x                           L3,M3,N3,B3,Bmat3,
     x                           L1,M1,N1,B1,Bmat1,
     x                           L2,M2,N2,B2,Bmat2,
     x                           gamA12,gamA13,gamA23,
     x                           gamB12,gamB13,gamB23,
     x                           Cmat,ZNUC,
     x                           vec_ans)

               call underflow(vec_ans)
               xgVPCg=xgVPCg+ZNUC*vec_ans

C              xgVECg1 (1,3)  (2,3)
C              --g(e1,p1) V^{eC}(e1) g(e2,p1)--
C
               call zero_gam(gamA12,gamA13,gamA23,
     *                       gamB12,gamB13,gamB23)
               gamA13=GAMMA1(IK)
               gamB23=GAMMA1(IL)

               call G3_MD_xgVxCg(I1,J1,K1,A1,Amat1,
     *                           I2,J2,K2,A2,Amat2,
     *                           I3,J3,K3,A3,Amat3,
     *                           L1,M1,N1,B1,Bmat1,
     *                           L2,M2,N2,B2,Bmat2,
     *                           L3,M3,N3,B3,Bmat3,
     *                           gamA12,gamA13,gamA23,
     *                           gamB12,gamB13,gamB23,
     *                           Cmat,ZNUC,
     *                           vec_ans)
               call underflow(vec_ans)
               xgVECg1=xgVECg1-ZNUC*vec_ans

C              xgVECg2 (2,3)  (2,3)
C              --g(e2,p1) V^{eC}(e1) g(e2,p1)--
               call zero_gam(gamA12,gamA13,gamA23,
     *                       gamB12,gamB13,gamB23)
               gamA23=GAMMA1(IK)
               gamB23=GAMMA1(IL)

               call G3_MD_xgVxCg(I1,J1,K1,A1,Amat1,
     *                           I2,J2,K2,A2,Amat2,
     *                           I3,J3,K3,A3,Amat3,
     *                           L1,M1,N1,B1,Bmat1,
     *                           L2,M2,N2,B2,Bmat2,
     *                           L3,M3,N3,B3,Bmat3,
     *                           gamA12,gamA13,gamA23,
     *                           gamB12,gamB13,gamB23,
     *                           Cmat,ZNUC,
     *                           vec_ans)
               call underflow(vec_ans)

               xgVECg2=xgVECg2-ZNUC*vec_ans
C End of loop over NAT (classical nuclei)
            END DO

C%%%%%%%%%%--End of 2 Gamma VEC Integrals--%%%%%%%%%%

C%%%%%%%%%%--2 Gamma VEP Integrals--%%%%%%%%%%

C           xgVEPg1  (1,2) (2,3)
C           --g(e1,p1) V^{ep}(e1,p1) g(e2,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23)
            gamA12=GAMMA1(IK)
            gamB23=GAMMA1(IL)

            call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                        I3,J3,K3,A3,Amat3,
     *                        I2,J2,K2,A2,Amat2,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L2,M2,N2,B2,Bmat2,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        xgVEPg1)
            call underflow(xgVEPg1)
            xgVEPg1=-xgVEPg1

C           xgVEPg2  (2,3) (2,3)
C           --g(e2,p1) V^{ep}(e1,p1) g(e2,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23)
            gamA23=GAMMA1(IK)
            gamB23=GAMMA1(IL)

            call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                        I3,J3,K3,A3,Amat3,
     *                        I2,J2,K2,A2,Amat2,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L3,M3,N3,B3,Bmat3,
     *                        L2,M2,N2,B2,Bmat2,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        xgVEPg2)
            call underflow(xgVEPg2)
            xgVEPg2=-xgVEPg2

C%%%%%%%%%%--End of 2 Gamma VEP Integrals--%%%%%%%%%%

C%%%%%%%%%%--2 Gamma VEE Integrals--%%%%%%%%%%

C           xgVEEg1  (1,3) (1,3)
C           --g(e1,p1) V^{ee}(e1,e2) g(e1,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     *                    gamB12,gamB13,gamB23)
            gamA13=GAMMA1(IK)
            gamB13=GAMMA1(IL)

            call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        xgVEEg1)
            call underflow(xgVEEg1)

C           xgVEEg2  (1,3) (2,3)
C           --g(e1,p1) V^{ee}(e1,e2) g(e2,p1)--
            call zero_gam(gamA12,gamA13,gamA23,
     1                    gamB12,gamB13,gamB23)
            gamA13=GAMMA1(IK)
            gamB23=GAMMA1(IL)

            call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
     *                        I2,J2,K2,A2,Amat2,
     *                        I3,J3,K3,A3,Amat3,
     *                        L1,M1,N1,B1,Bmat1,
     *                        L2,M2,N2,B2,Bmat2,
     *                        L3,M3,N3,B3,Bmat3,
     *                        gamA12,gamA13,gamA23,
     *                        gamB12,gamB13,gamB23,
     *                        xgVEEg2)
            call underflow(xgVEEg2)

c           evaluate overlap
c           --g(re1,rp1)g(re2,rp1)--         
            call zero_gam(gamA12,gamA13,gamA23, 
     *            gamB12,gamB13,gamB23)
            gamA13 = gamma1(ik)
            gamB23 = gamma1(il)
c           call G3ovlap(I1,J1,K1,A1,Amat1,
c    *                   I2,J2,K2,A2,Amat2,
c    *                   I3,J3,K3,A3,Amat3,
c    *                   L1,M1,N1,B1,Bmat1,
c    *                   L2,M2,N2,B2,Bmat2,
c    *                   L3,M3,N3,B3,Bmat3,
c    *                   gamA12,gamA13,gamA23,
c    *                   gamB12,gamB13,gamB23,
c    *                   xgsg)
            call G3_MD_xggs(I1,J1,K1,A1,Amat1,
     *                      I2,J2,K2,A2,Amat2,
     *                      I3,J3,K3,A3,Amat3,
     *                      L1,M1,N1,B1,Bmat1,
     *                      L2,M2,N2,B2,Bmat2,
     *                      L3,M3,N3,B3,Bmat3,
     *                      gamA12,gamA13,gamA23,
     *                      gamB12,gamB13,gamB23,
     *                      xgsg)


C  Sum up 2 gamma loop terms
            gTPg =gTPg +bcoef1(ik)*bcoef1(iL)*xgTPg
            gTEg1=gTEg1+bcoef1(ik)*bcoef1(iL)*xgTEg1
            gTEg2=gTEg2+bcoef1(ik)*bcoef1(iL)*xgTEg2
            gTEg3=gTEg3+bcoef1(ik)*bcoef1(iL)*xgTEg3

            gVPCg =gVPCg +bcoef1(ik)*bcoef1(iL)*xgVPCg
            gVECg1=gVECg1+bcoef1(ik)*bcoef1(iL)*xgVECg1
            gVECg2=gVECg2+bcoef1(ik)*bcoef1(iL)*xgVECg2

            gVEPg1=gVEPg1+bcoef1(ik)*bcoef1(iL)*xgVEPg1
            gVEPg2=gVEPg2+bcoef1(ik)*bcoef1(iL)*xgVEPg2

            gVEEg1=gVEEg1+bcoef1(ik)*bcoef1(iL)*xgVEEg1
            gVEEg2=gVEEg2+bcoef1(ik)*bcoef1(iL)*xgVEEg2

            gsg=gsg+bcoef1(ik)*bcoef1(iL)*xgsg

C%%%%%%%%%%--End of 2 Gamma VEE Integrals--%%%%%%%%%%
      end do
         end do
C --- end of 2 gamma loop ---

C sum of terms: 
      ansE = two * (gTE + gVEC + gVEP)
     x     + four * gVEE * half
     x     + gTPg + gVPCg
     x     + gTEg1 + gTEG3
     x     + two * (gVEPg1 + gVECg1)
     x     + gTEg2 + gVEPg2 + gVECg2
     x     + two*gVEEg1*half
     x     + two*gVEEg2*half

      ansS=gsg

      return
      end

C======================================================================
      subroutine loop_size(L1,L2,nproc,myrank,istart,iend)

C======================================================================
      implicit none
      integer L1,L2,nproc,myrank,istart,iend
      integer iwork1,iwork2

      iwork1 = (L2 - L1 + 1) / nproc
      iwork2 = MOD(L2 - L1 + 1, nproc)
      istart = myrank * iwork1 + L1 + MIN(myrank, iwork2)
      iend   = istart + iwork1 - 1
      IF (iwork2 > myrank) iend = iend + 1


      return
      end

