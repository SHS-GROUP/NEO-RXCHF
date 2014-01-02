!======================================================================
      subroutine chop_GAM4(Nchunks,ne,np,ngee,ng2,ng4)

!======================================================================
      implicit none
! Input Variables
      integer Nchunks
      integer ng4,ng2,ngee,ne,np
      
! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4
      integer,allocatable :: loop_map(:,:)
      double precision,allocatable :: GAM_ee(:)
      double precision,allocatable :: GAM_2s(:)
      double precision,allocatable :: GAM_4(:)


!----READ-GAM_ee-AND-GAM_2s-INTO-MEMORY-------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      allocate( GAM_ee(ngee),stat=istat )
      write(*,*) 'allocate GAM_ee: ',istat

      if(allocated(GAM_2s)) deallocate(GAM_2s)
      allocate( GAM_2s(ng2),stat=istat )
      write(*,*) 'allocate GAM_2s: ',istat

      call read_GAM_ee(ne,ngee,GAM_ee) 
      call read_GAM_2s(ne,np,ng2,GAM_2s) 
!----READ-GAM_ee-AND-GAM_2s-INTO-MEMORY-------------------------------)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(

      do ichunk=1,Nchunks

         call loop_size(1,ng4,Nchunks,ichunk-1,istart,iend)
         write(*,*)'after call loop size'
         write(*,*)'NG4=',ng4
         write(*,*)'Nchunks=',Nchunks
         write(*,*)'ichunk=',ichunk
         write(*,*)'istart=',istart
         write(*,*)'iend=',iend

! Segment of ng4:
         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )
         write(*,*) 'allocate loop_map: ',istat

         if(allocated(GAM_4)) deallocate(GAM_4)
         allocate( GAM_4(ng4_seg),stat=istat )
         write(*,*) 'allocate GAM_4: ',istat

! Nested loop compression for this chunk:
! There should be a better way to do this
         Loopi=0
         imas=0
         do ip=1,np
         do jp=1,np
            do ie1=1,ne
            do je1=1,ne
               do ie2=1,ne
               do je2=1,ne
                  do ie3=1,ne
                  do je3=1,ne
                    do ie4=1,ne
                    do je4=1,ne

                         imas=imas+1 ! imas is master_index
                         if(imas.ge.istart.and.imas.le.iend) then
                            Loopi=Loopi+1
                            loop_map(Loopi,1)=je4
                            loop_map(Loopi,2)=ie4
                            loop_map(Loopi,3)=je3
                            loop_map(Loopi,4)=ie3
                            loop_map(Loopi,5)=je2
                            loop_map(Loopi,6)=ie2
                            loop_map(Loopi,7)=je1
                            loop_map(Loopi,8)=ie1
                            loop_map(Loopi,9)=jp
                            loop_map(Loopi,10)=ip
                         end if

                     end do
                     end do
                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

         write(*,*)'done with nested loop compression'

         call thread_gam4(ne,np,ngee,ng2,ng4,ng4_seg,istart,iend,
     x                    loop_map,GAM_ee,GAM_2s,GAM_4)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      if(allocated(GAM_2s)) deallocate(GAM_2s)
      if(allocated(loop_map)) deallocate(loop_map)
      if(allocated(GAM_4)) deallocate(GAM_4)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)


      return
      end
!======================================================================
!     subroutine GAM4_OMP_MD(ne,np,ng4,ng2,ngee)
      subroutine thread_gam4(ne,np,ngee,ng2,ng4,ng4_seg,istart,iend,
     x                       loop_map,GAM_ee,GAM_2s,GAM_4)
!
!======================================================================
      implicit none
      include 'omp_lib.h'
!     include 'mpif.h'

      integer istart,iend,imap
!     integer npe   ! Number of primitive electronic basis functions
      integer ne    ! Number of contracted electronic basis functions
      integer np    ! Number of nuclear basis functions
      integer ngee  ! Number of contracted 2-electron integrals
      integer ng2   ! Number of contracted 3-particle integrals
      integer ng4   ! Number of contracted 5-particle integrals
      integer ng4_seg ! dimension of chunk of contracted 5-particle integrals

      integer loop_map(ng4_seg,10)
      double precision GAM_ee(ngee)   ! Array storage of 2e-integrals
      double precision GAM_2s(ng2)    ! Array storage of 3-part overlaps
      double precision GAM_4(ng4_seg) ! Array storage of 5-particle ints

!--------------------------------(
! MPI-Related Local variables
!     integer :: myid,ierr,nprocs
!     integer :: loopi,iLP
!     integer :: istart,iend
!     integer :: loop_map(ng4,10)
!--------------------------------)
!---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

! Local variables
      integer ia
      integer ip
      integer jp
      integer ie1
      integer ie2
      integer ie3
      integer ie4
      integer je1
      integer je2
      integer je3
      integer je4
!CWS_int_stat
      integer iii
      integer jjj
      integer icount1
      integer icount2
      double precision threshold1
      double precision threshold2
!CWS_int_stat

      double precision ans

!---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

         write(*,*)
         write(*,*)'**************************************'
         write(*,*)'    Computing GAM_4 Integrals    '
         write(*,*)
         write(*,*)'nebf      =',ne
         write(*,*)'npbf      =',np
         write(*,*)'Total ng4 =',ng4
         write(*,*)'ISTART    =',istart
         write(*,*)'IEND      =',iend
         write(*,*)' Available processors: ',omp_get_num_procs()
         write(*,*)' Available threads     ',omp_get_max_threads()
         write(*,*)' Threads in use        ',omp_get_num_threads()
         write(*,*)'**************************************'
         write(*,*)

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(ne,np)
!$ompx shared(ngee)
!$ompx shared(ng2)
!$ompx shared(ng4_seg)
!$ompx shared(GAM_2s)
!$ompx shared(GAM_ee)
!$ompx shared(GAM_4)
!$ompx shared(IFIL)
!$ompx private(iLp) 
!$ompx private(imap) 
!$ompx private(ia)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(ie3,je3) 
!$ompx private(ie4,je4) 
!$ompx private(ans)
!$ompx private(id)

      id= omp_get_thread_num()
      write(*,*)' Hello from process ',id
      if(id.eq.0) then
         write(*,*)'Threads in use', omp_get_num_threads()
      end if

!$omp do
      do iLP=istart,iend

         imap=iLp-istart+1
         je4=loop_map(imap,1) 
         ie4=loop_map(imap,2) 
         je3=loop_map(imap,3) 
         ie3=loop_map(imap,4) 
         je2=loop_map(imap,5) 
         ie2=loop_map(imap,6) 
         je1=loop_map(imap,7) 
         ie1=loop_map(imap,8) 
         jp =loop_map(imap,9) 
         ip =loop_map(imap,10)

!        call index_GAM_4PK2(ne,np,
!    x                      ip,jp,
!    x                      ie1,je1,
!    x                      ie2,je2,
!    x                      ie3,je3,
!    x                      ie4,je4,ia)


         call symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

!        GAM_4(ia)=ans
         GAM_4(imap)=ans

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

!---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime() - wtime
      write(*,*)'TIME TO CALCULATE GAM_4 INTEGRALS: ',wtime
!---OPENMP-TIMING------------------------------------------------------)

!---------------WRITE-GAM4-INTEGRALS-TO-DISK---------------------------(
      open(806,file='GAM_4.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      icount1=0
      threshold1=1.0d-06
      icount2=0
      threshold2=1.0d-10

      do iLP=istart,iend
  
            imap=iLp-istart+1
            je4=loop_map(imap,1) 
            ie4=loop_map(imap,2) 
            je3=loop_map(imap,3) 
            ie3=loop_map(imap,4) 
            je2=loop_map(imap,5) 
            ie2=loop_map(imap,6) 
            je1=loop_map(imap,7) 
            ie1=loop_map(imap,8) 
            jp =loop_map(imap,9) 
            ip =loop_map(imap,10)

            call index_GAM_4PK2(ne,np,
     x                         ip,jp,
     x                         ie1,je1,
     x                         ie2,je2,
     x                         ie3,je3,
     x                         ie4,je4,ia)

!           write(806,REC=ia) GAM_4(ia)
            write(806,REC=ia) GAM_4(imap)

            ans=GAM_4(imap)
            if(abs(ans).lt.threshold1) icount1=icount1+1
            if(abs(ans).lt.threshold2) icount2=icount2+1

         end do  ! End iLP Loop

      close(806)
!CWS_int_stat
!     write(*,*)'GAM4 integrals skipped=',jjj*ne*ne*ne*ne*ne*ne
!     write(*,*)'ngam4 integrals greater than 1d-06=',iii
!     write(*,*)'number ng4 integrals eq 0       =',ng4-iii
      write(*,*)'GAM4 INTEGRALS LESS THAN 1.0D-06: ',icount1
      write(*,*)'GAM4 INTEGRALS LESS THAN 1.0D-10: ',icount2
!CWS_int_stat
!     END IF  ! endif for myid
!---------------WRITE-GAM4-INTEGRALS-TO-DISK---------------------------)

     
      return
      end 
     
C======================================================================
      subroutine GAM4_OMP_MD(ne,np,ng4,ng2,ngee)
C
C======================================================================
      implicit none
      include 'omp_lib.h'
c     include 'mpif.h'

      integer npe   ! Number of primitive electronic basis functions
      integer ne    ! Number of contracted electronic basis functions
      integer np    ! Number of nuclear basis functions
      integer ngee  ! Number of contracted 2-electron integrals
      integer ng2   ! Number of contracted 3-particle integrals
      integer ng4   ! Number of contracted 5-particle integrals

      double precision GAM_ee(ngee)  ! Array storage of 2e-integrals
      double precision GAM_2s(ng2)   ! Array storage of 3-part overlaps
c Array storage of 5-particle ints
c     double precision GAM_4(ng4)    ! Array storage of 5-particle ints
      double precision GAM_4(ng4)
c     double precision XGAM_4(ng4)   ! Buffer Array

C--------------------------------(
C MPI-Related Local variables
c     integer myid,ierr,nprocs
c     integer loopi,iLP
c     integer istart,iend
c     integer loop_map(ng4,10)
C--------------------------------)
C---OPENMP-RELATED-VARIABLES-----(
      integer IFIL
      integer id
      integer loopi,iLP
      double precision wtime
c     integer loop_map(ng4,10)
      integer loop_map(ng4,10)
C---OPENMP-RELATED-VARIABLES-----)

C Local variables
      integer ia
      integer ip
      integer jp
      integer ie1
      integer ie2
      integer ie3
      integer ie4
      integer je1
      integer je2
      integer je3
      integer je4
CCWS_int_stat
      integer iii
      integer jjj
      integer icount1
      integer icount2
      double precision threshold1
      double precision threshold2
CCWS_int_stat

      double precision ans
c     double precision tolerance
c     logical prin

c     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
c     call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)

C---OPENMP-TIMING------------------------------------------------------(
C      wtime = omp_get_wtime()
C---OPENMP-TIMING------------------------------------------------------)

C-----------INITIALIZE-DATA-STRUCTURES-ON-MASTER-----------------------(
c     IF(myid .eq. 0) THEN

         write(*,*)
         write(*,*)'**************************************'
         write(*,*)'    Computing GAM_4 Integrals    '
         write(*,*)
         write(*,*)'nebf =',ne
         write(*,*)'npbf =',np
         write(*,*)'ng4  =',ng4
         write(*,*)' Available processors: ',omp_get_num_procs()
         write(*,*)' Available threads     ',omp_get_max_threads()
         write(*,*)' Threads in use        ',omp_get_num_threads()
         write(*,*)'**************************************'
         write(*,*)

C Nested loop compression:
         Loopi=0
         do ip=1,np
         do jp=1,np
            do ie1=1,ne
            do je1=1,ne
               do ie2=1,ne
               do je2=1,ne
                  do ie3=1,ne
                  do je3=1,ne
                    do ie4=1,ne
                    do je4=1,ne

                           Loopi=Loopi+1
                           loop_map(Loopi,1)=je4
                           loop_map(Loopi,2)=ie4
                           loop_map(Loopi,3)=je3
                           loop_map(Loopi,4)=ie3
                           loop_map(Loopi,5)=je2
                           loop_map(Loopi,6)=ie2
                           loop_map(Loopi,7)=je1
                           loop_map(Loopi,8)=ie1
                           loop_map(Loopi,9)=jp
                           loop_map(Loopi,10)=ip

c                          call index_GAM_4PK(ne,np,
c    x                                        ip,jp,
c    x                                        ie1,je1,
c    x                                        ie2,je2,
c    x                                        ie3,je3,
c    x                                        ie4,je4,ia)

c                          GAM_4(ia)=0.0d+00

                     end do
                     end do
                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do

         write(*,*)'done with nested loop compression'
C Read the 2-electron integrals from disk
         call read_GAM_ee(ne,ngee,GAM_ee) 

 
C Read the gam_2s integrals from disk
         call read_GAM_2s(ne,np,ng2,GAM_2s) 
c        write(*,*)'GAM_ee is:',GAM_ee
c        write(*,*)'GAM_2s is:',GAM_2s
 
CCWS_int_stat
c     tolerance=1.0d-02
c     iii=0
c     jjj=0
CCWS_int_stat
c     END IF  ! endif myid
C-----------INITIALIZE-DATA-STRUCTURES-ON-MASTER-----------------------)

C-------MPI-BROADCAST-DATA---------------------------------------------(
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(loop_map,ng4*10,mpi_integer,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(GAM_4,ng4,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(GAM_ee,ngee,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BCAST(GAM_2s,ng2,mpi_double_precision,0,
c    x     MPI_COMM_WORLD,ierr)

c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C-------MPI-BROADCAST-DATA---------------------------------------------)

C-----------CALCULATION-OF-GAM4-INTEGRALS-ON-ALL-PROCESSORS------------(
c     call loop_size(1,ng4,nprocs,myid,istart,iend)
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
c        IFIL=806
c        open(IFIL,file='GAM_4.ufm',form='unformatted',
c    x    status='unknown',access='direct',RECL=8)
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(ne,np)
!$ompx shared(ngee)
!$ompx shared(ng2)
!$ompx shared(ng4)
!$ompx shared(GAM_2s)
!$ompx shared(GAM_ee)
!$ompx shared(GAM_4)
!$ompx shared(IFIL)
!$ompx private(iLp) 
!$ompx private(ia)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(ie2,je2) 
!$ompx private(ie3,je3) 
!$ompx private(ie4,je4) 
!$ompx private(ans)
!$ompx private(id)
CCCCCCCCCCCCCCCCCCCCCCmpx reduction(+:GAM_4)

      id= omp_get_thread_num()
      write(*,*)' Hello from process ',id
      if(id.eq.0) then
         write(*,*)'Threads in use', omp_get_num_threads()
      end if

!$omp do
c     do iLP=istart,iend
c     do iLP=1,ng4
      do iLP=1,ng4
c        do ip=1,np
c        do jp=1,np
c           do ie1=1,ne
c           do je1=1,ne
c              do ie2=1,ne
c              do je2=1,ne
c                 do ie3=1,ne
c                 do je3=1,ne
c                   do ie4=1,ne
c                   do je4=1,ne
         je4=loop_map(iLP,1) 
         ie4=loop_map(iLP,2) 
         je3=loop_map(iLP,3) 
         ie3=loop_map(iLP,4) 
         je2=loop_map(iLP,5) 
         ie2=loop_map(iLP,6) 
         je1=loop_map(iLP,7) 
         ie1=loop_map(iLP,8) 
         jp =loop_map(iLP,9) 
         ip =loop_map(iLP,10)

         call index_GAM_4PK2(ne,np,
     x                      ip,jp,
     x                      ie1,je1,
     x                      ie2,je2,
     x                      ie3,je3,
     x                      ie4,je4,ia)


         call symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

c        gam_4(ie1,je1,ie2,je2,ie3,je3,ie4,je4)=ans
c        write(*,*)'IA=  ANS=',ia,ans
         GAM_4(ia)=ans
c        write(IFIL,REC=ia) ans

c                    end do
c                    end do
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
c        close(IFIL)
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

C---OPENMP-TIMING------------------------------------------------------(
C      wtime = omp_get_wtime() - wtime
C      write(*,*)'TIME TO CALCULATE GAM_4 INTEGRALS: ',wtime
C---OPENMP-TIMING------------------------------------------------------)


C-----------CALCULATION-OF-GAM4-INTEGRALS-ON-ALL-PROCESSORS------------)

c     if(myid.eq.0) then
c        write(*,*)'MYID=',myid
c        write(*,*)'GAM4 is:',GAM_4
c     end if
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     if(myid.eq.1) then
c        write(*,*)'MYID=',myid
c        write(*,*)'GAM4 is:',GAM_4
c     end if

C-----------GATHER-GAM_4-DATA-FROM-ALL-PROCESSORS----------------------(
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     call MPI_REDUCE(GAM_4,XGAM_4,ng4,
c    x                mpi_double_precision,mpi_sum,0,
c    x                MPI_COMM_WORLD,ierr)
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C-----------GATHER-GAM_4-DATA-FROM-ALL-PROCESSORS----------------------)

C---------------WRITE-GAM4-INTEGRALS-TO-DISK---------------------------(
c     IF(myid .eq. 0) THEN

c        GAM_4=XGAM_4

c        write(*,*)'AFTER MPI REDUCE: MYID=',myid
c        write(*,*)'GAM4 is now:',GAM_4

c     goto 800
         open(806,file='GAM_4.ufm',form='unformatted',
     x    status='unknown',access='direct',RECL=8)

         icount1=0
         threshold1=1.0d-06
         icount2=0
         threshold2=1.0d-10
c        do iLP=1,ng4
         do ip=1,np
         do jp=1,np
            do ie1=1,ne
            do je1=1,ne
               do ie2=1,ne
               do je2=1,ne
                  do ie3=1,ne
                  do je3=1,ne
                    do ie4=1,ne
                    do je4=1,ne
  
c           je4=loop_map(iLP,1) 
c           ie4=loop_map(iLP,2) 
c           je3=loop_map(iLP,3) 
c           ie3=loop_map(iLP,4) 
c           je2=loop_map(iLP,5) 
c           ie2=loop_map(iLP,6) 
c           je1=loop_map(iLP,7) 
c           ie1=loop_map(iLP,8) 
c           jp =loop_map(iLP,9) 
c           ip =loop_map(iLP,10)

            call index_GAM_4PK2(ne,np,
     x                         ip,jp,
     x                         ie1,je1,
     x                         ie2,je2,
     x                         ie3,je3,
     x                         ie4,je4,ia)

            ans = GAM_4(ia)

            write(806,REC=ia) ans

            if(ans.lt.threshold1) icount1=icount1+1
            if(ans.lt.threshold2) icount2=icount2+1

                     end do
                     end do
                  end do
                  end do
               end do
               end do
            end do
            end do
         end do
         end do
c        end do  ! End iLP Loop

         close(806)
c 800 continue
CCWS_int_stat
c     write(*,*)'GAM4 integrals skipped=',jjj*ne*ne*ne*ne*ne*ne
c     write(*,*)'ngam4 integrals greater than 1d-06=',iii
c     write(*,*)'number ng4 integrals eq 0       =',ng4-iii
      write(*,*)'GAM4 INTEGRALS LESS THAN 1.0D-06: ',icount1
      write(*,*)'GAM4 INTEGRALS LESS THAN 1.0D-10: ',icount2
CCWS_int_stat
c     END IF  ! endif for myid
C---------------WRITE-GAM4-INTEGRALS-TO-DISK---------------------------)


CDEBUG
c        write(*,*)'>>>> GAM_4 Integral Debug <<<<'
c        ip = 1
c        jp = 1
c        ie1= 1
c        je2= 1
c        ie2= 1
c        je2= 1
c        ie3= 2
c        je3= 2
c        ie4= 2
c        je4= 2

c     call index_GAM_4PK(ne,np,
c    x                   ip,jp,
c    x                   ie1,je1,
c    x                   ie2,je2,
c    x                   ie3,je3,
c    x                   ie4,je4,ia)

c        call symm_gam4(ne,np,GAM_2s,GAM_ee,
c    x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

c        write(*,*)'ia=',ia
c        write(*,*)'ans=',ans
c        write(*,*)
c        write(*,*)'>>>> end of GAM_4 Integral Debug <<<<'
CDEBUG
     
      return
      end 
     
C======================================================================
      SUBROUTINE symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     1                     ip,jp,ii,i,jj,j,kk,k,ll,l,ans)
c     SUBROUTINE symm_gam4(ne,np,sgm2,Vee,
c    1                     ii,i,jj,j,kk,k,ll,l,ans)
c     SUBROUTINE symm_omg4_3(NE,NG,NVEE,sgm2,Vee,
c    1                     ii,i,jj,j,kk,k,ll,l,ans)
C
C======================================================================
C     ii  i   jj  j   kk  k   ll  l   ip  jp
C     ie1,je1,ie2,je2,ie3,je3,ie4,je4,ip1,jp1,ans)
C 
      implicit none
c     include 'param.h'
      double precision two,four
      parameter(two=2.0d+00,four=4.0d+00)
C input
      integer ne
      integer np
      integer ng2
      integer ngee
c     double precision GAM_2s(np,np,ne,ne,ne,ne)
c     double precision GAM_ee(ne,ne,ne,ne)
      double precision GAM_ee(ngee)
CCWS-IO
      double precision GAM_2s(ng2)
c     double precision GAM_2s(1)
      integer ip
      integer jp
      integer ii,i
      integer jj,j
      integer kk,k
      integer ll,l

c     logical prin
C output
      double precision ans
C local
      double precision x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
      double precision x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24
      double precision  EIGHT


      EIGHT = FOUR * TWO

      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                         ip,jp,ii,jj,kk,ll,i,j,k,l,x1)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                         ip,jp,ii,jj,kk,ll,j,i,k,l,x2)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                         ip,jp,ii,jj,kk,ll,k,i,j,l,x3)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                         ip,jp,ii,jj,kk,ll,i,k,j,l,x4)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                         ip,jp,ii,jj,kk,ll,j,k,i,l,x5)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                         ip,jp,ii,jj,kk,ll,k,j,i,l,x6)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                         ip,jp,ii,jj,kk,ll,l,j,i,k,x7)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                         ip,jp,ii,jj,kk,ll,j,l,i,k,x8)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                         ip,jp,ii,jj,kk,ll,i,l,j,k,x9)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,l,i,j,k,x10)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,j,i,l,k,x11)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,i,j,l,k,x12)

      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,i,k,l,j,x13)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,k,i,l,j,x14)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,l,i,k,j,x15)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,i,l,k,j,x16)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,k,l,i,j,x17)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,l,k,i,j,x18)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,l,k,j,i,x19)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,k,l,j,i,x20)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,j,l,k,i,x21)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,l,j,k,i,x22)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,k,j,l,i,x23)
      call i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                        ip,jp,ii,jj,kk,ll,j,k,l,i,x24)

c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,i,j,k,l,x1)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,j,i,k,l,x2)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,k,i,j,l,x3)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,i,k,j,l,x4)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,j,k,i,l,x5)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,k,j,i,l,x6)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,l,j,i,k,x7)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,j,l,i,k,x8)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,i,l,j,k,x9)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,l,i,j,k,x10)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,j,i,l,k,x11)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,i,j,l,k,x12)

c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,i,k,l,j,x13)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,k,i,l,j,x14)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,l,i,k,j,x15)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,i,l,k,j,x16)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,k,l,i,j,x17)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,l,k,i,j,x18)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,l,k,j,i,x19)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,k,l,j,i,x20)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,j,l,k,i,x21)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,l,j,k,i,x22)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,k,j,l,i,x23)
c     call eazy10_3(NE,np,GAM_2s,GAM_ee,ip,jp,ii,jj,kk,ll,j,k,l,i,x24)

!-----Before modification 09-09-2011(
      ans = (FOUR*FOUR*x1)
     1    - (EIGHT    *x2)
     2    + (FOUR     *x3)
     3    - (EIGHT    *x4)
     4    + (FOUR     *x5)
     5    - (EIGHT    *x6)
     6    + (FOUR     *x7)
     7    - (TWO      *x8)
     8    + (FOUR     *x9)
     9    - (TWO      *x10)
     1    + (FOUR     *x11)
     2    - (EIGHT    *x12)
     3    + (FOUR     *x13)
     4    - (TWO      *x14)
     5    + (FOUR     *x15)
     6    - (EIGHT    *x16)
     7    + (FOUR     *x17)
     8    - (TWO      *x18)
     9    + (FOUR     *x19)
     1    - (TWO      *x20)
     2    + (FOUR     *x21)
     3    - (EIGHT    *x22)
     4    + (FOUR     *x23)
     5    - (TWO      *x24)
!-----Before modification 09-09-2011)
!-----After modification 09-09-2011(
!     ans = (FOUR*FOUR*x1)
!    x    - (EIGHT    *x2)
!    x    + (FOUR     *x3)
!    x    - (EIGHT    *x4)
!    x    + (FOUR     *x5)
!    x    - (EIGHT    *x6)
!    x    + (FOUR     *x7)
!    x    - (TWO      *x8)
!    x    + (FOUR     *x9)
!    x    - (TWO      *x10)
!    x    - (FOUR     *x11)
!    x    - (EIGHT    *x12)
!    x    + (FOUR     *x13)
!    x    - (TWO      *x14)
!    x    + (FOUR     *x15)
!    x    - (EIGHT    *x16)
!    x    - (FOUR     *x17)
!    x    - (TWO      *x18)
!    x    - (FOUR     *x19)
!    x    - (TWO      *x20)
!    x    + (FOUR     *x21)
!    x    - (EIGHT    *x22)
!    x    + (FOUR     *x23)
!    x    - (TWO      *x24)
!-----After modification 09-09-2011)

CDEBUG
c        write(*,*)'In symm_gam4 routine:'
c        write(*,*)
c       write(*,*)'x1=',x1
c       write(*,*)'x2=',x2
c       write(*,*)'x3=',x3
c       write(*,*)'x4=',x4
c       write(*,*)'x5=',x5
c       write(*,*)'x6=',x6
c       write(*,*)'x7=',x7
c       write(*,*)'x8=',x8
c       write(*,*)'x9=',x9
c       write(*,*)'x10=',x10
c       write(*,*)'x11=',x11
c       write(*,*)'x12=',x12
c       write(*,*)'x13=',x13
c       write(*,*)'x14=',x14
c       write(*,*)'x15=',x15
c       write(*,*)'x16=',x16
c       write(*,*)'x17=',x17
c       write(*,*)'x18=',x18
c       write(*,*)'x19=',x19
c       write(*,*)'x20=',x20
c       write(*,*)'x21=',x21
c       write(*,*)'x22=',x22
c       write(*,*)'x23=',x23
c       write(*,*)'x24=',x24
c        write(*,*)'ans=',ans

c        write(*,*)'>>>> end of GAM_4 Integral Debug <<<<'
CDEBUG

      return
      END

C======================================================================
      SUBROUTINE i10(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     1             ip,jp,ie1,ie2,ie3,ie4,je1,je2,je3,je4,ans)
c     SUBROUTINE eazy10_3(NE,NG2,NVEE,sgmat2,Vee,
c    1             ie1,ie2,ie3,ie4,je1,je2,je3,je4,ans)
C======================================================================
C Physicist to chemist notation conversion
      implicit none
      integer ne
      integer np
      integer ng2
      integer ngee
c     double precision GAM_2s(np,np,ne,ne,ne,ne)
c     double precision GAM_ee(ne,ne,ne,ne)
CCWS-IO
      double precision GAM_2s(ng2)
c     double precision GAM_2s(1)
      double precision GAM_ee(ngee)

c     logical prin

      integer ip,jp
      integer ie1,je1
      integer ie2,je2
      integer ie3,je3
      integer ie4,je4
      double precision ans


c     call calc_GAM_4_integral(NE,NG2,NVEE,ie1,je1,ie2,je2,
c    1                      ie3,je3,ie4,je4,sgmat2,Vee,ans)

      call calc_GAM_4_integral(ne,np,ng2,ngee,
     x                         GAM_2s,GAM_ee,ip,jp,
     x                         ie1,je1,ie2,je2,
     x                         ie3,je3,ie4,je4,
     x                         ans)
C

      return
      END


C======================================================================
      subroutine calc_GAM_4_integral(ne,np,ng2,ngee,
     x                               GAM_2s,GAM_ee,ip,jp,
     x                               ie1,je1,ie2,je2,
     x                               ie3,je3,ie4,je4,
     x                               ans)
C
C======================================================================
      implicit none

      double precision half
      double precision six
      parameter(half=0.5d+00,six=6.0d+00)

      integer ne
      integer np
      integer ng2
      integer ngee
      integer ip
      integer jp
      integer ie1
      integer ie2
      integer ie3
      integer ie4
      integer je1
      integer je2
      integer je3
      integer je4

c     double precision GAM_ee(ne,ne,ne,ne)
c     double precision GAM_2s(np,np,ne,ne,ne,ne)
CCWS-IO
      double precision GAM_2s(ng2)
c     double precision GAM_2s(1)
      double precision GAM_ee(ngee)
      double precision ans

c     logical prin

C Local variables
      double precision x1,x2,x3,x4,x5,x6

      integer iee34 
      integer iee24
      integer iee23
      integer iee14
      integer iee13
      integer iee12

      integer igs12
      integer igs13
      integer igs14
      integer igs23
      integer igs24
      integer igs34

      double precision ee34
      double precision ee24
      double precision ee23
      double precision ee14
      double precision ee13
      double precision ee12

      double precision gs12
      double precision gs13
      double precision gs14
      double precision gs23
      double precision gs24
      double precision gs34

c     open(804,file='GAM_2s.ufm',form='unformatted',
c    x status='unknown',access='direct',RECL=8)

c     x1=half*gam_ee(ie3,je3,ie4,je4)*GAM_2s(ie1,je1,ie2,je2)
c     x2=half*gam_ee(ie2,je2,ie4,je4)*GAM_2s(ie1,je1,ie3,je3)
c     x3=half*gam_ee(ie2,je2,ie3,je3)*GAM_2s(ie1,je1,ie4,je4)
c     x4=half*gam_ee(ie1,je1,ie4,je4)*GAM_2s(ie2,je2,ie3,je3)
c     x5=half*gam_ee(ie1,je1,ie3,je3)*GAM_2s(ie2,je2,ie4,je4)
c     x6=half*gam_ee(ie1,je1,ie2,je2)*GAM_2s(ie3,je3,ie4,je4)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C GAM_1(je1,ie1,jp,ip)    |  ==>  array(i1,i2,i3,i4)
C GAM_ee(je2,ie2,je1,ie1) |

cc    ee34=gam_ee(ie3,je3,ie4,je4) 
cc    ee24=gam_ee(ie2,je2,ie4,je4) 
cc    ee23=gam_ee(ie2,je2,ie3,je3) 
cc    ee14=gam_ee(ie1,je1,ie4,je4) 
cc    ee13=gam_ee(ie1,je1,ie3,je3) 
cc    ee12=gam_ee(ie1,je1,ie2,je2) 

cc    call pack_4D(ne,ne,ne,ie3,je3,ie4,je4,iee34)
cc    call pack_4D(ne,ne,ne,ie2,je2,ie4,je4,iee24)
cc    call pack_4D(ne,ne,ne,ie2,je2,ie3,je3,iee23)
cc    call pack_4D(ne,ne,ne,ie1,je1,ie4,je4,iee14)
cc    call pack_4D(ne,ne,ne,ie1,je1,ie3,je3,iee13)
cc    call pack_4D(ne,ne,ne,ie1,je1,ie2,je2,iee12)

      call pack_4D(ne,ne,ne,je4,ie4,je3,ie3,iee34)
      call pack_4D(ne,ne,ne,je4,ie4,je2,ie2,iee24)
      call pack_4D(ne,ne,ne,je3,ie3,je2,ie2,iee23)
      call pack_4D(ne,ne,ne,je4,ie4,je1,ie1,iee14)
      call pack_4D(ne,ne,ne,je3,ie3,je1,ie1,iee13)
      call pack_4D(ne,ne,ne,je2,ie2,je1,ie1,iee12)


      ee34=gam_ee(iee34) 
      ee24=gam_ee(iee24) 
      ee23=gam_ee(iee23) 
      ee14=gam_ee(iee14) 
      ee13=gam_ee(iee13) 
      ee12=gam_ee(iee12) 

      call underflow(ee34)
      call underflow(ee24)
      call underflow(ee23)
      call underflow(ee14)
      call underflow(ee13)
      call underflow(ee12)

cc    call underflow(gam_ee(ie3,je3,ie4,je4))
cc    call underflow(gam_ee(ie2,je2,ie4,je4))
cc    call underflow(gam_ee(ie2,je2,ie3,je3))
cc    call underflow(gam_ee(ie1,je1,ie4,je4))
cc    call underflow(gam_ee(ie1,je1,ie3,je3))
cc    call underflow(gam_ee(ie1,je1,ie2,je2))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

cc    gs12=GAM_2s(ip,jp,ie1,je1,ie2,je2) 
cc    gs13=GAM_2s(ip,jp,ie1,je1,ie3,je3) 
cc    gs14=GAM_2s(ip,jp,ie1,je1,ie4,je4) 
cc    gs23=GAM_2s(ip,jp,ie2,je2,ie3,je3) 
cc    gs24=GAM_2s(ip,jp,ie2,je2,ie4,je4) 
cc    gs34=GAM_2s(ip,jp,ie3,je3,ie4,je4) 


      call index_GAM_2PK(ne,np,ip,jp,ie1,je1,ie2,je2,igs12)
      call index_GAM_2PK(ne,np,ip,jp,ie1,je1,ie3,je3,igs13)
      call index_GAM_2PK(ne,np,ip,jp,ie1,je1,ie4,je4,igs14)
      call index_GAM_2PK(ne,np,ip,jp,ie2,je2,ie3,je3,igs23)
      call index_GAM_2PK(ne,np,ip,jp,ie2,je2,ie4,je4,igs24)
      call index_GAM_2PK(ne,np,ip,jp,ie3,je3,ie4,je4,igs34)

c     read(804,REC=igs12) gs12
c     read(804,REC=igs13) gs13
c     read(804,REC=igs14) gs14
c     read(804,REC=igs23) gs23
c     read(804,REC=igs24) gs24
c     read(804,REC=igs34) gs34

      gs12=GAM_2s(igs12) 
      gs13=GAM_2s(igs13) 
      gs14=GAM_2s(igs14) 
      gs23=GAM_2s(igs23) 
      gs24=GAM_2s(igs24) 
      gs34=GAM_2s(igs34) 

      call underflow(gs12)
      call underflow(gs13)
      call underflow(gs14)
      call underflow(gs23)
      call underflow(gs24)
      call underflow(gs34)

c     call underflow(GAM_2s(ip,jp,ie1,je1,ie2,je2))
c     call underflow(GAM_2s(ip,jp,ie1,je1,ie3,je3))
c     call underflow(GAM_2s(ip,jp,ie1,je1,ie4,je4))
c     call underflow(GAM_2s(ip,jp,ie2,je2,ie3,je3))
c     call underflow(GAM_2s(ip,jp,ie2,je2,ie4,je4))
c     call underflow(GAM_2s(ip,jp,ie3,je3,ie4,je4))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     if(prin) then
c     write(*,*)'xxx_gam_ee=',gam_ee(ie3,je3,ie4,je4) 
c     write(*,*)'xxx_gam_ee=',gam_ee(ie2,je2,ie4,je4) 
c     write(*,*)'xxx_gam_ee=',gam_ee(ie2,je2,ie3,je3) 
c     write(*,*)'xxx_gam_ee=',gam_ee(ie1,je1,ie4,je4) 
c     write(*,*)'xxx_gam_ee=',gam_ee(ie1,je1,ie3,je3) 
c     write(*,*)'xxx_gam_ee=',gam_ee(ie1,je1,ie2,je2) 

c     write(*,*)'yyy_GAM_2s=',GAM_2s(ip,jp,ie1,je1,ie2,je2) 
c     write(*,*)'yyy_GAM_2s=',GAM_2s(ip,jp,ie1,je1,ie3,je3) 
c     write(*,*)'yyy_GAM_2s=',GAM_2s(ip,jp,ie1,je1,ie4,je4) 
c     write(*,*)'yyy_GAM_2s=',GAM_2s(ip,jp,ie2,je2,ie3,je3) 
c     write(*,*)'yyy_GAM_2s=',GAM_2s(ip,jp,ie2,je2,ie4,je4) 
c     write(*,*)'yyy_GAM_2s=',GAM_2s(ip,jp,ie3,je3,ie4,je4) 

c        write(*,*)'In symm_gam4 routine:'
c        write(*,*)
c     write(*,*)'==========='
c     write(*,*)'ee34=',ee34 
c     write(*,*)'ee24=',ee24 
c     write(*,*)'ee23=',ee23 
c     write(*,*)'ee14=',ee14 
c     write(*,*)'ee13=',ee13 
c     write(*,*)'ee12=',ee12 
c     write(*,*)
c     write(*,*)'gs12=',gs12 
c     write(*,*)'gs13=',gs13 
c     write(*,*)'gs14=',gs14 
c     write(*,*)'gs23=',gs23 
c     write(*,*)'gs24=',gs24 
c     write(*,*)'gs34=',gs34 
c     write(*,*)'==========='
c     end if
c     call flshbf(6)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     x1=half*gam_ee(ie3,je3,ie4,je4)*GAM_2s(ip,jp,ie1,je1,ie2,je2)
c     x2=half*gam_ee(ie2,je2,ie4,je4)*GAM_2s(ip,jp,ie1,je1,ie3,je3)
c     x3=half*gam_ee(ie2,je2,ie3,je3)*GAM_2s(ip,jp,ie1,je1,ie4,je4)
c     x4=half*gam_ee(ie1,je1,ie4,je4)*GAM_2s(ip,jp,ie2,je2,ie3,je3)
c     x5=half*gam_ee(ie1,je1,ie3,je3)*GAM_2s(ip,jp,ie2,je2,ie4,je4)
c     x6=half*gam_ee(ie1,je1,ie2,je2)*GAM_2s(ip,jp,ie3,je3,ie4,je4)

      x1=half*ee34*gs12
      x2=half*ee24*gs13
      x3=half*ee23*gs14
      x4=half*ee14*gs23
      x5=half*ee13*gs24
      x6=half*ee12*gs34

      ans=(x1+x2+x3+x4+x5+x6)/six

c     close(804)

     
      return
      end 

C=======================================================================
      subroutine read_GAM_2s(ne,np,ng2,GAM_2s)
C=======================================================================
      implicit none
C Input Variables
      integer ne,np,ng2
C Variables Returned
      double precision GAM_2s(ng2)
C Local Variables
      integer ia
      integer ip,jp
      integer ie1,je1
      integer ie2,je2
      double precision ans

      open(804,file='GAM_2s.ufm',form='unformatted',
     x status='unknown',access='direct',RECL=8)

      do ip=1,np
       do jp=1,np
        do ie1=1,ne
         do je1=1,ne
          do ie2=1,ne
           do je2=1,ne

              call index_GAM_2PK(ne,np,ip,jp,ie1,je1,ie2,je2,ia)
              read(804,REC=ia) ans
              GAM_2s(ia)=ans

           end do
          end do
         end do
        end do
       end do
      end do

      close(804)


      return
      end
     
