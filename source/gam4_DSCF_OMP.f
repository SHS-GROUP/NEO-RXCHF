!======================================================================
      subroutine GAM4_DSCF(NG4CHK,ne,np,ngee,ng2,ng4,
     x                     DE,DP,focke,fockp,E_gam4)

!======================================================================
      implicit none
! Input Variables
      integer NG4CHK
      integer ng4,ng2,ngee,ne,np
      double precision E_gam4
      double precision DE(ne,ne)
      double precision DP(np,np)
      double precision focke(ne,ne)
      double precision fockp(np,np)
      
! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4
      integer,allocatable :: loop_map(:,:)
      double precision xfocke(ne,ne)
      double precision xfockp(np,np)
      double precision,allocatable :: GAM_ee(:)
      double precision,allocatable :: GAM_2s(:)
!     double precision,allocatable :: GAM_4(:)


!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_gam4=0.0d+00
      xfocke=0.0d+00
      xfockp=0.0d+00
!     call zero_out(nebf,xfocke)
!     call zero_out(npbf,xfockp)
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)


!----READ-GAM_ee-AND-GAM_2s-INTO-MEMORY-------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      allocate( GAM_ee(ngee),stat=istat )
!     write(*,*) 'allocate GAM_ee: ',istat

      if(allocated(GAM_2s)) deallocate(GAM_2s)
      allocate( GAM_2s(ng2),stat=istat )
!     write(*,*) 'allocate GAM_2s: ',istat

      call read_GAM_ee(ne,ngee,GAM_ee) 
      call read_GAM_2s(ne,np,ng2,GAM_2s) 
!----READ-GAM_ee-AND-GAM_2s-INTO-MEMORY-------------------------------)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(

      do ichunk=1,NG4CHK

         call loop_size(1,ng4,NG4CHK,ichunk-1,istart,iend)
!        write(*,*)'after call loop size'
!        write(*,*)'NG4=',ng4
!        write(*,*)'Nchunks=',NG4CHK
!        write(*,*)'ichunk=',ichunk
!        write(*,*)'istart=',istart
!        write(*,*)'iend=',iend

! Segment of ng4:
         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )
!        write(*,*) 'allocate loop_map: ',istat

!        if(allocated(GAM_4)) deallocate(GAM_4)
!        allocate( GAM_4(ng4_seg),stat=istat )
!        write(*,*) 'allocate GAM_4: ',istat

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

!        write(*,*)'done with nested loop compression'

!        call thread_gam4_DSCF(ne,np,ngee,ng2,ng4,ng4_seg,istart,iend,
!    x                 loop_map,GAM_ee,GAM_2s,focke,fockp,E_GAM4)
         call thread_gam4_DSCF(ne,np,ngee,ng2,ng4,ng4_seg,
     x                         istart,iend,loop_map,GAM_ee,GAM_2s,
     x                         DE,DP,xfocke,xfockp,E_GAM4)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      if(allocated(GAM_2s)) deallocate(GAM_2s)
      if(allocated(loop_map)) deallocate(loop_map)
!     if(allocated(GAM_4)) deallocate(GAM_4)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
      call add2fock(ne,xfocke,focke)
      call add2fock(np,xfockp,fockp)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!======================================================================
      subroutine thread_gam4_DSCF(ne,np,ngee,ng2,ng4,ng4_seg,
     x                            istart,iend,loop_map,GAM_ee,GAM_2s,
     x                            DE,DP,xfocke,xfockp,E_GAM4)
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
!     double precision GAM_4(ng4_seg) ! Array storage of 5-particle ints
      double precision xfocke(ne,ne)
      double precision xfockp(np,np)
      double precision DE(ne,ne)
      double precision DP(np,np)
      double precision E_GAM4
      double precision four
!     parameter(four=4.0d+00)

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

      four=4.0d+00

!---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

!        write(*,*)
!        write(*,*)'**************************************'
!        write(*,*)'    Computing GAM_4 Integrals    '
!        write(*,*)
!        write(*,*)'nebf      =',ne
!        write(*,*)'npbf      =',np
!        write(*,*)'Total ng4 =',ng4
!        write(*,*)'ISTART    =',istart
!        write(*,*)'IEND      =',iend
!        write(*,*)' Available processors: ',omp_get_num_procs()
!        write(*,*)' Available threads     ',omp_get_max_threads()
c        write(*,*)' Threads in use        ',omp_get_num_threads()
!        write(*,*)'**************************************'
!        write(*,*)

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(ne,np)
!$ompx shared(ngee)
!$ompx shared(DE)
!$ompx shared(DP)
!$ompx shared(ng2)
!$ompx shared(ng4_seg)
!$ompx shared(istart,iend)
!$ompx shared(four)
!$ompx shared(GAM_2s)
!$ompx shared(GAM_ee)
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
!$ompx reduction(+:E_GAM4)
!$ompx reduction(+:xfocke)
!$ompx reduction(+:xfockp)
!!!!!!!!!!!!!!mpx shared(GAM_4)
!!!mpx shared(xfocke)
!!!mpx shared(xfockp)

      id= omp_get_thread_num()
c     write(*,*)' Hello from process ',id
      if(id.eq.0) then
c        write(*,*)'Threads in use', omp_get_num_threads()
      end if

!$omp do SCHEDULE(RUNTIME)
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
!        GAM_4(imap)=ans
         ans=ans/(four*four)

         xfocke(ie1,je1)=xfocke(ie1,je1)+
     *       DP(ip,jp)*DE(ie2,je2)*DE(ie3,je3)*DE(ie4,je4)*four*ans

         xfockp(ip,jp)=xfockp(ip,jp)+
     *       DE(ie1,je1)*DE(ie2,je2)*DE(ie3,je3)*DE(ie4,je4)*ans

         E_gam4=E_gam4+
     x         DE(ie1,je1)*DE(ie2,je2)*
     x         DE(ie3,je3)*DE(ie4,je4)*DP(ip,jp)*ans


      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

!---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime() - wtime
!     write(*,*)'TIME TO CALCULATE GAM_4 INTEGRALS: ',wtime
!     write(*,*)
!---OPENMP-TIMING------------------------------------------------------)

     
      return
      end 
     
