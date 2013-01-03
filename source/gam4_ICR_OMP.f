!======================================================================
      subroutine GAM4_ICR(Nchunks,ne,np,ngee,ng2,ng4,
     x                    GAM_2s,GAM_4)

!======================================================================
      implicit none
      include 'omp_lib.h'
! Input Variables
      integer Nchunks
      integer ng4,ng2,ngee,ne,np
      double precision GAM_4(ng4)
      double precision GAM_2s(ng2)
      
! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4
      integer,allocatable :: loop_map(:,:)
      double precision,allocatable :: GAM_ee(:)
!     double precision,allocatable :: GAM_2s(:)
      double precision wtime,wtime2


         write(*,1000) ng4,nchunks,omp_get_num_procs(),
     xomp_get_max_threads(),1
         wtime = omp_get_wtime()

!----READ-GAM_ee-AND-GAM_2s-INTO-MEMORY-------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      allocate( GAM_ee(ngee),stat=istat )
!     write(*,*) 'allocate GAM_ee: ',istat
      call read_GAM_ee(ne,ngee,GAM_ee) 

!     if(.NOT.LG2IC1) then
!        if(allocated(GAM_2s)) deallocate(GAM_2s)
!        allocate( GAM_2s(ng2),stat=istat )
!        write(*,*) 'allocate GAM_2s: ',istat
!        call read_GAM_2s(ne,np,ng2,GAM_2s) 
!     end if
!----READ-GAM_ee-AND-GAM_2s-INTO-MEMORY-------------------------------)

!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------(

      do ichunk=1,Nchunks

         wtime2 = omp_get_wtime()

         call loop_size(1,ng4,Nchunks,ichunk-1,istart,iend)

         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )
!        write(*,*) 'allocate loop_map: ',istat

! Nested loop compression for this chunk:
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

         call thread_gam4_IC(ne,np,ngee,ng2,ng4,ng4_seg,istart,iend,
     x                       loop_map,GAM_ee,GAM_2s,GAM_4)

         wtime2 = omp_get_wtime() - wtime2
         write(*,2000)ichunk,wtime2

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(GAM_ee)) deallocate(GAM_ee)
!     if(allocated(GAM_2s)) deallocate(GAM_2s)
      if(allocated(loop_map)) deallocate(loop_map)
!     if(allocated(GAM_4)) deallocate(GAM_4)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

      wtime = omp_get_wtime() - wtime
      write(*,3000)wtime


 1000 FORMAT(/6X,'+---------------------------------------------+',/,
     x        6X,'|     CALCULATING 5-PARTICLE INTEGRALS        |',/,
     x        6X,'|            --IN-CORE APPROACH--             |',/,
     x        6X,'+---------------------------------------------+',/,
     x        8X,'                          ',/,
     x        8X,'   NUMBER OF 5-PARTICLE INTEGRALS: ',1X,I12/
     x        8X,'  NUMBER OF BLOCKS (USER DEFINED): ',1X,I12/
     x        8X,'                          ',/,
     x        8X,'  COMPUTATIONAL RESOURCES:',/,
     x        8X,'  ------------------------',/,
     x        8X,'     CORES PER NODE:',1X,I3/
     x        8X,'          AVAILABLE:',1X,I3/
     x        8X,'    NUMBER OF NODES:',1X,I3/)
                       
 2000 FORMAT(8X,'    TIME TO EVALUATE BLOCK ',1X,I4,1X,F10.2)

 3000 FORMAT(/8X,'  TIMING SUMMARY FOR 5-PARTICLE INTEGRALS:',/,
     x        8X,'  ----------------------------------------',/,
     x        8X,'    TIME TO EVALUATE ALL INTEGRALS:',1X,F12.4)


      return
      end
!======================================================================
      subroutine thread_gam4_IC(ne,np,ngee,ng2,ng4,ng4_seg,istart,iend,
     x                          loop_map,GAM_ee,GAM_2s,GAM_4)
!
!======================================================================
      implicit none
      include 'omp_lib.h'

      integer istart,iend,imap
      integer ne    ! Number of contracted electronic basis functions
      integer np    ! Number of nuclear basis functions
      integer ngee  ! Number of contracted 2-electron integrals
      integer ng2   ! Number of contracted 3-particle integrals
      integer ng4   ! Number of contracted 5-particle integrals
      integer ng4_seg ! dimension of chunk of contracted 5-particle integrals

      integer loop_map(ng4_seg,10)
      double precision GAM_ee(ngee)   ! Array storage of 2e-integrals
      double precision GAM_2s(ng2)    ! Array storage of 3-part overlaps
      double precision GAM_4(ng4)     ! Array storage of 5-particle ints

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
!     integer iii
!     integer jjj
!     integer icount1
!     integer icount2
!     double precision threshold1
!     double precision threshold2
!CWS_int_stat

      double precision four
!     parameter(four=4.0d+00)
      double precision ans

      four=4.0d+00

!---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(ne,np)
!$ompx shared(ngee)
!$ompx shared(ng2)
!$ompx shared(ng4_seg)
!$ompx shared(four)
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

         call index_GAM_4PK2(ne,np,
     x                       ip,jp,
     x                       ie1,je1,
     x                       ie2,je2,
     x                       ie3,je3,
     x                       ie4,je4,ia)


         call symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                  ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4,ans)

         ans=ans/(four*four)
         GAM_4(ia)=ans

      end do
!$omp end do
!$omp end parallel      
!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

!     icount1=0
!     threshold1=1.0d-06
!     icount2=0
!     threshold2=1.0d-10
!     ans=GAM_4(imap)
!     if(abs(ans).lt.threshold1) icount1=icount1+1
!     if(abs(ans).lt.threshold2) icount2=icount2+1
!     write(*,*)'GAM4 integrals skipped=',jjj*ne*ne*ne*ne*ne*ne
!     write(*,*)'ngam4 integrals greater than 1d-06=',iii
!     write(*,*)'number ng4 integrals eq 0       =',ng4-iii
!     write(*,*)'GAM4 INTEGRALS LESS THAN 1.0D-06: ',icount1
!     write(*,*)'GAM4 INTEGRALS LESS THAN 1.0D-10: ',icount2

     
      return
      end 
!=======================================================================
      subroutine EG4ICR(Nchunks,nebf,npbf,ng4,
     x                  DE,DP,GM4ICR,focke,fockp,E_gam4)

!=======================================================================
      implicit none
! Input Variables
      integer Nchunks
      integer ng4,nebf,npbf
      double precision GM4ICR(ng4)
      double precision DE(nebf,nebf),DP(npbf,npbf)
! Variables Returned
      double precision focke(nebf,nebf),fockp(npbf,npbf)
      double precision E_gam4

! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4
      integer,allocatable :: loop_map(:,:)
      double precision xfocke(nebf,nebf),xfockp(npbf,npbf)
     
         
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------(
      E_gam4=0.0d+00
      xfocke=0.0d+00
      xfockp=0.0d+00
!---------INITIALIZE-FOCK-DATA-STRUCTURES-----------------------------)
!-----CHOP-UP-THE-CALCULATION-OF-GAM_3--------------------------------(
      do ichunk=1,Nchunks
      
         call loop_size(1,ng4,Nchunks,ichunk-1,istart,iend)
!        write(*,*)'after call loop size'
!        write(*,*)'NG3=',ng3
!        write(*,*)'Nchunks=',Nchunks
!        write(*,*)'ichunk=',ichunk
!        write(*,*)'istart=',istart
!        write(*,*)'iend=',iend

! Segment of ng3:
         ng4_seg=1+iend-istart

         if(allocated(loop_map)) deallocate(loop_map)
         allocate( loop_map(ng4_seg,10),stat=istat )
!        write(*,*) 'allocate loop_map: ',istat

! Nested loop compression for this chunk:
         Loopi=0
         imas=0
         do ip=1,npbf
         do jp=1,npbf
            do ie1=1,nebf
            do je1=1,nebf
               do ie2=1,nebf
               do je2=1,nebf
                  do ie3=1,nebf
                  do je3=1,nebf
                    do ie4=1,nebf
                    do je4=1,nebf

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

         call thread_EG4ICR(istart,iend,ng4_seg,ng4,nebf,npbf,
     x                      loop_map,DE,DP,GM4ICR,
     x                      xfocke,xfockp,E_gam4)

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)
!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(loop_map)) deallocate(loop_map)
!     if(allocated(GAM_3)) deallocate(GAM_3)
!-----CLEAN-UP-AND-RETURN---------------------------------------------)

!--Update-the-full-Fock-matrices--------------------------------------(
      call add2fock(nebf,xfocke,focke)
      call add2fock(npbf,xfockp,fockp)
!--Update-the-full-Fock-matrices--------------------------------------)


      return
      end
!======================================================================
      subroutine thread_EG4ICR(istart,iend,ng4_seg,ng4,nebf,npbf,
     x                         loop_map,DE,DP,GM4ICR,
     x                         xfocke,xfockp,E_gam4)

!======================================================================
      implicit none
      include 'omp_lib.h'

! Input Variables
      integer istart,iend,ng4_seg
      integer nebf   ! Number contracted electronic basis functions
      integer npbf   ! Number nuclear basis functions
      integer ng4

      integer loop_map(ng4_seg,10)

      double precision DE(nebf,nebf),DP(npbf,npbf)
      double precision GM4ICR(ng4)

! Variables Returned
      double precision xfocke(nebf,nebf),xfockp(npbf,npbf)
      double precision E_gam4

! Local Variables
      integer ip,jp
      integer ie1,je1  !
      integer ie2,je2  ! Contracted elec basis function indices
      integer ie3,je3  !
      integer ie4,je4  !
      integer imap,ia
      double precision ans
      double precision four
!     parameter(four=4.0d+00)

!---OPENMP-RELATED-VARIABLES-----(
!     integer IFIL
!     integer id
      integer loopi,iLP
!     double precision wtime
!---OPENMP-RELATED-VARIABLES-----)

      four=4.0d+00

!---OPENMP-TIMING------------------------------------------------------(
!     wtime = omp_get_wtime()
!---OPENMP-TIMING------------------------------------------------------)

!--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(istart,iend)
!$ompx shared(nebf,npbf,ng4_seg)
!$ompx shared(ng4)
!$ompx shared(four)
!$ompx shared(DE,DP)
!$ompx shared(GM4ICR)
!$ompx private(iLp) 
!$ompx private(imap)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1)
!$ompx private(ie2,je2)
!$ompx private(ie3,je3)
!$ompx private(ie4,je4)
!$ompx private(ia)
!$ompx private(ans)
!$ompx reduction(+:E_gam4)
!$ompx reduction(+:xfocke)
!$ompx reduction(+:xfockp)

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

         call index_GAM_4PK2(nebf,npbf,
     x                       ip,jp,
     x                       ie1,je1,
     x                       ie2,je2,
     x                       ie3,je3,
     x                       ie4,je4,ia)


         ans=GM4ICR(ia)
!        ans=ans/(four*four)

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
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)


      return
      end     
