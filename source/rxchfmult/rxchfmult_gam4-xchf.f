!======================================================================
      subroutine RXCHFmult_GAM4_XCHF(Nchunks,nebf,nebfBE,npbf,
     x                               ngee,ng2,ng4,
     x                               GAM_2s,GAM_4)

C Reads GAMee over all-electron basis and stores over restricted basis
C Allows all subsequent calls to be to original XCHF routines
!======================================================================
      implicit none
      include 'omp_lib.h'
! Input Variables
      integer Nchunks
      integer ng4,ng2,ngee,nebf,nebfBE,npbf
      double precision GAM_4(ng4)
      double precision GAM_2s(ng2)
      
! Local Variables
      integer istat,ichunk,istart,iend,ng4_seg
      integer Loopi,imas
      integer ngeeBE,ia,ib
      integer ip,jp,ie1,je1,ie2,je2,ie3,je3,ie4,je4
      integer,allocatable :: loop_map(:,:)
      double precision,allocatable :: GAM_ee(:)
      double precision,allocatable :: GAM_eeBE(:)
      double precision wtime,wtime2


      write(*,1000) ng4,nchunks,omp_get_num_procs(),
     x              omp_get_max_threads(),1
      wtime = omp_get_wtime()

C Read in GAMee over all-electron basis set
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      allocate( GAM_ee(ngee),stat=istat )
      call read_GAM_ee(ne,ngee,GAM_ee) 

C Store GAMee over reduced basis set
      ngeeBE=nebfBE*nebfBE*nebfBE*nebfBE
      if(allocated(GAM_eeBE)) deallocate(GAM_eeBE)
      allocate( GAM_eeBE(ngeeBE),stat=istat )

      do ie1=1,nebfBE
      do je1=1,nebfBE
        do ie2=1,nebfBE
        do je2=1,nebfBE

C Find place in all-electron array
          call pack_4D(nebf,nebf,nebf,je2,ie2,je1,ie1,ia)

C Find new place in restricted array
          call pack_4D(nebfBE,nebfBE,nebfBE,je2,ie2,je1,ie1,ib)

          GAM_eeBE(ib)=GAM_ee(ia)

        end do
        end do
      end do
      end do

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
         do ip=1,npbf
         do jp=1,npbf
            do ie1=1,nebfBE
            do je1=1,nebfBE
               do ie2=1,nebfBE
               do je2=1,nebfBE
                  do ie3=1,nebfBE
                  do je3=1,nebfBE
                    do ie4=1,nebfBE
                    do je4=1,nebfBE

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

         call thread_gam4_IC(nebfBE,npbf,ngeeBE,
     x                       ng2,ng4,ng4_seg,istart,iend,
     x                       loop_map,GAM_eeBE,GAM_2s,GAM_4)

         wtime2 = omp_get_wtime() - wtime2
         write(*,2000)ichunk,wtime2

      end do !end loop over chunks
!-----CHOP-UP-THE-CALCULATION-OF-GAM_4--------------------------------)

!-----CLEAN-UP-AND-RETURN---------------------------------------------(
      if(allocated(GAM_eeBE)) deallocate(GAM_eeBE)
      if(allocated(GAM_ee)) deallocate(GAM_ee)
      if(allocated(loop_map)) deallocate(loop_map)
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

