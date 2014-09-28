!=======================================================================
      subroutine RXCHF_GAM2_pass_MPI(nproc,rank,
     x                               ng2_int,ng2_xchf,
     x                               ng2loc_int,ng2loc_xchf,
     x                               nebf,nebfBE,npbf,
     x                               mpistart_int,mpiend_int,
     x                               XGM2_1,XGM2s,
     x                               GM2_1,GM2s)

! Pass XCHF integrals calculated for ng2_int/ng2loc_int range (XGM2*)
! to arrays for ng2_xchf/ng2loc_xchf range (GM2*)
!
! This erroneously assumes that all of the required integrals are stored
! contiguously at the start of the XGM2* arrays
! This is unfortunately not the case!
!=======================================================================
      implicit none
      include 'mpif.h'

! Input Variables
      integer nproc,rank
      integer ng2_int,ng2_xchf
      integer ng2loc_int,ng2loc_xchf
      integer nebf,nebfBE,npbf
      integer mpistart_int,mpiend_int   ! Indices wrt ng2_int/ng2loc_int
      double precision XGM2_1(ng2loc_int),XGM2s(ng2loc_int)

! Output Variables
      double precision GM2_1(ng2loc_xchf),GM2s(ng2loc_xchf)

! Local Variables
      integer i,j
      integer nints,indshift
      integer nints_recv,indshift_recv
      integer currstart,currend
      integer startrank,endrank
      integer mpistart_xchf,mpiend_xchf ! Indices wrt ng2_xchf/ng2loc_xchf

#if MPI32
      integer*4 ierr
#else
      integer ierr
#endif

! Get target start/end indices
      call get_mpi_range(ng2_xchf,nproc,rank,mpistart_xchf,mpiend_xchf)
      if(rank.eq.(nproc-1)) mpiend_xchf=ng2_xchf

! Iterate over each chunk of target integrals
      do i=0,nproc-1

C ARS( testing
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C )

! Broadcast start/end indices of this chunk
        if(rank.eq.i) then
         currstart=mpistart_xchf
         currend=mpiend_xchf
         write(*,*) "i,mpistart_int,mpiend_int,currstart,currend:",
     x               i,mpistart_int,mpiend_int,currstart,currend
        end if
        call MPI_BCAST(currstart,1,MPI_INTEGER8,i,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(currend,1,MPI_INTEGER8,i,MPI_COMM_WORLD,ierr)

! Get start/end rank of processes with integrals for this chunk
        call get_mpi_proc(ng2_int,nproc,currstart,startrank)
        call get_mpi_proc(ng2_int,nproc,currend,endrank)
        if(rank.eq.i) write(*,*) "startrank,endrank:",startrank,endrank

! If some part of integrals are on this process, pass
        if (.not.((mpistart_int.gt.currend).or.
     x            (mpiend_int.lt.currstart))) then

         if (mpistart_int.le.currstart) then
          nints=min(mpiend_int,currend)-currstart+1
          indshift=0
C ARS( testing
          write(*,*) "rank,nints1,indshift1:",rank,nints,indshift
C )

          if(rank.eq.i) then
           do j=1,nints
             GM2_1(indshift+j)=XGM2_1(currstart-mpistart_int+j)
             GM2s(indshift+j)=XGM2s(currstart-mpistart_int+j)
           end do
           write(*,*) "passed ints1 in memory:",rank
          else
           write(*,*) "proc ",rank," sending to ",i
           write(*,*) "sending nints:",rank,nints
           call MPI_SEND(nints,1,MPI_INTEGER8,i,rank,
     x                   MPI_COMM_WORLD,ierr)
           write(*,*) "sent"
           write(*,*) "sending indshift:",rank,indshift
           call MPI_SEND(indshift,1,MPI_INTEGER8,i,rank+nproc,
     x                   MPI_COMM_WORLD,ierr)
           write(*,*) "sent"
           write(*,*) "sending xgm2_1:",rank
           call MPI_SEND(XGM2_1(currstart-mpistart_int+1),nints,
     x                   MPI_DOUBLE_PRECISION,
     x                   i,rank+2*nproc,
     x                   MPI_COMM_WORLD,ierr)
           write(*,*) "sent"
           write(*,*) "sending xgm2s:",rank
           call MPI_SEND(XGM2s(currstart-mpistart_int+1),nints,
     x                   MPI_DOUBLE_PRECISION,
     x                   i,rank+3*nproc,
     x                   MPI_COMM_WORLD,ierr)
           write(*,*) "sent"
          end if

         else if (mpiend_int.ge.currend) then ! if we branch here, mpistart_int > currstart
          nints=currend-mpistart_int+1
          indshift=mpistart_int-currstart
C ARS( testing
          write(*,*) "rank,nints2,indshift2:",rank,nints,indshift
C )

          if(rank.eq.i) then
           do j=1,nints
             GM2_1(indshift+j)=XGM2_1(j)
             GM2s(indshift+j)=XGM2s(j)
           end do
           write(*,*) "passed ints2 in memory:",rank
          else
           write(*,*) "proc ",rank," sending to ",i
           write(*,*) "sending nints:",rank,nints
           call MPI_SEND(nints,1,MPI_INTEGER8,i,rank,
     x                   MPI_COMM_WORLD,ierr)
           write(*,*) "sent"
           write(*,*) "sending indshift:",rank,indshift
           call MPI_SEND(indshift,1,MPI_INTEGER8,i,rank+nproc,
     x                   MPI_COMM_WORLD,ierr)
           write(*,*) "sent"
           write(*,*) "sending xgm2_1:",rank
           call MPI_SEND(XGM2_1(1),nints,
     x                   MPI_DOUBLE_PRECISION,
     x                   i,rank+2*nproc,
     x                   MPI_COMM_WORLD,ierr)
           write(*,*) "sent"
           write(*,*) "sending xgm2s:",rank
           call MPI_SEND(XGM2s(1),nints,
     x                   MPI_DOUBLE_PRECISION,
     x                   i,rank+3*nproc,
     x                   MPI_COMM_WORLD,ierr)
          end if

         end if

        end if

        if(rank.eq.i) then
         do j=startrank,endrank
           if(j.ne.i) then
            write(*,*) "proc ",i," receiving from ",j
            write(*,*) "receiving nints:",rank
            call MPI_RECV(nints_recv,1,MPI_INTEGER8,j,j,
     x                    MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            write(*,*) "received:",nints_recv
            write(*,*) "receiving indshift:",rank
            call MPI_RECV(indshift_recv,1,MPI_INTEGER8,j,j+nproc,
     x                    MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            write(*,*) "received:",indshift_recv
            write(*,*) "receiving gm2_1:",rank
            call MPI_RECV(GM2_1(indshift_recv+1),nints_recv,
     x                    MPI_DOUBLE_PRECISION,
     x                    j,j+2*nproc,
     x                    MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            write(*,*) "received"
            write(*,*) "receiving gm2s:",rank
            call MPI_RECV(GM2s(indshift_recv+1),nints_recv,
     x                    MPI_DOUBLE_PRECISION,
     x                    j,j+3*nproc,
     x                    MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            write(*,*) "received"
           end if
         end do
        end if

      end do

      return
      end

