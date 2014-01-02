      subroutine RXCHF_MPI_TESTALL(rank,nreqs,reqs)
      implicit none
      include 'mpif.h'

      integer    rank,nreqs
      integer*4  reqs(nreqs)

      integer    i,j
      integer*4  ierr
      integer*4  stat1(MPI_STATUS_SIZE)
      logical    flag
      integer*4, allocatable :: reqsloc(:)

      if(allocated(reqsloc)) deallocate(reqsloc)
      allocate(reqsloc(nreqs))
      do i=1,nreqs
        do j=1,nreqs
          reqsloc(j)=reqs(j)
        end do
        write(*,*) "loop:",rank,i
        write(*,*) "reqs1:",reqs(1),reqs(2),reqs(3)
        write(*,*) "reqsloc1:",rank,i,reqsloc(1),reqsloc(2),reqsloc(3)
      write(*,*) "allocated1:",rank,i,allocated(reqsloc)
        if (reqsloc(i).ne.MPI_REQUEST_NULL) then
         call MPI_TEST(reqsloc(i),flag,stat1,ierr)
C         call MPI_WAIT(reqsloc(i),stat1,ierr)
         write(*,*) "i,flag:",rank,i,flag
         write(*,*) "reqs2:",reqs(1),reqs(2),reqs(3)
         write(*,*) "reqsloc2:",rank,i,reqsloc(1),reqsloc(2),reqsloc(3)
      write(*,*) "allocated2:",rank,i,allocated(reqsloc)
        end if
      end do

      write(*,*) "after all:",rank

      if(allocated(reqsloc)) deallocate(reqsloc)

      return
      end

      subroutine RXCHF_MPI_WAITALL(rank,nreqs,reqs)
      implicit none
      include 'mpif.h'

      integer    rank,nreqs
      integer*4  reqs(nreqs)

      integer    i,j
      integer*4  ierr
      integer*4  stat1(MPI_STATUS_SIZE)
      integer*4, allocatable :: reqsloc(:)

      do i=1,nreqs
        if(allocated(reqsloc)) deallocate(reqsloc)
        allocate(reqsloc(nreqs))
        do j=1,nreqs
          reqsloc(j)=reqs(j)
        end do

C        write(*,*) "loop:",rank,i
C        write(*,*) "reqs1:",reqs(1),reqs(2),reqs(3)
C        write(*,*) "reqsloc1:",rank,i,reqsloc(1),reqsloc(2),reqsloc(3)
C      write(*,*) "allocated1:",rank,i,allocated(reqsloc)
        if (reqsloc(i).ne.MPI_REQUEST_NULL) then
         call MPI_WAIT(reqsloc(i),stat1,ierr)
C         call MPI_WAIT(reqsloc(i),stat1,ierr)
C         write(*,*) "i,flag:",rank,i,flag
C         write(*,*) "reqs2:",reqs(1),reqs(2),reqs(3)
C         write(*,*) "reqsloc2:",rank,i,reqsloc(1),reqsloc(2),reqsloc(3)
C      write(*,*) "allocated2:",rank,i,allocated(reqsloc)
        end if
      end do

C      write(*,*) "after all:",rank

      if(allocated(reqsloc)) deallocate(reqsloc)

      return
      end

      subroutine RXCHF_MPI_WAITALL2(rank,nreqs,reqs)
      implicit none
      include 'mpif.h'

      integer    rank,nreqs
      integer*4  reqs(nreqs)

      integer      i,j
      integer*4    req
      integer*4    ierr
      integer*4    stat1(MPI_STATUS_SIZE)
      integer      unitno
      character*3  istring

      write(istring,'(I3.3)') rank
      unitno=11+rank
      write(*,*) "nreqs:",rank,nreqs

      open(unit=unitno,file="reqs-"//istring//".ufm",form="unformatted")
      do i=1,nreqs
        write(unitno) reqs(i)
      end do
      close(unitno)

      open(unit=unitno,file="reqs-"//istring//".ufm",form="unformatted")
      do i=1,nreqs
        read(unitno) req
        if (req.ne.MPI_REQUEST_NULL) then
         call MPI_WAIT(req,stat1,ierr)
        end if
      end do
      close(unitno)

      return
      end

