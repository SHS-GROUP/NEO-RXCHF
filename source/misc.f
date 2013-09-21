!=======================================================================
      subroutine mpi_set(nproc,myid,ierr)
! Setup for MPI
!=======================================================================
      implicit none
!     include 'mpif.h'
! Variables Returned
      integer nproc,myid,ierr

!     call MPI_INIT(ierr)
!     call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
!     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

      return
      end
!=======================================================================
      subroutine mpi_final(ierr)

!=======================================================================
      implicit none
!     include 'mpif.h'
! Variables Returned
      integer ierr

!     call MPI_FINALIZE(ierr)

      return
      end
!======================================================================
      subroutine make_KPE(nebf,npebf,AMPEB2C,KPESTR,KPEEND)
 
!     Create a map between contracted EBF index and:
!     1) Beginning primitive index of the contracted shell:  KPESTR()
!     2) Ending primitive index of the contracted shell: KPEEND()
!======================================================================
      implicit none
! Input Variables
      integer nebf,npebf
      integer AMPEB2C(npebf)
! Variables Returned
      integer KPESTR(nebf)
      integer KPEEND(nebf)
! Local Variables
      integer ichk,ist,iend
      integer iep,iec


      ichk=1
      ist=1
      do iep=1,npebf
         iec=AMPEB2C(iep)
         if(iec.ne.ichk) then
            ichk=ichk+1
            iend=iep-1
            KPESTR(iec-1)=ist
            KPEEND(iec-1)=iend
            ist=iep
         end if
         if(iep.eq.npebf) then
            iend=iep
            KPESTR(iec)=ist
            KPEEND(iec)=iend
         end if
      end do


      return
      end
!=======================================================================
      subroutine NUCNORM(npbf,AGNBFCC,NUCEX,NUCAM,NUCBFC)
! Normalize the contraction coefficients of nuclear basis set
!=======================================================================
      implicit none
! Input Variables
      integer npbf
      integer NUCAM(npbf,3)  ! Angular mom for quantum nuclei
      double precision NUCEX(npbf) ! Exponents: nuc basis
      double precision NUCBFC(npbf,3) ! Basis centers: nuc basis
! Variables Returned
      double precision AGNBFCC(npbf) ! Map prim index to contract coef
! Local Variables
      integer ip
      integer I1,J1,K1
      double precision ans,A1,Amat1(3)

      do ip=1,npbf

         A1=NUCEX(ip)
         I1=NUCAM(ip,1)
         J1=NUCAM(ip,2)
         K1=NUCAM(ip,3)
!        Amat1(1)=NUCBFC(ip,1)
!        Amat1(2)=NUCBFC(ip,2)
!        Amat1(3)=NUCBFC(ip,3)
         Amat1(1)=0.0d+00
         Amat1(2)=0.0d+00
         Amat1(3)=0.0d+00

         call gfovlap(I1,J1,K1,A1,Amat1,
     2                I1,J1,K1,A1,Amat1,
     3                ans)

         ans=1.0d+00/sqrt(ans)
         AGNBFCC(ip)=AGNBFCC(ip)*ans
      end do

      return
      end

      subroutine get_mpi_range(n,nproc,rank,istart,iend)
      implicit none
      integer n,nproc,rank
      integer istart,iend

      istart=rank*(n/nproc)+1
      iend=(rank+1)*(n/nproc)

      return
      end

      subroutine copy_arr(n,arr1,arr)
! Copy n entries of arr1 into n-dimensional arr
      implicit none
      integer n
      double precision arr1(n), arr(n)
      integer i

      do i=1,n
        arr(i)=arr1(i)
      end do
      
      return
      end

