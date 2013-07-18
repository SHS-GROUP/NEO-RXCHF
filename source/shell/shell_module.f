        module shell
! pindex: primitive index for the functions in the current shell.
! atm   : the index number of the atom, for which the current shell is
!         centered.
! nfunc : # of the primitive functions in the current shell
!        (s:1,p:3,d:6,f:10...)
! expt  : the Gaussian exponent
! coord : the Cartesian coordinates for the atom(center)
! ang   : ang. momentum number for the primtives in the current shell.
! peindex : primitive index for the functions belong to this shell


        type :: eshell
          integer :: atm,nfunc
          double precision :: expt
          double precision :: coord(3)
          integer,allocatable,dimension(:,:):: ang 
          integer,allocatable,dimension(:) :: peindex
        end type eshell

        type :: pshell
          integer :: atm,nfunc
          double precision :: expt
          double precision :: coord(3)
          integer,allocatable,dimension(:,:):: ang 
          integer,allocatable,dimension(:) :: peindex
        end type pshell

        double precision, parameter :: eps=1.0d-15

        end module shell


