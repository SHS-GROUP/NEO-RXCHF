!======================================================================
      subroutine RXCHFmult_contr_mat(nebf,nebfBE,mat,matBE)
!
! Store restricted basis subset of matrix 
!            mat(nebf,nebf) : matrix over all bfs
!      matBE(nebfBE,nebfBE) : matrix over restricted bfs
! where each dim of mat is ordered as
!      1,...,nebfBE,nebfBE+1,...,nebf
!======================================================================
      implicit none

! Input Variables
      integer nebf
      integer nebfBE
      double precision mat(nebf,nebf)

! Variables Returned
      double precision matBE(nebfBE,nebfBE)

! Local variables
      integer i,j

      matBE=0.0d+00
      do i=1,nebfBE
      do j=1,nebfBE
        matBE(j,i)=mat(j,i)
      end do
      end do

      return
      end

!======================================================================
      subroutine process_elec_guess(nebf,nebfBE,elindBE,noccb)

! This processes the regular and special electronic guess MOs that must
! be provided by the user for restricted basis set runs
! This routine should not be called if LREORDGS=False
! The guesses are processed in this order:
!  1) reorders electronic guess MOs as per previous basis set reordering
!  2) truncates MO coeffs of spec elec guess that are not in the restr
!     basis set, orthogonalizes, then forms orth compl
!  3) projects reg elec guess onto orth compl from (2)
!
! Input:
!   integer nebf:         dimension of total space (all-AO basis)
!   integer nebfBE:       dimension of restricted space (restr basis)
!   integer elindBE:      see driver documentation
!   integer noccb:        number of occ special MOs
!   file    guessCAE.inp: regular MOs according to original basis set specification
!   file    guessCBE.inp: special MOs according to original basis set specification
!
! Output:
!   file    FinalCAE.dat: regular MOs according to new basis set specification
!   file    FinalCBE.dat: special MOs according to new basis set specification
!======================================================================
      implicit none

! Input variables
      integer nebf,nebfBE,noccb
      integer elindBE(nebfBE)

! Local variables
      integer i,j,k
      integer contrind,currcontrind
      logical laddbasis
      logical debug
      double precision ovlap
      double precision CA0(nebf,nebf)      ! reg  coeffs as given in input file
      double precision CB0(nebf,nebf)      ! spec coeffs as given in input file
      double precision CA(nebf,nebf-noccb) ! new regular coeffs
      double precision CB(nebfBE,noccb)    ! new special coeffs (occ only)
      double precision scrA(nebf,nebf)
      double precision scrB(nebf,nebf)
      double precision ocompl(nebf,nebf-noccb)
      double precision projvec(nebf-noccb,nebf-noccb)
      double precision CBpk(nebfBE,nebfBE)
      double precision S(nebf,nebf)
      double precision SBE(nebfBE,nebfBE)
      double precision vecA(nebf,nebf)
      double precision vecB(nebfBE,nebfBE)
      double precision junk(nebf,nebf)
      double precision zeroarr(nebf)
      double precision zero,one
      parameter(zero=0.0d+00,one=1.0d+00)

      debug=.true.

      CA=zero
      CB=zero
      scrA=zero
      scrB=zero
      ocompl=zero
      CBpk=zero
      projvec=zero
      vecA=zero
      vecB=zero
      junk=zero
      zeroarr=zero

! Read overlap matrix
      S=zero
      call read_elec_ovlap(nebf,S)
      call RXCHFmult_contr_mat(nebf,nebfBE,S,SBE)

! Read and reorder regular/special electron guesses

      call RXCHFmult_read_CAE(nebf,1,junk,CA0)
      call RXCHFmult_read_CBE(nebf,1,junk,CB0)

      write(*,*)
      write(*,*) "INITIAL REGULAR ELECTRONIC GUESS ORBITALS"
      write(*,*)
      call PREVNU(CA0,zeroarr,nebf,nebf,nebf)
      write(*,*)
      write(*,*) "INITIAL SPECIAL ELECTRONIC GUESS ORBITALS"
      write(*,*)
      call PREVNU(CB0,zeroarr,nebf,nebf,nebf)
      write(*,*)

      currcontrind=1
C First add bfs that are also in special electronic set
      do i=1,nebfBE
        contrind=elindBE(i)
        do j=1,nebf
          scrA(currcontrind,j)=CA0(contrind,j)
          scrB(currcontrind,j)=CB0(contrind,j)
        end do
        currcontrind=currcontrind+1
      end do
C Add remaining bfs for CA
C Omitting this for CB effectively truncates it to restr basis set
      do i=1,nebf
        laddbasis=.true.
        do j=1,nebfBE
          contrind=elindBE(j)
          if (i.eq.contrind) laddbasis=.false.
        end do
        if (laddbasis) then
         do j=1,nebf
           scrA(currcontrind,j)=CA0(i,j)
         end do
         currcontrind=currcontrind+1
        end if
      end do

! Repack special guess into condensed array
      do i=1,noccb
        do j=1,nebfBE
          CB(j,i)=scrB(j,i)
        end do
      end do

      if(debug) then
       write(*,*)
       write(*,*) "Reordered regular orbs before projection:"
       call PREVNU(scrA,zeroarr,nebf,nebf,nebf)
       write(*,*)
       write(*,*) "Reordered special orbs before orthonormalization:"
       call PREVNU(CB,zeroarr,noccb,nebfBE,nebfBE)
       write(*,*)
      end if

! Orthonormalize special MOs which are now nonorthogonal due to truncation
      call RXCHF_symmorth(nebfBE,noccb,CB,SBE)

! Final special electron MOs are now in CB - form orth compl to their span next
      do i=1,noccb
        do j=1,nebfBE
          CBpk(j,i)=CB(j,i)
        end do
      end do
      call calc_orth_compl(nebf,nebfBE,noccb,CBpk,S,ocompl)

      if(debug) then
       write(*,*)
       write(*,*) "Orthogonal complement of CB (cols):"
       call PREVNU(ocompl,zeroarr,nebf-noccb,nebf,nebf)
       write(*,*)
      end if

! Project and orthonormalize regular guess vectors in orth compl space
      call RXCHFmult_OCBSE_transVt(nebf,nebf-noccb,ocompl,
     x                             S,scrA,projvec)
! Project vectors back to all-AO space
!  - will be o-normalized and still span orth compl
      call dgemm('n','n',nebf,nebf-noccb,nebf-noccb,one,
     x           ocompl,nebf,projvec,nebf-noccb,zero,CA,nebf)

      if(debug) then
       write(*,*)
       write(*,*) "Checking orthormality of the special orbs:"
       do i=1,noccb
         do j=i,noccb
           call moovlap(nebfBE,CB(:,i),CB(:,j),SBE,ovlap)
           write(*,*) "   i,j,<spec_i|spec_j>:",i,j,ovlap
         end do
       end do
       write(*,*)
       write(*,*)
       write(*,*) "Checking orthormality of the regular orbs:"
       do i=1,nebf-noccb
         do j=i,nebf-noccb
           call moovlap(nebf,CA(:,i),CA(:,j),S,ovlap)
           write(*,*) "   i,j,<reg_i|reg_j>:",i,j,ovlap
         end do
       end do
       write(*,*)
       scrB=zero
       do i=1,noccb
         do j=1,nebfBE
           scrB(j,i)=CB(j,i)
         end do
       end do
       write(*,*) "Checking orthogonality of the special/regular orbs:"
       do i=1,noccb
         do j=1,nebf-noccb
           call moovlap(nebf,scrB(:,i),CA(:,j),S,ovlap)
           write(*,*) "   i,j,<spec_i|reg_j>:",i,j,ovlap
         end do
       end do
       write(*,*)
      end if

! Write out guess and also save on disk for reading from SCF driver
      do i=1,nebf-noccb
        do j=1,nebf
          vecA(j,i)=CA(j,i)
        end do
      end do
      do i=1,noccb
        do j=1,nebfBE
          vecB(j,i)=CB(j,i)
        end do
      end do

      write(*,*)
      write(*,*) "REORDERED REGULAR ELECTRONIC GUESS ORBITALS"
      write(*,*)
      call PREVNU(vecA,zeroarr,nebf,nebf,nebf)
      write(*,*)
      write(*,*) "writing to guessCAE-alt.inp ..."
      write(*,*)
      call write_MOs(880,nebf,vecA)

      write(*,*) "REORDERED SPECIAL ELECTRONIC GUESS ORBITALS"
      write(*,*)
      call PREVNU(vecB,zeroarr,nebfBE,nebfBE,nebfBE)
      write(*,*)
      write(*,*) "writing to guessCBE-alt.inp ..."
      write(*,*)
      call write_MOs(881,nebfBE,vecB)

      return
      end

!======================================================================
      subroutine calc_orth_compl(N,m,nvec,vec,S,orthvec)

! Calculates orthogonal complement of the (N-nvec)-dim subspace formed by
! the span of the first nvec vectors in {vec} in the global space given
! in all-AO coords
!  - calculated by finding nullspace of A^t*S
!  - A: matrix whose cols are vec given in R^N coords
!  - S: AO overlap matrix
!
! Input:
!   integer N:    dimension of total space
!                 (e.g. dimension of all-AO basis)
!   integer m:    dimension of subspace in whose coords the {vec} are given
!                 (e.g. dimension of restr-AO basis)
!   integer nvec: dimension of subspace to which orth compl will be found
!                 (e.g. # occ special vectors)
!   real*8  vec:  (m x m) array whose first nvec cols are vectors to which
!                 orth vecs will be found
!                 (e.g. all special vectors given in restr-AO basis coords)
!   real*8  S:    (m x m) matrix of overlap between bfs of total space
!                 (e.g. all-AO overlap matrix)
!
! Output:
!  real*8 orthvec:  (N x (N-vec)) array of vectors spanning orthogonal complement
!======================================================================

      implicit none

! Input variables
      integer N,m,nvec
      double precision vec(m,m)
      double precision S(N,N)

! Output variables
      double precision orthvec(N,N-Nvec)

! Local variables
      integer i,j,k
      integer istat
      integer dimnull
      logical debug
      double precision ovlap
      double precision A(N,nvec)
      double precision AtS(nvec,N)
      double precision aux(nvec,N)
      double precision svals(N),workq(1)
      double precision vt(N,N)
      double precision nullbas(N,N)
      double precision tvec(nvec)
      double precision zeroarr(N)
      double precision, allocatable :: work(:)
      double precision zero,one,tol
      parameter(zero=0.0d+00,one=1.0d+00,tol=1.0d-12)

      debug=.true.

      orthvec=zero
      A=zero
      AtS=zero
      aux=zero
      svals=zero
      vt=zero
      nullbas=zero
      tvec=zero
      zeroarr=zero
      workq=zero

! Construct A (cols are first nvec entries of {vec})
      do i=1,nvec
        do j=1,m
          A(j,i)=vec(j,i)
        end do
      end do

      if(debug) then
       write(*,*)
       write(*,*) "A matrix:"
       call PREVNU(A,zeroarr,nvec,N,N)
       write(*,*)
      end if

! Form A^t * S
      call dgemm('t','n',nvec,N,N,one,A,N,S,N,zero,AtS,nvec)

      if(debug) then
       write(*,*)
       write(*,*) "A^t * S matrix:"
       call PREVNU(AtS,zeroarr,N,nvec,nvec)
       write(*,*)
      end if

! Compute null space using SVD

! Query work array size for SVD and allocate work array
      istat=0
      aux(:,:)=AtS(:,:)
      call dgesvd("n","a",nvec,N,aux,nvec,
     x            svals,zero,1,vt,N,workq,-1,istat)
      if (istat.ne.0) then
       write(*,*) "Error in dgesvd query"
      end if
      if(allocated(work)) deallocate(work)
      allocate(work(int(workq(1))))

! Compute SVD
      istat=0
      aux(:,:)=AtS(:,:)
      call dgesvd("n","a",nvec,N,aux,nvec,
     x            svals,zero,1,vt,N,work,int(workq(1)),istat)
      if (istat.ne.0) then
       write(*,*) "Error in dgesvd"
      end if
      if(allocated(work)) deallocate(work)

      if(debug) then
       write(*,*)
       write(*,*) "singular values:"
       do i=1,N
         write(*,*) "i,s(i):",i,svals(i)
       end do
       write(*,*)
       write(*,*) "V^t matrix:"
       call PREVNU(vt,zeroarr,N,N,N)
       write(*,*)
      end if

! Find zero singular values and store corresponding rows of vt 
! in temporary array containing nullspace basis
      dimnull=0
      nullbas=zero
      do i=1,N
        if(abs(svals(i)).lt.tol) then
         dimnull=dimnull+1
         do j=1,N
           nullbas(j,dimnull)=vt(i,j)
         end do
        end if
      end do

! Dimension of null basis should be N-nvec
      if(dimnull.ne.(N-nvec)) then
       write(*,*) "ERROR: Dimension of null basis is incorrect"
       write(*,*) "       dimnull:",dimnull
       write(*,*) "       N-nvec :",N-nvec
      end if

! GS-orthogonalize basis
      orthvec(:,1)=nullbas(:,1)
      call moovlap(N,orthvec(:,1),orthvec(:,1),S,ovlap)
      do k=1,N
        orthvec(k,1)=orthvec(k,1)/dsqrt(ovlap)
      end do

      do i=2,dimnull
        orthvec(:,i)=nullbas(:,i)
        do j=i-1,1,-1
          call moovlap(N,orthvec(:,i),orthvec(:,j),S,ovlap)
          do k=1,N
            orthvec(k,i)=orthvec(k,i)-ovlap*orthvec(k,j)
          end do
        end do
        call moovlap(N,orthvec(:,i),orthvec(:,i),S,ovlap)
        do k=1,N
          orthvec(k,i)=orthvec(k,i)/dsqrt(ovlap)
        end do
      end do

      if(debug) then
       write(*,*)
       write(*,*) "GS-orthogonalized basis for null space (cols):"
       call PREVNU(orthvec,zeroarr,N-Nvec,N,N)
       write(*,*)
      end if

! Check nullity (i.e. check orthogonality to original set)
      do i=1,dimnull
        call dgemv('n',nvec,N,one,AtS,nvec,orthvec(:,i),1,zero,tvec,1)
        do j=1,nvec
          if (abs(tvec(j)).gt.tol) then
           write(*,*)
           write(*,*) " ***** WARNING ***** "
           write(*,*) "Computed null space vector is not correct!"
           write(*,*) "i:",i
           write(*,*) "x_i:",orthvec(:,i)
           write(*,*) "A^tS * x_i",tvec
           write(*,*)
          end if
        end do
      end do

! For overkill, manually check orthogonality of all computed vectors to each
! other and to the original set
      do i=1,dimnull
        do j=1,nvec
          call moovlap(N,orthvec(:,i),A(:,j),S,ovlap)
          if(abs(ovlap).gt.tol) then
           write(*,*)
           write(*,*) " ***** WARNING ***** "
           write(*,*) "Computed orth compl vector i is not"
           write(*,*) "orthogonal to initial set vector j!"
           write(*,*) "i:",i
           write(*,*) "x_i:",orthvec(:,i)
           write(*,*) "j:",j
           write(*,*) "v_j:",A(:,j)
           write(*,*) "<x_i | v_j>",ovlap
           write(*,*)
          end if
        end do
        do j=i+1,dimnull
          call moovlap(N,orthvec(:,i),orthvec(:,j),S,ovlap)
          if(abs(ovlap).gt.tol) then
           write(*,*)
           write(*,*) " ***** WARNING ***** "
           write(*,*) "Computed orth compl vector i is not"
           write(*,*) "orthogonal to orth compl vector j!"
           write(*,*) "i:",i
           write(*,*) "x_i:",orthvec(:,i)
           write(*,*) "j:",j
           write(*,*) "x_j:",orthvec(:,j)
           write(*,*) "<x_i | x_j>",ovlap
           write(*,*)
          end if
        end do
      end do

      return
      end

!======================================================================
      subroutine RXCHFmult_intersection(dimtot,nvec,vecA,ncanon,Sao,
     x                                  dimint,basint)
!
! Calculates an orthonormal basis for the intersection A \cap B where
!     A is spanned by the nvec columns of vecA
! and
!     B is spanned by the first ncanon canonical vectors of dim dimtot
!
! The intersection is calculated by forming a matrix [ vecA | {e_i} ]
! and calculating the null space using singular value decomposition
! (right singular vectors corresponding to zero singular values
! form an o-normal basis of the null space / intersection space)
!
! dimtot : dimension of basis in which vecA / canonical vectors are given
!   nvec : number of columns in vecA
!   vecA : (dimtot x nvec) matrix with columns corresponding to vectors
! ncanon : number of canonical vectors of dimension dimtot
!    Sao : overlap matrix in special electron AO basis
! dimint : dimension of intersection
! basint : orthonormal basis of intersection flattened to dim ncanon
!           - dimtot vectors are output
!           - the first dimint vectors are relevant
!======================================================================
      implicit none

! Input Variables
      integer          dimtot
      integer          nvec
      integer          ncanon
      double precision vecA(dimtot,nvec)
      double precision Sao(ncanon,ncanon)

! Variables Returned
      integer          dimint
      double precision basint(ncanon,dimtot) ! adjust on exit to dimint

! Local variables
      logical debug
      integer i,j,k
      integer m,n
      integer currind
      integer istat
      double precision ovlap
      double precision svals(max(dimtot,nvec+ncanon)),workq(1)
      double precision mat(dimtot,nvec+ncanon)
      double precision aux(dimtot,nvec+ncanon)
      double precision u(dimtot,dimtot)
      double precision s(dimtot,nvec+ncanon)
      double precision vt(nvec+ncanon,nvec+ncanon)
      double precision nullbas(nvec+ncanon,nvec+ncanon)
      double precision testvec(dimtot)
      double precision testmat(dimtot,nvec+ncanon)
      double precision zeroarr(dimtot)
      double precision, allocatable :: work(:)
      double precision, parameter   :: zero=0.0d+00, one=1.0d+00
      double precision, parameter   :: tol=1.0d-12

      debug=.true.

! Initialize
      svals=zero
      u=zero
      s=zero
      vt=zero
      aux=zero
      zeroarr=zero
      workq=zero

! Fill first nvec columns of mat with vecA
      do i=1,nvec
        do j=1,dimtot
          mat(j,i)=vecA(j,i)
        end do
      end do

! Fill next ncanon columns of mat with canonical vectors
      do i=1,ncanon
        do j=1,dimtot
          if (j.eq.i) then
           mat(j,i+nvec)=one
          else
           mat(j,i+nvec)=zero
          end if
        end do
      end do

! Query work array size for SVD and allocate work array
      m=dimtot
      n=nvec+ncanon

! GESVD
      aux(:,:)=mat(:,:)
      istat=0
      call dgesvd("A","A",m,n,mat,m,svals,u,m,vt,n,workq,-1,istat)
      if (istat.ne.0) then
       write(*,*) "Error in dgesvd query"
      end if
      if(allocated(work)) deallocate(work)
      allocate(work(int(workq(1))))

! Compute SVD
      aux(:,:)=mat(:,:)
      istat=0
      call dgesvd("A","A",m,n,aux,m,svals,u,m,vt,n,
     x            work,int(workq(1)),istat)
      if (istat.ne.0) then
       write(*,*) "Error in dgesvd"
      end if
      if(allocated(work)) deallocate(work)

      do i=1,min(m,n)
        s(i,i)=svals(i)
      end do

      if (debug) then
       write(*,*)
       write(*,*) "Input matrix:"
       write(*,*)
       call PREVNU(mat,zeroarr,n,m,m)

       write(*,*)
       write(*,*) "U matrix:"
       write(*,*)
       call PREVNU(u,zeroarr,m,m,m)

       write(*,*)
       write(*,*) "Sigma matrix:"
       write(*,*)
       call PREVNU(s,zeroarr,n,m,m)

       write(*,*)
       write(*,*) "V^t matrix:"
       write(*,*)
       call PREVNU(vt,zeroarr,n,n,n)
      end if

! Find zero singular values and store corresponding rows of vt 
! in temporary array containing nullspace basis
      dimint=0
      nullbas=zero
      do i=1,max(dimtot,nvec+ncanon)
        if(abs(svals(i)).lt.tol) then
         dimint=dimint+1
         do j=1,nvec+ncanon
           nullbas(j,dimint)=vt(i,j)
         end do
        end if
      end do

      if (debug) then
       write(*,*)
       write(*,*) "Product matrix:"
       call dgemm('n','n',m,n,n,one,s,m,vt,n,zero,testmat,m)
       call dgemm('n','n',m,n,m,one,u,m,testmat,m,zero,aux,m)
       call PREVNU(aux,zeroarr,n,m,m)
       write(*,*)
      end if

      if (debug) then
       write(*,*)
       write(*,*) "----------------------------------"
       write(*,*) " Dimension of intersection space:"
       write(*,'(2X,A,1X,I3)') "Max possible (nebfBE)      =",ncanon
       write(*,'(2X,A,1X,I3)') "Actual (after computation) =",dimint
       write(*,*) "----------------------------------"
       write(*,*)
      end if

      do i=1,dimint
        call dgemv('n',m,n,one,mat,m,nullbas(:,i),1,zero,testvec,1)
        do j=1,m
          if (abs(testvec(j)).gt.tol) then
           write(*,*)
           write(*,*) " ***** WARNING ***** "
           write(*,*) "Computed null space vector is not correct!"
           write(*,*) "i:",i
           write(*,*) "x_i:",nullbas(:,i)
           write(*,*) "Ax_i",testvec
           write(*,*)
          end if
        end do
      end do

! GS-orthogonalize basis
      basint(:,1)=nullbas(nvec+1:n,1)
      call moovlap(ncanon,basint(:,1),basint(:,1),Sao,ovlap)
      do k=1,ncanon
        basint(k,1)=basint(k,1)/dsqrt(ovlap)
      end do

      do i=2,dimint
        basint(:,i)=nullbas(nvec+1:n,i)
        do j=i-1,1,-1
          call moovlap(ncanon,basint(:,i),basint(:,j),Sao,ovlap)
          do k=1,ncanon
            basint(k,i)=basint(k,i)-ovlap*basint(k,j)
          end do
        end do
        call moovlap(ncanon,basint(:,i),basint(:,i),Sao,ovlap)
        do k=1,ncanon
          basint(k,i)=basint(k,i)/dsqrt(ovlap)
        end do
      end do

      return
      end

!======================================================================
      subroutine RXCHF_loworth(m,n,mat)

! Performs Lowdin orthogonalization for mat (replaced by symmetric mat)
! Computes SVD A = U * S * V^t and replaces A with U * V^t (uses GESVD)
! Assumes A(m x n) where the n columns of A are m-dim vecs to be orth
! Assumes an identity metric - use RXCHF_symmorth for orth wrt S_AO
!
!      m : rows of mat
!      n : cols of mat
!    mat : matrix whose columns will be orthogonalized
!======================================================================
      implicit none

! Input variables
      integer          :: m,n

! Input/Output variables
      double precision :: mat(m,n)

! Local variables
      integer          :: i,j
      integer          :: istat
      integer          :: matrank
      double precision :: svals(max(m,n)),workq(1)
      double precision :: u(m,m)
      double precision :: vt(n,n)
      double precision :: aux(m,n)

      double precision, allocatable :: work(:),red(:,:)
      double precision, parameter   :: zero=0.0d+00
      double precision, parameter   :: one=1.0d+00
      double precision, parameter   :: tol=1.0d-12

      svals=zero
      u=zero
      vt=zero
      aux=zero
      workq=zero

! Query work array size for SVD and allocate work array
      aux(:,:)=mat(:,:)
      istat=0
      call dgesvd("A","A",m,n,mat,m,svals,u,m,vt,n,workq,-1,istat)
      if (istat.ne.0) then
       write(*,*) "Error in dgesvd query"
      end if
      if(allocated(work)) deallocate(work)
      allocate(work(int(workq(1))))

! Compute SVD
      aux(:,:)=mat(:,:)
      istat=0
      call dgesvd("A","A",m,n,aux,m,svals,u,m,vt,n,
     x            work,int(workq(1)),istat)
      if (istat.ne.0) then
       write(*,*) "Error in dgesvd"
      end if
      if(allocated(work)) deallocate(work)

! Check that rank is n (since input columns should have been li)
      matrank=0
      do i=1,max(m,n)
        if(abs(svals(i)).ge.tol) then
         matrank=matrank+1
        end if
      end do
      if(matrank.ne.n) then
       write(*,*) "ERROR in RXCHF_loworth:"
       write(*,*) "   Calculated rank from SVD = ",matrank
       write(*,*) "       Number of input cols = ",n
       write(*,*) "Columns are possibly linearly dependent!"
       write(*,*) "svals:",svals
       write(*,*) "Exiting..."
       call abrt
      end if

! Form reduced form of either U or V^t and form orthogonalized mat
      mat=zero

      if (m.gt.n) then ! reduced form of U
       if(allocated(red)) deallocate(red)
       allocate(red(m,n))
       red=zero
       do i=1,n
       do j=1,m
         red(j,i)=u(j,i)
       end do
       end do
       call dgemm('n','n',m,n,n,one,red,m,vt,n,zero,mat,m)
       if(allocated(red)) deallocate(red)

      else if(n.gt.m) then ! reduced form of V^t
       if(allocated(red)) deallocate(red)
       allocate(red(m,n))
       red=zero
       do i=1,n
       do j=1,m
         red(j,i)=vt(j,i)
       end do
       end do
       call dgemm('n','n',m,n,m,one,u,m,red,m,zero,mat,m)
       if(allocated(red)) deallocate(red)

      else
       call dgemm('n','n',n,n,n,one,u,n,vt,n,zero,mat,n)
      end if

      return
      end

!======================================================================
      subroutine RXCHF_symmorth(m,n,mat,S)

! Performs symmetric orthogonalization for mat (replaced by symmetric mat)
! Computes SVD A = U * S * V^t and replaces A with U * V^t (uses GESVD)
! Assumes A(m x n) where the n columns of A are m-dim vecs to be orth
! Assumes a non-identity metric - can use RXCHF_loworth for S=I
!
!      m : rows of mat
!      n : cols of mat
!    mat : matrix whose columns will be orthogonalized
!      S : overlap matrix
!======================================================================
      implicit none

! Input variables
      integer          :: m,n
      double precision :: S(m,m)

! Input/Output variables
      double precision :: mat(m,n)

! Local variables
      integer          :: i,j
      integer          :: istat
      logical          :: debug
      double precision :: ovlap
      double precision :: Sevals(m),workq(1)
      double precision :: Sevecs(m,m)
      double precision :: Shalf(m,m)
      double precision :: Sminhalf(m,m)
      double precision :: aux(m,m)
      double precision :: transvec(m,n)
      double precision :: zeroarr(m)

      double precision, allocatable :: work(:)
      double precision, parameter   :: zero=0.0d+00
      double precision, parameter   :: one=1.0d+00
      double precision, parameter   :: tol=1.0d-12

      debug=.true.

      Sevals=zero
      Sevecs=zero
      Shalf=zero
      Sminhalf=zero
      aux=zero
      transvec=zero
      zeroarr=zero
      workq=zero

      if(debug) then
       write(*,*)
       write(*,*) "Matrix whose columns will be orthogonalized"
       call PREVNU(mat,zeroarr,n,m,m)
       write(*,*)
       do i=1,n
         do j=i,n
           call moovlap(m,mat(:,i),mat(:,j),S,ovlap)
           write(*,*) "i,j,<x_i|x_j>:",i,j,ovlap
         end do
       end do
       write(*,*)
      end if

      if(n.eq.1) then

! If one vector, then just normalize
       call moovlap(m,mat(:,1),mat(:,1),S,ovlap)
       do j=1,m
         mat(j,1)=mat(j,1)/dsqrt(ovlap)
       end do

      else

! Query work array size for diagonalization and allocate work array
       Sevecs(:,:)=S(:,:)
       istat=0
       call dsyev("v","l",m,Sevecs,m,Sevals,workq,-1,istat)
       if (istat.ne.0) then
        write(*,*) "Error in dsyev query"
       end if
       if(allocated(work)) deallocate(work)
       allocate(work(int(workq(1))))

! Diagonalize
       Sevecs(:,:)=S(:,:)
       istat=0
       call dsyev("v","l",m,Sevecs,m,Sevals,work,int(workq(1)),istat)
       if (istat.ne.0) then
        write(*,*) "Error in dsyev"
       end if
       if(allocated(work)) deallocate(work)

! Form S^(1/2) and S^(-1/2)
       do i=1,m
         Shalf(i,i)=dsqrt(Sevals(i))
         Sminhalf(i,i)=one/dsqrt(Sevals(i))
       end do

       aux=zero
       call dgemm('n','n',m,m,m,one,Sevecs,m,Shalf,m,zero,aux,m)
       Shalf=zero
       call dgemm('n','t',m,m,m,one,aux,m,Sevecs,m,zero,Shalf,m)
       if(debug) then
        write(*,*)
        write(*,*) "S^(1/2)"
        call PREVNU(Shalf,zeroarr,m,m,m)
        write(*,*)
       end if

       aux=zero
       call dgemm('n','n',m,m,m,one,Sevecs,m,Sminhalf,m,zero,aux,m)
       Sminhalf=zero
       call dgemm('n','t',m,m,m,one,aux,m,Sevecs,m,zero,Sminhalf,m)
       if(debug) then
        write(*,*)
        write(*,*) "S^(-1/2)"
        call PREVNU(Sminhalf,zeroarr,m,m,m)
        write(*,*)
       end if

! Transform vectors to orthogonal basis
       call dgemm('n','n',m,n,m,one,Shalf,m,mat,m,zero,transvec,m)
       if(debug) then
        write(*,*)
        write(*,*) "Transformed vectors:"
        call PREVNU(transvec,zeroarr,n,m,m)
        write(*,*)
       end if

! Lowdin orthogonalization (metric is now identity)
       call RXCHF_loworth(m,n,transvec)
       if(debug) then
        write(*,*)
        write(*,*) "Lowdin-orthogonalized transformed vectors:"
        call PREVNU(transvec,zeroarr,n,m,m)
        write(*,*)
       end if

! Back-transform from orthogonal basis
       mat=zero
       call dgemm('n','n',m,n,m,one,Sminhalf,m,transvec,m,zero,mat,m)
       if(debug) then
        write(*,*)
        write(*,*) "Lowdin-orthogonalized back-transformed vectors:"
        call PREVNU(mat,zeroarr,n,m,m)
        write(*,*)
       end if

! Normalize
       do i=1,n
         call moovlap(m,mat(:,i),mat(:,i),S,ovlap)
         do j=1,m
           mat(j,i)=mat(j,i)/dsqrt(ovlap)
         end do
       end do
       if(debug) then
        write(*,*)
        write(*,*) "Normalized vectors:"
        call PREVNU(mat,zeroarr,n,m,m)
        write(*,*)
        do i=1,n
          do j=i,n
            call moovlap(m,mat(:,i),mat(:,j),S,ovlap)
            write(*,*) "i,j,<x_i|x_j>:",i,j,ovlap
          end do
        end do
        write(*,*)
       end if

! Check orthogonality of all computed vectors to each other
       do i=1,n
         do j=i+1,n
           call moovlap(m,mat(:,i),mat(:,j),S,ovlap)
           if(abs(ovlap).gt.tol) then
            write(*,*)
            write(*,*) " ***** WARNING ***** "
            write(*,*) "Computed symm orth vector i is not"
            write(*,*) "orthogonal to symm orth vector j!"
            write(*,*) "i:",i
            write(*,*) "x_i:",mat(:,i)
            write(*,*) "j:",j
            write(*,*) "x_j:",mat(:,j)
            write(*,*) "<x_i | x_j>",ovlap
            write(*,*)
           end if
         end do
       end do

      end if

      return
      end

C======================================================================
      subroutine RXCHFmult_GAM_2PK(nebf,nebf2,npbf,
     x                             ip,jp,
     x                             ie1,je1,
     x                             ie2,je2,ia)
C
C Analog of index_GAM_2PK for different special electron basis
C======================================================================
      implicit none

      integer ia
      integer nebf
      integer nebf2
      integer npbf

      integer ip
      integer jp
      integer ie1
      integer je1
      integer ie2
      integer je2

C GAM_2(je2,ie2,je1,ie1,jp,ip)
C   ie1,je1=1,...,nebf
C   ie2,je2=1,...,nebf2
C     ip,jp=1,...,npbf

C  Array decalred as A(N1,N2,...)
C then A(i1,i2,...) = 
C    the (i1-1) + N1*((i2-1) + N2*(...))th element

        ia = 1 + (je2-1) + 
     x   nebf2*( (ie2-1) + 
     x   nebf2*( (je1-1) + 
     x    nebf*( (ie1-1) + 
     x    nebf*( (jp -1) +
     x    npbf*(  ip -1) ) ) ) )

      return
      end

C======================================================================
      subroutine RXCHFmult_GAM_3PK(nebf,nebf2,npbf,
     x                             ip,jp,
     x                             ie1,je1,
     x                             ie2,je2,
     x                             ie3,je3,ia)
C
C GAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
C======================================================================
      implicit none

      integer ia
      integer nebf
      integer nebf2
      integer npbf

      integer ip
      integer jp
      integer ie1
      integer ie2
      integer ie3
      integer je1
      integer je2
      integer je3

C GAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
C   ie1,je1=1,...,nebf
C   ie2,je2=1,...,nebf2
C   ie3,je3=1,...,nebf2
C     ip,jp=1,...,npbf

C  Array decalred as A(N1,N2,...)
C then A(i1,i2,...) = 
C    the (i1-1) + N1*((i2-1) + N2*(...))th element

      ia = 1 + (je3-1) + 
     x nebf2*( (ie3-1) + 
     x nebf2*( (je2-1) + 
     x nebf2*( (ie2-1) + 
     x nebf2*( (je1-1) + 
     x  nebf*( (ie1-1) + 
     x  nebf*( (jp -1) +
     x  npbf*(  ip -1) ) ) ) ) ) )

      return
      end

C======================================================================
      subroutine RXCHFmult_GAM_4PK(nebf,nebf2,npbf,
     x                             ip,jp,
     x                             ie1,je1,
     x                             ie2,je2,
     x                             ie3,je3,
     x                             ie4,je4,ia)
C
C GAM_4(je4,ie4,je3,ie3,je2,ie2,je1,ie1,jp,ip)
C======================================================================
      implicit none

      integer ia
      integer nebf
      integer nebf2
      integer npbf

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

C GAM_4(je4,ie4,je3,ie3,je2,ie2,je1,ie1,jp,ip)
C   ie1,je1=1,...,nebf
C   ie2,je2=1,...,nebf2
C   ie3,je3=1,...,nebf2
C   ie4,je4=1,...,nebf2
C     ip,jp=1,...,npbf

C  Array decalred as A(N1,N2,...)
C then A(i1,i2,...) = 
C    the (i1-1) + N1*((i2-1) + N2*(...))th element

      ia = 1 + (je4-1) + 
     x nebf2*( (ie4-1) + 
     x nebf2*( (je3-1) + 
     x nebf2*( (ie3-1) + 
     x nebf2*( (je2-1) + 
     x nebf2*( (ie2-1) + 
     x nebf2*( (je1-1) + 
     x  nebf*( (ie1-1) + 
     x  nebf*( (jp -1) +
     x  npbf*(  ip -1) ) ) ) ) ) ) ) )

      return
      end

