!======================================================================
      subroutine DIIS_driver(nebf,nocca,noccb,nstore,iter,
     x                       vecAE,vecBE,Selec,
     x                       errvec,vecDIIS,err)
!
! Invokes DIIS procedure of Ionova & Carter (JCP 1995: 102, 1251)
!
!======================================================================
      implicit none

! Input variables
      integer nebf
      integer nocca,noccb
      integer nstore
      integer iter
      double precision vecAE(nebf,nebf)  ! regular electrons curr it
      double precision vecBE(nebf,nebf)  ! special electrons curr it
      double precision Selec(nebf,nebf)  ! overlap matrix in AO basis

! Input/output variables
      double precision errvec(nebf,nebf,nstore)  ! stored as (AO ind,MO ind,iter) where
      double precision vecDIIS(nebf,nebf,nstore) ! MOs stored in order {occA,occB,virt}
      double precision err                       ! Matrix norm of current errvec

! Local variables
      integer i,j,k,ia
      integer nvirt
      logical debug
      double precision vecAEocc(nebf,nocca)             !
      double precision vecBEocc(nebf,noccb)             ! curr it partition
      double precision vecEvirt(nebf,nebf-nocca-noccb)  !

      double precision vecAE0occ(nebf,nocca)            !
      double precision vecBE0occ(nebf,noccb)            ! prev it partition
      double precision vecE0virt(nebf,nebf-nocca-noccb) !

      double precision errvecnorm
      double precision errveccurr(nebf,nebf)
      double precision errvecocc(nebf,nocca+noccb,nstore)

      double precision diiscoeffs(nstore)

      double precision newvec(nebf,nebf)

      double precision zero1(nebf)
      double precision zero2(nocca)
      double precision zero3(noccb)
      double precision zero4(nebf-nocca-noccb)
      double precision zero
      parameter(zero=0.0d+00)

      debug=.true.

      zero1=zero
      zero2=zero
      zero3=zero
      zero4=zero

      if(debug) then
       write(*,*) "Initial vecAE:"
       call PREVNU(vecAE,zero1,nebf,nebf,nebf)
       write(*,*) "Initial vecBE:"
       call PREVNU(vecBE,zero1,nebf,nebf,nebf)
      end if

! Auxiliary storage
      nvirt=nebf-nocca-noccb
      do i=1,nocca
        vecAEocc(:,i)=vecAE(:,i)
        vecAE0occ(:,i)=vecDIIS(:,i,nstore)
      end do
      do i=1,noccb
        vecBEocc(:,i)=vecBE(:,i)
        vecBE0occ(:,i)=vecDIIS(:,i+nocca,nstore)
      end do
      do i=1,nvirt
        vecEvirt(:,i)=vecBE(:,i+noccb)
        vecE0virt(:,i)=vecDIIS(:,i+nocca+noccb,nstore)
      end do

      if(debug) then
       write(*,*) "vecAE0occ before swap:"
       call PREVNU(vecAE0occ,zero2,nocca,nebf,nebf)
       write(*,*) "vecBE0occ before swap:"
       call PREVNU(vecBE0occ,zero3,noccb,nebf,nebf)
       write(*,*) "vecE0virt before swap:"
       call PREVNU(vecE0virt,zero4,nvirt,nebf,nebf)
       write(*,*) "vecAEocc before swap:"
       call PREVNU(vecAEocc,zero2,nocca,nebf,nebf)
       write(*,*) "vecBEocc before swap:"
       call PREVNU(vecBEocc,zero3,noccb,nebf,nebf)
       write(*,*) "vecEvirt before swap:"
       call PREVNU(vecEvirt,zero4,nvirt,nebf,nebf)
      end if

! Phase/align vecs from this iteration to the previous one
      call DIIS_swap(nebf,nocca,Selec,vecAE0occ,vecAEocc)
      call DIIS_swap(nebf,noccb,Selec,vecBE0occ,vecBEocc)
      call DIIS_swap(nebf,nvirt,Selec,vecE0virt,vecEvirt)
      if(debug) then
       write(*,*) "vecAE0occ after swap:"
       call PREVNU(vecAE0occ,zero2,nocca,nebf,nebf)
       write(*,*) "vecBE0occ after swap:"
       call PREVNU(vecBE0occ,zero3,noccb,nebf,nebf)
       write(*,*) "vecE0virt after swap:"
       call PREVNU(vecE0virt,zero4,nvirt,nebf,nebf)
       write(*,*) "vecAEocc after swap:"
       call PREVNU(vecAEocc,zero2,nocca,nebf,nebf)
       write(*,*) "vecBEocc after swap:"
       call PREVNU(vecBEocc,zero3,noccb,nebf,nebf)
       write(*,*) "vecEvirt after swap:"
       call PREVNU(vecEvirt,zero4,nvirt,nebf,nebf)
      end if

! Shift stored vecs from prev iterations and store latest
      do i=1,nstore-1
        vecDIIS(:,:,i)=vecDIIS(:,:,i+1)
      end do
      do i=1,nocca
        vecDIIS(:,i,nstore)=vecAEocc(:,i)
      end do
      do i=1,noccb
        vecDIIS(:,i+nocca,nstore)=vecBEocc(:,i)
      end do
      do i=1,nvirt
        vecDIIS(:,i+nocca+noccb,nstore)=vecEvirt(:,i)
      end do

! Get error vector of this iteration
      call DIIS_errvec(nebf,Selec,
     x                 vecDIIS(:,:,nstore-1),
     x                 vecDIIS(:,:,nstore),
     x                 errveccurr,err)

      if(debug) then
       write(*,*) "Current error vector:"
       call PREVNU(errveccurr,zero1,nebf,nebf,nebf)
      end if

! Shift stored errorvecs from prev iterations and store latest
      do i=1,nstore-1
        errvec(:,:,i)=errvec(:,:,i+1)
      end do
      do i=1,nocca
        errvec(:,i,nstore)=errveccurr(:,i)
      end do
      do i=1,noccb
        errvec(:,i+nocca,nstore)=errveccurr(:,i+nocca)
      end do
      do i=1,nvirt
        errvec(:,i+nocca+noccb,nstore)=errveccurr(:,i+nocca+noccb)
      end do

! If completed nstore iterations, start DIIS rotations
! In order to apply procedure for 1<iter<nstore, we resort to unfortunate hacks
      if(iter.ge.2) then

       errvecocc=zero
       diiscoeffs=zero
       newvec=zero

! Solve for mixing coefficients using only occ part of errvec
       ia=0
       do i=max(1,nstore-iter+1),nstore
         ia=ia+1
         do j=1,nocca+noccb
           do k=1,nebf
             errvecocc(k,j,ia)=errvec(k,j,i)
!!!!!!! screen on linear dependence (see Pulay open shell diis)
           end do
         end do
       end do
       call DIIS_solve(nebf,nocca,noccb,min(iter,nstore),
     x                 errvecocc(:,:,1:min(iter,nstore)),
     x                 diiscoeffs(1:min(iter,nstore)))

! Rotate orbitals using mixing coefficients and all of errvec
       call DIIS_rotate(nebf,min(iter,nstore),
     x                  vecDIIS(:,:,max(1,nstore-iter+1):nstore),
     x                  errvec(:,:,max(1,nstore-iter+1):nstore),
     x                  diiscoeffs(1:min(iter,nstore)),newvec)

! Pass new orbitals
       do i=1,nocca
         vecAE(:,i)=newvec(:,i)
       end do
       do i=1,noccb
         vecBE(:,i)=newvec(:,i+nocca)
       end do
       do i=1,nvirt
         vecAE(:,i+nocca)=newvec(:,i+nocca+noccb) ! Stores same set of virtuals
         vecBE(:,i+noccb)=newvec(:,i+nocca+noccb) ! in vecAE and vecBE
       end do

      end if

      if(debug) then
       write(*,*) "Final vecAE:"
       call PREVNU(vecAE,zero1,nebf,nebf,nebf)
       write(*,*) "Final vecBE:"
       call PREVNU(vecBE,zero1,nebf,nebf,nebf)
      end if

      return
      end
!======================================================================
      subroutine DIIS_swap(nebf,nvec,Selec,vec0,vec)
!
! Interchanges columns of vec to be in maximum correspondence with vec0
! Phasing is also performed if necessary
!
!======================================================================
      implicit none

! Input variables
      integer nebf,nvec
      double precision vec0(nebf,nebf)  ! reference orbitals (prev it)
      double precision Selec(nebf,nebf) ! overlap matrix in AO basis

! Input/output variables
      double precision vec(nebf,nebf)   ! target orbitals (curr it)

! Local variables
      integer i,j,k,l
      integer maxind
      double precision ovlap,maxovlap
      double precision tempvec(nebf)
      double precision zero,smallovlap
      parameter(zero=0.0d+00,smallovlap=1.0d-01)

      do i=1,nvec-1 ! iterate over bfs of previous it
        maxind=0
        maxovlap=zero
        tempvec=zero
        do j=i,nvec ! find bf of current it most like orb i
          ovlap=zero
          do k=1,nebf
          do l=1,nebf
            ovlap=ovlap+vec0(k,i)*vec(l,j)*Selec(k,l)
          end do
          end do
          if(abs(ovlap).gt.abs(maxovlap)) then
           maxind=j
           maxovlap=ovlap
          end if
        end do
        if(maxind.ne.i) then ! swap orbs if necessary
         do k=1,nebf
           tempvec(k)=vec(k,maxind)
           vec(k,maxind)=vec(k,i)
           vec(k,i)=tempvec(k)
         end do
        end if
        if(maxovlap.lt.zero) then ! phase if necessary
         do k=1,nebf
           vec(k,i)=-vec(k,i)
         end do
        end if
      end do

! Phase last orbital if needed
      ovlap=zero
      do k=1,nebf
      do l=1,nebf
        ovlap=ovlap+vec0(k,nvec)*vec(l,nvec)*Selec(k,l)
      end do
      end do
      if(ovlap.lt.zero) then
       do k=1,nebf
         vec(k,nvec)=-vec(k,nvec)
       end do
      end if

! Raise alarm if the last orbital is not close to the prev it
      if(abs(ovlap).lt.smallovlap) then
       write(*,*) "---WARNING---"
       write(*,*) "Last orbital after swapping is very different:"
       write(*,*) "   overlap:",ovlap
      end if

      return
      end
!======================================================================
      subroutine DIIS_errvec(nebf,Selec,vec0,vec,errvec,errnorm)
!
! Calculates error vector for current iteration as the principal
! logarithm (power series formula) of
!     A = vec0^t * Selec * (vec - vec0)
! which is first checked for convergence by calculating the norm of I-A
!
!======================================================================
      implicit none

! Input variables
      integer nebf
      double precision Selec(nebf,nebf) ! overlap matrix in AO basis
      double precision vec0(nebf,nebf)  ! MO coeffs (prev it)
      double precision vec(nebf,nebf)   ! MO coeffs (curr it)

! Output variables
      double precision errvec(nebf,nebf)
      double precision errnorm

! Local variables
      integer i,j,k
      logical convgd
      logical debug
      double precision coeff,matdiff
      double precision Ck(nebf,nebf),Ckt(nebf,nebf),Cn(nebf,nebf)
      double precision Cpow(nebf,nebf),aux(nebf,nebf),aux2(nebf,nebf)
      double precision errvec0(nebf,nebf)
      double precision zero1(nebf)
      double precision zeromat(nebf,nebf)

      integer maxpower ! maximum length of power series to eval logarithm
      parameter(maxpower=100)
      double precision zero,one,tol
      parameter(zero=0.0d+00,one=1.0d+00,tol=1.0d-06)

      debug=.true.

      if(debug) then
       write(*,*) "vec0:"
       call PREVNU(vec0,zero1,nebf,nebf,nebf)
       write(*,*) "vec:"
       call PREVNU(vec,zero1,nebf,nebf,nebf)
      end if

      Ck=zero
      Ckt=zero
      Cn=zero
      Cpow=zero
      aux=zero
      aux2=zero
      errvec=zero
      zero1=zero
      zeromat=zero

! Construct necessary matrices
      do i=1,nebf
      do j=1,nebf
        Ck(j,i)=vec0(j,i)
        Ckt(j,i)=vec0(i,j)
        Cn(j,i)=vec(j,i)
        aux(j,i)=Cn(j,i)-Ck(j,i)
      end do
      end do

      call RXCHF_matmult(nebf,nebf,nebf,nebf,
     x                   Selec,aux,aux2)
      call RXCHF_matmult(nebf,nebf,nebf,nebf,
     x                   Ckt,aux2,Cpow)

      if(debug) then
       write(*,*) "Arg-I of logarithm:"
       call PREVNU(Cpow,zero1,nebf,nebf,nebf)
      end if

! Check if norm of Cpow < 1 (or else power series not convergent)
      call DIIS_diffmatnorm(nebf,nebf,Cpow,zeromat,errnorm)
      if(debug) then
       write(*,*) "||Cpow|| = ",errnorm
      end if
      if(errnorm.ge.one) then
       write(*,*) "---  ERROR  ---"
       write(*,*) "Norm of logarithm argument is:",errnorm
       write(*,*) "This must be < 1 for principal logarithm to be"
       write(*,*) "defined as a power series! Exiting..."
       call ABRT
      end if

      errvec(:,:)=Cpow(:,:) ! first term in series

!      if(debug) then
!       write(*,*) "Error vector with series of length:",1
!       call PREVNU(errvec,zero1,nebf,nebf,nebf)
!      end if

      convgd=.false.
      i=2 ! exponent in power series expansion

      aux(:,:)=Cpow(:,:)

      do while (.not.convgd)

        errvec0(:,:)=errvec(:,:)

        if(mod(i,2).eq.0) then
         coeff=-one/dble(i)
        else
         coeff=one/dble(i)
        end if

! Form matrix polynomial
        call RXCHF_matmult(nebf,nebf,nebf,nebf,
     x                     Cpow,aux,aux2)
        aux(:,:)=aux2(:,:)

! Update error vector
        do j=1,nebf
        do k=1,nebf
          errvec(k,j)=errvec(k,j)+coeff*aux(k,j)
        end do
        end do

! Evaluate convergence using matrix norm of difference matrix
        call DIIS_diffmatnorm(nebf,nebf,errvec,errvec0,matdiff)
        if(matdiff.lt.tol) convgd=.true.

!        if(debug) then
!         write(*,*) "Error vector with series of length:",i
!         call PREVNU(errvec,zero1,nebf,nebf,nebf)
!         write(*,*) "has difference:",matdiff
!        end if

        i=i+1
        if((i.gt.maxpower).and.(.not.convgd)) then
         write(*,*) "--- WARNING ---"
         write(*,*) "Trouble converging logarithm after"
         write(*,*) "expanding power series to",maxpower," terms"
         write(*,*) "  latest difference (matrix norm):",matdiff
         convgd=.true.
        end if

      end do

! Evaluate DIIS error as current error vector norm
      call DIIS_diffmatnorm(nebf,nebf,errvec,zeromat,errnorm)

      return
      end
!======================================================================
      subroutine DIIS_solve(nebf,nocca,noccb,nstore,errvec,diiscoeffs)
!
! Solves linear system of equations (involving B matrix) to determine
! optimal mixing coefficients (mixing of error vectors of prev iters)
! Overlap of error matrices is taken to be the Frobenius inner product
!
!======================================================================
      implicit none

! Input variables
      integer nebf,nstore
      integer nocca,noccb
      double precision errvec(nebf,nocca+noccb,nstore)

! Output variables
      double precision diiscoeffs(nstore)

! Local variables
      integer i,j,k
      integer ierr,lwork
      double precision matnorm,matcond,matprod
      double precision Bmat(nstore+1,nstore+1)
      double precision Bmat0(nstore+1,nstore+1)
C      double precision Bmatfac(nstore+1,nstore+1)
      double precision aux(nstore+1,1)
      double precision c(nstore+1,1)
      double precision zero1(nstore+1)
      double precision zero2(nocca+noccb)
      double precision zeromat(nstore+1,nstore+1)
      double precision workq(1),ipiv(nstore+1),ipiv0(nstore+1)
      double precision, allocatable :: work(:)
      logical debug

      double precision zero,one
      parameter(zero=0.0d+00,one=1.0d+00)

      debug=.true.

      Bmat=zero
      Bmat0=zero
C      Bmatfac=zero
      c=zero
      zeromat=zero
      zero1=zero
      zero2=zero

      if(debug) then
       do i=1,nstore
         write(*,*) "Error vector at position:",i
         call PREVNU(errvec(:,:,i),zero2,nocca+noccb,nebf,nebf)
       end do
      end if

! Get B11 for normalization
      call DIIS_matprod(nebf,nocca+noccb,
     x                  errvec(:,:,1),errvec(:,:,1),matnorm)
      matnorm=one ! do not normalize for now

! Query workspace for matrix routines
      ierr=0
      workq=zero
      ipiv=zero
      call dsytrf("L",nstore+1,Bmat,nstore+1,ipiv,workq,-1,ierr)
      if(ierr.ne.0) then
       write(*,*) "Error in dsytrf query:",ierr
       return
      end if
      lwork=int(workq(1))
      if(allocated(work)) deallocate(work)
      allocate(work(lwork))

! Form B and c
      Bmat=zero
      do i=1,nstore
        Bmat(i+1,1)=-one/matnorm
        Bmat(1,i+1)=-one/matnorm
      end do
      do i=2,nstore+1
      do j=2,nstore+1
        call DIIS_matprod(nebf,nocca+noccb,
     x                    errvec(:,:,i-1),errvec(:,:,j-1),matprod)
        Bmat(i,j)=matprod/matnorm
      end do
      end do

      c(1,1)=-one/matnorm

      if(debug) then
       write(*,*) "Bmat:"
       call PREVNU(Bmat,zero1,nstore+1,nstore+1,nstore+1)
      end if

      if(debug) then
       write(*,*) "c:"
       call PREVNU(c,0,1,nstore+1,nstore+1)
      end if

      Bmat0(:,:)=Bmat(:,:)

! Factorize B (B overwritten with factor necessary for subsequent calls)
      ierr=0
      work=zero
      ipiv=zero
      call dsytrf("L",nstore+1,Bmat,nstore+1,ipiv,work,lwork,ierr)
      if(ierr.ne.0) then
       write(*,*) "Error in dsytrf:",ierr
       return
      end if

C      Bmatfac(:,:)=Bmat(:,:)
C      ipiv0(:)=ipiv(:)

C! Get condition number of B
C      call DIIS_diffmatnorm(nstore+1,nstore+1,Bmat0,zeromat,matnorm)
C      ierr=0
C      work=zero
C      matcond=zero
C      call dsycon("L",nstore+1,Bmat,nstore+1,ipiv,
C     x            matnorm,matcond,work,lwork,ierr)
C      if(ierr.ne.0) then
C       write(*,*) "Error in dsycon:",ierr
C       return
C      end if
C      write(*,*) "Reciprocal condition number of B:",matcond

! Solve Bx=c (c overwritten with solution after sytrs call)
      ierr=0
C      Bmat(:,:)=Bmatfac(:,:)
C      ipiv(:)=ipiv0(:)
      call dsytrs("L",nstore+1,1,Bmat,nstore+1,ipiv,c,nstore+1,ierr)
      if(ierr.ne.0) then
       write(*,*) "Error in dsytrs:",ierr
       return
      end if

      if(debug) then
       write(*,*) "x:"
       call PREVNU(c,0,1,nstore+1,nstore+1)
      end if

! Check c=Bx
      if(debug) then
        call RXCHF_matmult(nstore+1,nstore+1,nstore+1,1,
     x                     Bmat0,c,aux)
       write(*,*) "Bx:"
       call PREVNU(aux,0,1,nstore+1,nstore+1)
      end if

! Pass mixing coefficients
      do i=1,nstore
        diiscoeffs(i)=c(i+1,1)
      end do

      if(allocated(work)) deallocate(work)

      return
      end
!======================================================================
      subroutine DIIS_rotate(nebf,nstore,vecDIIS,errvec,diiscoeffs,vec)
!
! Determines optimal vec at current iteration by mixing vecs from prev
! iterations with DIIS coefficients
!
!======================================================================
      implicit none

! Input variables
      integer nebf,nstore
      integer nocca,noccb
      double precision vecDIIS(nebf,nebf,nstore) ! stored as (AO ind,MO ind,iter) where
      double precision errvec(nebf,nebf,nstore)  ! MOs stored in order {occA,occB,virt}
      double precision diiscoeffs(nstore)

! Output variables
      double precision vec(nebf,nebf)

! Local variables
      integer i,j,k,l
      integer maxind
      logical convgd
      logical debug
      double precision maxcoeff,errnorm
      double precision coeff,matdiff
      double precision Cpow(nebf,nebf),delmat(nebf,nebf)
      double precision Cex(nebf,nebf),Cex0(nebf,nebf)
      double precision aux(nebf,nebf),aux2(nebf,nebf)
      double precision oldvec(nebf,nebf)
      double precision newerrvec(nebf,nebf),zeromat(nebf,nebf)
      double precision zero1(nebf)

      integer maxpower ! maximum length of power series to eval exponential
      parameter(maxpower=100)
      double precision zero,tol
      parameter(zero=0.0d+00,tol=1.0d-06)

      debug=.true.

      Cpow=zero
      Cex=zero
      Cex0=zero
      aux=zero
      aux2=zero
      oldvec=zero
      newerrvec=zero
      vec=zero
      zero1=zero
      zeromat=zero

! Determine which of the nstore saved vectors has the largest DIIS coeff
      maxind=0
      maxcoeff=zero
      do i=1,nstore
        if(abs(diiscoeffs(i)).gt.maxcoeff) then
         maxind=i
         maxcoeff=abs(diiscoeffs(i))
        end if
      end do

! Construct necessary matrices
      do k=1,nstore
        if(k.ne.maxind) then

! Calculate del_ki from {del_j}
         delmat=zero
         if(k.gt.maxind) then
          do l=maxind+1,k
            do i=1,nebf
            do j=1,nebf
              delmat(j,i)=delmat(j,i)+errvec(j,i,l)
            end do
            end do
          end do
         else
          do l=k+1,maxind
            do i=1,nebf
            do j=1,nebf
              delmat(j,i)=delmat(j,i)-errvec(j,i,l)
            end do
            end do
          end do
         end if

! Weight {del_ki} with {tau_i}
         do i=1,nebf
         do j=1,nebf
           Cpow(j,i)=Cpow(j,i)+diiscoeffs(k)*delmat(j,i)
         end do
         end do

        end if
      end do

      if(debug) then
       write(*,*) "Arg of exponential:"
       call PREVNU(Cpow,zero1,nebf,nebf,nebf)
      end if

      do i=1,nebf
        Cex(i,i)=1.0d+00
        do j=1,nebf
          Cex(j,i)=Cex(j,i)+Cpow(j,i) ! first term in series (1+x)
        end do
      end do

!      if(debug) then
!       write(*,*) "Exponential with series of length:",1
!       call PREVNU(Cex,zero1,nebf,nebf,nebf)
!      end if

      convgd=.false.
      i=2 ! exponent in power series expansion

      aux(:,:)=Cpow(:,:)

      coeff=1.0d+00

      do while (.not.convgd)

        Cex0(:,:)=Cex(:,:)

        coeff=coeff*dble(i)

! Form matrix polynomial
        call RXCHF_matmult(nebf,nebf,nebf,nebf,
     x                     Cpow,aux,aux2)
        aux(:,:)=aux2(:,:)

! Update Cex
        do j=1,nebf
        do k=1,nebf
          Cex(k,j)=Cex(k,j)+1.0d+00/coeff*aux(k,j)
        end do
        end do

! Evaluate convergence using matrix norm of difference matrix
        call DIIS_diffmatnorm(nebf,nebf,Cex,Cex0,matdiff)
        if(matdiff.lt.tol) convgd=.true.

!        if(debug) then
!         write(*,*) "Exponential with series of length:",i
!         call PREVNU(Cex,zero1,nebf,nebf,nebf)
!         write(*,*) "has difference:",matdiff
!        end if

        i=i+1
        if((i.gt.maxpower).and.(.not.convgd)) then
         write(*,*) "--- WARNING ---"
         write(*,*) "Trouble converging exponential after"
         write(*,*) "expanding power series to",maxpower," terms"
         write(*,*) "  latest difference (matrix norm):",matdiff
         convgd=.true.
        end if

      end do

! Calculate final orbital by rotating C_maxind
      oldvec(:,:)=vecDIIS(:,:,maxind)
      call RXCHF_matmult(nebf,nebf,nebf,nebf,
     x                   oldvec,Cex,vec)
      if(debug) then
       write(*,*) "Pre-rotated MOs (occA,occB,virt) at curr it"
       call PREVNU(oldvec,zero1,nebf,nebf,nebf)
       write(*,*) "Rotated using relative iteration:",maxind
       write(*,*) "Rotated MOs (occA,occB,virt) at curr it"
       call PREVNU(vec,zero1,nebf,nebf,nebf)
      end if

! Update vecDIIS (last entry for curr it)
      vecDIIS(:,:,nstore)=vec(:,:)

C! Update errvec (last entry for curr it)
C      do k=1,nstore
C      do i=1,nebf
C      do j=1,nebf
C        newerrvec(j,i)=newerrvec(j,i)+diiscoeffs(k)*errvec(j,i,k)
C      end do
C      end do
C      end do
C      errvec(:,:,nstore)=newerrvec(:,:)
C
C      if(debug) then
C       write(*,*) "Updated error vector:"
C       call PREVNU(newerrvec,zero1,nebf,nebf,nebf)
C      end if
C
C! Evaluate DIIS error as current error vector norm
C      call DIIS_diffmatnorm(nebf,nebf,newerrvec,zeromat,errnorm)
C      if(debug) then
C       write(*,*) "Updated error vector norm:",errnorm
C      end if

      return
      end
!======================================================================
      subroutine DIIS_diffmatnorm(m,n,mat1,mat2,diffnorm)
!
! Calculates matrix norm of mat1-mat2 using Frobenius norm
!
!======================================================================
      implicit none

! Input variables
      integer m,n
      double precision mat1(m,n),mat2(m,n)

! Output variables
      double precision diffnorm

! Local variables
      integer i,j
      double precision diffmat(m,n)
      double precision zero
      parameter(zero=0.0d+00)

      diffnorm=zero

      do i=1,n
      do j=1,m
        diffmat(j,i)=mat1(j,i)-mat2(j,i)
      end do
      end do

      do i=1,n
        do j=1,m
          diffnorm=diffnorm+diffmat(j,i)*diffmat(j,i)
        end do
      end do

      diffnorm=dsqrt(diffnorm)

      return
      end
!======================================================================
      subroutine DIIS_matprod(m,n,mat1,mat2,innerprod)
!
! Calculates inner product <mat1,mat2> using Frobenius inner product
! <A|B> = tr(A^t B) which induces the Frobenius matrix norm
!
!======================================================================
      implicit none

! Input variables
      integer m,n
      double precision mat1(m,n),mat2(m,n)

! Output variables
      double precision innerprod

! Local variables
      integer i,j
      double precision mat1t(n,m)
      double precision aux(n,n)
      double precision zero
      parameter(zero=0.0d+00)

      innerprod=zero
      mat1t=zero
      aux=zero

      do i=1,n
      do j=1,m
        mat1t(i,j)=mat1(j,i)
      end do
      end do

      call RXCHF_matmult(n,m,m,n,mat1t,mat2,aux)

      do i=1,n
        innerprod=innerprod+aux(i,i)
      end do

      return
      end

