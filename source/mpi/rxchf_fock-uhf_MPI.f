C======================================================================
      subroutine RXCHF_fock_uhf_MPI(nproc,rank,
     x                             LCMF,nebf,nebf2,nelec,ngee,
     x                             DE,GAM_ecore,GAM_ee,
     x                             focke,E_total,E_ecore,E_ee,
     x                              DE_alp,DE_beta,FAE_alp,FAE_beta)

C UHF Fock procedure (not NEO-HF)
C======================================================================
      implicit none

C Input variables
      integer           nproc,rank
      logical           LCMF
      integer           nebf
      integer           nebf2
      integer           ngee
      integer           nelec
      double precision  DE(nebf,nebf)
      double precision  DE_alp(nebf,nebf)
      double precision  DE_beta(nebf,nebf)
      double precision  FAE_alp(nebf,nebf)
      double precision  FAE_beta(nebf,nebf)
      double precision  GAM_ecore(nebf2)
      double precision  GAM_ee(ngee)

C Output variables
      double precision  focke(nebf,nebf)
      double precision  E_total
      double precision  E_ecore, E_ecore_tmp
      double precision  E_ee, E_ee_tmp

C Local variables
      integer           nebflt
      double precision  zero
      parameter(zero=0.0d+00)


C Initialize values
      E_total=zero
      E_ecore=zero
      E_ee=zero
      E_ecore_tmp=zero
      E_ee_tmp=zero
      focke=zero
      FAE_alp=0
      FAE_beta=0

C One-electron energy
      call E_from_GAM_ecore(nebf,nebf2,GAM_ecore,DE_alp,
     x                      FAE_alp,E_ecore_tmp)
      E_ecore = E_ecore + E_ecore_tmp

      call E_from_GAM_ecore(nebf,nebf2,GAM_ecore,DE_beta,
     x                       FAE_beta,E_ecore_tmp)
      E_ecore = E_ecore + E_ecore_tmp

C Two-electron energy
      if(nelec.gt.1) then
         call E_from_GAM_ee_uhf(nebf,ngee,GAM_ee,DE,FAE_alp,
     x                          E_ee_tmp,De_alp)
      E_ee = E_ee + E_ee_tmp
         call E_from_GAM_ee_uhf(nebf,ngee,GAM_ee,DE,FAE_beta,
     x                          E_ee_tmp,De_beta)
      E_ee = E_ee + E_ee_tmp
      end if

      E_total=E_ecore+E_ee

      if((rank.eq.0).and.(LCMF)) then

       call UFM_sym_check2(nebf,focke)

       nebflt=nebf*(nebf+1)/2


      end if
      return
      end


!======================================================================
      subroutine E_from_GAM_ee_uhf(nebf,ngee,
     * GAM_ee,De_total,focke,E_ee,De_spin)
!
!======================================================================
      implicit none

      double precision zero,two,half
      parameter(zero=0.0d+00,two=2.0d+00,half=5.0d-01)
      integer nebf
      integer ngee
      double precision GAM_ee(ngee)
      double precision DE_total(nebf,nebf)
      double precision DE_spin(nebf,nebf)
      double precision focke(nebf,nebf)
      double precision E_ee

! Local variables
      integer nebfLT
      integer ia1,ia2
      integer ie1,je1,ie2,je2
      double precision vee1
      double precision vee2
      double precision val1, val2
      double precision xfocke(nebf,nebf)

      E_ee=zero

      call zero_out(nebf,xfocke)

      do ie2=1,nebf
         do je2=1,nebf

            do ie1=1,nebf
               do je1=1,nebf

                  call pack_4D(nebf,nebf,nebf,
     x                         je2,ie2,je1,ie1,ia1)

                  call pack_4D(nebf,nebf,nebf,
     x                         je1,ie2,je2,ie1,ia2)

                  vee1=GAM_ee(ia1)
                  vee2=GAM_ee(ia2)
                  val1=vee1*DE_TOTAL(ie2,je2)
                  val2=vee2*DE_SPIN(ie2,je2)
                  xfocke(ie1,je1)=xfocke(ie1,je1)+(val1-val2)
                  E_ee=E_ee+half*DE_spin(ie1,je1)*(val1-val2)
               end do
            end do

         end do
      end do

!     nebfLT=(nebf+nebf*nebf)/2
!     write(*,*)'IN E_from_GAM_ee: EE Contribution to Elec FOCK Matrix:'
!     call print_my_fock(nebf,nebfLT,xfocke)

!  Update the full electronic Fock matrix
      call add2fock(nebf,xfocke,focke)

      return
      end
!=======================================================================
      subroutine RXCHFmult_construct_DE_uhf(nelec,nebf,CE,DE)

! Form electronic density matrix for UHF
!=======================================================================
      implicit none

! Input Variables
      integer nelec,nebf
      double precision CE(nebf,nebf)

! Variables Returned
      double precision DE(nebf,nebf)

! Local Variables
      integer i,j,k
      integer nocc
      double precision coeff
      double precision one,two
      parameter(one=1.0d+00,two=2.0d+00)

      coeff=one
      nocc=nelec
      
      DO i=1,nebf
         DO j=1,nebf
          DE(j,i)=0.0d+00
          do k=1,nocc
            DE(j,i) = DE(j,i) + coeff*CE(j,k)*CE(i,k)
          end do

         END DO
      END DO


      return
      end
!======================================================================
      subroutine RXCHFmult_guess_elec_uhf(nae,nbe,nebf,xxse,GAM_ecore,
     x                            DAE,DBE,CAE,CBE,
     x                          nae_alp,nae_beta,DAE_alp,DAE_beta,
     x                           CAE_alp,CAE_beta)
 
!     Diagonalize the core electron Hamiltonian
!     to construct initial regular and special electronic guess density
!======================================================================
      implicit none
! Input Variables
      integer nebf
      integer nae
c@@@krb
      integer nae_alp,nae_beta
      integer nbe
      double precision xxse(nebf,nebf)
      double precision GAM_ecore(nebf,nebf)
! Variables Returned
      double precision DAE(nebf,nebf)
c@@@krb
      double precision DAE_alp(nebf,nebf)
      double precision DAE_beta(nebf,nebf)
      double precision DBE(nebf,nebf)
      double precision CAE(nebf,nebf)
c@@@krb
      double precision CAE_alp(nebf,nebf)
      double precision CAE_beta(nebf,nebf)
      double precision CBE(nebf,nebf)
! Local variables
      integer i,j
      integer nocca,noccb
      double precision C(nebf,nebf)
      double precision EVF(nebf)
      double precision zero
      parameter(zero=0.0d+00)

c@@@KB commenting this out and just assuming N_alp > N_beta
c      if (nae.gt.1) then
c       nocca=nae/2
c      else
c       nocca=nae
c      end if
      nocca=nae_alp

      if (nbe.gt.1) then
       noccb=nbe/2
      else
       noccb=nbe
      end if

      DAE=zero
      DBE=zero
      DAE_alp  = zero
      DAE_beta = zero
      CAE=zero
      CAE_alp  = zero
      CAE_beta = zero
      CBE=zero

      call UROOTHAN(C,EVF,xxse,GAM_ecore,nebf)

! Store first nocca evectors as occ alpha reg elec vectors
      do i=1,nae_alp
        do j=1,nebf
          CAE_alp(j,i)=C(j,i)
        end do
      end do
! Store first nocca (or nocca-1) evectors as occ alpha reg elec vectors
c@@@krb eventually change this so they aren't always identical
      do i=1,nae_alp
        do j=1,nebf
          CAE_beta(j,i)=C(j,i)
        end do
      end do
! Store next noccb evectors as occ spec elec vectors
      do i=1,noccb
        do j=1,nebf
          CBE(j,i)=C(j,i+nocca)
        end do
      end do

! Store nebf-nocca-noccb remaining evectors as virt reg and spec elec vectors
      do i=nocca+noccb+1,nebf
        do j=1,nebf
          CAE_alp(j,i-noccb)=C(j,i)
c          CAE_beta(j,i-noccb-(nae_alp-nae_beta))=C(j,i)
          CAE_beta(j,i-noccb)=C(j,i)
          CBE(j,i-nocca)=C(j,i)
        end do
      end do

      call RXCHFmult_construct_DE_uhf(nae_alp,nebf,CAE_alp,DAE_alp)
      call RXCHFmult_construct_DE_uhf(nae_beta,nebf,CAE_beta,DAE_beta)
      
      
      call RXCHFmult_construct_DE(NBE,nebf,CBE,DBE)
      do i = 1, nebf
        do j = 1, nebf
          DAE(i,j) = DAE_alp(i,j)+DAE_beta(i,j)
        enddo
      enddo
      return
      end