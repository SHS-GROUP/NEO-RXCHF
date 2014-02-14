C======================================================================
      subroutine RXCHF_fock_hf_MPI(nproc,rank,
     x                             LCMF,nebf,nebf2,nelec,ngee,
     x                             DE,GAM_ecore,GAM_ee,
     x                             focke,E_total,E_ecore,E_ee)

C HF Fock procedure (not NEO-HF)
C Exact duplicate of ../rxchfmult/rxchfmult_fock-hf.f except with
C rank passed to not have horrid amounts of output when LCMF is on
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
      double precision  GAM_ecore(nebf2)
      double precision  GAM_ee(ngee)

C Output variables
      double precision  focke(nebf,nebf)
      double precision  E_total
      double precision  E_ecore
      double precision  E_ee

C Local variables
      integer           nebflt
      double precision  zero
      parameter(zero=0.0d+00)


C Initialize values
      E_total=zero
      E_ecore=zero
      E_ee=zero
      focke=zero

C One-electron energy
      call E_from_GAM_ecore(nebf,nebf2,GAM_ecore,DE,focke,E_ecore)

C Two-electron energy
      if(nelec.gt.1) then
         call E_from_GAM_ee(nebf,ngee,GAM_ee,DE,focke,E_ee)
      end if

      E_total=E_ecore+E_ee

      if((rank.eq.0).and.(LCMF)) then

       call UFM_sym_check2(nebf,focke)

       nebflt=nebf*(nebf+1)/2
       write(*,*) "FAE HF:"
       call prt_lower_triangle(nebf,nebflt,focke)
       write(*,*)

      end if

      return
      end

