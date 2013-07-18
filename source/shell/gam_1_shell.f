C======================================================================
      subroutine GAM1_Shell(nebf,npebf,npbf,ng1,ng1prm,nat,ngtg1,
     x                       pmass,cat,zan,bcoef1,gamma1,
     x                       AMPEB2C,AGEBFCC,AGNBFCC,
     x                       ELCEX,NUCEX,ELCAM,NUCAM,ELCBFC,NUCBFC,
     $                       neshell,nnshell,esh,psh)

C======================================================================
      use shell
      implicit none
      include 'omp_lib.h'
c     include 'mpif.h'

C modified from gam_1_OMP.f
C g2_md_xggvec and
C g2_md_xggvee are done separately in the shell
C loops 

C Input Variables
      integer nebf   ! Number of contracted elec basis functions
      integer npbf   ! Number of nuclear basis functions
      integer npebf  ! Number of primitive elec basis functions
      integer ng1    ! Number of contracted GAM_1 integrals
      integer ng1prm ! Number primitive gamma1 integrals
      integer nat    ! Number of atoms
      integer ngtg1  ! Number BGammas
      integer :: neshell,nnshell  ! # of elec,nuc shells
C-------Basis Set Info-------(
      integer ELCAM(npebf,3)  ! Angular mom for electrons
      integer NUCAM(npbf,3)   ! Angular mom for quantum nuclei
      double precision ELCEX(npebf) ! Exponents: elec basis
      double precision NUCEX(npbf)  ! Exponents: nuc basis
      double precision ELCBFC(npebf,3) ! Basis centers: elec basis
      double precision NUCBFC(npbf,3)  ! basis centers: nuc basis
      integer AMPEB2C(npebf) ! Map primitive index to contracted
      double precision AGEBFCC(npebf) ! Map prim index to contract coef
      double precision AGNBFCC(npbf)  ! Nuclear contract coef
C-------Basis Set Info-------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)
!-------SHELL----------------------------------------------
      type(eshell) :: esh(neshell)
      type(pshell) :: psh(nnshell)
      integer :: ig,jg,scnt,ii,jj
!shell quartet info
      integer :: qshell(4,nnshell**2*neshell**2)
!elec/nuc shell idices
      integer :: ise,jse,isn,jsn
!b/gam's
      double precision :: gami,gamj,bi,bj
!elec shell pair quantities
      double precision :: KABe(neshell,neshell,3)
      double precision :: Pmate(neshell,neshell,3)
      double precision :: Qmate(neshell,neshell,3)
      double precision :: qABe(neshell,neshell)
      double precision :: Pe(neshell,neshell)
!Nuc shell pair quantities
      double precision :: KABn(nnshell,nnshell,3)
      double precision :: Pmatn(nnshell,nnshell,3)
      double precision :: Qmatn(nnshell,nnshell,3)
      double precision :: qABn(nnshell,nnshell)
      double precision :: Pn(nnshell,nnshell)
!prefactors
      double precision :: c0,qP1P2,alpha,expnt
      double precision :: X_P1P2,Y_P1P2,Z_P1P2
!-------SHELL---------------------------------------------------

C Local variables
      integer ia   ! Packing index
      integer ip   ! Nuc basis index
      integer jp   ! Nuc basis index
      integer ie1  ! Primitive elec basis func index 
      integer je1  ! Primitive elec basis func index 
      integer iec1 ! Contracted elec basis func index 
      integer jec1 ! Contracted elec basis func index 
      double precision GAM_1(ng1)
c     double precision XGAM_1(ng1)
      double precision GAM_1S(ng1)
c     double precision XGAM_1S(ng1)
      double precision Cof_ie1
      double precision Cof_je1
      double precision Cof_ip
      double precision Cof_jp
      double precision ans
      double precision ansKE,ansVE
      double precision ansS
      double precision zero
      parameter(zero=0.0d+00)
C--------------------------------(
C Basis set-related local variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3

      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision A3,Amat3(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
      double precision B3,Bmat3(3)
C--------------------------------)
C---OPENMP-RELATED-VARIABLES-----(
      integer id
      integer loopi,iLP
      double precision wtime
      integer loop_map(ng1prm,4)
C---OPENMP-RELATED-VARIABLES-----)


C---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime()
C---OPENMP-TIMING------------------------------------------------------)

c     IF(myid .eq. 0) THEN
         write(*,*)
         write(*,*)'**************************************'
         write(*,*)'       GAM_1 Integrals    '
         write(*,*)
         write(*,*)'METHOD = MD'
         write(*,*)'npebf  =',npebf
         write(*,*)'nebf   =',nebf
         write(*,*)'npbf   =',npbf
         write(*,*)'ng1    =',ng1
         write(*,*)' Available processors: ',omp_get_num_procs()
         write(*,*)' Available threads     ',omp_get_max_threads()
         write(*,*)' Threads in use        ',omp_get_num_threads()
         write(*,*)'**************************************'
         write(*,*)


C  zero out the 1-D arrays to hold the contracted integrals
            GAM_1=zero
            GAM_1s=zero
C Compress nested loops
         Loopi=0
         do ip=1,npbf
            do jp=1,npbf
               do ie1=1,npebf
                  do je1=1,npebf

                     Loopi=Loopi+1
                     loop_map(Loopi,1)=je1
                     loop_map(Loopi,2)=ie1
                     loop_map(Loopi,3)=jp
                     loop_map(Loopi,4)=ip

                  end do
               end do
            end do
         end do

C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------(
!$omp parallel 
!$ompx shared(loop_map)
!$ompx shared(ELCEX,ELCAM,ELCBFC,NUCEX,NUCAM,NUCBFC) 
!$ompx shared(AMPEB2C,AGEBFCC,AGNBFCC)
!$ompx shared(nat,ngtg1,pmass,cat,zan,bcoef1,gamma1)
!$ompx shared(nebf,npebf,npbf)
!$ompx shared(ng1)
!$ompx shared(ng1prm)
!$ompx shared(GAM_1,GAM_1s)
!$ompx private(iLp) 
!$ompx private(ia)
!$ompx private(ip,jp) 
!$ompx private(ie1,je1) 
!$ompx private(iec1,jec1)
!$ompx private(A1,I1,J1,K1,Amat1)
!$ompx private(A2,I2,J2,K2,Amat2)
!$ompx private(B1,L1,M1,N1,Bmat1)
!$ompx private(B2,L2,M2,N2,Bmat2)
!$ompx private(Cof_ie1,Cof_je1)
!$ompx private(Cof_ip,Cof_jp)
!$ompx private(ansE,ansS)
!$ompx private(id)
CCCCCCCCCCCCCCCCCC!$ompx reduction(+:GAM_1,GAM_1s)

      id= omp_get_thread_num()
      write(*,*)' Hello from process ',id
      if(id.eq.0) then
         write(*,*)'Threads in use', omp_get_num_threads()
      end if

!$omp do
      do iLP=1,ng1prm
C  Map loop indices
         je1=loop_map(iLP,1)
         ie1=loop_map(iLP,2)
         jp =loop_map(iLP,3)
         ip =loop_map(iLP,4)

C Get Basis set info:
         A1=ELCEX(ie1)
         I1=ELCAM(ie1,1)
         J1=ELCAM(ie1,2)
         K1=ELCAM(ie1,3)
         Amat1(1)=ELCBFC(ie1,1)
         Amat1(2)=ELCBFC(ie1,2)
         Amat1(3)=ELCBFC(ie1,3)

         A2=NUCEX(ip)
         I2=NUCAM(ip,1)
         J2=NUCAM(ip,2)
         K2=NUCAM(ip,3)
         Amat2(1)=NUCBFC(ip,1)
         Amat2(2)=NUCBFC(ip,2)
         Amat2(3)=NUCBFC(ip,3)

         B1=ELCEX(je1)
         L1=ELCAM(je1,1)
         M1=ELCAM(je1,2)
         N1=ELCAM(je1,3)
         Bmat1(1)=ELCBFC(je1,1)
         Bmat1(2)=ELCBFC(je1,2)
         Bmat1(3)=ELCBFC(je1,3)

         B2=NUCEX(jp)
         L2=NUCAM(jp,1)
         M2=NUCAM(jp,2)
         N2=NUCAM(jp,3)
         Bmat2(1)=NUCBFC(jp,1)
         Bmat2(2)=NUCBFC(jp,2)
         Bmat2(3)=NUCBFC(jp,3)

C  Calculate primitive integrals 
C  (Only the kinetic energy and overlap
C   Integrals)

         call xcalc_GAM1_KES(I1,J1,K1,A1,Amat1,
     x                      I2,J2,K2,A2,Amat2,
     x                      L1,M1,N1,B1,Bmat1,
     x                      L2,M2,N2,B2,Bmat2,
     x                      nat,ngtg1,
     x                      pmass,cat,zan,
     x                      bcoef1,gamma1,
     x                      ansKE,ansS)
                           

C  Map from primitive BF indices to contracted indices
         iec1=AMPEB2C(ie1)
         jec1=AMPEB2C(je1)

C  Get primitive Electron Basis Function Contraction Coefficients
         Cof_ie1=AGEBFCC(ie1)
         Cof_je1=AGEBFCC(je1)

C  Get Nuclear Basis Function Contraction Coefficients
         Cof_ip=AGNBFCC(ip)
         Cof_jp=AGNBFCC(jp)


C  Map the 4-index contracted integral to 1-D:
         call pack_4D(nebf,nebf,npbf,
     x                jec1,iec1,jp,ip,ia)

         GAM_1(ia)=GAM_1(ia)+ansKE*
     x             Cof_ip*Cof_jp*Cof_ie1*Cof_je1

         GAM_1s(ia)=GAM_1s(ia)+ansS*
     x              Cof_ip*Cof_jp*Cof_ie1*Cof_je1


      end do

!$omp end do
!$omp end parallel
C--------------%%%--PARALLEL--LOOPS--%%%-------------------------------)

C Now loop over the SHELL quartets and b/gam's.
C Formerly, g2_md_xggvec and g2_md_xggvee were used in one-integral-at
C -a-time fashion to calculate VCE,VCP, and VMX integrals.

!scnt:shell quartet count
      scnt=0
      do ise=1,neshell
       do jse=1,neshell
        do isn=1,nnshell
         do jsn=1,nnshell
            scnt=scnt+1
            qshell(1,scnt)=ise
            qshell(2,scnt)=jse
            qshell(3,scnt)=isn
            qshell(4,scnt)=jsn
         enddo
        enddo
       enddo
      enddo            

! shell pair quantities calculated before entering the gamma loops
! ELC SHELL PAIR
        do ii=1,neshell
          do jj=1,neshell
          call gauss_prod(esh(ii)%expt, esh(jj)%expt,
     $                    esh(ii)%coord,esh(jj)%coord,
     $                    KABe(ii,jj,:),qABe(ii,jj),Pe(ii,jj),
     $                    Pmate(ii,jj,:),Qmate(ii,jj,:))
          enddo
        enddo
! NUC SHELL PAIR
        do ii=1,nnshell
          do jj=1,nnshell
          call gauss_prod(psh(ii)%expt, psh(jj)%expt,
     $                    psh(ii)%coord,psh(jj)%coord,
     $                    KABn(ii,jj,:),qABn(ii,jj),Pn(ii,jj),
     $                    Pmatn(ii,jj,:),Qmatn(ii,jj,:))
          enddo
        enddo

!scnt=nnshell**2*neshell**2
      DO ig=1,ngtg1
        gami=gamma1(ig) 
          bi=bcoef1(ig)
        do ii=1,scnt 
         ise=qshell(1,ii); jse=qshell(2,ii)
         isn=qshell(3,ii); jsn=qshell(4,ii)

         call g2_vee_shell(Pe(ise,jse),Pn(isn,jsn),
     $              Pmate(ise,jse,:),Pmatn(isn,jsn,:),gami,
     $              ng1,GAM_1,ise,jse,isn,jsn,esh,psh,neshell,nnshell,
     $              AMPEB2C,AGEBFCC,AGNBFCC,npebf,npbf,
     $              bi)


!         g2_vec_shell contains atom loops inside.
!         Calculates 1e-CL Nuc attraction
!         and 1p-CL Nuc repulsion
         call g2_vec_shell(Pe(ise,jse),Pn(isn,jsn),
     $           Pmate(ise,jse,:),Pmatn(isn,jsn,:),gami,
     $           KABe(ise,jse,:),KABn(isn,jsn,:),qABe(ise,jse),
     $           qABn(isn,jsn),Qmate(ise,jse,:),Qmatn(isn,jsn,:),
     $           ng1,GAM_1,ise,jse,isn,jsn,esh,psh,neshell,nnshell,
     $           AMPEB2C,AGEBFCC,AGNBFCC,npebf,npbf,
     $           cat,zan,nat,bi)

        enddo ! ii loop

      DO jg=1,ngtg1
         gamj=gamma1(jg)
           bj=bocef1(jg)

        do jj=1,scnt 
         ise=qshell(1,jj); jse=qshell(2,jj)
         isn=qshell(3,jj); jsn=qshell(4,jj)

         call g2_vee_shell(Pe(ise,jse),Pn(isn,jsn),
     $                Pmate(ise,jse,:),Pmatn(isn,jsn,:),gami+gamj,
     $                ng1,GAM_1,ise,jse,isn,jsn,esh,psh,neshell,nnshell,
     $                AMPEB2C,AGEBFCC,AGNBFCC,npebf,npbf,
     $                bi,bj)

         call g2_vec_shell(Pe(ise,jse),Pn(isn,jsn),
     $              Pmate(ise,jse,:),Pmatn(isn,jsn,:),gami+gamj,
     $              ng1,GAM_1,ise,jse,isn,jsn,esh,psh,neshell,nnshell,
     $              AMPEB2C,AGEBFCC,AGNBFCC,npebf,npbf,
     $              cat,zan,nat,
     $              bi,bj)
                                 
                                 
                                 
        enddo ! jj loop

      ENDDO  ! jg loop
      ENDDO  ! ig loop



C---OPENMP-TIMING------------------------------------------------------(
      wtime = omp_get_wtime() - wtime
      write(*,*)'TIME TO CALCULATE GAM_1 INTEGRALS: ',wtime
C---OPENMP-TIMING------------------------------------------------------)

C  Open files to store integrals
         open(801,file='GAM_1.ufm',form='unformatted',
     x    status='unknown',access='direct',RECL=8)

         open(802,file='GAM_1s.ufm',form='unformatted',
     x    status='unknown',access='direct',RECL=8)

c        GAM_1=XGAM_1
c        GAM_1s=XGAM_1s

         do ip=1,npbf
            do jp=1,npbf
               do iec1=1,nebf
                  do jec1=1,nebf

                     call pack_4D(nebf,nebf,npbf,
     x                            jec1,iec1,jp,ip,ia)

C  Write to file
                     write(801,REC=ia) GAM_1(ia)
                     write(802,REC=ia) GAM_1s(ia)


                  end do
               end do
            end do
         end do

         close(801)
         close(802)

      return
      end

C======================================================================
      subroutine xcalc_GAM1_KES(I1,J1,K1,A1,Amat1,
     x                         I2,J2,K2,A2,Amat2,
     x                         L1,M1,N1,B1,Bmat1,
     x                         L2,M2,N2,B2,Bmat2,
     x                         nat,ngtg1,
     x                         pmass,cat,zan,
     x                         bcoef1,gamma1,
     x                         ansKE,ansS)

C======================================================================
      implicit none
C Input Variables
      integer nat,ngtg1
C--------------------------------(
C Basis set-related variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3

      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision A3,Amat3(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
      double precision B3,Bmat3(3)
C--------------------------------)
      double precision pmass    ! Mass of nonelectron quantum particle 
      double precision zan(nat) ! Classical nuclear charges
      double precision cat(3,nat) ! XYZ Coordinates of atoms
      double precision bcoef1(ngtg1)
      double precision gamma1(ngtg1)
C Variables Returned
      double precision ansKE
      double precision ansS
C Local Variables
      integer iii
      integer K,L
      logical debug

      double precision ans
      double precision ZNUC
      double precision Cmat(3)
      double precision OXX
      double precision OGX
      double precision OGG
      double precision TEXG  
      double precision TPXG  
      double precision TEGX  
      double precision TPGX  
      double precision TEGG  
      double precision TPGG  
      double precision GV 
      double precision BCFV 
      double precision GW 
      double precision BCFW 
      double precision ogxix,ogxiy,ogxiz 
      double precision tegxi,tpgxi 
      double precision texgi,tpxgi 
      double precision oggix,oggiy,oggiz 
      double precision teggi,tpggi
      double precision dum1,dum2,dum3,sval

      double precision S   
      double precision TE  
      double precision TP  
      double precision VCE 
      double precision VCP 
      double precision VMX 

      double precision ZERO,one 
      parameter(ZERO=0.0d+00,one=1.0d+00)


C     debug=.true.
      debug=.false.

      OGX=ZERO
      OGG=ZERO
      TEXG=ZERO
      TPXG=ZERO
      TEGX=ZERO
      TPGX=ZERO
      TEGG=ZERO
      TPGG=ZERO
C
      DO K=1,NGTG1
           GV=GAMMA1(K)
         BCFV=BCOEF1(K)

         CALL PGIKE1(I1,J1,K1,A1,Amat1,
     1               I2,J2,K2,A2,Amat2,
     2               L1,M1,N1,B1,Bmat1,
     3               L2,M2,N2,B2,Bmat2,
     4               gv,zero,one,tegxi,sval)
c    5               ogxix,ogxiy,ogxiz)
c        write(*,*)'electronic kinetic gx integral = ',tegxi
C
         CALL PGIKE1(I2,J2,K2,A2,Amat2,
     1               I1,J1,K1,A1,Amat1,
     2               L2,M2,N2,B2,Bmat2,
     3               L1,M1,N1,B1,Bmat1,
     4               gv,zero,pmass,tpgxi,sval)
c    5               ogxix,ogxiy,ogxiz)
C
c        OGX=OGX+BCFV*OGXIX*OGXIY*OGXIZ
         tegx=tegx+bcfv*tegxi
         tpgx=tpgx+bcfv*tpgxi

         CALL PGIKE1(I1,J1,K1,A1,Amat1,
     1               I2,J2,K2,A2,Amat2,
     2               L1,M1,N1,B1,Bmat1,
     3               L2,M2,N2,B2,Bmat2,
     4               zero,gv,one,texgi,sval)
c    5               dum1,dum2,dum3)
c        write(*,*)'electronic kinetic xg integral = ',texgi
C
         CALL PGIKE1(I2,J2,K2,A2,Amat2,
     1               I1,J1,K1,A1,Amat1,
     2               L2,M2,N2,B2,Bmat2,
     3               L1,M1,N1,B1,Bmat1,
     4               zero,gv,pmass,tpxgi,sval)
c    5               dum1,dum2,dum3)
C
         texg=texg+bcfv*texgi
         tpxg=tpxg+bcfv*tpxgi
C
         DO L=1,NGTG1
            GW=GAMMA1(L)
            BCFW=BCOEF1(L)
C
            CALL PGIKE1(I1,J1,K1,A1,Amat1,
     1                  I2,J2,K2,A2,Amat2,
     2                  L1,M1,N1,B1,Bmat1,
     3                  L2,M2,N2,B2,Bmat2,
     4                  gv,gw,one,teggi,sval)
c    5                  oggix,oggiy,oggiz)
c        write(*,*)'electronic kinetic gg integral = ',teggi
c        write(*,*)
C
            CALL PGIKE1(I2,J2,K2,A2,Amat2,
     1                  I1,J1,K1,A1,Amat1,
     2                  L2,M2,N2,B2,Bmat2,
     3                  L1,M1,N1,B1,Bmat1,
     4                  gv,gw,pmass,tpggi,sval)
c    5                  oggix,oggiy,oggiz)
C
c           OGG=OGG+BCFV*BCFW*OGGIX*OGGIY*OGGIZ
            tegg=tegg+bcfv*bcfw*teggi
            tpgg=tpgg+bcfv*bcfw*tpggi

         END DO
      END DO

      TE=TEGX+TEXG+TEGG
      TP=TPGX+TPXG+TPGG

cko      V terms are  calculated in the shell loops
cko      VCE=2.0d+00*VCEGX+VCEGG
cko      VCP=2.0d+00*VCPGX+VCPGG
cko      VMX=2.0d+00*VMGX+VMGG
cko      ansE=TE+TP-VCE+VCP-VMX

      ansKE=TE+TP

CCCCC  Overlap

      OXX=ZERO
      OGX=ZERO
      OGG=ZERO

C
      DO K=1,NGTG1
           GV=GAMMA1(K)
         BCFV=BCOEF1(K)

         call G2_MD_xggs(I1,J1,K1,A1,Amat1,
     1                   I2,J2,K2,A2,Amat2,
     2                   L1,M1,N1,B1,Bmat1,
     3                   L2,M2,N2,B2,Bmat2,
     4                   GV,zero,OGXIX)


         OGX=OGX+BCFV*OGXIX
C
         DO L=1,NGTG1
            GW=GAMMA1(L)
            BCFW=BCOEF1(L)
            call G2_MD_xggs(I1,J1,K1,A1,Amat1,
     1                      I2,J2,K2,A2,Amat2,
     2                      L1,M1,N1,B1,Bmat1,
     3                      L2,M2,N2,B2,Bmat2,
     4                      GV,GW,OGGIX)
C
            OGG=OGG+BCFV*BCFW*OGGIX
         end do
      end do

      S=OXX+2.0d+00*OGX+OGG
      ansS=s

      return
      end



