      subroutine g2_vec_shell(P1,P2,Pmat1,Pmat2,gamm,
     $                       KABe,KABn,qABe,qABn,Qmate,Qmatn,
     $                       ng1,GAM_1,ise,jse,isn,jsn,
     $                       esh,psh,nesh,nnsh,
     $                       AMPEB2C,AGEBFCC,AGNBFCC,npebf,npbf, 
     $                       ,cat,zan,natom,bi,bj)

      use shell
      implicit none
!     Ref: Persson and Taylor 1997 (TCA)
!     For a given shell quartet:
!      calculates either nuc-nuc repulsion with GTG
!      or e-(QM Nuc) attraction with GTG.
!     P1,P2,Pmat1,and Pmat2 are the shell pair quantities.
!     ise,jse,isn,jsn are the individual shell indices.
!     bi,bj are the bcoef's of GTG, and the latter is present only for
!       the nested gamma loop.
!     cat: nuclear Cartesian coordinates
!     zan: nuclear charge 
      integer, intent(in):: ng1,ise,jse,isn,jsn,npebf,npbf,nesh,nnsh
      integer, intent(in):: AMPEB2C(npebf),natom
      double precision, intent(in)::AGEBFCC(npebf),AGNBFCC(npbf)
      double precision, intent(in):: P1,P2,gamm,Pmat1(3),Pmat2(3),bi
      double precision, intent(in):: KABe(3),KABn(3),Qmate(3),Qmatn(3)
      double precision, intent(in):: qABe,qABn
      double precision, intent(in):: cat(3,natom),zan(natom)
      double precision, intent(in),optional ::bj 
      type(eshell),intent(in) ::esh(nesh)
      type(pshell),intent(in) ::psh(nnsh)
      integer, intent(inout)::GAM_1(ng1)
      
!! Local variables
      double precision :: c1,qP1P2,alpha,expnt,cmat(3),znuc
      double precision :: DC1,DC2,RCSx,RCSy,RCSz
      integer :: i,j,iep,jep,inp,jnp
      integer :: tlim,ulim,vlim      ! total(particles 1&2)  max ang momentum
      integer :: tlim2,ulim2,vlim2   ! max ang momentum (only for particle2)
      integer :: nef1,nef2,nnf1,nnf2 ! # of primitives of each shell
! tlim (max of t1+t2), 0<=t1<=(n1+n1'). 1:EShell Pair/2:NShell Pair
! angmx_e has maximum ang momentum number for the cartesian components of a given shell pair
! angmx_e(xyz(3):shell_index(2)). Similary for angmx_n (nuc shell pair)
      integer,dimension(3,2) :: angmx_e,angmx_n 
      double precision, allocatable,dimension(:,:,:) ::dKpp 
      double precision, allocatable,dimension(:,:,:) ::Rtuv
!  primitive info and Hermite coeff's
      integer :: I1,J1,K1,I2,J2,K2
      integer :: L1,M1,N1,L2,M2,N2,nloop,nloop2,ii,jj,kk
! Hermite index mapping in a 1D loop
      integer, allocatable :: hidx1(:,:),hidx2(:,:)
      double precision :: A1,B1,A2,B2
      double precision, allocatable,dimension(:) ::Et1,Eu1,Ev1
      double precision, allocatable,dimension(:) ::Et2,Eu2,Ev2
!Vec Herimte integrals for given primitives
      double precision :: ans
! The intermediate array arrI is I(k,g,h) in the McMurchie-Davidson paper
      double precision, allocatable,:: arrI(:,:,:)

!  Particle #1: Elec
!  Particle #2: QM Nuc
!  1/R1C operator(attraction) used first, and then 1/R2C (repulsion) operator 

      tlim=0;tlim2=0
      ulim=0;ulim2=0
      vlim=0;vlim2=0
     
      call ang_max_pair(ise,jse,isn,jsn,esh,psh,nesh,nnsh,
     $             angmx_e,angmx_n)


      do i=1,2
         tlim=tlim+angmx_e(1,i)+angmx_n(1,i)
         ulim=ulim+angmx_e(2,i)+angmx_n(2,i)
         vlim=vlim+angmx_e(3,i)+angmx_n(3,i)
      enddo

      do i=1,2
         tlim2=tlim2+angmx_n(1,i)
         ulim2=ulim2+angmx_n(2,i)
         vlim2=vlim2+angmx_n(3,i)
      enddo

! Now set up (N,L,M)k in the McMurchie-Davidson paper using tlim2, ulim2, and vlim2.
       nloop2=(tlim2+1)*(ulim2+1)*(vlim2+1)
       if(allocated(hidx2))deallocate(hidx2)
       allocate(hidx2(3,nloop2))

       nloop2=0
       do ii=0,tlim2
        do jj=0,ulim2
         do kk=0,vlim2
          nloop2=nloop2+1
          hidx2(1,nloop2)=ii
          hidx2(2,nloop2)=jj
          hidx2(3,nloop2)=kk
         enddo
        enddo
       enddo

      if(allocated(dKpp)deallocate(dKpp)
      allocate(dKpp(0:tlim,0:ulim,0:vlim))

! Derivatives of Kpp':exp(-qR2pp'), Eq.80
      qP1P2=P1*P2*gamm/(P1*P2+P1*gamm+P2*gamm)
      call deriv_kpp(tlim,ulim,vlim,Pmat1,Pmat2,qP1P2,dKpp)

      
      if(allocated(Rtuv)deallocate(Rtuv)
      allocate(Rtuv(0:tlim,0:ulim,0:vlim))

! # of the primitives in each shell (e.g. maximum S:1 P:3 D:6...)
      nef1=esh(ise)%nfunc
      nef2=esh(jse)%nfunc
      nnf1=psh(isn)%nfunc
      nnf2=psh(jsn)%nfunc

      if(allocated(arrI))deallocate(arrI)
      allocate(arrI(nloop2,nef1,nef2))
      
        A1=esh(ise)%expt
        B1=esh(jse)%expt
        A2=psh(isn)%expt
        B2=psh(jsn)%expt


! ATOM LOOP
      DO i=1,natom
         cmat=cat(:,i) 
         znuc=zan(i)
         call g2_vec_prefactors(P1,P2,
     $                    Pmat1,Pmat2,cmat,gamm,
     $                    c1,qP1P2,alpha,expnt,
     $                    DC1,DC2,
     $                    RCSx,RCSy,RCSz)

         call RTUV_sh(tlim,ulim,vlim,tlim+ulim+vlim,
     $                expnt,alpha,RCSx,RCSy,RCSz,Rtuv)

! shell member loop (g,h loop in the MD paper)
        do iep=1,nef1
          do jep=1,nef2
            I1=esh(ise)%ang(iep,1);L1=esh(jse)%ang(jep,1)
            J1=esh(ise)%ang(iep,2);M1=esh(jse)%ang(jep,2)
            K1=esh(ise)%ang(iep,3);N1=esh(jse)%ang(jep,3)
            if(allocated(Et1))deallocate(Et1)
            allocate(Et1(0:I1+L1))
            if(allocated(Eu1))deallocate(Eu1)
            allocate(Eu1(0:J1+M1))
            if(allocated(Ev1))deallocate(Ev1)
            allocate(Ev1(0:K1+N1))
            call ghec(I1,L1,KABe(1),A1,B1,P1,qABe,Qmate(1),Et1)
            call ghec(J1,M1,KABe(2),A1,B1,P1,qABe,Qmate(2),Eu1)
            call ghec(K1,N1,KABe(3),A1,B1,P1,qABe,Qmate(3),Ev1)

            nloop=(I1+L1+1)*(J1+M1+1)*(K1+N1+1)
            if(allocated(hidx1))deallocate(hidx1)
            allocate(hidx1(3,nloop))

            nloop=0
            do ii=0,I1+L1
             do jj=0,J1+M1
              do kk=0,K1+N1
                 nloop=nloop+1
                 hidx1(1,nloop)=ii
                 hidx1(2,nloop)=jj
                 hidx1(3,nloop)=kk
              enddo
             enddo
            enddo

            arrI=0.0d0
            do ii=1,nloop2 
             do jj=1,nloop
             call dherm_vec_sh(hidx1(1,jj),hidx1(2,jj),hidx1(3,jj),
     $                         hidx2(1,ii),hidx2(2,ii),hidx2(3,ii),
     $                         prefactors,....,ans)
              arrI(ii,iep,jep)=arrI(ii,iep,jep)+c1*ans
     $                       *Et1(hidx1(1,jj))
     $                       *Eu1(hidx1(2,jj))
     $                       *Ev1(hidx1(3,jj))
               enddo ! nloop 
            enddo ! nloop2
          enddo ! jep loop
        enddo ! iep loop

        do inp=1,nnf1
         do jnp=1,nnf2
            I2=psh(isn)%ang(inp,1);L2=psh(jsn)%ang(jnp,1)
            J2=psh(isn)%ang(inp,2);M2=psh(jsn)%ang(jnp,2)
            K2=psh(isn)%ang(inp,3);N2=psh(jsn)%ang(jnp,3)
            if(allocated(Et2))deallocate(Et2)
            allocate(Et2(0:I2+L2))
            if(allocated(Eu2))deallocate(Eu2)
            allocate(Eu2(0:J2+M2))
            if(allocated(Ev2))deallocate(Ev2)
            allocate(Ev2(0:K2+N2))

            call ghec(I2,L2,KABn(1),A2,B2,P2,qABn,Qmatn(1),Et2)
            call ghec(J2,M2,KABn(2),A2,B2,P2,qABn,Qmatn(2),Eu2)
            call ghec(K2,N2,KABn(3),A2,B2,P2,qABn,Qmatn(3),Ev2)

            do iep=1,nef1
             do jep=1,nef2
               do ii=1,nloop2 

             enddo ! jep loop
            enddo ! iep loop
         
         enddo ! jnp loop
        enddo ! inp loop


      ENDDO ! atom loop
