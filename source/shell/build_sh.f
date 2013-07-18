      subroutine build_shell(npebf,nshell,ncenter,ectrfst,ectrlst,
     $                         ELCAM,ELCEX,ELCBFC,
     $                         nssh,npsh,ndsh,nfsh,esh)
! Construct NSHELL shells
        use shell
        implicit none
        integer :: npebf,nshell,ncenter
! indices of starting and ending primitives in each center
        integer,dimension(ncenter) ::ectrfst,ectrlst 
! # of shells in each center
        integer,dimension(ncenter) ::nssh,npsh,ndsh,nfsh
        integer :: ELCAM(npebf,3)
        double precision ::ELCEX(npebf),ELCBFC(npebf,3)
        type(eshell) ::esh(nshell)
! Local variables
        integer :: iic,istrt,iend,ishell,scount,ik,i,j
        integer :: loopi,loopf,itms,itmp,itmd,itmf,ijprod
        integer,allocatable :: shfunc(:)
        integer :: psh_prim(3),dsh_prim(6),fsh_prim(10)

! ishell should be equal to nshell at the end
       ishell=0

       DO iic=1,ncenter
        itms=0;itmp=0;itmd=0;itmf=0
        istrt=ectrfst(iic)
        iend =ectrlst(iic)
cko S shell (consists of one primitive)
        do i=istrt,iend
         if(sum(ELCAM(i,:))==0)then
          itms=itms+1
          ishell=ishell+1
          esh(ishell)%nfunc=1
          esh(ishell)%expt=ELCEX(i) 
          esh(ishell)%atm=iic
          esh(ishell)%coord=ELCBFC(i,:)
          allocate(esh(ishell)%ang(1,3))
          esh(ishell)%ang=0
          allocate(esh(ishell)%peindex(1)) 
          esh(ishell)%peindex(1)=i
         endif
        enddo

cko P shell (up to three primitives)
         do i=istrt,iend
          if(sum(ELCAM(i,:))==1)itmp=itmp+1
         enddo
         allocate(shfunc(istrt+itms:istrt+itms+itmp-1))

         loopi=istrt+itms
         loopf=istrt+itms+itmp-1
         shfunc=1
cko loop i goes through only the unique exponents representing each
cko shell.
         do i=loopi,loopi+npsh(iic)-1
           psh_prim(1)=i
! scount counts # of components in a given shell
           scount=1
           do j=i+1,loopf
             ijprod=shfunc(i)*shfunc(j)
             if(abs(ELCEX(i)-ELCEX(j))<eps .and. ijprod==1) then
                shfunc(j)=0
                scount=scount+1
                psh_prim(scount)=j
             endif
           enddo
           ishell=ishell+1
           esh(ishell)%nfunc=scount
           esh(ishell)%expt=ELCEX(i)
           esh(ishell)%atm=iic
           esh(ishell)%coord=ELCBFC(i,:)
           allocate(esh(ishell)%ang(scount,3))
           allocate(esh(ishell)%peindex(scount))
           do ik=1,scount
            esh(ishell)%ang(ik,:)=ELCAM(psh_prim(ik),:)
            esh(ishell)%peindex(ik)=psh_prim(ik)
           enddo 
         enddo
         deallocate(shfunc)

cko D shell
        do i=istrt,iend
         if(sum(ELCAM(i,:))==2)itmd=itmd+1
        enddo

        allocate(shfunc(istrt+itms+itmp:istrt+itms
     $                  +itmp+itmd-1))
        loopi=istrt+itms+itmp
        loopf=istrt+itms+itmp+itmd-1
        shfunc=1
        do i=loopi,loopi+ndsh(iic)-1
           dsh_prim(1)=i
           scount=1
           do j=i+1,loopf
             ijprod=shfunc(i)*shfunc(j)
             if(abs(ELCEX(i)-ELCEX(j))<eps .and. ijprod==1) then
                shfunc(j)=0
                scount=scount+1
                dsh_prim(scount)=j
             endif
           enddo
           ishell=ishell+1
           esh(ishell)%nfunc=scount
           esh(ishell)%expt=ELCEX(i)
           esh(ishell)%atm=iic
           esh(ishell)%coord=ELCBFC(i,:)
           allocate(esh(ishell)%ang(scount,3))
           allocate(esh(ishell)%peindex(scount))
           do ik=1,scount
            esh(ishell)%ang(ik,:)=ELCAM(dsh_prim(ik),:)
            esh(ishell)%peindex(ik)=dsh_prim(ik)
           enddo
        enddo 
        deallocate(shfunc)

cko F shell
        do i=istrt,iend
         if(sum(ELCAM(i,:))==3)itmf=itmf+1
        enddo

        allocate(shfunc(istrt+itms+itmp+itmd:istrt+itms
     $                  +itmp+itmd+itmf-1))
        loopi=istrt+itms+itmp+itmd
        loopf=istrt+itms+itmp+itmd+itmf-1
        shfunc=1

        do i=loopi,loopi+nfsh(iic)-1
           fsh_prim(1)=i
           scount=1
           do j=i+1,loopf
             ijprod=shfunc(i)*shfunc(j)
             if(abs(ELCEX(i)-ELCEX(j))<eps .and. ijprod==1) then
                shfunc(j)=0
                scount=scount+1
                fsh_prim(scount)=j
             endif
           enddo
           ishell=ishell+1
           esh(ishell)%nfunc=scount
           esh(ishell)%expt=ELCEX(i)
           esh(ishell)%atm=iic
           esh(ishell)%coord=ELCBFC(i,:)
           allocate(esh(ishell)%ang(scount,3))
           allocate(esh(ishell)%peindex(scount))
           do ik=1,scount
            esh(ishell)%ang(ik,:)=ELCAM(fsh_prim(ik),:)
            esh(ishell)%peindex(ik)=fsh_prim(ik)
           enddo
        enddo 
        deallocate(shfunc)
        print *,'building shells',nshell,ishell
        if(nshell /=ishell)then
          print *,' # of shell mismatch, in,out',nshell,ishell
          stop
        endif
       ENDDO   

      end subroutine
