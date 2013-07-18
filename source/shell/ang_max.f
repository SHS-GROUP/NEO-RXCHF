      subroutine ang_max_pair(ise,jse,isn,jsn,esh,psh,nesh,nnsh,
     $                        angmx_e,angmx_n)
      use shell
      implicit none
      integer, intent(in) :: ise,jse,isn,jsn,nesh,nnsh
      type(eshell),intent(in) ::esh(nesh)
      type(pshell),intent(in) ::psh(nnsh)
      integer,dimension(3,2),intent(out) :: angmx_e,angmx_n 
! Local vars
      integer :: i,j,nf,tmp_mx,ss,si     
      angmx_e=0
      angmx_n=0

! 2 shells of elec

      DO ss=1,2
        if(ss==1)then
         si=ise  ! si: actual shell index
        else
         si=jse
        endif

        nf=esh(si)%nfunc
        do j=1,3
          tmp_mx=0
         do i=1,nf
            if(esh(si)%ang(i,j)>tmp_mx) tmp_mx=esh(si)%ang(i,j)
         enddo                            
         angmx_e(j,ss)=tmp_mx         
        enddo
      ENDDO
 
! 2 shells of NUC

      DO ss=1,2
        if(ss==1)then
         si=isn  
        else
         si=jsn
        endif

        nf=psh(si)%nfunc
        do j=1,3
          tmp_mx=0
         do i=1,nf
            if(psh(si)%ang(i,j)>tmp_mx) tmp_mx=psh(si)%ang(i,j)
         enddo                            
         angmx_n(j,ss)=tmp_mx         
        enddo
      ENDDO

      end subroutine
