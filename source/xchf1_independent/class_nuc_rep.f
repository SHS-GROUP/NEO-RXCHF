C=======================================================================
      subroutine class_nuc_rep(nat,zan,cat)

C=======================================================================
      implicit none
C Input Variables
      integer nat
      double precision cat(3,nat)
      double precision zan(nat)
C Variables Returned
      double precision E_nuc_rep
C Local Variables
      integer nst
      integer i,j,k
      integer ni
      double precision R2
      double precision zero
      parameter(zero=0.0d+00)


      E_nuc_rep = zero

      if(nat.lt.2) goto 10

      nst = 2
      do i = nst,nat
         ni = i-1
         do j = nst-1,ni
            R2 = zero
            do k = 1,3
               R2 = R2+(cat(k,i)-cat(k,j))**2
            end do
            if(R2.ne.zero) E_nuc_rep = E_nuc_rep+zan(i)*zan(j)/SQRT(R2)
         end do
      end do


  10  CONTINUE
      open(800,file='ENUCRP.dat',status='unknown')
      write(800,*) E_nuc_rep
      close(800)


      return
      end
