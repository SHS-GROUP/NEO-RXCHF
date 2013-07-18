
      subroutine RTUV_sh(tLim,uLim,vLim,nLim,
     x                expt,alpha,XPC,YPC,ZPC,RTUV)
C
C=======================================================================
      implicit none

C Input variables
      integer tLim
      integer uLim
      integer vLim
      integer nLim
      double precision expt 
      double precision alpha 
      double precision XPC 
      double precision YPC 
      double precision ZPC 

C Variables Returned
      double precision RTUV(0:tLim,0:uLim,0:vLim)

C Local Variables
      integer t,u,v,n
      double precision Fmx
      double precision zero
      parameter(zero=0.0d+00)
      double precision two
      parameter(two=2.0d+00)
c     double precision F(0:nLim)
      double precision R(0:tLim+1,0:uLim+1,0:vLim+1,0:nLim+1)

C     Initialization
      R=zero

C     ===R(0,0,0,n)===
      do n=0,nLim
c         call iboys2(n,expt,Fmx)
         call boys_interpol(n,expt,Fmx)
         R(0,0,0,n)=(-two * alpha)**dble(n) * Fmx
      end do  

C     ===R(1,0,0,n)===
      t=1
      do n=0,nLim
         R(t,0,0,n)=XPC*R(t-1,0,0,n+1)
      end do  

C     ===R(t>=2,0,0,n)===
      do t=2,tLim
         do n=0,nLim
            R(t,0,0,n)=dble(t-1)*R(t-2,0,0,n+1) + XPC*R(t-1,0,0,n+1)
         end do  
      end do  

CCCCCCCCCCCCCC

C     ===R(t,u=1,0,n)===
      u=1
c     v=0
      do t=0,tLim
         do n=0,nLim
            R(t,u,0,n)=YPC*R(t,u-1,0,n+1)
         end do  
      end do  

C     ===R(t,u>=2,0,n)===
      do u=2,uLim
         do t=0,tLim
            do n=0,nLim
               R(t,u,0,n)=dble(u-1)*R(t,u-2,0,n+1) + YPC*R(t,u-1,0,n+1)
            end do  
         end do  
      end do
 
CCCCCCCCCCCCCC

C     ===R(t,u,v=1,n)===
      v=1
      do t=0,tLim
         do u=0,uLim
            do n=0,nLim
               R(t,u,v,n)=ZPC*R(t,u,v-1,n+1)
            end do  
         end do  
      end do  

C     ===R(t,u,v>=2,n)===
      do v=2,vLim
         do t=0,tLim
            do u=0,uLim
               do n=0,nLim
                R(t,u,v,n)=dble(v-1)*R(t,u,v-2,n+1) + ZPC*R(t,u,v-1,n+1)
               end do  
            end do  
         end do
      end do
 
CCCCCCCCCCCCCC

C  All we need:
C      RTUV_OUT=R(tLim,uLim,vLim,0)

C  Output
       do t=0,tLim
        do u=0,uLim
         do v=0,vLim
            RTUV(t,u,v)=R(t,u,v,0)
         enddo
        enddo
       enddo
C Print output for testing:
c     do t=0,tLim
c        do u=0,uLim
c           do v=0,vLim
c              do n=0,nLim
c               write(*,*)'R(',t,u,v,')=',RTUV_OUT(t,u,v)
c              end do  
c           end do  
c        end do
c     end do


      return
