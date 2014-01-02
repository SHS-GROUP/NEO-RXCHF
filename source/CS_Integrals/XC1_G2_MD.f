C======================================================================
c     subroutine cws_gam1_xggs(I1,J1,K1,A1,Amat1,
      subroutine G2_MD_xggs(I1,J1,K1,A1,Amat1,
     1                      I2,J2,K2,A2,Amat2,
     2                      L1,M1,N1,B1,Bmat1,
     3                      L2,M2,N2,B2,Bmat2,
     4                      gamma1,gamma2,xggs)

C Test xggs integral using analytical Persson + Taylor
C expressions (uses s-p-d-functions)
C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2

      double precision A1
      double precision A2
      double precision B1
      double precision B2

c     double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2
c     double precision Cmat(3)
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Bmat1(3)
      double precision Bmat2(3)

C Variables returned
      double precision xggs

C Local variables
      integer t1,u1,v1
      integer t2,u2,v2
      integer Nt,Nu,Nv

      double precision c1,c2,c3
      double precision dgx,dgy,dgz
      double precision ans

      double precision one,two,three,pi,half
      parameter(two=2.0d+00,three=3.0d+00,pi=3.14159265358979d+00)
      parameter(half=0.5d+00,one=1.0d+00)

      double precision KAB1(3)
      double precision KAB2(3)

      double precision qAB1
      double precision qAB2

      double precision P1
      double precision P2
      double precision qP1P2

      double precision Pmat1(3)
      double precision Pmat2(3)

      double precision Qmat1(3)
      double precision Qmat2(3)

c     double precision Et1(0:I1,0:L1,0:I1+L1)
c     double precision Eu1(0:J1,0:M1,0:J1+M1)
c     double precision Ev1(0:K1,0:N1,0:K1+N1)

c     double precision Et2(0:I2,0:L2,0:I2+L2)
c     double precision Eu2(0:J2,0:M2,0:J2+M2)
c     double precision Ev2(0:K2,0:N2,0:K2+N2)

      double precision Et1(0:I1+L1)
      double precision Eu1(0:J1+M1)
      double precision Ev1(0:K1+N1)

      double precision Et2(0:I2+L2)
      double precision Eu2(0:J2+M2)
      double precision Ev2(0:K2+N2)

      double precision GAM


      GAM=gamma1+gamma2

C Calculate overlap distributions:  
C Gaussian product of exponentials:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
C Hermite expansion coef.
      call ghec(I1,L1,KAB1(1),A1,B1,P1,qAB1,Qmat1(1),Et1)
      call ghec(J1,M1,KAB1(2),A1,B1,P1,qAB1,Qmat1(2),Eu1)
      call ghec(K1,N1,KAB1(3),A1,B1,P1,qAB1,Qmat1(3),Ev1)

      call ghec(I2,L2,KAB2(1),A2,B2,P2,qAB2,Qmat2(1),Et2)
      call ghec(J2,M2,KAB2(2),A2,B2,P2,qAB2,Qmat2(2),Eu2)
      call ghec(K2,N2,KAB2(3),A2,B2,P2,qAB2,Qmat2(3),Ev2)

c     write(*,*)
c     write(*,*)'Et1',Et1 
c     write(*,*)'Eu1',Eu1  
c     write(*,*)'Ev1',Ev1  
c     write(*,*)'Et2',Et2  
c     write(*,*)'Eu2',Eu2  
c     write(*,*)'Ev2',Ev2  
C Evaluate integral over hermite gaussians

      qP1P2=(p1*p2*GAM)/(p1*p2+p1*GAM+p2*GAM)
      c1=(pi*pi/(P1*P2+P1*GAM+P2*GAM))**(three/two)

      xggs=0.0d+00

      do t1=0,I1+L1
      do u1=0,J1+M1
      do v1=0,K1+N1

      do t2=0,I2+L2
      do u2=0,J2+M2
      do v2=0,K2+N2

         c2=dble((-1)**(t2+u2+v2))

         Nt=t1+t2
         Nu=u1+u2
         Nv=v1+v2

c        call gdel_x(N,x,alp,Ax,gd0)
         call gdel_x(Nt,Pmat1(1),qP1P2,Pmat2(1),dgx)
         call gdel_x(Nu,Pmat1(2),qP1P2,Pmat2(2),dgy)
         call gdel_x(Nv,Pmat1(3),qP1P2,Pmat2(3),dgz)

         ans=c1*c2*
     x       Et1(t1)*Et2(t2) *
     x       Eu1(u1)*Eu2(u2) *
     x       Ev1(v1)*Ev2(v2) *
     x       dgx*dgy*dgz

         xggs=xggs+ans
         
c            write(*,*)
c            write(*,*)'t1=',t1
c            write(*,*)'u1=',u1
c            write(*,*)'v1=',v1
c            write(*,*)
c            write(*,*)'t2=',t2
c            write(*,*)'u2=',u2
c            write(*,*)'v2=',v2
c            write(*,*)
c            write(*,*)'KAB1(1)=',KAB1(1)
c            write(*,*)'KAB1(2)=',KAB1(2)
c            write(*,*)'KAB1(3)=',KAB1(3)
c            write(*,*)
c            write(*,*)'KAB2(1)=',KAB2(1)
c            write(*,*)'KAB2(2)=',KAB2(2)
c            write(*,*)'KAB2(3)=',KAB2(3)
c            write(*,*)
c            write(*,*)'dgx=',dgx
c            write(*,*)'dgy=',dgy
c            write(*,*)'dgz=',dgz
c            write(*,*)
c            write(*,*)'c1=',c1
c            write(*,*)'c2=',c2
c            write(*,*)
c            write(*,*)'Et1=',Et1(t1)
c            write(*,*)'Eu1=',Eu1(u1)
c            write(*,*)'Ev1=',Ev1(v1)
c            write(*,*)
c            write(*,*)'Et2=',Et2(t2)
c            write(*,*)'Eu2=',Eu2(u2)
c            write(*,*)'Ev2=',Ev2(v2)
c            write(*,*)

      end do
      end do
      end do

      end do
      end do
      end do



      return
      end


C======================================================================
c     subroutine cws_gam1_xggvec(I1,J1,K1,A1,Amat1,
      subroutine G2_MD_xggvec(I1,J1,K1,A1,Amat1,
     1                        I2,J2,K2,A2,Amat2,
     2                        L1,M1,N1,B1,Bmat1,
     3                        L2,M2,N2,B2,Bmat2,
     4                        gamma1,gamma2,Cmat,ZNUC,
     5                        xggvec)

C Test xggvec integral using analytical Persson + Taylor
C expressions (uses s-p-d-functions)
C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2

      double precision A1
      double precision A2
      double precision B1
      double precision B2

      double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2
      double precision Cmat(3)
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Bmat1(3)
      double precision Bmat2(3)

C Variables returned
      double precision xggvec

C Local variables
      integer t1,u1,v1
      integer t2,u2,v2
      integer Nt,Nu,Nv

      double precision c1
      double precision c2
      double precision DC1
      double precision DC2
      double precision GAM
      double precision ans
      double precision Dgf

      double precision one,two,three,pi,half
      parameter(two=2.0d+00,three=3.0d+00,pi=3.14159265358979d+00)
      parameter(half=0.5d+00,one=1.0d+00)

      double precision KAB1(3)
      double precision KAB2(3)

      double precision qAB1
      double precision qAB2

      double precision P1
      double precision P2
      double precision qP1P2

      double precision RCSx
      double precision RCSy
      double precision RCSz
      double precision alpha
      double precision expT

      double precision Pmat1(3)
      double precision Pmat2(3)

      double precision Qmat1(3)
      double precision Qmat2(3)

      double precision Et1(0:I1+L1)
      double precision Eu1(0:J1+M1)
      double precision Ev1(0:K1+N1)

      double precision Et2(0:I2+L2)
      double precision Eu2(0:J2+M2)
      double precision Ev2(0:K2+N2)


      GAM=gamma1+gamma2

C Calculate overlap distributions:  
C Gaussian product of exponentials:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
C Hermite expansion coef.
      call ghec(I1,L1,KAB1(1),A1,B1,P1,qAB1,Qmat1(1),Et1)
      call ghec(J1,M1,KAB1(2),A1,B1,P1,qAB1,Qmat1(2),Eu1)
      call ghec(K1,N1,KAB1(3),A1,B1,P1,qAB1,Qmat1(3),Ev1)

      call ghec(I2,L2,KAB2(1),A2,B2,P2,qAB2,Qmat2(1),Et2)
      call ghec(J2,M2,KAB2(2),A2,B2,P2,qAB2,Qmat2(2),Eu2)
      call ghec(K2,N2,KAB2(3),A2,B2,P2,qAB2,Qmat2(3),Ev2)

C Get prefactors for calculating integral
c     call gam1_vec_prefactors(P1,P2,Pmat1,Pmat2,Cmat,GAM,
      call g2_vec_prefactors(P1,P2,Pmat1,Pmat2,Cmat,GAM,
     x                       c1,qP1P2,alpha,expT,
     x                       DC1,DC2,
     x                       RCSx,RCSy,RCSz)

c     write(*,*)
c     write(*,*)'C1=',c1
c     write(*,*)'qP1P2=',qP1P2
c     write(*,*)'alpha=',alpha
c     write(*,*)'expT=',expT
c     write(*,*)'RCSx=',RCSx
c     write(*,*)'RCSy=',RCSy
c     write(*,*)'RCSz=',RCSz
c     write(*,*)
c     write(*,*)'KAB1=',KAB1
c     write(*,*)'KAB2=',KAB2
c     write(*,*)
c     write(*,*)'Et1',Et1 
c     write(*,*)'Eu1',Eu1  
c     write(*,*)'Ev1',Ev1  
c     write(*,*)'Et2',Et2  
c     write(*,*)'Eu2',Eu2  
c     write(*,*)'Ev2',Ev2  
C Evaluate integral over hermite gaussians

c     write(*,*)
c     write(*,*)'I1+L1=',I1+L1
c     write(*,*)'J1+M1=',J1+M1
c     write(*,*)'K1+N1=',K1+N1
c     write(*,*)'I2+L2=',I2+L2
c     write(*,*)'J2+M2=',J2+M2
c     write(*,*)'K2+N2=',K2+N2



      xggvec=0.0d+00

      do t1=0,I1+L1
        do u1=0,J1+M1
          do v1=0,K1+N1

            do t2=0,I2+L2
              do u2=0,J2+M2
                do v2=0,K2+N2

c                  c2=(-one)**dble(t2+u2+v2)

c                  Nt=t1+t2
c                  Nu=u1+u2
c                  Nv=v1+v2


c                  call gam1_dhermite_vec(t1,u1,v1,t2,u2,v2,
                   call g2_dhermite_vec(t1,u1,v1,t2,u2,v2,
     x                                  RCSx,RCSy,RCSz,
     x                                  DC1,DC2,
     x                                  qP1P2,expT,alpha,
     x                                  Pmat1,Pmat2,Dgf)


                   ans=c1*
     x                 Et1(t1)*Et2(t2) *
     x                 Eu1(u1)*Eu2(u2) *
     x                 Ev1(v1)*Ev2(v2) *
     x                 Dgf

c     write(*,*)'Et1',Et1(t1) 
c     write(*,*)'Eu1',Eu1(u1)
c     write(*,*)'Ev1',Ev1(v1)
c     write(*,*)'Et2',Et2(t2)
c     write(*,*)'Eu2',Eu2(u2)
c     write(*,*)'Ev2',Ev2(v2)
c     write(*,*)'c1=',c1
c     write(*,*)'Dgf=',Dgf
                       xggvec=xggvec+ans
         

                end do
              end do
            end do

          end do
        end do
      end do



      return
      end


C=======================================================================
c     subroutine gam1_vec_prefactors(P1,P2,Pmat1,Pmat2,Cmat,GAM,
      subroutine g2_vec_prefactors(P1,P2,Pmat1,Pmat2,Cmat,GAM,
     x                             c1,qP1P2,alpha,expT,
     x                             DC1,DC2,
     x                             RCSx,RCSy,RCSz)
C
C=======================================================================
      implicit none

C Input variables
      double precision P1
      double precision P2
      double precision GAM
      double precision Pmat1(3)
      double precision Pmat2(3)
      double precision Cmat(3)

C Variables returned
      double precision c1
      double precision qP1P2
      double precision alpha
      double precision expT
      double precision RCSx
      double precision RCSy
      double precision RCSz
      double precision DC1
      double precision DC2

C Local variables
      double precision s
      double precision Smat(3)
      double precision xx
      double precision R2CS
      double precision pi
      double precision five
      double precision two
      double precision three
      parameter(pi=3.14159265358979d+00,two=2.0d+00,five=5.0d+00)
      parameter(three=3.0d+00)


      s = ( p1*p2 + p1*gam + p2*gam ) / ( p2+gam )

      alpha=s

      c1 = ( pi / ( p2+gam ) )**(three/two) * two*pi/s

      xx= p1*p2 + p1*gam + p2*gam

      Smat(1)= ( ( p1*p2+p1*gam ) * Pmat1(1) + p2*gam*Pmat2(1) ) / xx
      Smat(2)= ( ( p1*p2+p1*gam ) * Pmat1(2) + p2*gam*Pmat2(2) ) / xx
      Smat(3)= ( ( p1*p2+p1*gam ) * Pmat1(3) + p2*gam*Pmat2(3) ) / xx

      RCSx= Cmat(1) - Smat(1) 
      RCSy= Cmat(2) - Smat(2) 
      RCSz= Cmat(3) - Smat(3) 

      R2CS=RCSx**2 + RCSy**2 + RCSz**2

      expT = alpha*R2CS
      
      qP1P2 = ( p1*p2*GAM) / xx

C Differentiation constants:
      DC1 = ( p1*p2 + p1*gam ) / xx
      DC2 = ( p2*gam ) / xx

      

c     write(*,*)
c     write(*,*)'xxx'
c     write(*,*)'p1=',p1
c     write(*,*)'p2=',p2
c     write(*,*)'gam=',gam
c     write(*,*)



      return
      end

C=======================================================================
c     subroutine gam1_dhermite_vec(t1,u1,v1,t2,u2,v2,
      subroutine g2_dhermite_vec(t1,u1,v1,t2,u2,v2,
     x                           RCSx,RCSy,RCSz,
     x                           DC1,DC2,
     x                           qP1P2,expT,alpha,
     x                           Pmat1,Pmat2,Dgf)
C
C=======================================================================
         implicit none

C Input variables
         integer t1,u1,v1
         integer t2,u2,v2

         double precision Pmat1(3)
         double precision Pmat2(3)
         double precision qP1P2
         double precision expT
         double precision alpha
         double precision RCSx
         double precision RCSy
         double precision RCSz
         double precision DC1
         double precision DC2

C Variables returned
         double precision Dgf

C Local variables
         integer xb1,yb1,zb1
         integer xb2,yb2,zb2
         integer ax1,ay1,az1
         integer ax2,ay2,az2
         integer nLim

         double precision A1,A2,A3
         double precision B1,B2,B3
         double precision C1,C2,C3
         double precision dgx,dgy,dgz
         double precision dgxyz
         double precision dfxyz
         double precision zero
         double precision one
         parameter(zero=0.0d+00)
         parameter(one=1.0d+00)

c        write(*,*)
c        write(*,*)'t1=',t1
c        write(*,*)'u1=',u1
c        write(*,*)'v1=',v1
c        write(*,*)'t2=',t2
c        write(*,*)'u2=',u2
c        write(*,*)'v2=',v2

         Dgf=zero

         do xb1=0,t1
            ax1=t1-xb1
            call cbinom(t1,xb1,A1)

            do yb1=0,u1
               ay1=u1-yb1
               call cbinom(u1,yb1,A2)

               do zb1=0,v1
                  az1=v1-zb1
                  call cbinom(v1,zb1,A3)

                  do xb2=0,t2
                     ax2=t2-xb2
                     call cbinom(t2,xb2,B1)

                     do yb2=0,u2
                        ay2=u2-yb2
                        call cbinom(u2,yb2,B2)

                        do zb2=0,v2
                           az2=v2-zb2
                           call cbinom(v2,zb2,B3)

                       C1=(-one)**(ax2+ay2+az2)
                       C2=(-DC1)**(xb1+yb1+zb1)
                       C3=(-DC2)**(xb2+yb2+zb2)

                       call gdel_x(ax1+ax2,Pmat1(1),qP1P2,Pmat2(1),dgx)
                       call gdel_x(ay1+ay2,Pmat1(2),qP1P2,Pmat2(2),dgy)
                       call gdel_x(az1+az2,Pmat1(3),qP1P2,Pmat2(3),dgz)
                       dgxyz=dgx*dgy*dgz

c                          call gdel_x(ax1,Pmat1(1),qP1P2,Pmat2(1),dgx1)
c                          call gdel_x(ay1,Pmat1(2),qP1P2,Pmat2(2),dgy1)
c                          call gdel_x(az1,Pmat1(3),qP1P2,Pmat2(3),dgz1)

c                          call gdel_x(ax2,Pmat2(1),qP1P2,Pmat1(1),dgx2)
c                          call gdel_x(ay2,Pmat2(2),qP1P2,Pmat1(2),dgy2)
c                          call gdel_x(az2,Pmat2(3),qP1P2,Pmat1(3),dgz2)

c                          dgxyz=dgx1*dgy1*dgz1*dgx2*dgy2*dgz2

c                          call gdel(argx,argy,argz,qP1P2,Pmat1,Pmat2,dgxyz)

                       nLim=xb1+yb1+zb1+xb2+yb2+zb2
                       call RTUV(xb1+xb2,yb1+yb2,zb1+zb2,nLim,
     x                           expT,alpha,RCSx,RCSy,RCSz,dfxyz)
c                      nLim=xb2+yb2+zb2
c                      call RTUV(xb2,yb2,zb2,nLim,
c    x                           expT,alpha,RCSx,RCSy,RCSz,dfxyz1)

c                      write(*,*)
c                      write(*,*)'A1=',A1
c                      write(*,*)'A2=',A2
c                      write(*,*)'A3=',A3
c                      write(*,*)'B1=',B1
c                      write(*,*)'B2=',B2
c                      write(*,*)'B3=',B3
c                      write(*,*)'C1=',C1
c                      write(*,*)'C2=',C2
c                      write(*,*)'C3=',C3
c                      write(*,*)'DC1=',DC1
c                      write(*,*)'DC2=',DC2
c                      write(*,*)'dgxyz=',dgxyz
c                      write(*,*)'dfxyz=',dfxyz
c                      write(*,*)
                       Dgf=Dgf+A1*A2*A3*
     x                         B1*B2*B3*
     x                         C1*C2*C3*
     x                         dgxyz*dfxyz


                        end do 
                     end do 
                  end do
               end do 
            end do 
         end do


c        write(*,*)
c        write(*,*)'argx=',argx
c        write(*,*)'argy=',argy
c        write(*,*)'argz=',argz
c        write(*,*)'qP1P2=',qP1P2
c        write(*,*)'bin_tx1=',bin_tx1
c        write(*,*)'bin_uy1=',bin_uy1
c        write(*,*)'bin_vz1=',bin_vz1
c        write(*,*)'bin_tx2=',bin_tx2
c        write(*,*)'bin_uy2=',bin_uy2
c        write(*,*)'bin_vz2=',bin_vz2
c        write(*,*)'Dg_sign=',Dg_sign
c        write(*,*)'Df_sign=',Df_sign
c        write(*,*)'dgx=',dgx
c        write(*,*)'dgy=',dgy
c        write(*,*)'dgz=',dgz
c        write(*,*)'dgxyz=',dgxyz
c        write(*,*)'dfxyz=',dfxyz
c        write(*,*)



      return
      end

CCCCCCCCC

C======================================================================
c     subroutine cws_gam1_xggvee(I1,J1,K1,A1,Amat1,
      subroutine G2_MD_xggvee(I1,J1,K1,A1,Amat1,
     1                        I2,J2,K2,A2,Amat2,
     2                        L1,M1,N1,B1,Bmat1,
     3                        L2,M2,N2,B2,Bmat2,
     4                        gamma1,gamma2,xggvee)

C Test xggvep integral using analytical Persson + Taylor
C expressions (uses s-p-d-functions)
C======================================================================
      implicit none

C Input variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2

      double precision A1
      double precision A2
      double precision B1
      double precision B2

c     double precision ZNUC
      double precision GAMMA1
      double precision GAMMA2
c     double precision Cmat(3)
      double precision Amat1(3)
      double precision Amat2(3)
      double precision Bmat1(3)
      double precision Bmat2(3)

C Variables returned
      double precision xggvee

C Local variables
      integer t1,u1,v1
      integer t2,u2,v2
      integer Nt,Nu,Nv

      double precision c1
      double precision c2
      double precision GAM
      double precision ans
      double precision Dgf

      double precision one,two,three,pi,half
      parameter(two=2.0d+00,three=3.0d+00,pi=3.14159265358979d+00)
      parameter(half=0.5d+00,one=1.0d+00)

      double precision KAB1(3)
      double precision KAB2(3)

      double precision qAB1
      double precision qAB2

      double precision P1
      double precision P2
      double precision qP1P2

      double precision X_P1P2
      double precision Y_P1P2
      double precision Z_P1P2
      double precision alpha
      double precision expT

      double precision Pmat1(3)
      double precision Pmat2(3)

      double precision Qmat1(3)
      double precision Qmat2(3)

      double precision Et1(0:I1+L1)
      double precision Eu1(0:J1+M1)
      double precision Ev1(0:K1+N1)

      double precision Et2(0:I2+L2)
      double precision Eu2(0:J2+M2)
      double precision Ev2(0:K2+N2)


      GAM=gamma1+gamma2

C Calculate overlap distributions:  
C Gaussian product of exponentials:
      call Gauss_prod(A1,B1,Amat1,Bmat1,KAB1,qAB1,P1,Pmat1,Qmat1)
      call Gauss_prod(A2,B2,Amat2,Bmat2,KAB2,qAB2,P2,Pmat2,Qmat2)
C Hermite expansion coef.
      call ghec(I1,L1,KAB1(1),A1,B1,P1,qAB1,Qmat1(1),Et1)
      call ghec(J1,M1,KAB1(2),A1,B1,P1,qAB1,Qmat1(2),Eu1)
      call ghec(K1,N1,KAB1(3),A1,B1,P1,qAB1,Qmat1(3),Ev1)

      call ghec(I2,L2,KAB2(1),A2,B2,P2,qAB2,Qmat2(1),Et2)
      call ghec(J2,M2,KAB2(2),A2,B2,P2,qAB2,Qmat2(2),Eu2)
      call ghec(K2,N2,KAB2(3),A2,B2,P2,qAB2,Qmat2(3),Ev2)

C Get prefactors for calculating integral
c     call gam1_ee_prefactors(P1,P2,Pmat1,Pmat2,GAM,
      call g2_ee_prefactors(P1,P2,Pmat1,Pmat2,GAM,
     x                      c1,qP1P2,alpha,expT,
     x                      X_P1P2,Y_P1P2,Z_P1P2)

c     write(*,*)
c     write(*,*)'C1=',c1
c     write(*,*)'qP1P2=',qP1P2
c     write(*,*)'alpha=',alpha
c     write(*,*)'expT=',expT
c     write(*,*)'X_P1P2=',X_P1P2
c     write(*,*)'Y_P1P2=',Y_P1P2
c     write(*,*)'Z_P1P2=',Z_P1P2
c     write(*,*)
c     write(*,*)'KAB1=',KAB1
c     write(*,*)'KAB2=',KAB2
c     write(*,*)
c     write(*,*)'Et1',Et1 
c     write(*,*)'Eu1',Eu1  
c     write(*,*)'Ev1',Ev1  
c     write(*,*)'Et2',Et2  
c     write(*,*)'Eu2',Eu2  
c     write(*,*)'Ev2',Ev2  
C Evaluate integral over hermite gaussians


      xggvee=0.0d+00

      do t1=0,I1+L1
      do u1=0,J1+M1
      do v1=0,K1+N1

      do t2=0,I2+L2
      do u2=0,J2+M2
      do v2=0,K2+N2

         c2=dble((-1)**(t2+u2+v2))

         Nt=t1+t2
         Nu=u1+u2
         Nv=v1+v2


c        call gam1_dhermite_vee(Nt,Nu,Nv,
         call g2_dhermite_vee(Nt,Nu,Nv,
     x                        X_P1P2,Y_P1P2,Z_P1P2,
     x                        qP1P2,expT,alpha,
     x                        Pmat1,Pmat2,Dgf)


         ans=c1*c2*
     x   Et1(t1)*Et2(t2) *
     x   Eu1(u1)*Eu2(u2) *
     x   Ev1(v1)*Ev2(v2) *
     x   Dgf

         xggvee=xggvee+ans
         

      end do
      end do
      end do

      end do
      end do
      end do



      return
      end


C=======================================================================
c     subroutine gam1_ee_prefactors(P1,P2,Pmat1,Pmat2,GAM,
      subroutine g2_ee_prefactors(P1,P2,Pmat1,Pmat2,GAM,
     x                            c0,qP1P2,alpha,expT,
     x                            X_P1P2,Y_P1P2,Z_P1P2)
C
C=======================================================================
      implicit none

C Input variables
      double precision P1
      double precision P2
      double precision GAM
      double precision Pmat1(3)
      double precision Pmat2(3)

C Variables returned
      double precision C0
      double precision qP1P2
      double precision alpha
      double precision expT
      double precision X_P1P2
      double precision Y_P1P2
      double precision Z_P1P2

C Local variables
      double precision xx
      double precision R2P1P2
      double precision pi
      double precision five
      double precision two
      parameter(pi=3.14159265358979d+00,two=2.0d+00,five=5.0d+00)


c     write(*,*)
c     write(*,*)'xxx'
c     write(*,*)'p1=',p1
c     write(*,*)'p2=',p2
c     write(*,*)'gam=',gam
c     write(*,*)

      xx = p1*p2+p1*gam+p2*gam

      c0 = ( two * (pi)**(five/two) ) / ( xx * sqrt(p1+p2) )

c     qP1P2 = ( p1*p2*GAM) / xx
      qP1P2 = (p1*p2*GAM)/(p1*p2+p1*GAM+p2*GAM) 

      alpha = ( p1*p1*p2*p2 ) / ( xx * (p1 + p2) )

      X_P1P2=Pmat1(1)-Pmat2(1)
      Y_P1P2=Pmat1(2)-Pmat2(2)
      Z_P1P2=Pmat1(3)-Pmat2(3)

      R2P1P2=X_P1P2**2 + Y_P1P2**2 + Z_P1P2**2

      expT = alpha*R2P1P2


      return
      end

C=======================================================================
c     subroutine gam1_dhermite_vee(Nt,Nu,Nv,
      subroutine g2_dhermite_vee(Nt,Nu,Nv,
     x                           X_P1P2,Y_P1P2,Z_P1P2,
     x                           qP1P2,expT,alpha,
     x                           Pmat1,Pmat2,Dgf)
C
C=======================================================================
         implicit none

C Input variables
         integer Nt,Nu,Nv

         double precision Pmat1(3)
         double precision Pmat2(3)
         double precision qP1P2
         double precision expT
         double precision alpha
         double precision X_P1P2
         double precision Y_P1P2
         double precision Z_P1P2

C Variables returned
         double precision Dgf

C Local variables
         integer xb1,yb1,zb1
         integer nLim
         integer argx,argy,argz

         double precision bin_tx1
         double precision bin_uy1
         double precision bin_vz1
         double precision dgx,dgy,dgz
         double precision dgxyz
         double precision dfxyz
         double precision zero
         double precision one
         parameter(zero=0.0d+00)
         parameter(one=1.0d+00)


         Dgf=zero

         do xb1=0,Nt
            argx=Nt-xb1
            call cbinom(Nt,xb1,bin_tx1)

         do yb1=0,Nu
            argy=Nu-yb1
            call cbinom(Nu,yb1,bin_uy1)

         do zb1=0,Nv
            argz=Nv-zb1
            call cbinom(Nv,zb1,bin_vz1)

            call gdel_x(argx,Pmat1(1),qP1P2,Pmat2(1),dgx)
            call gdel_x(argy,Pmat1(2),qP1P2,Pmat2(2),dgy)
            call gdel_x(argz,Pmat1(3),qP1P2,Pmat2(3),dgz)
            dgxyz=dgx*dgy*dgz
c           call gdel(argx,argy,argz,qP1P2,Pmat1,Pmat2,dgxyz)

            nLim=xb1+yb1+zb1

            call RTUV(xb1,yb1,zb1,nLim,
     x                expT,alpha,X_P1P2,Y_P1P2,Z_P1P2,dfxyz)


            Dgf=Dgf+bin_tx1*bin_uy1*bin_vz1
     x             *dgxyz*dfxyz


c        write(*,*)
c        write(*,*)'argx=',argx
c        write(*,*)'argy=',argy
c        write(*,*)'argz=',argz
c        write(*,*)'qP1P2=',qP1P2
c        write(*,*)'bin_tx1=',bin_tx1
c        write(*,*)'bin_uy1=',bin_uy1
c        write(*,*)'bin_vz1=',bin_vz1
c        write(*,*)'bin_tx2=',bin_tx2
c        write(*,*)'bin_uy2=',bin_uy2
c        write(*,*)'bin_vz2=',bin_vz2
c        write(*,*)'Dg_sign=',Dg_sign
c        write(*,*)'Df_sign=',Df_sign
c        write(*,*)'dgx=',dgx
c        write(*,*)'dgy=',dgy
c        write(*,*)'dgz=',dgz
c        write(*,*)'dgxyz=',dgxyz
c        write(*,*)'dfxyz=',dfxyz
c        write(*,*)



         end do
         end do
         end do

      return
      end

CCCCCCCCC
