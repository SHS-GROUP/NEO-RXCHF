C G3VeC_AUX_g13g23V1.f
C G3VpC_AUX_g13g23V3.f 
C G3Vep_AUX_g13g23V13.f
C G3Vee_AUX_g13V12.f  
C G3Vee_AUX_g23V12.f   

C=======================================================================
      subroutine G3VeC_AUX_g13g23V1(I1,J1,K1,A1,Amat1,
     *                              I2,J2,K2,A2,Amat2,
     *                              I3,J3,K3,A3,Amat3,
     *                              L1,M1,N1,B1,Bmat1,
     *                              L2,M2,N2,B2,Bmat2,
     *                              L3,M3,N3,B3,Bmat3,
     *                              gamA13,gamA23,
     *                              gamB13,gamB23,
     *                              Cmat,ZNUC,xgVECg)

C Evaluates the following VeC integrals:
C xgVEC    ::  gA(2,3) V(1,Rc)
C xgVECg1  ::  gA(1,3) V(1,Rc) gB(2,3)
C xgVECg2  ::  gA(2,3) V(1,Rc) gB(2,3)
C
C calls G2_MD_xggvec to evaluate auxilliary integrals
C
C Note:  Proton is assumed to be particle 3
C=======================================================================
      implicit none

C Input Variables
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
      double precision gamA13
      double precision gamA23
      double precision gamB13
      double precision gamB23
      double precision Cmat(3)
      double precision ZNUC

C Variables Returned
      double precision xgVECg

C Local Variables
      integer t1,t2,xb 
      integer u1,u2,yb 
      integer v1,v2,zb 
      double precision bin_t1,bin_t2,bin_xb
      double precision bin_u1,bin_u2,bin_yb
      double precision bin_v1,bin_v2,bin_zb
      double precision gam13
      double precision gam23
      double precision p2
      double precision Pmat2(3)
      double precision KAB2
      double precision XPA2,YPA2,ZPA2
      double precision XPB2,YPB2,ZPB2
      double precision C23
      double precision aC23
      double precision a
      double precision CHI
      double precision ans
      double precision G2VeC
      double precision G3VeC



      gam13=gamA13+gamB13
      gam23=gamA23+gamB23


C  Get basis function contraction values for particle 2:
      call contract_BF(A2,B2,Amat2,Bmat2,p2,Pmat2,KAB2,
     x                 XPA2,YPA2,ZPA2,XPB2,YPB2,ZPB2)

      C23=gam23/(p2+gam23)
      aC23=p2*C23

      a=p2+gam23

      G3VeC=0.0d+00

      do t1=0,I2
       call cbinom(I2,t1,bin_t1)
       do t2=0,L2
        call cbinom(L2,t2,bin_t2)
        do xb=0,t1+t2
         call cbinom(t1+t2,xb,bin_xb)

         do u1=0,J2
          call cbinom(J2,u1,bin_u1)
          do u2=0,M2
           call cbinom(M2,u2,bin_u2)
           do yb=0,u1+u2
            call cbinom(u1+u2,yb,bin_yb)

            do v1=0,K2
             call cbinom(K2,v1,bin_v1)
             do v2=0,N2
              call cbinom(N2,v2,bin_v2)
              do zb=0,v1+v2
               call cbinom(v1+v2,zb,bin_zb)

               call evaluate_CHI(xb,yb,zb,a,CHI)

C Calculate the "auxiliary" 2-particle integral 
               call aux_G2VeC(I1,J1,K1,A1,Amat1,
     x                        I3,J3,K3,A3,Amat3,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L3,M3,N3,B3,Bmat3,
     x                        aC23,Pmat2,
     x                        t1,t2,xb,
     x                        u1,u2,yb,
     x                        v1,v2,zb,
     x                        gam13,Cmat,ZNUC,G2VeC)

               
               ans=bin_t1*bin_t2*bin_xb
     x            *bin_u1*bin_u2*bin_yb
     x            *bin_v1*bin_v2*bin_zb
     x            *XPA2**(I2-t1) * YPA2**(J2-u1) * ZPA2**(K2-v1)
     x            *XPB2**(L2-t2) * YPB2**(M2-u2) * ZPB2**(N2-v2)
     x            *C23**(t1+t2-xb)
     x            *C23**(u1+u2-yb)
     x            *C23**(v1+v2-zb)
     x            *CHI*G2VeC


               G3VeC=G3VeC+ans

              end do
             end do
            end do

           end do
          end do
         end do

        end do
       end do
      end do

      G3VeC=G3VeC*KAB2
c     xgVECg=-1.0d+00*G3VeC
      xgVECg=G3VeC


      return
      end

C=======================================================================
      subroutine aux_G2VeC(I1,J1,K1,A1,Amat1,
     x                     I3,J3,K3,A3,Amat3,
     x                     L1,M1,N1,B1,Bmat1,
     x                     L3,M3,N3,B3,Bmat3,
     x                     aC23,Pmat2,
     x                     t1,t2,xb,
     x                     u1,u2,yb,
     x                     v1,v2,zb,
     x                     gam13,Cmat,ZNUC,G2VeC)

C=======================================================================
      implicit none

C Input Variables
      integer I1,J1,K1
      integer I3,J3,K3
      integer L1,M1,N1
      integer L3,M3,N3
      integer t1,t2,xb
      integer u1,u2,yb
      integer v1,v2,zb
       
      double precision A1,Amat1(3)
      double precision A3,Amat3(3)
      double precision B1,Bmat1(3)
      double precision B3,Bmat3(3)
      double precision aC23,Pmat2(3)
      double precision gam13
      double precision Cmat(3)
      double precision ZNUC 

C Variables Returned
      double precision G2VeC

C Local Variables
      integer t3,u3,v3
      integer t4,u4,v4
      integer IX3,JX3,KX3
      integer LX3,MX3,NX3

      double precision bin_t3,bin_u3,bin_v3
      double precision bin_t4,bin_u4,bin_v4
      double precision P3
      double precision Pmat3(3)
      double precision KAB3
      double precision XPA3,YPA3,ZPA3
      double precision XPB3,YPB3,ZPB3
      double precision ans
      double precision xG2VeC
      double precision zero
      parameter(zero=0.0d+00)

C  Get basis function contraction values for particle 3:
      call contract_BF(A3,B3,Amat3,Bmat3,p3,Pmat3,KAB3,
     x                 XPA3,YPA3,ZPA3,XPB3,YPB3,ZPB3)

      IX3=t1+t2-xb
      JX3=u1+u2-yb
      KX3=v1+v2-zb

      G2VeC=zero

      do t3=0,I3
       call cbinom(I3,t3,bin_t3)
       do t4=0,L3
        call cbinom(L3,t4,bin_t4)
        do u3=0,J3
         call cbinom(J3,u3,bin_u3)
         do u4=0,M3
          call cbinom(M3,u4,bin_u4)
          do v3=0,K3
           call cbinom(K3,v3,bin_v3)
           do v4=0,N3
            call cbinom(N3,v4,bin_v4)


            LX3=t3+t4
            MX3=u3+u4
            NX3=v3+v4



c     write(*,*)'I1=',I1
c     write(*,*)'J1=',J1
c     write(*,*)'K1=',K1
c     write(*,*)'A1=',A1
c     write(*,*)'Amat1=',Amat1
c     write(*,*)'IX3=',IX3
c     write(*,*)'JX3=',JX3
c     write(*,*)'KX3=',KX3
c     write(*,*)'ac23=',ac23
c     write(*,*)'Pmat2=',Pmat2
c     write(*,*)'L1=',L1
c     write(*,*)'M1=',M1
c     write(*,*)'N1=',N1
c     write(*,*)'B1=',B1
c     write(*,*)'Bmat1=',Bmat1
c     write(*,*)'LX3=',LX3
c     write(*,*)'MX3=',MX3
c     write(*,*)'NX3=',NX3
c     write(*,*)'p3=',p3
c     write(*,*)'Pmat3=',Pmat3
c     write(*,*)'Cmat=',Cmat
c           call cws_gam1_xggvec(I1,J1,K1,A1,Amat1,
c           call G2_MD_xggvec(I1,J1,K1,A1,Amat1,
c    x                        IX3,JX3,KX3,aC23,Pmat2,
c    x                        L1,M1,N1,B1,Bmat1,
c    x                        LX3,MX3,NX3,p3,Pmat3,
c    x                        gam13,zero,Cmat,ZNUC,
c    x                        xG2VeC)
            call pgivec(I1,J1,K1,A1,Amat1,
     x                  IX3,JX3,KX3,aC23,Pmat2,
     x                  L1,M1,N1,B1,Bmat1,
     x                  LX3,MX3,NX3,p3,Pmat3,
     x                  gam13,zero,Cmat,xG2VeC)
c           xG2VeC=ZNUC*xG2VeC


            ans=bin_t3*bin_t4
     x         *bin_u3*bin_u4
     x         *bin_v3*bin_v4
     x         *XPA3**(I3-t3)
     x         *XPB3**(L3-t4)
     x         *YPA3**(J3-u3)
     x         *YPB3**(M3-u4)
     x         *ZPA3**(K3-v3)
     x         *ZPB3**(N3-v4)
     x         *xG2VeC

            G2VeC=G2VeC+ans

           end do
          end do
         end do
        end do
       end do
      end do

      G2VeC=G2VeC*KAB3


      return
      end



C=======================================================================
      subroutine G3VpC_AUX_g13g23V3(I1,J1,K1,A1,Amat1,
     *                              I2,J2,K2,A2,Amat2,
     *                              I3,J3,K3,A3,Amat3,
     *                              L1,M1,N1,B1,Bmat1,
     *                              L2,M2,N2,B2,Bmat2,
     *                              L3,M3,N3,B3,Bmat3,
     *                              gamA13,gamA23,
     *                              gamB13,gamB23,
     *                              Cmat,ZNUC,xgVPCg)

C Evaluates the following VpC integrals
C xgVpCg   ::  gA(1,3) V(3,Rc) gB(2,3)
C
C calls G2_MD_xggvec to evaluate auxilliary integrals
C
C Note:  Proton is assumed to be particle 3
C=======================================================================
      implicit none

C Input Variables
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
      double precision gamA13
      double precision gamA23
      double precision gamB13
      double precision gamB23
      double precision Cmat(3)
      double precision ZNUC

C Variables Returned
      double precision xgVPCg

C Local Variables
      integer t1,t2,xb 
      integer u1,u2,yb 
      integer v1,v2,zb 
      double precision bin_t1,bin_t2,bin_xb
      double precision bin_u1,bin_u2,bin_yb
      double precision bin_v1,bin_v2,bin_zb
      double precision gam13
      double precision gam23
      double precision p2
      double precision Pmat2(3)
      double precision KAB2
      double precision XPA2,YPA2,ZPA2
      double precision XPB2,YPB2,ZPB2
      double precision C23
      double precision aC23
      double precision a
      double precision CHI
      double precision ans
      double precision G2VpC
      double precision G3VpC
CCWS-BEDUG
      double precision ans1
      double precision ans2
      double precision ans3
CCWS-BEDUG



      gam13=gamA13+gamB13
      gam23=gamA23+gamB23


C  Get basis function contraction values for particle 2:
      call contract_BF(A2,B2,Amat2,Bmat2,p2,Pmat2,KAB2,
     x                 XPA2,YPA2,ZPA2,XPB2,YPB2,ZPB2)

      C23=gam23/(p2+gam23)
      aC23=p2*C23

      a=p2+gam23

      G3VpC=0.0d+00
      xgVPCg=0.0d+00

CCWS-DEBUG
c           call gfovlap(I1,J1,K1,a1,Amat1,
c    2                   L1,M1,N1,b1,Bmat1,ans1)

c           call gfovlap(I2,J2,K2,a2,Amat2,
c    2                   L2,M2,N2,b2,Bmat2,ans2)

c           call gfvec(I3,J3,K3,a3,Amat3,
c    2                 L3,M3,N3,b3,Bmat3,
c    4                 Cmat,ans3)
c           write(*,*)'SIMPLE ANS=',ans1*ans2*ans3*znuc

CCWS-DEBUG
      do t1=0,I2
       call cbinom(I2,t1,bin_t1)
       do t2=0,L2
        call cbinom(L2,t2,bin_t2)
        do xb=0,t1+t2
         call cbinom(t1+t2,xb,bin_xb)

         do u1=0,J2
          call cbinom(J2,u1,bin_u1)
          do u2=0,M2
           call cbinom(M2,u2,bin_u2)
           do yb=0,u1+u2
            call cbinom(u1+u2,yb,bin_yb)

            do v1=0,K2
             call cbinom(K2,v1,bin_v1)
             do v2=0,N2
              call cbinom(N2,v2,bin_v2)
              do zb=0,v1+v2
               call cbinom(v1+v2,zb,bin_zb)

               call evaluate_CHI(xb,yb,zb,a,CHI)

C Calculate the "auxiliary" 2-particle integral 
               call aux_G2VpC(I1,J1,K1,A1,Amat1,
     x                        I3,J3,K3,A3,Amat3,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L3,M3,N3,B3,Bmat3,
     x                        aC23,Pmat2,
     x                        t1,t2,xb,
     x                        u1,u2,yb,
     x                        v1,v2,zb,
     x                        gam13,Cmat,ZNUC,G2VpC)

               
               ans=bin_t1*bin_t2*bin_xb
     x            *bin_u1*bin_u2*bin_yb
     x            *bin_v1*bin_v2*bin_zb
     x            *XPA2**(I2-t1) * YPA2**(J2-u1) * ZPA2**(K2-v1)
     x            *XPB2**(L2-t2) * YPB2**(M2-u2) * ZPB2**(N2-v2)
     x            *C23**(t1+t2-xb)
     x            *C23**(u1+u2-yb)
     x            *C23**(v1+v2-zb)
     x            *CHI*G2VpC


               G3VpC=G3VpC+ans

              end do
             end do
            end do

           end do
          end do
         end do

        end do
       end do
      end do

      G3VpC=G3VpC*KAB2
      xgVPCg=G3VpC


      return
      end

C=======================================================================
      subroutine aux_G2VpC(I1,J1,K1,A1,Amat1,
     x                     I3,J3,K3,A3,Amat3,
     x                     L1,M1,N1,B1,Bmat1,
     x                     L3,M3,N3,B3,Bmat3,
     x                     aC23,Pmat2,
     x                     t1,t2,xb,
     x                     u1,u2,yb,
     x                     v1,v2,zb,
     x                     gam13,Cmat,ZNUC,G2VpC)

C=======================================================================
      implicit none

C Input Variables
      integer I1,J1,K1
      integer I3,J3,K3
      integer L1,M1,N1
      integer L3,M3,N3
      integer t1,t2,xb
      integer u1,u2,yb
      integer v1,v2,zb
       
      double precision A1,Amat1(3)
      double precision A3,Amat3(3)
      double precision B1,Bmat1(3)
      double precision B3,Bmat3(3)
      double precision aC23,Pmat2(3)
      double precision gam13
      double precision Cmat(3)
      double precision ZNUC 

C Variables Returned
      double precision G2VpC


C Local Variables
      integer t3,u3,v3
      integer t4,u4,v4
      integer IX3,JX3,KX3
      integer LX3,MX3,NX3

      double precision bin_t3,bin_u3,bin_v3
      double precision bin_t4,bin_u4,bin_v4
      double precision P3
      double precision Pmat3(3)
      double precision KAB3
      double precision XPA3,YPA3,ZPA3
      double precision XPB3,YPB3,ZPB3
      double precision ans
CCWS-BEDUG
      double precision ans1
      double precision ans2
CCWS-BEDUG
      double precision xG2VpC
      double precision zero
      parameter(zero=0.0d+00)

C  Get basis function contraction values for particle 3:
      call contract_BF(A3,B3,Amat3,Bmat3,p3,Pmat3,KAB3,
     x                 XPA3,YPA3,ZPA3,XPB3,YPB3,ZPB3)

      IX3=t1+t2-xb
      JX3=u1+u2-yb
      KX3=v1+v2-zb

      G2VpC=zero
      xG2VpC=zero

      do t3=0,I3
       call cbinom(I3,t3,bin_t3)
       do t4=0,L3
        call cbinom(L3,t4,bin_t4)
        do u3=0,J3
         call cbinom(J3,u3,bin_u3)
         do u4=0,M3
          call cbinom(M3,u4,bin_u4)
          do v3=0,K3
           call cbinom(K3,v3,bin_v3)
           do v4=0,N3
            call cbinom(N3,v4,bin_v4)


            LX3=t3+t4
            MX3=u3+u4
            NX3=v3+v4


c           call cws_gam1_xggvec(I1,J1,K1,A1,Amat1,
c    x                           IX3,JX3,KX3,aC23,Pmat2,
c    x                           L1,M1,N1,B1,Bmat1,
c    x                           LX3,MX3,NX3,p3,Pmat3,
c    x                           gam13,zero,Cmat,ZNUC,
c    x                           xG2VeC)

c           call cws_gam1_xggvec(IX3,JX3,KX3,aC23,Pmat2,
c           call G2_MD_xggvec(IX3,JX3,KX3,aC23,Pmat2,
c    x                        I1,J1,K1,A1,Amat1,
c    x                        LX3,MX3,NX3,p3,Pmat3,
c    x                        L1,M1,N1,B1,Bmat1,
c    x                        gam13,zero,Cmat,ZNUC,
c    x                        xG2VpC)
c           xG2VpC=xG2VpC*ZNUC

            call pgivec(IX3,JX3,KX3,aC23,Pmat2,
     x                  I1,J1,K1,A1,Amat1,
     x                  LX3,MX3,NX3,p3,Pmat3,
     x                  L1,M1,N1,B1,Bmat1,
     x                  gam13,zero,Cmat,xG2VpC)
c           xG2VeC=ZNUC*xG2VeC
CCWS-DEBUG
c           call gfovlap(I1,J1,K1,a1,Amat1,
c    2                   L1,M1,N1,b1,Bmat1,ans1)

c           call gfvec(IX3,JX3,KX3,ac23,Pmat2,
c    2                 LX3,MX3,NX3,p3,Pmat3,
c    4                 Cmat,ans2)
c           xG2VpC=ans1*ans2

CCWS-DEBUG


            ans=bin_t3*bin_t4
     x         *bin_u3*bin_u4
     x         *bin_v3*bin_v4
     x         *XPA3**(I3-t3)
     x         *XPB3**(L3-t4)
     x         *YPA3**(J3-u3)
     x         *YPB3**(M3-u4)
     x         *ZPA3**(K3-v3)
     x         *ZPB3**(N3-v4)
     x         *xG2VpC

            G2VpC=G2VpC+ans

           end do
          end do
         end do
        end do
       end do
      end do

      G2VpC=G2VpC*KAB3


      return
      end



C=======================================================================
      subroutine G3Vep_AUX_g13g23V13(I1,J1,K1,A1,Amat1,
     *                               I2,J2,K2,A2,Amat2,
     *                               I3,J3,K3,A3,Amat3,
     *                               L1,M1,N1,B1,Bmat1,
     *                               L2,M2,N2,B2,Bmat2,
     *                               L3,M3,N3,B3,Bmat3,
     *                               gamA13,gamA23,
     *                               gamB13,gamB23,
     *                               xgVEPg)

C Evaluates the following Vep integrals:
C xgVEP    ::  gA(2,3) V(1,3)
C xgVEPg1  ::  gA(1,3) V(1,3) gB(2,3)
C xgVEPg2  ::  gA(2,3) V(1,3) gB(2,3)
C
C calls G2_MD_xggvee to evaluate auxilliary integrals
C
C Note:  Proton is assumed to be particle 3
C=======================================================================
      implicit none

C Input Variables
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
      double precision gamA13
      double precision gamA23
      double precision gamB13
      double precision gamB23

C Variables Returned
      double precision xgVEPg

C Local Variables
      integer t1,t2,xb 
      integer u1,u2,yb 
      integer v1,v2,zb 
      double precision bin_t1,bin_t2,bin_xb
      double precision bin_u1,bin_u2,bin_yb
      double precision bin_v1,bin_v2,bin_zb
      double precision gam13
      double precision gam23
      double precision p2
      double precision Pmat2(3)
      double precision KAB2
      double precision XPA2,YPA2,ZPA2
      double precision XPB2,YPB2,ZPB2
      double precision C23
      double precision aC23
      double precision a
      double precision G3Vep
      double precision CHI
      double precision ans
      double precision G2Vep



      gam13=gamA13+gamB13
      gam23=gamA23+gamB23


C  Get basis function contraction values for particle 2:
      call contract_BF(A2,B2,Amat2,Bmat2,p2,Pmat2,KAB2,
     x                 XPA2,YPA2,ZPA2,XPB2,YPB2,ZPB2)

      C23=gam23/(p2+gam23)
      aC23=p2*C23

      a=p2+gam23

      G3Vep=0.0d+00

      do t1=0,I2
       call cbinom(I2,t1,bin_t1)
       do t2=0,L2
        call cbinom(L2,t2,bin_t2)
        do xb=0,t1+t2
         call cbinom(t1+t2,xb,bin_xb)

         do u1=0,J2
          call cbinom(J2,u1,bin_u1)
          do u2=0,M2
           call cbinom(M2,u2,bin_u2)
           do yb=0,u1+u2
            call cbinom(u1+u2,yb,bin_yb)

            do v1=0,K2
             call cbinom(K2,v1,bin_v1)
             do v2=0,N2
              call cbinom(N2,v2,bin_v2)
              do zb=0,v1+v2
               call cbinom(v1+v2,zb,bin_zb)

               call evaluate_CHI(xb,yb,zb,a,CHI)

C Calculate the "auxiliary" 2-particle integral 
               call aux_G2Vep(I1,J1,K1,A1,Amat1,
     x                        I3,J3,K3,A3,Amat3,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L3,M3,N3,B3,Bmat3,
     x                        aC23,Pmat2,
     x                        t1,t2,xb,
     x                        u1,u2,yb,
     x                        v1,v2,zb,
     x                        gam13,G2Vep)

               
               ans=bin_t1*bin_t2*bin_xb
     x            *bin_u1*bin_u2*bin_yb
     x            *bin_v1*bin_v2*bin_zb
     x            *XPA2**(I2-t1) * YPA2**(J2-u1) * ZPA2**(K2-v1)
     x            *XPB2**(L2-t2) * YPB2**(M2-u2) * ZPB2**(N2-v2)
     x            *C23**(t1+t2-xb)
     x            *C23**(u1+u2-yb)
     x            *C23**(v1+v2-zb)
     x            *CHI*G2Vep




               G3Vep=G3Vep+ans

              end do
             end do
            end do

           end do
          end do
         end do

        end do
       end do
      end do

      G3Vep=G3Vep*KAB2
      xgVEPg=G3Vep


      return
      end

C=======================================================================
      subroutine aux_G2Vep(I1,J1,K1,A1,Amat1,
     x                     I3,J3,K3,A3,Amat3,
     x                     L1,M1,N1,B1,Bmat1,
     x                     L3,M3,N3,B3,Bmat3,
     x                     aC23,Pmat2,
     x                     t1,t2,xb,
     x                     u1,u2,yb,
     x                     v1,v2,zb,
     x                     gam13,G2Vep)

C=======================================================================
      implicit none

C Input Variables
      integer I1,J1,K1
      integer I3,J3,K3
      integer L1,M1,N1
      integer L3,M3,N3
      integer t1,t2,xb
      integer u1,u2,yb
      integer v1,v2,zb
       
      double precision A1,Amat1(3)
      double precision A3,Amat3(3)
      double precision B1,Bmat1(3)
      double precision B3,Bmat3(3)
      double precision aC23,Pmat2(3)
      double precision gam13

C Variables Returned
      double precision G2Vep


C Local Variables
      integer t3,u3,v3
      integer t4,u4,v4
      integer IX3,JX3,KX3
      integer LX3,MX3,NX3

      double precision bin_t3,bin_u3,bin_v3
      double precision bin_t4,bin_u4,bin_v4
      double precision P3
      double precision Pmat3(3)
      double precision KAB3
      double precision XPA3,YPA3,ZPA3
      double precision XPB3,YPB3,ZPB3
      double precision ans
      double precision xG2Vep
      double precision zero
      parameter(zero=0.0d+00)

C  Get basis function contraction values for particle 3:
      call contract_BF(A3,B3,Amat3,Bmat3,p3,Pmat3,KAB3,
     x                 XPA3,YPA3,ZPA3,XPB3,YPB3,ZPB3)

      IX3=t1+t2-xb
      JX3=u1+u2-yb
      KX3=v1+v2-zb

      G2Vep=zero

      do t3=0,I3
       call cbinom(I3,t3,bin_t3)
       do t4=0,L3
        call cbinom(L3,t4,bin_t4)
        do u3=0,J3
         call cbinom(J3,u3,bin_u3)
         do u4=0,M3
          call cbinom(M3,u4,bin_u4)
          do v3=0,K3
           call cbinom(K3,v3,bin_v3)
           do v4=0,N3
            call cbinom(N3,v4,bin_v4)


            LX3=t3+t4
            MX3=u3+u4
            NX3=v3+v4

c           call cws_gam1_xggvee(I1,J1,K1,A1,Amat1,
c           call G2_MD_xggvee(I1,J1,K1,A1,Amat1,
c    1                        IX3,JX3,KX3,aC23,Pmat2,
c    2                        L1,M1,N1,B1,Bmat1,
c    3                        LX3,MX3,NX3,p3,Pmat3,
c    4                        gam13,zero,xG2Vep)

            call pgivee(I1,J1,K1,A1,Amat1,
     1                  IX3,JX3,KX3,aC23,Pmat2,
     2                  L1,M1,N1,B1,Bmat1,
     3                  LX3,MX3,NX3,p3,Pmat3,
     4                  gam13,zero,xG2Vep)


            ans=bin_t3*bin_t4
     x         *bin_u3*bin_u4
     x         *bin_v3*bin_v4
     x         *XPA3**(I3-t3)
     x         *XPB3**(L3-t4)
     x         *YPA3**(J3-u3)
     x         *YPB3**(M3-u4)
     x         *ZPA3**(K3-v3)
     x         *ZPB3**(N3-v4)
     x         *xG2Vep

            G2Vep=G2Vep+ans

           end do
          end do
         end do
        end do
       end do
      end do

      G2Vep=G2Vep*KAB3


      return
      end



C=======================================================================
      subroutine G3Vee_AUX_g13V12(I1,J1,K1,A1,Amat1,
     *                            I2,J2,K2,A2,Amat2,
     *                            I3,J3,K3,A3,Amat3,
     *                            L1,M1,N1,B1,Bmat1,
     *                            L2,M2,N2,B2,Bmat2,
     *                            L3,M3,N3,B3,Bmat3,
     *                            gamA13,gamA23,
     *                            gamB13,gamB23,
     *                            xgVEE)

C Evaluates the following Vee integrals:
C xgVEEg1  ::  gA(1,3) V(1,2) gB(1,3) 
C
C calls G2_MD_xggvee to evaluate auxilliary integrals
C
C Note:  Proton is assumed to be particle 3
C=======================================================================
      implicit none

C Input Variables
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
      double precision gamA13
      double precision gamA23
      double precision gamB13
      double precision gamB23

C Variables Returned
      double precision xgVEE

C Local Variables
      integer t1,t2,xb 
      integer u1,u2,yb 
      integer v1,v2,zb 
      double precision bin_t1,bin_t2,bin_xb
      double precision bin_u1,bin_u2,bin_yb
      double precision bin_v1,bin_v2,bin_zb
      double precision gam13
      double precision gam23
      double precision p3
      double precision Pmat3(3)
      double precision KAB3
      double precision XPA3,YPA3,ZPA3
      double precision XPB3,YPB3,ZPB3
      double precision C13
      double precision aC13
      double precision a
      double precision G3Vee
      double precision CHI
      double precision ans
      double precision G2Vee



      gam13=gamA13+gamB13
      gam23=gamA23+gamB23


C  Get basis function contraction values for particle 3:
      call contract_BF(A3,B3,Amat3,Bmat3,p3,Pmat3,KAB3,
     x                 XPA3,YPA3,ZPA3,XPB3,YPB3,ZPB3)

      C13=gam13/(p3+gam13)
      aC13=p3*C13

      a=p3+gam13

      G3Vee=0.0d+00

      do t1=0,I3
       call cbinom(I3,t1,bin_t1)
       do t2=0,L3
        call cbinom(L3,t2,bin_t2)
        do xb=0,t1+t2
         call cbinom(t1+t2,xb,bin_xb)

         do u1=0,J3
          call cbinom(J3,u1,bin_u1)
          do u2=0,M3
           call cbinom(M3,u2,bin_u2)
           do yb=0,u1+u2
            call cbinom(u1+u2,yb,bin_yb)

            do v1=0,K3
             call cbinom(K3,v1,bin_v1)
             do v2=0,N3
              call cbinom(N3,v2,bin_v2)
              do zb=0,v1+v2
               call cbinom(v1+v2,zb,bin_zb)

               call evaluate_CHI(xb,yb,zb,a,CHI)

C Calculate the "auxiliary" 2-particle integral 
               call aux_G2Vee(I1,J1,K1,A1,Amat1,
     x                        I2,J2,K2,A2,Amat2,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L2,M2,N2,B2,Bmat2,
     x                        aC13,Pmat3,
     x                        t1,t2,xb,
     x                        u1,u2,yb,
     x                        v1,v2,zb,
     x                        G2Vee)

               
               ans=bin_t1*bin_t2*bin_xb
     x            *bin_u1*bin_u2*bin_yb
     x            *bin_v1*bin_v2*bin_zb
     x            *XPA3**(I3-t1) * YPA3**(J3-u1) * ZPA3**(K3-v1)
     x            *XPB3**(L3-t2) * YPB3**(M3-u2) * ZPB3**(N3-v2)
     x            *C13**(t1+t2-xb)
     x            *C13**(u1+u2-yb)
     x            *C13**(v1+v2-zb)
     x            *CHI*G2Vee


               G3Vee=G3Vee+ans
c     write(*,*)'---------------'
c     write(*,*)'bin_t1=',bin_t1
c     write(*,*)'bin_t2=',bin_t2
c     write(*,*)'bin_u1=',bin_u1
c     write(*,*)'bin_u2=',bin_u2
c     write(*,*)'bin_v1=',bin_v1
c     write(*,*)'bin_v2=',bin_v2
c     write(*,*)'CHI=',CHI
c     write(*,*)'XPA3=',XPA3**(I3-t1)
c     write(*,*)'XPB3=',XPB3**(l3-t2)
c     write(*,*)'YPA3=',YPA3**(j3-u1)
c     write(*,*)'YPB3=',YPB3**(m3-u2)
c     write(*,*)'ZPA3=',ZPA3**(k3-v1)
c     write(*,*)'ZPB3=',ZPB3**(n3-v2)
c     write(*,*)'V1=',V1
c     write(*,*)'V2=',V2
c     write(*,*)'zb=',zb
c     write(*,*)'C13=',C13
c     write(*,*)'aC13=',aC13
c     write(*,*)'C13t=',C13**(t1+t2-xb)
c     write(*,*)'C13u=',C13**(u1+u2-yb)
c     write(*,*)'C13v=',C13**(v1+v2-zb)
c     write(*,*)'G2Vee=',G2Vee
c     write(*,*)'KAB3=',KAB3
c     write(*,*)'---------------'

              end do
             end do
            end do

           end do
          end do
         end do

        end do
       end do
      end do

      G3Vee=G3Vee*KAB3
      xgVEE=G3Vee


      return
      end

C=======================================================================
      subroutine aux_G2Vee(I1,J1,K1,A1,Amat1,
     x                     I2,J2,K2,A2,Amat2,
     x                     L1,M1,N1,B1,Bmat1,
     x                     L2,M2,N2,B2,Bmat2,
     x                     aC13,Pmat3,
     x                     t1,t2,xb,
     x                     u1,u2,yb,
     x                     v1,v2,zb,
     x                     G2Vee)

C=======================================================================
      implicit none

C Input Variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2
      integer t1,t2,xb
      integer u1,u2,yb
      integer v1,v2,zb
       
      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
      double precision aC13,Pmat3(3)

C Variables Returned
      double precision G2Vee

C Local Variables
      integer t3,u3,v3
      integer t4,u4,v4
      integer IX1,JX1,KX1
      integer LX1,MX1,NX1

      double precision bin_t3,bin_u3,bin_v3
      double precision bin_t4,bin_u4,bin_v4
      double precision P1
      double precision Pmat1(3)
      double precision KAB1
      double precision XPA1,YPA1,ZPA1
      double precision XPB1,YPB1,ZPB1
      double precision ans
      double precision xG2Vee
      double precision zero
      parameter(zero=0.0d+00)

C  Get basis function contraction values for particle 1:
      call contract_BF(A1,B1,Amat1,Bmat1,P1,Pmat1,KAB1,
     x                 XPA1,YPA1,ZPA1,XPB1,YPB1,ZPB1)

      IX1=t1+t2-xb
      JX1=u1+u2-yb
      KX1=v1+v2-zb

      G2Vee=zero

      do t3=0,I1
       call cbinom(I1,t3,bin_t3)
       do t4=0,L1
        call cbinom(L1,t4,bin_t4)
        do u3=0,J1
         call cbinom(J1,u3,bin_u3)
         do u4=0,M1
          call cbinom(M1,u4,bin_u4)
          do v3=0,K1
           call cbinom(K1,v3,bin_v3)
           do v4=0,N1
            call cbinom(N1,v4,bin_v4)


            LX1=t3+t4
            MX1=u3+u4
            NX1=v3+v4

c           call G2_MD_xggvee(IX1,JX1,KX1,aC13,Pmat3,
c    x                        I2,J2,K2,A2,Amat2,
c    3                        LX1,MX1,NX1,P1,Pmat1,
c    2                        L2,M2,N2,B2,Bmat2,
c    4                        zero,zero,xG2Vee)

            call pgivee(IX1,JX1,KX1,aC13,Pmat3, 
     1                  I2,J2,K2,A2,Amat2, 
     2                  LX1,MX1,NX1,P1,Pmat1,
     3                  L2,M2,N2,B2,Bmat2, 
     4                  zero,zero,xG2Vee)  

c           call gfvee(IX1,JX1,KX1,aC13,Pmat3,
c    1                 I2,J2,K2,A2,Amat2,
c    2                 LX1,MX1,NX1,P1,Pmat1,
c    3                 L2,M2,N2,B2,Bmat2,
c    4                 xG2Vee)


            ans=bin_t3*bin_t4
     x         *bin_u3*bin_u4
     x         *bin_v3*bin_v4
     x         *XPA1**(I1-t3)
     x         *YPA1**(J1-u3)
     x         *ZPA1**(K1-v3)
     x         *XPB1**(L1-t4)
     x         *YPB1**(M1-u4)
     x         *ZPB1**(N1-v4)
     x         *xG2Vee

            G2Vee=G2Vee+ans

c     write(*,*)'-X-X-X-X-X-X-X-'
c     write(*,*)'IX1=',IX1
c     write(*,*)'JX1=',JX1
c     write(*,*)'KX1=',KX1
c     write(*,*)'aC13=',aC13
c     write(*,*)'Pmat3=',Pmat3
c     write(*,*)'I2=',I2
c     write(*,*)'J2=',J2
c     write(*,*)'K2=',K2
c     write(*,*)'A2=',A2
c     write(*,*)'Amat2=',Amat2
c     write(*,*)'LX1=',LX1
c     write(*,*)'MX1=',MX1
c     write(*,*)'NX1=',NX1
c     write(*,*)'P1=',P1
c     write(*,*)'Pmat1=',Pmat1
c     write(*,*)'L2=',L2
c     write(*,*)'M2=',M2
c     write(*,*)'N2=',N2
c     write(*,*)'B2=',B2
c     write(*,*)'Bmat2=',Bmat2
c     write(*,*)'bin_t3=',bin_t3
c     write(*,*)'bin_t4=',bin_t4
c     write(*,*)'bin_u3=',bin_u3
c     write(*,*)'bin_u4=',bin_u4
c     write(*,*)'bin_v3=',bin_v3
c     write(*,*)'bin_v4=',bin_v4
c     write(*,*)'XPA1=',XPA1**(I1-t3)
c     write(*,*)'XPB1=',XPB1**(l1-t4)
c     write(*,*)'YPA1=',YPA1**(j1-u3)
c     write(*,*)'YPB1=',YPB1**(m1-u4)
c     write(*,*)'ZPA1=',ZPA1**(k1-v3)
c     write(*,*)'ZPB1=',ZPB1**(n1-v4)
c     write(*,*)'V3=',V3
c     write(*,*)'V4=',V4
c     write(*,*)'zb=',zb
c     write(*,*)'C13=',C13
c     write(*,*)'aC13=',aC13
c     write(*,*)'xG2Vee=',xG2Vee
c     write(*,*)'KAB1=',KAB1
c     write(*,*)'-X-X-X-X-X-X-X-'

           end do
          end do
         end do
        end do
       end do
      end do

      G2Vee=G2Vee*KAB1


      return
      end



C=======================================================================
      subroutine G3Vee_AUX_g23V12(I1,J1,K1,A1,Amat1,
     *                            I2,J2,K2,A2,Amat2,
     *                            I3,J3,K3,A3,Amat3,
     *                            L1,M1,N1,B1,Bmat1,
     *                            L2,M2,N2,B2,Bmat2,
     *                            L3,M3,N3,B3,Bmat3,
     *                            gamA13,gamA23,
     *                            gamB13,gamB23,
     *                            xgVEE)

C Evaluates the following Vee integrals:
C xgVEE  ::  gA(2,3) V(1,2) 
C
C calls pgivee --or-- G2_MD_xggvee to evaluate auxilliary integrals
C
C Note:  Proton is assumed to be particle 3
C=======================================================================
      implicit none

C Input Variables
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
      double precision gamA13
      double precision gamA23
      double precision gamB13
      double precision gamB23

C Variables Returned
      double precision xgVEE

C Local Variables
      integer t1,t2,xb 
      integer u1,u2,yb 
      integer v1,v2,zb 
      double precision bin_t1,bin_t2,bin_xb
      double precision bin_u1,bin_u2,bin_yb
      double precision bin_v1,bin_v2,bin_zb
      double precision gam13
      double precision gam23
      double precision p3
      double precision Pmat3(3)
      double precision KAB3
      double precision XPA3,YPA3,ZPA3
      double precision XPB3,YPB3,ZPB3
      double precision C23
      double precision aC23
      double precision a
      double precision G3Vee
      double precision CHI
      double precision ans
      double precision G2Vee



      gam13=gamA13+gamB13
      gam23=gamA23+gamB23


C  Get basis function contraction values for particle 3:
      call contract_BF(A3,B3,Amat3,Bmat3,p3,Pmat3,KAB3,
     x                 XPA3,YPA3,ZPA3,XPB3,YPB3,ZPB3)

      C23=gam23/(p3+gam23)
      aC23=p3*C23

      a=p3+gam23

      G3Vee=0.0d+00

      do t1=0,I3
       call cbinom(I3,t1,bin_t1)
       do t2=0,L3
        call cbinom(L3,t2,bin_t2)
        do xb=0,t1+t2
         call cbinom(t1+t2,xb,bin_xb)

         do u1=0,J3
          call cbinom(J3,u1,bin_u1)
          do u2=0,M3
           call cbinom(M3,u2,bin_u2)
           do yb=0,u1+u2
            call cbinom(u1+u2,yb,bin_yb)

            do v1=0,K3
             call cbinom(K3,v1,bin_v1)
             do v2=0,N3
              call cbinom(N3,v2,bin_v2)
              do zb=0,v1+v2
               call cbinom(v1+v2,zb,bin_zb)

               call evaluate_CHI(xb,yb,zb,a,CHI)

C Calculate the "auxiliary" 2-particle integral 
               call aux_G2Vee2(I1,J1,K1,A1,Amat1,
     x                         I2,J2,K2,A2,Amat2,
     x                         L1,M1,N1,B1,Bmat1,
     x                         L2,M2,N2,B2,Bmat2,
     x                         aC23,Pmat3,
     x                         t1,t2,xb,
     x                         u1,u2,yb,
     x                         v1,v2,zb,
     x                         G2Vee)

               
               ans=bin_t1*bin_t2*bin_xb
     x            *bin_u1*bin_u2*bin_yb
     x            *bin_v1*bin_v2*bin_zb
     x            *XPA3**(I3-t1) * YPA3**(J3-u1) * ZPA3**(K3-v1)
     x            *XPB3**(L3-t2) * YPB3**(M3-u2) * ZPB3**(N3-v2)
     x            *C23**(t1+t2-xb)
     x            *C23**(u1+u2-yb)
     x            *C23**(v1+v2-zb)
     x            *CHI*G2Vee


               G3Vee=G3Vee+ans

              end do
             end do
            end do

           end do
          end do
         end do

        end do
       end do
      end do

      G3Vee=G3Vee*KAB3
      xgVEE=G3Vee


      return
      end

C=======================================================================
      subroutine aux_G2Vee2(I1,J1,K1,A1,Amat1,
     x                      I2,J2,K2,A2,Amat2,
     x                      L1,M1,N1,B1,Bmat1,
     x                      L2,M2,N2,B2,Bmat2,
     x                      aC23,Pmat3,
     x                      t1,t2,xb,
     x                      u1,u2,yb,
     x                      v1,v2,zb,
     x                      G2Vee)

C=======================================================================
      implicit none

C Input Variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer L1,M1,N1
      integer L2,M2,N2
      integer t1,t2,xb
      integer u1,u2,yb
      integer v1,v2,zb
       
      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
      double precision aC23,Pmat3(3)

C Variables Returned
      double precision G2Vee

C Local Variables
      integer t3,u3,v3
      integer t4,u4,v4
      integer IX2,JX2,KX2
      integer LX2,MX2,NX2

      double precision bin_t3,bin_u3,bin_v3
      double precision bin_t4,bin_u4,bin_v4
      double precision P2
      double precision Pmat2(3)
      double precision KAB2
      double precision XPA2,YPA2,ZPA2
      double precision XPB2,YPB2,ZPB2
      double precision ans
      double precision xG2Vee
      double precision zero
      parameter(zero=0.0d+00)

C  Get basis function contraction values for particle 1:
      call contract_BF(A2,B2,Amat2,Bmat2,P2,Pmat2,KAB2,
     x                 XPA2,YPA2,ZPA2,XPB2,YPB2,ZPB2)

      IX2=t1+t2-xb
      JX2=u1+u2-yb
      KX2=v1+v2-zb

      G2Vee=zero

      do t3=0,I2
       call cbinom(I2,t3,bin_t3)
       do t4=0,L2
        call cbinom(L2,t4,bin_t4)
        do u3=0,J2
         call cbinom(J2,u3,bin_u3)
         do u4=0,M2
          call cbinom(M2,u4,bin_u4)
          do v3=0,K2
           call cbinom(K2,v3,bin_v3)
           do v4=0,N2
            call cbinom(N2,v4,bin_v4)


            LX2=t3+t4
            MX2=u3+u4
            NX2=v3+v4

c           call cws_gam1_xggvee(I1,J1,K1,A1,Amat1,
c           call G2_MD_xggvee(I1,J1,K1,A1,Amat1,
c    x                        IX2,JX2,KX2,aC23,Pmat3,
c    2                        L1,M1,N1,B1,Bmat1,
c    3                        LX2,MX2,NX2,P2,Pmat2,
c    4                        zero,zero,xG2Vee)

            call pgivee(I1,J1,K1,A1,Amat1, 
     x                  IX2,JX2,KX2,aC23,Pmat3, 
     x                  L1,M1,N1,B1,Bmat1, 
     x                  LX2,MX2,NX2,P2,Pmat2,
     x                  zero,zero,xG2Vee)  

            ans=bin_t3*bin_t4
     x         *bin_u3*bin_u4
     x         *bin_v3*bin_v4
     x         *XPA2**(I2-t3)
     x         *YPA2**(J2-u3)
     x         *ZPA2**(K2-v3)
     x         *XPB2**(L2-t4)
     x         *YPB2**(M2-u4)
     x         *ZPB2**(N2-v4)
     x         *xG2Vee

            G2Vee=G2Vee+ans

           end do
          end do
         end do
        end do
       end do
      end do

      G2Vee=G2Vee*KAB2


      return
      end


