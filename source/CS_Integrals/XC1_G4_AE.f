C G4Vep_AUX_g14g24V34.f
C G4Vee_AUX_g14g34V12.f    
C G4Vee_AUXDFT_g23g14V12.f 
C=======================================================================
      subroutine G4Vep_AUX_g14g24V34(I1,J1,K1,A1,Amat1,
     *                               I2,J2,K2,A2,Amat2,
     *                               I3,J3,K3,A3,Amat3,
     *                               I4,J4,K4,A4,Amat4,
     *                               L1,M1,N1,B1,Bmat1,
     *                               L2,M2,N2,B2,Bmat2,
     *                               L3,M3,N3,B3,Bmat3,
     *                               L4,M4,N4,B4,Bmat4,
     *                               gamA14,gamA24,
     *                               gamB14,gamB24,
     *                               xgVEPg)

C Evaluates the following Vep integral:
C xgVEPg  ::  gA(1,4) V(3,4) gB(2,4)
C
C Calls aux_G3Vep 
C Uses G3Vep_AUX_g13g23V13 routine to evaluate auxiliary 3-particle ints
C Note:  Proton is assumed to be particle 4
C=======================================================================
      implicit none

C Input Variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4
      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision A3,Amat3(3)
      double precision A4,Amat4(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
      double precision B3,Bmat3(3)
      double precision B4,Bmat4(3)
      double precision gamA14
      double precision gamA24
      double precision gamB14
      double precision gamB24

C Variables Returned
      double precision xgVEPg

C Local Variables
      integer t1,t2,xb 
      integer u1,u2,yb 
      integer v1,v2,zb 
      double precision bin_t1,bin_t2,bin_xb
      double precision bin_u1,bin_u2,bin_yb
      double precision bin_v1,bin_v2,bin_zb
      double precision gam14
      double precision gam24
      double precision p2
      double precision Pmat2(3)
      double precision KAB2
      double precision XPA2,YPA2,ZPA2
      double precision XPB2,YPB2,ZPB2
      double precision C24
      double precision aC24
      double precision a
      double precision G4Vep
      double precision CHI
      double precision ans
      double precision G3Vep

      double precision ans1,ans2,ans3,ans4



      gam14=gamA14+gamB14
      gam24=gamA24+gamB24


C  Get basis function contraction values for particle 2:
      call contract_BF(A2,B2,Amat2,Bmat2,p2,Pmat2,KAB2,
     x                 XPA2,YPA2,ZPA2,XPB2,YPB2,ZPB2)

      C24=gam24/(p2+gam24)
      aC24=p2*C24

      a=p2+gam24

      G4Vep=0.0d+00

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

C Calculate the "auxiliary" 3-particle integral 
               call aux_G3Vep(I1,J1,K1,A1,Amat1,
     x                        I3,J3,K3,A3,Amat3,
     x                        I4,J4,K4,A4,Amat4,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L3,M3,N3,B3,Bmat3,
     x                        L4,M4,N4,B4,Bmat4,
     x                        aC24,Pmat2,
     x                        t1,t2,xb,
     x                        u1,u2,yb,
     x                        v1,v2,zb,
     x                        gam14,G3Vep)

               
               ans=bin_t1*bin_t2*bin_xb
     x            *bin_u1*bin_u2*bin_yb
     x            *bin_v1*bin_v2*bin_zb
     x            *XPA2**(I2-t1) * YPA2**(J2-u1) * ZPA2**(K2-v1)
     x            *XPB2**(L2-t2) * YPB2**(M2-u2) * ZPB2**(N2-v2)
     x            *C24**(t1+t2-xb)
     x            *C24**(u1+u2-yb)
     x            *C24**(v1+v2-zb)
     x            *CHI*G3Vep




               G4Vep=G4Vep+ans

              end do
             end do
            end do

           end do
          end do
         end do

        end do
       end do
      end do

      G4Vep=G4Vep*KAB2
      xgVEPg=G4Vep

CCWS-DEBUG(
c           call gfovlap(I1,J1,K1,a1,Amat1,
c    2                   L1,M1,N1,b1,Bmat1,ans1)

c           call gfovlap(I2,J2,K2,a2,Amat2,
c    2                   L2,M2,N2,b2,Bmat2,ans2)

cc          call gfovlap(I3,J3,K3,a3,Amat3,
cc   2                   L3,M3,N3,b3,Bmat3,ans3)

c           call gfvee(I3,J3,K3,a3,Amat3,
c    2                 I4,J4,K4,a4,Amat4,
c    2                 L3,M3,N3,b3,Bmat3,
c    2                 L4,M4,N4,b4,Bmat4,
c    4                 ans4)
c           write(*,*)'      ANS1=',ans1
c           write(*,*)'      ANS2=',ans2
c           write(*,*)'      ANS4=',ans4*-1.0d+00
c           write(*,*)'SIMPLE ANS=',ans1*ans2*ans4*-1.0d+00

CCWS-DEBUG)


      return
      end

C=======================================================================
      subroutine aux_G3Vep(I1,J1,K1,A1,Amat1,
     x                     I3,J3,K3,A3,Amat3,
     x                     I4,J4,K4,A4,Amat4,
     x                     L1,M1,N1,B1,Bmat1,
     x                     L3,M3,N3,B3,Bmat3,
     x                     L4,M4,N4,B4,Bmat4,
     x                     aC24,Pmat2,
     x                     t1,t2,xb,
     x                     u1,u2,yb,
     x                     v1,v2,zb,
     x                     gam14,G3Vep)

C=======================================================================
      implicit none

C Input Variables
      integer I1,J1,K1
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L3,M3,N3
      integer L4,M4,N4
      integer t1,t2,xb
      integer u1,u2,yb
      integer v1,v2,zb
       
      double precision A1,Amat1(3)
      double precision A3,Amat3(3)
      double precision A4,Amat4(3)
      double precision B1,Bmat1(3)
      double precision B3,Bmat3(3)
      double precision B4,Bmat4(3)
      double precision aC24,Pmat2(3)
      double precision gam14

C Variables Returned
      double precision G3Vep


C Local Variables
      integer t3,u3,v3
      integer t4,u4,v4
      integer IX4,JX4,KX4
      integer LX4,MX4,NX4

      double precision bin_t3,bin_u3,bin_v3
      double precision bin_t4,bin_u4,bin_v4
      double precision P4
      double precision Pmat4(3)
      double precision KAB4
      double precision XPA4,YPA4,ZPA4
      double precision XPB4,YPB4,ZPB4
      double precision ans
      double precision xG3Vep
      double precision zero
      parameter(zero=0.0d+00)

C  Get basis function contraction values for particle 3:
      call contract_BF(A4,B4,Amat4,Bmat4,p4,Pmat4,KAB4,
     x                 XPA4,YPA4,ZPA4,XPB4,YPB4,ZPB4)

      IX4=t1+t2-xb
      JX4=u1+u2-yb
      KX4=v1+v2-zb

      G3Vep=zero

      do t3=0,I4
       call cbinom(I4,t3,bin_t3)
       do t4=0,L4
        call cbinom(L4,t4,bin_t4)
        do u3=0,J4
         call cbinom(J4,u3,bin_u3)
         do u4=0,M4
          call cbinom(M4,u4,bin_u4)
          do v3=0,K4
           call cbinom(K4,v3,bin_v3)
           do v4=0,N4
            call cbinom(N4,v4,bin_v4)


            LX4=t3+t4
            MX4=u3+u4
            NX4=v3+v4

C 2 ways to calculate this auxilliary integral---MCM or another AUX:

C MCM:
c           call cws_gam2_xgVeeg(I3,J3,K3,A3,Amat3,
c           call G3_MD_xgVeeg(I3,J3,K3,A3,Amat3,
c    x                        IX4,JX4,KX4,aC24,Pmat2,
c    x                        I1,J1,K1,A1,Amat1,
c    x                        L3,M3,N3,B3,Bmat3,
c    x                        LX4,MX4,NX4,P4,Pmat4,
c    x                        L1,M1,N1,B1,Bmat1,
c    x                        zero,zero,gam14,
c    x                        zero,zero,zero,
c    x                        xG3Vep)


C AUX:
            call G3Vep_AUX_g13g23V13(I3,J3,K3,A3,Amat3,
     x                               I1,J1,K1,A1,Amat1,
     x                               IX4,JX4,KX4,aC24,Pmat2,
     x                               L3,M3,N3,B3,Bmat3,
     x                               L1,M1,N1,B1,Bmat1,
     x                               LX4,MX4,NX4,P4,Pmat4,
     x                               zero,gam14,
     x                               zero,zero,
     x                               xG3Vep)

c           call G3Vep-AUX-g13g23V13(I3,J3,K3,A3,Amat3,
c    x                               I1,J1,K1,A1,Amat1,
c    x                            t1+t2-xb,u1+u2-yb,v1+v2-zb,aC24,Pmat2,
c    x                               L3,M3,N3,B3,Bmat3,
c    x                               L1,M1,N1,B1,Bmat1,
c    x                               t3+t4,u3+u4,v3+v4,P4,Pmat4,
c    x                               zero,gam23,
c    x                               zero,zero,
c    x                               xG3Vep)


            ans=bin_t3*bin_t4
     x         *bin_u3*bin_u4
     x         *bin_v3*bin_v4
     x         *XPA4**(I4-t3)
     x         *XPB4**(L4-t4)
     x         *YPA4**(J4-u3)
     x         *YPB4**(M4-u4)
     x         *ZPA4**(K4-v3)
     x         *ZPB4**(N4-v4)
     x         *xG3Vep

            G3Vep=G3Vep+ans

           end do
          end do
         end do
        end do
       end do
      end do

      G3Vep=G3Vep*KAB4


      return
      end



C=======================================================================
c     subroutine cws_G4_xgVeeg_AUX(I1,J1,K1,A1,Amat1,
      subroutine G4Vee_AUX_g14g34V12(I1,J1,K1,A1,Amat1,
     *                               I2,J2,K2,A2,Amat2,
     *                               I3,J3,K3,A3,Amat3,
     *                               I4,J4,K4,A4,Amat4,
     *                               L1,M1,N1,B1,Bmat1,
     *                               L2,M2,N2,B2,Bmat2,
     *                               L3,M3,N3,B3,Bmat3,
     *                               L4,M4,N4,B4,Bmat4,
     *                               gam14,gam34,
     *                               xgVeeg)

C Evaluates the following Vep integral:
C xgVEEg2 ::  gA(1,4) V(1,2) gB(3,4)
C
C Calls aux_G3Vee 
C Uses either G3_MD_xgVeeg or G3Vee_AUX_g13V12 to evaluate 
C  3-particle auxiliary integrals
C
C Note:  Proton is assumed to be particle 4
C=======================================================================
      implicit none

C Input Variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4

      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4

      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision A3,Amat3(3)
      double precision A4,Amat4(3)

      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
      double precision B3,Bmat3(3)
      double precision B4,Bmat4(3)

      double precision gam14
      double precision gam34

C Variables Returned
      double precision xgVeeg

C Local Variables
      integer t1,u1,v1
      integer t2,u2,v2
      integer xb,yb,zb

      double precision PI
      parameter(PI=3.14159265358979d+00)

      double precision P3
      double precision a
      double precision ans
      double precision CHI
      double precision Pmat3(3)
      double precision KAB3
      double precision XPA3,YPA3,ZPA3
      double precision XPB3,YPB3,ZPB3
      double precision C34
      double precision aC34
      double precision G4Vee
      double precision G3Vee
      double precision bin_t1,bin_t2
      double precision bin_u1,bin_u2
      double precision bin_v1,bin_v2
      double precision bin_xb
      double precision bin_yb
      double precision bin_zb
      
c     write(*,*)' IN cws_G4_xgVeeg_AUX'

C  Get basis function contraction values for particle 3:
      call contract_BF(A3,B3,Amat3,Bmat3,p3,Pmat3,KAB3,
     x                 XPA3,YPA3,ZPA3,XPB3,YPB3,ZPB3)

      C34=gam34/(p3+gam34)
      aC34=p3*C34

      a=p3+gam34

      G4Vee=0.0d+00

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

c     write(*,*)'CHECK CHI =',CHI-(PI/(p3+gam34))**(3.0d+00/2.0d+00)

C Calculate the "auxiliary" 3-particle integral 
               call aux_G3Vee(I1,J1,K1,A1,Amat1,
     x                        I2,J2,K2,A2,Amat2,
     x                        I4,J4,K4,A4,Amat4,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L2,M2,N2,B2,Bmat2,
     x                        L4,M4,N4,B4,Bmat4,
     x                        aC34,Pmat3,
     x                        t1,t2,xb, 
     x                        u1,u2,yb, 
     x                        v1,v2,zb,
     x                        gam14,G3Vee) 

               ans=bin_t1*bin_t2*bin_xb
     x            *bin_u1*bin_u2*bin_yb
     x            *bin_v1*bin_v2*bin_zb
     x            *XPA3**(I3-t1) * YPA3**(J3-u1) * ZPA3**(K3-v1)
     x            *XPB3**(L3-t2) * YPB3**(M3-u2) * ZPB3**(N3-v2)
     x            *C34**(t1+t2-xb)
     x            *C34**(u1+u2-yb)
     x            *C34**(v1+v2-zb)
     x            *CHI*G3Vee
c    x            *CHI_xb*CHI_yb*CHI_zb

c                 write(*,*)
c                 write(*,*)'bin_t1=',bin_t1
c                 write(*,*)'bin_t2=',bin_t2
c                 write(*,*)'bin_xb=',bin_xb
c                 write(*,*)'bin_u1=',bin_u1
c                 write(*,*)'bin_u2=',bin_u2
c                 write(*,*)'bin_yb=',bin_yb
c                 write(*,*)'bin_v1=',bin_v1
c                 write(*,*)'bin_v2=',bin_v2
c                 write(*,*)'bin_zb=',bin_zb
c                 write(*,*)'XPA3=',XPA3**(I3-t1)
c                 write(*,*)'YPA3=',YPA3**(j3-u1)
c                 write(*,*)'ZPA3=',ZPA3**(k3-v1)
c                 write(*,*)'XPB3=',XPB3**(l3-t2)
c                 write(*,*)'YPB3=',YPB3**(m3-u2)
c                 write(*,*)'ZPB3=',ZPB3**(n3-v2)
c                 write(*,*)'C34=',C34**(t1+t2-xb)
c                 write(*,*)'C34=',C34**(u1+u2-yb)
c                 write(*,*)'C34=',C34**(v1+v2-zb)
c                 write(*,*)'CHI=',CHI
c                 write(*,*)'G3Vee=',G3Vee
c                 write(*,*)'KAB3=',KAB3
c                 write(*,*)
    

               G4Vee=G4Vee+ans

              end do
             end do
            end do

           end do
          end do
         end do

        end do
       end do
      end do

      G4Vee=G4Vee*KAB3
      xgVeeg=G4Vee


      return 
      end


C=======================================================================
      subroutine aux_G3Vee(I1,J1,K1,A1,Amat1,
     x                     I2,J2,K2,A2,Amat2,
     x                     I4,J4,K4,A4,Amat4,
     x                     L1,M1,N1,B1,Bmat1,
     x                     L2,M2,N2,B2,Bmat2,
     x                     L4,M4,N4,B4,Bmat4,
     x                     aC34,Pmat3,
     x                     t1,t2,xb, 
     x                     u1,u2,yb, 
     x                     v1,v2,zb,
     x                     gam14,G3Vee) 

C=======================================================================
      implicit none

C Input Variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L4,M4,N4
      integer t1,t2,xb
      integer u1,u2,yb
      integer v1,v2,zb

      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision A4,Amat4(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
      double precision B4,Bmat4(3)
      double precision aC34,Pmat3(3)
      double precision gam14

C Variables Returned
      double precision G3Vee

C Local Variables
      integer t3,u3,v3
      integer t4,u4,v4
      integer I3,J3,K3
      integer L3,M3,N3

      double precision bin_t3,bin_u3,bin_v3
      double precision bin_t4,bin_u4,bin_v4
      double precision xG3Vee
      double precision gamA12
      double precision gamA13
      double precision gamA23
      double precision gamB12
      double precision gamB13
      double precision gamB23
      double precision P4
      double precision Pmat4(3)
      double precision KAB4
      double precision XPA4,YPA4,ZPA4
      double precision XPB4,YPB4,ZPB4
      double precision zero
      parameter(zero=0.0d+00)
      double precision ans


C  Get basis function contraction values for particle 4:
      call contract_BF(A4,B4,Amat4,Bmat4,p4,Pmat4,KAB4,
     x                 XPA4,YPA4,ZPA4,XPB4,YPB4,ZPB4)

      G3Vee=zero

      do t3=0,I4
       call cbinom(I4,t3,bin_t3)
       do t4=0,L4
        call cbinom(L4,t4,bin_t4)
        do u3=0,J4
         call cbinom(J4,u3,bin_u3)
         do u4=0,M4
          call cbinom(M4,u4,bin_u4)
          do v3=0,K4
           call cbinom(K4,v3,bin_v3)
           do v4=0,N4
            call cbinom(N4,v4,bin_v4)

            gamA12=zero
            gamA13=gam14
            gamA23=zero
            gamB12=zero
            gamB13=zero
            gamB23=zero

            I3=t1+t2-xb
            J3=u1+u2-yb
            K3=v1+v2-zb

            L3=t3+t4
            M3=u3+u4
            N3=v3+v4

c                 write(*,*)
c                 write(*,*)'========================================='
c                 write(*,*)'I1=',I1
c                 write(*,*)'J1=',J1
c                 write(*,*)'K1=',K1
c                 write(*,*)'I2=',I2
c                 write(*,*)'J2=',J2
c                 write(*,*)'K2=',K2
c                 write(*,*)'I3=',I3
c                 write(*,*)'J3=',J3
c                 write(*,*)'K3=',K3
c                 write(*,*)'L1=',L1
c                 write(*,*)'M1=',M1
c                 write(*,*)'N1=',N1
c                 write(*,*)'L2=',L2
c                 write(*,*)'M2=',M2
c                 write(*,*)'N2=',N2
c                 write(*,*)'L3=',L3
c                 write(*,*)'M3=',M3
c                 write(*,*)'N3=',N3
c                 write(*,*)'A1=',A1
c                 write(*,*)'A2=',A2
c                 write(*,*)'aC34=',aC34
c                 write(*,*)'B1=',B1
c                 write(*,*)'B2=',B2
c                 write(*,*)'P4=',P4
c                 write(*,*)'Amat1=',Amat1(1)+Amat1(2)+Amat1(3)
c                 write(*,*)'Amat2=',Amat2(1)+Amat2(2)+Amat2(3)
c                 write(*,*)'Pmat3=',Pmat3(1)+Pmat3(2)+Pmat3(3)
c                 write(*,*)'Bmat1=',Bmat1(1)+Bmat1(2)+Bmat1(3)
c                 write(*,*)'Bmat2=',Bmat2(1)+Bmat2(2)+Bmat2(3)
c                 write(*,*)'Pmat4=',Pmat4(1)+Pmat4(2)+Pmat4(3)
c                 write(*,*)'gamA12=',gamA12
c                 write(*,*)'gamA13=',gamA13
c                 write(*,*)'gamA23=',gamA23
c                 write(*,*)'gamB12=',gamB12
c                 write(*,*)'gamB13=',gamB13
c                 write(*,*)'gamB23=',gamB23
c                 write(*,*)

            xG3Vee=zero

c           call cws_gam2_xgVeeg(I1,J1,K1,A1,Amat1,
c           call G3_MD_xgVeeg(I1,J1,K1,A1,Amat1,
c    *                        I2,J2,K2,A2,Amat2,
c    *                        I3,J3,K3,aC34,Pmat3,
c    *                        L1,M1,N1,B1,Bmat1,
c    *                        L2,M2,N2,B2,Bmat2,
c    *                        L3,M3,N3,P4,Pmat4,
c    *                        gamA12,gamA13,gamA23,
c    *                        gamB12,gamB13,gamB23,
c    *                        xG3Vee)

            call G3Vee_AUX_g13V12(I1,J1,K1,A1,Amat1,
     *                            I2,J2,K2,A2,Amat2,
     *                            I3,J3,K3,aC34,Pmat3,
     *                            L1,M1,N1,B1,Bmat1,
     *                            L2,M2,N2,B2,Bmat2,
     *                            L3,M3,N3,P4,Pmat4,
     *                            gamA13,gamA23,
     *                            gamB13,gamB23,
     *                            xG3Vee)


            ans=bin_t3*bin_t4
     x         *bin_u3*bin_u4
     x         *bin_v3*bin_v4
     x         *XPA4**(I4-t3) 
     x         *XPB4**(L4-t4) 
     x         *YPA4**(J4-u3) 
     x         *YPB4**(M4-u4) 
     x         *ZPA4**(K4-v3) 
     x         *ZPB4**(N4-v4) 
     x         *xG3Vee

            G3Vee=G3Vee+ans

c                 write(*,*)
c                 write(*,*)'bin_t3=',bin_t3
c                 write(*,*)'bin_t4=',bin_t4
c                 write(*,*)'bin_u3=',bin_u3
c                 write(*,*)'bin_u4=',bin_u4
c                 write(*,*)'bin_v3=',bin_v3
c                 write(*,*)'bin_v4=',bin_v4
c                 write(*,*)'XPA4=',XPA4**(I4-t3)
c                 write(*,*)'XPB4=',XPB4**(l4-t4)
c                 write(*,*)'YPA4=',YPA4**(j4-u3)
c                 write(*,*)'YPB4=',YPB4**(m4-u4)
c                 write(*,*)'ZPA3=',ZPA4**(k4-v3)
c                 write(*,*)'ZPB3=',ZPB4**(n4-v4)
c                 write(*,*)'xG3Vee=',xG3Vee
c                 write(*,*)'KAB4=',KAB4
c                 write(*,*)
           end do
          end do
         end do
        end do
       end do
      end do

      G3Vee=G3Vee*KAB4


      return
      end


C=======================================================================
      subroutine G4Vee_AUXDFT_g23g14V12(I1,J1,K1,A1,Amat1,
     *                                  I2,J2,K2,A2,Amat2,
     *                                  I3,J3,K3,A3,Amat3,
     *                                  I4,J4,K4,A4,Amat4,
     *                                  L1,M1,N1,B1,Bmat1,
     *                                  L2,M2,N2,B2,Bmat2,
     *                                  L3,M3,N3,B3,Bmat3,
     *                                  L4,M4,N4,B4,Bmat4,
     *                                  gamA23,gamA14,
     *                                  gamB23,gamB14,
     *                                  xgVEEg)

C Evaluates the following Vee integral:
C xgVEEg  ::  g(2,3) V(1,2) g(1,4)
C
C Calls: aux123_G3Vee ===> G3Vee_AUX_g23V12 
C Note:  Protons are assumed to be particles 3,4
C=======================================================================
      implicit none

C Input Variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer I4,J4,K4
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer L4,M4,N4
      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision A3,Amat3(3)
      double precision A4,Amat4(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
      double precision B3,Bmat3(3)
      double precision B4,Bmat4(3)
      double precision gamA14
      double precision gamA23
      double precision gamB14
      double precision gamB23

C Variables Returned
      double precision xgVEEg

C Local Variables
      integer t1,t2,xb 
      integer u1,u2,yb 
      integer v1,v2,zb 
      double precision bin_t1,bin_t2,bin_xb
      double precision bin_u1,bin_u2,bin_yb
      double precision bin_v1,bin_v2,bin_zb
      double precision gam14
      double precision gam23
      double precision p4
      double precision Pmat4(3)
      double precision KAB4
      double precision XPA4,YPA4,ZPA4
      double precision XPB4,YPB4,ZPB4
      double precision C14
      double precision aC14
      double precision a
      double precision G4Vee
      double precision CHI
      double precision ans
      double precision G3Vee



      gam14=gamA14+gamB14
      gam23=gamA23+gamB23


C  Get basis function contraction values for particle 4:
      call contract_BF(A4,B4,Amat4,Bmat4,p4,Pmat4,KAB4,
     x                 XPA4,YPA4,ZPA4,XPB4,YPB4,ZPB4)

      C14=gam14/(p4+gam14)
      aC14=p4*C14

      a=p4+gam14

      G4Vee=0.0d+00

      do t1=0,I4
       call cbinom(I4,t1,bin_t1)
       do t2=0,L4
        call cbinom(L4,t2,bin_t2)
        do xb=0,t1+t2
         call cbinom(t1+t2,xb,bin_xb)

         do u1=0,J4
          call cbinom(J4,u1,bin_u1)
          do u2=0,M4
           call cbinom(M4,u2,bin_u2)
           do yb=0,u1+u2
            call cbinom(u1+u2,yb,bin_yb)

            do v1=0,K4
             call cbinom(K4,v1,bin_v1)
             do v2=0,N4
              call cbinom(N4,v2,bin_v2)
              do zb=0,v1+v2
               call cbinom(v1+v2,zb,bin_zb)

               call evaluate_CHI(xb,yb,zb,a,CHI)

C Calculate the "auxiliary" 3-particle integral 
               call aux123_G3Vee(I1,J1,K1,A1,Amat1,
     x                           I2,J2,K2,A2,Amat2,
     x                           I3,J3,K3,A3,Amat3,
     x                           L1,M1,N1,B1,Bmat1,
     x                           L2,M2,N2,B2,Bmat2,
     x                           L3,M3,N3,B3,Bmat3,
     x                           aC14,Pmat4,
     x                           t1,t2,xb,
     x                           u1,u2,yb,
     x                           v1,v2,zb,
     x                           gam23,G3Vee)

               
               ans=bin_t1*bin_t2*bin_xb
     x            *bin_u1*bin_u2*bin_yb
     x            *bin_v1*bin_v2*bin_zb
     x            *XPA4**(I4-t1) * YPA4**(J4-u1) * ZPA4**(K4-v1)
     x            *XPB4**(L4-t2) * YPB4**(M4-u2) * ZPB4**(N4-v2)
     x            *C14**(t1+t2-xb)
     x            *C14**(u1+u2-yb)
     x            *C14**(v1+v2-zb)
     x            *CHI*G3Vee




               G4Vee=G4Vee+ans

              end do
             end do
            end do

           end do
          end do
         end do

        end do
       end do
      end do

      G4Vee=G4Vee*KAB4
      xgVEEg=G4Vee


      return
      end

C=======================================================================
      subroutine aux123_G3Vee(I1,J1,K1,A1,Amat1,
     x                        I2,J2,K2,A2,Amat2,
     x                        I3,J3,K3,A3,Amat3,
     x                        L1,M1,N1,B1,Bmat1,
     x                        L2,M2,N2,B2,Bmat2,
     x                        L3,M3,N3,B3,Bmat3,
     x                        aC14,Pmat4,
     x                        t1,t2,xb,
     x                        u1,u2,yb,
     x                        v1,v2,zb,
     x                        gam23,G3Vee)

C=======================================================================
      implicit none

C Input Variables
      integer I1,J1,K1
      integer I2,J2,K2
      integer I3,J3,K3
      integer L1,M1,N1
      integer L2,M2,N2
      integer L3,M3,N3
      integer t1,t2,xb
      integer u1,u2,yb
      integer v1,v2,zb
       
      double precision A1,Amat1(3)
      double precision A2,Amat2(3)
      double precision A3,Amat3(3)
      double precision B1,Bmat1(3)
      double precision B2,Bmat2(3)
      double precision B3,Bmat3(3)
      double precision aC14,Pmat4(3)
      double precision gam23

C Variables Returned
      double precision G3Vee

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
      double precision xG3Vee
      double precision zero
      parameter(zero=0.0d+00)

C  Get basis function contraction values for particle 1:
      call contract_BF(A1,B1,Amat1,Bmat1,p1,Pmat1,KAB1,
     x                 XPA1,YPA1,ZPA1,XPB1,YPB1,ZPB1)

      IX1=t1+t2-xb
      JX1=u1+u2-yb
      KX1=v1+v2-zb

      G3Vee=zero

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

C calculate this auxilliary integral:

C AUX:
      call G3Vee_AUX_g23V12(IX1,JX1,KX1,aC14,Pmat4,
     *                      I2,J2,K2,A2,Amat2,
     *                      I3,J3,K3,A3,Amat3,
     *                      LX1,MX1,NX1,P1,Pmat1,
     *                      L2,M2,N2,B2,Bmat2,
     *                      L3,M3,N3,B3,Bmat3,
     *                      zero,gam23,
     *                      zero,zero,
     *                      xg3Vee)

c     call G3Vee_AUX_g23V12(t1+t2-xb,u1+u2-yb,v1+v2-zb,aC14,Pmat4,
c    *                      I2,J2,K2,A2,Amat2,
c    *                      I3,J3,K3,A3,Amat3,
c    *                      t3+t4,u3+u4,v3+v4,P1,Pmat1,
c    *                      L2,M2,N2,B2,Bmat2,
c    *                      L3,M3,N3,B3,Bmat3,
c    *                      zero,gamA23,
c    *                      zero,zero,
c    *                      xgVEE)



            ans=bin_t3*bin_t4
     x         *bin_u3*bin_u4
     x         *bin_v3*bin_v4
     x         *XPA1**(I1-t3)
     x         *XPB1**(L1-t4)
     x         *YPA1**(J1-u3)
     x         *YPB1**(M1-u4)
     x         *ZPA1**(K1-v3)
     x         *ZPB1**(N1-v4)
     x         *xG3Vee

            G3Vee=G3Vee+ans

           end do
          end do
         end do
        end do
       end do
      end do

      G3Vee=G3Vee*KAB1


      return
      end

