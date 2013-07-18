      subroutine g2_vee_shell(P1,P2,Pmat1,Pmat2,gamm,
     $                       ng1,GAM_1,ise,jse,isn,jsn,
     $                       esh,psh,nesh,nnsh,
     $                       AMPEB2C,AGEBFCC,AGNBFCC,npebf,npbf, 
     $                       bi,bj)
      use shell
      
!     For a given shell quartet:
!      calculates either e-e repulsion with GTG
!      or e-(QM Nuc) attraction with GTG.
!     P1,P2,Pmat1,and Pmat2 are the shell pair quantities.
!     ise,jse,isn,jsn are the individual shell indices.
!     bi,bj are the bcoef's of GTG, and the latter is present only for
!     the nested gamma loop
      integer, intent(in):: ng1,ise,jse,isn,jsn,npebf,npbf,nesh,nnsh
      integer, intent(in):: AMPEB2C(npebf)
      double precision, intent(in)::AGEBFCC(npebf),AGNBFCC(npbf)
      double precision, intent(in):: P1,P2,gamm,Pmat1(3),Pmat2(3),bi
      double precision, intent(in),optional ::bj 
      type(eshell),intent(in) ::esh(nesh)
      type(pshell),intent(in) ::psh(nnsh)
      integer, intent(inout)::GAM_1(ng1)
      
!Local variables
      double precision :: c0,qP1P2,alpha,expnt
      double precision :: X_P1P2,Y_P1P2,Z_P1P2
      
      
         call g2_ee_prefactors(P1,P2,
     $                    Pmat1,Pmat2,gamm,
     $                    c0,qP1P2,alpha,expnt,
     $                    X_P1P2,Y_P1P2,Z_P1P2) 

