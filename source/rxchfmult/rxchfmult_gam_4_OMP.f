C======================================================================
      subroutine RXCHFmult_symm_gam4(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                     ip,jp,ii,i,jj,j,kk,k,ll,l,ans)
C
C======================================================================
C     ii  i   jj  j   kk  k   ll  l   ip  jp
C     ie1,je1,ie2,je2,ie3,je3,ie4,je4,ip1,jp1,ans)
C 
      implicit none
c     include 'param.h'
      double precision two,four
      parameter(two=2.0d+00,four=4.0d+00)
C input
      integer ne
      integer np
      integer ng2
      integer ngee
      double precision GAM_ee(ngee)
      double precision GAM_2s(ng2)
      integer ip
      integer jp
      integer ii,i
      integer jj,j
      integer kk,k
      integer ll,l

      double precision ans
      double precision x1,x2,x3,x4,x5,x6


      call RXCHFmult_i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,   ! 1 2 3 4
     x                                   ip,jp,ii,jj,kk,ll,i,j,k,l,x1)
      call RXCHFmult_i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                                   ip,jp,ii,jj,kk,ll,i,j,l,k,x2)
      call RXCHFmult_i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                                   ip,jp,ii,jj,kk,ll,i,k,j,l,x3)
      call RXCHFmult_i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                                   ip,jp,ii,jj,kk,ll,i,k,l,j,x4)
      call RXCHFmult_i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                                   ip,jp,ii,jj,kk,ll,i,l,j,k,x5)
      call RXCHFmult_i10(NE,np,ng2,ngee,GAM_2s,GAM_ee,
     x                                   ip,jp,ii,jj,kk,ll,i,l,k,j,x6)

      ans = x1
     x    - x2/TWO
     x    - x3/TWO
     x    + x4/FOUR
     x    + x5/FOUR
     x    - x6/TWO

      return
      end

C======================================================================
      SUBROUTINE RXCHFmult_i10(ne,np,ng2,ngee,GAM_2s,GAM_ee,
     x                     ip,jp,ie1,ie2,ie3,ie4,je1,je2,je3,je4,ans)
C======================================================================
C Physicist to chemist notation conversion
      implicit none
      integer ne
      integer np
      integer ng2
      integer ngee
c     double precision GAM_2s(np,np,ne,ne,ne,ne)
c     double precision GAM_ee(ne,ne,ne,ne)
CCWS-IO
      double precision GAM_2s(ng2)
c     double precision GAM_2s(1)
      double precision GAM_ee(ngee)

c     logical prin

      integer ip,jp
      integer ie1,je1
      integer ie2,je2
      integer ie3,je3
      integer ie4,je4
      double precision ans

      call RXCHFmult_calc_GAM_4_integral(ne,np,ng2,ngee,
     x                                   GAM_2s,GAM_ee,ip,jp,
     x                                   ie1,je1,ie2,je2,
     x                                   ie3,je3,ie4,je4,
     x                                   ans)
C

      return
      END


C======================================================================
      subroutine RXCHFmult_calc_GAM_4_integral(ne,np,ng2,ngee,
     x                                         GAM_2s,GAM_ee,ip,jp,
     x                                         ie1,je1,ie2,je2,
     x                                         ie3,je3,ie4,je4,
     x                                         ans)
C
C======================================================================
      implicit none

      double precision half,three
      parameter(half=0.5d+00,three=3.0d+00)

      integer ne
      integer np
      integer ng2
      integer ngee
      integer ip
      integer jp
      integer ie1
      integer ie2
      integer ie3
      integer ie4
      integer je1
      integer je2
      integer je3
      integer je4

      double precision GAM_2s(ng2)
      double precision GAM_ee(ngee)
      double precision ans

C Local variables
      double precision x1,x2,x3

      integer iee14
      integer iee13
      integer iee12

      integer igs23
      integer igs24
      integer igs34

      double precision ee14
      double precision ee13
      double precision ee12

      double precision gs23
      double precision gs24
      double precision gs34

      call pack_4D(ne,ne,ne,je4,ie4,je1,ie1,iee14)
      call pack_4D(ne,ne,ne,je3,ie3,je1,ie1,iee13)
      call pack_4D(ne,ne,ne,je2,ie2,je1,ie1,iee12)


      ee14=gam_ee(iee14) 
      ee13=gam_ee(iee13) 
      ee12=gam_ee(iee12) 

      call underflow(ee14)
      call underflow(ee13)
      call underflow(ee12)

      call index_GAM_2PK(ne,np,ip,jp,ie2,je2,ie3,je3,igs23)
      call index_GAM_2PK(ne,np,ip,jp,ie2,je2,ie4,je4,igs24)
      call index_GAM_2PK(ne,np,ip,jp,ie3,je3,ie4,je4,igs34)

      gs23=GAM_2s(igs23) 
      gs24=GAM_2s(igs24) 
      gs34=GAM_2s(igs34) 

      call underflow(gs23)
      call underflow(gs24)
      call underflow(gs34)

      x1=ee14*gs23
      x2=ee13*gs24
      x3=ee12*gs34

      ans=(x1+x2+x3)/three

      return
      end 

