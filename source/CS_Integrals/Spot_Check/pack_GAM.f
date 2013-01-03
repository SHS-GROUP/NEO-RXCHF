C GAM_1(je1,ie1,jp,ip)
C GAM_2(je2,ie2,je1,ie1,jp,ip)
C GAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
C GAM_4(je4,ie4,je3,ie3,je2,ie2,je1,ie1,jp,ip)
C
C  Array decalred as A(N1,N2,...)
C then A(i1,i2,...) = 
C    the (i1-1) + N1*((i2-1) + N2*(...))th element

C======================================================================
      subroutine pack_2D(n1,i1,i2,ia)
C
C
C======================================================================

      integer ia
      integer n1
      integer i1,i2

      ia = 1 + (i1-1) + 
     x    n1*( (i2-1) ) 


      return
      end


C======================================================================
      subroutine pack_4D(n1,n2,n3,
     x                   i1,i2,i3,i4,ia)
C
c     subroutine index_GAM_1PK(nebf,npbf,
c    x                         ip,jp,
c    x                         ie1,je1,ia)
C
C GAM_ee:
c     call pack_4D(nebf,nebf,nebf,
c    x             je2,ie2,je1,ie1,ia)
C GAM_ep/GAM_1:
c     call pack_4D(nebf,nebf,npbf,
c    x             je1,ie1,jp,ip,ia)
C
C   i1  i2  i3  i4
C   je1 ie1 jp  ip
C   je2 ie2 je1 ie1
C
C GAM_1(je1,ie1,jp,ip)    |  ==>  array(i1,i2,i3,i4)
C GAM_ee(je2,ie2,je1,ie1) |
C
c     ia = 1 + (je1-1) + 
c    x  nebf*( (ie1-1) + 
c    x  nebf*( (jp -1) +
c    x  npbf*(  ip -1) ) )
C
C======================================================================

      integer ia

      integer n1,n2,n3
      integer i1,i2,i3,i4

      ia = 1 + (i1-1) + 
     x    n1*( (i2-1) + 
     x    n2*( (i3-1) +
     x    n3*(  i4-1) ) )


      return
      end

C======================================================================
      subroutine index_GAM_2PK(nebf,npbf,
     x                         ip,jp,
     x                         ie1,je1,
     x                         ie2,je2,ia)
C
C======================================================================

      integer ia
      integer nebf
      integer npbf

      integer ip
      integer jp
      integer ie1
      integer je1
      integer ie2
      integer je2

C GAM_2(je2,ie2,je1,ie1,jp,ip)

C  Array decalred as A(N1,N2,...)
C then A(i1,i2,...) = 
C    the (i1-1) + N1*((i2-1) + N2*(...))th element

      ia = 1 + (je2-1) + 
     x  nebf*( (ie2-1) + 
     x  nebf*( (je1-1) + 
     x  nebf*( (ie1-1) + 
     x  nebf*( (jp -1) +
     x  npbf*(  ip -1) ) ) ) )


      return
      end

C======================================================================
      subroutine index_GAM_3PK(nebf,npbf,
     x                         ip,jp,
     x                         ie1,je1,
     x                         ie2,je2,
     x                         ie3,je3,ia)
C
C GAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)
C======================================================================

      integer ia
      integer nebf
      integer npbf

      integer ip
      integer jp
      integer ie1
      integer ie2
      integer ie3
      integer je1
      integer je2
      integer je3

C GAM_3(je3,ie3,je2,ie2,je1,ie1,jp,ip)

C  Array decalred as A(N1,N2,...)
C then A(i1,i2,...) = 
C    the (i1-1) + N1*((i2-1) + N2*(...))th element

      ia = 1 + (je3-1) + 
     x  nebf*( (ie3-1) + 
     x  nebf*( (je2-1) + 
     x  nebf*( (ie2-1) + 
     x  nebf*( (je1-1) + 
     x  nebf*( (ie1-1) + 
     x  nebf*( (jp -1) +
     x  npbf*(  ip -1) ) ) ) ) ) )


      return
      end

C======================================================================
      subroutine index_GAM_4PK(nebf,npbf,
     x                         ip,jp,
     x                         ie1,je1,
     x                         ie2,je2,
     x                         ie3,je3,
     x                         ie4,je4,ia)
C
C======================================================================

      integer ia
      integer nebf
      integer npbf

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

C GAM_4(je4,ie4,je3,ie3,je2,ie2,je1,ie1,jp,ip)

C  Array decalred as A(N1,N2,...)
C then A(i1,i2,...) = 
C    the (i1-1) + N1*((i2-1) + N2*(...))th element

      ia = 1 + (je4-1) + 
     x  nebf*( (ie4-1) + 
     x  nebf*( (je3-1) + 
     x  nebf*( (ie3-1) + 
     x  nebf*( (je2-1) + 
     x  nebf*( (ie2-1) + 
     x  nebf*( (je1-1) + 
     x  nebf*( (ie1-1) + 
     x  nebf*( (jp -1) +
     x  npbf*(  ip -1) ) ) ) ) ) ) ) )


      return
      end

C======================================================================
      subroutine index_GAM_4PK2(nebf,npbf,
     x                          ip,jp,
     x                          ie1,je1,
     x                          ie2,je2,
     x                          ie3,je3,
     x                          ie4,je4,ia_8)
C
C======================================================================

      integer*8 ia_8
      integer nebf
      integer npbf

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

C GAM_4(je4,ie4,je3,ie3,je2,ie2,je1,ie1,jp,ip)

C  Array decalred as A(N1,N2,...)
C then A(i1,i2,...) = 
C    the (i1-1) + N1*((i2-1) + N2*(...))th element

      ia_8 = 1 + (je4-1) + 
     x  nebf*( (ie4-1) + 
     x  nebf*( (je3-1) + 
     x  nebf*( (ie3-1) + 
     x  nebf*( (je2-1) + 
     x  nebf*( (ie2-1) + 
     x  nebf*( (je1-1) + 
     x  nebf*( (ie1-1) + 
     x  nebf*( (jp -1) +
     x  npbf*(  ip -1) ) ) ) ) ) ) ) )


      return
      end

