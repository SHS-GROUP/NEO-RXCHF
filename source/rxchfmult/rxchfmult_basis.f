C======================================================================
      subroutine RXCHFmult_GAM_2PK(nebf,nebf2,npbf,
     x                             ip,jp,
     x                             ie1,je1,
     x                             ie2,je2,ia)
C
C Analog of index_GAM_2PK for different special electron basis
C======================================================================
      implicit none

      integer ia
      integer nebf,nebf2
      integer npbf

      integer ip
      integer jp
      integer ie1
      integer je1
      integer ie2
      integer je2

C GAM_2(je2,ie2,je1,ie1,jp,ip)
C   ie1,je1=1,...,nebf
C   ie2,je2=1,...,nebf2
C     ip,jp=1,...,npbf

C  Array decalred as A(N1,N2,...)
C then A(i1,i2,...) = 
C    the (i1-1) + N1*((i2-1) + N2*(...))th element

        ia = 1 + (je2-1) + 
     x   nebf2*( (ie2-1) + 
     x   nebf2*( (je1-1) + 
     x    nebf*( (ie1-1) + 
     x    nebf*( (jp -1) +
     x    npbf*(  ip -1) ) ) ) )

      return
      end
C======================================================================
      subroutine RXCHFmult_GAM_3PK(nebf,nebf2,npbf,
     x                             ip,jp,
     x                             ie1,je1,
     x                             ie2,je2,
     x                             ie3,je3,ia)
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
C   ie1,je1=1,...,nebf
C   ie2,je2=1,...,nebf2
C   ie3,je3=1,...,nebf2
C     ip,jp=1,...,npbf

C  Array decalred as A(N1,N2,...)
C then A(i1,i2,...) = 
C    the (i1-1) + N1*((i2-1) + N2*(...))th element

      ia = 1 + (je3-1) + 
     x nebf2*( (ie3-1) + 
     x nebf2*( (je2-1) + 
     x nebf2*( (ie2-1) + 
     x nebf2*( (je1-1) + 
     x  nebf*( (ie1-1) + 
     x  nebf*( (jp -1) +
     x  npbf*(  ip -1) ) ) ) ) ) )

      return
      end
C======================================================================
      subroutine RXCHFmult_GAM_4PK(nebf,nebf2,npbf,
     x                             ip,jp,
     x                             ie1,je1,
     x                             ie2,je2,
     x                             ie3,je3,
     x                             ie4,je4,ia)
C
C GAM_4(je4,ie4,je3,ie3,je2,ie2,je1,ie1,jp,ip)
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
C   ie1,je1=1,...,nebf
C   ie2,je2=1,...,nebf2
C   ie3,je3=1,...,nebf2
C   ie4,je4=1,...,nebf2
C     ip,jp=1,...,npbf

C  Array decalred as A(N1,N2,...)
C then A(i1,i2,...) = 
C    the (i1-1) + N1*((i2-1) + N2*(...))th element

      ia = 1 + (je4-1) + 
     x nebf2*( (ie4-1) + 
     x nebf2*( (je3-1) + 
     x nebf2*( (ie3-1) + 
     x nebf2*( (je2-1) + 
     x nebf2*( (ie2-1) + 
     x nebf2*( (je1-1) + 
     x  nebf*( (ie1-1) + 
     x  nebf*( (jp -1) +
     x  npbf*(  ip -1) ) ) ) ) ) ) ) )

      return
      end

