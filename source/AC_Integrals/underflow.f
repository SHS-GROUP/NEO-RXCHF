C**********************************************************
      SUBROUTINE underflow(x)
C**********************************************************
      implicit none
      double precision x
      double precision TINY
      parameter (TINY = 1.0d-100)

      if(dabs(x) .lt. TINY) x = 0.0d0

      END

