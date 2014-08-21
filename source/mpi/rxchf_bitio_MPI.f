C=======================================================================
      subroutine bitpack(ip,jp,i1,j1,i2,j2,i3,j3,i4,j4,packp,packe)

C Packs two proton indices into packp as [ip jp] and
C packs eight electron indices into packe as [i1 j1 i2 j2 i3 j3 i4 j4]
C=======================================================================
      implicit none

! Input variables
      integer ip,jp
      integer i1,j1
      integer i2,j2
      integer i3,j3
      integer i4,j4

! Output variables
      integer(kind=4) packp
      integer(kind=8) packe

! Pack both proton indices into one 4-byte integer as [ip jp]
      packp=jp
      packp=packp+ISHFT(ip,8)

! Pack all electron indices into one 8-byte integer as [i1 j1 i2 j2 i3 j3 i4 j4]
      packe=j4
      packe=packe+ISHFT(i4,8)
      packe=packe+ISHFT(j3,16)
      packe=packe+ISHFT(i3,24)
      packe=packe+ISHFT(j2,32)
      packe=packe+ISHFT(i2,40)
      packe=packe+ISHFT(j1,48)
      packe=packe+ISHFT(i1,56)

      return
      end
C=======================================================================
      subroutine bitunpack(packp,packe,ip,jp,i1,j1,i2,j2,i3,j3,i4,j4)

C Unpacks two proton indices from packp as [ip jp] and
C unpacks eight electron indices from packe [i1 j1 i2 j2 i3 j3 i4 j4]
C=======================================================================
      implicit none

! Input variables
      integer(kind=4) packp
      integer(kind=8) packe

! Output variables
      integer ip,jp
      integer i1,j1
      integer i2,j2
      integer i3,j3
      integer i4,j4

! Unpack both proton indices from one 4-byte integer as [ip jp]
      jp=IAND(packp,255)
      ip=ISHFT(packp,-8)

! Unpack all electron indices from one 8-byte integer as [i1 j1 i2 j2 i3 j3 i4 j4]
      j4=IAND(packe,255)
      i4=IAND(ISHFT(packe,-8),255)
      j3=IAND(ISHFT(packe,-16),255)
      i3=IAND(ISHFT(packe,-24),255)
      j2=IAND(ISHFT(packe,-32),255)
      i2=IAND(ISHFT(packe,-40),255)
      j1=IAND(ISHFT(packe,-48),255)
      i1=ISHFT(packe,-56)

      return
      end
C=======================================================================
      subroutine writebitint(nproc,rank,n,namelen,fname,pinds,einds,arr)

C Writes GAM4 integrals to file as
C   entry 1: number of nontrivial integrals
C   entry 2: packed proton index array (in chunks of size nchunk)
C   entry 3: packed electron index array (in chunks of size nchunk)
C   entry 4: nontrivial integral array (in chunks of size nchunk)
C=======================================================================
      implicit none

! Input variables
      integer                nproc,rank
      integer                n,namelen
      character(len=namelen) fname
      integer(kind=4)        pinds(n)
      integer(kind=8)        einds(n)
      real*8                 arr(n)

! Local variables
      integer                            :: unitno,i
      integer, parameter                 :: nchunk=100000000

C Write to file
      unitno=200+rank
      open(unit=unitno,file=fname,form="unformatted")

! Write number of nontrivial integrals
      write(unitno) n

! Write packed proton indices
      do i=1,n/nchunk
        write(unitno) pinds((i-1)*nchunk+1:i*nchunk)
      end do
      if(mod(n,nchunk).ne.0) then
        write(unitno) pinds((n/nchunk)*nchunk+1:n)
      end if

! Write packed electron indices
      do i=1,n/nchunk
        write(unitno) einds((i-1)*nchunk+1:i*nchunk)
      end do
      if(mod(n,nchunk).ne.0) then
        write(unitno) einds((n/nchunk)*nchunk+1:n)
      end if

! Write integrals
      do i=1,n/nchunk
        write(unitno) arr((i-1)*nchunk+1:i*nchunk)
      end do
      if(mod(n,nchunk).ne.0) then
        write(unitno) arr((n/nchunk)*nchunk+1:n)
      end if

      close(unitno)

      return
      end

