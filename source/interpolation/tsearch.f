        subroutine tsearch(t,ia,delta)
!search for the x index in the F table
! t: Fm(t) input t value
! ix: closest (0.1 interval) value for t, which will be used in ftable(m,ia)
! delta: t-x(ia) delta_x in the Taylor expansion
        implicit none
        double precision, intent(in) :: t 
        integer, intent(out) :: ia
        double precision, intent(out) :: delta 

        double precision :: xk,half,xk1
        double precision, parameter::eps=1.0d-15
        integer :: jx,k
! jx and ia are the table indices.

        delta=-1.0d0

        jx=int(t)*10
        half=int(t)+0.5d0
        if(t < half)then
          do k=0,4
            xk=int(t)+dble(k)*0.1d0
            xk1=xk+0.1d0
            if(abs(t-xk) < eps)then
              ia=jx+k
              delta=0.0d0
              exit
            elseif(abs(t-xk1) < eps)then
              ia=jx+k+1
              delta=0.0d0
              exit
            elseif(t>xk .and. t<xk1 )then
              ia=jx+k 
              exit
            endif
          enddo
        elseif (t > half .or. abs(t-half) <eps)then
          do k=5,9
            xk=int(t)+dble(k)*0.1d0
            xk1=xk+0.1d0
            if(abs(t-xk) < eps)then
              ia=jx+k
              delta=0.0d0
              exit
            elseif(abs(t-xk1) < eps)then
              ia=jx+k+1
              delta=0.0d0
              exit
            elseif(t>xk .and. t<xk1 )then
              ia=jx+k 
              exit
            endif
          enddo
        endif

        if(delta <0.0d0)delta=t-ia*0.1d0
                
        end
