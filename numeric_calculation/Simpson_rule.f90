      subroutine fx(x,val)
      double precision x,val
      val=sqrt(1-x*x)
      return
      end

      subroutine simpson_rule(N,integral_val,x_min,x_max)
      integer N,i
      double precision integral_val,x_val,x_min,x_max,y(1000100),delta

      if (N.lt.1) then
        stop
      endif

      x_val=x_min
      delta=(x_max-x_min) / (2*N)
      i=1
      do 
        call fx(x_val,y(i))
        i=i+1
        x_val=x_val+delta
        if (x_val.ge.x_max) then
          exit
        endif
      enddo

      integral_val=0d0
      i=1
      do
        integral_val=integral_val+4*y(i)
        i=i+2
        if (i.gt.(2*N-1)) then
          exit
        endif
      enddo

      i=2
      do
        integral_val=integral_val+2*y(i)
        i=i+2
        if (i.gt.(2*N-1)) then
          exit
        endif
      enddo

      integral_val=integral_val+y(0)+y(2*N)
      integral_val=integral_val+delta/3d0
      write(*,*) "--------------------------------------"
      write(*,*) "N=",N,", integral_value=",integral_val
      write(*,*) "--------------------------------------"


      return
      end


