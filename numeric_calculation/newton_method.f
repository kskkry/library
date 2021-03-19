!!下のfxを改変する
      subroutine newton_method(x_min,x_max) 
      double precision x_left,x_right,x_left_val,x_right_val
      double precision x_min,x_max
      double precision tmp,tmp_val,x_mid,x_mid_val

      x_left=x_min
      x_right=x_max
      call fx(x_left,x_left_val)
      call fx(x_right,x_right_val)

      if((x_left_val.ge.0).and.(x_right_val.ge.0)) then
        write(*,*) "RANGE ERROR"
      endif
      if((x_left_val.le.0).and.(x_right_val.le.0)) then
        write(*,*) "RANGE ERROR"
      endif

      !二分法の処理
      !fx<0のときのxをx_leftとしfx>0のときのxをx_rightとする
      if ((x_left_val.ge.0).and.(x_right_val.le.0)) then
        tmp=x_left
        x_left=x_right
        x_right=tmp
      endif

      do 
        x_mid=(x_left+x_right)/2
        call fx(x_left,x_left_val)
        call fx(x_right,x_right_val)
        call fx(x_mid,x_mid_val)

        if (abs(x_left-x_right)<1e-10) then
          write(*,*) "optimal value=",x_mid
          exit
        endif

        if (x_mid_val.ge.0) then
          x_right=x_mid
        else
          x_left=x_mid
        endif
      enddo
      stop
      end
      

      subroutine fx(x,val)
      double precision x,val
      val=x**3 - 16*(x**2) + 65*x -50
      return
      end