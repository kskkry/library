        subroutine least_square_method(n,a,b,x,y)
        integer n
        double precision a,b,x(100),y(100)
        double precision x_mean,y_mean
        double precision x_sum,y_sum,xx_sum,xy_sum

        x_sum=0.0d0
        y_sum=0.0d0
        xx_sum=0.0d0
        xy_sum=0.0d0

        do i=1,n
          x_sum=x_sum+x(i)
          y_sum=y_sum+y(i)
          xx_sum=xx_sum + x(i)*x(i)
          xy_sum=xy_sum + x(i)*y(i)
        enddo

        x_mean=x_sum / dble(n)
        y_mean=y_sum / dble(n)

        a=(n*xy_sum - x_sum*y_sum) / (n*xx_sum - x_sum**2)
        b=(xx_sum*y_sum - xy_sum*x_sum) / (n*xx_sum - x_sum**2)
        return
        end

        subroutine compute_correlation_coefficient(n,x,y,r_val)
          integer n,i
          double precision x(100),y(100),r_val
          double precision x_sum,y_sum,x_avg,y_avg
          !共分散とx,yの分散
          double precision x_var,y_var,xy_var

          !compute mean of x,y
          x_sum=0.0d0
          y_sum=0.0d0
          do i=1,n
            x_sum = x_sum+x(i)
            y_sum = y_sum+y(i)
          enddo
          x_avg = x_sum / dble(N)
          y_avg = y_sum / dble(N)

          !compute var of x,y,xy
          x_var=0.0d0
          y_var=0.0d0
          xy_var=0.0d0
          do i=1,n
            x_var=x_var + (x(i)-x_avg)**2 / dble(N)
            y_var=y_var + (y(i)-y_avg)**2 / dble(N)
            xy_var=xy_var + (x(i)-x_avg)*(y(i)-y_avg) / dble(N)
          enddo

          r_val = xy_var / sqrt(x_var*y_var)
        return 
        end


