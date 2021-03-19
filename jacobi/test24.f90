      Program main
      implicit none
      integer N,i,j,k,p,q,stop_loop_flag
      double precision a(100,100),eigenvec(100),eigenval(100)
      double precision g_mx(100,100),transform_g(100,100)
      double precision ag_prod(100,100),gag_prod(100,100)
      double precision T(100,100),tmp_T(100,100),inv_T(100,100)
      double precision AT_prod(100,100),lambda_mx(100,100)
      read(*,*) N

      !参考
      !http://nkl.cc.u-tokyo.ac.jp/16n/Eigen.pdf
      !http://www.na.scitec.kobe-u.ac.jp/~yamamoto/lectures/appliedmathematics2004/chapter11.PDF
      do i=1,N
        read(*,*) (a(i,j),j=1,N)
      enddo

      do i=1,1
        call get_given_matrix(N,a,p,q,g_mx,stop_loop_flag)
        write(*,*) "---------------------------------"
        write(*,*) "i=",i
        write(*,*) "matrix a"
        do j=1,N
          write(*,*) (a(j,k),k=1,N)
        enddo
        call get_transform_mx(N,g_mx,transform_g)
        call compute_mx_product(N,a,g_mx,ag_prod)
        call compute_mx_product(N,transform_g,ag_prod,gag_prod)
        a=gag_prod
        write(*,*) "---------------------------------"
        
        if (i.eq.1) then
          T=gag_prod
        else 
          tmp_T=T
          call compute_mx_product(N,tmp_T,g_mx,T)
        endif

        !ループを抜けるか判定するフラグ
        stop_loop_flag=0
        !call do_break_loop(N,a,stop_loop_flag)
        if (stop_loop_flag.eq.1) then 
          exit
        endif
      enddo
      call compute_inv_matrix(N,T,inv_T)
      call compute_mx_product(N,a,T,AT_prod)
      call compute_mx_product(N,inv_T,AT_prod,lambda_mx)


      write(*,*) "-------------------"
      write(*,*) "T"
      do i=1,N
        write(*,*) (T(i,j),j=1,N)
      enddo
      write(*,*) "-------------------"


      write(*,*) "-------------------"
      write(*,*) "inv T"
      do i=1,N
        write(*,*) (inv_T(i,j),j=1,N)
      enddo
      write(*,*) "-------------------"

 
      write(*,*) "-------------------"
      write(*,*) "Optimized matrix"
      do i=1,N
        write(*,*) (lambda_mx(i,j),j=1,N)
      enddo
      write(*,*) "-------------------"

      stop
      end

      subroutine compute_inv_matrix(n,a,inv_a)
      integer n,i,j,k
      double precision a(100,100),tmp_a(100,100),inv_a(100,100),tmp
      !参考
      !https://thira.plavox.info/blog/2008/06/_c.html
      !https://risalc.info/src/Gauss-Jordan-elimination.html

      tmp_a=a
      !初期化
      do i=1,n
        do j=1,n
          if (i.eq.j) then
            inv_a(i,j)=1d0
          else
            inv_a(i,j)=0d0
          endif
        enddo
      enddo
        
      do i=1,n
        tmp=tmp_a(i,i)
        do j=1,n
          tmp_a(i,j)=tmp_a(i,j)/tmp
          inv_a(i,j)=inv_a(i,j)/tmp
        enddo

        do j=1,n
          if (i.ne.j) then
            tmp=tmp_a(j,i)
            !後から値が変化してしまうので先にtmpに格納する
            do k=1,n
              tmp_a(j,k)=tmp_a(j,k)-tmp_a(i,k)*tmp
              inv_a(j,k)=inv_a(j,k)-inv_a(i,k)*tmp
            enddo
          endif
        enddo
      enddo
      return
      end 



      subroutine do_break_loop(n,mx,stop_loop_flag)
      integer i,j,n,stop_loop_flag
      double precision max_val,min_val,threshold,mx(100,100)
      stop_loop_flag=0
      max_val=mx(1,1)
      min_val=mx(1,1)
      do i=1,n
        do j=1,n
          if (mx(i,j).gt.max_val) then
            max_val=mx(i,j)
          endif
          if (mx(i,j).lt.min_val) then
            min_val=mx(i,j)
          endif
        enddo
      enddo

      !閾値
      threshold=1e-16
      if ((max_val-min_val).lt.threshold) then 
        stop_loop_flag=1
      endif

      return
      end

      subroutine get_transform_mx(n,mx,t_mx)
      integer n,i,j
      double precision mx(100,100),t_mx(100,100)
      do i=1,n
        do j=1,n
          t_mx(i,j)=mx(j,i)
        enddo
      enddo

      return 
      end


      subroutine compute_mx_product(n,mx1,mx2,prod)
      integer n,i,j,k
      double precision mx1(100,100),mx2(100,100),prod(100,100)
      do i=1,n 
        do j=1,n
          prod(i,j)=0d0
          do k=1,n
            prod(i,j)=prod(i,j) + mx1(i,k)*mx2(k,j)
          enddo
        enddo 
      enddo

      return
      end

      subroutine get_given_matrix(n,mx,p,q,g_mx,stop_loop_flag)
      integer n,p,q,i,j,stop_loop_flag
      double precision max_val,mx(100,100),g_mx(100,100)
      double precision theta,threshold

      max_val=mx(1,2)
      p=1
      q=2
      do i=1,n
        do j=1,n
          if (i.eq.j) then 
            cycle
          endif
          if (mx(i,j).gt.max_val) then
            max_val=mx(i,j)
            p=i
            q=j
          endif
        enddo
      enddo

      threshold=1e-10
      if (abs(mx(p,p)-mx(q,q)).lt.threshold) then
        stop_loop_flag=1
      endif
      
      !given変換のためのthetaの値を得る
      call get_theta_val(mx,p,q,theta)

      !given変換のための行列の初期化
      do i=1,n
        do j=1,n
          if (i.eq.j) then 
            g_mx(i,j)=1.d0
          else
            g_mx(i,j)=0.d0
          endif
        enddo
      enddo

      !thetaの弧度法に注意
      !theta=theta*180/(4*atan(1d0))
      g_mx(p,p)=cos(theta)
      g_mx(q,q)=cos(theta)
      g_mx(p,q)=(-1)*sin(theta)
      g_mx(q,p)=sin(theta)
      write(*,*) "-------------------------------"
      write(*,*) "theta=",theta
      write(*,*) "-------------------------------"


      return 
      end

      subroutine get_theta_val(mx,p,q,theta)
      double precision mx(100,100),theta
      integer p,q
      theta = 0.5*atan((2)*mx(p,q) / (mx(p,p) - mx(q,q)))
      return
      end


