      Program main
      implicit none
      integer N,i,j,k,p,q,stop_loop_flag
      double precision a(100,100)
      double precision g_mx(100,100),inv_g_mx(100,100)
      double precision gag_mx(100,100),ag_mx(100,100)
      read(*,*) N

      !参考
      !http://nkl.cc.u-tokyo.ac.jp/16n/Eigen.pdf
      !http://www.na.scitec.kobe-u.ac.jp/~yamamoto/lectures/appliedmathematics2004/chapter11.PDF
      do i=1,N
        read(*,*) (a(i,j),j=1,N)
      enddo

      do i=1,5000
        !ループを抜けるか判定するフラグ
        stop_loop_flag=0
        call update_mx_matrix(N,a,p,q,g_mx,stop_loop_flag,i)
        if (stop_loop_flag.eq.1) then 
          exit
        endif
        write(*,*) "---------------------------------"
        write(*,*) "i=",i
        write(*,*) "matrix a"
        do j=1,N
          write(*,*) (a(j,k),k=1,N)
        enddo
        write(*,*) "---------------------------------"
      enddo

      !非対角成分の微小量の項を0にする
      do i=1,N
        do j=1,N
          if (abs(a(i,j)).lt.1e-5) then
            a(i,j)=0.0d0
          endif
        enddo
      enddo

      write(*,*)"-------------------------------"
      write(*,*)"eigenvalue"
      do i=1,N
        write(*,*) (a(i,j),j=1,N)
      enddo
      write(*,*)"-------------------------------"

      write(*,*)"-------------------------------"
      write(*,*)"eigenvector"
      do i=1,N
        write(*,*) (g_mx(i,j),j=1,N)
      enddo
      write(*,*)"-------------------------------"
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
            inv_a(i,j)=1.d0
          else
            inv_a(i,j)=0.d0
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

      subroutine get_transform_mx(n,mx,t_mx)
      !転置行列を得る
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

      subroutine update_mx_matrix(n,mx,p,q,g_mx,stop_loop_flag,iteration)
      integer n,p,q,i,j,k,stop_loop_flag,iteration
      double precision max_val,mx(100,100)
      double precision theta,sin_val,cos_val,tan_val,threshold
      double precision a_pk,a_qk,a_kp,a_kq
      double precision g_mx(100,100),new_g_mx(100,100),tmp_g_mx(100,100)
      stop_loop_flag=0
      max_val=mx(1,2)
      p=1
      q=2
      !行列G
      do i=1,n
        do j=1,n
          if (i.eq.j) then 
            cycle
          endif
          if (abs(mx(i,j)).gt.max_val) then
            max_val=abs(mx(i,j))
            p=i
            q=j
          endif
        enddo
      enddo

      threshold=1e-16
      if (abs(mx(p,q)).le.threshold) then
        stop_loop_flag=1
        write(*,*) "p=",p,", q=",q,", abs mx(p,q)=",mx(p,q)
        !stop
      endif

      if (abs(mx(p,p)-mx(q,q)).lt.threshold) then
        theta=atan(1.d0)
      else 
        tan_val=2.d0*mx(p,q) / (mx(p,p)-mx(q,q))
        theta=0.5d0*atan(tan_val)
      endif
      
      !thetaの弧度法に注意
      !theta=theta*180/(4*atan(1d0))
      sin_val=sin(theta)
      cos_val=cos(theta)
      write(*,*) "theta=",theta

      do k=1,n
        a_pk=mx(p,k)
        a_qk=mx(q,k)
        mx(p,k)=a_pk*cos_val + a_qk*sin_val
        mx(q,k)=(-1)*a_pk*sin_val + a_qk*cos_val
      enddo

      do k=1,n
        a_kp=mx(k,p)
        a_kq=mx(k,q)
        mx(k,p)=a_kp*cos_val + a_kq*sin_val
        mx(k,q)=(-1)*a_kp*sin_val + a_kq*cos_val
      enddo

      !行列gの更新
      do i=1,n
        do j=1,n
          if (i.eq.j) then
            tmp_g_mx(i,j)=1.0d0
          else
            tmp_g_mx(i,j)=0.d0
          endif
        enddo
      enddo

      tmp_g_mx(p,p)=cos_val
      tmp_g_mx(q,q)=cos_val
      tmp_g_mx(p,q)=sin_val
      tmp_g_mx(q,p)=(-1)*sin_val

      if (iteration.eq.1) then
        g_mx=tmp_g_mx
      else
        call compute_mx_product(n,g_mx,tmp_g_mx,new_g_mx)
        g_mx=new_g_mx
      endif
      return 
      end

      subroutine get_theta_val(mx,p,q,theta)
      double precision mx(100,100),theta
      integer p,q
      theta = 0.5*atan((2)*mx(p,q) / (mx(p,p) - mx(q,q)))
      return
      end


