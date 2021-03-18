      subroutine compute_inv_matrix(n,a,inv_a)
      integer n,i,j,k
      double precision a(100,100),tmp_a(100,100),inv_a(100,100),tmp
      !掃き出し法
      !参考
      !https://thira.plavox.info/blog/2008/06/_c.html
      !https://risalc.info/src/Gauss-Jordan-elimination.html
      
      !input:
      !n: integer 
      !a: (n,n) matrix
      !inv_a: (n,n) matrix
      
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
            !後から値が変化することを避けるため先にtmpに格納する
            do k=1,n
              tmp_a(j,k)=tmp_a(j,k)-tmp_a(i,k)*tmp
              inv_a(j,k)=inv_a(j,k)-inv_a(i,k)*tmp
            enddo
          endif
        enddo
      enddo
      return
      end
