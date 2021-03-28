      module sample
        implicit none
        double precision, parameter :: pi=acos(-1.0d0)

        contains
        subroutine compute_area(r)
            double precision r
            write(*,*) r*r*pi
        return
        end
      end module

      program  main
          use sample
          implicit none

          call compute_area(dble(1)) !=> 3.1415926535897931
          
      end program  main