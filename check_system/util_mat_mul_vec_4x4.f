*CMZ :          28/12/2021  10.32.29  by  Michael Scheer
*CMZ : 00.00/02 21/07/2004  15.43.47  by  Michael Scheer
*-- Author :    Michael Scheer   21/07/2004
      subroutine util_mat_mul_vec_4x4(a,vin,w)

      implicit none

      double precision a(4,4),v(4),w(4),vin(4)

      v=vin

      w(1)=a(1,1)*v(1)+a(1,2)*v(2)+a(1,3)*v(3)+a(1,4)*v(4)
      w(2)=a(2,1)*v(1)+a(2,2)*v(2)+a(2,3)*v(3)+a(2,4)*v(4)
      w(3)=a(3,1)*v(1)+a(3,2)*v(2)+a(3,3)*v(3)+a(3,4)*v(4)
      w(4)=a(4,1)*v(1)+a(4,2)*v(2)+a(4,3)*v(3)+a(4,4)*v(4)

      return
      end
