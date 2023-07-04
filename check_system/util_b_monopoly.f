*CMZ :          23/10/2017  09.35.02  by  Michael Scheer
*-- Author : Michael Scheer

c +PATCH,//UTIL/FOR
c +DECK,util_b_monopoly.

      subroutine util_b_monopoly(nq,qgrid,x,y,z,bx,by,bz,ifail)

      implicit none

c Nq charges qi=qgrid(4,i) at xi=qgrid(1), yi=qgrid(2,i), zi=qgrid(3,i)
c
c ri = (x-xi,y-yi,z-zi)
c
c Bx(x,y,z) = Sum_n_1_N-1 {  Qi / ri**3 * (x-xi) } + Qn / rn**3 * (x-xn)
c By(x,y,z) = Sum_n_1_N-1 {  Qi / ri**3 * (y-yi) } + Qn / rn**3 * (y-yn)
c Bz(x,y,z) = Sum_n_1_N-1 {  Qi / ri**3 * (z-zi) } + Qn / rn**3 * (z-zn)
c
c Qn = Sum_n_1_N-1 { -Qi }
c

      double precision qgrid(4,nq),x,y,z,bx,by,bz,qr3,dx,dy,dz

      integer ifail, nq, iq

      bx=0.0d0
      by=0.0d0
      bz=0.0d0
      ifail=0

      do iq=1,nq
        dx=x-qgrid(1,iq)
        dy=y-qgrid(2,iq)
        dz=z-qgrid(3,iq)
        qr3=qgrid(4,iq)/sqrt(dx**2+dy**2+dz**2)**3
        if (qr3.eq.0.0d0) then
          ifail=1
          cycle
        endif
        bx=bx+qr3*dx
        by=by+qr3*dy
        bz=bz+qr3*dz
      enddo

      end
