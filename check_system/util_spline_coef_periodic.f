*CMZ :          26/08/2020  14.16.55  by  Michael Scheer
*CMZ : 00.00/11 26/05/2011  12.50.27  by  Michael Scheer
*CMZ : 00.00/08 10/12/2010  14.03.25  by  Michael Scheer
*CMZ : 00.00/07 10/12/2010  13.21.53  by  Michael Scheer
*-- Author :    Michael Scheer   10/12/2010
      subroutine util_spline_coef_periodic(x,y,n,ypp,aa,bb,cc,c,cn,ifail)

c--- calculates spline coefficients for periodic function
c--- the interval must be closed, i.e. x(n)-x(1)=periodelength and y(n)=y(1)

c--   input:

c-       n: number of x,y-values, must be at least five
c-       x: array of x-values
c-       y: array of y-values

c--   ouput:

c-       ypp:   spline-coefficients
c-     ifail:   error status

c--   workingspace: aa(n),bb(n),cc(n),c(n),cn(n)

      implicit none

      integer n,j,n1,n2,ifail
      real*8  x(n),y(n),ypp(n),aa(n),bb(n),cc(n),c(n),cn(n),ymax

      if(n.lt.5) then
        ifail=-1
        return
      endif

      ymax=0.0d0
      do j=1,n
        if (abs(y(j)).gt.ymax) ymax=abs(y(j))
      enddo

      if(abs(y(n)-y(1))/ymax.gt.1.0d-9) then
        ifail=-2
        return
      else
        y(1)=(y(n)+y(1))/2.0d0
        y(n)=y(1)
      endif

      cn=0.0d0 !letzte Spalte der Matrix

c      y(n)=y(1)
c      ypp(n)=ypp(1)

c n=5
c      aa(1)*ypp(4)+bb(1)*ypp(1)+cc(1)*ypp(2)=c(1)
c      aa(2)*ypp(1)+bb(2)*ypp(2)+cc(2)*ypp(3)=c(2)
c      aa(3)*ypp(2)+bb(3)*ypp(3)+cc(3)*ypp(4)=c(3)
c      aa(4)*ypp(3)+bb(4)*ypp(4)+cc(4)*ypp(1)=c(4)

      n1=n-1
      n2=n-2

      aa(1)=(x(n)-x(n1))/6.d0
      bb(1)=(x(2)-x(1)+(x(n)-x(n1)))/3.d0
      cc(1)=(x(2)-x(1))/6.d0
      c(1)=(y(2)-y(1))/(x(2)-x(1))
     &    -(y(n)-y(n1))/(x(n)-x(n1))

      do j=2,n1
          aa(j)=(x(j  )-x(j-1))/6.d0
          bb(j)=(x(j+1)-x(j-1))/3.d0
          cc(j)=(x(j+1)-x(j  ))/6.d0
          c(j)=(y(j+1)-y(j  ))/(x(j+1)-x(j  ))
     &          -(y(j  )-y(j-1))/(x(j  )-x(j-1))
      enddo !j

c Auf Dreiecksmatrix bringen

c Oberste Zeile

      cc(1)=cc(1)/bb(1)
      aa(1)=aa(1)/bb(1)
      c(1) = c(1)/bb(1)
      bb(1)=1.0d0

c 2. Zeile, d.h. die erste reguläre, cn(j) ist die letzte Spalte der Matrix

      bb(2)=bb(2)/aa(2)
      cc(2)=cc(2)/aa(2)
      c(2) = c(2)/aa(2)
      aa(2)=1.0d0

      aa(2)=0.0d0
      bb(2)=bb(2)-cc(1)
      cn(2)=-aa(1)
      c(2) = c(2)-c(1)

      cc(2)=cc(2)/bb(2)
      cn(2)=cn(2)/bb(2)
      c(2) = c(2)/bb(2)

      bb(2)=1.0d0

c Nun die hoeheren bis n-3

      do j=3,n-3

        bb(j)=bb(j)/aa(j)
        cc(j)=cc(j)/aa(j)
        c(j) = c(j)/aa(j)

        aa(j)=0.0d0
        bb(j)=bb(j)-cc(j-1)
        cn(j)=-cn(j-1)
        c(j) = c(j)-c(j-1)

        cc(j)=cc(j)/bb(j)
        cn(j)=cn(j)/bb(j)
        c(j) = c(j)/bb(j)

        bb(j)=1.0d0

      enddo

c vorletzte zeile

      bb(n2)=bb(n2)/aa(n2)
      cc(n2)=cc(n2)/aa(n2)
      c(n2) = c(n2)/aa(n2)
      aa(n2)=1.0d0

      aa(n2)=0.0d0
      bb(n2)=bb(n2)-cc(n2-1)
      cc(n2)=cc(n2)-cn(n2-1)
      c(n2) = c(n2)-c(n2-1)

      cc(n2)=cc(n2)/bb(n2)
      c(n2) = c(n2)/bb(n2)

      bb(n2)=1.0d0

c Letzte Zeile, benutze ypp als Puffer

      ypp=0.0d0
      ypp(n2)=aa(n1)/cc(n1)
      ypp(n1)=bb(n1)/cc(n1)
      c(n1)=c(n1)/cc(n1)
      cc(n1)=1.0d0
      ypp(1)=cc(n1)

c Oberste Zeile abziehen

      ypp(1)=0.0d0
      ypp(2)=-cc(1)
      ypp(n1)=ypp(n1)-aa(1)
      c(n1)=c(n1)-c(1)

      do j=2,n2

        c(n1)=c(n1)/ypp(j)
        ypp=ypp/ypp(j)

        ypp(j)=ypp(j)-bb(j)
        ypp(j+1)=ypp(j+1)-cc(j)
        ypp(n1)=ypp(n1)-cn(j)
        c(n1)=c(n1)-c(j)

      enddo

      c(n1)=c(n1)/ypp(n1)

c Ernten

      ypp(n1)=c(n1)

c Vorletzte Zeile

      bb(n2)=bb(n2)/cc(n2)
      c(n2)=c(n2)/cc(n2)
      cc(n2)=1.0d0

      c(n2)=c(n2)-c(n1)
      cc(n2)=0.0d0

      c(n2)=c(n2)/bb(n2)
      bb(n2)=1.0d0

      ypp(n2)=c(n2)

c Letzte Spale nullen und normieren, letzte und vorletzte sind bereits fertig

      bb(1)=bb(1)/aa(1)
      cc(1)=cc(1)/aa(1)
      c(1)=c(1)/aa(1)
      aa(1)=1.0d0

      c(1)=c(1)-c(n1)
      aa(1)=0.0d0

      cc(1)=cc(1)/bb(1)
      c(1)=c(1)/bb(1)
      bb(1)=1.0d0

c Reguläre Zeilen
      do j=2,n-3

        bb(j)=bb(j)/cn(j)
        cc(j)=cc(j)/cn(j)
        c(j)=c(j)/cn(j)
        cn(j)=1.0d0

        c(j)=c(j)-c(n1)
        cn(j)=0.0d0

        cc(j)=cc(j)/bb(j)
        c(j)=c(j)/bb(j)
        bb(j)=1.0d0

      enddo

      do j=n-3,2,-1
         ypp(j)=c(j)-cc(j)*ypp(j+1)
      enddo

      ypp(1)=(c(1)-cc(1)*ypp(2)-aa(1)*ypp(n1))/bb(1)
      ypp(n)=ypp(1)

      ifail=0

      return
      end
