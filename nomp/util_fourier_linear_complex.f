*CMZ :  4.01/05 28/03/2024  11.30.43  by  Michael Scheer
*CMZ : 00.00/02 22/04/99  17.27.52  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_fourier_linear_complex(nx,x,y,nom,om,ft,ifail)

c---  calculates integral(f(x)*sin(omega*x)) and integral(f(x)*cos(omega*x))

      implicit none

      integer nx,nom,iom,ix,ifail,ic

      real*8, parameter :: pi1=3.141592653589793d0

      complex*16 :: y(nx),fw,ci=(0.0d0,1.0d0),ex0oma,ft(nom),cn,dexpomx,dexpom,ex00,
     &  ddexpomx

      real*8 x(nx),om(nom),const,a0,yr(nx),oma,dom,dx

      ifail=0
      const=1.0d0/sqrt(2.0d0*pi1)

      dx=x(2)-x(1)
      dom=om(2)-om(1)

      dexpomx=exp(-ci*dx*om(1))
      ddexpomx=exp(-ci*dx*dom)

      do ic=1,2

        if (ic.eq.1) then
          yr(1:nx)=dreal(y(1:nx))
        else
          yr(1:nx)=dimag(y(1:nx))
        endif

        ex00=exp(-ci*x(1)*om(1))
        dexpom=exp(-ci*x(1)*dom)

        do iom=1,nom

          oma=om(iom)

          if (abs(oma).lt.1.0d-9) then
            a0=0.0d0
            do ix=1,nx-1
              a0=a0+(yr(ix)+yr(ix+1))/2.0d0*dx
            enddo
            if (ic.eq.1) then
              ft(iom)=dcmplx(a0*const,0.0d0)
            else
              ft(iom)=ft(iom)+dcmplx(0.0d0,a0*const)
            endif
            ex00=ex00*dexpom
            dexpomx=dexpomx*ddexpomx
            cycle
          endif

          ex0oma=ex00

          cn=yr(1)*ex0oma/2.0d0
          fw=cn*dx
          ex0oma=ex0oma*dexpomx

          do ix=2,nx-1
            cn=yr(ix)*ex0oma
            fw=fw+cn*dx
            ex0oma=ex0oma*dexpomx
          enddo !ix

          cn=yr(nx)*ex0oma/2.0d0
          fw=fw+cn*dx

          if (ic.eq.1) then
            ft(iom)=fw*const
          else
            ft(iom)=ft(iom)+ci*fw*const
          endif

          ex00=ex00*dexpom
          dexpomx=dexpomx*ddexpomx

        enddo !oma

      enddo !ic

      return
      end
