*CMZ :  4.01/05 24/03/2024  13.11.02  by  Michael Scheer
*CMZ : 00.00/02 22/04/99  17.27.52  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_fourier_spline(nx,X,Y,Nom,OM,ftcos,ftsin,IFAIL)

C---  CALCULATES INTEGRAL(F(x)*SIN(omega*x)) AND INTEGRAL(F(x)*COS(omega*x))

      IMPLICIT NONE

      INTEGER NX,NOM,IOM,Ix,IFAIL

      real*8, parameter :: pi1=3.141592653589793D0

      complex*16 :: fw,ci=(0.0d0,1.0d0),ex0oma,ex1oma,ex0x1oma,ft(nom),cn,cd,
     &  ex00,dexpom,dexpomx,ddexpomx

      REAL*8 X(NX),om(nom),yp(nx),ypp(nx),const,a0,y(nx),ftcos(nom),ftsin(nom),
     &  oma,oma2,oma3,oma4,x0,x02,x03,x1,x12,x13,y0,y1,ypp0,ypp1,dx,dom

      IFAIL=0
      const=1.0d0/sqrt(2.0d0*pi1)


      call util_coef_spline(nx,x,y,0.0d0,0.0d0,yp,ypp,ifail)
      if (ifail.ne.0) return

      dx=x(2)-x(1)
      dom=om(2)-om(1)

      dexpom=exp(ci*x(1)*dom)
      ex00=exp(ci*x(1)*om(1))
      dexpomx=exp(ci*dx*om(1))
      ddexpomx=exp(ci*dx*dom)

      DO IOM=1,nom

        oma=om(iom)

        if (abs(oma).lt.1.0d-9) then

          a0=0.0d0
          do ix=1,nx-1
            a0=a0
     &        +(X(ix+1)-X(ix))*0.5D0
     &        *(y(ix)+y(ix+1))
     &        -(X(ix+1)-X(ix))**3/24.D0
     &        *(ypp(ix)+ypp(ix+1))
          enddo

          ft(iom)=dcmplx(a0*const,0.0d0)

          ex00=ex00*dexpom
          dexpomx=dexpomx*ddexpomx
          cycle

        endif

        oma2=oma*oma
        oma3=oma2*oma
        oma4=oma3*oma

        fw=dcmplx(0.0d0,0.0d0)
        ex0oma=ex00

        do ix=1,nx-1

          x0=x(ix)
          x02=x0*x0
          x03=x02*x0
          x1=x(ix+1)
          x12=x1*x1
          x13=x12*x1
          y0=y(ix)
          y1=y(ix+1)
          ypp0=ypp(ix)
          ypp1=ypp(ix+1)

          ex1oma=ex0oma*dexpomx
          ex0x1oma=ex0oma*ex1oma

          cn=(-(ex0oma*(2.0d0*x12*ypp1*ci*oma2
     &      +x12*ypp0*ci*oma2-4.0d0*x1*x0*ypp1*ci*oma2-2.0d0*x1*x0*ypp0*ci*oma2-
     &      6.0d0*x1*y1*oma3+6.0d0*x1*ypp1*oma+2.0d0*x02*ypp1*ci*oma2+
     &      x02*ypp0*ci*oma2+6.0d0*x0*y1*oma3-6.0d0*x0*ypp1*oma+6.0d0*y1*ci
     &      *oma2-6.0d0*y0*ci*oma2-6.0d0*ypp1*ci+6.0d0*ypp0*ci)+
     &      ex1oma*(x12*ypp1*ci*oma2+2.0d0*
     &      x12*ypp0*ci*oma2-2.0d0*x1*x0*ypp1*ci*oma2-4.0d0*x1*x0*ypp0*ci*oma2+
     &      6.0d0*x1*y0*oma3-6.0d0*x1*ypp0*oma+
     &      x02*ypp1*ci*oma2+
     &      2.0d0*x02*ypp0*ci*oma2-6.0d0*x0*y0
     &      *oma3+6.0d0*x0*ypp0*oma-6.0d0*y1*ci*oma2+6.0d0*y0
     &      *ci*oma2+6.0d0*ypp1*ci-6.0d0*ypp0*ci))*ci)

          cd=6.0d0*ex0x1oma*(x1-x0)*oma4

          if (abs(cd).ge.1.0d-12) fw=fw+cn/cd

          ex0oma=ex0oma*dexpomx

        enddo !ix

        ft(iom)=fw*const

        ex00=ex00*dexpom
        dexpomx=dexpomx*ddexpomx

      enddo !oma

      ftcos=dreal(ft)
      ftsin=dimag(ft)

      return
      end
