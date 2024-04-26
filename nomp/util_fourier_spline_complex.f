*CMZ :  4.01/05 28/03/2024  08.39.34  by  Michael Scheer
*CMZ : 00.00/02 22/04/99  17.27.52  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_fourier_spline_complex(nx,X,Y,Nom,OM,ft,IFAIL)

C---  CALCULATES INTEGRAL(F(x)*SIN(omega*x)) AND INTEGRAL(F(x)*COS(omega*x))

      IMPLICIT NONE

      INTEGER NX,NOM,IOM,Ix,IFAIL,ic

      real*8, parameter :: pi1=3.141592653589793D0

      complex*16 :: y(nx),fw,ci=(0.0d0,1.0d0),ex0oma,ex1oma,ex0x1oma,ft(nom),cn,cd,
     &  ex00,dexpom,dexpomx,ddexpomx

      REAL*8 X(NX),om(nom),yp(nx),ypp(nx),const,a0,yr(nx),
     &  oma,oma2,oma3,oma4,x0,x02,x03,x1,x12,x13,y0,y1,ypp0,ypp1,dx,dom

      IFAIL=0
      const=1.0d0/sqrt(2.0d0*pi1)

      dx=x(2)-x(1)
      dom=om(2)-om(1)

      dexpom=exp(ci*x(1)*dom)
      ddexpomx=exp(ci*dx*dom)

      do ic=1,2

        if (ic.eq.1) then
          yr(1:nx)=dreal(y(1:nx))
        else
          yr(1:nx)=dimag(y(1:nx))
        endif

        call util_coef_spline(nx,x,yr,0.0d0,0.0d0,yp,ypp,ifail)
        if (ifail.ne.0) return

        ex00=exp(ci*x(1)*om(1))
        dexpomx=exp(ci*dx*om(1))

        DO IOM=1,nom

          oma=om(iom)

          if (abs(oma).lt.1.0d-9) then
            a0=0.0d0
            do ix=1,nx-1
              a0=a0
     &          +(X(ix+1)-X(ix))*0.5D0
     &          *(Yr(ix)+Yr(ix+1))
     &          -(X(ix+1)-X(ix))**3/24.D0
     &          *(ypp(ix)+ypp(ix+1))
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
            y0=yr(ix)
            y1=yr(ix+1)
            ypp0=ypp(ix)
            ypp1=ypp(ix+1)

            ex1oma=ex0oma*dexpomx
            ex0x1oma=ex0oma*ex1oma

            cn=(-(ex0oma*(2.0d0*x12*ypp1*ci*oma2
     &        +x12*ypp0*ci*oma2-4.0d0*x1*x0*ypp1*ci*oma2-2.0d0*x1*x0*ypp0*ci*oma2-
     &        6.0d0*x1*y1*oma3+6.0d0*x1*ypp1*oma+2.0d0*x02*ypp1*ci*oma2+
     &        x02*ypp0*ci*oma2+6.0d0*x0*y1*oma3-6.0d0*x0*ypp1*oma+6.0d0*y1*ci
     &        *oma2-6.0d0*y0*ci*oma2-6.0d0*ypp1*ci+6.0d0*ypp0*ci)+
     &        ex1oma*(x12*ypp1*ci*oma2+2.0d0*
     &        x12*ypp0*ci*oma2-2.0d0*x1*x0*ypp1*ci*oma2-4.0d0*x1*x0*ypp0*ci*oma2+
     &        6.0d0*x1*y0*oma3-6.0d0*x1*ypp0*oma+
     &        x02*ypp1*ci*oma2+
     &        2.0d0*x02*ypp0*ci*oma2-6.0d0*x0*y0
     &        *oma3+6.0d0*x0*ypp0*oma-6.0d0*y1*ci*oma2+6.0d0*y0
     &        *ci*oma2+6.0d0*ypp1*ci-6.0d0*ypp0*ci))*ci)

            cd=6.0d0*ex0x1oma*(x1-x0)*oma4

            if (abs(cd).ge.1.0d-12) fw=fw+cn/cd

            ex0oma=ex0oma*dexpomx

          enddo !ix

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
