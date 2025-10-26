*CMZ :          26/09/2025  13.16.41  by  Michael Scheer
*-- Author :    Michael Scheer   26/09/2025
      subroutine util_e_to_stokes(e,specnor,s)

      implicit none

      complex(8), dimension(4,3), parameter ::
     &  vstokes=reshape([
     &  ( 0.0000000000000000d0,  0.0000000000000000d0),
     &  ( 0.0000000000000000d0,  0.0000000000000000d0),
     &  ( 0.0000000000000000d0,  0.0000000000000000d0),
     &  ( 0.0000000000000000d0,  0.0000000000000000d0),
     &  ( 0.0000000000000000d0,  0.0000000000000000d0),
     &  ( 0.0000000000000000d0, -0.70710678118654746d0),
     &  ( 0.0000000000000000d0, -0.70710678118654746d0),
     &  ( 0.70710678118654746d0, 0.0000000000000000d0),
     &  (-0.70710678118654746d0,-0.70710678118654746d0),
     &  ( 0.70710678118654746d0, 0.0000000000000000d0),
     &  (-0.70710678118654746d0, 0.0000000000000000d0),
     &  (-0.70710678118654746d0, 0.0000000000000000d0)
     &  ],[4,3])


      complex(8) e(3),apolh,apolr,apoll,apol45
      real(8) s(4),s1,s2,s3,s4,specnor

      apolh=
     &  e(1)*conjg(vstokes(1,1))
     &  +e(2)*conjg(vstokes(1,2))
     &  +e(3)*conjg(vstokes(1,3))

      apolr=
     &  e(1)*conjg(vstokes(2,1))
     &  +e(2)*conjg(vstokes(2,2))
     &  +e(3)*conjg(vstokes(2,3))

      apoll=
     &  e(1)*conjg(vstokes(3,1))
     &  +e(2)*conjg(vstokes(3,2))
     &  +e(3)*conjg(vstokes(3,3))

      apol45=
     &  e(1)*conjg(vstokes(4,1))
     &  +e(2)*conjg(vstokes(4,2))
     &  +e(3)*conjg(vstokes(4,3))

      s1=dreal(apolr*conjg(apolr)+apoll*conjg(apoll))
      s2=dreal(-s1+2.0d0*apolh*conjg(apolh))
      s3=dreal(2.0d0*apol45*conjg(apol45)-s1)
      s4=dreal(apolr*conjg(apolr)-apoll*conjg(apoll))

      s=[s1,s2,s3,s4]*specnor

      end
