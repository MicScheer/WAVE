*CMZ : 00.00/20 11/10/2016  15.15.18  by  Michael Scheer
*CMZ : 00.00/07 29/06/2010  15.40.05  by  Michael Scheer
*CMZ : 00.00/02 29/07/99  16.24.48  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_LINEAR_FIT_USER(NARG,NFUN,NPAR,NDIMPOI,IPOI,FUNDATA,T)

c +PATCH,//UTIL/FOR
c +DECK,util_fit_2d-parabola_user.

C     THIS SUBROUTINE IS CALLED BY UTIL_LINEAR_FIT:
C
C     THE USER HAS TO PROVIDE HERE THE TERMS T(IFUN,IPAR)
C     WHERE T(IFUN,IPAR) IS THE DERIVATIVE OF THE iTH FITTED
C     FUNCTION TO THE iparTH FITPARAMETER FOR THE DATA-POINT ipoi.
C
C     EXAMPLE:
C
C     BX(X,Y,Z) IS FIRST FUNCTION OF ARGUMENTS X,Y,Z
C     BY(X,Y,Z) IS SECOND FUNCTION OF ARGUMENTS X,Y,Z
C     BZ(X,Y,Z) IS THIRD FUNCTION OF ARGUMENTS X,Y,Z
C     I.E. NARG=3 AND NFUN=3
C
C     P1...P4 ARE PARAMETER TO BE FITTED SIMULTANIOUSLY FOR THE
C     FUNCTIONS BX,BY AND BZ
C     I.E. NPAR=4
C
C     BX=P1*X(I)*Y(I)*SIN(K*Z(I))+P2*X(I)*Y(I)*SIN(3*K*Z(I))
C     BY=P3*X(I)*Y(I)*COS(K*Z(I))+P4*X(I)*Y(I)*COS(3*K*Z(I))
C       BZ=(3*COS(K*Z)*P1*y+COS(3*K*Z)*P2*y
C          -3*SIN(K*Z)*P3*x-SIN(3*K*Z)*P4*x)/(3*K)
C
C
C          T(1,1)=dBX/dP1=X(I)*Y(I)*SIN(K*Z(I))
C          T(1,2)=dBX/dP2=X(I)*Y(I)*SIN(3*K*Z(I))
C          T(1,3)=dBX/dP3=0.
C          T(1,4)=dBX/dP4=0.
C          T(2,1)=dBY/dP1=0.
C          T(2,1)=dBY/dP2=0.
C          T(2,3)=dBY/dP3=X(I)*Y(I)*COS(K*Z(I)
C          T(2,4)=dBY/dP4=X(I)*Y(I)*COS(3*K*Z(I)
C          T(3,1)=dBZ/dP1=(COS(K*Z)*Y)/K
C          T(3,2)=dBZ/dP2=(COS(3*K*Z)*Y)/(3*K)
C          T(3,3)=dBZ/dP3=(-SIN(K*Z)*X)/K
C          T(3,4)=dBZ/dP4=(-SIN(3*K*Z)*X)/(3*K)
C

      IMPLICIT NONE

        INTEGER NARG,NFUN,NPAR,IPOI,NDIMPOI

      DOUBLE PRECISION FUNDATA,T
      DIMENSION FUNDATA(NARG+NFUN,NDIMPOI),T(NFUN,NPAR)

c      f(x,y)=p1+p2*x+p3*y+p4*x**2+p5*y**2+p6*x*y

        t(1,1)=1.0d0
        t(1,2)=fundata(1,ipoi)
        t(1,3)=fundata(2,ipoi)
        t(1,4)=fundata(1,ipoi)**2
        t(1,5)=fundata(2,ipoi)**2
        t(1,6)=fundata(1,ipoi)*fundata(2,ipoi)

      RETURN
      END
