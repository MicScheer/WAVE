*CMZ : 00.00/07 22/03/2010  15.28.00  by  Michael Scheer
*CMZ : 00.00/01 20/06/95  10.09.17  by  Michael Scheer
*-- Author :    Michael Scheer   26/01/95

      SUBROUTINE util_min_parabel(NDIM,X,F,XMN,FMN,WSX,WSF,JFAIL)

C--- TO FIND MINIMUM OF ARRAY FUNCTION F(X)

C     INPUT : F(NDIM)   ARRAY OF FUNCTION
C             X(NDIM)   ARRAY OF ARGUMENTS
C             WSX(NDIM) WORKINGSPACE
C             WSF(NDIM) WORKINGSPACE

C     OUTPUT:  XMN ARGUMENT WHERE FUNCTION REACHES EXTREMUM
C              FMN MINIMUM OF FUNCTION
C              JFAIL FLAG: =0, IF OK, =1 ELSE

      IMPLICIT NONE

      INTEGER NDIM,I,IFAIL,JFAIL

      REAL*8 X(NDIM),F(NDIM),XMN,FMN,WSX(NDIM),WSF(NDIM)
      REAL*8 XDUM(3),FDUM(3),A(3),FP(3),FMIN,XMIN

      FMIN=1.0D30
      DO I=1,NDIM
          WSX(I)=X(I)
          WSF(I)=F(I)
          IF (WSF(I).LT.FMIN) THEN
              FMIN=WSF(I)
              XMIN=WSX(I)
          ENDIF
      ENDDO

      CALL UTIL_SORT_FUNC(NDIM,WSF,WSX)

      DO I=1,3
         XDUM(I)=WSX(I)
         FDUM(I)=WSF(I)
      ENDDO

      CALL UTIL_PARABEL(XDUM,FDUM,A,FP,XMN,FMN,IFAIL)

      IF (IFAIL.NE.0.OR.FMN.LT.FMIN) THEN
          JFAIL=1
c          WRITE(6,*)
c          WRITE(6,*)'*** WARNING IN UTIL_MIN_PARABEL: SEARCH FAILED ***'
c          WRITE(6,*)'*** MINIMUM OF ARRAY TAKEN ***'
c          WRITE(6,*)
          FMN=FMIN
          XMN=XMIN
          RETURN
      ENDIF

      JFAIL=0

      RETURN
      END
