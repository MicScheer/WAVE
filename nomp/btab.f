*CMZ :  4.00/16 22/07/2022  10.31.35  by  Michael Scheer
*CMZ :  4.00/15 28/03/2022  12.44.37  by  Michael Scheer
*CMZ :  4.00/11 17/05/2021  11.34.14  by  Michael Scheer
*CMZ :  4.00/07 09/07/2020  16.52.53  by  Michael Scheer
*CMZ :  4.00/03 07/05/2019  14.24.38  by  Michael Scheer
*CMZ :  3.05/23 23/11/2018  18.01.47  by  Michael Scheer
*CMZ :  3.04/00 19/01/2018  13.23.56  by  Michael Scheer
*CMZ :  3.03/02 16/02/2017  13.03.18  by  Michael Scheer
*CMZ :  3.01/04 09/05/2014  14.31.02  by  Michael Scheer
*CMZ :  3.01/03 20/03/2014  13.03.02  by  Michael Scheer
*CMZ :  2.63/02 24/01/2008  15.24.28  by  Michael Scheer
*CMZ :  2.52/11 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.52/09 29/10/2004  12.30.03  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.13/10 25/03/2000  14.36.30  by  Michael Scheer
*CMZ :  2.13/11 22/03/2000  13.00.11  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.25.04  by  Michael Scheer
*CMZ :  2.13/00 03/12/99  16.05.55  by  Michael Scheer
*CMZ :  2.11/00 10/05/99  17.40.30  by  Michael Scheer
*CMZ :  1.04/00 25/11/98  14.08.28  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.07.21  by  Michael Scheer
*CMZ : 00.02/05 03/03/97  12.25.12  by  Michael Scheer
*CMZ : 00.01/12 15/10/96  12.15.47  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.56.48  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.48.39  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.35  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BTAB(XIN,Y,Z,BX,BY,BZ,AX,AY,AZ)
*KEEP,gplhint.
*KEND.

C     BTAB LIEST B-FELD TABELLE UND BERECHNET DURCH SPLINE-INTERPOLATION B-FELD

      use fbtabzymod

      IMPLICIT NONE

      CHARACTER(60) BTABCOM

      INTEGER ISYM,ICAL,I,NPOINT,IWARNY,IWARN,IMONO,ifaili,ieof,ifaile

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEND.

      DOUBLE PRECISION XA(NBTABP),BYA(NBTABP),Y2A(NBTABP)
      DOUBLE PRECISION XSCALE,BYSCALE,X,Y,Z,BX,BY,BZ,XIN,TOTLEN,TOTLEN2
      DOUBLE PRECISION AX,AY,AZ,apl,aph,x0l,x0h

      COMMON/BTABC/XA,BYA,Y2A

      DATA ISYM/0/,ICAL/0/
      DATA IWARN/0/

      IF (ICAL.NE.1) THEN

        IWARNY=0

        OPEN (UNIT=LUNTB,FILE = FILETB,STATUS = 'OLD',FORM = 'FORMATTED')

        if (irbtab.gt.0.or.irbtabzy.gt.0.or.irbtabxyz.gt.0) then
          call util_skip_comment_end(luntb,ieof)
          READ(LUNTB,'(1A60)') BTABCOM
          call util_skip_comment_end(luntb,ieof)
          READ(LUNTB,*) XSCALE,BYSCALE
          call util_skip_comment_end(luntb,ieof)
          READ(LUNTB,*) NPOINT
          call util_skip_comment_end(luntb,ieof)
        else
          npoint=0
          btabcom=filetb(1:60)
1         continue
          call util_skip_comment_end(luntb,ieof)
          if (ieof.ne.0) goto 9
          READ(LUNTB,*) x,by
          npoint=npoint+1
          goto 1
9         rewind(luntb)
          xscale=1.0d0
          if (irbtab.eq.-2) xscale=0.001d0
          byscale=1.0d0
        endif

        IF (NPOINT.GT.NBTABP.OR.-2*NPOINT-1.GT.NBTABP) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN BTAB ***'
          WRITE(LUNGFO,*)'DIMENSION EXCEEDED, INCREASE NBTABP***'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BTAB ***'
          WRITE(6,*)'DIMENSION EXCEEDED, INCREASE NBTABP***'
          WRITE(6,*)
          STOP
        ENDIF

        IF (ABS(NPOINT).LT.2) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN BTAB ***'
          WRITE(LUNGFO,*)
     &      'LESS THAN TWO POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BTAB ***'
          WRITE(6,*)'LESS THAN TWO POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(6,*)
          STOP
        ENDIF

        IF (NPOINT.GT.0.AND.NPOINT.LT.3) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN BTAB ***'
          WRITE(LUNGFO,*)
     &      'LESS THAN THREE POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BTAB ***'
          WRITE(6,*)'LESS THAN THREE POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(6,*)
        ENDIF

        IF (NPOINT.LT.0) THEN
          ISYM=1
          NPOINT=-NPOINT
        ENDIF

        DO I=1,NPOINT
          call util_skip_comment_end(luntb,ieof)
          READ(LUNTB,*)XA(I),BYA(I)
          XA(I)=XA(I)*XSCALE-XSHBTAB
          BYA(I)=BYA(I)*BYSCALE
c          write(6,*)i,xa(i),bya(i)
        END DO

        call util_parabola_to_zero(xa(1:2),bya(1:2),apl,x0l,ifaili)
        call util_parabola_to_zero(xa(npoint-1:npoint),bya(npoint-1:npoint),
     &    aph,x0h,ifaile)

        CLOSE(LUNTB)

        call util_sort_func(npoint,xa,bya)

        CALL UTIL_CHECK_MONOTON(NPOINT,XA,IMONO)
        IF (ABS(IMONO).NE.2) THEN
          PRINT *,'*** ERROR IN BTAB: Field data not monoton'
          WRITE(LUNGFO,*)'*** ERROR IN BTAB: FIELD DATA NOT MONOTON'
          STOP '*** Program WAVE aborted'
        ENDIF

C--- SYMMETRISIEREN

        IF (ISYM.EQ.1) THEN

          IF(XA(1).LT.0.) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN BTAB ***'
            WRITE(LUNGFO,*)'FIRST X-VALUE ON DATA FILE LOWER ZERO'
            WRITE(LUNGFO,*)'SYMMETRY OPTION REQUIRES POSITIV X-VALUES'
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN BTAB ***'
            WRITE(6,*)'FIRST X-VALUE ON DATA FILE LOWER ZERO'
            WRITE(6,*)'SYMMETRY OPTION REQUIRES POSITIV X-VALUES'
            WRITE(6,*)
            STOP
          ENDIF !XA

          IF(XA(1).NE.0.) THEN

            DO I=1,NPOINT        !SHIFT
              XA(I+NPOINT)= XA(I)
              BYA(I+NPOINT)=BYA(I)
            END DO

            DO I=1,NPOINT
              XA(I)=-XA(2*NPOINT-I+1)
              BYA(I)=BYA(2*NPOINT-I+1)
            END DO

            NPOINT=2*NPOINT

          ELSE

            DO I=1,NPOINT        !SHIFT

              XA(2*NPOINT-I)= XA(NPOINT-I+1)
              BYA(2*NPOINT-I)=BYA(NPOINT-I+1)

            END DO

            DO I=1,NPOINT-1
              XA(I) =-XA(2*NPOINT-I)
              BYA(I)=BYA(2*NPOINT-I)
            END DO

            NPOINT=2*NPOINT-1

          ENDIF

        ENDIF

        TOTLEN=DABS(XA(NPOINT)-XA(1))
        TOTLEN2=TOTLEN/2.D0
        DEVLEN=TOTLEN
        DEVLEN2=TOTLEN2

        WRITE (LUNGFO,*)
        WRITE (LUNGFO,*)'     SR BTAB: Magnetic field data read from file'
        WRITE (LUNGFO,*)'     ',FILETB
        WRITE (LUNGFO,*)
        WRITE (LUNGFO,*)'     BTAB comment:',BTABCOM
        WRITE (LUNGFO,*)
        WRITE (LUNGFO,*)'     Length of device:     ',TOTLEN
        WRITE (LUNGFO,*)'     Half length of device:',TOTLEN2
        WRITE (LUNGFO,*)

C060793  CALL SPLINETB(XA,BYA,NPOINT,2.D30,2.D30,Y2A)
        CALL SPLINETB(XA,BYA,NPOINT,0.D0,0.D0,Y2A)

        IF (XSTART.EQ.9999.) XSTART=XA(1)
        IF ( XSTOP.EQ.9999.)  XSTOP=XA(NPOINT)

        nfbtabc=npoint
        ical=1
      ENDIF !(ICAL)

      IF(Y.NE.0..AND.IWARNY.EQ.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING SR BTAB ***'
        WRITE(LUNGFO,*)'Y-COORDINATE OF ELECTRON NOT ZERO'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING SR BTAB ***'
        WRITE(6,*)'Y-COORDINATE OF ELECTRON NOT ZERO'
        WRITE(6,*)
        IWARNY=1
      ENDIF

C22.3.93 X=DMOD(XIN,TOTLEN) !2.12.91
C22.3.93 IF (X.GT.TOTLEN2) THEN
C22.3.93     X=X-TOTLEN
C22.3.93 ELSEIF (X.LT.-TOTLEN2) THEN
C22.3.93     X=X+TOTLEN
C22.3.93 ENDIF

C22.3.93 IF(X.LT.XA(1).OR.X.GT.XA(NPOINT)) THEN
C22.3.93     STOP '*** S/R BTAB: X MORE THAN ONE STEP OUT OF TABLE ***'
C22.3.93 ENDIF


C22.3.93 --------------------------------------------------------

      IF (XIN.EQ.9999.) THEN
        X=XSTART
      ELSE
        X=XIN
      ENDIF

      IF(iperiodg.eq.0.and.IWARN.EQ.0.AND.
     &    (X.LT.XA(1)-1./MYINUM.OR.X.GT.XA(NPOINT)+1./MYINUM)) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN BTAB ***'
        WRITE(LUNGFO,*)'X MORE THAN ONE STEP OUT OF TABLE'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'X, XMIN, XMAX:'
        WRITE(LUNGFO,*)X,XA(1),XA(NPOINT)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN BTAB ***'
        WRITE(6,*)'X MORE THAN ONE STEP OUT OF TABLE'
        WRITE(6,*)
        WRITE(6,*)'X, XMIN, XMAX:'
        WRITE(6,*)X,XA(1),XA(NPOINT)
        print*,'Field and vector potential smoothed to zero'
        IWARN=1
      ENDIF

      IF(X.LT.XA(1)) THEN
        BX=0.0D0
        BZ=0.0D0
        if (x.lt.x0l.or.ifaili.ne.0) then
          if (ifaili.eq.-1) then
            print*,"*** WARNING IN BTAB: Failed to find extrapolation parabola at intrance"
            ifaili=-2
          endif
          BY=0.0D0
        else
          by=apl*(x-x0l)**2
        endif
        AX=0.5*BY*Z
        AY=0.0
        AZ=-0.5*BY*X
        RETURN
      ENDIF !X.LT.XA(1)

      IF(X.GT.XA(NPOINT)) THEN
        BX=0.0D0
        BZ=0.0D0
        if (x.gt.x0h.or.ifaile.ne.0) then
          if (ifaile.eq.-1) then
            print*,"*** WARNING IN BTAB: Failed to find extrapolation parabola at exit"
            ifaile=-2
          endif
          BY=0.0D0
        else
          by=aph*(x-x0h)**2
        endif
        AX=0.5*BY*Z
        AY=0.0
        AZ=-0.5*BY*X

        RETURN
      ENDIF !X.GT.XA(NPOINT)

C22.3.93 --------------------------------------------------------

c      write(6,*)ical,npoint,x,by
      CALL SPLINTB(NPOINT,X,BY)

      BX=0.
      BZ=0.

      AX=0.5*BY*Z
      AY=0.0
      AZ=-0.5*BY*X

      RETURN
      END
