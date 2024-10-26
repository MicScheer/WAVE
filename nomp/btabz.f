*CMZ :  4.00/15 28/03/2022  12.42.45  by  Michael Scheer
*CMZ :  4.00/11 17/05/2021  11.34.33  by  Michael Scheer
*CMZ :  3.01/04 21/03/2014  14.02.03  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.54/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.52/09 29/10/2004  12.30.12  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.25.04  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.07.21  by  Michael Scheer
*CMZ : 00.01/08 31/05/95  13.36.58  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.57.38  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.48.44  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.59  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BTABZ(XIN,Y,Z,BX,BY,BZ,AX,AY,AZ,XS,XE)
*KEEP,gplhint.
*KEND.

C     BTAB-Version for horizontal field Bz

      use fbtabzymod

      IMPLICIT NONE

      CHARACTER(60) BTABCOM

      INTEGER :: ISYM,ICAL,I,NPOINT,ILOW,IHIGH,IWARN,ieof=0

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
      DOUBLE PRECISION AL(3),AH(3),AX,AY,AZ
      DOUBLE PRECISION XS,XE,ZDUM

      COMMON/BTABCz/XA,BYA,Y2A


      DATA ILOW/0/
      DATA IHIGH/0/
      DATA ISYM/0/,ICAL/0/

      IF (ICAL.NE.1) THEN

        ZDUM=Z

        OPEN (UNIT=LUNTBz,FILE = FILETBz,STATUS = 'OLD',FORM = 'FORMATTED')

        if (irbtabzy.gt.0.or.irbtabxyz.gt.0) then
          READ(LUNTBZ,'(1A60)') BTABCOM
          READ(LUNTBZ,*) XSCALE,BYSCALE
          READ(LUNTBZ,*) NPOINT
        else
          btabcom=trim(filetbz)
          xscale=1.0d0
          if (irbtabzy.eq.-2.or.irbtabxyz.eq.-2) xscale=0.001d0
          byscale=1.0d0
          npoint=0
          do while (ieof.eq.0)
            call util_skip_comment_end(luntbz,ieof)
            read(luntbz,*)x,by
            if (ieof.ne.0) then
              rewind(luntbz)
              exit
            endif
            npoint=npoint+1
          enddo
        endif

        IF (NPOINT.GT.NBTABP.OR.-2*NPOINT-1.GT.NBTABP) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN BTABZ ***'
          WRITE(LUNGFO,*)'DIMENSION EXCEEDED, INCREASE NBTABP***'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BTABZ ***'
          WRITE(6,*)'DIMENSION EXCEEDED, INCREASE NBTABP***'
          WRITE(6,*)
          STOP
        ENDIF

        IF (ABS(NPOINT).LT.2) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN BTABZ ***'
          WRITE(LUNGFO,*)
     &      'LESS THAN TWO POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BTABZ ***'
          WRITE(6,*)'LESS THAN TWO POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(6,*)
        ENDIF

        IF (NPOINT.GT.0.AND.NPOINT.LT.3) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN BTABZ ***'
          WRITE(LUNGFO,*)
     &      'LESS THAN THREE POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BTABZ ***'
          WRITE(6,*)'LESS THAN THREE POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(6,*)
        ENDIF

        IF (NPOINT.LT.0) THEN
          ISYM=1
          NPOINT=-NPOINT
        ENDIF

        DO I=1,NPOINT
          READ(LUNTBZ,*)XA(I),BYA(I)
          XA(I)=XA(I)*XSCALE-XSHBTAB
          BYA(I)=BYA(I)*BYSCALE
        END DO

        CLOSE(LUNTBZ)

        call util_sort_func(npoint,xa,bya)

C--- SYMMETRISIEREN

        IF (ISYM.EQ.1) THEN

          IF(XA(1).LT.0.) THEN
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)'*** ERROR IN BTABZ ***'
            WRITE(LUNGFO,*)'FIRST X-VALUE ON DATA FILE LOWER ZERO'
            WRITE(LUNGFO,*)'SYMMETRY OPTION REQUIRES POSITIV X-VALUES'
            WRITE(LUNGFO,*)
            WRITE(LUNGFO,*)
            WRITE(6,*)
            WRITE(6,*)'*** ERROR IN BTABZ ***'
            WRITE(6,*)'FIRST X-VALUE ON DATA FILE LOWER ZERO'
            WRITE(6,*)'SYMMETRY OPTION REQUIRES POSITIV X-VALUES'
            WRITE(6,*)
            STOP
          ENDIF !XA

          IF(XA(1).NE.0) THEN

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
        WRITE (LUNGFO,*)'     SR BTABZ: Magnetic field data read from file'
        WRITE (LUNGFO,*)'     ',FILETBz
        WRITE (LUNGFO,*)
        WRITE (LUNGFO,*)'     BTABZ comment:',BTABCOM
        WRITE (LUNGFO,*)
        WRITE (LUNGFO,*)'     Length of device:     ',TOTLEN
        WRITE (LUNGFO,*)'     Half length of device:',TOTLEN2
        WRITE (LUNGFO,*)

C060793  CALL SPLINETB(XA,BYA,NPOINT,2.D30,2.D30,Y2A)
        CALL SPLINETB(XA,BYA,NPOINT,0.D0,0.D0,Y2A)

        IF (XS.EQ.9999.)  XSTART=XA(1)
        IF (XE.EQ.9999.)  XSTOP=XA(NPOINT)

        ICAL=1

        nfbtabcz=npoint

      ENDIF !(ICAL)

      IF(Y.NE.0.) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN BTABZ ***'
        WRITE(LUNGFO,*)'Y-COORDINATE OF ELECTRON NOT ZERO'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN BTABZ ***'
        WRITE(6,*)'Y-COORDINATE OF ELECTRON NOT ZERO'
        WRITE(6,*)
        STOP
      ENDIF

      IF (XIN.EQ.9999.) THEN
        X=XSTART
      ELSE
        X=XIN
      ENDIF

      IF(X.LT.XA(1)-1./MYINUM.OR.X.GT.XA(NPOINT)+1./MYINUM) THEN
        IF (IWARN.EQ.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING IN BTABZ ***'
          WRITE(LUNGFO,*)'X MORE THAN ONE STEP OUT OF TABLE'
          WRITE(LUNGFO,*)'X:',X
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING IN BTABZ ***'
          WRITE(6,*)'X MORE THAN ONE STEP OUT OF TABLE'
          WRITE(6,*)'X:',X
          WRITE(6,*)
          IWARN=1
        ENDIF
        BX=0.0D0
        BY=0.0D0
        BZ=0.0D0
        AX=0.0D0
        AY=0.0D0
        AZ=0.0D0
        RETURN
      ENDIF

      IF(X.LT.XA(1)) THEN
        IF (ILOW.NE.1) THEN
          CALL PARABEL(XA(1),BYA(1),XA(2),BYA(2),XA(3),BYA(3),AL)
          ILOW=1
        ENDIF   ! ILOW

        Bx=0.
        By=0.
        Bz=AL(1)+AL(2)*X+AL(3)*X*X

        AX=-0.5*Bz*y
        Ay= 0.5*Bz*x
        Az= 0.0

        RETURN
      ENDIF !X.LT.XA(1)

      IF(X.GT.XA(NPOINT)) THEN

        IF (IHIGH.NE.1) THEN

          CALL PARABEL(XA(NPOINT-2),BYA(NPOINT-2)
     &      ,XA(NPOINT-1),BYA(NPOINT-1)
     &      ,XA(NPOINT)  ,BYA(NPOINT)  ,AH)

          IHIGH=1
        ENDIF   ! IHIGH

        By=0.
        Bx=0.
        Bz=AH(1)+AH(2)*X+AH(3)*X*X

        AX=-0.5*Bz*y
        Ay= 0.5*Bz*x
        Az= 0.0

        RETURN

      ENDIF !X.GT.XA(NPOINT)

      CALL SPLINTBz(NPOINT,X,Bz)

      BX=0.
      By=0.

      AX=-0.5*Bz*y
      Ay= 0.5*Bz*x
      Az= 0.0

      RETURN
      END
