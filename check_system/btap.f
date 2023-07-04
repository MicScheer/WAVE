*CMZ :  4.00/07 10/07/2020  09.20.19  by  Michael Scheer
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
      SUBROUTINE BTAP(X,Y,Z,BX,BY,BZ,AX,AY,AZ)
*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*KEND.

      IMPLICIT NONE

      INTEGER ISYM,ICAL,I,NPOINT,IWARNY,IWARN,IMONO,ifaili,ieof,ifaile,lunt

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEND.

      DOUBLE PRECISION XA(NBTABP),BYA(NBTABP),Y2A(NBTABP)
      DOUBLE PRECISION XSCALE,BYSCALE,X,Y,Z,BX,BY,BZ,TOTLEN,TOTLEN2
      DOUBLE PRECISION AX,AY,AZ,apl,aph,x0l,x0h,xx,bb

      CHARACTER(60) BTAPCOM

      DATA ISYM/0/,ICAL/0/
      DATA IWARN/0/

      save

      IF (ICAL.NE.1) THEN

        ical=1
        IWARNY=0

        OPEN (newUNIT=lunt,FILE = FILEFTV,STATUS = 'OLD',FORM = 'FORMATTED')

        call util_skip_comment_end(lunt,ieof)
        READ(lunt,'(1A60)') BTAPCOM
        call util_skip_comment_end(lunt,ieof)
        READ(lunt,*) XSCALE,BYSCALE

        npoint=0

1       continue
        call util_skip_comment_end(lunt,ieof)
        if (ieof.ne.0) goto 9
        READ(lunt,*,end=9) xx,bb
        npoint=npoint+1
        IF (NPOINT.GT.NBTABP) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN BTAP ***'
          WRITE(LUNGFO,*)'DIMENSION EXCEEDED, INCREASE NBTABP***'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BTAP ***'
          WRITE(6,*)'DIMENSION EXCEEDED, INCREASE NBTABP***'
          WRITE(6,*)
          STOP
        ENDIF

        xa(npoint)=xx
        bya(npoint)=bb
        goto 1
9       close(lunt)

        IF (ABS(NPOINT).LT.2) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN BTAP ***'
          WRITE(LUNGFO,*)
     &      'LESS THAN TWO POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BTAP ***'
          WRITE(6,*)'LESS THAN TWO POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(6,*)
          STOP
        ENDIF

        IF (NPOINT.GT.0.AND.NPOINT.LT.3) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN BTAP ***'
          WRITE(LUNGFO,*)
     &      'LESS THAN THREE POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN BTAP ***'
          WRITE(6,*)'LESS THAN THREE POINTS ON DATA FILE OF MAGNETIC FIELD'
          WRITE(6,*)
        ENDIF

        call util_parabola_to_zero(xa(1:2),bya(1:2),apl,x0l,ifaili)
        call util_parabola_to_zero(xa(npoint-1:npoint),bya(npoint-1:npoint),
     &    aph,x0h,ifaile)

        call util_sort_func(npoint,xa,bya)

        CALL UTIL_CHECK_MONOTON(NPOINT,XA,IMONO)
        IF (ABS(IMONO).NE.2) THEN
          PRINT *,'*** ERROR IN BTAP: Field data not monoton'
          WRITE(LUNGFO,*)'*** ERROR IN BTAP: FIELD DATA NOT MONOTON'
          STOP '*** Program WAVE aborted'
        ENDIF

        WRITE (LUNGFO,*)
        WRITE (LUNGFO,*)'     Subroutine BTAP: Magnetic field data read from file'
        WRITE (LUNGFO,*)'     ',fileftv
        WRITE (LUNGFO,*)
        WRITE (LUNGFO,*)'     BTAP comment:',BTAPCOM
        WRITE (LUNGFO,*)

        CALL SPLINETB(XA,BYA,NPOINT,0.D0,0.D0,Y2A)

      ENDIF !(ICAL)

      IF(Y.NE.0..AND.IWARNY.EQ.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING SR BTAP ***'
        WRITE(LUNGFO,*)'Y-COORDINATE OF ELECTRON NOT ZERO'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING SR BTAP ***'
        WRITE(6,*)'Y-COORDINATE OF ELECTRON NOT ZERO'
        WRITE(6,*)
        IWARNY=1
      ENDIF

      IF(IWARN.EQ.0.AND.
     &    (X.LT.XA(1)-1./MYINUM.OR.X.GT.XA(NPOINT)+1./MYINUM)) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN BTAP ***'
        WRITE(LUNGFO,*)'X MORE THAN ONE STEP OUT OF TABLE'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'X, XMIN, XMAX:'
        WRITE(LUNGFO,*)X,XA(1),XA(NPOINT)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN BTAP ***'
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
            print*,"*** WARNING IN BTAP: Failed to find extrapolation parabola at intrance"
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
            print*,"*** WARNING IN BTAP: Failed to find extrapolation parabola at exit"
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
      CALL SPLINTap(NPOINT,xa,bya,y2a,X,BY)

      BX=0.
      BZ=0.

      AX=0.5*BY*Z
      AY=0.0
      AZ=-0.5*BY*X

      RETURN
      END
