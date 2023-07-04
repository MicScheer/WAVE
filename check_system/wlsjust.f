*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.28.35  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.04.09  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.20  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WLSJUST(XI,YI)
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

C UTILITY ROUTINE THAT ESTIMATES A NEW VALUE FOR X TO
C TO BRING Y=F(X) TO ZERO. THE PREVIOUS AND CURRENT VALUES
C WRITTEN TO THE FILE FILEJ, I.E. THE FILE IS UPDATED BY THE
C CALLS TO THE ROUTINE

      IMPLICIT NONE

      INTEGER ICAL
      DOUBLE PRECISION XI,YI,X0,Y0,X1,Y1,X,Y,AAA,BBB

      DATA ICAL/0/

      DATA X0,Y0,X1,Y1,X/5*0.0/

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEND.

      IF (IJUST.EQ.1) THEN

          IF (ICAL.NE.0) THEN

         X0=XI
         Y0=YI

         OPEN(UNIT=LUNJ,FILE=FILEJ,STATUS='NEW',FORM='FORMATTED')
         WRITE(LUNJ,*)X0,Y0
         CLOSE(LUNJ)

          ELSE

         ICAL=1

          ENDIF

      ELSE IF (IJUST.EQ.2) THEN

          IF (ICAL.NE.0) THEN

         Y1=YI

         WRITE(LUNJ,*)X0,Y0
         WRITE(LUNJ,*)X1,Y1

         CLOSE(LUNJ)

          ELSE

         OPEN(UNIT=LUNJ,FILE=FILEJ,STATUS='OLD',FORM='FORMATTED')
         READ(LUNJ,*)X0,Y0
         REWIND (LUNJ)

         X1=X0*1.10
         XI=X1

         ICAL=1

          ENDIF

      ELSE IF (IJUST.EQ.3) THEN

          IF (ICAL.NE.0) THEN

         Y=YI

         WRITE(LUNJ,*)X1,Y1
         WRITE(LUNJ,*)X,Y

         CLOSE(LUNJ)

          ELSE

         OPEN(UNIT=LUNJ,FILE=FILEJ,STATUS='OLD',FORM='FORMATTED')
         READ(LUNJ,*)X0,Y0
         READ(LUNJ,*)X1,Y1
         REWIND (LUNJ)

         AAA=(Y1-Y0)/(X1-X0)
         BBB=Y1-AAA*X1
         X=-BBB/AAA

         XI=X

         ICAL=1

          ENDIF

      ELSE IF (IJUST.EQ.4) THEN

          IF (ICAL.EQ.0) THEN

         OPEN(UNIT=LUNJ,FILE=FILEJ,STATUS='OLD',FORM='FORMATTED')
         READ(LUNJ,*)X0,Y0
         READ(LUNJ,*)X1,Y1

         XI=X1

         CLOSE(LUNJ)

         ICAL =1

          ENDIF

      ELSE

          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** ERROR IN WLSJUST ***'
          WRITE(LUNGFO,*)'IJUST WRONG. CHECK NAMELIST CONTRL'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN WLSJUST ***'
          WRITE(6,*)'IJUST WRONG. CHECK NAMELIST CONTRL'
          WRITE(6,*)
          STOP

      ENDIF


      RETURN
      END
