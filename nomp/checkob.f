*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.17/00 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.44  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.12/03 22/07/99  10.49.26  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.40.33  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.49.00  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.51  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE CHECKOB
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

*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- CHECKS CONSISTENCE OF OBSERVATION POINTS AND COLLIMATOR

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,colli.
      include 'colli.cmn'
*KEND.

      INTEGER IOBSV

      DOUBLE PRECISION XOB,YOB,ZOB
      DOUBLE PRECISION UP1,DOWN1,RIGHT1,LEFT1,UP2,DOWN2,RIGHT2,LEFT2
      DOUBLE PRECISION YSLOPEU,YSLOPED,ZSLOPEL,ZSLOPER,DLEN,DIST

      UP1=CY1+HIG1/2.D0
      DOWN1=CY1-HIG1/2.D0
      RIGHT1=CZ1+WID1/2.D0
      LEFT1=CZ1-WID1/2.D0

      UP2=CY2+HIG2/2.D0
      DOWN2=CY2-HIG2/2.D0
      RIGHT2=CZ2+WID2/2.D0
      LEFT2=CZ2-WID2/2.D0

      DLEN=CX2-CX1

      IF (DLEN.EQ.0.0) THEN
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** WARNING SR CHECKOB ***'
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'COLLIMATOR HAS ZERO LENGTH, OBSERVER POSITIONS ARE NOT CHECKED, BE CAREFUL !'
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(6,*)'*** WARNING SR CHECKOB ***'
              WRITE(6,*)
              WRITE(6,*)'COLLIMATOR HAS ZERO LENGTH, OBSERVER POSITIONS ARE NOT CHECKED, BE CAREFUL !'
              WRITE(6,*)
              WRITE(6,*)
            RETURN
      ENDIF

      YSLOPEU=(UP2-DOWN1)/DLEN
      YSLOPED=(DOWN2-UP1)/DLEN
      ZSLOPEL=(LEFT2-RIGHT1)/DLEN
      ZSLOPER=(RIGHT2-LEFT1)/DLEN

      DO IOBSV=1,NOBSV
          XOB=OBSV(1,IOBSV)
          YOB=OBSV(2,IOBSV)
          ZOB=OBSV(3,IOBSV)
          DIST=XOB-CX2

          IF(
     &         XOB.LT.CX2
     &         .OR.
     &         YOB.GT.UP2+YSLOPEU*DIST
     &         .OR.
     &         YOB.LT.DOWN2+YSLOPED*DIST
     &         .OR.
     &         ZOB.GT.RIGHT2+ZSLOPER*DIST
     &         .OR.
     &         ZOB.LT.LEFT2+ZSLOPEL*DIST
     &          ) THEN
              WRITE(6,*)
              WRITE(6,*)'*** WARNING IN SR CHECKOB ***'
              WRITE(6,*)
     &        'OBSERVATION POINTS AND COLLIMATOR INCONSISTENT'
              WRITE(6,*)'CHECK RESULTS CAREFULLY'
              WRITE(6,*)
              WRITE(6,*)'OBSERVATION POINT (X,Y,Z):'
              WRITE(6,*)XOB,YOB,ZOB
              WRITE(6,*)
              WRITE(6,*)'CHECK ALSO DEFINITION OF COLLIMATOR'
              WRITE(6,*)
     &'IN NAMELIST COLLIN; MAYBE DEFAULTS NOT USEFUL'
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'*** WARNING IN SR CHECKOB ***'
              WRITE(LUNGFO,*)
     &        'OBSERVATION POINTS AND COLLIMATOR INCONSISTENT'
              WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'OBSERVATION POINT (X,Y,Z):'
              WRITE(LUNGFO,*)XOB,YOB,ZOB
              WRITE(LUNGFO,*)'(REMAINING POINTS NOT CHECKED)'
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)
              WRITE(LUNGFO,*)'CHECK ALSO DEFINITION OF COLLIMATOR'
              WRITE(LUNGFO,*)
     &'IN NAMELIST COLLIN; MAYBE DEFAULTS NOT USEFUL'
            RETURN
         ENDIF
      ENDDO

      RETURN
      END
