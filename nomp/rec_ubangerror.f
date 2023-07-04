*CMZ :  4.00/15 27/04/2022  09.02.34  by  Michael Scheer
*CMZ :  3.03/02 20/01/2016  12.04.25  by  Michael Scheer
*CMZ :  3.02/03 03/11/2014  12.08.30  by  Michael Scheer
*CMZ :  2.44/00 30/10/2002  13.51.40  by  Michael Scheer
*CMZ :  2.42/04 28/10/2002  11.43.42  by  Michael Scheer
*CMZ :  2.42/00 09/09/2002  18.30.24  by  Michael Scheer
*-- Author :    Michael Scheer   04/09/2002

      SUBROUTINE REC_UBANGERROR(N,BC,NURANMOD,UBANGERR,ubansig,USIGOFFY)
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

      INTEGER N,NURANMOD,I,NPOLTOT,J
      DOUBLE PRECISION BC(N),UBANGERR,BDIFF,BCORRI,BCORRE,UBCSUM2
     &                  ,USIGOFFY,GRARAD,ANGERR,ubansig
      REAL xran(1),XWALK,XCORR,rr(2)

      PARAMETER (GRARAD=0.0174532925199D0)

      IF (N.LE.6) THEN
          RETURN
      ELSEIF (MOD(N,3).NE.0) THEN
          WRITE(6,*)
          WRITE(6,*)'*** ERROR IN REC_UBCERROR: ZERO MAGNET STRENGTH  ***'
          STOP
      ENDIF

      IF (BC(8).EQ.0.D0) THEN
          WRITE(6,*)
          WRITE(6,*)'*** WARNING IN REC_UBCERROR: ZERO MAGNET STRENGTH  ***'
          WRITE(6,*)'*** LEAVING...NO ERRORS GENERATED'
          WRITE(6,*)
          RETURN
      ENDIF

      DO I=1,6
          BC(I)=0.D0
      ENDDO !I

      DO I=1,6
          BC(N+1-I)=0.D0
      ENDDO !I

        DO I=6+1,N-6,3
123       CALL RNORML(XRAN,1,rr)
          if (abs(xran(1)).gt.ubansig) goto 123
           XWALK=XWALK+(xran(1)-XCORR)
           IF (ABS(XWALK).GT.USIGOFFY) THEN
             XCORR=XWALK/2.D0
           ENDIF
           J=(I+2)*2-1
           ANGERR=(xran(1)-XCORR)*UBANGERR*GRARAD
           BC(I)=BC(I)*ANGERR
           BC(I+1)=BC(I+1)*ANGERR
           BC(I+2)=BC(I+2)*ANGERR
        ENDDO

      BDIFF=0.D0
      DO I=6+1,N-6,3
                BDIFF=BDIFF+BC(I)
      ENDDO !I

      DO I=1,3
                BC(I)=+BDIFF/4.D0
      ENDDO

      DO I=4,6
                BC(I)=-BDIFF*3.D0/4.D0
      ENDDO

      DO I=1,3
                BC(N+1-I)=+BDIFF/4.D0
      ENDDO

      DO I=4,6
                BC(N+1-I)=-BDIFF*3.D0/4.D0
      ENDDO

C--- COMPENSATE VERTICAL OFFSET WITH ENDPOLE

      IF (NURANMOD.EQ.4) THEN

          UBCSUM2=0.D0
        NPOLTOT=N/3
        DO I=1,NPOLTOT
          J=-1+I*3
          UBCSUM2=UBCSUM2+(NPOLTOT+1-I)*BC(J)
        ENDDO !I

        IF (UBCSUM2.NE.0.D0) THEN
                    BCORRI=1.D0-UBCSUM2/(NPOLTOT-1)/BC(2)
                    BCORRE=1.D0+UBCSUM2/(NPOLTOT-1)/BC(N+1-2)
        ELSE
                    BCORRI=1.D0
                    BCORRE=1.D0
        ENDIF

        DO I=1,3
                    BC(I)=BC(I)*BCORRI
        ENDDO


        DO I=1,3
                    BC(N+1-I)=BC(N+1-I)*BCORRE
        ENDDO

      ENDIF   !NURANMOD.EQ.4

      RETURN
      END
