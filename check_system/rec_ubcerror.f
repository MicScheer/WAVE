*CMZ :  4.00/15 27/04/2022  09.00.10  by  Michael Scheer
*CMZ :  3.02/03 03/11/2014  12.00.05  by  Michael Scheer
*CMZ :  2.42/00 05/09/2002  18.59.21  by  Michael Scheer
*-- Author :    Michael Scheer   04/09/2002
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

      SUBROUTINE REC_UBCERROR(N,BC,NURANMOD,BCRAN,BCRANSIG)

      INTEGER N,NURANMOD,I,NPOLTOT
      DOUBLE PRECISION BC(N),BCRAN,BERR,BDIFF,BCORRI,BCORRE,UBCSUM2,BCSIGN
      REAL BCRANSIG,xran(1),rr(2)

      IF (N.LE.6) THEN
          RETURN
      ELSE IF (MOD(N,3).NE.0) THEN
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

      DO I=6+1,N-6,3
1         CALL RNORML(xran(1),1,rr)
          IF (ABS(xran(1)).GT.BCRANSIG) GOTO 1
          BERR=(1.+BCRAN*xran(1))
          BC(I)=BC(I)*BERR
          BC(I+1)=BC(I+1)*BERR
          BC(I+2)=BC(I+2)*BERR
      ENDDO !I

C--- COMPENSATE THE KICK

      BDIFF=0.D0
      DO I=6+1+1,N-6,6
            BDIFF=BDIFF-BC(I)+BC(I+3)
      ENDDO !I

      BCORRI=1.D0-BDIFF/BC(8)
      BCORRE=1.D0+BDIFF/BC(8)

      DO I=1,6
            BC(I)=BC(I)*BCORRI
            BC(N+1-I)=BC(N+1-I)*BCORRE
      ENDDO

      IF (NURANMOD.EQ.4) THEN

C--- COMPENSATE OFFSET WITH BY ENDPOLE

           UBCSUM2=0.D0
           NPOLTOT=N/3

           DO I=1,NPOLTOT
                    IF (I.EQ.1) THEN
                        BCSIGN=1.D0
                    ELSEIF(I.EQ.2) THEN
                        BCSIGN=-1.D0
                    ELSEIF(I.EQ.NPOLTOT-1) THEN
                        BCSIGN=-1.D0
                    ELSEIF(I.EQ.NPOLTOT) THEN
                        BCSIGN=1.D0
                    ELSE
                        BCSIGN=(-1.D0)**I
                    ENDIF
                    UBCSUM2=UBCSUM2+(NPOLTOT+1-I)*BCSIGN*BC(I*3-1)
                ENDDO !I

                BCORRI=1.D0-UBCSUM2/(NPOLTOT-2)/BC(8)*2.D0

                DO I=1,6
                    BC(I)=BC(I)*BCORRI
                ENDDO

                BCORRE=BCORRI

                DO I=1,6
                    BC(N+1-I)=BC(N+1-I)*BCORRE
                ENDDO

            ENDIF   !NURANMOD.EQ.4

      RETURN
      END
