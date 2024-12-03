*CMZ :  4.00/14 22/12/2021  18.07.25  by  Michael Scheer
*CMZ :  4.00/07 28/04/2020  21.30.22  by  Michael Scheer
*CMZ :  3.08/01 02/04/2019  15.33.15  by  Michael Scheer
*CMZ :  3.07/01 21/03/2019  15.31.19  by  Michael Scheer
*CMZ :  3.05/14 28/09/2018  13.00.20  by  Michael Scheer
*CMZ :  3.05/01 08/05/2018  16.09.00  by  Michael Scheer
*CMZ :  3.02/06 15/04/2015  11.35.50  by  Michael Scheer
*CMZ :  3.02/03 06/11/2014  14.24.54  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  14.04.56  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.55.57  by  Michael Scheer
*CMZ :  2.50/00 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.49/00 22/03/2004  14.04.26  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.45  by  Michael Scheer
*CMZ :  2.12/00 27/05/99  10.08.55  by  Michael Scheer
*CMZ :  2.10/01 24/02/99  10.20.40  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.18.17  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.54.05  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.44  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE SOUINT(ISOUR,IBUFF)

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
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,reargf90u.
      include 'reargf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEND.
      use ompmod
C--- EVALUATE INTEGRALES FOR A SINGLE SOURCE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,datetime.
      include 'datetime.cmn'
*KEND.

      INTEGER ISOUR,IBUFF,IOBSV,JX10,JDX10,IX10,ICAL,NNBUFF,ICYCLE

      INTEGER NTUPP
      PARAMETER (NTUPP=22)
      CHARACTER(4) CTUP(NTUPP)

      DATA ICAL/0/

      data ctup /'t','x','y','z','rx','ry','rz','rt','p','rea','ima','roi'
     &            ,'iob','ie','yob','zob','betn','dtom','emod','dmod'
     &            ,'spec','te'/

      if (iomp.eq.1.and.ipin.ne.0) then
        !print*,"******************* OMP ****************"
        call souint_omp(isour,ibuff)
        return
      endif

      IF (ICAL.EQ.0.AND.IBUFF.EQ.1)  THEN
        IX10=1
        JDX10=NBUFF*NOBSV/10
        JX10=JDX10
        NNBUFF=1
        IF (IWFILINT.EQ.-ISOUR)
     &    CALL hbookm(NIDSOURCE,'RADIATION INTEGRAL',NTUPP
     &    ,'//WAVE',nlpoi/jwfilint+2*jwfilint,CTUP)
      ENDIF

      IF (JDX10.LT.1) JDX10=1

      IF (ISOUR.EQ.1.AND.IBUFF.EQ.1.AND.NOBSV.GT.1) THEN
        WRITE(6,*)' '
        WRITE(6,*)
     &    '      counting from 1 to 10 for first source to show progress:'
        WRITE(6,*)' '
      ENDIF

C--- LOOP OVER ALL OBSERVATION POINTS

      DO IOBSV=1,NOBSV

C- CALCULATE FREQUENCE INDEPENDENT PARTS OF INTEGRANTES,STORE RESULT IN ARRAYS

        CALL REARG(ISOUR,IOBSV)  !REAL PARTS OF INTEGRANTS

C--- INTEGRATION FOR ALL FREQUENCES

        CALL ARGSUM(ISOUR,IOBSV,IBUFF)

        IF(IWFILINT.EQ.ISOUR.AND.IOBSV.EQ.1) CALL WFILINT

        NNBUFF=NNBUFF+1

        IF (ICAL.EQ.0.AND.NNBUFF.EQ.JX10.AND.IX10.LE.10.AND.NOBSV.GT.1) THEN
          JX10=JX10+JDX10
          CALL date_and_time(dtday,dttime,dtzone,idatetime)
          WRITE(6,*)' ',IX10,' ',dttime(1:2),':',dttime(3:4),':',dttime(5:6)
          IX10=IX10+1
        ENDIF

      ENDDO !LOOP OVER ALL OBSERVATION POINTS


      IF (IBUFF.EQ.NBUFF) THEN
        IF (IWFILINT.EQ.-ISOUR) THEN
          CALL MHROUT(NIDSOURCE,ICYCLE,' ')
          CALL hdeletm(NIDSOURCE)
        ENDIF
        ICAL=1
      ENDIF

      RETURN
      END
