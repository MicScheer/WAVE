*CMZ :  3.02/00 28/08/2014  08.42.42  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.69/00 26/10/2012  11.36.41  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  11.44.11  by  Michael Scheer
*CMZ :  2.62/03 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.62/02 16/07/2007  07.15.46  by  Michael Scheer
*CMZ :  2.34/09 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.45  by  Michael Scheer
*CMZ :  1.04/03 11/12/98  16.34.58  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.47.16  by  Michael Scheer
*CMZ : 00.02/04 24/02/97  12.37.49  by  Michael Scheer
*CMZ : 00.02/00 19/11/96  14.57.13  by  Michael Scheer
*CMZ : 00.01/08 22/06/95  17.29.50  by  Michael Scheer
*CMZ : 00.01/06 01/02/95  16.35.43  by  Michael Scheer
*CMZ : 00.01/04 29/11/94  10.17.51  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.07.25  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.52.56  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.26  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE PINCENIN
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
*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.

C--- INITIALIZE PINCEN

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,trackf90.
      include 'trackf90.cmn'
*KEND.

      integer i

      IF (IPINCIRC.NE.0.OR.IF1DIM.EQ.2) THEN

        PINW=2.D0*PINR
        PINH=2.D0*PINR

      ENDIF !IPINCIRC

      IF (PINCEN(2).EQ.9999.) THEN
        IF     (IPBRILL.EQ.0) THEN
          PINCEN(2)=0.0
        ELSE IF (IPBRILL.EQ.1) THEN
          PINCEN(2)=PINH/2.D0
        ELSE IF (IPBRILL.EQ.2) THEN
          PINCEN(2)=PINH/2.D0
        ELSE IF (IPBRILL.EQ.3) THEN
          PINCEN(2)=-PINH/2.D0
        ELSE IF (IPBRILL.EQ.4) THEN
          PINCEN(2)=-PINH/2.D0
        ENDIF
      ELSE IF (PINCEN(2).EQ.-9999.) THEN
        PINCEN(2)=YSTART+VYIN/VXIN*(PINCEN(1)-XSTART)
      ELSE IF (PINCEN(2).EQ.-8888.) THEN
        PINCEN(2)=yoffstr+yslopetr*PINCEN(1)
       else IF (PINCEN(2).EQ.-9000.) THEN
         pincen(2)=0.0d0
         do i=1,nco
           PINCEN(2)=pincen(2)+
     &       wtra(2,1,i)+wtra(2,2,i)/wtra(1,2,i)*(PINCEN(1)-wtra(1,1,i))
         enddo
         pincen(2)=pincen(2)/nco
      ENDIF

      IF (PINCEN(3).EQ.9999.) THEN
        IF     (IPBRILL.EQ.0) THEN
          PINCEN(3)=0.0
        ELSE IF (IPBRILL.EQ.1) THEN
          PINCEN(3)=PINW/2.D0
        ELSE IF (IPBRILL.EQ.2) THEN
          PINCEN(3)=-PINW/2.D0
        ELSE IF (IPBRILL.EQ.3) THEN
          PINCEN(3)=-PINW/2.D0
        ELSE IF (IPBRILL.EQ.4) THEN
          PINCEN(3)=PINW/2.D0
        ENDIF
      ELSE IF (PINCEN(3).EQ.-9999.) THEN
        PINCEN(3)=ZSTART+VZIN/VXIN*(PINCEN(1)-XSTART)
      ELSE IF (PINCEN(3).EQ.-8888.) THEN
        PINCEN(3)=zoffstr+zslopetr*PINCEN(1)
      else IF (PINCEN(3).EQ.-9000.) THEN
        pincen(3)=0.0d0
        do i=1,nco
          PINCEN(3)=pincen(3)+
     &      wtra(3,1,i)+wtra(3,2,i)/wtra(1,2,i)*(PINCEN(1)-wtra(1,1,i))
        enddo
        pincen(3)=pincen(3)/nco
      ENDIF

      RETURN
      END
