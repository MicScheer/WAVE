*CMZ :  3.00/00 11/03/2013  10.40.59  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/08 12/11/2009  16.27.11  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.53.11  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.27  by  Michael Scheer
*-- Author : Michael Scheer
C****************************************************************************
      SUBROUTINE POWSPL(IWALL,IMODE,IPOL,ISTART,N,YP1,YPN)
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

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,spect.
      include 'spect.cmn'
*KEND.


      INTEGER IPOL,IMODE,IWALL,K,N,I,ISTART,IS

       DOUBLE PRECISION QN,SIG,P,UN,YP1,YPN

      ALLOCATE(POWSPLU(N))

      IS=ISTART-1

      IF (N.LT.3) THEN
          POWS2(1)=YP1
          POWS2(N)=YPN
          GOTO 9999
      ENDIF


      IF (YP1.GT..99D30) THEN
        POWS2(1)=0.
        POWSPLU(1)=0.
      ELSE
        POWS2(1)=-0.5D0
        POWSPLU(1)=(3.D0/(RADPOW(IWALL,IPOL,IS+2)
     &  -RADPOW(IWALL,IPOL,IS+1)))
     &  *((RADPOW(IWALL+2+2*IMODE,IPOL,IS+2)
     &   -RADPOW(IWALL+2+2*IMODE,IPOL,IS+1))
     &  /(RADPOW(IWALL,IPOL,IS+2)
     &   -RADPOW(IWALL,IPOL,IS+1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(RADPOW(IWALL,IPOL,IS+I)-RADPOW(IWALL,IPOL,IS+I-1))
     &  /(RADPOW(IWALL,IPOL,IS+I+1)-RADPOW(IWALL,IPOL,IS+I-1))
        P=SIG*POWS2(I-1)+2.D0
        POWS2(I)=(SIG-1.)/P
        POWSPLU(I)=(
     &       6.D0*(
     &                (RADPOW(IWALL+2+2*IMODE,IPOL,IS+I+1)
     &                -RADPOW(IWALL+2+2*IMODE,IPOL,IS+I))
     &                /
     &                 (RADPOW(IWALL,IPOL,IS+I+1)-RADPOW(IWALL,IPOL,IS+I))
     &                -
     &                 (RADPOW(IWALL+2+2*IMODE,IPOL,IS+I)
     &                 -RADPOW(IWALL+2+2*IMODE,IPOL,IS+I-1))
     &                /
     &                 (RADPOW(IWALL,IPOL,IS+I)-RADPOW(IWALL,IPOL,IS+I-1))
     &            )
     &      /
     &       (RADPOW(IWALL,IPOL,IS+I+1)-RADPOW(IWALL,IPOL,IS+I-1))
     &      -
     &       SIG*POWSPLU(I-1)
     &      )/P
11    CONTINUE
      IF (YPN.GT..99D30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5D0
        UN=(3.D0/(RADPOW(IWALL,IPOL,IS+N)-RADPOW(IWALL,IPOL,IS+N-1)))
     &     *(YPN-(RADPOW(IWALL+2+2*IMODE,IPOL,IS+N)
     &     -RADPOW(IWALL+2+2*IMODE,IPOL,IS+N-1))
     &     /(RADPOW(IWALL,IPOL,IS+N)-
     &     RADPOW(IWALL,IPOL,IS+N-1)))
      ENDIF
      POWS2(N)=(UN-QN*POWSPLU(N-1))/(QN*POWS2(N-1)+1.D0)
      DO 12 K=N-1,1,-1
        POWS2(K)=POWS2(K)*POWS2(K+1)+POWSPLU(K)
12    CONTINUE

C-----------------------------------------------------------------

9999  CONTINUE

      DEALLOCATE(POWSPLU)

      RETURN
      END
