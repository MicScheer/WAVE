*CMZ :  2.41/08 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  17.44.31  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  1.00/00 30/06/97  11.38.41  by  Michael Scheer
*CMZ : 00.02/04 10/02/97  14.07.52  by  Michael Scheer
*CMZ : 00.02/03 04/02/97  16.50.14  by  Michael Scheer
*-- Author :    Michael Scheer   22/01/97
      SUBROUTINE VPOT2DH
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
     &(NPOI,IFAIL)

*KEEP,bpoly2dhf90u.
      include 'bpoly2dhf90u.cmn'
*KEND.

C---  TO FIT POTENTIAL WITH TRANSVERSAL POLYNOMIAL AND
C     LONGITUDINAL SIN/COS-LIKE ANSATZ
C     OF A MAGNETIC FIELD B=(BX,BY,BZ)=-GRAD(V)
C

C--- INPUT:

C     NPOI  : NUMBER OF DATA POINTS X,Y,Z,BX,BY,BZ

C--- OUTPUT:

C     IFAIL : FAILURE FLAG

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpoly2dh.
      include 'bpoly2dh.cmn'
*KEND.

      INTEGER NPOI


      INTEGER IFAIL,NARG,NFUN,IPOI,IPAR,NPARP
      PARAMETER (NPARP=NPARTOTP,NFUN=3,NARG=3)

      DOUBLE PRECISION WS(NARG+NFUN)
      DOUBLE PRECISION PARAM(NPARP),A(NPARP,NPARP),T(NFUN,NPARP)

C NOTE CHANGE COORDINATE SYSTEMS!

      ALLOCATE(FUNDATA(NARG+NFUN,NPOI))
      DO IPOI=1,NPOI
         fundata(1,IPOI)=-z(IPOI)
         fundata(2,IPOI)=y(IPOI)
         fundata(3,IPOI)=x(IPOI)
         fundata(4,IPOI)=-BZ(IPOI)
         fundata(5,IPOI)=BY(IPOI)
         fundata(6,IPOI)=BX(IPOI)
      ENDDO

      CALL UTIL_LINEAR_FIT
     &  (IFAIL,NPARTOT,PARAM,NPOI,NPOI,NARG,NFUN,A,T,FUNDATA,WS)
      IF (IFAIL.NE.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING IN VPOT2DH: FIT FAILED'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'IFAIL:',IFAIL
            WRITE(LUNGFO,*) '(maybe no transversal field gradient'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING IN VPOT2DH: FIT FAILED'
          WRITE(6,*)
          WRITE(6,*)'IFAIL:',IFAIL
            WRITE(6,*) '(maybe no transversal field gradient'
          WRITE(6,*)
      ENDIF

      DO IPAR=1,NPARTOT
          PAR2DH(IPAR)=PARAM(IPAR)
      ENDDO

      DEALLOCATE(FUNDATA)
      RETURN
      END
