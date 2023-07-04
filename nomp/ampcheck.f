*CMZ :  2.70/12 01/03/2013  15.45.11  by  Michael Scheer
*CMZ :  2.70/00 08/11/2012  11.33.18  by  Michael Scheer
*CMZ :  2.69/02 08/11/2012  10.17.41  by  Michael Scheer
*CMZ :  2.66/07 24/02/2010  17.18.39  by  Michael Scheer
*CMZ :  2.37/07 11/12/2001  17.50.58  by  Michael Scheer
*CMZ :  2.36/00 08/11/2001  14.16.25  by  Michael Scheer
*CMZ :  2.15/00 02/11/2001  14.21.12  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.16.24  by  Michael Scheer
*CMZ :  2.00/02 12/01/99  15.53.04  by  Michael Scheer
*CMZ :  2.00/00 07/01/99  11.21.39  by  Michael Scheer
*-- Author :    Michael Scheer   06/01/99

      SUBROUTINE AMPCHECK(
     &  NSOURCE,NOBSV,NFREQ,IFREQ2P,
     &  NOBSVZ,NOBSVY,MOBSVZ,MOBSVY,
     &  MEDGEZ,MEDGEY,MMEDGEZ,MMEDGEY,
     &  PINW,PINH,PINR,IPIN,IF1DIM,IPINCIRC,AMPFREQ,IAMPREP,
     &  ibunch,iubunch,bunchlen,
     &  IERROR)

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

      INTEGER ICAL,IERROR ,ICHECK,IAMPREP

      INTEGER
     &  NSOURCE,NOBSV,NFREQ,IFREQ2P,
     &  NOBSVZ,NOBSVY,MOBSVZ,MOBSVY,
     &  MEDGEZ,MEDGEY,MMEDGEZ,MMEDGEY,
     &  IF1DIM,IPIN,IPINO,IPINCIRC,
     &  ibunch,iubunch

      DOUBLE PRECISION
     &  PINW,PINH,PINR,AMPFREQ,AMPFREQO,bunchlen

      INTEGER
     &  NSOURCEO,NOBSVO,NFREQO,IFREQ2PO,
     &  NOBSVZO,NOBSVYO,MOBSVZO,MOBSVYO,
     &  MEDGEZO,MEDGEYO,MMEDGEZO,MMEDGEYO,
     &  IF1DIMO,IPINCIRCO,
     &  iubuncho,ibuncho

      DOUBLE PRECISION
     &  PINWO,PINHO,PINRO,RCHECK,ACHECK,bunchleno


      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN

        NSOURCEO=NSOURCE
        NOBSVO=NOBSV
        NFREQO=NFREQ
        IFREQ2PO=IFREQ2P
        NOBSVZO=NOBSVZ
        NOBSVYO=NOBSVY
        MOBSVZO=MOBSVZ
        MOBSVYO=MOBSVY
        MEDGEZO=MEDGEZ
        MEDGEYO=MEDGEY
        MMEDGEZO=MMEDGEZ
        MMEDGEYO=MMEDGEY
        IF1DIMO=IF1DIM
        IPINCIRCO=IPINCIRC
        IPINO=IPIN
        ibuncho=ibunch
        iubuncho=iubunch
        bunchleno=bunchlen

        PINWO=PINW
        PINHO=PINH
        PINRO=PINR

        AMPFREQO=AMPFREQ

        ICAL=1
        IERROR=0

      ELSE  !ICAL


        ICHECK=
     &    +ABS(       NSOURCEO-NSOURCE
     &    )+ABS(      NOBSVO-NOBSV
     &    )+ABS(      NFREQO-NFREQ
     &    )+ABS(      IFREQ2PO-IFREQ2P
     &    )+ABS(      NOBSVZO-NOBSVZ
     &    )+ABS(      NOBSVYO-NOBSVY
     &    )+ABS(      MOBSVZO-MOBSVZ
     &    )+ABS(      MOBSVYO-MOBSVY
     &    )+ABS(      MEDGEZO-MEDGEZ
     &    )+ABS(      MEDGEYO-MEDGEY
     &    )+ABS(      MMEDGEZO-MMEDGEZ
     &    )+ABS(      MMEDGEYO-MMEDGEY
     &    )+ABS(      IF1DIMO-IF1DIM
     &    )+ABS(      IPINO-IPIN
     &    )+ABS(      IPINCIRCO-IPINCIRC)

        RCHECK=
     &    +DABS(      PINWO-PINW
     &    )+DABS(     PINHO-PINH
     &    )+DABS(     PINRO-PINR)

        IF (IAMPREP.GE.0.and.ampfreqo.ne.0d0) THEN
          ACHECK=DABS(AMPFREQO-AMPFREQ)/AMPFREQO
        ELSE
          ACHECK=0.0D0
        ENDIF

        IF (ICHECK.NE.0.OR.RCHECK.NE.0.D0.OR.ACHECK.GT.1.D-6) IERROR=-1

        if (ibuncho.ne.0.and.bunchleno.le.0.0d0) then
          write(16,*)' '
          write(16,*)'      *** Warning in AMPCHECK:'
          write(16,*)' '
          write(16,*)
     &      '      IBUNCH has not been zero on the file read. '
          write(16,*)'      Be careful, especially if BUNCHLEN has not been zero.'
          write(16,*)'      IBUNCH, IUBUNCH and BUNCHLEN read from file:'
          write(16,*)'      ',ibunch, iubunch, bunchlen
          write(16,*)' '
          write(6,*)' '
          write(6,*)'      *** Warning in AMPCHECK:'
          write(6,*)' '
          write(6,*)
     &      '      IBUNCH has not been zero on the file read. '
          write(6,*)'      Be careful, especially if BUNCHLEN has not been zero.'
          write(6,*)'      IBUNCH, IUBUNCH and BUNCHLEN read from file:'
          write(6,*)'      ',ibunch, iubunch, bunchlen
        endif

      ENDIF !ICAL

      RETURN
      END
