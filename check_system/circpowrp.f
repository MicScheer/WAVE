*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.34/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.31/00 23/04/2001  18.27.11  by  Michael Scheer
*CMZ :  2.20/09 23/03/2001  11.01.07  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  14.08.07  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.08.20  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.22.27  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.04.41  by  Michael Scheer
*CMZ :  1.00/00 31/07/97  17.36.12  by  Michael Scheer
*CMZ : 00.02/05 18/03/97  15.48.44  by  Michael Scheer
*CMZ : 00.02/04 26/02/97  12.07.43  by  Michael Scheer
*CMZ : 00.02/00 10/12/96  18.07.52  by  Michael Scheer
*CMZ : 00.01/09 01/09/95  12.58.16  by  Michael Scheer
*CMZ : 00.01/02 24/11/94  15.49.44  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.38  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.04  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE CIRCPOWRP(ISOUR)
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
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- INTEGRATES THE POWER DENSITY INSIDE THE PINHOLE

      IMPLICIT NONE


*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.


      INTEGER ISOUR,IY,IZ,IOBSV,ICAL,IR,IP,NR,MR,MP,KDUM
      INTEGER IWBLEN
      INTEGER IWRPHIS,IWSOUR
      INTEGER IW_BLENO

      DOUBLE PRECISION DSUM,RPHI,SUM

      DOUBLE PRECISION SUMY(NDOBSVZP)
      DOUBLE PRECISION R(NDOBSVZP),PHI(NDOBSVYP)
     &        ,FPHI(NDOBSVYP)
     &        ,SZ(NDOBSVZP),SY(NDOBSVYP)
     &        ,FZ(NDOBSVZP),FY(NDOBSVYP),DPHI,DR,X,Y


      DATA ICAL/0/

      IWBLEN=0
      IW_BLENO=IW_BLEN

C--- INTEGRATION OVER PHI

      IF (ICAL.EQ.0) THEN

          DR=DMIN1(OBSVDZ,OBSVDY)
          MR=NINT(PINR/DR)+1
          DR=PINR/(MR-1)
          MEDGER=MIN( MEDGEZ, MEDGEY)
          NR=MR+MEDGER
          MP=MOBSVY

          IF (MR.LT.1) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** ERROR IN CIRCPOWRP ***'
           WRITE(LUNGFO,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
           WRITE(LUNGFO,*)'INCREASE PARAMETER MPINZ IN NAMELIST PINHOLE'
           WRITE(LUNGFO,*)
           WRITE(6,*)
           WRITE(6,*)'*** ERROR IN CIRCPOWRP ***'
           WRITE(6,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
           WRITE(6,*)'INCREASE PARAMETER MPINZ IN NAMELIST PINHOLE'
           WRITE(6,*)
          ENDIF
          IF (MP.LT.4) THEN
           WRITE(LUNGFO,*)
           WRITE(LUNGFO,*)'*** ERROR IN CIRCPOWRP ***'
           WRITE(LUNGFO,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
           WRITE(LUNGFO,*)'INCREASE PARAMETER MPINY IN NAMELIST PINHOLE'
           WRITE(LUNGFO,*)
           WRITE(6,*)
           WRITE(6,*)'*** ERROR IN CIRCPOWRP ***'
           WRITE(6,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
           WRITE(6,*)'INCREASE PARAMETER MPINY IN NAMELIST PINHOLE'
           WRITE(6,*)
          ENDIF

          DPHI=2.D0*PI1/(MP-1)
          DO IP=1,MP
               PHI(IP)=(IP-1)*DPHI
          ENDDO

          DO IR=1,NR
               R(IR)=(IR-1)*DR
          ENDDO

      DO IR=2,NR
      DO IP=1,MP
          XC(IP+(IR-1)*NOBSVY)=R(IR)*DCOS(PHI(IP))+PINCEN(3)
          YC(IP+(IR-1)*NOBSVY)=R(IR)*DSIN(PHI(IP))+PINCEN(2)
      ENDDO !IP
      ENDDO !IR

          ICAL=1

      ENDIF !ICAL

C--- INTERPOLATION OF INTENSITY FOR CIRCULAR GRID

      DO IR=2,NR
          DO IP=1,MP

            X=XC(IP+(IR-1)*NOBSVY)
            Y=YC(IP+(IR-1)*NOBSVY)

         DO IY=1,NOBSVY
         DO IZ=1,NOBSVZ
             IOBSV=(IY-1)*NOBSVZ+IZ
             FZ(IZ)=SPECPOW(ISOUR+(IOBSV-1)*NSOURCE)
         ENDDO !IZ

       CALL FSPLINDX(OBSVDZ,FZ,NOBSVZ,0.D0,0.D0,SZ)
         CALL SPLINZY(NOBSVZ,X,FY(IY),OBSVZ,FZ,SZ,KDUM)

         ENDDO !IY

         CALL FSPLINDX(OBSVDY,FY,NOBSVY,0.D0,0.D0,SY)
         CALL SPLINZY(NOBSVY,Y,FPHIR(IP+(IR-1)*NOBSVY),OBSVY,FY,SY,KDUM)
         IF(IWRPHIS.NE.ISOUR.AND.
     &              FPHIR(IP+(IR-1)*NOBSVY).LT.0.0) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** WARNING SR CIRCPOWRP ***'
             WRITE(LUNGFO,*)
     &              'SPLINE INTERPOLATION FOR OPTION IRPHI FAILED, RESULTS NOT RELIABLE'
             WRITE(LUNGFO,*)'SOURCE POINT:',ISOUR
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             IWRPHIS=ISOUR
           IW_BLEN=1
         ENDIF !IWPHIR

          ENDDO !IP
      ENDDO !IR

C--- DO THE INTEGRATION OF FPHIR OVER PHI AND R

      SUMY(1)=0.0
      DO IR=2,NR  !FIRST RADIUS IS ZERO

          DO IP=1,MP
         FPHI(IP)=FPHIR(IP+(IR-1)*NOBSVY)
          ENDDO   !IP

          CALL FSPLPER(DPHI,FPHI,MP,SY)

          SUMY(IR)=0.D0
          RPHI=R(IR)*DPHI
          DO IP=1,MP-1

         DSUM=
     &       RPHI*0.5D0*(FPHI(IP)+FPHI(IP+1))
     &          -RPHI**3/24.D0*(SY(IP)+SY(IP+1))

         IF(IWSOUR.NE.ISOUR.AND.DSUM.LT.0.0) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** WARNING SR CIRCPOWRP ***'
             WRITE(LUNGFO,*)
     &              'SPLINE INTEGRATION FAILED, RESULTS NOT RELIABLE'
             WRITE(LUNGFO,*)'SOURCE POINT:',ISOUR
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             IWSOUR=ISOUR
           IW_BLEN=1
         ENDIF !IWSOUR

         SUMY(IR)=SUMY(IR)+DSUM

          ENDDO   !IP

      ENDDO !IR

      CALL FSPLINDX(DR,SUMY,NR,0.D0,0.D0,SZ)

          SUM=0.0
          DO IR=1,MR-1

         DSUM=
     &          DR*0.5D0
     &          *(SUMY(IR)+SUMY(IR+1))
     &          -DR**3/24.D0
     &          *(SZ(IR)+SZ(IR+1))

         IF(IWSOUR.NE.ISOUR.AND.DSUM.LT.0.0) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)'*** WARNING SR CIRCPOWRP ***'
             WRITE(LUNGFO,*)
     &              'SPLINE INTEGRATION FAILED, RESULTS NOT RELIABLE'
             WRITE(LUNGFO,*)'SOURCE POINT:',ISOUR
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
             IWSOUR=ISOUR
           IW_BLEN=1
         ENDIF !IWSOUR

         SUM=SUM+DSUM

          ENDDO

        SPECPOWVH(ISOUR)=SUM

      IW_BLEN=IW_BLENO

      RETURN
      END
