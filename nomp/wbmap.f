*CMZ :  4.00/15 07/03/2022  16.15.45  by  Michael Scheer
*CMZ :  4.00/14 22/12/2021  16.48.08  by  Michael Scheer
*CMZ :  4.00/04 17/05/2019  14.22.20  by  Michael Scheer
*CMZ :  3.05/05 13/07/2018  09.24.44  by  Michael Scheer
*CMZ :  3.03/04 19/10/2017  14.59.02  by  Michael Scheer
*CMZ :  3.02/03 23/10/2014  13.43.13  by  Michael Scheer
*CMZ :  3.00/00 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  14.04.56  by  Michael Scheer
*CMZ :  2.68/05 01/10/2012  13.51.04  by  Michael Scheer
*CMZ :  2.68/02 02/07/2012  12.58.24  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.55.58  by  Michael Scheer
*CMZ :  2.45/03 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.44/01 10/12/2002  17.49.10  by  Michael Scheer
*CMZ :  2.44/00 07/11/2002  15.12.01  by  Michael Scheer
*CMZ :  2.42/04 14/09/2002  07.16.48  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.02  by  Michael Scheer
*CMZ :  2.39/01 15/01/2002  16.47.24  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  1.03/06 06/08/98  18.01.11  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.28  by  Michael Scheer
*CMZ : 00.02/03 23/01/97  17.29.29  by  Michael Scheer
*CMZ : 00.02/00 25/11/96  09.51.41  by  Michael Scheer
*-- Author :    Michael Scheer   28/09/95
      SUBROUTINE WBMAP
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

C--- TO WRITE 3D FIELD MAP TO FILE


      IMPLICIT NONE

      EXTERNAL DCOSD,DSIND
      DOUBLE PRECISION DCOSD,DSIND

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,bpoly3d.
      include 'bpoly3d.cmn'
*KEEP,bpoly2dh.
      include 'bpoly2dh.cmn'
*KEEP,bpoly3dg.
      include 'bpoly3dg.cmn'
*KEEP,bmap.
      include 'bmap.cmn'
*KEEP,datetime.
      include 'datetime.cmn'
*KEND.

      DOUBLE PRECISION DX,DY,DZ,BX,BY,BZ,AX,AY,AZ,X,Y,Z,R,PHI

      INTEGER IPOI,IX,IY,IZ

      INTEGER NTUP_P,ICYCLE
      PARAMETER (NTUP_P=6)
      REAL*8 TUP_D(NTUP_P)
      CHARACTER(3) CHTAGS_D(NTUP_P)
      character(2048) cline

      data chtags_d/'x','y','z','bx','by','bz'/
      DATA DX,DY,DZ/0.D0,0.D0,0.D0/

      IF (IWBMAP.EQ.4) THEN
        CALL WBMAP4
        RETURN
      ELSE IF (IWBMAP.EQ.5) THEN
        CALL WBMAP5
      ELSE IF (IWBMAP.EQ.7) THEN
        CALL WBMAP_FOR_SPECTRA
        return
      ENDIF

      IF (IWBMAP.EQ.3) THEN
        OPEN(UNIT=LUNBMAP,FILE=FILEBMAP,STATUS='NEW'
     &    ,FORM='UNFORMATTED')
      ELSE IF (IWBMAP.eq.6) THEN
        open(unit=lunbmap,file=filebmap,status='new',form='formatted')
        call date_and_time(dtday,dttime,dtzone,idatetime)

        write(cline,*)'! WAVE: x y z Bx By Bz with x as long. beam axis'
        write(lunbmap,'(a)')cline(2:len_trim(cline))
        write(cline,*)'@ date (yyyy.month.day) and time = ',
     &    dtday(1:4),'.',dtday(5:6),'.',dtday(7:8),' ',
     &    dttime(1:2),':',dttime(3:4),':',dttime(5:6)
        write(lunbmap,'(a)')cline(2:len_trim(cline))
        write(cline,*)'@ run =  ',icode
        write(lunbmap,'(a)')cline(2:len_trim(cline))
        write(cline,*)'@ comment = ',code
        write(lunbmap,'(a)')cline(2:len_trim(cline))
        write(cline,*)'@ scaling = 1.0 1.0 1.0 1.0 1.0 1.0'
        write(lunbmap,'(a)')cline(2:len_trim(cline))
        write(cline,*)'@ offset = 0.0, 0.0, 0.0 0.0 0.0 0.0'
        write(lunbmap,'(a)')cline(2:len_trim(cline))

      else if (iwbmap.gt.0) then
        open(unit=lunbmap,file=filebmap,status='new'
     &    ,form='formatted')
      ENDIF

      IF (IBMRADIAL.NE.0) THEN
        IF (XMAPMN.EQ.9999.) THEN
          WRITE(6,*)'*** WARNING WBMAP: DEFAULT FOR XMAPMN NOT USEFUL'
          WRITE(6,*)'(IBMRADIL IS NOT ZERO, SEE NAMELIST WBMAPN)'
          WRITE(LUNGFO,*)'*** WARNING WBMAP: DEFAULT FOR XMAPMN NOT USEFUL'
          WRITE(LUNGFO,*)'(IBMRADIL IS NOT ZERO, SEE NAMELIST WBMAPN)'
        ENDIF
        IF (XMAPMX.EQ.9999.) THEN
          WRITE(6,*)'*** WARNING WBMAP: DEFAULT FOR XMAPMX NOT USEFUL'
          WRITE(6,*)'(IBMRADIL IS NOT ZERO, SEE NAMELIST WBMAPN)'
          WRITE(LUNGFO,*)'*** WARNING WBMAP: DEFAULT FOR XMAPMX NOT USEFUL'
          WRITE(LUNGFO,*)'(IBMRADIL IS NOT ZERO, SEE NAMELIST WBMAPN)'
        ENDIF
        IF (YMAPMN.EQ.9999.) THEN
          WRITE(6,*)'*** WARNING WBMAP: DEFAULT FOR YMAPMN NOT USEFUL'
          WRITE(6,*)'(IBMRADIL IS NOT ZERO, SEE NAMELIST WBMAPN)'
          WRITE(LUNGFO,*)'*** WARNING WBMAP: DEFAULT FOR YMAPMN NOT USEFUL'
          WRITE(LUNGFO,*)'(IBMRADIL IS NOT ZERO, SEE NAMELIST WBMAPN)'
        ENDIF
        IF (YMAPMX.EQ.9999.) THEN
          WRITE(6,*)'*** WARNING WBMAP: DEFAULT FOR YMAPMX NOT USEFUL'
          WRITE(6,*)'(IBMRADIL IS NOT ZERO, SEE NAMELIST WBMAPN)'
          WRITE(LUNGFO,*)'*** WARNING WBMAP: DEFAULT FOR YMAPMX NOT USEFUL'
          WRITE(LUNGFO,*)'(IBMRADIL IS NOT ZERO, SEE NAMELIST WBMAPN)'
        ENDIF
        IF (XMAPMX.LT.0.D0.OR.XMAPMN.LT.0.D0) THEN
          WRITE(6,*)'*** ERROR WBMAP: XMAPMX.OR.XMAPMN .LT. 0'
          WRITE(6,*)'(SEE NAMELIST WBMAPN)'
          WRITE(LUNGFO,*)'*** ERROR WBMAP: XMAPMX.OR.XMAPMN .LT. 0'
          WRITE(LUNGFO,*)'(SEE NAMELIST WBMAPN)'
          STOP '*** PROGRAM WAVE ABORTED ***'
        ENDIF
      ENDIF !IBMRADIAL

      IF (XMAPMN.EQ.9999.) XMAPMN=XSTART
      IF (XMAPMX.EQ.9999.) XMAPMX=XSTOP

      IF (XMAPMX.LT.XMAPMN) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN WBMAP: XMAPMX.LT.XMAPMN'
        WRITE(LUNGFO,*)'CHECK NAMELIST WBMAP IN WAVE.IN'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN WBMAP: XMAPMX.LT.XMAPMN'
        WRITE(6,*)'CHECK NAMELIST WBMAP IN WAVE.IN'
        WRITE(6,*)
        STOP '*** PROGRAM WAVE ABORTED ***'
      ENDIF

      IF (YMAPMX.LT.YMAPMN) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN WBMAP: YMAPMX.LT.YMAPMN'
        WRITE(LUNGFO,*)'CHECK NAMELIST WBMAP IN WAVE.IN'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN WBMAP: YMAPMX.LT.YMAPMN'
        WRITE(6,*)'CHECK NAMELIST WBMAP IN WAVE.IN'
        WRITE(6,*)
        STOP '*** PROGRAM WAVE ABORTED ***'
      ENDIF

      IF (ZMAPMX.LT.ZMAPMN) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN WBMAP: ZMAPMX.LT.ZMAPMN'
        WRITE(LUNGFO,*)'CHECK NAMELIST WBMAP IN WAVE.IN'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN WBMAP: ZMAPMX.LT.ZMAPMN'
        WRITE(6,*)'CHECK NAMELIST WBMAP IN WAVE.IN'
        WRITE(6,*)
        STOP '*** PROGRAM WAVE ABORTED ***'
      ENDIF

      IF (NMAPX.EQ.-9999) NMAPX=NINT((XSTOP-XSTART)*MYINUM)+1

      IF (NMAPX.GT.1)   DX=(XMAPMX-XMAPMN)/(NMAPX-1)
      IF (NMAPY.GT.1) DY=(YMAPMX-YMAPMN)/(NMAPY-1)
      IF (NMAPZ.GT.1) DZ=(ZMAPMX-ZMAPMN)/(NMAPZ-1)

      IF (IHBPOLY3D.NE.0.OR.IHBPOLY2DH.NE.0)
     &  CALL hbookm(NIDBPOLY-1,'WBMAP',NTUP_P,'//WAVE',nmapx*nmapy*nmapz,
     &  CHTAGS_D)

      IF (IWBMAP.NE.3.and.iwbmap.ne.6) THEN

        IPOI=0
        DO IZ=1,NMAPZ
          DO IY=1,NMAPY
            DO IX=1,NMAPX
              IPOI=IPOI+1
              IF (IBMRADIAL.EQ.0) THEN
                X=XMAPMN+(IX-1)*DX
                Y=YMAPMN+(IY-1)*DY
                Z=ZMAPMN+(IZ-1)*DZ
              ELSE IF (IBMRADIAL.EQ.1) THEN
                R=XMAPMN+(IX-1)*DX
                PHI=YMAPMN+(IY-1)*DY
                X=BMRADX0+R*DCOSD(PHI)
                Y=BMRADY0+R*DSIND(PHI)
                Z=ZMAPMN+(IZ-1)*DZ
              ELSE
                R=ZMAPMN+(IZ-1)*DZ
                PHI=YMAPMN+(IY-1)*DY
                Z=BMRADX0+R*DCOSD(PHI)
                Y=BMRADY0+R*DSIND(PHI)
                X=XMAPMN+(IX-1)*DX
              ENDIF
              CALL MYBFELD(X,Y,Z,BX,BY,BZ,AX,AY,AZ)
              IF (IBMAPX.EQ.0) THEN
                BX=-9999.D0
              ENDIF
              IF (IBMAPY.EQ.0) THEN
                BY=-9999.D0
              ENDIF
              IF (IBMAPZ.EQ.0) THEN
                BZ=-9999.D0
              ENDIF

              IF (IWBMAP.EQ.1) THEN
                WRITE(LUNBMAP,*)X,IPOI
                WRITE(LUNBMAP,*)Y
                WRITE(LUNBMAP,*)Z
                WRITE(LUNBMAP,*)BX
                WRITE(LUNBMAP,*)BY
                WRITE(LUNBMAP,*)BZ
                WRITE(LUNBMAP,*)
              ELSE IF (IWBMAP.EQ.2) THEN
                WRITE(LUNBMAP,'(6(1PD21.12))')X,Y,Z,BX,BY,BZ
              ENDIF

              IF (IHBPOLY3D.NE.0.OR.IHBPOLY2DH.NE.0) THEN
                TUP_D(1)=X
                TUP_D(2)=Y
                TUP_D(3)=Z
                TUP_D(4)=BX
                TUP_D(5)=BY
                TUP_D(6)=BZ
                CALL hfm(NIDBPOLY-1,TUP_D)
              ENDIF

            ENDDO
          ENDDO
        ENDDO

      ELSE IF (IWBMAP.EQ.6) THEN

        IPOI=0

        DO IX=1,NMAPX
          DO IY=1,NMAPY
            DO IZ=1,NMAPZ

              IPOI=IPOI+1

              IF (IBMRADIAL.EQ.0) THEN
                X=XMAPMN+(IX-1)*DX
                Y=YMAPMN+(IY-1)*DY
                Z=ZMAPMN+(IZ-1)*DZ
              ELSE IF (IBMRADIAL.EQ.1) THEN
                R=XMAPMN+(IX-1)*DX
                PHI=YMAPMN+(IY-1)*DY
                X=BMRADX0+R*DCOSD(PHI)
                Y=BMRADY0+R*DSIND(PHI)
                Z=ZMAPMN+(IZ-1)*DZ
              ELSE
                R=ZMAPMN+(IZ-1)*DZ
                PHI=YMAPMN+(IY-1)*DY
                Z=BMRADX0+R*DCOSD(PHI)
                Y=BMRADY0+R*DSIND(PHI)
                X=XMAPMN+(IX-1)*DX
              ENDIF
              CALL MYBFELD(X,Y,Z,BX,BY,BZ,AX,AY,AZ)
              IF (IBMAPX.EQ.0) THEN
                BX=-9999.D0
              ENDIF
              IF (IBMAPY.EQ.0) THEN
                BY=-9999.D0
              ENDIF
              IF (IBMAPZ.EQ.0) THEN
                BZ=-9999.D0
              ENDIF

              WRITE(LUNBMAP,'(6(1PE21.12))')X,Y,Z,BX,BY,BZ

              IF (IHBPOLY3D.NE.0.OR.IHBPOLY2DH.NE.0) THEN
                TUP_D(1)=X
                TUP_D(2)=Y
                TUP_D(3)=Z
                TUP_D(4)=BX
                TUP_D(5)=BY
                TUP_D(6)=BZ
                CALL hfm(NIDBPOLY-1,TUP_D)
              ENDIF

            ENDDO
          ENDDO
        ENDDO

      ELSE IF (IWBMAP.EQ.3) THEN

        WRITE(LUNBMAP)ICODE,CODE
        WRITE(LUNBMAP)NMAPX,XMAPMN,XMAPMX
        WRITE(LUNBMAP)NMAPY,YMAPMN,YMAPMX
        WRITE(LUNBMAP)NMAPZ,ZMAPMN,ZMAPMX

        IPOI=0

        DO IX=1,NMAPX
          DO IY=1,NMAPY
            DO IZ=1,NMAPZ

              IPOI=IPOI+1

              IF (IBMRADIAL.EQ.0) THEN
                X=XMAPMN+(IX-1)*DX
                Y=YMAPMN+(IY-1)*DY
                Z=ZMAPMN+(IZ-1)*DZ
              ELSE IF (IBMRADIAL.EQ.1) THEN
                R=XMAPMN+(IX-1)*DX
                PHI=YMAPMN+(IY-1)*DY
                X=BMRADX0+R*DCOSD(PHI)
                Y=BMRADY0+R*DSIND(PHI)
                Z=ZMAPMN+(IZ-1)*DZ
              ELSE
                R=ZMAPMN+(IZ-1)*DZ
                PHI=YMAPMN+(IY-1)*DY
                Z=BMRADX0+R*DCOSD(PHI)
                Y=BMRADY0+R*DSIND(PHI)
                X=XMAPMN+(IX-1)*DX
              ENDIF
              CALL MYBFELD(X,Y,Z,BX,BY,BZ,AX,AY,AZ)
              IF (IBMAPX.EQ.0) THEN
                BX=-9999.D0
              ENDIF
              IF (IBMAPY.EQ.0) THEN
                BY=-9999.D0
              ENDIF
              IF (IBMAPZ.EQ.0) THEN
                BZ=-9999.D0
              ENDIF

              WRITE(LUNBMAP)SNGL(BX),SNGL(BY),SNGL(BZ)

              IF (IHBPOLY3D.NE.0.OR.IHBPOLY2DH.NE.0) THEN
                TUP_D(1)=X
                TUP_D(2)=Y
                TUP_D(3)=Z
                TUP_D(4)=BX
                TUP_D(5)=BY
                TUP_D(6)=BZ
                CALL hfm(NIDBPOLY-1,TUP_D)
              ENDIF

            ENDDO
          ENDDO
        ENDDO

      ENDIF !IWBMAP.EQ.3

      IF (IWBMAP.GT.0) CLOSE (LUNBMAP)

      IF (IWBMAP.GT.0) THEN

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'     SR WBMAP: FIELD MAP WRITTEN TO FILE'
        WRITE(LUNGFO,*)'     ',FILEBMAP
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'     SR WBMAP: FIELD MAP WRITTEN TO FILE'
        WRITE(6,*)'     ',FILEBMAP
        WRITE(6,*)

      ENDIF

      IF (IHBPOLY3D.NE.0.OR.IHBPOLY2DH.NE.0) THEN
        CALL MHROUT(NIDBPOLY-1,ICYCLE,' ')
        CALL hdeletm(NIDBPOLY-1)
      ENDIF

      RETURN
      END
