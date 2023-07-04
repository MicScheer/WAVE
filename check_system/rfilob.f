*CMZ :  4.00/11 09/06/2021  14.32.44  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.67/02 19/04/2012  08.52.32  by  Michael Scheer
*CMZ :  2.17/00 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  11.19.57  by  Michael Scheer
*CMZ : 00.02/00 19/11/96  14.56.18  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.17.17  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.53.48  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.41  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE RFILOB
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

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

C--- SUBROUTINE READS OBSERVATION POINTS FOR WHICH PHOTON FLUX IS CALCULATED

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEND.

      INTEGER IOB,IC,ipos(2,8),nwords,istat,ieof
      character(256) cline

      OBSVDZ=0.
      OBSVDY=0.

C--- ONE SINGLE OBSERVATION POINT

      IF (IRFILOB.EQ.0) THEN
        NOBSV=1
        NOBSVZ=1
        NOBSVY=1
        OBSV(1,1)=OBS1X
        OBSV(2,1)=OBS1Y
        OBSV(3,1)=OBS1Z
        RETURN
      ENDIF !IRFILOB

      OPEN (UNIT=LUNOB,FILE=FILEOB,STATUS='OLD',FORM='FORMATTED',ERR=999)

      call util_skip_comment_end(lunob,ieof)

      read(lunob,'(a)') cline
      call util_string_split(cline,8,nwords,ipos,istat)

      IF(IHPIN.NE.0) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'*** ERROR IN SR RFILOB ***'
         WRITE(LUNGFO,*)'CONTROL FLAG IHPIN IS NOT COMPATIBLE WITH'
         WRITE(LUNGFO,*)'IRFILOB'
         WRITE(6,*)
         WRITE(6,*)'*** ERROR IN SR RFILOB ***'
         WRITE(6,*)'CONTROL FLAG IHPIN IS NOT COMPATIBLE WITH'
         WRITE(6,*)'IRFILOB'
      ENDIF

      if (nwords.eq.1) then
        READ(LUNOB,*,ERR=99) NOBSV
      else
        nobsv=0
        rewind(lunob)
        do while (.true.)
          call util_skip_comment_end(lunob,ieof)
          if (ieof.ne.0) exit
          read(lunob,'(a)',end=9,err=9) cline
          nobsv=nobsv+1
        enddo
9       rewind(lunob)
      endif

      IF(NOBSV.GT.NDOBSV) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** ERROR IN SR RFILOB ***'
        WRITE(LUNGFO,*)'TOO MANY OBSERVATION POINTS ON FILEOB'
        WRITE(LUNGFO,*)'INCREASE PARAMETER NDOBSVP IN CMPARA.CMN'
        WRITE(6,*)
        WRITE(6,*)'*** ERROR IN SR RFILOB ***'
        WRITE(6,*)'TOO MANY OBSERVATION POINTS ON FILEOB'
        WRITE(6,*)'INCREASE PARAMETER NDOBSVP IN CMPARA.CMN'
        STOP
      ENDIF

      IF (IOBSV_A.NE.NOBSV) THEN
          IF (IOBSV_A.NE.0) DEALLOCATE(OBSV)
          ALLOCATE(OBSV(3,NOBSV))
          IOBSV_A=NOBSV
      ENDIF !(IOBSV_A.LT.NOBSV)
      IF (IOBSVZ_A.NE.NOBSVZ) THEN
          IF (IOBSVZ_A.NE.0) DEALLOCATE(OBSVZ)
          ALLOCATE(OBSVZ(NOBSVZ))
          IOBSVZ_A=NOBSVZ
      ENDIF !(IOBSVY_A.LT.NOBSVY)
      IF (IOBSVY_A.NE.NOBSVY) THEN
          IF (IOBSVY_A.NE.0) DEALLOCATE(OBSVY)
          ALLOCATE(OBSVY(NOBSVY))
          IOBSVY_A=NOBSVY
      ENDIF !(IOBSV_A.LT.NOBSV)

      DO IOB=1,NOBSV
          call util_skip_comment_end(lunob,ieof)
          READ(LUNOB,*,ERR=99) (OBSV(IC,IOB),IC=1,3)
      ENDDO
      CLOSE(LUNOB)

      RETURN

999   CONTINUE
      WRITE(6,*)
      WRITE(6,*)'*** ERROR IN SR RFILOB ***'
      WRITE(6,*)'FILE OPENING ERROR '
      WRITE(6,*)'FILE, UNIT:'
      WRITE(6,*)FILEOB, LUNOB
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'*** ERROR IN SR RFILOB ***'
      WRITE(LUNGFO,*)'FILE OPENING ERROR '
      WRITE(LUNGFO,*)'FILE, UNIT:'
      WRITE(LUNGFO,*)FILEOB, LUNOB
      STOP
99    CONTINUE
      WRITE(6,*)
      WRITE(6,*)'*** ERROR IN SR RFILOB ***'
      WRITE(6,*)'FILE READING ERROR '
      WRITE(6,*)'FILE, UNIT:'
      WRITE(6,*)FILEOB, LUNOB
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'*** ERROR IN SR RFILOB ***'
      WRITE(LUNGFO,*)'FILE READING ERROR '
      WRITE(LUNGFO,*)'FILE, UNIT:'
      WRITE(LUNGFO,*)FILEOB, LUNOB
      STOP
      END
