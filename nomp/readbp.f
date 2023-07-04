*CMZ :  3.05/10 13/08/2018  14.40.26  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.48/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  16.03.23  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  13.48.26  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  17.11.07  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.53.27  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.09  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE READBP
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

C     READBP LIEST PANDIRA OUTPUT IN ARRAYS EIN

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEND.

*KEEP,pandir.
      include 'pandir.cmn'
*KEND.

      INTEGER LTOPR, KTOPR,LLO,LU,LO,KL,KR,KMIN,KTOP,LMIN,LTOP,LL
     &         ,LDUM,KDUM,IREAD,L,K,IL,IK
      REAL X,Y,BX,BY

C      DATA BPAN/0./

      NXPAN=NXPANP
      NYPAN=NYPANP


      OPEN (
     1 UNIT=LUNP,
     1   FILE = FILEP,
     1   STATUS = 'OLD',
     1   FORM = 'FORMATTED')

      READ(LUNP,'(1A32)') PANTIT
      READ(LUNP,*) KMIN,KTOP,LMIN,LTOP

      IF (KTOP-KMIN.LT.1 .OR. 2*KTOP-1.GT.NXPAN .OR.
     &      LTOP-LMIN.LT.1 .OR. 2*LTOP-1.GT.NYPAN) THEN
         WRITE(6,*)
         WRITE(6,*) '*** S/R READBP:'
         WRITE(6,*)
     & 'KMIN.LT.3 .OR. 2*KTOP-1.GT.NXPAN .OR. LMIN.LT.3 .OR. 2*LTOP-1.GT.NYPAN'
         WRITE(16,*)
         WRITE(16,*) '*** S/R READBP:'
         WRITE(16,*)
     & 'KMIN.LT.3 .OR. 2*KTOP-1.GT.NXPAN .OR. LMIN.LT.3 .OR. 2*LTOP-1.GT.NYPAN'
         STOP
      ENDIF

C--- H-SHAPED MAGNET SYMMETRIE

      DO IL=1,2*LTOP-1
         DO IK=1,2*KTOP-1
         LPANB(IK,IL)=-1
         ENDDO
      ENDDO

      K=0
      L=1
      DO IREAD=1,(KTOP-KMIN+1)*(LTOP-LMIN+1)
         READ (LUNP,*,END=900)
     &            KDUM,LDUM,X,Y,BX,BY

C--- SCALIEREN (17.07.91)
C        BX=BX*BPANSC
C        BY=BY*BPANSC
C---
              K=K+1
         IF(K.GT.KTOP) THEN
             K=1
             L=L+1
         ENDIF

         KR=KTOP-1+K
         KL=KTOP+1-K
         LO=LTOP-1+L
         LU=LTOP+1-L

         BXPAN(KR,LO)=BX
         BYPAN(KR,LO)=BY
         XPAN(KR,LO)=X
         YPAN(KR,LO)=Y
         LPANB(KR,LO)=1

         BXPAN(KL,LO)=-BX  !VORZEICHEN KORRIGIERT 3.4.91
         BYPAN(KL,LO)=BY
         XPAN(KL,LO)=-X
         YPAN(KL,LO)=Y
         LPANB(KL,LO)=1

         BXPAN(KR,LU)=-BX  !VORZEICHEN KORRIGIERT 3.4.91
         BYPAN(KR,LU)=BY
         XPAN(KR,LU)=X
         YPAN(KR,LU)=-Y
         LPANB(KR,LU)=1

         BXPAN(KL,LU)=BX
         BYPAN(KL,LU)=BY
         XPAN(KL,LU)=-X
         YPAN(KL,LU)=-Y
         LPANB(KL,LU)=1
      ENDDO

900   CLOSE(LUNP)

      KTOP=2*KTOP-1
      LTOP=2*LTOP-1

C--- NORMIEREN

CN    DO IL=1,LTOP
C        IL=LTOP/2+1 !CN
C        DO IK=1,KTOP
C           BPAN=AMAX1(BPAN,
C     & SQRT(BYPAN(IK,IL)*BYPAN(IK,IL)+BXPAN(IK,IL)*BXPAN(IK,IL)))
C        ENDDO
CN    ENDDO
C     BP=BPANMX/BPAN
C     DO IL=1,LTOP
C        DO IK=1,KTOP
C        BXPAN(IK,IL)=BXPAN(IK,IL)*BP
C        BYPAN(IK,IL)=BYPAN(IK,IL)*BP
C        ENDDO
C     ENDDO

      WRITE(6,*)
      WRITE(6,*)'*** S/R READBP: B-FELD VON PANDIRA GELESEN ***'
      WRITE(6,*)
      WRITE(16,*)
      WRITE(16,*)'*** S/R READBP: B-FELD VON PANDIRA GELESEN ***'
      WRITE(16,*)
      WRITE(16,*)'PANDIRA-TITLE:',PANTIT
      WRITE(16,*)
     &'VIERER-SYMMETRIE VORAUSGESETZT, KTOP,LTOP VERDOPPELT'
      WRITE(16,*)'IREAD,KMIN,KTOP,LMIN,LTOP:',IREAD,KMIN,KTOP,LMIN,LTOP
      WRITE(16,*)'BX(KMIN,LMIN), BY(KMIN,LMIN):',
     &              BXPAN(KMIN,LMIN),BYPAN(KMIN,LMIN)

C-------12.12.90, REDUKTION VON BXPAN,BYPAN-> BXPANR,BYPANR

      DO 100 L=1,LTOP
          LLO=L
          DO 101 K=1,KTOP
                IF (LPANB(K,L).LT.0) GOTO 100
101       CONTINUE
            GOTO 200
100   CONTINUE
200   CONTINUE

      DO 300 L=LLO,LTOP-LLO+1
          LL=L-LLO+1
          DO 301 K=1,KTOP
         BXPANR(K,LL)=BXPAN(K,L)
         BYPANR(K,LL)=BYPAN(K,L)
         XPANR(K,LL)=XPAN(K,L)
         YPANR(K,LL)=YPAN(K,L)
301       CONTINUE
300   CONTINUE

      KTOPR=KTOP
      LTOPR=LTOP-2*LLO+2

      IF (LTOPR.LT.3) THEN
         WRITE(6,*)
         WRITE(6,*) '*** S/R READBP: LTOPR.LT.3 ***'
         WRITE(6,*)
         WRITE(16,*)
         WRITE(16,*) '*** S/R READBP: LTOPR.LT.3 ***'
         WRITE(16,*)
         STOP
      ENDIF


      WRITE(6,*)
      WRITE(6,*)'*** S/R REABDP: BXPAN,BYPAN REDUZIERT ***'
      WRITE(6,*)'KTOPR,LTOPR: ',KTOPR,LTOPR
      WRITE(16,*)
      WRITE(16,*)'*** S/R REABDP: BXPAN,BYPAN REDUZIERT ***'
      WRITE(16,*)'KTOPR,LTOPR: ',KTOPR,LTOPR

C     DO L=1,LTOPR
C     DO K=1,KTOPR
C        WRITE(88,1000)K,XPANR(K,L),YPANR(K,L),BXPANR(K,L),BYPANR(K,L)
C     ENDDO
C     WRITE(88,*)
C     WRITE(88,*)
C     WRITE(88,*)
C     ENDDO
C1000    FORMAT(I5,4F10.6)

C     WRITE(16,*)
C     WRITE(16,*)'SCALENFAKTOR FUER PANDIRA B-FELD:',BPANSC
C     WRITE(16,*)

      RETURN
      END
