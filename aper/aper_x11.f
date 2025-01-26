*CMZ :          26/01/2025  19.11.35  by  Michael Scheer
*CMZ :  1.00/01 25/01/2025  14.19.36  by  Michael Scheer
*CMZ :  1.00/00 14/01/98  15.26.48  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.13  by  Michael Scheer
*-- Author : Michael Scheer
      PROGRAM APER_X11

      implicit none

C    PROGRAMM ZEICHNET APERTUR UND LICHTSTRAHLEN

      integer, PARAMETER :: NWORDS=200000
      real rpaw
      COMMON/PAWC/RPAW(NWORDS)

      integer, PARAMETER :: LINDIM=100

      character XLABFR(5),YLABFR(5),XLABAX(5)

      REAL*4 :: XLABF(1)=0.0,YLABF(1)=0.0,XLABA(1)=0.0,y0axis,x0axis,phi,rg,xfin,yfin

      EQUIVALENCE (XLABFR(1),XLABF(1)),(YLABFR(1),YLABF(1)),
     &  (XLABAX(1),XLABA(1))

*KEEP,aperseq.
      integer modeli,idev,modus,modeb,modeq,nangle

      real HEIGLI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY,
     &  blimit,wlspos,phimin,phimax

      COMMON /MAGNET/ HEIGLI,MODELI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY

      CHARACTER*80 DFILE

      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
     &  ,PHIMIN(10),PHIMAX(10),NANGLE
*KEEP,magseq.
      integer, parameter :: nmagsp=1000
      integer nmags
      real xplmags(2,5,nmagsp),widmag
      common/magseqc/nmags,xplmags,widmag

*KEND.

      integer modax,moday,luni,i
      real sheaxd,sheaxu,theaxd,theaxu,dirax,ymaxfr,theayr,theayl,theax,sheayr,sheayl,diray,
     &  sheax,xmint,xmaxt,xmaxfr,xlaxis,ylaxis

      character :: w,cblank=' '

      DATA SHEAXD/ -0.3/
      DATA SHEAXU/-0.3/
      DATA THEAXD/ 0.1/
      DATA THEAXU/-0.1/

      DATA YMINFR/0./          ! UNTER PLOTRAND AUF PAPIER (CM)
      DATA YMAXFR/20./      ! OBERER PLOTRAND AUF PAPIER (CM)

      DATA SHEAYL/0.3/
      DATA SHEAYR/+0.3/
      DATA THEAYL/-0.1/
      DATA THEAYR/+0.1/

      DATA DIRAX/0./
      DATA SHEAX/-0.3/
      DATA THEAX/0.2/
      DATA MODAX/ 9 /

      DATA DIRAY/0./
      DATA MODAY/ 0 /

      idev=1

C--- RAHMEN

      XMINFR=0.          ! LINKER PLOTRAND AUF PAPIER (CM)
      XMAXFR=20.         ! RECHTER PLOTRAND AUF PAPIER (CM)
      SMIAXX=-150.00 !COMMON      ! SKALENMINIMUM
      SMAAXX=1000.00        ! SKALEMMAXIMUM

      LUNI=99

      OPEN(UNIT=luni,FILE='aper.inp',STATUS='UNKNOWN')

C--- SKALA

      READ (LUNI,*) SMIAXX,SMAAXX
      READ (LUNI,*) XMINT,XMAXT

      SMIAXY=  -SMAAXX/2./(XMAXFR-XMINFR)*(YMAXFR-YMINFR)
      SMAAXY=   SMAAXX/2./(XMAXFR-XMINFR)*(YMAXFR-YMINFR)
      SMIAXY=  -SMAAXX/2./SQRT(2.)
      SMAAXY=   SMAAXX/2./SQRT(2.)

      SMIAXY=-30. ! HIER FESTE Y-SCALA WAEHLEN
      SMAAXY=+70.

      READ (LUNI,*) SMIAXY,SMAAXY

C--- z-width of magnets

      READ (LUNI,*)widmag

C--- X-POSITION DES WLS

      READ (LUNI,*)WLSPOS

      CLOSE(99)

C--- DATENFILE WAEHLEN

      DFILE='wave_track.dat'

C--- INITIALISIERUNG

      call mshplt_init(20,-20.,-20.,0,0,800,800,'aper.eps',cblank,cblank,0.0)

C--- RAHMEN

      XLAXIS=XMAXFR-XMINFR
      YLAXIS=YMAXFR-YMINFR

      CALL mPLFRA(SMIAXX,SMAAXX,SMIAXY,SMAAXY,' ')
c      call mshplt_set_text_angle(90.0)
c      CALL mshplt_text_NDC(-0.1,0.9,'cm')
c      call mshplt_set_text_angle(0.0)
c      CALL mshplt_text_NDC(0.95,-0.5,'cm')
      call mshplt_reset_clipping
      call mplax('cm','cm')

C--- ACHSE

      X0AXIS=XMINFR
      Y0AXIS=(YMAXFR-YMINFR)/2.+YMINFR

C--- MAGNETE ZEICHNEN
      CALL magseq
c      call lattice

C--- MODUS WAEHLEN

101   continue
      MODUS=0
      MODEQ=0
      MODEB=0
      BLIMIT=0

      WRITE(6,*) 'MODUS:'
      WRITE(6,*) '0 = STRAHLEN VON QUELLPUNKTEN ZEICHNEN'
      WRITE(6,*) '1 = STRAHLENKEGEL VOM PUNKT MIT MAX. B-FELD (OEFFNUNGSWINKEL 2* 3 MRAD)'
      WRITE(6,*) '2 = STRAHLEN AUS DEM BEREICH OBER- ODER UNTERHALB |B|limit'
      WRITE(6,*) '3 = STRAHLEN AUS WINKELBEREICH (Z.B. BEAMLINE-APERTUR)'
      READ(5,*) MODUS

      IF (MODUS.NE.0.AND.MODUS.NE.1.AND.MODUS.NE.2.AND.MODUS.NE.3) GOTO 101

      IF (MODUS.EQ.0) THEN
          WRITE(6,*)
          WRITE(6,*) 'Q-MODUS:'
202       WRITE(6,*)
          WRITE(6,*) '0 : LICHT VOM GESAMTEN QUELLPUNKT'
          WRITE(6,*) '1 : LICHT VOM ZENTRUM DES QUELLPUNKTES'
          READ(5,*)MODEQ
        ENDIF !MODUS

      IF (MODEQ.NE.0.AND.MODEQ.NE.1) GOTO 202

      IF (MODUS.EQ.2) THEN
102       WRITE(6,*) 'B-MODUS:'
          WRITE(6,*) '0 : |B| > ( |B|max * |B|limit)'
          WRITE(6,*) '1 :  B  > ( |B|max * |B|limit)'
          WRITE(6,*) '2 : -B  > (-|B|min * |B|limit)'
          WRITE(6,*)
          READ (5,*) MODEB
          IF (MODEB.NE.0.AND.MODEB.NE.1.AND.MODEB.NE.2) GOTO 102

            WRITE (6,*) 'MINIMALES B-FELD, DAS BERUECKSICHTIGT WERDEN SOLL'
          WRITE(6,*) 'EINGABE RELATIV D.H. AUF BMAX NORMIERT:'
          READ (5,*) BLIMIT
      ENDIF

      IF (MODUS.EQ.3) THEN
         WRITE(6,*)
         WRITE(6,*)'ANZAHL DER WINKELBEREICHE (MAX. 10):'
         READ(5,*) NANGLE
         IF (NANGLE.GT.10)
     &          STOP 'FEHLER: NUR 10 WINKELBEREICHE ERLAUBT'
         RG=180./3.141593
         DO I=1,NANGLE
            WRITE(6,*)'MINIMALER UND MAXIMALER WINKEL (mrad):'
            READ(5,*)PHIMIN(I),PHIMAX(I)
            PHIMIN(I)=PHIMIN(I)/1000.*RG
            PHIMAX(I)=PHIMAX(I)/1000.*RG
            IF (PHIMIN(I).GT.PHIMAX(I)) THEN
               PHI=PHIMIN(I)
               PHIMIN(I)=PHIMAX(I)
               PHIMAX(I)=PHI
            ENDIF
         ENDDO
      ENDIF

C--- SYNCHROTRONSTRAHLUNG ZEICHNEN
      call light(XMINT,XMAXT)

C--- TERMINIEREN

      CALL mUWK(0,0)
      W='W'
      WRITE(6,*)
      WRITE(6,*)'WEITER ODER STOP (W/*) [W]:'
      READ(5,'(A1)') W
      IF (W.EQ.'W' .OR. W.EQ.'w' .OR. W.EQ.' ') GOTO 101

      call mplend

      XFIN=99999.
      YFIN=99999.

      end

      include 'mshplt.f'
*CMZ :  1.00/01 23/01/2025  16.31.42  by  Michael Scheer
*CMZ :  1.00/00 14/01/98  12.20.43  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.14  by  Michael Scheer
*-- Author :
C**********************************************************************
      SUBROUTINE LATTICE
C*********************************************************************
      PARAMETER (LINDIM=100)
      character XLABFR(5),YLABFR(5),XLABAX(5)
      REAL*4 XLABF(1),YLABF(1),XLABA(1)
      EQUIVALENCE (XLABFR(1),XLABF(1)),(YLABFR(1),YLABF(1)),
     &  (XLABAX(1),XLABA(1))
      DIMENSION XLINE(LINDIM),YLINE(LINDIM)

*KEEP,aperseq.
      integer modeli,idev,modus,modeb,modeq,nangle

      real HEIGLI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY,
     &  blimit,wlspos,phimin,phimax

      COMMON /MAGNET/ HEIGLI,MODELI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY

      CHARACTER*80 DFILE

      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
     &  ,PHIMIN(10),PHIMAX(10),NANGLE
*KEND.

      DATA IBESSY/-2/   !0 KEINE MAGNETE
      !1  BESSY I
      !2  BESSY II DOWN STREAM VOM WLS
      !-2 BESSY II UP AND DOWN STREAM VOM WLS

      HEIGLI=0.
      MODELI=2

C--- MAGNETE   DRAMAG(X0,Y0,XLAENGE,YLAENGE,APERTUR,ABLENKWINKEL)

      IF (IBESSY.EQ.2.OR.IBESSY.EQ.-2) THEN

C        call DRAMAG(290.0,0.,31.,50.15,7.2,0.) ! QUADRUPOL TYP B
C        call DRAMAG(333.0,0.,18.,54.5,7.,0.)  ! SEXTUPOL  TYP B
C        call DRAMAG(363.0,0.,49.,50.15,7.2,0.) ! QUADRUPOL TYP A
C        call DRAMAG(424.0,0.,18.,54.5,7.,0.)  ! SEXTUPOL  TYP B
C        call DRAMAG(454.0,0.,31.,50.15,7.2,0.) ! QUADRUPOL TYP B
C        call DIPOL (513.,0.0,50.,20.0,420.,90.,-12.)  ! DIPOL

        call DIPOL (475.394,0.0,50.,20.0,435.897,90.,-11.25)  ! DIPOL
        call DIPOL (931.166,-82.120,50.,20.0,435.897,78.75,-11.25)  ! DIPOL
        call DIPOL (1891.343,-470.773,50.,20.0,435.897,67.5,-11.25)  ! DIPOL
        call DIPOL (2281.000,-721.059,50.,20.0,435.897,56.25,-11.25)  ! DIPOL

        XLINE(1)=0.0
        XLINE(2)=475.394
        YLINE(1)=0.0
        YLINE(2)=0.0
        call mpl(2,XLINE,YLINE)
        call mshplt_arc( 475.394,-435.897,435.897,435.897,78.75,90.00)
        XLINE(1)=560.429
        XLINE(2)=931.166
        YLINE(1)=-8.376
        YLINE(2)=-82.12
        call mpl(2,XLINE,YLINE)
        call mshplt_arc( 846.127,- 509.641 ,435.897,435.897,67.50,78.75)
        XLINE(1)=1012.937
        XLINE(2)=1891.343
        YLINE(1)=-106.925
        YLINE(2)=-470.773
        call mpl(2,xline,yline)
        call mshplt_arc(1724.533,- 873.489 ,435.897,435.897,56.25,67.50)
        XLINE(1)=1966.704
        XLINE(2)=2281.0
        YLINE(1)=-511.054
        YLINE(2)=-721.059
        call mpl(2,xline,yline)
        call mshplt_arc(2038.828,-1083.494 ,435.897,435.897,45.00,56.25)

        IF (IBESSY.EQ.-2) THEN
C        call DRAMAG(-290.0,0.,-31.,50.15,7.2,0.) ! QUADRUPOL TYP B
C        call DRAMAG(-333.0,0.,-18.,54.5,7.,0.)  ! SEXTUPOL  TYP B
C        call DRAMAG(-363.0,0.,-49.,50.15,7.2,0.) ! QUADRUPOL TYP A
C        call DRAMAG(-424.0,0.,-18.,54.5,7.,0.)  ! SEXTUPOL  TYP B
C        call DRAMAG(-454.0,0.,-31.,50.15,7.2,0.) ! QUADRUPOL TYP B
          call DIPOL (-513.,0.0,50.,20.0,420.,90.,12.)  ! DIPOL
        ENDIF

      ELSE IF (IBESSY.EQ.1) THEN

C        call DRAMAG(-100.,0.,84.4,35.,8.0,0.) ! CAVITY
C        call DRAMAG(119.6,0.,32.4,35.,8.0,0.) ! STROMMONITOR ?
C        call DRAMAG(0.,0.,109.6,35.,8.0,0.)   ! WLS
C        call DRAMAG(180.0,0.,42.0,35.,8.0,0.) ! QUADRUPOL
C        call DRAMAG(250.0,0.,42.0,35.,8.0,0.) ! QUADRUPOL
C        call DRAMAG(359.0,0.,123.,35.,8.0,30.)! DIPOL
C        call DRAMAG(-359.0,0.,-123.,35.,8.0,-30.)! DIPOL

        call DIPOL( 332.06,0.0,50.,8.0,177.9,90.,-30.)  ! DIPOL

        call mshplt_arc(0.,0.,25.,42.5,0.,360.00)
        XLINE(1)=-332.06
        XLINE(2)= 332.06
        YLINE(1)=0.0
        YLINE(2)=0.0
        call mpl(2,xline,yline)
        call mshplt_arc(332.06,-177.9,177.9,177.9,60.,90.00)
        XLINE(1)=421.01
        XLINE(2)=594.215
        YLINE(1)=-23.834
        YLINE(2)=-123.834
        call mpl(2,xline,yline)

      ENDIF

      call muwk(0,0)

      RETURN
      END
*CMZ :          26/01/2025  13.25.00  by  Michael Scheer
*CMZ :  1.00/01 23/01/2025  16.31.42  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.14  by  Michael Scheer
*-- Author :
C**************************************************************************
      SUBROUTINE DRAMAG(X0,Y0,XLEN,YLEN,AP,WIN)
C**************************************************************************

C     SR ZEICHNET MAGNET

      PARAMETER (LINDIM=100)
      LOGICAL*1 XLABFR(5),YLABFR(5),XLABAX(5),YLABAX(5)
      REAL*4 XLABF(1),YLABF(1),XLABA(1),YLABA(1)
      EQUIVALENCE (XLABFR(1),XLABF(1)),(YLABFR(1),YLABF(1)),
     &              (XLABAX(1),XLABA(1))
*KEEP,aperseq.
      integer modeli,idev,modus,modeb,modeq,nangle

      real HEIGLI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY,
     &  blimit,wlspos,phimin,phimax

      COMMON /MAGNET/ HEIGLI,MODELI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY

      CHARACTER*80 DFILE

      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
     &  ,PHIMIN(10),PHIMAX(10),NANGLE
*KEND.
      DIMENSION XLINE(LINDIM),YLINE(LINDIM)
c      COMMON /MAGNET/ HEIGLI,MODELI,FACTX,FACTY,
c     &                  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY
c      CHARACTER*80 DFILE
c      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ

      WINB=WIN*1.74533E-2    ! BOGENMASS
      COSW=COS(WINB)
      SINW=SIN(WINB)

      NPOI=7

      XLINE(1)=X0 !X0 IST KOORDINATE DER VORDEREN STIRNFLAECHE
      XLINE(2)=X0+XLEN/2.
      XLINE(3)=XLINE(2)+XLEN/2.*COSW
      XLINE(4)=XLINE(3)
      XLINE(5)=XLINE(2)
      XLINE(6)=XLINE(1)
      XLINE(7)=XLINE(1)

      YLINE(1)=Y0-YLEN/2.
      YLINE(2)=YLINE(1)
      YLINE(3)=YLINE(2)-XLEN/2.*SINW
      YLINE(4)=YLINE(3)+YLEN/2.-AP/2.
      YLINE(5)=YLINE(2)+YLEN/2.-AP/2.
      YLINE(6)=YLINE(1)+YLEN/2.-AP/2.
      YLINE(7)=YLINE(1)

CX    XMINLI=AMIN(XLINE,NPOI)
CX    YMINLI=AMIN(YLINE,NPOI)
      XMINLI=1.D30
      YMINLI=1.D30
      DO III=1,NPOI
         IF (XMINLI.GT.XLINE(III)) XMINLI=XLINE(III)
         IF (YMINLI.GT.YLINE(III)) YMINLI=YLINE(III)
      ENDDO

      X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
      Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

      IF (IDEV.EQ.1.OR.IDEV.EQ.72) THEN
         CALL HPLINE(XLINE,YLINE,NPOI,' ')
         IF (IDEV.EQ.1.OR.IDEV.EQ.72) CALL muwk(0,1)
      ENDIF

      YLINE(1)=YLINE(1)+YLEN/2.+AP/2.
      YLINE(2)=YLINE(2)+YLEN/2.+AP/2.
      YLINE(3)=YLINE(3)+YLEN/2.+AP/2.
      YLINE(4)=YLINE(4)+YLEN/2.+AP/2.
      YLINE(5)=YLINE(5)+YLEN/2.+AP/2.
      YLINE(6)=YLINE(6)+YLEN/2.+AP/2.
      YLINE(7)=YLINE(1)

CX    XMINLI=AMIN(XLINE,NPOI)
CX    YMINLI=AMIN(YLINE,NPOI)
      XMINLI=1.D30
      YMINLI=1.D30
      DO III=1,NPOI
         IF (XMINLI.GT.XLINE(III)) XMINLI=XLINE(III)
         IF (YMINLI.GT.YLINE(III)) YMINLI=YLINE(III)
      ENDDO
      X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
      Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

      IF (IDEV.EQ.1.OR.IDEV.EQ.72) THEN
         CALL HPLINE(XLINE,YLINE,NPOI,' ')
         IF (IDEV.EQ.1) CALL muwk(0,1)
      ENDIF

      RETURN
      END
*CMZ :  1.00/01 23/01/2025  16.31.42  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.14  by  Michael Scheer
*-- Author :
C***********************************************************************
      SUBROUTINE DIPOL(X0,Y0,WIDTH,APER,RADIUS,ANGLE,BENDING)
C***********************************************************************

C DRAWS A CIRCLE

      PARAMETER (LINDIM=100)
      LOGICAL*1 XLABFR(5),YLABFR(5),XLABAX(5)
      REAL*4 XLABF(1),YLABF(1),XLABA(1)
      EQUIVALENCE (XLABFR(1),XLABF(1)),(YLABFR(1),YLABF(1)),
     &              (XLABAX(1),XLABA(1))
      DIMENSION XLINE(LINDIM),YLINE(LINDIM)

*KEEP,aperseq.
      integer modeli,idev,modus,modeb,modeq,nangle

      real HEIGLI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY,
     &  blimit,wlspos,phimin,phimax

      COMMON /MAGNET/ HEIGLI,MODELI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY

      CHARACTER*80 DFILE

      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
     &  ,PHIMIN(10),PHIMAX(10),NANGLE
*KEND.

      COSW=COSD(ANGLE)
      SINW=SIND(ANGLE)

      XCEN=X0-RADIUS*COSW
      YCEN=Y0-RADIUS*SINW

      XLINE(1)=X0+APER/2.*COSW   !X0 IST KOORDINATE DER VORDEREN STIRNFLAECHE
      YLINE(1)=Y0+APER/2.*SINW
      XLINE(2)=X0+WIDTH/2.*COSW
      YLINE(2)=Y0+WIDTH/2.*SINW

      NPOI=50
      MPOI=(NPOI-4)/2
      DO I=1,MPOI-1
          XLINE(I+2)=XCEN+(RADIUS+WIDTH/2.)*COSD(ANGLE+I*BENDING/MPOI)
          YLINE(I+2)=YCEN+(RADIUS+WIDTH/2.)*SIND(ANGLE+I*BENDING/MPOI)
      ENDDO

      XLINE(2+MPOI)=XCEN+(RADIUS+WIDTH/2.)*COSD(ANGLE+BENDING)
      YLINE(2+MPOI)=YCEN+(RADIUS+WIDTH/2.)*SIND(ANGLE+BENDING)
      XLINE(2+MPOI+1)=XCEN+(RADIUS+APER/2.)*COSD(ANGLE+BENDING)
      YLINE(2+MPOI+1)=YCEN+(RADIUS+APER/2.)*SIND(ANGLE+BENDING)

      DO I=1,MPOI-1
          XLINE(I+2+MPOI+1)=XCEN+(RADIUS+APER/2.)
     &                        *COSD(ANGLE+BENDING-I*BENDING/MPOI)
          YLINE(I+2+MPOI+1)=YCEN+(RADIUS+APER/2.)
     &                        *SIND(ANGLE+BENDING-I*BENDING/MPOI)
      ENDDO

      XLINE(2*MPOI+2)=XLINE(1)
      YLINE(2*MPOI+2)=YLINE(1)

      NPOI=2*MPOI+2
CX    XMINLI=AMIN(XLINE,NPOI)
CX    YMINLI=AMIN(YLINE,NPOI)
      XMINLI=1.0e30
      YMINLI=1.0e30
      DO III=1,NPOI
         IF (XMINLI.GT.XLINE(III)) XMINLI=XLINE(III)
         IF (YMINLI.GT.YLINE(III)) YMINLI=YLINE(III)
      ENDDO
      X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
      Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

      CALL mshplt_pline(npoi,XLINE,YLINE)
c      CALL muwk(0,1)

      XLINE(1)=X0-APER/2.*COSW   !X0 IST KOORDINATE DER VORDEREN STIRNFLAECHE
      YLINE(1)=Y0-APER/2.*SINW
      XLINE(2)=X0-WIDTH/2.*COSW
      YLINE(2)=Y0-WIDTH/2.*SINW

      NPOI=50
      MPOI=(NPOI-4)/2
      DO I=1,MPOI-1
          XLINE(I+2)=XCEN+(RADIUS-WIDTH/2.)*COSD(ANGLE+I*BENDING/MPOI)
          YLINE(I+2)=YCEN+(RADIUS-WIDTH/2.)*SIND(ANGLE+I*BENDING/MPOI)
      ENDDO

      XLINE(2+MPOI)=XCEN+(RADIUS-WIDTH/2.)*COSD(ANGLE+BENDING)
      YLINE(2+MPOI)=YCEN+(RADIUS-WIDTH/2.)*SIND(ANGLE+BENDING)
      XLINE(2+MPOI+1)=XCEN+(RADIUS-APER/2.)*COSD(ANGLE+BENDING)
      YLINE(2+MPOI+1)=YCEN+(RADIUS-APER/2.)*SIND(ANGLE+BENDING)

      DO I=1,MPOI-1
          XLINE(I+2+MPOI+1)=XCEN+(RADIUS-APER/2.)
     &                        *COSD(ANGLE+BENDING-I*BENDING/MPOI)
          YLINE(I+2+MPOI+1)=YCEN+(RADIUS-APER/2.)
     &                        *SIND(ANGLE+BENDING-I*BENDING/MPOI)
      ENDDO

      XLINE(2*MPOI+2)=XLINE(1)
      YLINE(2*MPOI+2)=YLINE(1)

      NPOI=2*MPOI+2
CX    XMINLI=AMIN(XLINE,NPOI)
CX    YMINLI=AMIN(YLINE,NPOI)
      XMINLI=1.0e30
      YMINLI=1.0e30
      DO III=1,NPOI
         IF (XMINLI.GT.XLINE(III)) XMINLI=XLINE(III)
         IF (YMINLI.GT.YLINE(III)) YMINLI=YLINE(III)
      ENDDO
      X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
      Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

      CALL mshplt_pline(npoi,XLINE,YLINE)
      CALL muwk(0,1)

      RETURN
      END
*CMZ :          26/01/2025  19.17.01  by  Michael Scheer
*CMZ :  1.00/01 24/01/2025  12.52.37  by  Michael Scheer
*CMZ :  1.00/00 14/01/98  14.34.50  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.14  by  Michael Scheer
*-- Author :
C***********************************************************************
      SUBROUTINE LIGHT(XMINT,XMAXT)
C***********************************************************************

C     LIEST DATENFILE MIT TRAJEKTORIE UND B-FELD UND ZEICHNET LICHTKEGEL

*KEEP,SOURCE.

      PARAMETER (NDIMS=10)
      DIMENSION XISOUR(NDIMS),XESOUR(NDIMS),ICOLS(NDIMS)

      COMMON/SOURC/XISOUR,XESOUR,NSOURCE,ICOLS
*KEND.

      PARAMETER (LINDIM=100)
      LOGICAL*1 XLABFR(5),YLABFR(5)
      REAL*4 XLABF(1),YLABF(1)
      EQUIVALENCE (XLABFR(1),XLABF(1)),(YLABFR(1),YLABF(1))
      DIMENSION IISOUR(NDIMS),IESOUR(NDIMS)
      DIMENSION IC(NDIMS)
      DIMENSION XLINE(LINDIM),YLINE(LINDIM)
*KEEP,aperseq.
      integer modeli,idev,modus,modeb,modeq,nangle

      real HEIGLI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY,
     &  blimit,wlspos,phimin,phimax

      COMMON /MAGNET/ HEIGLI,MODELI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY

      CHARACTER*80 DFILE

      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
     &  ,PHIMIN(10),PHIMAX(10),NANGLE
*KEND.
c      COMMON /MAGNET/ HEIGLI,MODELI,FACTX,FACTY,
c     &                  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY

      REAL*8 DET0,DET1,DET2,DET3,A2,A1,A0,DX0,DXL,DXH,DYL,
     &      DY0,DYH,DXB0,DB0,DBL,DBH,DXL2,DX02,DXH2
      PARAMETER (MAXPOI=100000)
      REAL*4 X(MAXPOI),Y(MAXPOI),B(MAXPOI)

      COMMON/LICHT/X,Y,B,IANZ

c      CHARACTER*80 DFILE
c      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
c     &               ,PHIMIN(10),PHIMAX(10),NANGLE

C--- WLS ABMESSUNGEN

      DATA XLWLS/260./
      DATA YLWLS/50./

      RCOLOUR=1.
      CALL mshplt_igset('PLCI',RCOLOUR)

      ICOLOUR=1
      DO ISOUR=1,NDIMS
          ICOLOUR=ICOLOUR+1
          IF (ICOLOUR.GT.7) ICOLOUR=1
          ICOLS(ISOUR)=ICOLOUR
      ENDDO

C--- ABSORBER

      OPEN(UNIT=47,FILE='aper.in',STATUS='OLD',FORM='FORMATTED')
      !READ(47,*) iabsor
      iabsor=0
      do while (.true.)
        READ(47,*,end=9,iostat=istat) IDIST,APOSX,APOSY,ALENDUM,AWINDUM,ITRANS
        if (istat.ne.0) cycle
        iabsor=iabsor+1
        IF (IDIST.LT.0) THEN
          ALEN=SQRT((APOSX-ALENDUM)**2+(APOSY-AWINDUM)**2)
          IF (ABS(ALENDUM-APOSX).GT.1.E-10)  THEN
            AWIN=ATAND((AWINDUM-APOSY)/(ALENDUM-APOSX))
          ELSE
            IF (AWINDUM.GT.APOSY) THEN
              AWIN=90.
            ELSE
              AWIN=-90.
            ENDIF
          ENDIF
        ELSE
          ALEN=ALENDUM
          AWIN=AWINDUM
        ENDIF

        IF (ITRANS.GE.0) THEN
          CALL ABSORBIN(APOSX,APOSY,ALEN,AWIN)
          CALL muwk(0,1)
        ELSE
          CALL mshplt_set_line_style(2)
          NPOI=2
          XLINE(1)=APOSX
          XLINE(2)=XLINE(1)+COSD(AWIN)*ALEN
          YLINE(1)=APOSY
          YLINE(2)=YLINE(1)+SIND(AWIN)*ALEN
          CALL mshplt_pline(2,XLINE,YLINE)
          CALL muwk(0,1)
          CALL mshplt_set_line_style(1)
        ENDIF

      ENDDO
    9 CLOSE(47)

C--- WLS ZEICHNEN

c     CALL DRAMAG(WLSPOS-125.,0.,90.,50.,0.,0.) ! RANDPOL
c     CALL DRAMAG(WLSPOS- 15.,0.,30.,50.,0.,0.) ! HAUPTPOL
c     CALL DRAMAG(WLSPOS+ 35.,0.,90.,50.,0.,0.) ! RANDPOL

C--- TRAJEKTORIE LESEN

      CALL README(ICODE)

C--- TRAJEKTORIE ZEICHNEN

CX       XMINLI=AMIN(X,IANZ)
CX       YMINLI=AMIN(Y,IANZ)
      XMINLI=1.e30
      YMINLI=1.e30

      DO III=1,IANZ
         IF (XMINLI.GT.X(III)) XMINLI=X(III)
         IF (YMINLI.GT.Y(III)) YMINLI=Y(III)
      ENDDO

      X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
      Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

      call mshplt_pline(ianz,x,y)
      CALL muwk(0,1)

C--- STRAHLEN ZEICHNEN

      NSOURCE=1
      XISOUR(1)=-1.E30
      XESOUR(1)=1.E30

      IF (MODUS.EQ.0) THEN

C QUELLPUNKTE LESEN

        CALL READL0(ICODEL,NSOURCE,XISOUR,XESOUR,IISOUR,IESOUR,NDIMS,GAMMA1)

        if (nsource.eq.0) then
          nsource=1
          xisour(1)=xminli
          xesour(1)=xmaxli
        endif

        IF (ICODEL.NE.ICODE) THEN
          WRITE(6,*) '*** APER: FILES FOR TRAJEKTORY AND SOURCE POINTS NOT CONSISTENT'
          STOP
        ENDIF !ICODEL

        BLIMIT=0.0
        MODUS=2
        MODEB=0

      ENDIF !MODUS

      IF (MODUS.EQ.1)  THEN

C--- LICHTKEGEL VOM ORT DES MAX. FELDES

        BMAX1=-1.E30
        BMAX2=-1.E30
        DO I=1,IANZ
          IF (ABS(B(I)).GT.BMAX1) THEN
            BMAX1=ABS(B(I))
            IBMAX1=I
          ENDIF
          IF (ABS(B(I)).GE.BMAX2) THEN
            BMAX2=ABS(B(I))
            IBMAX2=I
          ENDIF
        END DO

        IBMAX=(IBMAX1+IBMAX2)/2
        BMAX=ABS(B(IBMAX))
        IF (BMAX.NE.BMAX1) STOP '*** BMAX NICHT GEFUNDEN ***'

C---BFELD INTERPOLIEREN (PARABEL DURCH DREI PUNKTE)

        X0=X(IBMAX)
        Y0=Y(IBMAX)
        B0=B(IBMAX)

        XL=X(IBMAX-1)
        XH=X(IBMAX+1)
        YL=Y(IBMAX-1)
        YH=Y(IBMAX+1)
        BL=B(IBMAX-1)
        BH=B(IBMAX+1)

C SKIP PARABEL IF IN CONSTANT FIELD REGION

        IF (ABS(BMAX-ABS(BL)).GT.0.001*BMAX
     &      .AND.ABS(BMAX-ABS(BH)).GT.0.001*BMAX) THEN

          DXL=XL
          DX0=X0
          DXH=XH
          DYL=YL
          DY0=Y0
          DYH=YH
          DBL=BL
          DB0=B0
          DBH=BH

          DXL2=DXL*DXL
          DX02=DX0*DX0
          DXH2=DXH*DXH

          DET0=(DX02*DXH-DXH2*DX0)-(DXL2*DXH-DXH2*DXL)+(DXL2*DX0-DX02*DXL)
          DET1=(DB0 *DXH-DBH *DX0)-(DBL *DXH-DBH *DXL)+(DBL *DX0-DB0 *DXL)
          DET2=(DX02*DBH-DXH2*DB0)-(DXL2*DBH-DXH2*DBL)+(DXL2*DB0-DX02*DBL)
          DET3=DBL*(DX02*DXH-DXH2*DX0)-DB0*(DXL2*DXH-DXH2*DXL)
     &      +DBH*(DXL2*DX0-DX02*DXL)

          A2=DET1/DET0
          A1=DET2/DET0
          A0=DET3/DET0

C--- INTERPOLATIONSPOLYNOM LAUTET JETZT B(X)=A2*X*X+A1*X+A0

          DXB0=X0
          IF(A2.NE.0.0)DXB0=-A1/(2.D0*A2)   !MAXIMUM DES POLYNOMS NEU BERECHEN

C---KURVE INTERPOLIEREN (PARABEL DURCH DREI PUNKTE)

          DET1=(DY0 *DXH-DYH *DX0)-(DYL *DXH-DYH *DXL)+(DYL *DX0-DY0 *DXL)
          DET2=(DX02*DYH-DXH2*DY0)-(DXL2*DYH-DXH2*DYL)+(DXL2*DY0-DX02*DYL)
          DET3=DYL*(DX02*DXH-DXH2*DX0)-DY0*(DXL2*DXH-DXH2*DXL)
     &      +DYH*(DXL2*DX0-DX02*DXL)

          A2=DET1/DET0
          A1=DET2/DET0
          A0=DET3/DET0

C--- INTERPOLATIONSPOLYNOM LAUTET JETZT Y(X)=A2*X*X+A1*X+A0

          Y0=sngl(A2*DXB0*DXB0+A1*DXB0+A0)     ! Y0 NEU BERECHEN
          SLOPE0=sngl(2.0*A2*DXB0+A1)

        ELSE   !SKIP PLATEAU

          SLOPE0=(YH-YL)/(XH-XL)

        ENDIF !SKIP PLATEAU


        XE=SMAAXX
        YE=Y0+SLOPE0*(XE-X0)
        CALL ABSORB (X0,Y0,XE,YE) !CHECK ABSORBER
        IF (YE.GT.SMAAXY) THEN
          YE=SMAAXY
          XE=(YE-Y0)/SLOPE0+X0
        ELSE IF (YE.LT.SMIAXY) THEN
          YE=SMIAXY
          XE=(YE-Y0)/SLOPE0+X0
        ENDIF
        NPOI=2
        XLINE(1)=X0
        YLINE(1)=Y0
        XLINE(2)=XE
        YLINE(2)=YE
CX    XMINLI=AMIN(XLINE,NPOI)
CX    YMINLI=AMIN(YLINE,NPOI)
        XMINLI=1.e30
        YMINLI=1.e30
        DO III=1,NPOI
          IF (XMINLI.GT.XLINE(III)) XMINLI=XLINE(III)
          IF (YMINLI.GT.YLINE(III)) YMINLI=YLINE(III)
        ENDDO
        X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
        Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

        IF (IDEV.EQ.1.OR.IDEV.EQ.72) THEN
          IF (X0.GE.XMINT.AND.X0.LE.XMAXT)
     &      CALL MSPLINE(XLINE,YLINE,NPOI,' ')
          IF (IDEV.EQ.1) CALL muwk(0,1)
        ENDIF
C     SLOPE1=TAN(ATAN(SLOPE0)-ATAN(GAMMA1))
        SLOPE1=TAN(ATAN(SLOPE0)-0.003)
        XE=SMAAXX
        YE=Y0+SLOPE1*(XE-X0)
        CALL ABSORB (X0,Y0,XE,YE) !CHECK ABSORBER
        IF (YE.GT.SMAAXY) THEN
          YE=SMAAXY
          XE=(YE-Y0)/SLOPE1+X0
        ELSE IF (YE.LT.SMIAXY) THEN
          YE=SMIAXY
          XE=(YE-Y0)/SLOPE1+X0
        ENDIF
        NPOI=2
        XLINE(1)=X0
        YLINE(1)=Y0
        XLINE(2)=XE
        YLINE(2)=YE
CX    XMINLI=AMIN(XLINE,NPOI)
CX    YMINLI=AMIN(YLINE,NPOI)
        XMINLI=1.e30
        YMINLI=1.e30
        DO III=1,NPOI
          IF (XMINLI.GT.XLINE(III)) XMINLI=XLINE(III)
          IF (YMINLI.GT.YLINE(III)) YMINLI=YLINE(III)
        ENDDO
        X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
        Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

        IF (IDEV.EQ.1.OR.IDEV.EQ.72) THEN
          IF (X0.GE.XMINT.AND.X0.LE.XMAXT)
     &      CALL MSPLINE(XLINE,YLINE,NPOI,' ')
          IF (IDEV.EQ.1) CALL muwk(0,1)
        ENDIF
C     SLOPE2=TAN(ATAN(SLOPE0)+ATAN(GAMMA1))
        SLOPE2=TAN(ATAN(SLOPE0)+0.003)
        XE=SMAAXX
        YE=Y0+SLOPE2*(XE-X0)
        CALL ABSORB (X0,Y0,XE,YE) !CHECK ABSORBER
        IF (YE.GT.SMAAXY) THEN
          YE=SMAAXY
          XE=(YE-Y0)/SLOPE2+X0
        ELSE IF (YE.LT.SMIAXY) THEN
          YE=SMIAXY
          XE=(YE-Y0)/SLOPE2+X0
        ENDIF
        NPOI=2
        XLINE(1)=X0
        YLINE(1)=Y0
        XLINE(2)=XE
        YLINE(2)=YE
CX    XMINLI=AMIN(XLINE,NPOI)
CX    YMINLI=AMIN(YLINE,NPOI)
        XMINLI=1.e30
        YMINLI=1.e30
        DO III=1,NPOI
          IF (XMINLI.GT.XLINE(III)) XMINLI=XLINE(III)
          IF (YMINLI.GT.YLINE(III)) YMINLI=YLINE(III)
        ENDDO
        X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
        Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

        IF (IDEV.EQ.1.OR.IDEV.EQ.72) THEN
          IF (X0.GE.XMINT.AND.X0.LE.XMAXT)
     &      CALL MSPLINE(XLINE,YLINE,NPOI,' ')
          IF (IDEV.EQ.1) CALL muwk(0,1)
        ENDIF

      ELSE IF (MODUS.EQ.2) THEN

        BMAX=-1.E30
        BLOW= 1.E30
        DO I=1,IANZ
          IF (B(I).GT.BMAX) THEN
            BMAX=B(I)
            IBMAX=I
          ENDIF
          IF (B(I).LT.BLOW) THEN
            BLOW=B(I)
            IBLOW=I
          ENDIF
        END DO

C---------------------------------------------------
        IF (MODEB.EQ.0) THEN
          BLEVEL=AMAX1(ABS(BLOW),ABS(BMAX))*ABS(BLIMIT)
        ELSE IF (MODEB.EQ.1) THEN
          BLEVEL=BMAX*ABS(BLIMIT)
        ELSE IF (MODEB.EQ.2) THEN
          BLEVEL=-BLOW*ABS(BLIMIT)
        ENDIF

C FIND MIDDLE OF SOURCE (DEFINED BY ANGLE)

        IF (MODEQ.EQ.1) THEN
          DO ISOUR=1,NSOURCE
            SLOPEI= (Y(IISOUR(ISOUR))-Y(IISOUR(ISOUR)-1))/
     &        (X(IISOUR(ISOUR))-X(IISOUR(ISOUR)-1))
            SLOPEE= (Y(IESOUR(ISOUR))-Y(IESOUR(ISOUR)-1))/
     &        (X(IESOUR(ISOUR))-X(IESOUR(ISOUR)-1))
            SLOPEM=(SLOPEI+SLOPEE)/2.
            DIFFM=1.E30
            DO J=IISOUR(ISOUR)+1,IESOUR(ISOUR)-1
              SLOPE=(Y(J+1)-Y(J-1))/(X(J+1)-X(J-1))
              IF (ABS(SLOPE-SLOPEM).LT.DIFFM) THEN
                DIFFM=ABS(SLOPE-SLOPEM)
                IC(ISOUR)=J
              ENDIF   !SLOPE
            ENDDO !J
          ENDDO !ISOUR
        ENDIF !MODEQ

        DO 500 I=2,IANZ-1

          IF (MODEB.EQ.0) THEN
            BCUR=ABS(B(I))
          ELSE IF (MODEB.EQ.1) THEN
            BCUR=B(I)
          ELSE IF (MODEB.EQ.2) THEN
            BCUR=-B(I)
          ENDIF

          X0=X(I)
          Y0=Y(I)
          SLOPE=(Y(I+1)-Y(I-1))/(X(I+1)-X(I-1))
          IFLAGS=0
          IF(MODEQ.EQ.0) THEN
            JFLAGS=1
          ELSE    !MODEQ
            JFLAGS=0
          ENDIF   !MODEQ

          DO ISOUR=1,NSOURCE
            IF (X0.GE.XISOUR(ISOUR).AND.X0.LE.XESOUR(ISOUR)) THEN
              IFLAGS=1
            ENDIF
            IF (IFLAGS.EQ.1.AND.MODEQ.EQ.1) THEN
              IF (I.EQ.IC(ISOUR)) JFLAGS=1
            ENDIF !MODEQ
          ENDDO !ISOUR

          IF (BCUR.GE.BLEVEL.AND.IFLAGS.EQ.1.AND.JFLAGS.EQ.1) THEN
            XE=SMAAXX
            YE=Y0+SLOPE*(XE-X0)
            CALL ABSORB (X0,Y0,XE,YE) !CHECK ABSORBER
            IF (YE.GT.SMAAXY) THEN
              YE=SMAAXY
              XE=(YE-Y0)/SLOPE+X0
            ELSE IF (YE.LT.SMIAXY) THEN
              YE=SMIAXY
              XE=(YE-Y0)/SLOPE+X0
            ENDIF
            NPOI=2
            XLINE(1)=X0
            YLINE(1)=Y0
            XLINE(2)=XE
            YLINE(2)=YE
CX       XMINLI=AMIN(XLINE,NPOI)
CX       YMINLI=AMIN(YLINE,NPOI)
            XMINLI=1.e30
            YMINLI=1.e30
            DO III=1,NPOI
              IF (XMINLI.GT.XLINE(III)) XMINLI=XLINE(III)
              IF (YMINLI.GT.YLINE(III)) YMINLI=YLINE(III)
            ENDDO
            X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
            Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

            IF (IDEV.EQ.1.OR.IDEV.EQ.72) THEN
              IF (X0.GE.XMINT.AND.X0.LE.XMAXT)
     &          CALL MSPLINE(XLINE,YLINE,NPOI,' ')
              IF (IDEV.EQ.1) CALL muwk(0,1)
            ENDIF
            IF (MODEQ.EQ.1) THEN

              SLOPEC=SLOPE
              SLOPE=(SLOPEC+GAMMA1)/(1.-SLOPEC*GAMMA1)
              XE=SMAAXX
              YE=Y0+SLOPE*(XE-X0)
              CALL ABSORB (X0,Y0,XE,YE) !CHECK ABSORBER
              IF (YE.GT.SMAAXY) THEN
                YE=SMAAXY
                XE=(YE-Y0)/SLOPE+X0
              ELSE IF (YE.LT.SMIAXY) THEN
                YE=SMIAXY
                XE=(YE-Y0)/SLOPE+X0
              ENDIF
              NPOI=2
              XLINE(1)=X0
              YLINE(1)=Y0
              XLINE(2)=XE
              YLINE(2)=YE
CX          XMINLI=AMIN(XLINE,NPOI)
CX          YMINLI=AMIN(YLINE,NPOI)
              XMINLI=1.e30
              YMINLI=1.e30
              DO III=1,NPOI
                IF (XMINLI.GT.XLINE(III)) XMINLI=XLINE(III)
                IF (YMINLI.GT.YLINE(III)) YMINLI=YLINE(III)
              ENDDO
              X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
              Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

              IF (IDEV.EQ.1.OR.IDEV.EQ.72) THEN
                IF (X0.GE.XMINT.AND.X0.LE.XMAXT)
     &            CALL MSPLINE(XLINE,YLINE,NPOI,' ')
                IF (IDEV.EQ.1) CALL muwk(0,1)
              ENDIF
              SLOPE=(SLOPEC-GAMMA1)/(1.+SLOPEC*GAMMA1)
              XE=SMAAXX
              YE=Y0+SLOPE*(XE-X0)
              CALL ABSORB (X0,Y0,XE,YE) !CHECK ABSORBER
              IF (YE.GT.SMAAXY) THEN
                YE=SMAAXY
                XE=(YE-Y0)/SLOPE+X0
              ELSE IF (YE.LT.SMIAXY) THEN
                YE=SMIAXY
                XE=(YE-Y0)/SLOPE+X0
              ENDIF
              NPOI=2
              XLINE(1)=X0
              YLINE(1)=Y0
              XLINE(2)=XE
              YLINE(2)=YE
CX          XMINLI=AMIN(XLINE,NPOI)
CX          YMINLI=AMIN(YLINE,NPOI)
              XMINLI=1.e30
              YMINLI=1.e30
              DO III=1,NPOI
                IF (XMINLI.GT.XLINE(III)) XMINLI=XLINE(III)
                IF (YMINLI.GT.YLINE(III)) YMINLI=YLINE(III)
              ENDDO
              X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
              Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

              IF (IDEV.EQ.1.OR.IDEV.EQ.72) THEN
                IF (X0.GE.XMINT.AND.X0.LE.XMAXT)
     &            CALL MSPLINE(XLINE,YLINE,NPOI,' ')
                IF (IDEV.EQ.1) CALL muwk(0,1)
              ENDIF
            ENDIF !MODEQ

          ENDIF  !BCUR
C---------------------------------------------------
500     CONTINUE

      ELSEIF (MODUS.EQ.3) THEN   ! ZEICHNE STRAHL, FALLS IN WINKELBEREICH

        DO I=2,IANZ-1
          X0=X(I)
          Y0=Y(I)
          SLOPE=(Y(I+1)-Y(I-1))/(X(I+1)-X(I-1))

          DO IPHI=1,NANGLE
            IF (     SLOPE.GE.TAND(PHIMIN(IPHI))
     &          .AND.SLOPE.LE.TAND(PHIMAX(IPHI))) THEN
              XE=SMAAXX
              YE=Y0+SLOPE*(XE-X0)
              CALL ABSORB (X0,Y0,XE,YE) !CHECK ABSORBER
              IF (YE.GT.SMAAXY) THEN
                YE=SMAAXY
                XE=(YE-Y0)/SLOPE+X0
              ELSE IF (YE.LT.SMIAXY) THEN
                YE=SMIAXY
                XE=(YE-Y0)/SLOPE+X0
              ENDIF
              NPOI=2
              XLINE(1)=X0
              XLINE(2)=XE
              YLINE(1)=Y0
              YLINE(2)=YE
              IF (IDEV.EQ.1.OR.IDEV.EQ.72) THEN
                IF (X0.GE.XMINT.AND.X0.LE.XMAXT)
     &            CALL MSPLINE(XLINE,YLINE,NPOI,' ')
                IF (IDEV.EQ.1) CALL muwk(0,1)
              ENDIF
            ENDIF
          ENDDO

        ENDDO

      ENDIF   !MODUS

      CALL mshplt_igset('PLCI',1.)

      RETURN
      END
*CMZ :          26/01/2025  17.44.03  by  Michael Scheer
*CMZ :  1.00/01 23/01/2025  16.31.42  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.14  by  Michael Scheer
*-- Author :
C************************************************************************
      SUBROUTINE README(ICODE)
C*************************************************************************
      PARAMETER (MAXPOI=100000)
      REAL*4 X(MAXPOI),Y(MAXPOI),B(MAXPOI)

      COMMON/LICHT/X,Y,B,IANZ

*KEEP,aperseq.
      integer modeli,idev,modus,modeb,modeq,nangle

      real HEIGLI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY,
     &  blimit,wlspos,phimin,phimax

      COMMON /MAGNET/ HEIGLI,MODELI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY

      CHARACTER*80 DFILE

      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
     &  ,PHIMIN(10),PHIMAX(10),NANGLE
*KEND.
c      CHARACTER*80 DFILE,COMMENT1
c      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
c     &               ,PHIMIN(10),PHIMAX(10),NANGLE

CW    IF(DFILE.NE.'TRACK.DAT') THEN
CW       OPEN(UNIT=46,FILE=DFILE,STATUS='OLD',FORM='UNFORMATTED')
CW       READ(46)IANZ
CW       IF (IANZ.GT.MAXPOI) STOP '*** SR READ IANZ.GT.MAXPOI ***'
CW       DO I=1,IANZ
CW          READ(46) X(I),DUMMY,Y(I),B(I)
CW       ENDDO
CW    ELSE
CW       OPEN(UNIT=46,FILE=DFILE,STATUS='OLD',FORM='FORMATTED')
CW
CW          I=0
CW1357         I=I+1
CW          IF (IANZ+1.GT.MAXPOI) STOP '*** SR READ IANZ.GT.MAXPOI ***'
CW          READ(46,*,END=579) X(I),DUMMY,Y(I),DUMMY,DUMMY,DUMMY,
CW     &                                DUMMY,B(I),DUMMY,DUMMY,DUMMY,DUMMY
CW          IANZ=IANZ+1
CW          GOTO 1357
CW579       CONTINUE
CW    ENDIF
CW
CW    CLOSE(46)

      OPEN(UNIT=46,FILE=DFILE,STATUS='OLD',FORM='FORMATTED',iostat=istat)

      READ(46,*,iostat=istat) ICODE
      READ(46,'(A80)',iostat=istat) COMMENT1
      do k=1,maxpoi
        READ(46,*,END=579) X(K),D,Y(K),D,D,D,D,B(K)
      enddo

579   IANZ=K-1

      CLOSE(46)

      DO I=1,IANZ
         X(I)=X(I)*100.+WLSPOS
c        Y(I)=-Y(I)*100.
         Y(I)=+Y(I)*100.
      ENDDO

      RETURN
      END
*CMZ :          26/01/2025  16.20.04  by  Michael Scheer
*CMZ :  1.00/01 23/01/2025  16.31.42  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.14  by  Michael Scheer
*-- Author :
C**********************************************************************
      SUBROUTINE ABSORBIN(XA0,YA0,ALENA,WINA)
C**********************************************************************

C INITIALIESIERT ABSORBER FUER SR ABSORB
CXC XA0,YA0,ALEN,WINA = ZENTRUM, LAENGE UND NEIGUNG (GRAD) DES ABSORBERS
C XA0,YA0,ALEN,WINA = ANFANG, LAENGE UND NEIGUNG (GRAD) DES ABSORBERS

C     IMPLICIT NONE
      INTEGER MABSO
      PARAMETER (MABSO=100)
      INTEGER ICALL,IHIT(MABSO),NABSO
      REAL*4 XA0,YA0,ALENA,WINA
      REAL*4 X0ABSO(MABSO),Y0ABSO(MABSO),XABSO(MABSO,2),YABSO(MABSO,2),
     &          SLOPEA(MABSO),OFFA(MABSO),XHIT(MABSO),YHIT(MABSO),
     &          ALEN(MABSO),WIN(MABSO)

      COMMON/ABSORBER/X0ABSO,Y0ABSO,XABSO,YABSO,SLOPEA,OFFA,XHIT,YHIT,IHIT,
     &  NABSO
*KEEP,aperseq.
      integer modeli,idev,modus,modeb,modeq,nangle

      real HEIGLI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY,
     &  blimit,wlspos,phimin,phimax

      COMMON /MAGNET/ HEIGLI,MODELI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY

      CHARACTER*80 DFILE

      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
     &  ,PHIMIN(10),PHIMAX(10),NANGLE
*KEND.
c      CHARACTER*80 DFILE
c      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
c     &               ,PHIMIN(10),PHIMAX(10),NANGLE

      DATA ICALL/0/

      ICALL=ICALL+1
      NABSO=ICALL

      IF(NABSO.GT.MABSO) STOP '*** ERROR SR ABSORBIN: DIMENSION FALSCH ***'

      X0ABSO(NABSO)=XA0
      Y0ABSO(NABSO)=YA0
      ALEN(NABSO)=ALENA
      WIN(NABSO)=WINA

CX    XABSO(NABSO,1)=XA0-COSD(WINA)*ALENA/2.
CX    XABSO(NABSO,2)=XA0+COSD(WINA)*ALENA/2.
CX    YABSO(NABSO,1)=YA0-SIND(WINA)*ALENA/2.
CX    YABSO(NABSO,2)=YA0+SIND(WINA)*ALENA/2.

      XABSO(NABSO,1)=XA0
      XABSO(NABSO,2)=XA0+COSD(WINA)*ALENA
      YABSO(NABSO,1)=YA0
      YABSO(NABSO,2)=YA0+SIND(WINA)*ALENA

      SLOPEA(NABSO)=9999.
      OFFA(NABSO)=-9999.
      IF(ABS(XABSO(NABSO,1)-XABSO(NABSO,2)).GT.1.E-5) THEN
          SLOPEA(NABSO)=(YABSO(NABSO,2)-YABSO(NABSO,1))/
     &                     (XABSO(NABSO,2)-XABSO(NABSO,1))
          OFFA(NABSO)=YABSO(NABSO,1)-SLOPEA(NABSO)*XABSO(NABSO,1)
      ENDIF

C--- ABSORBER ZEICHEN
      CALL DRABSO(XABSO(NABSO,1),XABSO(NABSO,2),
     &          YABSO(NABSO,1),YABSO(NABSO,2))

      RETURN
      END
*CMZ :  1.00/01 23/01/2025  16.31.42  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.14  by  Michael Scheer
*-- Author :
C**********************************************************************
      SUBROUTINE ABSORB(X0,Y0,XE,YE)
C**********************************************************************

C BERECHNET ENDPUNKT EINES LICHTSTRAHLES UNTER BERUECKSICHTIGUNG VON ABSORBERN
C X0,Y0,XE,YE,SLOPE0,OFF0 = ANFANG,ENDE,STEIGUNG,Y-ABSCHNITT DES LICHTSTRAHLES
C X0ABSO,Y0ABSO,ALEN,WIN = ZENTRUM, LAENGE UND WINKEL(GRAD) DES ABSORBERS
C XABSO,YABSO,SLOPEA,OFFA = ANFANG,ENDE,STEIGUNG,Y-ABSCHNITT DES ABSORBERS
C XHIT,YHIT,IHIT = TREFFERPUNKT UND TREFFER-FLAGGE
C NABSO = ANZAHLE DER ABSORBER

C ROUTINE SETZT VORAUS, DASS LICHTSTRAHL VON LINKS NACH RECHTS GEHT

      IMPLICIT NONE
      INTEGER MABSO
      PARAMETER (MABSO=100)
      INTEGER I,IHIT(MABSO),NABSO
      REAL*4 X0ABSO(MABSO),Y0ABSO(MABSO),XABSO(MABSO,2),YABSO(MABSO,2),
     &          SLOPEA(MABSO),OFFA(MABSO),XHIT(MABSO),YHIT(MABSO),
     &          ALEN(MABSO),WIN(MABSO),SLOPE0,OFF0,X0,Y0,XE,YE

      COMMON/ABSORBER/X0ABSO,Y0ABSO,XABSO,YABSO,SLOPEA,OFFA,XHIT,YHIT,IHIT,
     &              NABSO

      IF (NABSO.GT.MABSO) STOP '*** SR ABSOR: DIMENSION FALSCH ***'

C--- LICHTSTRAHL
      SLOPE0=(YE-Y0)/(XE-X0)
      OFF0=YE-SLOPE0*XE

C--- CHECK ABSORBER
      DO I=1,NABSO
          XHIT(I)=1.E30
          YHIT(I)=1.E30
          IHIT(I)=0
          IF (SLOPEA(I).NE.SLOPE0.AND.SLOPEA(I).NE.9999.) THEN
         XHIT(I)=(OFFA(I)-OFF0)/(SLOPE0-SLOPEA(I))
         YHIT(I)=SLOPE0*XHIT(I)+OFF0
          ELSE IF (SLOPEA(I).NE.SLOPE0.AND.SLOPEA(I).EQ.9999.) THEN
         XHIT(I)=X0ABSO(I)
         YHIT(I)=SLOPE0*XHIT(I)+OFF0
          ELSE IF (SLOPEA(I).EQ.SLOPE0) THEN
          IF (OFFA(I).EQ.OFF0) THEN
             IF (X0.GT.XABSO(I,1).AND.X0.LT.XABSO(I,2)) THEN
            X0=XABSO(I,1)
            XE=XABSO(I,2)
             ELSEIF (X0.LT.XABSO(I,1)
     &                 .AND.XE.GT.XABSO(I,1).AND.XE.LT.XABSO(I,2)) THEN
            XE=XABSO(I,1)
             ENDIF
          ENDIF
          GOTO 90
          END IF

        IF((XHIT(I)-XABSO(I,1))*(XHIT(I)-XABSO(I,2)).LT.1.E-5
     &       .AND.
     &      (YHIT(I)-YABSO(I,1))*(YHIT(I)-YABSO(I,2)).LT.1.E-5)
     &       IHIT(I)=1
          IF (IHIT(I).EQ.1.AND.XHIT(I).GE.X0.AND.XHIT(I).LT.XE) THEN
           XE=XHIT(I)
           YE=YHIT(I)
        END IF
90    CONTINUE
        END DO

      RETURN
      END
*CMZ :          26/01/2025  13.05.07  by  Michael Scheer
*CMZ :  1.00/01 23/01/2025  16.31.42  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.14  by  Michael Scheer
*-- Author :
C**************************************************************************
      SUBROUTINE DRABSO(X1,X2,Y1,Y2)
C**************************************************************************

C     SR ZEICHNET ABSORBER

      PARAMETER (LINDIM=100)
      LOGICAL*1 XLABFR(5),YLABFR(5),XLABAX(5),YLABAX(5)
      REAL*4 XLABF(1),YLABF(1),XLABA(1),YLABA(1)
      EQUIVALENCE (XLABFR(1),XLABF(1)),(YLABFR(1),YLABF(1)),
     &              (XLABAX(1),XLABA(1))
      DIMENSION XLINE(LINDIM),YLINE(LINDIM)
*KEEP,aperseq.
      integer modeli,idev,modus,modeb,modeq,nangle

      real HEIGLI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY,
     &  blimit,wlspos,phimin,phimax

      COMMON /MAGNET/ HEIGLI,MODELI,FACTX,FACTY,
     &  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY

      CHARACTER*80 DFILE

      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
     &  ,PHIMIN(10),PHIMAX(10),NANGLE
*KEND.
c      CHARACTER*80 DFILE
c      COMMON/CONTRL/IDEV,DFILE,MODUS,BLIMIT,WLSPOS,MODEB,MODEQ
c     &               ,PHIMIN(10),PHIMAX(10),NANGLE
c      COMMON /MAGNET/ HEIGLI,MODELI,FACTX,FACTY,
c     &                  XMINFR,YMINFR,SMIAXX,SMAAXX,SMAAXY,SMIAXY


      NPOI=2

      XLINE(1)=X1
      XLINE(2)=X2

      YLINE(1)=Y1
      YLINE(2)=Y2

      IF(X1.GT.X2) THEN
         XLINE(2)=X1
         XLINE(1)=X2
         YLINE(2)=Y1
         YLINE(1)=Y2
      ENDIF

      SLOPE0=1.E15
      IF(X1.NE.X2) SLOPE0=(Y2-Y1)/(X2-X1)
c050893  XLINE(2)=SMAAXX
c050893  YLINE(2)=SLOPE0*(XLINE(2)-XLINE(1))+YLINE(1)

      IF (YLINE(2).GT.SMAAXY.AND.SLOPE0.NE.0.0) THEN
            YLINE(2)=SMAAXY
          XLINE(2)=(YLINE(2)-YLINE(1))/SLOPE0+XLINE(1)
        ELSE IF (YLINE(2).LT.SMIAXY.AND.SLOPE0.NE.0.0) THEN
            YLINE(2)=SMIAXY
          XLINE(2)=(YLINE(2)-YLINE(1))/SLOPE0+XLINE(1)
      ENDIF

CX    XMINLI=AMIN(XLINE,NPOI)
CX    YMINLI=AMIN(YLINE,NPOI)
      XMINLI=1.D30
      YMINLI=1.D30
      DO III=1,NPOI
         IF (XMINLI.GT.XLINE(III)) XMINLI=XLINE(III)
         IF (YMINLI.GT.YLINE(III)) YMINLI=YLINE(III)
      ENDDO
      X0LINE=XMINFR+(XMINLI-SMIAXX)*FACTX
      Y0LINE=YMINFR+(YMINLI-SMIAXY)*FACTY

      IF (IDEV.EQ.1.OR.IDEV.EQ.72) THEN
         CALL HPLINE(XLINE,YLINE,NPOI,' ')
         IF (IDEV.EQ.1.OR.IDEV.EQ.72) CALL muwk(0,1)
      ENDIF
      RETURN
      END
*CMZ :  1.00/01 23/01/2025  16.31.42  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.14  by  Michael Scheer
*-- Author :
C************************************************************
      SUBROUTINE PAWLINE(X,Y,N)
      REAL*4 X(N),Y(N)

C--- SCHREIBT X,Y-KOORDINATEN PUNKTWEISE AUF DATENFILE FUER PAW

      DO I=1,N-1
          WRITE(55,*)X(I)/100.,Y(I)/100.
          WRITE(55,*)X(I+1)/100.,Y(I+1)/100.
      ENDDO
      RETURN
      END
*CMZ :          26/01/2025  15.19.57  by  Michael Scheer
*CMZ :  1.00/01 23/01/2025  16.31.42  by  Michael Scheer
*CMZ :  1.00/00 14/01/98  14.16.45  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.14  by  Michael Scheer
*-- Author :
C************************************************************
      SUBROUTINE READL0
     &  (ICODE,NSOURCE,XISOUR,XESOUR,IISOUR,IESOUR,NDIM,GAMMA1)

      DIMENSION XISOUR(NDIM),XESOUR(NDIM),IISOUR(NDIM),IESOUR(NDIM)

      OPEN(UNIT=49,FILE='wave_l0.dat',STATUS='OLD',FORM='FORMATTED',iostat=istat)

      if (istat.eq.0) then

        READ(49,*) ICODE
        READ(49,*) GAMMA
        GAMMA1=1./GAMMA
        READ(49,*) NSOURCE
        READ(49,*)DUMMY,DUMMY,DUMMY
        READ(49,*)DUMMY,DUMMY,DUMMY,DUMMY,DUMMY
        READ(49,*)DUMMY,DUMMY,DUMMY,DUMMY,DUMMY
        DO ISOUR=1,NSOURCE
          READ(49,*)IISOUR(ISOUR),IESOUR(ISOUR)
          READ(49,*)XISOUR(ISOUR)
     &      ,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY
     &      ,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY
          READ(49,*)XESOUR(ISOUR)
     &      ,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY
     &      ,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY
          XISOUR(ISOUR)=XISOUR(ISOUR)*100.
          XESOUR(ISOUR)=XESOUR(ISOUR)*100.
        ENDDO   !ISOUR

        CLOSE(49)
      else
        nsource=0
      endif

      RETURN
      END
*CMZ :          26/01/2025  16.56.17  by  Michael Scheer
*CMZ :  1.00/01 23/01/2025  16.31.42  by  Michael Scheer
*-- Author :    Michael Scheer   22/01/98

      SUBROUTINE MSPLINE(X,Y,N,CHOPT)

*KEEP,SOURCE.

      PARAMETER (NDIMS=10)
      DIMENSION XISOUR(NDIMS),XESOUR(NDIMS),ICOLS(NDIMS)

      COMMON/SOURC/XISOUR,XESOUR,NSOURCE,ICOLS
*KEND.

      CHARACTER CHOPT
      DIMENSION X(N),Y(N)

      R=1.
      DO ISOUR=1,NSOURCE
          IF (X(1).GE.XISOUR(ISOUR).AND.X(1).LE.XESOUR(ISOUR)) THEN
         R=ICOLS(ISOUR)
         GOTO 90
          ELSE
         R=1.
          ENDIF
      ENDDO

90    CALL mshplt_igset('PLCI',R)
      CALL HPLINE(X,Y,N,chopt)

      CALL mshplt_igset('PLCI',1.)

      call mshplt_flush_buff

      RETURN
      END
*CMZ :  1.00/01 25/01/2025  19.25.27  by  Michael Scheer
*CMZ :  1.00/00 14/01/98  12.20.43  by  Michael Scheer
*CMZ :  0.00/00 13/01/98  18.06.14  by  Michael Scheer
*-- Author :
C**********************************************************************
      SUBROUTINE magseq
C*********************************************************************

      implicit none

*KEEP,magseq.
      integer, parameter :: nmagsp=1000
      integer nmags
      real xplmags(2,5,nmagsp),widmag
      common/magseqc/nmags,xplmags,widmag

*KEND.

      double precision, parameter :: grarad1=0.017453292519943295
      integer, parameter :: linelenp=256

      double precision xzmags(2,5,nmagsp),x,y,z,xmin,xmax,ymin,ymax,zmin,zmax,
     &  scale(6),offset(6),vx,vy,vz,phi,center(3),vin(3),vout(3),vrot(3),zmean,xlen,
     &  rho,angle,fringe,xkz,grad,perlen,dum

      integer lunmg,nlines,icomm,ipos(2,100),nwords,istat,i,im,iop,ianf,iend,l,nfour,nper
      integer :: koper(nmagsp),noper,kmags(nmagsp)

      real xpl(5),zpl(5)

      character(linelenp) cline,cword,clow
      character(linelenp) cbuff(nmagsp),cnam(nmagsp),ctyp(nmagsp)
      character c1

      open(newunit=lunmg,file='magseq.in')

      nlines=0
      icomm=0

      do while (.true.)
        read(lunmg,'(a)',end=9,err=9) cline
        cline=trim(adjustl(cline))
        c1=cline(1:1)
        if (c1.eq.'{') then
          icomm=1
          cycle
        else if (c1.eq.'}') then
          icomm=0
          cycle
        endif
        if (icomm.eq.1.or.c1.eq.''.or.c1.eq.' '.or.c1.eq.'!'.or.c1.eq.'#'.or.c1.eq.'%'.or.c1.eq.'*') cycle
        nlines=nlines+1
        if (nlines.gt.nmagsp) then
          print*,"*** Dimension exceeded in subroutine magseq!"
          print*,"*** Must not exceed ",nmagsp
          print*,"*** Will ignore the follow these lines!"
          nlines=nlines-1
          exit
        endif
        cbuff(nlines)=cline
      enddo

9     close(lunmg)

      noper=0
      nmags=0

      do i=1,nlines
        cline=cbuff(i)
        call util_string_split(cline,linelenp,nwords,ipos,istat)
        cword=cline(ipos(1,1):ipos(2,1))
        clow=cword
        call util_lower_case(clow)
        if (clow.eq.'rotate'.or.clow.eq.'shift') then
          noper=noper+1
          koper(noper)=i
        else
          nmags=nmags+1
          cnam(nmags)=cword
          if (cnam(nmags)(1:1).eq."'".or.cnam(nmags)(1:1).eq.'"') then
            l=len_trim(cnam(nmags))
            cnam(nmags)=cnam(nmags)(2:l-1)
          endif
          ctyp(nmags)=cline(ipos(1,2):ipos(2,2))
          if (ctyp(nmags)(1:1).eq."'".or.ctyp(nmags)(1:1).eq.'"') then
            l=len_trim(ctyp(nmags))
            ctyp(nmags)=ctyp(nmags)(2:l-1)
          endif
          kmags(nmags)=i
        endif
      enddo

      do im=1,nmags
        cline=cbuff(kmags(im))
        call util_string_split(cline,linelenp,nwords,ipos,istat)
        if (ctyp(im).eq.'MAP') then
          xmin=1.0d30
          xmax=-1.0d30
          ymin=1.0d30
          ymax=-1.0d30
          zmin=1.0d30
          zmax=-1.0d30
          scale=1.0d0
          offset=0.0d0
          open(newunit=lunmg,file=cline(ipos(1,3):ipos(2,3)),status='old')
          do while (.true.)
            read(lunmg,'(a)',end=99) cline
            cline=trim(adjustl(cline))
            c1=cline(1:1)
            if (c1.eq.''.or.c1.eq.' '.or.c1.eq.'!'.or.c1.eq.'#'.or.c1.eq.'%'.or.c1.eq.'*') cycle
            call util_string_substring(cline,'scaling',ianf,iend,istat)
            if (istat.eq.0) then
              do i=1,linelenp
                if (cline(i:i).eq.'=') exit
              enddo
              read(cline(i+1:linelenp),*) scale
            endif
            call util_string_substring(cline,'offset',ianf,iend,istat)
            if (istat.eq.0) then
              do i=1,linelenp
                if (cline(i:i).eq.'=') exit
              enddo
              read(cline(i+1:linelenp),*) offset
            endif
            if (c1.eq.'@') cycle
            read(cline,*) x,y,z
            x=scale(1)*x+offset(1)
            y=scale(2)*y+offset(2)
            z=scale(3)*z+offset(3)
            if (x.lt.xmin) xmin=x*100.
            if (x.gt.xmax) xmax=x*100.
            if (y.lt.ymin) ymin=y*100.
            if (y.gt.ymax) ymax=y*100.
            if (z.lt.zmin) zmin=z*100.
            if (z.gt.zmax) zmax=z*100.
          enddo
 99       close(lunmg)
          zmean=(zmax+zmin)/2.0d0
          zmin=zmean-widmag/2.
          zmax=zmean+widmag/2.
          xzmags(1,1,im)=xmin
          xzmags(2,1,im)=zmin
          xzmags(1,2,im)=xmax
          xzmags(2,2,im)=zmin
          xzmags(1,3,im)=xmax
          xzmags(2,3,im)=zmax
          xzmags(1,4,im)=xmin
          xzmags(2,4,im)=zmax
          xzmags(1,5,im)=xmin
          xzmags(2,5,im)=zmin
        else if (ctyp(im).eq.'FOUR') then
          open(newunit=lunmg,file=cline(ipos(1,3):ipos(2,3)),status='old')
          read(lunmg,'(a)',end=99) cline
          read(lunmg,*) xlen
          close(lunmg)
          xlen=xlen*100.
          xmin=-xlen/2.0d0
          xmax=xlen/2.0d0
          zmean=0.0d0
          zmin=zmean-widmag/2.
          zmax=zmean+widmag/2.
          xzmags(1,1,im)=xmin
          xzmags(2,1,im)=zmin
          xzmags(1,2,im)=xmax
          xzmags(2,2,im)=zmin
          xzmags(1,3,im)=xmax
          xzmags(2,3,im)=zmax
          xzmags(1,4,im)=xmin
          xzmags(2,4,im)=zmax
          xzmags(1,5,im)=xmin
          xzmags(2,5,im)=zmin
        else if (ctyp(im).eq.'DIF'.or.ctyp(im).eq.'DHF') then
          read(cline,*,iostat=istat) clow,cword,angle,rho,x,fringe,xkz,nfour,phi
          if (istat.ne.0) then
            read(cline,*,iostat=istat) clow,cword,angle,rho,x,fringe,xkz,nfour
            phi=0.0d0
          endif
          x=x*100.
          rho=rho*100.
          xlen=rho*angle
          xmin=x-xlen/2.0d0
          xmax=x+xlen/2.0d0
          zmean=0.0d0
          zmin=zmean-widmag/2.
          zmax=zmean+widmag/2.
          xzmags(1,1,im)=xmin
          xzmags(2,1,im)=zmin
          xzmags(1,2,im)=xmax
          xzmags(2,2,im)=zmin
          xzmags(1,3,im)=xmax
          xzmags(2,3,im)=zmax
          xzmags(1,4,im)=xmin
          xzmags(2,4,im)=zmax
          xzmags(1,5,im)=xmin
          xzmags(2,5,im)=zmin
          if (phi.ne.0.0d0) then
            center=[x,0.0d0,0.0d0]
            vrot=[0.0d0,1.0d0,0.0d0]
            phi=phi*grarad1
            do i=1,5
              vin=[xzmags(1,i,im),0.0d0,xzmags(2,i,im)]
              call util_rotate(center,vrot,phi,vin,vout,istat)
              xzmags(1,i,im)=vout(1)
              xzmags(2,i,im)=vout(3)
            enddo
          endif
          xzmags(1,:,im)=xzmags(1,:,im)+x
        else if (ctyp(im).eq.'UE') then
          read(cline,*,iostat=istat) clow,cword,dum,dum,dum,dum,x,perlen,nper,dum,dum,dum,dum,z,phi
          x=x*100.0
          perlen=perlen*100.0
          z=z*100.
          xlen=nper*perlen
          zmean=z
          zmin=zmean-widmag/2.
          zmax=zmean+widmag/2.
          xzmags(1,1,im)=xmin
          xzmags(2,1,im)=zmin
          xzmags(1,2,im)=xmax
          xzmags(2,2,im)=zmin
          xzmags(1,3,im)=xmax
          xzmags(2,3,im)=zmax
          xzmags(1,4,im)=xmin
          xzmags(2,4,im)=zmax
          xzmags(1,5,im)=xmin
          xzmags(2,5,im)=zmin
          if (phi.ne.0.0d0) then
            center=[x,0.0d0,0.0d0]
            vrot=[0.0d0,1.0d0,0.0d0]
            phi=phi*grarad1
            do i=1,5
              vin=[xzmags(1,i,im),0.0d0,xzmags(2,i,im)]
              call util_rotate(center,vrot,phi,vin,vout,istat)
              xzmags(1,i,im)=vout(1)
              xzmags(2,i,im)=vout(3)
            enddo
          endif
          xzmags(1,:,im)=xzmags(1,:,im)+x
          xzmags(2,:,im)=xzmags(2,:,im)+z
        else if (ctyp(im).eq.'QP'.or.ctyp(im).eq.'SX'.or.ctyp(im).eq.'QF') then
          read(cline,*,iostat=istat) clow,cword,xlen,grad,x,z,phi
          if (istat.ne.0) then
            read(cline,*,iostat=istat) clow,cword,xlen,grad,x,z
            phi=0.0d0
          endif
          xlen=xlen*100.0
          x=x*100.0
          z=z*100.0
          xmin=x-xlen/2.0d0
          xmax=x+xlen/2.0d0
          zmean=z
          zmin=zmean-widmag/2.
          zmax=zmean+widmag/2.
          xzmags(1,1,im)=xmin
          xzmags(2,1,im)=zmin
          xzmags(1,2,im)=xmax
          xzmags(2,2,im)=zmin
          xzmags(1,3,im)=xmax
          xzmags(2,3,im)=zmax
          xzmags(1,4,im)=xmin
          xzmags(2,4,im)=zmax
          xzmags(1,5,im)=xmin
          xzmags(2,5,im)=zmin
        else if (ctyp(im).eq.'DI'.or.ctyp(im).eq.'DH') then
          read(cline,*,iostat=istat) clow,cword,angle,rho,x,fringe,phi
          if (istat.ne.0) then
            read(cline,*,iostat=istat) clow,cword,angle,rho,x,fringe
            phi=0.0d0
          endif
          rho=rho*100.
          x=x*100.
          xlen=rho*angle
          xmin=x-xlen/2.0d0
          xmax=x+xlen/2.0d0
          zmean=0.0d0
          zmin=zmean-widmag/2.
          zmax=zmean+widmag/2.
          xzmags(1,1,im)=xmin
          xzmags(2,1,im)=zmin
          xzmags(1,2,im)=xmax
          xzmags(2,2,im)=zmin
          xzmags(1,3,im)=xmax
          xzmags(2,3,im)=zmax
          xzmags(1,4,im)=xmin
          xzmags(2,4,im)=zmax
          xzmags(1,5,im)=xmin
          xzmags(2,5,im)=zmin
        endif

        do iop=1,noper
          read(cbuff(koper(iop)),*) clow,cword
          if (cword.ne.cnam(im)) cycle
          call util_lower_case(clow)
          if (clow.eq.'shift') then
            read(cbuff(koper(iop)),*) clow,cword,x,y,z
            xzmags(1,1:5,im)=xzmags(1,1:5,im)+x*100.
            xzmags(2,1:5,im)=xzmags(2,1:5,im)+z*100.
          endif
          if (clow.eq.'rotate') then
            read(cbuff(koper(iop)),*) clow,cword,x,y,z,vx,vy,vz,phi
            center=[x,y,z]
            vrot=[vx,vy,vz]
            phi=phi*grarad1
            do i=1,5
              vin=[xzmags(1,i,im),0.0d0,xzmags(2,i,im)]
              call util_rotate(center,vrot,phi,vin,vout,istat)
              xzmags(1,i,im)=vout(1)
              xzmags(2,i,im)=vout(3)
            enddo
          endif
        enddo

        xpl(1:5)=sngl(xzmags(1,1:5,im))
        zpl(1:5)=sngl(xzmags(2,1:5,im))
        call mshplt_pline(5,xpl,zpl)

      enddo !magnets

      return
      end
*CMZ :          26/01/2025  13.00.20  by  Michael Scheer
*-- Author :    Michael Scheer   26/01/2025
      subroutine hpline(x,y,n,copt)

      real x(*),y(*)
      character(*) copt
      character(32) c

      c=copt

      call mshplt_pline(n,x,y)

      return
      end
*CMZ : 00.00/20 06/12/2016  18.30.31  by  Michael Scheer
*CMZ : 00.00/16 19/03/2014  12.14.29  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.29.49  by  Michael Scheer
*CMZ : 00.00/06 08/03/2007  14.02.27  by  Michael Scheer
*CMZ : 00.00/05 07/03/2007  12.58.44  by  Michael Scheer
*-- Author :    Michael Scheer   07/03/2007
      subroutine util_string_substring(cline,substring,ianf,iend,istat)

c Input:
c      cline, substring

c If substring is passed as variable, full length of substring is checked,
c i.e. pending invisible characters are tested as well
c Probably, you want to test trim(substring)!

c Output:
c      ianf,iend: start and end position of substring, 0 if not found
c      istat: error, i.e. string not found

c Evtl. besser FORTRAN-functions scan oder index benutzen

      implicit none

        integer ilenl,ilens,istat,ianf,iend,i

        character(*) cline,substring

        istat=-1
        ianf=0
        iend=0

        ilenl=len(cline)
        ilens=len(substring)

        if (ilens.gt.ilenl) return

        do i=1,ilenl-ilens+1
          if (cline(i:i+ilens-1).eq.substring) then
            ianf=i
            iend=ianf+ilens-1
            istat=0
            return
          endif
        enddo

      return
      end
*CMZ :          14/10/2021  05.12.31  by  Michael Scheer
*-- Author :    Michael Scheer   06/04/2018
      subroutine util_rotate(cen,vrot,phi,vin,vout,istat)

c +PATCH,//UTIL/FOR
c +DECK,util_rotate.
      implicit none

      double precision cen(3),phi,vin(3),vout(3),vrot(3),vlen,o(3),
     &  r(3),s,c,c1,rm(3,3)

      integer istat

      istat=0

      vlen=norm2(vrot)
      if (vlen.eq.0.0d0) then
        vout=vin
        istat=1
        return
      endif

      o=vrot/vlen
      s=sin(phi)
      c=cos(phi)
      c1=1.0d0-c

      rm(1,1)=o(1)*o(1)*c1+c
      rm(2,2)=o(2)*o(2)*c1+c
      rm(3,3)=o(3)*o(3)*c1+c

      rm(1,2)=o(1)*o(2)*c1-o(3)*s
      rm(1,3)=o(1)*o(3)*c1+o(2)*s

      rm(2,1)=o(1)*o(2)*c1+o(3)*s
      rm(2,3)=o(2)*o(3)*c1-o(1)*s

      rm(3,1)=o(1)*o(3)*c1-o(2)*s
      rm(3,2)=o(2)*o(3)*c1+o(1)*s

      r=vin-cen

      vout(1)=rm(1,1)*r(1)+rm(1,2)*r(2)+rm(1,3)*r(3)+cen(1)
      vout(2)=rm(2,1)*r(1)+rm(2,2)*r(2)+rm(2,3)*r(3)+cen(2)
      vout(3)=rm(3,1)*r(1)+rm(3,2)*r(2)+rm(3,3)*r(3)+cen(3)

      return
      end
*CMZ :          08/03/2018  14.14.52  by  Michael Scheer
*CMZ : 00.00/16 13/10/2014  09.07.28  by  Michael Scheer
*-- Author :    Michael Scheer   13/10/2014
      subroutine util_lower_case(cline)

      implicit none

      character(*) cline

      integer i,ic1
      character c1
      equivalence (ic1,c1)

      ic1=0

      do i=1,len_trim(cline)
        c1=cline(i:i)
        if (ic1.ge.65.and.ic1.le.90) ic1=ic1+32
        cline(i:i)=c1
      enddo

      return
      end
