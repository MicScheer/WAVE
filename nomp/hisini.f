*CMZ :  4.01/02 14/05/2023  13.16.19  by  Michael Scheer
*CMZ :  4.01/00 09/01/2023  20.53.38  by  Michael Scheer
*CMZ :  4.00/15 08/04/2022  17.00.19  by  Michael Scheer
*CMZ :  4.00/14 31/12/2021  08.37.37  by  Michael Scheer
*CMZ :  4.00/13 20/12/2021  16.34.10  by  Michael Scheer
*CMZ :  4.00/11 22/04/2021  19.01.32  by  Michael Scheer
*CMZ :  4.00/06 09/12/2019  14.36.37  by  Michael Scheer
*CMZ :  4.00/04 21/08/2019  14.35.31  by  Michael Scheer
*CMZ :  3.05/10 13/08/2018  14.40.26  by  Michael Scheer
*CMZ :  3.03/02 26/02/2016  10.43.50  by  Michael Scheer
*CMZ :  3.03/00 17/08/2015  10.03.05  by  Michael Scheer
*CMZ :  3.02/06 08/06/2015  13.41.37  by  Michael Scheer
*CMZ :  3.02/05 07/04/2015  15.57.59  by  Michael Scheer
*CMZ :  3.02/04 18/03/2015  09.26.30  by  Michael Scheer
*CMZ :  3.02/02 16/10/2014  10.56.27  by  Michael Scheer
*CMZ :  3.02/00 24/09/2014  11.52.00  by  Michael Scheer
*CMZ :  3.01/07 23/06/2014  13.59.37  by  Michael Scheer
*CMZ :  3.01/06 20/06/2014  16.49.38  by  Michael Scheer
*CMZ :  3.01/05 13/06/2014  11.46.51  by  Michael Scheer
*CMZ :  3.00/00 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.70/08 14/01/2013  17.46.40  by  Michael Scheer
*CMZ :  2.70/07 14/01/2013  16.55.40  by  Michael Scheer
*CMZ :  2.70/06 03/01/2013  13.55.49  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  15.10.55  by  Michael Scheer
*CMZ :  2.70/03 14/12/2012  16.42.06  by  Michael Scheer
*CMZ :  2.70/00 08/11/2012  12.18.50  by  Michael Scheer
*CMZ :  2.68/05 16/10/2012  17.25.52  by  Michael Scheer
*CMZ :  2.68/02 26/06/2012  11.55.39  by  Michael Scheer
*CMZ :  2.68/01 29/05/2012  16.43.49  by  Michael Scheer
*CMZ :  2.67/06 23/05/2012  15.50.53  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.51.59  by  Michael Scheer
*CMZ :  2.67/02 09/05/2012  09.56.14  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  13.01.28  by  Michael Scheer
*CMZ :  2.66/12 24/06/2010  12.50.52  by  Michael Scheer
*CMZ :  2.66/07 10/03/2010  09.23.32  by  Michael Scheer
*CMZ :  2.66/06 24/11/2009  09.57.50  by  Michael Scheer
*CMZ :  2.65/02 23/09/2009  08.34.24  by  Michael Scheer
*CMZ :  2.61/06 12/04/2007  13.36.22  by  Michael Scheer
*CMZ :  2.61/03 27/03/2007  13.00.05  by  Michael Scheer
*CMZ :  2.52/16 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.48/04 17/03/2004  13.52.09  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  2.36/00 06/11/2001  16.31.58  by  Michael Scheer
*CMZ :  2.17/00 02/11/2000  16.12.50  by  Michael Scheer
*CMZ :  2.16/08 20/10/2000  11.44.49  by  Michael Scheer
*CMZ :  2.16/04 26/06/2000  16.57.45  by  Michael Scheer
*CMZ :  2.13/11 20/03/2000  12.32.30  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.26.25  by  Michael Scheer
*CMZ :  2.12/00 27/05/99  10.11.30  by  Michael Scheer
*CMZ :  2.10/01 24/02/99  12.04.52  by  Michael Scheer
*CMZ :  1.04/00 16/10/98  10.54.26  by  Michael Scheer
*CMZ :  1.03/06 29/09/98  14.00.01  by  Michael Scheer
*CMZ :  1.03/03 27/03/98  13.37.41  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  10.31.28  by  Michael Scheer
*CMZ : 00.02/00 25/11/96  14.46.50  by  Michael Scheer
*CMZ : 00.01/08 18/07/95  11.55.34  by  Michael Scheer
*CMZ : 00.01/07 10/03/95  13.17.24  by  Michael Scheer
*CMZ : 00.01/04 30/01/95  11.29.57  by  Michael Scheer
*CMZ : 00.00/05 29/04/94  19.50.41  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  18.06.58  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.28  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE HISINI
*KEEP,gplhint.
*KEND.

*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

      use bunchmod
      use mhbook_mod
      use clustermod

C--- INITIALIZES HBOOK

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,klotz.
      include 'klotz.cmn'
*KEEP,phasef90.
      include 'phasef90.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,photon.
      include 'photon.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,ntuple2.
      include 'ntuple2.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEND.

      INTEGER ISTAT,ICYCLE,ISCRATCH,LREC

      INTEGER IGL,I,J,iopen,last,iline,iblank
      BYTE I1(96)

      CHARACTER(32) LINE_1,LINE_2,LINE_3
      CHARACTER(1) LINE32_1(32),LINE32_2(32),LINE32_3(32)
      CHARACTER(96) C96,EMPTY,cwavein(5000)
      CHARACTER(1) C1(96)
      CHARACTER(4) CHTAGS(41),chtagsg1(100) !if you change the dimension of chtags here, check all routines with NIDBUNCH
      character(4) chtagsus(ntupusp)

      EQUIVALENCE(I1,C1,C96)
      EQUIVALENCE(LINE_1,LINE32_1)
      EQUIVALENCE(LINE_2,LINE32_2)
      EQUIVALENCE(LINE_3,LINE32_3)

      DATA EMPTY/' '/
      data iopen/0/,icycle/0/

      COMMON/CLINEc/LINE_1,LINE_2,LINE_3

      data chtags/'ibun','iel','nel'
     &  ,'xbun','xi','yi','zi','ypi','zpi'
     &  ,'xe','ye','ze','ype','zpe','ei','ee','xobs','yobs','zobs',
     &  'iene','ener','spec','s0','s1','s2','s3','pow','isou','te',
     &  'exr','exi','eyr','eyi','ezr','ezi',
     &  'bxr','bxi','byr','byi','bzr','bzi'
     &  /

      data chtagsus /
     &  'u1  ',
     &  'u2  ',
     &  'u3  ',
     &  'u4  ',
     &  'u5  ',
     &  'u6  ',
     &  'u7  ',
     &  'u8  ',
     &  'u9  ',
     &  'u10 '
     &  /

      character(4) chtagsip(ntupip) /
     &  'ebea',
     &  'curr',
     &  'ipin',
     &  'icir',
     &  'pinx',
     &  'piny',
     &  'pinz',
     &  'pinw',
     &  'pinh',
     &  'pinr',
     &  'mpiz',
     &  'mpiy',
     &  'mpir',
     &  'mpip',
     &  'icbr',
     &  'selx',
     &  'sely',
     &  'selz',
     &  'phax',
     &  'sigz',
     &  'sigy',
     &  'sgzp',
     &  'sgyp',
     &  'espr',
     &  'if2p',
     &  'nfrq',
     &  'flow',
     &  'fhig',
     &  'ispe',
     &  'ispm',
     &  'idip',
     &  'nlpo',
     &  'bw',
     &  'ibun',
     &  'nbun',
     &  'neib',
     &  'iamp',
     &  'ielo',
     &  'ifol',
     &  'iefo',
     &  'irun',
     &  'idat',
     &  'itim',
     &  'vers',
     &  'isto',
     &  'noby',
     &  'nobz',
     &  'wal1',
     &  'wal2',
     &  'xabs',
     &  'zab1',
     &  'zab2',
     &  'ibri',
     &  'ibet',
     &  'kamp'
     &  /

      IF (IHISINI_C.gt.0) RETURN

      VMX=0.

      IDCODE=1
      IDINPUT=12+ABS(IHINPUT)
      IDOUTP=16

      IDTRACKZ=1100  !TRAJECTORY (Z-COORDINATE)
      IDTRACKY=1200  !TRAJECTORY (Y-COORDINATE)

      IDTRCKSZ=1100  !TRAJECTORY OF SOURCE POINTS (Z-COORDINATE)
      IDTRCKSY=1200  !TRAJECTORY OF SOURCE POINTS (Y-COORDINATE)

      IDBX=2100   !MAGNETIC FIELD (X-COMPONENT)
      IDBY=2200   !MAGNETIC FIELD (Y-COMPONENT)
      IDBZ=2300   !MAGNETIC FIELD (Z-COMPONENT)

      IDSBX=2100  !MAGNETIC FIELD OF SOURCE (X-COMPONENT)
      IDSBY=2200  !MAGNETIC FIELD OF SOURCE (Y-COMPONENT)
      IDSBZ=2300  !MAGNETIC FIELD OF SOURCE (Z-COMPONENT)

      IDAX=4100   !MAGNETIC VECTOR POTENTIAL (X-COMPONENT)
      IDAY=4200   !MAGNETIC VECTOR POTENTIAL (Y-COMPONENT)
      IDAZ=4300   !MAGNETIC VECTOR POTENTIAL (Z-COMPONENT)

      IDBETAH= 3100  !HORIZONTAL BETA-FUNCTION, NUMERICALLY
      IDBETAV= 3200  !VERTICAL BETA-FUNCTION, NUMERICALLY
      IDBETHPC=3110  !HORIZONTAL BETA-FUNCTION, ANALYTICALLY WITH
C            ! BETA=(BETA0+X**2/BETA0)*(COS(THETA))**2
      IDBETHP=3120  !HORIZONTAL BETA-FUNCTION, ANALYTICALLY WITH
C           ! BETA=(BETA0+X**2/BETA0)
      IDBETAPH=3300  !DERIVATIVE OF HORIZONTAL BETA-FUNCTION
      IDBPHPC= 3310  !DERIVATIVE OF
C            !HORIZONTAL BETA-FUNCTION, ANALYTICALLY WITH
C            ! BETA=(BETA0+X**2/BETA0)*(COS(THETA))**2
      IDBPHP= 3320  !DERIVATIVE OF
C            !HORIZONTAL BETA-FUNCTION, ANALYTICALLY WITH
C            ! BETA=(BETA0+X**2/BETA0)
      IDBETAPV=3400  !DERIVATIVE OF VERTICAL BETA-FUNCTION
      IDETA=   3450  !DISPERSION
      IDHS=    3800  !H(s), CHROMATIC VARIABLE
      IDHX=    3801  !H(x), CHROMATIC VARIABLE
      IDHSR3=  3802  !H(s)/abs(rho)**3, NUMERICALLY
      IDHXR3=  3803  !H(x)/abs(rho)**3, ANALYTICALLY WITH PARABOLIC ANSATZ

      IDSPEC=11000   !INTENSITY IN PINHOLE
      IDSPECH=12000  !INTENSITY IN ORBIT PLANE
      IDSPECV=13000  !INTENSITY VERTICAL TO ORBIT PLANE

      IDSPECF= 21000 !FOLDED INTENSITY IN PINHOLE
      IDSPECFH=22000 !FOLDED INTENSITY IN ORBIT PLANE
      IDSPECFV=23000 !FOLDED INTENSITY VERTICAL TO ORBIT PLANE
      IDFOLFNZ=24000 !HORIZONTAL FOLDING FUNCTION
      IDFOLFNY=25000 !VERTICAL FOLDING FUNCTION

      IDSTOKES0=31000   !STOKES S0 IN PINHOLE
      IDSTOKES1=31100   !STOKES S1 IN PINHOLE
      IDSTOKES2=31200   !STOKES S2 IN PINHOLE
      IDSTOKES3=31300   !STOKES S3 IN PINHOLE

      IDSTOKFS0=41000   !STOKES S0 IN PINHOLE (FOLDED)
      IDSTOKFS1=41100   !STOKES S1 IN PINHOLE (FOLDED)
      IDSTOKFS2=41200   !STOKES S2 IN PINHOLE (FOLDED)
      IDSTOKFS3=41300   !STOKES S3 IN PINHOLE (FOLDED)


      IDFREQ=48000  !SPECTRUM
      IDFREQF=49000  !FOLDED SPECTRUM

      IDFRS0=50000   !SPECTRUM OF S0
      IDFRS1=51000   !SPECTRUM OF S1
      IDFRS2=52000   !SPECTRUM OF S2
      IDFRS3=53000   !SPECTRUM OF S3
      IDFRP =54000   !SPECTRUM OF P
      IDFRP1=54500   !SPECTRUM OF P1
      IDFRP2=55000   !SPECTRUM OF P2
      IDFRP3=56000   !SPECTRUM OF P3
      IDFRP23=57000   !SPECTRUM OF P23
      IDFRG3=58000   !SPECTRUM OF S0*P3**2
      IDFRG23=59000   !SPECTRUM OF S0*P23**2

      IDFRS0F=60000  !FOLDED SPECTRUM OF S0
      IDFRS1F=61000  !FOLDED SPECTRUM OF S1
      IDFRS2F=62000  !FOLDED SPECTRUM OF S2
      IDFRS3F=63000  !FOLDED SPECTRUM OF S3
      IDFRPF =64000  !FOLDED SPECTRUM OF P
      IDFRP1F=64500  !FOLDED SPECTRUM OF P1
      IDFRP2F=65000  !FOLDED SPECTRUM OF P2
      IDFRP3F=66000  !FOLDED SPECTRUM OF P3
      IDFRP23F=67000  !FOLDED SPECTRUM OF P23
      IDFRG3F=68000  !FOLDED SPECTRUM OF S0*P3**2
      IDFRG23F=69000  !FOLDED SPECTRUM OF S0*P23**2

      IDFRS0E=70000 !SPECTRUM OF S0, E-FOLDED
      IDFRS1E=71000 !SPECTRUM OF S1, E-FOLDED
      IDFRS2E=72000 !SPECTRUM OF S2, E-FOLDED
      IDFRS3E=73000 !SPECTRUM OF S3, E-FOLDED
      IDFRPE =74000 !SPECTRUM OF P, E-FOLDED
      IDFRP1E=74500 !SPECTRUM OF P1, E-FOLDED
      IDFRP2E=75000 !SPECTRUM OF P2, E-FOLDED
      IDFRP3E=76000 !SPECTRUM OF P3, E-FOLDED
      IDFRP23E=77000 !SPECTRUM OF P23, E-FOLDED
      IDFRG3E=78000 !SPECTRUM OF S0*P3**2, E-FOLDED
      IDFRG23E=79000 !SPECTRUM OF S0*P23**2, E-FOLDED

      IDFRS0EF=80000  !FOLDED SPECTRUM OF S0, E-FOLDED
      IDFRS1EF=81000  !FOLDED SPECTRUM OF S1, E-FOLDED
      IDFRS2EF=82000  !FOLDED SPECTRUM OF S2, E-FOLDED
      IDFRS3EF=83000  !FOLDED SPECTRUM OF S3, E-FOLDED
      IDFRPEF =84000  !FOLDED SPECTRUM OF P , E-FOLDED
      IDFRP1EF=84500  !FOLDED SPECTRUM OF P1, E-FOLDED
      IDFRP2EF=85000  !FOLDED SPECTRUM OF P2, E-FOLDED
      IDFRP3EF=86000  !FOLDED SPECTRUM OF P3, E-FOLDED
      IDFRP23EF=87000  !FOLDED SPECTRUM OF P23, E-FOLDED
      IDFRG3EF=88000  !FOLDED SPECTRUM OF S0*P3**2, E-FOLDED
      IDFRG23EF=89000  !FOLDED SPECTRUM OF S0*P23**2, E-FOLDED

      ICFREQ=148000  !SPECTRUM
      ICFREQF=149000  !FOLDED SPECTRUM

      ICFRS0=150000   !CENTRAL SPECTRUM OF S0
      ICFRS1=151000   !CENTRAL SPECTRUM OF S1
      ICFRS2=152000   !CENTRAL SPECTRUM OF S2
      ICFRS3=153000   !CENTRAL SPECTRUM OF S3
      ICFRP =154000   !CENTRAL SPECTRUM OF P
      ICFRP1=154500   !CENTRAL SPECTRUM OF P1
      ICFRP2=155000   !CENTRAL SPECTRUM OF P2
      ICFRP3=156000   !CENTRAL SPECTRUM OF P3
      ICFRP23=157000   !CENTRAL SPECTRUM OF P23
      ICFRG3=158000   !CENTRAL SPECTRUM OF S0*P3**2
      ICFRG23=159000   !CENTRAL SPECTRUM OF S0*P23**2

      ICFRS0F=160000  !FOLDED CENTRAL SPECTRUM OF S0
      ICFRS1F=161000  !FOLDED CENTRAL SPECTRUM OF S1
      ICFRS2F=162000  !FOLDED CENTRAL SPECTRUM OF S2
      ICFRS3F=163000  !FOLDED CENTRAL SPECTRUM OF S3
      ICFRP1F=164500  !FOLDED CENTRAL SPECTRUM OF P1
      ICFRPF =164000  !FOLDED CENTRAL SPECTRUM OF P
      ICFRP2F=165000  !FOLDED CENTRAL SPECTRUM OF P2
      ICFRP3F=166000  !FOLDED CENTRAL SPECTRUM OF P3
      ICFRP23F=167000  !FOLDED CENTRAL SPECTRUM OF P23
      ICFRG3F=168000  !FOLDED CENTRAL SPECTRUM OF S0*P3**2
      ICFRG23F=169000  !FOLDED CENTRAL SPECTRUM OF S0*P23**2

      ICFRS0E=170000 !CENTRAL SPECTRUM OF S0, E-FOLDED
      ICFRS1E=171000 !CENTRAL SPECTRUM OF S1, E-FOLDED
      ICFRS2E=172000 !CENTRAL SPECTRUM OF S2, E-FOLDED
      ICFRS3E=173000 !CENTRAL SPECTRUM OF S3, E-FOLDED
      ICFRPE =174000 !CENTRAL SPECTRUM OF P, E-FOLDED
      ICFRP1E=174500 !CENTRAL SPECTRUM OF P, E-FOLDED
      ICFRP2E=175000 !CENTRAL SPECTRUM OF P2, E-FOLDED
      ICFRP3E=176000 !CENTRAL SPECTRUM OF P3, E-FOLDED
      ICFRP23E=177000 !CENTRAL SPECTRUM OF P23, E-FOLDED
      ICFRG3E=178000 !CENTRAL SPECTRUM OF S0*P3**2, E-FOLDED
      ICFRG23E=179000 !CENTRAL SPECTRUM OF S0*P23**2, E-FOLDED

      ICFRS0EF=180000  !FOLDED CENTRAL SPECTRUM OF S0, E-FOLDED
      ICFRS1EF=181000  !FOLDED CENTRAL SPECTRUM OF S1, E-FOLDED
      ICFRS2EF=182000  !FOLDED CENTRAL SPECTRUM OF S2, E-FOLDED
      ICFRS3EF=183000  !FOLDED CENTRAL SPECTRUM OF S3, E-FOLDED
      ICFRPEF =184000  !FOLDED CENTRAL SPECTRUM OF P , E-FOLDED
      ICFRP1EF=184500  !FOLDED CENTRAL SPECTRUM OF P1, E-FOLDED
      ICFRP2EF=185000  !FOLDED CENTRAL SPECTRUM OF P2, E-FOLDED
      ICFRP3EF=186000  !FOLDED CENTRAL SPECTRUM OF P3, E-FOLDED
      ICFRP23EF=187000  !FOLDED CENTRAL SPECTRUM OF P23, E-FOLDED
      ICFRG3EF=188000  !FOLDED CENTRAL SPECTRUM OF S0*P3**2, E-FOLDED
      ICFRG23EF=189000  !FOLDED CENTRAL SPECTRUM OF S0*P23**2, E-FOLDED

      ICFRB0=190000   !CENTRAL BRILLIANCE
      ICFRB1=191000   !CENTRAL BRILLIANCE
      ICFRB2=192000   !CENTRAL BRILLIANCE
      ICFRB3=193000   !CENTRAL BRILLIANCE

      ICFRB0E=194000   !CENTRAL BRILLIANCE
      ICFRB1E=195000   !CENTRAL BRILLIANCE
      ICFRB2E=196000   !CENTRAL BRILLIANCE
      ICFRB3E=197000   !CENTRAL BRILLIANCE

      ICFRB0F=190500   !CENTRAL BRILLIANCE
      ICFRB1F=191500   !CENTRAL BRILLIANCE
      ICFRB2F=192500   !CENTRAL BRILLIANCE
      ICFRB3F=193500   !CENTRAL BRILLIANCE

      ICFRB0EF=194500   !CENTRAL BRILLIANCE
      ICFRB1EF=195500   !CENTRAL BRILLIANCE
      ICFRB2EF=196500   !CENTRAL BRILLIANCE
      ICFRB3EF=197500   !CENTRAL BRILLIANCE

      NIDTRACK=10   !NTUPLE CONTAINS ARRAY WTRA
      NIDMINI=11   !NTUPLE CONTAINS MINI-TRAJECTORIES
      NIDSOURCE=12   !NTUPLE CONTAINS INTEGRANTS OF SOURCES
      NIDTRCKG=20  !NTUPLE CONTAINS ARRAY WTRA
      NIDPHANA=21  !NTUPLE PHASE ERROR ANALYSIS
      nidbunch=30  !ntuple of bunch electrons

      NIDTRAC=500   !X,XP,Y,YP AT START AND STOP OF TRAJECTORIES IN SR OPTI
      NIDFREQ=3500  !SPECTRUM FOR SINGLE OBSERVATION POINTS, NTUPLE
      NIDFREQP=3600  !SPECTRUM FOR PINHOLE, NTUPLE
      NIDSPEC=3700  !SPECTRUM-ARRAY FOR ALL OBSERVATION POINTS, NTUPLE
      NIDSPECRPHI=5700  !SPECTRUM-ARRAY FOR ALL OBSERVATION POINTS, NTUPLE
      NIDSPECF=3701  !FOLDED SPECTRUM-ARRAY FOR ALL OBSERVATION POINTS, NTUPLE
      NIDPOW=2000  !NTUPLE FOR POWER-DENSITY
      NIDSTOK=4700  !STOKES-ARRAYS FOR ALL OBSERVATION POINTS, NTUPLE
      NIDSTOKF=4701  !FOLDED STOKES-ARRAYS FOR ALL OBSERVATION POINTS, NTUPLE
      NIDSTOKE=4702  !E-FOLDED STOKES-ARRAYS FOR ALL OBSERVATION POINTS, NTUPLE
      NIDSTOKEF=4703  !FOLDED AND E-FOLDED STOKES-ARRAYS FOR ALL OBSERVATION POINTS, NTUPLE

      IDPOWER=40000  !POWER DENSITY ALONG BEAMLINE

      NIDFIL=5000    !PHOTO ABSORPTION COEFFICIENTS OF ABSORBING FILTER
      NIDMUEN=5001   !ENERGY ABSORPTION COEFFICIENTS
      NIDPHASE=6000   !PHASE
      IDPHASE=206000 !PHASE
      NIDGEO=6100   !PHASE SPACE FROM GEOMETRICAL OPTIC
      IDSEL=206100   !PHASE SPACE RATIO

C*** THIS 600 IS USED IN SR BMESSINI AS 599 (600-1), CHANGE ACCORDINGLY
C    IF NECESSARY
      NIDBPOLY=600 !field map and fitted field of SR BPOLY3DFIT

C---    THE PARAMETER NDPAWCP INQUIRES THE WORK SPACE FOR HBOOK
C     HISTOGRAM 1 IS USED TO IDENTIFY THE HISTOGRAM FILE, IT CONTAINS
C     THE RUN NUMBER OF WAVE

c      IF (IHLIMIT_C.EQ.0..and.iroottrees.ge.0) then
      IF (IHLIMIT_C.EQ.0) THEN
        CALL hlimitm(mhbookp)
        if (icluster.lt.0) then
          open(newunit=nscr7777,file='msh_hbook_ntuple_7777.scr')
          open(newunit=nscr3600,file='msh_hbook_ntuple_3600.scr')
          open(newunit=nscr3601,file='msh_hbook_ntuple_3601.scr')
          open(newunit=nscr3700,file='msh_hbook_ntuple_3700.scr')
          open(newunit=nscr4600,file='msh_hbook_ntuple_4600.scr')
          open(newunit=nscr4700,file='msh_hbook_ntuple_4700.scr')
          open(newunit=nscr30,file='msh_hbook_ntuple_30.scr')
        endif
        IHLIMIT_C=1
      ENDIF

      if (iopen.eq.0) then
        if (iroottrees.ge.0) then
          CALL hropenm(LUNHB,'WAVE',FILEHB,'NQ',LREC,ISTAT)
          CALL htitlem(CODE)
        endif
        last=len_trim(CODE)
        call hbook1m(IDCODE,CODE(1:last)//'',2,0.5,2.5,0.)
        CALL hfillm(IDCODE,1.,0.,0.0d0)
        CALL hfillm(IDCODE,2.,0.,ABS(dble(ICODE)))
        CALL MHROUT(IDCODE,ICYCLE,' ')
        iopen=1
      endif

      if (ihisini_c.eq.-2) then
        idg1=777777
        nidphotons=7777
        chtagsg1(1)='trck'
        chtagsg1(2)='step'
        chtagsg1(3)='x'
        chtagsg1(4)='y'
        chtagsg1(5)='z'
        chtagsg1(6)='pxe'
        chtagsg1(7)='pye'
        chtagsg1(8)='pze'
        chtagsg1(9)='ee'
        chtagsg1(10)='pxg'
        chtagsg1(11)='pyg'
        chtagsg1(12)='pzg'
        chtagsg1(13)='eg'
        chtagsg1(14)='ec'
        call hbookm(nidphotons,
     &    'quantum fluctuation photons',14,'//WAVE',2**16,chtagsg1)
        ihisini_c=-1
        return
      endif

      if (ihisini_c.gt.0) return

      nid222=222
      CALL hbookm(nid222,'Input Parameters',ntupip,'//WAVE',1024,CHTAGSip)
      nid223=223
      CALL hbookm(nid223,'Variables of user common block',ntupusp,
     &  '//WAVE',1024,CHTAGSus)

      if (iroottrees.ge.0) then
        IF (IHINPUT.NE.0) then
          CALL hbntm(IDINPUT,'INPUT FILE WAVE.IN',' ')
        endif
        IF (ABS(IHINPUT).EQ.1) THEN
          CALL hbnamcm(IDINPUT,'LINE',
     &      LINE_1,'LINE_1:C*32')
        ELSE IF (ABS(IHINPUT).EQ.2) THEN
          CALL hbnamcm(IDINPUT,'LINE',
     &      LINE_1,'LINE_1:C*32,LINE_2:C*32')
        ELSE IF (ABS(IHINPUT).EQ.3) THEN
          CALL hbnamcm(IDINPUT,'LINE',
     &      LINE_1,'LINE_1:C*32,LINE_2:C*32,LINE_3:C*32')
        ELSE
          GOTO 999
        ENDIF
      endif
C        OPEN(UNIT=LUNGFI,FILE=FILEI,STATUS='OLD',READONLY)

      REWIND(LUNGFI)

      iline=0
100   C96=EMPTY

      IF (IHINPUT.ne.0) THEN

        READ (LUNGFI,'(A96)',END=90) C96

C- SCAN LINES FOR COMMENTS ETC.

        IGL=96
        DO I=1,96
          IF (C1(I).EQ.'=') THEN
            IGL=I
            GOTO 9
          ENDIF
        ENDDO

9       CONTINUE

        DO I=1,IGL-1
          IF (C1(I).EQ.'') GOTO 15
          IF (C1(I).EQ.'!'.and.ihinput.lt.0) GOTO 100
        ENDDO

15      CONTINUE

C- REPLACE TAB BY ONE BLANK

        iblank=1
        DO I=1,96
          IF (I1(I).EQ.9) C1(I)=' '
          if (c1(i).ne.' ') iblank=0
        ENDDO

        if (iblank.eq.1) goto 100

C- SUPRESS MULTIPLE BLANKS

        J=0
        C1(96)=' '
        DO I=1,96-1
          IF (C1(I).NE.' '.OR.C1(I+1).NE.' ') THEN
            J=J+1
            C1(J)=C1(I)
          ENDIF
        ENDDO

        do i=j+1,96
          C1(i)=' '
        enddo

C- FILL CWN BUFFER

20      CONTINUE

        DO I=1,32
          LINE32_1(I)=C1(I)
          LINE32_2(I)=C1(I+32)
          LINE32_3(I)=C1(I+64)
        ENDDO

        if(iroottrees.ge.0) CALL hfntm(IDINPUT)

        iline=iline+1

        if (abs(ihinput).eq.1) then
          cwavein(iline)=trim(line_1)
        else if (abs(ihinput).eq.2) then
          cwavein(iline)=trim(line_1) // trim(line_2)
        else if (abs(ihinput).eq.3) then
          cwavein(iline)=trim(line_1) // trim(line_2) // trim(line_3)
        endif

        GOTO 100

90      CONTINUE

c        CALL MHROUT(IDINPUT,ICYCLE,' ')

        if (iroottrees.ge.0) then
          DO ISCRATCH=1,ICYCLE-1
            CALL HSCRm(IDINPUT,ISCRATCH,' ')
          ENDDO
          CALL hdeletnomh(IDINPUT)
        endif

        CLOSE(LUNGFI)

        write(luns_mh(kfile_mh),'(a)')" ! ------------------------------------------------"
        write(luns_mh(kfile_mh),'(a)')"                   15                    4    ! id and kind of histogram of Ntuple"
        write(luns_mh(kfile_mh),'(a)')"          18    ! length of title"
        write(luns_mh(kfile_mh),'(a)')"INPUT FILE WAVE.IN"
        if (abs(ihinput).eq.1) then
          write(luns_mh(kfile_mh),'(a)')"                    6           1          11    ! length of pathname, number of variables, and length of variable names"
          write(luns_mh(kfile_mh),'(a)')"//WAVE"
          write(luns_mh(kfile_mh),'(a)')"LINE_1:C*32"
        else if (abs(ihinput).eq.2) then
          write(luns_mh(kfile_mh),'(a)')"                    6           2          11    ! length of pathname, number of variables, and length of variable names"
          write(luns_mh(kfile_mh),'(a)')"//WAVE"
          write(luns_mh(kfile_mh),'(a)')"LINE_1:C*32"
          write(luns_mh(kfile_mh),'(a)')"LINE_2:C*32"
        else if (abs(ihinput).eq.3) then
          write(luns_mh(kfile_mh),'(a)')"                    6           3          11    ! length of pathname, number of variables, and length of variable names"
          write(luns_mh(kfile_mh),'(a)')"//WAVE"
          write(luns_mh(kfile_mh),'(a)')"LINE_1:C*32"
          write(luns_mh(kfile_mh),'(a)')"LINE_2:C*32"
          write(luns_mh(kfile_mh),'(a)')"LINE_3:C*32"
        endif
        write(luns_mh(kfile_mh),*)iline," ! number of entries"

        do i=1,iline
          if (abs(ihinput).eq.1) then
            write(luns_mh(kfile_mh),'(a)') cwavein(i)(1:32)
          else if (abs(ihinput).eq.2) then
            write(luns_mh(kfile_mh),'(a)') cwavein(i)(1:64)
          else if (abs(ihinput).eq.3) then
            write(luns_mh(kfile_mh),'(a)') cwavein(i)(1:96)
          endif
        enddo

      ENDIF !IHINPUT.GT.0


999   CONTINUE

      nidbunch=30
      if (ihbunch.ne.0.and.ispec.ne.0) then
        call hbookm(nidbunch,'beam electrons',41,'//WAVE',nbunch*neinbunch*nintfreq,
     &    chtags)
      endif

      IF (XSTARTH.EQ.9999.) XSTARTH=XSTART
      IF (XSTOPH.EQ.9999.)  XSTOPH=XSTOP

      IHISINI_C=1

      RETURN
      END
