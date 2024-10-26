*CMZ :  4.01/03 01/06/2023  08.07.36  by  Michael Scheer
*CMZ :  4.00/14 09/02/2022  16.42.19  by  Michael Scheer
*CMZ :  4.00/07 30/03/2020  16.55.30  by  Michael Scheer
*CMZ :  3.07/00 15/03/2019  12.15.25  by  Michael Scheer
*CMZ :  3.05/02 15/05/2018  08.43.18  by  Michael Scheer
*CMZ :  3.05/01 09/05/2018  09.07.17  by  Michael Scheer
*CMZ :  3.02/06 15/04/2015  11.35.50  by  Michael Scheer
*CMZ :  3.02/03 06/11/2014  14.24.54  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  14.04.56  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.55.57  by  Michael Scheer
*CMZ :  2.50/00 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.49/00 22/03/2004  14.04.26  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.45  by  Michael Scheer
*CMZ :  2.12/00 27/05/99  10.08.55  by  Michael Scheer
*CMZ :  2.10/01 24/02/99  10.20.40  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.18.17  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.54.05  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.44  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE SOUINT_omp(ISOUR,IBUFF)

*KEEP,gplhint.
*KEND.

*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,reargf90u.
      include 'reargf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.
      use ompmod

C--- EVALUATE INTEGRALES FOR A SINGLE SOURCE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,datetime.
      include 'datetime.cmn'
*KEND.

c      double precision reromp(8,ndwsou)
      double precision, dimension (:,:), allocatable :: reromp

      INTEGER ISOUR,IBUFF,IOBSV,JX10,JDX10,IX10,ICAL,NNBUFF,ICYCLE,
     &  ith,kfreq,k,i

      INTEGER NTUPP
      PARAMETER (NTUPP=22)
      CHARACTER(4) CTUP(NTUPP)

      DATA ICAL/0/

      data ctup /'t','x','y','z','rx','ry','rz','rt','p','rea','ima','roi'
     &            ,'iob','ie','yob','zob','betn','dtom','emod','dmod'
     &            ,'spec','te'/

      save ical

      allocate(reromp(11,ndwsou))

      IF (ICAL.EQ.0.AND.IBUFF.EQ.1)  THEN
c        if (ndwsou.gt.100000) then
c          write(lungfo,*)
c          write(lungfo,*)"      *** Warning in souint_omp: Too many steps in source point for ISPECMODE = 3 ***"
c          write(lungfo,*)"      *** Reduce MYINUM or NLPOI or change ISPECMODE or do not use OMP, i.e. set MTREADS = 0 ***"
c          write(lungfo,*)
c          print*
c          print*,"      *** Warning in souint_omp: Too many steps in source point for ISPECMODE = 3 ***"
c          print*,"      *** Reduce MYINUM or NLPOI or change ISPECMODE or do not use OMP, i.e. set MTREADS = 0 ***"
c          print*
c        endif
        reaima=0.0d0
        IX10=1
        JDX10=NBUFF*NOBSV/10
        JX10=JDX10
        NNBUFF=1
        IF (IWFILINT.EQ.-ISOUR)
     &    CALL hbookm(NIDSOURCE,'RADIATION INTEGRAL$',NTUPP
     &    ,'//WAVE',1024,CTUP)
      ENDIF

      if (isour.ne.isouro) then
        afreq=(0.0d0,0.0d0)
        if (istokes.ne.0) stokes=0.0d0
        specpow=0.0d0
      endif

      IF (JDX10.LT.1) JDX10=1

      IF (ISOUR.EQ.1.AND.IBUFF.EQ.1.AND.NOBSV.GT.1) THEN
        WRITE(6,*)' '
        WRITE(6,*)
     &    ' counting from 1 to 10 for first source to show progress:'
        WRITE(6,*)' '
      ENDIF

C--- LOOP OVER ALL OBSERVATION POINTS

      iobsv=1

C- CALCULATE FREQUENCE INDEPENDENT PARTS OF INTEGRANTES,STORE RESULT IN ARRAYS

      CALL REARG_omp(ISOUR,IOBSV,ndwsou,reromp)  !REAL PARTS OF INTEGRANTS

C--- INTEGRATION FOR ALL FREQUENCES

      CALL ARGSUM_omp(ISOUR,IOBSV,IBUFF,ndwsou,reromp)

      IF(IWFILINT.EQ.ISOUR.AND.IOBSV.EQ.1) CALL WFILINT

      NNBUFF=NNBUFF+1

      IF (ICAL.EQ.0.AND.NNBUFF.EQ.JX10.AND.IX10.LE.10.AND.NOBSV.GT.1) THEN
        JX10=JX10+JDX10
        CALL date_and_time(dtday,dttime,dtzone,idatetime)
        WRITE(6,*)' ',IX10,' ',dttime(1:2),':',dttime(3:4),':',dttime(5:6)
        IX10=IX10+1
      ENDIF

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& SHARED(ical,isour,ibuff,iwfilint,nobsv,nargum,tbuff,afreq)
!$OMP& SHARED(nnbuff,ix10,jx10,jdx10,ndwsou)

!$OMP DO

      DO IOBSV=2,NOBSV-1

C- CALCULATE FREQUENCE INDEPENDENT PARTS OF INTEGRANTES,STORE RESULT IN ARRAYS

        CALL REARG_omp(ISOUR,IOBSV,ndwsou,reromp)  !REAL PARTS OF INTEGRANTS
C        if (iobsv.eq.1) print*,ith,iobsv,reromp(3,1),
C     &    reromp(3,nargum(iobsv,isour))

C--- INTEGRATION FOR ALL FREQUENCES

        CALL ARGSUM_omp(ISOUR,IOBSV,IBUFF,ndwsou,reromp)

        IF(IWFILINT.EQ.ISOUR.AND.IOBSV.EQ.1) CALL WFILINT

        NNBUFF=NNBUFF+1

        IF (ICAL.EQ.0.AND.NNBUFF.EQ.JX10.AND.IX10.LE.10.AND.NOBSV.GT.1) THEN
          JX10=JX10+JDX10
          CALL date_and_time(dtday,dttime,dtzone,idatetime)
          WRITE(6,*)' ',IX10,' ',dttime(1:2),':',dttime(3:4),':',dttime(5:6)
          IX10=IX10+1
        ENDIF

      ENDDO !LOOP OVER ALL OBSERVATION POINTS

!$OMP END DO
!$OMP END PARALLEL

      iobsv=nobsv

C- CALCULATE FREQUENCE INDEPENDENT PARTS OF INTEGRANTES,STORE RESULT IN ARRAYS

      CALL REARG_omp(ISOUR,IOBSV,ndwsou,reromp)  !REAL PARTS OF INTEGRANTS

C--- INTEGRATION FOR ALL FREQUENCES

      CALL ARGSUM_omp(ISOUR,IOBSV,IBUFF,ndwsou,reromp)

      IF(IWFILINT.EQ.ISOUR.AND.IOBSV.EQ.1) CALL WFILINT

      NNBUFF=NNBUFF+1

      IF (ICAL.EQ.0.AND.NNBUFF.EQ.JX10.AND.IX10.LE.10.AND.NOBSV.GT.1) THEN
        JX10=JX10+JDX10
        CALL date_and_time(dtday,dttime,dtzone,idatetime)
        WRITE(6,*)' ',IX10,' ',dttime(1:2),':',dttime(3:4),':',dttime(5:6)
        IX10=IX10+1
      ENDIF

      IF (IBUFF.EQ.NBUFF) THEN
        IF (IWFILINT.EQ.-ISOUR) THEN
          CALL MHROUT(NIDSOURCE,ICYCLE,' ')
          CALL hdeletm(NIDSOURCE)
        ENDIF
        ICAL=1
      ENDIF

      deallocate(reromp)

      RETURN
      END
