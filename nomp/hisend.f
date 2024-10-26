*CMZ :  4.01/00 09/01/2023  20.52.32  by  Michael Scheer
*CMZ :  4.00/14 22/12/2021  13.31.15  by  Michael Scheer
*CMZ :  4.00/13 20/12/2021  16.32.19  by  Michael Scheer
*CMZ :  4.00/11 22/04/2021  18.57.39  by  Michael Scheer
*CMZ :  4.00/06 09/12/2019  14.35.09  by  Michael Scheer
*CMZ :  4.00/04 23/08/2019  15.23.40  by  Michael Scheer
*CMZ :  3.03/00 17/08/2015  10.09.12  by  Michael Scheer
*CMZ :  3.02/05 08/06/2015  13.40.53  by  Michael Scheer
*CMZ :  3.02/04 18/03/2015  09.33.03  by  Michael Scheer
*CMZ :  3.02/00 24/09/2014  11.52.00  by  Michael Scheer
*CMZ :  3.01/07 23/06/2014  13.59.37  by  Michael Scheer
*CMZ :  3.01/05 13/06/2014  11.42.17  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.13.36  by  Michael Scheer
*CMZ :  2.70/06 03/01/2013  13.47.56  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  15.31.59  by  Michael Scheer
*CMZ :  2.70/00 08/11/2012  12.29.49  by  Michael Scheer
*CMZ :  2.68/05 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.68/02 26/06/2012  11.55.39  by  Michael Scheer
*CMZ :  2.67/01 01/03/2012  16.47.26  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  13.32.15  by  Michael Scheer
*CMZ :  2.63/05 03/08/2009  08.58.29  by  Michael Scheer
*CMZ :  2.41/13 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.41/10 14/08/2002  17.34.01  by  Michael Scheer
*CMZ :  1.00/00 24/09/97  17.11.48  by  Michael Scheer
*CMZ : 00.01/09 06/10/95  10.20.31  by  Michael Scheer
*CMZ : 00.01/08 18/07/95  10.46.26  by  Michael Scheer
*CMZ : 00.01/04 30/01/95  16.00.05  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.29  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE HISEND
*KEEP,gplhint.
*KEND.

C--- TERMINATES HBOOK AND WRITE HISTOS TO FILE

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,reargf90u.
      include 'reargf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEEP,afreqf90u.
      include 'afreqf90u.cmn'
*KEEP,amplif90u.
      include 'amplif90u.cmn'
*KEND.

      use bunchmod
      use mhbook_mod

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,optic.
      include 'optic.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,specdip.
      include 'specdip.cmn'
*KEEP,colli.
      include 'colli.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,wusem.
      include 'wusem.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,whbook.
      include 'whbook.cmn'
*KEEP,pawcmn.
*KEEP,phasef90.
      include 'phasef90.cmn'
*KEEP,ntuple2.
      include 'ntuple2.cmn'
*KEEP,datetime.
      include 'datetime.cmn'
*KEEP,wvers.
      include 'wvers.cmn'
*KEND.

      INTEGER ICYCLE,ISCRATCH

      INTEGER I
      integer :: iline=0
      BYTE I1(96)

      CHARACTER(32) LINE_1,LINE_2,LINE_3
      CHARACTER(1) LINE32_1(32),LINE32_2(32),LINE32_3(32)
      CHARACTER(96) C96,EMPTY
      CHARACTER(1) C1(96)

      EQUIVALENCE(I1,C1,C96)
      EQUIVALENCE(LINE_1,LINE32_1)
      EQUIVALENCE(LINE_2,LINE32_2)
      EQUIVALENCE(LINE_3,LINE32_3)

      COMMON/CLINEc/LINE_1,LINE_2,LINE_3

      integer jday,jtim,ius

      DATA EMPTY/' '/

      fillp=-9999.

      call date_and_time(dtday,dttime,dtzone,idatetime)
      read(dtday(3:8),*)jday
      read(dttime(1:6),*)jtim

      fillp(1)=dmyenergy
      fillp(2)=dmycur
      if (ispec.ne.0) then
        fillp(3)=ipin
        fillp(4)=ipincirc
        fillp(5)=pincen(1)
        fillp(6)=pincen(2)
        fillp(7)=pincen(3)
        fillp(8)=pinw
        fillp(9)=pinh
        fillp(10)=pinr
        if (ipin.eq.3) then
          fillp(11)=mpinzorig
          fillp(12)=mpinyorig
        else
          fillp(11)=mpinz
          fillp(12)=mpiny
        endif
        fillp(13)=mpinr
        fillp(14)=mpinphi
        fillp(15)=icbrill
        fillp(16)=obsv(1,icbrill)
        fillp(17)=obsv(2,icbrill)
        fillp(18)=obsv(3,icbrill)
      endif
      fillp(19)=phcenx
      fillp(20)=bsigz(1)
      fillp(21)=bsigy(1)
      fillp(22)=bsigzp(1)
      fillp(23)=bsigyp(1)
      fillp(24)=espread
      if (ispec.ne.0) then
        fillp(25)=ifreq2p
        fillp(26)=nintfreq
        fillp(27)=freqlow
        fillp(28)=freqhig
      endif
      fillp(29)=ispec
      if (ispec.ne.0) then
        fillp(30)=ispecmode
        fillp(31)=ispecdip
        fillp(32)=nlpoi
        fillp(33)=banwid
      endif
      fillp(34)=ibunch
      fillp(35)=nbunch
      fillp(36)=neinbunch
      if (ispec.ne.0) then
        fillp(37)=iampli
      endif
      fillp(38)=ieneloss
      fillp(39)=ifold
      fillp(40)=iefold
      fillp(41)=icode
      fillp(42)=jday
      fillp(43)=jtim
      fillp(44)=wversion
      fillp(45)=istokes
      fillp(46)=nobsvy
      fillp(47)=nobsvz
      fillp(48)=wall(1)
      fillp(49)=wall(2)
      fillp(50)=xabsorb
      fillp(51)=zabsorb(1)
      fillp(52)=zabsorb(2)
      if (ispec.ne.0) then
        fillp(53)=ibrill
      endif
      fillp(54)=ihbeta
      fillp(55)=kampli

      call hfm(nid222,fillp)

      do ius=1,ntupusp
        fillu(ius)=user(ius)
      enddo

      call hfm(nid223,fillu)

      CALL MHROUT(0,ICYCLE,' ')

      if (iroottrees.lt.0) return

      IF (IHOUTP.NE.0) THEN

        CALL hbntm(IDOUTP,'OUTPUT FILE WAVE.OUT',' ')
        CALL hbnamcm(IDOUTP,'LINE',LINE_1
     &    ,'LINE_1:C*32,LINE_2:C*32,LINE_3:C*32')

        REWIND(LUNGFO)

100     C96=EMPTY
        READ (LUNGFO,'(A96)',END=90) C96

C- FILL CWN BUFFER

        DO I=1,32
          LINE32_1(I)=C1(I)
          LINE32_2(I)=C1(I+32)
          LINE32_3(I)=C1(I+64)
        ENDDO

        CALL hfntm(IDOUTP)
        iline=iline+1
        GOTO 100

90      CONTINUE


        DO ISCRATCH=1,ICYCLE-1
          CALL HSCRm(IDOUTP,ISCRATCH,' ')
        ENDDO

        CALL hdeletnomh(IDOUTP)

        rewind(lungfo)

        write(luns_mh(kfile_mh),'(a)')" ! ------------------------------------------------"
        write(luns_mh(kfile_mh),'(a)')"                   16                    4    ! id and kind of histogram of Ntuple"
        write(luns_mh(kfile_mh),'(a)')"          20     ! length of title"
        write(luns_mh(kfile_mh),'(a)')"OUTPUT FILE WAVE.OUT"
        write(luns_mh(kfile_mh),'(a)')"                    6           3          11    ! length of pathname, number of variables, and length of variable names"
        write(luns_mh(kfile_mh),'(a)')"//WAVE"
        write(luns_mh(kfile_mh),'(a)')"LINE_1:C*32"
        write(luns_mh(kfile_mh),'(a)')"LINE_2:C*32"
        write(luns_mh(kfile_mh),'(a)')"LINE_3:C*32"
        write(luns_mh(kfile_mh),*)iline," ! number of entries"

        do i=1,iline
          read (lungfo,'(a96)') c96
          write(luns_mh(kfile_mh),'(a)')c96
        enddo

        close(lungfo)

        open(unit=lungfo,file=filegfo,status='old',access='append')

      ENDIF !IHOUTP

      IF (IHINDEX.NE.0) then
        open(unit=99,file='wave_index.his')
        CALL HOUTPUm(99)
        CALL HLDIRm(' ','IT')
        close(99)
        call houtpum(6)
c obsolete          call system('mv fort.99 wave_index.his')
        print*,' '
        print*,'      Index of histograms written to wave_index.his'
        print*,' '
      endif

      if (iroottrees.ge.0) then
        CALL hrendm('WAVE')
      endif

      CLOSE(LUNHB)

        call mh_end

      RETURN
      END
