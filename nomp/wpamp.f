*CMZ :          02/05/2024  11.53.53  by  Michael Scheer
*CMZ :  4.01/03 12/06/2023  11.10.19  by  Michael Scheer
*CMZ :  4.00/17 15/11/2022  10.11.12  by  Michael Scheer
*CMZ :  4.00/15 14/03/2022  09.02.26  by  Michael Scheer
*CMZ :  4.00/06 05/12/2019  13.19.05  by  Michael Scheer
*CMZ :  4.00/04 05/08/2019  14.06.21  by  Michael Scheer
*CMZ :  4.00/02 12/04/2019  14.50.24  by  Michael Scheer
*CMZ :  3.08/01 02/04/2019  12.40.27  by  Michael Scheer
*CMZ :  3.07/01 29/03/2019  15.44.29  by  Michael Scheer
*CMZ :  3.03/02 02/03/2016  10.29.09  by  Michael Scheer
*CMZ :  3.02/00 28/08/2014  08.52.10  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.66/11 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.66/10 04/05/2010  11.49.38  by  Michael Scheer
*CMZ :  2.66/09 29/04/2010  11.46.31  by  Michael Scheer
*-- Author :    Michael Scheer   17/03/2010
      subroutine wpamp
*KEEP,gplhint.
*KEND.

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
      use clustermod
      !use waveenv

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
*KEEP,primkin.
      include 'primkin.cmn'
*KEEP,strings.
      include 'strings.cmn'
*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.

      double precision rmax2,rabs2,rea1,rea2,wpspecnoro,specnoro,bunnoro,
     &  s1,s2,s3,s4,pow,powt,buno,dw

      real*8 corrins !NIDBUNCH

      integer kfreq,iobsv,jcode,jfreq,jobsv,lastch,ifrmx,iobmx,njobs,ijob,
     &  ifirstch,jliobfr,islash,i,iexist,lun,lun98

      character(2048) chfile
      character(256) cstat,cjobnum
      character cslash

      equivalence (cslash,islash)

      wpspecnor=0.0d0


      if (nsource.ne.1) then
        write(6,*)
     &    ' '
        write(6,*)
     &    '*** Error in WPAMP: More then one source not allowd for WPAMP ***'
        write(6,*)
     &    '          *** Program WAVE aborted ***'
        write(lungfo,*)
     &    ' '
        write(lungfo,*)
     &    '*** Error in WPAMP: More then one source not allowd for WPAMP ***'
        write(lungfo,*)
     &    '          *** Program WAVE aborted ***'
        stop
      endif

      if (icluster.gt.0) then

        write(lungfo,*)
        write(lungfo,*)'     WPAMP: Processing spawned runs'
        write(lungfo,*)'     (If some information is missing, please check wave.out of spawned runs, i.e. ....stage.1/wave.out etc.)'
        write(lungfo,*)

        call util_string_trim(trim(chwavedir),ifirstch,lastch)

        njobs=0
        reaima=0.0d0
        spec=0.0d0
        specpow=0.0d0

        corrins=dble(nwinstances)/dble(nwgood)
        if (corrins.ne.1.0d0) then
          print*,""
          print*,"*** Warning in WPAMP: NOT ALL INSTANCES HAVE FINISHED CORRECTLY: BE CAREFUL, ESPECIALLY WITH NORMALIZATION!!"
          print*,""
          write(lungfo,*)""
          write(lungfo,*)"*** Warning in WPAMP: NOT ALL INSTANCES HAVE FINISHED CORRECTLY: BE CAREFUL, ESPECIALLY WITH NORMALIZATION!!"
          write(lungfo,*)""
        endif

        if (istokes.ne.0) stokes=0.0d0

        do njobs=1,nwinstances

          if (iwstat(njobs).ne.0) cycle

          write(cjobnum,*)njobs

          chfile=
     &      chwavedir(ifirstch:lastch)//
     &      '/.stage.'//trim(adjustl(cjobnum))//
     &      '/wave_cluster.dat'

          open(newunit=lun,file=chfile(1:len_trim(chfile)),status='unknown',
     &      err=99)

          read(lun,*)jcode

          if (ibunch.ne.0) then

            if (nbunch.eq.1.and.neinbunch.ne.1) then

              do kfreq=1,nfreq
                do iobsv=1,nobsv

                  iobfr=iobsv+nobsv*(kfreq-1)

                  read(lun,*)
     &              jfreq,jobsv

                  if (jfreq.ne.kfreq.or.jobsv.ne.iobsv) then
                    write(lungfo,*)
     &                '*** Error in WPAMP: Bad photon enery or observation point ***'
                    write(6,*)
     &                '*** Error in WPAMP: Bad photon enery or observation point ***'
                    stop '*** Program WAVE aborted ***'
                  endif

                  do i=1,10
                    read(lun,*)rea1,rea2
                    reaima(i,1,iobfr)=rea1+reaima(i,1,iobfr)
                    reaima(i,2,iobfr)=rea2+reaima(i,2,iobfr)
                  enddo

                enddo !iobsv
              enddo !kfreq

              read(lun,*) wpspecnor

              if (njobs.eq.1) then
                read(lun,*)ecsour(1:4,nsource)
                read(lun,*)ecmax(nsource),speccut
              endif

              read(lun,*)pow
              specpow(nsource)=specpow(nsource)+pow*corrins

              close(lun)

              dw=(wpspecnoro-wpspecnor)/(wpspecnor+wpspecnoro)

              if (njobs.gt.1.and.Abs(dw).gt.1.0e-12) then
                print*,
     &            '*** Warning in WPAMP: Different normalizations on files ***'
                print*,"Deviation:",sngl(dw)
              else
                wpspecnoro=wpspecnor
              endif

            else if (nbunch.gt.1) then

              do kfreq=1,nfreq
                do iobsv=1,nobsv

                  iliobfr=nsource+nsource*(iobsv-1+nobsv*(kfreq-1))
                  read(lun,*)jliobfr

                  if (jliobfr.ne.iliobfr) then
                    write(lungfo,*)
     &                '*** Error in WPAMP: Bad photon enery or observation point ***'
                    write(6,*)
     &                '*** Error in WPAMP: Bad photon enery or observation point ***'
                    stop '*** Program WAVE aborted ***'
                  endif

                  read(lun,*)rea1
                  spec(iliobfr)=spec(iliobfr)+rea1*corrins

                  if (istokes.ne.0) then
                    iobfr=iobsv+nobsv*(kfreq-1)
                    read(lun,*)s1,s2,s3,s4
                    stokes(1,iobfr)=stokes(1,iobfr)+s1*corrins
                    stokes(2,iobfr)=stokes(2,iobfr)+s2*corrins
                    stokes(3,iobfr)=stokes(3,iobfr)+s3*corrins
                    stokes(4,iobfr)=stokes(4,iobfr)+s4*corrins
                  endif

                enddo !iobsv
              enddo !kfreq

              do kfreq=1,nfreq
                do iobsv=1,nobsv
                  iobfr=iobsv+nobsv*(kfreq-1)
                  do i=1,10
                    read(lun,*)rea1,rea2
                    reaima(i,1,iobfr)=rea1+reaima(i,1,iobfr)
                    reaima(i,2,iobfr)=rea2+reaima(i,2,iobfr)
                  enddo
                enddo !iobsv
              enddo !kfreq

              read(lun,*)specnor,bunnor

              if (njobs.eq.1) then
                read(lun,*)ecsour(1:4,nsource)
                read(lun,*)ecmax(nsource),speccut
              endif

              read(lun,*)pow
              specpow(nsource)=specpow(nsource)+pow

              close(lun)

              if (njobs.gt.1.and.(specnor.ne.specnoro.or.bunnor.ne.bunnoro))
     &          stop
     &          '*** Error in WPAMP: Different normalizations on files ***'

              specnoro=specnor
              bunnoro=bunnor

            endif !neinbunch

          endif !(ibunch.ne.0) then

        enddo !njobs

        write(lungfo,*)

        if (nbunch.eq.1.and.neinbunch.ne.1) then
          reaima=reaima/sqrt(dble(nwgood))
        else if (nbunch.gt.1) then
          spec=spec/nwgood
          if (istokes.ne.0) stokes=stokes/nwgood
        endif

c 2.5.2024        call wpafreq

      else if (icluster.lt.0) then

        print*,""
        print*," Writing wave_cluster.dat"
        print*,""

        open(newunit=lun,file='wave_cluster.dat',status='unknown')

        write(lun,*)icode

        if (iclubun.eq.1.and.neinbunch.ne.1) then

          do kfreq=1,nfreq
            do iobsv=1,nobsv

              iobfr=iobsv+nobsv*(kfreq-1)

              write(lun,*)
     &          kfreq,iobsv

              do i=1,10
                write(lun,*) reaima(i,1,iobfr),reaima(i,2,iobfr)
              enddo

              rabs2=
     &          reaima(1,1,iobfr)**2+
     &          reaima(1,2,iobfr)**2+
     &          reaima(2,1,iobfr)**2+
     &          reaima(2,2,iobfr)**2+
     &          reaima(3,1,iobfr)**2+
     &          reaima(3,2,iobfr)**2

              if(rabs2.gt.rmax2) then
                rmax2=rabs2
                ifrmx=kfreq
                iobmx=iobsv
              endif

            enddo !iobsv
          enddo !kfreq

          if (rmax2.gt.0.0d0) then

            wpspecnor=spec(iobmx+nobsv*(ifrmx-1))/rmax2
            write(lun,*) wpspecnor

          else
            write(lun,*)0.0d0
          endif

        else if (iclubun.gt.1) then
          do kfreq=1,nfreq
            do iobsv=1,nobsv

              iliobfr=nsource+nsource*(iobsv-1+nobsv*(kfreq-1))
              write(lun,*)iliobfr
              write(lun,*)spec(iliobfr)

              if (istokes.ne.0) then
                iobfr=iobsv+nobsv*(kfreq-1)
                write(lun,*)stokes(1,iobfr),stokes(2,iobfr),
     &            stokes(3,iobfr),stokes(4,iobfr)
              endif

            enddo !iobsv
          enddo !kfreq

          do kfreq=1,nfreq
            do iobsv=1,nobsv
                iobfr=iobsv+nobsv*(kfreq-1)
                do i=1,10
                  write(lun,*)reaima(i,1:2,iobfr)
                enddo
            enddo !iobsv
          enddo !kfreq

          write(lun,*) specnor,bunnor

        endif !(nbunch.eq.1) then

        write(lun,*)ecsour(1:4,nsource)
        write(lun,*)ecmax(nsource),speccut
        write(lun,*)specpow(nsource)

        flush(lun)
        close(lun)

      endif !icluster

      return

99    write(lungfo,*)
     &  '*** Error in WPAMP: File not found'
      write(6,*)
     &  '*** Error in WPAMP: File not found'
      stop '*** Program WAVE aborted ***'

      end
