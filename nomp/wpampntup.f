*CMZ :          02/05/2024  12.36.06  by  Michael Scheer
*CMZ :  4.01/05 20/04/2024  09.47.23  by  Michael Scheer
*CMZ :  4.00/17 15/11/2022  10.12.04  by  Michael Scheer
*CMZ :  4.00/15 14/03/2022  09.02.26  by  Michael Scheer
*CMZ :  4.00/04 05/08/2019  15.46.18  by  Michael Scheer
*CMZ :  3.08/01 02/04/2019  12.38.49  by  Michael Scheer
*CMZ :  3.07/01 29/03/2019  12.42.19  by  Michael Scheer
*CMZ :  3.03/02 01/03/2016  20.59.13  by  Michael Scheer
*CMZ :  3.02/00 28/08/2014  08.52.10  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.66/11 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.66/10 04/05/2010  11.49.38  by  Michael Scheer
*CMZ :  2.66/09 29/04/2010  11.46.31  by  Michael Scheer
*-- Author :    Michael Scheer   17/03/2010
      subroutine wpampntup
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
*KEEP,photon.
      include 'photon.cmn'
*KEEP,strings.
      include 'strings.cmn'
*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.

      double precision rmax2,rabs2,rea1,rea2,wpspecnoro,specnoro,bunnoro,
     &  s1,s2,s3,s4

      real*8 count,fillb(100) !NIDBUNCH

      integer ifreq,iobsv,jcode,jfreq,jobsv,lastch,ifrmx,iobmx,njobs,ijob,
     &  ifirstch,isour,jliobfr,islash,i,iexist,lun99,lun98,iclo

      character(2048) chfile
      character(256) cstat,cjobnum
      character cslash

      equivalence (cslash,islash)

      iclo=icluster

      if (icluster.gt.0) then

        icluster=9999

        call util_string_trim(trim(chwavedir),ifirstch,lastch)

        do njobs=1,nwinstances

          if (iwstat(njobs).ne.0) cycle

          write(cjobnum,*)njobs

          chfile=
     &      chwavedir(ifirstch:lastch)//
     &      '/.stage.'//trim(adjustl(cjobnum))//
     &      '/wave_cluster.dat'

          open(newunit=lun99,file=chfile(1:len_trim(chfile)),status='old',
     &      err=99)

          read(lun99,*)jcode

          if (ibunch.ne.0) then

            if (ihphotons.ne.0) then

              chfile=
     &          chwavedir(ifirstch:lastch)//
     &          '/.stage.'//trim(adjustl(cjobnum))//
     &          '/msh_hbook_ntuple_7777.scr'

              call fexist(chfile,iexist)

              if (iexist.eq.1) then
                open(newunit=lun99,file=chfile)
13611           read(lun99,*,end=93611) fillb(1:14)
                call hfm(7777,fillb)
                goto 13611
93611           close(lun99)
              endif

            endif

            if (ihfreq.ne.0) then

              chfile=
     &          chwavedir(ifirstch:lastch)//
     &          '/.stage.'//trim(adjustl(cjobnum))//
     &          '/msh_hbook_ntuple_3601.scr'

              call fexist(chfile,iexist)

              if (iexist.eq.1) then
                open(newunit=lun99,file=chfile)
1361            read(lun99,*,end=9361) fillb(1:36)
                call hfm(3601,fillb)
                goto 1361
9361            close(lun99)
              endif

              chfile=
     &          chwavedir(ifirstch:lastch)//
     &          '/.stage.'//trim(adjustl(cjobnum))//
     &          '/msh_hbook_ntuple_3600.scr'

              call fexist(chfile,iexist)
              if (iexist.eq.1) then
                open(newunit=lun99,file=chfile)
136             read(lun99,*,end=936) fillb(1:2)
                call hfm(3600,fillb)
                goto 136
936             close(lun99)
              endif

              chfile=
     &          chwavedir(ifirstch:lastch)//
     &          '/.stage.'//trim(adjustl(cjobnum))//
     &          '/msh_hbook_ntuple_3700.scr'

              call fexist(chfile,iexist)
              if (iexist.eq.1) then
                open(newunit=lun99,file=chfile)
137             read(lun99,*,end=937) fillb(1:34)
                call hfm(3700,fillb)
                goto 137
937             close(lun99)
              endif

              if (istokes.ne.0) then

                chfile=
     &            chwavedir(ifirstch:lastch)//
     &            '/.stage.'//trim(adjustl(cjobnum))//
     &            '/msh_hbook_ntuple_4600.scr'

                call fexist(chfile,iexist)
                if (iexist.eq.1) then
                  open(newunit=lun99,file=chfile)
146               read(lun99,*,end=946) fillb(1:5)
                  call hfm(4600,fillb)
                  goto 146
946               close(lun99)
                endif

                chfile=
     &            chwavedir(ifirstch:lastch)//
     &            '/.stage.'//trim(adjustl(cjobnum))//
     &            '/msh_hbook_ntuple_4700.scr'

                call fexist(chfile,iexist)
                if (iexist.eq.1) then
                  open(newunit=lun99,file=chfile)
147               read(lun99,*,end=947) fillb(1:12)
                  call hfm(4700,fillb)
                  goto 147
947               close(lun99)
                endif

              endif !istokes

            endif !ihfreq

            if (ihbunch.ne.0) then

              chfile=
     &          chwavedir(ifirstch:lastch)//
     &          '/.stage.'//trim(adjustl(cjobnum))//
     &          '/msh_hbook_ntuple_30.scr'

              open(newunit=lun99,file=chfile)
              do i=1,nbunch*neinbunch/nwgood
                do iobsv=1,nobsv
                  do ifreq=1,nfreq
                    read(lun99,*,end=991) fillb(1:41)
                    if (iobsv.eq.1.and.ifreq.eq.1) count=count+1.0d0
                    if (nbunch.eq.1) then
                      fillb(2)=count
                    else
                      fillb(1)=count
                    endif
                    fillb(3)=count
                    call hfm(nidbunch,fillb)
                  enddo
                enddo
              enddo
991           close(lun99)

            endif !(ihbunch.ne.0) then

          endif !(ibunch.ne.0) then

        enddo !instances

      endif !icluster

      icluster=iclo
      return

99    write(lungfo,*)
     &    '*** Error in WPAMPNTUP: File not found'
      write(6,*)
     &  '*** Error in WPAMPNTUP: File not found'
      stop '*** Program WAVE aborted ***'

      close(lun98)
      end
