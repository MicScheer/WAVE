*CMZ :  2.66/14 20/08/2010  13.13.20  by  Michael Scheer
*CMZ :  2.66/09 22/03/2010  08.59.02  by  Michael Scheer
*CMZ :  2.66/03 29/10/2009  11.56.38  by  Michael Scheer
*CMZ :  2.66/01 15/10/2009  11.55.22  by  Michael Scheer
*CMZ :  2.66/00 13/10/2009  11.57.12  by  Michael Scheer
*CMZ :  2.65/03 30/09/2009  14.52.57  by  Michael Scheer
*-- Author :    Michael Scheer   30/09/2009
      program truebrill

      implicit none

      double precision fieldpar1min,fieldpar1max,dfieldpar1,fieldpar1
     &  ,fieldpar2min,fieldpar2max,dfieldpar2,fieldpar2,fluxmax
c     &  ,ahwpol
      double precision perlenmin,perlenmax,dperlen,perlen
     &  ,waveperlen,wavedevlen
      double precision wavefieldpar1,wavefieldpar2,freqmax,brillmax,
     &  park,wlen1,eharm1,dmygamma,dmyenergy,freqlow,freqhig,dfreq,
     &  sigrc,sigrpc,defl,dlamb

      integer
     &  nfieldpar1,ifieldpar1
     &  ,nfieldpar2,ifieldpar2
     &  ,nperlen,luntbout,luntbnam
     &  ,iperlen,
     &  lchprogmaster,lchprogwork,
     &  lchwavein,lchwaveout,
     &  kchfieldpar1,kchfieldpar2,lchfieldpar1,lchfieldpar2
     &  ,kchperlen,lchperlen
     &  ,kchwaveperlen,lchwaveperlen
c     &  ,kchahwpol,lchahwpol
     &  ,kchnperwave,lchnperwave
     &  ,kchwaveenergy,lchwaveenergy
     &  ,kchwavefieldpar1,lchwavefieldpar1
     &  ,kchwavefieldpar2,lchwavefieldpar2,
     &  kchfreqlow,lchfreqlow,
     &  kchfreqhig,lchfreqhig,
     &  kchsigrc,lchsigrc,
     &  kchsigrpc,lchsigrpc,
     &  iharmmin,iharmmax,idharm,iharm,nperwave,
     &  ifold,iefold

      integer util_igetlastchar
      external function util_igetlastchar

      character(512) chwavein,chprogmaster,chwaveout,chprogwork,
     &  chtbout,chtbnam,command,chprog,chprogout

      character(12) chfieldpar1,chfieldpar2
     &  ,chperlen
     &  ,chwaveperlen
c     &  ,chahwpol
     &  ,chnperwave,chsigrc,chsigrpc
     &  ,chwaveenergy,chwavefieldpar1,chwavefieldpar2,chfreqlow,chfreqhig
      character c1

      logical lexist

      namelist /truebrilln/
     &  wavedevlen,nperwave,
     &  fieldpar1min,fieldpar1max,dfieldpar1,
     &  fieldpar2min,fieldpar2max,dfieldpar2
     &  ,perlenmin,perlenmax,dperlen
     &  ,iharmmin,iharmmax,idharm
     &  ,chprog,chprogmaster,chprogwork,chwavein,chprogout
     &  ,dmyenergy,sigrc,sigrpc,
     &  ifold,iefold

*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,phycon1.
      include 'phycon1.cmn'
*KEND.

      luntbnam=15
      luntbout=16
      chtbnam='truebrill.nam'
      chtbout='truebrill.out'

      chwaveout='wave.in'

      open(unit=luntbnam,file=chtbnam,status='old')
      read(luntbnam,truebrilln)
      close(luntbnam)

      dmygamma=dmyenergy/emassg1

      lchprogmaster=util_igetlastchar(1,256,chprogmaster,c1)
      lchprogwork=util_igetlastchar(1,256,chprogwork,c1)

      lchwavein=util_igetlastchar(1,256,chwavein,c1)
      lchwaveout=util_igetlastchar(1,256,chwaveout,c1)

      if (dfieldpar1.eq.0.0d0) dfieldpar1=1.0d0
      if (dfieldpar2.eq.0.0d0) dfieldpar2=1.0d0

      nfieldpar1=nint((fieldpar1max-fieldpar1min)/dfieldpar1)+1
      nfieldpar2=nint((fieldpar2max-fieldpar2min)/dfieldpar2)+1

      nperlen=nint((perlenmax-perlenmin)/dperlen)+1

      call system('touch truebrill.cont')
      open(unit=luntbout,file=chtbout,status='unknown')

      do iperlen=1,nperlen

        perlen=perlenmin+dperlen*(iperlen-1)

        if (perlen.lt.0.0d0) then
          stop
     &'*** Error in TRUEBRILL:  Negative periodlength!'
        endif

        if (wavedevlen.lt.0.0d0.and.nperwave.gt.0) then
          waveperlen=perlen/1000.0d0
          wavedevlen=-nperwave*waveperlen
        else if (nperwave.lt.0.and.wavedevlen.gt.0.0d0) then
          waveperlen=perlen/1000.0d0
          nperwave=-nint(wavedevlen/waveperlen)
        else
          stop
     &'*** Error in TRUEBRILL: Either WAVEDEVLEN or nperwave must be negative'
        endif

        write(chwaveperlen,'(g12.5)')waveperlen
        lchwaveperlen=util_igetlastchar(1,12,chwaveperlen,c1)
        call util_string_igetfirstchar(kchwaveperlen,chwaveperlen,c1)

        write(chperlen,'(g12.5)')perlen
        print*,'Perlen: ',chperlen
        lchperlen=util_igetlastchar(1,12,chperlen,c1)
        call util_string_igetfirstchar(kchperlen,chperlen,c1)

        write(chnperwave,*)abs(nperwave)
        lchnperwave=util_igetlastchar(1,12,chnperwave,c1)
        call util_string_igetfirstchar(kchnperwave,chnperwave,c1)

        do ifieldpar1=1,nfieldpar1
          do ifieldpar2=1,nfieldpar2

          fieldpar1=fieldpar1min+dfieldpar1*(ifieldpar1-1)
          write(chfieldpar1,'(g12.5)')fieldpar1
          print*,'fieldpar1: ',chfieldpar1
          lchfieldpar1=util_igetlastchar(1,12,chfieldpar1,c1)
          call util_string_igetfirstchar(kchfieldpar1,chfieldpar1,c1)

          fieldpar2=fieldpar2min+dfieldpar2*(ifieldpar2-1)
          write(chfieldpar2,'(g12.5)')fieldpar2
          print*,'fieldpar2: ',chfieldpar2
          lchfieldpar2=util_igetlastchar(1,12,chfieldpar2,c1)
          call util_string_igetfirstchar(kchfieldpar2,chfieldpar2,c1)

          do iharm=iharmmin,iharmmax,idharm

            inquire(file='truebrill.cont',exist=lexist)
            if (lexist.eq..false.) then
              stop '*** Stopping due to missing file truebrill.cont'
            endif

            if (iharm.eq.iharmmin) then
              if (chprog.ne.'') then

                if (chprog.eq.'qnk') then
                  wavefieldpar1=fieldpar1
                  wavefieldpar2=fieldpar2
                  dlamb=perlen/10.0d0
                  DEFL=ECHARGE1*wavefieldpar1*
     &              DLAMB/100.0d0/(2.*PI1*EMASSKG1*CLIGHT1)
                  fieldpar1=defl
                else !qnk

                  command='sed ' //
     & '-e s/truebrillfieldpar1/' // chfieldpar1(kchfieldpar1:lchfieldpar1)// '/ ' //
     & '-e s/truebrillfieldpar2/' // chfieldpar2(kchfieldpar2:lchfieldpar2)// '/ ' //
     & '-e s/truebrillperlen/' // chperlen(kchfieldpar1:lchfieldpar1)// '/ ' //
     & chprogmaster(1:lchprogmaster) // ' > ' //
     & chprogwork(1:lchprogwork)

                  print*,command
                  call system(command)

                  call system(chprog)

                  open(unit=99,file=chprogout,status='old')
                  read(99,*)wavefieldpar1,wavefieldpar2
                  close(99)

                endif !(chprog.eq.'qnk') then

              else !(chprog.ne.'') then
                wavefieldpar1=fieldpar1
                wavefieldpar2=fieldpar2
              endif !(chprog.ne.'') then

            endif !(iharmmin) then

c            ahwpol=nint(1.5/waveperlen)*2+1

            park=echarge1*dabs(wavefieldpar1)*waveperlen/
     &        (2.0d0*pi1*emasskg1*clight1)
            wlen1=(1+park**2/2.0d0)/2.0d0/dmygamma**2*waveperlen*1.0d9
            eharm1=wtoe1/wlen1

c            write(chahwpol,'(g12.5)')ahwpol
c            lchahwpol=util_igetlastchar(1,12,chahwpol,c1)
c            call util_string_igetfirstchar(kchahwpol,chahwpol,c1)

            dfreq=eharm1/iharm/max(abs(nperwave),1)
            freqlow=iharm*(eharm1-dfreq)
            freqhig=iharm*(eharm1+dfreq)

            write(chwaveenergy,'(g12.5)')dmyenergy
            print*,'waveenergy: ',chwaveenergy
            lchwaveenergy=util_igetlastchar(1,12,chwaveenergy,c1)
            call util_string_igetfirstchar(kchwaveenergy,chwaveenergy,c1)

            write(chwavefieldpar1,'(g12.5)')wavefieldpar1
            print*,'wavefieldpar1: ',chwavefieldpar1
            lchwavefieldpar1=util_igetlastchar(1,12,chwavefieldpar1,c1)
            call util_string_igetfirstchar(kchwavefieldpar1,chwavefieldpar1,c1)

            write(chwavefieldpar2,'(g12.5)')wavefieldpar2
            print*,'wavefieldpar2: ',chwavefieldpar2
            lchwavefieldpar2=util_igetlastchar(1,12,chwavefieldpar2,c1)
            call util_string_igetfirstchar(kchwavefieldpar2,chwavefieldpar2,c1)

            write(chfreqlow,'(g12.5)')freqlow
            lchfreqlow=util_igetlastchar(1,12,chfreqlow,c1)
            call util_string_igetfirstchar(kchfreqlow,chfreqlow,c1)

            write(chfreqhig,'(g12.5)')freqhig
            lchfreqhig=util_igetlastchar(1,12,chfreqhig,c1)
            call util_string_igetfirstchar(kchfreqhig,chfreqhig,c1)

            write(chsigrc,'(g12.5)')sigrc
            lchsigrc=util_igetlastchar(1,12,chsigrc,c1)
            call util_string_igetfirstchar(kchsigrc,chsigrc,c1)

            write(chsigrpc,'(g12.5)')sigrpc
            lchsigrpc=util_igetlastchar(1,12,chsigrpc,c1)
            call util_string_igetfirstchar(kchsigrpc,chsigrpc,c1)

            command='sed ' //
     &        '-e s/truebrillenergy/' // chwaveenergy(kchwaveenergy:lchwaveenergy)// '/ ' //
     &        '-e s/truebrillfieldpar1/' // chwavefieldpar1(kchwavefieldpar1:lchwavefieldpar1)// '/ ' //
     &        '-e s/truebrillfieldpar2/' // chwavefieldpar2(kchwavefieldpar2:lchwavefieldpar2)// '/ ' //
     &        '-e s/truebrillperlen/' // chwaveperlen(kchwaveperlen:lchwaveperlen)// '/ ' //
     &        '-e s/truebrillnper/' // chnperwave(kchnperwave:lchnperwave)// '/ ' //
     &        '-e s/truebrillfreqlow/' // chfreqlow(kchfreqlow:lchfreqlow)// '/ ' //
     &        '-e s/truebrillfreqhig/' // chfreqhig(kchfreqhig:lchfreqhig)// '/ ' //
     &        '-e s/truebrillsigrc/' // chsigrc(kchsigrc:lchsigrc)// '/ ' //
     &        '-e s/truebrillsigrpc/' // chsigrpc(kchsigrpc:lchsigrpc)// '/ ' //
c     &        '-e s/tbahwpol/' // chahwpol(kchahwpol:lchahwpol)// '/ ' //
     &        chwavein(1:lchwavein) // ' > ' //
     &        chwaveout(1:lchwaveout)

            print*,command
            call system(command)

            call system('/scheer/wav/truebrill/wave_truebrill.sh')

            if (ifold.eq.0.and.iefold.eq.0) then
              open(unit=99,file='brill_brilliance.dat'
     &          ,status='old')
            else if (ifold.ne.0.and.iefold.eq.0) then
              open(unit=99,file='brill_brilliance_f.dat'
     &          ,status='old')
            else if (ifold.eq.0.and.iefold.ne.0) then
              open(unit=99,file='brill_brilliance_e.dat'
     &          ,status='old')
            else if (ifold.ne.0.and.iefold.ne.0) then
              open(unit=99,file='brill_brilliance_ef.dat'
     &          ,status='old')
            endif
            read(99,*)freqmax,brillmax
            close(99)

            if (ifold.eq.0.and.iefold.eq.0) then
              open(unit=99,file='brill_flux.dat'
     &          ,status='old')
            else if (ifold.ne.0.and.iefold.eq.0) then
              open(unit=99,file='brill_flux_f.dat'
     &          ,status='old')
            else if (ifold.eq.0.and.iefold.ne.0) then
              open(unit=99,file='brill_flux_e.dat'
     &          ,status='old')
            else if (ifold.ne.0.and.iefold.ne.0) then
              open(unit=99,file='brill_flux_ef.dat'
     &          ,status='old')
            endif
            read(99,*)freqmax,fluxmax
            close(99)

            if (chprog.eq.'qnk') then
              fieldpar2=fluxmax/pi1/alpha1/nperwave*echarge1
            endif

            write(luntbout,'(4g12.5,I4,4g12.5)')
     &        ,fieldpar1,fieldpar2
     &        ,wavefieldpar1,wavefieldpar2,iharm,waveperlen,
     &        freqmax,brillmax,fluxmax

          enddo !harm

        enddo !fieldpar2
        enddo !fieldpar1

      enddo !perlen

      close(luntbout)

      stop
      end
