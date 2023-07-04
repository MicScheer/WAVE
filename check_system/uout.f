*CMZ :  4.00/17 15/11/2022  10.13.07  by  Michael Scheer
*CMZ :  4.00/15 01/06/2022  16.47.26  by  Michael Scheer
*CMZ :  4.00/11 28/06/2021  10.33.06  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine uout

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.

      use sourcef90
      use observf90
      use afreqf90
      !use waveenv

      implicit none

      complex*16 amp0(3),amp(3),zexp

      double precision dtelec,dtpho,t0,perlen,dph,dobs(3),drn(3),cosang,
     &  dobsn(3),r0(3),r(3),dr(3),v(3),t,dt,dist,obs(3),om,fd,dist0

      double precision :: enemax=0.0d0,s0max=-1.0d30

      integer ifreq,iobsv,nper,i,lunio,iseed
      integer, parameter :: ndimsplit=100
      integer :: nwords=0, ipos(2,ndimsplit),istat=0,mode=0,io,ifr

      character(2048) cline

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,track0.
      include 'track0.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,observ.
      include 'observ.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,berror.
      include 'berror.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,waveenv.
      include 'waveenv.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

      if (iuout.eq.20220601) then
        enemax=spectotmx(3)
        s0max=spectotmx(4)/1.0d6
        open(newunit=lunio,file='uout.in',status='old',iostat=iseed)
        if (iseed.eq.0) then
          read(lunio,*) iseed
        else
          iseed=0
        endif
        close(lunio)
        open(newunit=lunio,file='serie_pherror.out',access='append')
        if (iberror.ne.0) then
          if (nberror.eq.0) then
            write(lunio,*) trim(chwstage)//'-1/2_1_-1/2',iseed,b0error,resrms*360.0d0,
     &        enemax,s0max
          else if (nberror.eq.1) then
            write(lunio,*) trim(chwstage)//'_1_',iseed,b0error,resrms*360.0d0,
     &        enemax,s0max
          else if (nberror.eq.2) then
            write(lunio,*) trim(chwstage)//'-1_1_',iseed,b0error,resrms*360.0d0,
     &        enemax,s0max
          endif
        else
          if (userchar(1).eq.'b0error') then
            write(lunio,*) trim(chwstage)//"_berr",iseed,'-1',user(1),resrms*360.0d0,
     &        enemax,s0max
          else
            write(lunio,*) trim(chwstage)//'_pherr',iseed,'-2',pherror,resrms*360.0d0,
     &        enemax,s0max
          endif
        endif
        flush(lunio)
        close(lunio)
        return
      endif

      if (iuout.eq.20220516) then
c        s0max=-1.0d30
c        do io=1,nobsv
c          do ifr=1,nfreq
c            if (spectot(io+nobsv*(ifr-1)).gt.s0max) then
c              enemax=freq(ifr)
c              s0max=spectot(io+nobsv*(ifr-1))/1.0d6
c            endif
c          enddo
c        enddo
        enemax=spectotmx(3)
        s0max=spectotmx(4)/1.0d6
        open(newunit=lunio,file='uout.in',status='old',iostat=iseed)
        if (iseed.eq.0) then
          read(lunio,*) iseed
        else
          iseed=0
        endif
        close(lunio)
        open(newunit=lunio,file='serie_b0error.out',access='append')
        write(lunio,'(2I10,4(1pe15.7))') icode,iseed,b0error,resrms*360.0d0,
     &    enemax,s0max
        flush(lunio)
        close(lunio)
        return
      endif

      if (iuout.eq.20220517) then
c        s0max=-1.0d30
c        do io=1,nobsv
c          do ifr=1,nfreq
c            if (spectot(io+nobsv*(ifr-1)).gt.s0max) then
c              enemax=freq(ifr)
c              s0max=spectot(io+nobsv*(ifr-1))/1.0d6
c            endif
c          enddo
c        enddo
        enemax=spectotmx(3)
        s0max=spectotmx(4)/1.0d6
        open(newunit=lunio,file='uout.in',status='old',iostat=iseed)
        if (iseed.eq.0) then
          read(lunio,'(a)') cline
          call util_string_split(cline,ndimsplit,nwords,ipos,istat)
        else
          cline='0 unknown'
        endif
        close(lunio)
        open(newunit=lunio,file='serie_amprep.out',access='append')
        write(lunio,'(I10," ",a," ",3(1pe15.7))') icode,trim(cline),
     &    pherror,enemax,s0max
        flush(lunio)
        close(lunio)
        return
      endif

      afreq=(0.0d0,0.0d0)
      dtelec=tftrack-t0track

      r0=[x0,y0,z0]
      dr=[xf0-x0,yf0-y0,zf0-z0]
      drn=dr/norm2(dr)

      perlen=norm2(dr)
      dtpho=perlen/clight1

      nper=100

      do iobsv=1,nobsv
        obs=obsv(1:3,iobsv)
        dist0=norm2(obs-r0)
        do ifreq=1,nfreq
          om=freq(ifreq)/hbarev1
          ifrob=iobsv+nobsv*(ifreq-1)
          amp=(0.0d0,0.0d0)
          amp0(1:3)=dcmplx(reaima(1:3,1,ifrob),reaima(1:3,2,ifrob))
          t=-dt
          do i=1-nper/2,nper-nper/2
            r=r0+i*dr
            dobs=obs-r
            dist=norm2(obs-r)
            dobsn=dobs/dist
            cosang=dot_product(drn,dobsn)
            dt=dtelec-dtpho*cosang
            t=t+dt
            dph=om*t
            zexp=cdexp(dcmplx(0.0d0,dph))
            amp=amp+amp0*zexp*dist0/dist
          enddo
          afreq(:,ifrob)=amp
          fd=sum(real(amp)**2+imag(amp)**2)*specnor
          write(66,*)sngl(obsv(:,iobsv)),sngl(freq(ifreq)),
     &      sngl(real(amp(3))),sngl(imag(amp(3))),sngl(fd)
        enddo !ifreq
      enddo !nobsv


      return
      end
