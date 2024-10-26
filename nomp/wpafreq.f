*CMZ :  4.00/04 02/08/2019  18.44.17  by  Michael Scheer
*CMZ :  3.08/01 31/03/2019  12.56.02  by  Michael Scheer
*CMZ :  3.07/01 29/03/2019  14.26.55  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.66/11 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.66/10 04/05/2010  12.39.03  by  Michael Scheer
*CMZ :  2.66/09 29/04/2010  11.46.31  by  Michael Scheer
*-- Author :    Michael Scheer   17/03/2010
      subroutine wpafreq
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
*KEND.

      complex*16 apolh,apolr,apol45,apoll,ax,ay,az
      double precision aza,aya,phi

      integer ifreq,iobsv,isour

      if (nbunch.eq.1.and.neinbunch.ne.1) then

        do iobsv=1,nobsv
          do ifreq=1,nfreq

            iobfr=iobsv+nobsv*(ifreq-1)

            ax=dcmplx(reaima(1,1,iobfr),reaima(1,2,iobfr))
            ay=dcmplx(reaima(2,1,iobfr),reaima(2,2,iobfr))
            az=dcmplx(reaima(3,1,iobfr),reaima(3,2,iobfr))

            spec(iobfr)=
     &        dreal(
     &        ax*conjg(ax)
     &        +ay*conjg(ay)
     &        +az*conjg(az)
     &        )*specnor

            if (istokes.ne.0) then

              apolh=
     &          ax*conjg(vstokes(1,1))
     &          +ay*conjg(vstokes(1,2))
     &          +az*conjg(vstokes(1,3))

              apolr=
     &          ax*conjg(vstokes(2,1))
     &          +ay*conjg(vstokes(2,2))
     &          +az*conjg(vstokes(2,3))

              apoll=
     &          ax*conjg(vstokes(3,1))
     &          +ay*conjg(vstokes(3,2))
     &          +az*conjg(vstokes(3,3))

              apol45=
     &          ax*conjg(vstokes(4,1))
     &          +ay*conjg(vstokes(4,2))
     &          +az*conjg(vstokes(4,3))

              stokes(1,iobfr)=(apolr*conjg(apolr)+apoll*conjg(apoll))
              stokes(2,iobfr)=(-stokes(1,iobfr)+2.*apolh*conjg(apolh))
              stokes(3,iobfr)=(2.*apol45*conjg(apol45)-stokes(1,iobfr))
              stokes(4,iobfr)=(apolr*conjg(apolr)-apoll*conjg(apoll))

              stokes(1:4,iobfr)=stokes(1:4,iobfr)*specnor

            endif !istokes

          enddo !ifreq
        enddo !iobsv

      else if (neinbunch.eq.1) then

        isour=1
        reaima=0.0d0

        if (specnor*bunnor.eq.0.0d0) then
          print*,'*** Warning in WPAFREQ: Zero normalization on cluster files'
          return
        endif

        wpspecnor=specnor*bunnor

        do iobsv=1,nobsv
          do ifreq=1,nfreq

            iliobfr=isour+nsource*(iobsv-1+nobsv*(ifreq-1))
            iobfr=iobsv+nobsv*(ifreq-1)

            if (istokes.eq.0) then

              reaima(3,1,iobfr)=sqrt(spec(iliobfr)/wpspecnor)

            else !istokes

              aza=sqrt((stokes(1,iobfr)+stokes(2,iobfr))/2.0d0/wpspecnor)
              aya=sqrt((stokes(1,iobfr)-stokes(2,iobfr))/2.0d0/wpspecnor)

              phi=0.0d0
              if (aza*aya.ne.0.0d0) then
                if (abs(stokes(3,iobfr)).gt.abs(stokes(2,iobfr))) then
                  phi=asin(stokes(3,iobfr))/2.0d0/aza/aya
                else
                  phi=acos(stokes(2,iobfr))/2.0d0/aza/aya
                endif
              endif

              reaima(2,1,iobfr)=aya*cos(phi)
              reaima(2,2,iobfr)=aya*sin(phi)
              reaima(3,1,iobfr)=aza
              reaima(3,2,iobfr)=0.0d0

            endif !istokes

          enddo !ifreq
        enddo !iobsv

      endif !nbunch

      return
      end
