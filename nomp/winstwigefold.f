*CMZ :          04/11/2025  14.11.47  by  Michael Scheer
*CMZ :  4.02/00 19/09/2025  13.23.19  by  Michael Scheer
*CMZ :  4.01/03 12/06/2023  11.06.51  by  Michael Scheer
*CMZ :  4.01/00 05/12/2022  09.54.57  by  Michael Scheer
*CMZ :  4.00/17 15/11/2022  10.06.37  by  Michael Scheer
*CMZ :  4.00/15 01/06/2022  16.29.17  by  Michael Scheer
*CMZ :  4.00/13 03/12/2021  16.25.25  by  Michael Scheer
*CMZ :  4.00/07 18/05/2020  09.47.26  by  Michael Scheer
*CMZ :  4.00/06 05/12/2019  13.13.01  by  Michael Scheer
*CMZ :  4.00/05 30/11/2019  15.57.36  by  Michael Scheer
*CMZ :  4.00/04 25/11/2019  15.41.55  by  Michael Scheer
*CMZ :  4.00/03 09/05/2019  10.59.58  by  Michael Scheer
*CMZ :  4.00/02 12/04/2019  15.05.49  by  Michael Scheer
*CMZ :  4.00/01 11/04/2019  14.40.23  by  Michael Scheer
*CMZ :  4.00/00 04/04/2019  12.29.10  by  Michael Scheer
*CMZ :  3.08/01 03/04/2019  11.56.15  by  Michael Scheer
*CMZ :  3.07/01 29/03/2019  14.35.23  by  Michael Scheer
*-- Author :    Michael Scheer   27/03/2019
      subroutine winstwigefold(kwigerr)

      use omp_lib
      use clustermod
      use bunchmod

      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,berror.
      include 'berror.cmn'
*KEEP,ampli.
      include 'ampli.cmn'
*KEEP,halbach.
      include 'halbach.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEEP,ellip.
      include 'ellip.cmn'
*KEEP,photon.
      include 'photon.cmn'
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEEP,phasef90.
      include 'phasef90.cmn'
*KEEP,random.
      include 'random.cmn'
*KEEP,waveenv.
      include 'waveenv.cmn'
*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEEP,wvers.
      include 'wvers.cmn'
*KEND.

      complex*16, dimension (:,:,:,:), allocatable :: esourz,esoury

      double precision, dimension (:,:,:,:,:,:,:), allocatable :: wigefold
      double precision, dimension (:,:), allocatable :: fdwigzy,fdwigtzty
      double precision, dimension (:), allocatable :: thez,they

      double precision, dimension (:), allocatable :: zw,yw

      double precision :: ebeam,ebeammin,ebeammax,debeam,deltae,ezr,ezi,eyr,eyi,wig,
     &  g(nwigefold+1),gsum,be(1000),bw(1000)

      integer isystem
      external isystem

      integer :: ibackspace=8,ic
      character cbs
      equivalence (ibackspace,cbs)

      integer iwrun,ipos(2,nwigefold),kwigerr,npola,kpola,iz,iy,itz,ity,ifrq,kfrq,iwcode,
     &  lz,ly,ltz,lty,lpola,isour,nsource

      integer :: iline=0,kempty,iwig,iwigdum
      integer :: ihtracko,ihfreqo,lunin,lunout,nlast,nl,lunrun,nfirst,ni,l2,l1,k2,k1,istat,ipid,
     &  iend,ianf,m1,m2,n1,n2,lunfis,ieof,lunwef

      character(2048) cline,cstage,cstage0,cins,clineb,cline1
      character(64) c64,cpid

      logical lexist

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

      if (kellip.ne.0.and.nharmell.ne.0) then
        write(6,*)
        write(6,*)'*** WARNING IN WINSTWIGEFOLD: KELLIP.NE.0.AND.NHARMELL.NE.0!'
        write(6,*)'*** This means that for all beam energies the undulator field is different!'
        write(6,*)'*** This is only useful for test purposes!'
        write(6,*)
      endif

      if (KHALBA.ne.0.and.NHHALBA.ne.0) then
        write(6,*)
        write(6,*)'*** WARNING IN WINSTWIGEFOLD: KHALBA.NE.0.AND.NHHALBA.NE.0!'
        write(6,*)'*** This means that for all beam energies the undulator field is different!'
        write(6,*)'*** This is only useful for test purposes!'
        write(6,*)
      endif

      if (KHALBASY.ne.0.and.NHHALBASY.ne.0) then
        write(6,*)
        write(6,*)'*** WARNING IN WINSTWIGEFOLD: KHALBASY.NE.0.AND.NHHALBASY.NE.0!'
        write(6,*)'*** This means that for all beam energies the undulator field is different!'
        write(6,*)'*** This is only useful for test purposes!'
        write(6,*)
      endif

      open(newunit=lunwef,file="wigner.wef")

      kwigerr=0

      ebeammin=dmyenergy*(1.0d0-nsige*espread)
      ebeammax=dmyenergy*(1.0d0+nsige*espread)
      debeam=(ebeammax-ebeammin)/dble(nwigefold-1)
      deltae=espread*dmyenergy

      nwigefold=nwigefold/2*2+1

      ihtracko=ihtrack
      ihfreqo=ihfreq

      ebeam=ebeammin

      gsum=0.0d0
      do iwig=1,nwigefold
        g(iwig)=exp(-((ebeammin+(iwig-1)*debeam-dmyenergy)/deltae)**2/2.0d0)/sqrt(twopi1)/deltae
        gsum=gsum+g(iwig)
      enddo
      if (nwigefold.le.3.or.nosplinewef.ne.0) then
        g=g/gsum
        gsum=1.0d0
      endif

      do iwigefold=1,nwigefold

        write(cline,*)'Winstwigefold:',iwigefold,' of ',nwigefold
        write(6,'(a)',advance='no') cline(1:len_trim(cline))
        do ic=1,len_trim(cline)
           write(6,'(a)',advance='no') cbs
        enddo

        be(iwigefold)=ebeam

        open(newunit=lunout,file="wave.tmp")
        open(newunit=lunin,file="wave.in",status='old')

        do while (.true.)

          read(lunin,'(a)',end=99)cline1

          iline=iline+1
          cline=cline1

          call util_string_trim(cline1,nfirst,nlast)

          if (nfirst.lt.0) then
            write(lunout,'(a)') trim(cline)
            cycle
          endif

          if (nfirst.lt.0.or.cline1(nfirst:nfirst).eq.'!') then
            write(lunout,'(a)') trim(cline)
            cycle
          endif

          call util_lower_case(cline1)

          c64="DMYENERGY"
          call util_string_substring_igncase(cline1,trim(c64),ianf,iend,istat)
          if (istat.eq.0) then
            write(c64,*) ebeam
            call util_string_trim(c64,ni,nl)
            cline="      DMYENERGY=" // c64(ni:nl)
          endif

          c64="NWIGEFOLD"
          call util_string_substring_igncase(cline1,trim(c64),ianf,iend,istat)
          if (istat.eq.0) then
            cline="      NWIGEFOLD=0"
          endif

          c64="IHTRACK"
          call util_string_substring_igncase(cline1,trim(c64),ianf,iend,istat)
          if (istat.eq.0) then
            cline="      IHTRACK=0"
          endif

          c64="IHFREQ"
          call util_string_substring_igncase(cline1,trim(c64),ianf,iend,istat)
          if (istat.eq.0) then
            cline="      IHFREQ=0"
          endif

          write(lunout,'(a)') trim(cline)

        enddo

99      close(lunout)
        close(lunin)

        call util_string_trim(chplatform,k1,k2)
        chplatform=chplatform(k1:k2)
        call util_upper_case(chplatform)

        call util_string_trim(chwavehome,l1,l2)

        cstage=trim(chwavedir)//"/.stage.wig"
        call util_string_trim(cstage,m1,m2)

        kempty=0
        call wave_delete_dir(trim(cstage),kempty,istat)
        call wave_make_dir(trim(cstage),istat)

        if (istat.ne.0) then
          print*,"*** Error winstwigefold: Can't create sub-directory"
          print*,trim(cstage)
        endif

        call wave_copy_file('wave.tmp',trim(cstage)//'/wave.in',istat)
        if (istat.ne.0) then
          print*,"*** Error waveinstances: Can't copy wave_ibunch.tmp to ",
     &      trim(cstage)//'/wave.in'
        endif

        if (kmagseq.ne.0) then
          call wave_copy_file(trim(filemg),trim(cstage),istat)
        endif

        if (irfilf.ne.0) then
          call wave_copy_file(trim(filef),trim(cstage),istat)
        endif

        if (irbtab.ne.0) then
          call wave_copy_file(trim(filetb),trim(cstage),istat)
        endif

        if (kbpolyh.ne.0) then
          call wave_copy_file("wave_bpolyharm_coef.dat",trim(cstage),istat)
        endif

        if (kbpoly3d.ne.0) then
          call wave_copy_file(trim(file3dfit),trim(cstage),istat)
        endif

        if (kbpoly2dh.ne.0) then
          call wave_copy_file(trim(file2dhfit),trim(cstage),istat)
        endif

        if (kbpharm.ne.0) then
          call wave_copy_file(trim(filephfit),trim(cstage),istat)
        endif

        if (kbpoly3d.ne.0) then
          call wave_copy_file(trim(file3dfit),trim(cstage),istat)
        endif

        if (irfilp.ne.0) then
          call wave_copy_file(trim(filep),trim(cstage),istat)
        endif

        if (kbgenesis.ne.0) then
          call wave_copy_file(trim(filegeni),trim(cstage),istat)
          call wave_copy_file(trim(filegenl),trim(cstage),istat)
        endif

        if (irfilb0.ne.0) then
          call wave_copy_file(trim(fileb0),trim(cstage),istat)
        endif

        if (imagspln.gt.0) then
          call wave_copy_file("magjob.dat",trim(cstage),istat)
        endif

        if (irfill0.ne.0) then
          call wave_copy_file(trim(filel0),trim(cstage),istat)
        endif

        if (ifreq2p.eq.0) then
          call wave_copy_file(trim(filefr),trim(cstage),istat)
        endif

        if (irfilsp0.ne.0) then
          call wave_copy_file(trim(filesp0),trim(cstage),istat)
        endif

        if (irfilsto.ne.0) then
          call wave_copy_file(trim(filesto),trim(cstage),istat)
        endif

        if (iampli.gt.0) then
          call wave_copy_file(trim(fileampli),trim(cstage),istat)
        endif

        if (ifilter.ne.0) then
          call wave_copy_file(trim(fileabs),trim(cstage),istat)
        endif

        if (ifilmul.ne.0) then
          open(newunit=lunfis,file=trim(fileam))
          do while (.true.)
            call util_skip_comment_end(lunfis,ieof)
            if (ieof.ne.0) exit
            read(lunfis,'(a)')cline
            call util_string_trim(cline,n1,n2)
            call wave_copy_file(cline(n1:n2),trim(cstage),istat)
          enddo
          close(lunfis)
        endif

        if (ieffi.ne.0) then
          call wave_copy_file(trim(fileff),trim(cstage),istat)
        endif

        if (iefield.ne.0) then
          call wave_copy_file("efield.dat",trim(cstage),istat)
        endif

        if (iubunch.eq.4) then
          call wave_copy_file("fourier-bunch.dat",trim(cstage),istat)
        endif

        if (nberror.ne.0.and.nberrmod.eq.-1) then
          call wave_copy_file("fibonacci_available_cut.dat",trim(cstage),istat)
          call wave_copy_file("fibonacci_used_cut.dat",trim(cstage),istat)
        endif

        if (ibmask.eq.100.or.jbmask.eq.100) then
          call wave_copy_file("wave.bmask",trim(cstage),istat)
        endif

        open(newunit=lunfis,file="wave_input-files.lis")
        do while (.true.)
          call util_skip_comment_end(lunfis,ieof)
          if (ieof.ne.0) exit
          read(lunfis,'(a)')cline
          call util_string_trim(cline,n1,n2)
          if (cline(n1:n2).eq.'wave.in') cycle
          call wave_copy_file(cline(n1:n2),trim(cstage),istat)
        enddo
        close(lunfis)

        if (chplatform.eq.'LINUX'.or.chplatform.eq.'CYGWIN') then

          cline="ln -s "
     &      // chwavehome(l1:l2) // "/bin/wave.exe "
     &      // chwavehome(l1:l2) // "/bin/wave_spawned.exe 2>/dev/null"

          istat=isystem(trim(cline))

          cline="cd " // trim(cstage) // " && "  // chwavehome(l1:l2)
     &      // "/bin/wave_spawned.exe > " //
     &      trim(cstage) // "/wave.log 2>&1 &"

        else if (chplatform.eq.'WINDOWS') then

          cline="cd " // trim(cstage) // " && START /b cmd "  // chwavehome(l1:l2)
     &      // "\bin\wave_spawned.exe > " //
     &      trim(cstage) // "/wave.log"

        endif !LINUX

c        print*,trim(cline)

        call util_break

        istat=isystem(trim(cline))

        iwrun=99
        do while (iwrun.ne.0)
          inquire(file=trim(cstage)//"/wave.status",exist=lexist)
          if (lexist.eqv..true.) then
            open(newunit=lunrun,file=trim(cstage)//chpathsep//"wave.status",status='old')
            read(lunrun,*)iwrun
            close(lunrun)
          endif
          call sleep(1)
        enddo

        open(newunit=lunrun,file=trim(cstage)//chpathsep//"wigner.wav",status='old')
        read(lunrun,'(a)') cline
        read(cline(2:),*) nsource,nfreq,iwigdum,mphasez,mphasey,nwigthetaz,nwigthetay

        if (iwigefold.eq.1) then

          write(lunwef,'(a)') trim(cline)

          allocate(thez(nwigthetaz),they(nwigthetay),
     &      zw(mphasez),yw(mphasey),
     &      fdwigzy(mphasez,mphasey),fdwigtzty(nwigthetaz,nwigthetay))

          if (iwigner.gt.0) then
            npola=4
          else
            npola=1
          endif

          allocate(wigefold(mphasez,mphasey,nwigthetaz,nwigthetay,nfreq,nwigefold,npola),
     &      esourz(mphasez,mphasey,nfreq,npola),
     &      esoury(mphasez,mphasey,nfreq,npola))

          wigefold=0.0d0

        endif

        do isour=1,nsource
          do ifrq=1,nfreq

            do kpola=1,npola

              do ity=1,nwigthetay
                do itz=1,nwigthetaz
                  do iy=1,mphasey
                    do iz=1,mphasez

                      read(lunrun,*) lpola,lz,ly,ltz,lty,kfrq,
     &                  freq(kfrq),zw(iz),yw(iy),thez(itz),they(ity),
     &                  ezr,ezi,eyr,eyi,
     &                  wig,fdwigzy(iz,iy),fdwigtzty(itz,ity)

                      esourz(iz,iy,ifrq,kpola)=dcmplx(ezr,ezi)
                      esoury(iz,iy,ifrq,kpola)=dcmplx(eyr,eyi)

                      if (nosplinewef.ne.0) then
                        wigefold(iz,iy,itz,ity,ifrq,1,kpola)=wigefold(iz,iy,itz,ity,ifrq,1,kpola)
     &                    +g(iwigefold)*debeam*wig
                      else
                        wigefold(iz,iy,itz,ity,ifrq,1,kpola)=wigefold(iz,iy,itz,ity,ifrq,1,kpola)
     &                    +g(iwigefold)*wig
                      endif
                    enddo
                  enddo
                enddo
              enddo

            enddo !kpola
          enddo !nfreq
        enddo !nsource

        close(lunrun)

        ebeam=ebeam+debeam

      enddo !iwigefold=1,nwigefold

      do ifrq=1,nfreq
        do kpola=1,npola
          do ity=1,nwigthetay
            do itz=1,nwigthetaz
              do iy=1,mphasey
                do iz=1,mphasez
                  if (nosplinewef.ne.0) then
                    write(lunwef,'(6I5,15(1pe15.5))') kpola,iz,iy,itz,ity,ifrq,
     &                freq(ifrq),zw(iz),yw(iy),thez(itz),they(ity),
     &                dreal(esourz(iz,iy,ifrq,kpola)),dimag(esourz(iz,iy,ifrq,kpola)),
     &                dreal(esoury(iz,iy,ifrq,kpola)),dimag(esoury(iz,iy,ifrq,kpola)),
     &                wigefold(iz,iy,itz,ity,ifrq,1,kpola),
     &                fdwigzy(iz,iy),fdwigtzty(itz,ity)
                  else
                    do iwig=1,nwigefold
                      bw(iwig)=wigefold(iz,iy,itz,ity,ifrq,iwig,kpola)
                    enddo
                    call util_integral_spline(be,bw,nwigefold,wig)
                    write(lunwef,'(6I5,15(1pe15.5))') kpola,iz,iy,itz,ity,ifrq,
     &                freq(ifrq),zw(iz),yw(iy),thez(itz),they(ity),
     &                dreal(esourz(iz,iy,ifrq,kpola)),dimag(esourz(iz,iy,ifrq,kpola)),
     &                dreal(esoury(iz,iy,ifrq,kpola)),dimag(esoury(iz,iy,ifrq,kpola)),
     &                wig,
     &                fdwigzy(iz,iy),fdwigtzty(itz,ity)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      close(lunwef)

      deallocate(thez,they,zw,yw,esourz,esoury,fdwigzy,fdwigtzty)

      return
      end
