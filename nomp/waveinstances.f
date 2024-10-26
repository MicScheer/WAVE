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
      subroutine waveinstances

      use omp_lib
      use clustermod
      !use waveenv
      use bunchmod
      !use waveenv

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
*KEEP,photon.
      include 'photon.cmn'
*KEEP,random.
      include 'random.cmn'
*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.

      integer lun0,lunin,lunclu,kins,ins,m1,m2,n1,n2,irun,istat,lunout,ianf,iend,ipid,
     &  lunfis,ieof,l1,l2,ialldone,masterpid,n1ins,n2ins,lpid,npids,i,idum,
     &  lunpid,lun10,i10,k10,ndim,kcount,lstat,kempty,k1,k2,ierr,kbuncherr,
     &  kstat,luni,lun,lunsi,iel1,iel2,lunspai,lunspao,iutil_fexist,nread,
     &  nfirst,nlast,icheckpid,lunbun,lunbu,ni,nl,nwords

      integer inspid(maxinstp),iwruns(0:maxinstp),ipos(2,maxinstp)
      integer :: iline=0

      real rn(1)

      logical lexist

      character(2048) cline,cstage,cstage0,cins,clineb,cline1
      character(64) c64,cpid

      print*,""
      print*,""

      WRITE(6,*)
      WRITE(6,*)'          *********************************************'
      WRITE(6,*)'          *          PROGRAM WAVE                     *'
      WRITE(6,*)'          *                                           *'
      WRITE(6,*)
*KEEP,wversion.
      include 'wversion.cmn'
*KEND.
      WRITE(6,*)'          *                                           *'
      WRITE(6,*)'          *          Michael Scheer                   *'
      WRITE(6,*)'          *              BESSY                        *'
      WRITE(6,*)'          *********************************************'
      WRITE(6,*)

*KEEP,gplhint.
*KEND.

      !call wavesystem

      kbuncherr=0

1     continue

      open(newunit=lunin,file="wave.in",status='old')
      read(lunin,contrl)
      read(lunin,randomn)
      read(lunin,cluster)
      read(lunin,b0scglobn)
      read(lunin,myfiles)
      read(lunin,bunchn)
      read(lunin,freqn)
      read(lunin,berrorn)
      read(lunin,photonn)
      close(lunin)

      if (kampli.ne.0.or.iundulator.eq.2) then
        !if (iundulator.eq.2) mthreads=-1
        ibunch=0
        if (mthreads.lt.0.or.mthreads.gt.OMP_GET_MAX_THREADS()) then
          mthreads=OMP_GET_MAX_THREADS()
        else if (mthreads.eq.0) then
          mthreads=1
        endif
        mampthreads=mthreads
        iamppin=ipin
        iamppincirc=ipincirc
      endif

      if (ibunch.ne.0) then
        if (nbunch.le.0) nbunch=1
        if (neinbunch.le.0) neinbunch=1
      endif

      if (kbuncherr.gt.0) then
        neinbunch=kbuncherr
      else if (kbuncherr.lt.0) then
        nbunch=kbuncherr
      endif

      nwinstances=1

      if (icluster.lt.0) then
        open(newunit=lunclu,file="wave.ins",status='old')
        read(lunclu,*)iwinstance,nwinstances
        close(lunclu)
      endif

      if (
     &    ipin.eq.0.and.ibunch.eq.0
     &    .or.
     &    ibunch.ne.0.and.neinbunch*nbunch.eq.1
     &    ) then
        mthreads=0
      endif

      if (mthreads.lt.0.or.mthreads.gt.OMP_GET_MAX_THREADS()) then
        mthreads=OMP_GET_MAX_THREADS()
      endif

      if (ibunch.ne.0.and.icluster.eq.0.and.mthreads.ne.0) then
        if (mthreads.lt.0) then
          mthreads=OMP_GET_MAX_THREADS()
        endif
        icluster=1
      endif

C--- RANDOM NUMBERS

      if (kbuncherr.eq.0) then
        open(newunit=lunin,file="WAVE_CODE.DAT")
        call util_skip_comment_end(lunin,ieof)
        if (ieof.ne.0) then
          irun=1
          rewind(lunin)
          write(lunin,*) irun
        else
          read(lunin,*)irun
          irun=irun+1
        endif
        close(lunin)
      endif

      irnsize=64

      if (irnmode.lt.0) then
        open(newunit=lun,file='wave.seeds',status='old')
        read(lun,*)irnsize
        do i=1,irnsize
          read(lun,*)idum,irnseed(i)
        enddo
        close(lun)
        call util_random_set_seed(irnsize,irnseed)
      else if (irnmode.eq.1) then
        call util_random_set_seed(irnsize,irnseed)
      else
        call util_random_init(irnsize,irnseed)
      endif

      if (icluster.lt.0) then
        irnsize=64
        irnseed=irnseed+(iwinstance-1)
        call util_random_set_seed(irnsize,irnseed)
      endif

      call util_random_get_seed(irnsize,irnseedi)

      open(newunit=lunsi,file='wave_start.seeds')
      write(lunsi,*)irnsize,irun
      do i=1,irnsize
        write(lunsi,*)i,irnseedi(i)
      enddo
      close(lunsi)

      nbuncho=nbunch
      neinbuncho=neinbunch

      if (ibunch.eq.0.and.icluster.ne.0) then
        print*,"*** WARNING in WAVEINSTANCES: ICLUSTER .ne. 0, but IBUNCH .eq. 0 ***"
        print*,"ICLUSTER set ICLUSTER=0"
        icluster=0
        return
      endif

      wppath=trim(chwavedir)
      call util_string_split_sep(wppath,maxinstp,nwords,ipos,chpathsep,istat)
      chwstage=wppath(ipos(1,nwords):ipos(2,nwords))

      open(newunit=lunclu,file="wave.pid",status='old')
      read(lunclu,*)lpid
      if (lpid.ne.kpid) then
        stop "*** ERROR  in waveinstances: Bad PID ***"
      endif
      close(lunclu)

      nwinstances=1
      nwgood=1

      if (icluster.lt.0) then
        open(newunit=lunclu,file="wave.ins",status='old')
        read(lunclu,*)iwinstance,nwinstances
        close(lunclu)
      endif

      if (icluster.ge.0) then
        nwinstances=max(icluster,mthreads,1)
        if (ibunch.ne.0) mthreads=0
        if (nwinstances.gt.maxinstp) then
          print*,"--- Warning in waveinstances: Number of instances limited to ",
     &      maxinstp
          nwinstances=maxinstp
        endif
        if (ibunch.ne.0) then
          if (nbunch.eq.1.and.neinbunch.ge.1) then
            if (neinbunch.lt.nwinstances) then
              nwinstances=neinbunch
            endif
            if (neinbunch/nwinstances.lt.2) then
              nwinstances=1
            endif
          else if (nbunch.gt.1.and.neinbunch.eq.1) then
            if (nbunch.lt.nwinstances) then
              nwinstances=nbunch
            endif
          endif
        endif !(ibunch.ne.0) then
        iwstat=1
        inspid=0
        masterpid=kpid
      else
        inspid(iwinstance)=kpid
        goto 9999
      endif

      if (ibunch.ne.0.and.iubunch.eq.3) then
        open(newunit=lunspai,file="wave_phasespace.dat")
        nread=0
        ieof=0
        do i=1,nbunch*neinbunch
          call util_skip_comment_end(lunspai,ieof)
          read(lunspai,'(a)',iostat=ierr) cline
          if (ierr.ne.0) then
            if (nread.eq.0) then
              print*,"*** Error in waveinstaces: no electrons from wave_phasespace.dat ***"
              stop '*** Program WAVE aborted ***'
            else if (neinbunch.eq.1) then
              print*,"*** Warning in waveinstaces: Could read only",
     &          nread," electrons from wave_phasespace.dat ***"
              print*,"Setting NBUNCH = ",(nread/nwinstances)*nwinstances
              NBUNCH=(nread/nwinstances)*nwinstances
              exit
            else if (nbunch.eq.1) then
              print*,"*** Warning in waveinstaces: Could read only",
     &          nread," electrons from wave_phasespace.dat ***"
              print*,"Setting NEINBUNCH = ",(nread/nwinstances)*nwinstances
              NEINBUNCH=(nread/nwinstances)*nwinstances
              exit
            else
              print*,"*** Error in waveinstaces: Could read only",
     &          nread," electrons from wave_phasespace.dat ***"
              stop '*** Program WAVE aborted ***'
            endif
          endif
          nread=nread+1
        enddo
        close(lunspai)
      endif !iubunch

      if (icluster.gt.0) then

        open(newunit=lunclu,file="wave.ins")
        write(lunclu,*)"0 ",nwinstances
        close(lunclu)

        ! Set up the working directories

        open(newunit=lunbun,file="wave_ibunch.tmp")
        open(newunit=lunout,file="wave.tmp")
        open(newunit=lunin,file="wave.in",status='old')

        do while (.true.)

          read(lunin,'(a)',end=99)cline1
          iline=iline+1
          cline=cline1

          !print*,iline,trim(cline)
          call util_string_trim(cline1,nfirst,nlast)
          !if (nfirst.lt.0) cycle
          if (nfirst.lt.0) then
            write(lunout,'(a)') trim(cline)
            write(lunbun,'(a)') trim(cline)
            cycle
          endif

          if (nfirst.lt.0.or.cline1(nfirst:nfirst).eq.'!') then
            write(lunout,'(a)') trim(cline)
            write(lunbun,'(a)') trim(cline)
            cycle
          endif

          call util_lower_case(cline1)

          c64="CODE="
          call util_string_substring_igncase(cline1,trim(c64),ianf,iend,istat)
          if (istat.eq.0) then
            write(lunout,'(a)') trim(cline)
            write(lunbun,'(a)') trim(cline)
            cycle
          endif

          c64="irnmode=0"
          call util_string_substring_igncase(cline1,trim(c64),ianf,iend,istat)
          if (istat.eq.0) then
            cline(ianf:iend)="IRNMODE=2"
          endif

          c64="iclubun="
          call util_string_substring_igncase(cline1,trim(c64),ianf,iend,istat)
          if (istat.eq.0) then
            write(c64,*) nbunch
            call util_string_trim(c64,ni,nl)
            cline="      ICLUBUN=" // c64(ni:nl)
          endif

          c64="icluster="
          call util_string_substring_igncase(cline1,trim(c64),ianf,iend,istat)
          if (istat.eq.0) then
            cline="      ICLUSTER=-1"
          endif

          c64="mthreads="
          call util_string_substring_igncase(cline1,trim(c64),ianf,iend,istat)
          if (istat.eq.0) then
            cline="      MTHREADS=0"
          endif

          if (nbunch.gt.1) then
            c64="nbunch="
            call util_string_substring_igncase(cline1,trim(c64),ianf,iend,kstat)
            c64="neinbunch="
            call util_string_substring_igncase(cline1,trim(c64),ianf,iend,istat)
            if (kstat.eq.0.and.istat.ne.0) then
                write(c64,*) nbunch/nwinstances
                write(lunout,'(a)')"      NBUNCH=" // trim(c64)
                write(lunbun,'(a)')"      NBUNCH=" // trim(c64)
                cycle
            endif
          else if (neinbunch.gt.1) then
            c64="neinbunch="
            call util_string_substring_igncase(cline1,trim(c64),ianf,iend,kstat)
            if (kstat.eq.0) then
                write(c64,*) neinbunch/nwinstances
                write(lunout,'(a)')"      NEINBUNCH=" // trim(c64)
                write(lunbun,'(a)')"      NEINBUNCH=" // trim(c64)
                cycle
            endif
          endif

          if (ibunch.eq.-1) then
            c64="ibunch="
            call util_string_substring_igncase(cline1,trim(c64),ianf,iend,kstat)
            if (kstat.eq.0) then
                write(c64,*) neinbunch/nwinstances
                write(lunout,'(a)')"      IBUNCH=1"
                write(lunbun,'(a)')"      IBUNCH=-1"
                cycle
            endif
          endif

          write(lunout,'(a)') trim(cline)
          write(lunbun,'(a)') trim(cline)

        enddo

99      close(lunout)
        close(lunbun)
        close(lunin)

        kins=1
        if (ibunch.ne.0.and.iubunch.eq.3) kins=0

        do ins=kins,nwinstances

          cstage=trim(chwavedir)//"/.stage."
          call util_string_trim(cstage,m1,m2)
          write(cins,*)ins
          call util_string_trim(cins,n1,n2)
          cstage=cstage(m1:m2)//cins(n1:n2)

          kempty=0
          call wave_delete_dir(trim(cstage),kempty,istat)
          call wave_make_dir(trim(cstage),istat)
          if (istat.ne.0) then
            print*,"*** Error waveinstances: Can't create sub-directory"
            print*,trim(cstage)
          endif
          if (kins.eq.0.and.ins.eq.0.or.kins.eq.1.and.ins.eq.1) then
            call wave_copy_file('wave_ibunch.tmp',trim(cstage)//'/wave.in',istat)
            if (istat.ne.0) then
              print*,"*** Error waveinstances: Can't copy wave_ibunch.tmp to ",
     &          trim(cstage)//'/wave.in'
            endif
          else
            call wave_copy_file('wave.tmp',trim(cstage)//'/wave.in',istat)
            if (istat.ne.0) then
              print*,"*** Error waveinstances: Can't copy wave.tmp to ",
     &          trim(cstage)//'/wave.in'
            endif
          endif

          open(newunit=lunin,file=trim(cstage)//"/WAVE_CODE.DAT")

          call util_random(1,rn)
          irun=rn(1)*10000000-1
          write(lunin,*)irun
          iwruns(ins)=irun
          close(lunin)

          open(newunit=lunclu,file=trim(cstage)//"/wave.ins")
          write(lunclu,*)ins,nwinstances
          close(lunclu)

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

        enddo !instances

        !Create instances of WAVE

        call util_string_trim(chwavehome,l1,l2)
        call util_string_trim(chwavedir,m1,m2)

        print*,""
        call zeit(6)

        print*,""
        print*,"     Spawning WAVE instances"
        print*,""

        if (iubunch.eq.3) then
          open(newunit=lunspai,file="wave_phasespace.dat",status='old')
        endif

        do ins=kins,nwinstances

          cstage=trim(chwavedir)//"/.stage."
          call util_string_trim(cstage,m1,m2)

          write(cins,*)ins
          call util_string_trim(cins,n1,n2)
          cstage=cstage(m1:m2)//cins(n1:n2)

          if (iubunch.eq.3) then
            open(newunit=lunspao,file=trim(cstage)//"/wave_phasespace.dat")
            iel1=nbunch*neinbunch/nwinstances*(ins-1)+1
            iel2=nbunch*neinbunch/nwinstances*ins
            if (ins.eq.0) then
              iel1=1
              iel2=1
            else if (ins.eq.1) then
              rewind(lunspai)
            endif
            nread=0
            do i=iel1,iel2
              call util_skip_comment_end(lunspai,ieof)
              read(lunspai,'(a)',iostat=ierr) cline
              if (ierr.ne.0) then
                if (nread.eq.0) then
                  print*,"*** Error in waveinstaces: no electrons from wave_phasespace.dat ***"
                  stop '*** Program WAVE aborted ***'
                else if (neinbunch.eq.1) then
                  print*,"*** Warning in waveinstaces: Could read only",
     &              nread," electrons from wave_phasespace.dat ***"
                  print*,"Setting NBUNCH = ",(nread/nwinstances)*nwinstances
                  kbuncherr=-(nread/nwinstances)*nwinstances
                  close(lunspai)
                  close(lunspao)
                  goto 1
                else if (nbunch.eq.1) then
                  print*,"*** Warning in waveinstaces: Could read only",
     &              nread," electrons from wave_phasespace.dat ***"
                  print*,"Setting NEINBUNCH = ",(nread/nwinstances)*nwinstances
                  close(lunspai)
                  close(lunspao)
                  kbuncherr=(nread/nwinstances)*nwinstances
                  goto 1
                else
                  print*,"*** Error in waveinstaces: Could read only",
     &              nread," electrons from wave_phasespace.dat ***"
                  stop '*** Program WAVE aborted ***'
                endif
              endif
              nread=nread+1
              write(lunspao,'(a)') trim(cline)
            enddo
            flush(lunspao)
            close(lunspao)
          endif !iubunch

          call util_string_trim(chplatform,k1,k2)
          chplatform=chplatform(k1:k2)
          call util_upper_case(chplatform)

          if (chplatform.eq.'LINUX'.or.chplatform.eq.'CYGWIN') then

            if (ins.gt.0.and.iubunch.eq.3) then
              call wave_copy_file(trim(cstage0)//'/wave_source.clu',
     &          trim(cstage)//'/wave_source.clu',istat)
              if (istat.ne.0) then
                print*,"*** Error waveinstances: Can't copy wave_source.clu to "
              endif
            endif

            cline="ln -s "
     &        // chwavehome(l1:l2) // "/bin/wave.exe "
     &        // chwavehome(l1:l2) // "/bin/wave_spawned.exe 2>/dev/null"

            istat=system(trim(cline))

            cline="cd " // trim(cstage) // " && "  // chwavehome(l1:l2)
     &        // "/bin/wave_spawned.exe > " //
     &        trim(cstage) // "/wave.log 2>&1 &"

            print*,trim(cline)
            istat=system(trim(cline))

            ipid=0
            do while(ipid.eq.0)
              inquire(file=trim(cstage)//"/wave.pid",exist=lexist)
              if (lexist.eqv..true.) then
                open(newunit=lunpid,file=trim(cstage)//"/wave.pid",status='old')
                read(lunpid,*)inspid(ins)
                close(lunpid)
                exit
              endif
              call sleep(1)
            enddo

          else if (chplatform.eq.'WINDOWS') then

            if (ins.gt.0.and.ibunch.eq.3) then
              call wave_copy_file(trim(cstage0)//'\wave_source.clu',
     &          trim(cstage)//'/wave_source.clu',istat)
              if (istat.ne.0) then
                print*,"*** Error waveinstances: Can't copy wave_source.clu to "
              endif
            endif

            cline="cd " // trim(cstage) // " && START /b cmd "  // chwavehome(l1:l2)
     &        // "\bin\wave_spawned.exe > " //
     &        trim(cstage) // "/wave.log"

            print*,trim(cline)
            istat=system(trim(cline))

          endif !LINUX

          if (ins.eq.0) then
            cstage0=cstage
            do while (.true.)
              call sleep(1)
              inquire(file=trim(cstage)//"/wave.status",exist=lexist)
              if (lexist.eqv..true.) then
                open(newunit=lun0,file=trim(cstage)//'/wave.status')
                read(lun0,*)istat
                close(lun0)
                if (istat.eq.0) exit
              endif
              print*,"Waiting for job on stage.0"
            enddo
          endif !(ins.gt.0) then

        enddo !instances

        if (iubunch.eq.3) close(lunspai)

        !Check instances

        istat=0
        kstat=0
        ialldone=0
        k10=0

        do while (ialldone.lt.nwinstances)

          call sleep(3)

          do ins=1,nwinstances

            if (iwstat(ins).le.0) then
              cycle
            endif

            call wave_check_pid(inspid(ins),icheckpid)

            if (icheckpid.eq.-1) then
              print*,"*** Warning in waveinstances: Could not check if subprocess is running."
              print*,"*** Check wave.log in subdirectories .stage..., if program hangs"
            endif

            write(cins,*) ins
            call util_string_trim(cins,n1,n2)

            cstage=trim(chwavedir)//"/.stage."
            call util_string_trim(cstage,m1,m2)
            cstage=cstage(m1:m2)//cins(n1:n2)

            inquire(file=trim(cstage)//"/wave.status",exist=lexist)
            if (lexist.eqv..false.) cycle

            open(newunit=lunin,file=trim(cstage)//"/wave.status")
            read(lunin,*)kstat
            iwstat(ins)=kstat
            close(lunin)

            if (iwstat(ins).eq.0) then
              print*,""
              print*,"--- WAVE instance",ins," has finished ---"
              ialldone=ialldone+1
              if (ialldone.eq.nwinstances) exit
            else if (iwstat(ins).ne.1.or.icheckpid.eq.0) then
              print*,""
              print*,"--- WAVE instance",ins," has probably crashed ---"
              print*,"*** Will be ignored, be careful ***"
              print*,""
              iwstat(ins)=-1
              ialldone=ialldone+1
              if (ialldone.eq.nwinstances) exit
            endif

          enddo !nwinstances

          if (ialldone.lt.nwinstances.and.kcount.eq.0) then
            print*,""
            CALL ZEIT(6)
            print*,"     Counting from 1 to 10 to show progress"
            print*,""
            kcount=1
          endif

          cline=trim(cstage)//"/wave.n10"
          inquire(file=trim(cline),exist=lexist)

          if ((lexist.eqv..true.) .and. kcount.gt.0.and.iwstat(ins).ne.-1) then
            open(newunit=lun10,
     &        file=trim(cstage)//"/wave.n10")
            read(lun10,'(a)',end=8)c64
            read(lun10,*,end=8)i10
 8          close(lun10)
            if (i10.gt.k10) then
              print*,i10,"  ",c64
              k10=i10
            endif
          endif

        enddo

        nwgood=0
        if (iwbunch.ne.0) then
          open(newunit=lunbu,file="wave_phasespace_bunch.dat",status='new')
          close(lunbu)
          open(newunit=lunbu,file="wave_phasespace_bunch.dat",access='append')
        endif
        do ins=1,nwinstances
          if (iwstat(ins).eq.0) then
            nwgood=nwgood+1
            if (iwbunch.ne.0) then
              cstage=trim(chwavedir)//"/.stage."
              call util_string_trim(cstage,m1,m2)
              write(cins,*)ins
              call util_string_trim(cins,n1,n2)
              cstage=cstage(m1:m2)//cins(n1:n2)
              cline=trim(cstage)//chpathsep//"wave_phasespace_bunch.dat"
              open(newunit=lunbun,file=trim(cline))
              do while (.true.)
                read(lunbun,'(a)',end=91)cline
                write(lunbu,*)trim(cline)
              enddo
91            close(lunbun)
            endif
          endif
        enddo
        if (iwbunch.ne.0) close(lunbu)

        print*,""

        if (nwgood.eq.0) then
          print*,"*** Error in waveinstances: All instances crashed ***"
          print*,"*** Check log files .stage.../wave.log ***"
          print*,"*** Check if needed files are not listed in wave_input-files.lis ***"
          stop "*** Program WAVE aborted ***"
        endif

      endif !icluster

9999  continue

      return
      end
