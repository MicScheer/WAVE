*CMZ :  2.66/07 19/02/2010  17.23.15  by  Michael Scheer
*-- Author :    Michael Scheer   17/02/2010
      program wave_remote_post

      implicit none

      integer nfilesp
      parameter (nfilesp=1000)

      double precision, dimension (:), allocatable ::
     &  e,fdmm,fdmrad,fdmms,fdmrads
      double precision ene

      integer job,njobstot,nene,iene,i,ifirst,ilast,istat,iffirst

      character c1
      character(10) stat(nfilesp),cjob
      character(256) file

      open(unit=20,file='wave-remote.jobs',status='old')

      njobstot=0

1     read(20,*,end=9) job

      njobstot=njobstot+1
      if (njobstot.gt.nfilesp) stop '*** Error: Too many files ***'
      goto 1

9     rewind(20)

      do job=1,njobstot
        read(20,*)i,stat(job)
        if (i.ne.job) stop '*** Error: Bad job order ***'
      enddo

      close(20)

      iffirst=0
      do job=1,njobstot
        if (stat(job).eq.'finished') then

          write(cjob,'(I10)')job
          call util_string_firstcharacter(cjob,ifirst,c1,istat)
          call util_string_lastcharacter(cjob,ilast,c1,istat)

          print*,'Reading output of job ',job

          file='wave_on-axis_flux-density-mm2.dat.'//cjob(ifirst:ilast)
          open(unit=20,file=file,status='old')

          file='wave_on-axis_flux-density-mrad2.dat.'//cjob(ifirst:ilast)
          open(unit=21,file=file,status='old')

          if (iffirst.eq.0) then
            nene=0
11          read(20,*,end=99)ene
            nene=nene+1
            goto 11
99          rewind(20)
            iffirst=1
            allocate(e(nene))
            allocate(fdmm(nene))
            allocate(fdmrad(nene))
            allocate(fdmms(nene))
            allocate(fdmrads(nene))
            fdmms=0.0d0
            fdmrads=0.0d0
          endif

          do iene=1,nene
            read(20,*)e(iene),fdmm(iene)
            read(21,*)e(iene),fdmrad(iene)
          enddo

          fdmms=fdmms+fdmm
          fdmrads=fdmrads+fdmrad

          close(20)
          close(21)

        endif
      enddo

      open(unit=20,file='wave_on-axis_flux-density-mm2.dat',status='unknown')
      open(unit=21,file='wave_on-axis_flux-density-mrad2.dat',status='unknown')

      do iene=1,nene
        write(20,*)e(iene),fdmms(iene)
        write(21,*)e(iene),fdmrads(iene)
      enddo

      close(20)
      close(21)

      stop
      end
