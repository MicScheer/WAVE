*CMZ :  2.68/03 14/08/2012  15.02.09  by  Michael Scheer
*CMZ :  2.68/02 26/06/2012  09.23.12  by  Michael Scheer
*-- Author :    Michael Scheer   25/06/2012
      program nsri

c +PATCH,//WAVE/MAIN
c +DECK,nsrwi.

c Program to read srwi-formatted spectra files and convert them to column
c format: y,z, Egamma, Fluxden
c
c WAVE has x as longitudinal direction, thus z is -x from file
c
      implicit none

      real y,z,egam,fluxden,egammin,egammax,zmin,zmax,ymin,ymax,dy,dz,
     &  degam

      integer last,narg,lastin,lastout,lenocc,negam,ny,nz,ianf,iend,istat,
     &  iy,iz,iegam
      character(512) filein,fileout,cline
      character(1024) command
      logical lexist

c      integer ic
c      character c1
c      equivalence (ic,c1)

      narg=iarg()
      if (narg.ne.2) then
        if (narg.ne.1) then
          print*,'Please enter name of input file:'
          read(5,'(a)')filein
        endif
        print*,'Please enter name of output file:'
        read(5,'(a)')fileout
      else
        call getarg(1,filein)
        call getarg(2,fileout)
      endif

      lastin=lenocc(filein)
      lastout=lenocc(fileout)

      inquire(file=filein(1:lastin),exist=lexist)
      if (lexist.eq..false.) then
        print*,'*** Input-File'
        print*,filein(1:lastin)
        print*,'*** not found ***'
        stop '*** nsrwi.exe aborted ***'
      endif

      inquire(file=fileout(1:lastout),exist=lexist)
      if (lexist.eq..true.) then
        command='mv ' // fileout(1:lastout) // ' ' // fileout(1:lastout) // '.bck'
        call system(command(1:lenocc(command)))
      endif

      open(unit=20,file=filein(1:lastin),status='old')

      open(unit=21,file=fileout(1:lastout),form='formatted',
     &  carriagecontrol='list',recl=512,
     &  status='unknown')

 1    read(20,'(a)'),cline

      last=lenocc(cline)

      if (last.le.1.or.
     &    cline(1:1).eq.'%'.or.
     &    cline(1:1).eq.'!'.or.
     &    cline(1:1).eq.'#'.or.
     &    cline(1:1).eq.'@') then

        write(6,*)cline(1:last)

        call util_string_substring(cline(1:last),'startX',ianf,iend,istat)
        if (istat.eq.0) then
          call util_string_substring(cline(1:last),'=',ianf,iend,istat)
          read(cline(ianf+1:last),*)zmin
          zmax=-zmin ! WAVE coord.-system
        endif

        call util_string_substring(cline(1:last),'endX',ianf,iend,istat)
        if (istat.eq.0) then
          call util_string_substring(cline(1:last),'=',ianf,iend,istat)
          read(cline(ianf+1:last),*)zmax
          zmin=-zmax ! WAVE coord.-system
        endif

        call util_string_substring(cline(1:last),'xPoints',ianf,iend,istat)
        if (istat.eq.0) then
          call util_string_substring(cline(1:last),'=',ianf,iend,istat)
          read(cline(ianf+1:last),*) nz
        endif

        call util_string_substring(cline(1:last),'startY',ianf,iend,istat)
        if (istat.eq.0) then
          call util_string_substring(cline(1:last),'=',ianf,iend,istat)
          read(cline(ianf+1:last),*)ymin
        endif

        call util_string_substring(cline(1:last),'endY',ianf,iend,istat)
        if (istat.eq.0) then
          call util_string_substring(cline(1:last),'=',ianf,iend,istat)
          read(cline(ianf+1:last),*)ymax
        endif

        call util_string_substring(cline(1:last),'yPoints',ianf,iend,istat)
        if (istat.eq.0) then
          call util_string_substring(cline(1:last),'=',ianf,iend,istat)
          read(cline(ianf+1:last),*) ny
        endif

        call util_string_substring(cline(1:last),'startPhotonEnergy',
     &    ianf,iend,istat)
        if (istat.eq.0) then
          call util_string_substring(cline(1:last),'=',ianf,iend,istat)
          read(cline(ianf+1:last),*) egammin
        endif

        call util_string_substring(cline(1:last),'endPhotonEnergy',
     &    ianf,iend,istat)
        if (istat.eq.0) then
          call util_string_substring(cline(1:last),'=',ianf,iend,istat)
          read(cline(ianf+1:last),*) egammax
        endif

        call util_string_substring(cline(1:last),'photonEnergyPoints',
     &    ianf,iend,istat)
        if (istat.eq.0) then
          call util_string_substring(cline(1:last),'=',ianf,iend,istat)
          read(cline(ianf+1:last),*) negam
        endif

        goto 1

      endif !classify line

      backspace(20)

      dy=0.0d0
      if (ny.gt.1) dy=(ymax-ymin)/(ny-1)
      dz=0.0d0
      if (nz.gt.1) dz=(zmax-zmin)/(nz-1)
      degam=0.0d0
      if (negam.gt.1) degam=(egammax-egammin)/(negam-1)

      do iy=1,ny
        y=ymin+(iy-1)*dy
        do iz=1,nz
        z=zmin+(iz-1)*dz
          do iegam=1,negam
            egam=egammin+(iegam-1)*degam
            read(20,*)fluxden
            write(21,*)z,y,egam,fluxden
          enddo
        enddo
      enddo

      close(20)
      close(21)

      stop
      end
