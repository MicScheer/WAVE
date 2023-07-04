*CMZ :  4.00/17 15/11/2022  10.46.52  by  Michael Scheer
*CMZ :  4.00/07 11/05/2020  15.42.59  by  Michael Scheer
*-- Author :    Michael Scheer   11/05/2020
      module waveenv

      integer, parameter :: maxinstp=100
      integer iwstat(maxinstp),iwinstance,nwinstances,nwgood
      character(2048) chuserhome,chwavehome,chwavedir,chuser

      end module waveenv

      use iso_c_binding
      use iso_fortran_env !, only: int64
      !use waveenv

      implicit none
*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.
      integer kpython

      character chpathsep
      character(32) chplatform, chrmtree, chrmdir,chrmfile,chmkdir,chcpfile,
     &  chcptree,chmvfile

      character(1024) chpythonpath,chhostname,chcompiler,chcwd,chwavepath

      common/platformc/ kpython,chpythonpath,
     &  chplatform, chrmtree, chrmdir,chrmfile,chmkdir,chcpfile,chcptree,
     &  chmvfile,chpathsep,chhostname,chcompiler,chcwd,chwavepath

      namelist/platformn/ chpythonpath,
     &  chplatform, chrmtree, chrmdir,chrmfile,chmkdir,chcpfile,chcptree,
     &  chmvfile,chpathsep,chwavepath

      logical lexist

      integer lunin,k1,k2,istat,luns
      character(1024) cline
      character(32) chos,chsys,chshell

      chpythonpath='python3'
      chplatform='UNKNOWN'
      chwavepath='auto'
      chrmtree='rm -r -f' ! system command to remove directory tree, even if not empty
      chrmdir='rmdir' ! system command to remove empty directory
      chmvfile='mv' ! system command to move  file
      chrmfile='rm' ! system command to remove file
      chmkdir='mkdir' ! system command to create directory
      chcptree='cp -a' ! system command to copy tree
      chcpfile='cp' ! system command to copy file
      chpathsep='/'     ! path seperator

c      inquire(file='wave_platform.nam',exist=lexist)
c      if (lexist.eqv..true.) then
c        open(newunit=lunin,file="wave_platform.nam",status='old')
c        read(lunin,platformn)
c        close(lunin)
c        call util_test_python3(chpythonpath,kpython)
c        if (kpython.eq.0) then
c          print*," "
c          print*,"--------------------------------------- ---------------------"
c          print*,"Python3 seems not to be available..., check wave_platform.nam"
c          print*,"-----------------------------------    ----------------------"
c          print*," "
c        endif
c        call util_string_trim(chwavepath,k1,k2)
c        cline=chwavepath(k1:k2)
c        if (cline.eq.'auto'.or.cline.eq.'AUTO'.or.cline.eq.'') then
c          chwavepath=trim(chwavehome) // trim(chpathsep) // 'bin'
c          call util_string_trim(chwavepath,k1,k2)
c          chwavepath=chwavepath(k1:k2)
c        endif
c      endif

      chcompiler=compiler_version()

      call get_environment_variable("OS", chos)
      call get_environment_variable("SHELL", chshell)

      if (len_trim(chos).gt.0) then
        chsys=chos
        if (chsys.eq.'Windows_NT') then
          chplatform='WINDOWS'
        endif
        istat = system('uname > .wavesystem')
        if (istat.eq.0) then
          open(newunit=luns,file='.wavesystem',status='old')
          read(luns,'(a)') chsys
          close(luns)
          if (chsys(1:5).eq.'MINGW') then
            chplatform='MINGW'
          endif
        endif !uname
      else
        istat = system('uname > .wavesystem')
        if (istat.ne.0) then
          print*,'Warning in wavesystem: Assumed LINUX, but uname did not work??'
        else
          open(newunit=luns,file='.wavesystem',status='old')
          read(luns,'(a)') chsys
          close(luns)
          if (chsys.eq.'Linux') then
            chplatform='LINUX'
          endif !(chsys.eq.'Linux') then
        endif !uname
      endif !len_trim(chos)

      call hostnm(chhostname)
      call getcwd(chcwd)
      call getlog(chuser)

      if (chplatform.ne.'LINUX' .and.
     &    chplatform.ne.'WINDOS' .and.
     &    chplatform.ne.'MINGW') then
        print*
        print*,'*** Warning in wavesystem: System seems not to be Windows, Linux, or MINGW ***'
        print*,'*** Provide wavesystem.nam or be aware of problems ***'
        print*
      endif !(chplatform.eq.'LINUX') then

      end
