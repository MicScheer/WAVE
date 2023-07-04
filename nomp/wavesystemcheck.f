*CMZ :  4.01/00 14/03/2023  11.10.53  by  Michael Scheer
*CMZ :  4.00/17 28/11/2022  15.21.16  by  Michael Scheer
*CMZ :  4.00/07 29/05/2020  12.52.31  by  Michael Scheer
*CMZ :  4.00/04 01/07/2019  14.04.23  by  Michael Scheer
*CMZ :  4.00/03 09/05/2019  10.59.58  by  Michael Scheer
*CMZ :  4.00/02 12/04/2019  15.05.49  by  Michael Scheer
*CMZ :  4.00/01 11/04/2019  14.40.23  by  Michael Scheer
*CMZ :  4.00/00 04/04/2019  12.29.10  by  Michael Scheer
*CMZ :  3.08/01 03/04/2019  11.56.15  by  Michael Scheer
*CMZ :  3.07/01 29/03/2019  14.35.23  by  Michael Scheer
*-- Author :    Michael Scheer   27/03/2019
      subroutine wavesystemcheck

      use iso_c_binding
      use iso_fortran_env !, only: int64
      !use waveenv

      implicit none

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
*KEEP,waveenv.
      include 'waveenv.cmn'
*KEND.

      integer :: nbacksl=0, nslash=0, inam=0

      integer lunin,iutil_fexist,k1,k2,istat,luns,ipos(2,1000),i,nwords
      character(1024) cline,chos,chsys,chshell,chpath,chtest,chpy
      character c1

      kpython=0

      chpythonhome='/usr/bin'
      chpythoncom='python3'
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

      chpy=trim(chpythonhome) // chpathsep // trim(chpythoncom)

      print*,""
      print*,"-- Checking wave_platform.nam --"
      print*,""

      if (iutil_fexist('wave_platform.nam').ne.0) then
        print*,"-- wave_platform.nam found --"
        print*,"-- Reading namelist PLATFORMN --"
        inam=1
        open(newunit=lunin,file="wave_platform.nam",status='old')
        read(lunin,platformn)
        close(lunin)
        call util_string_split_sep(chpythonpath,1000,nwords,ipos,chpathsep,istat)
        chpythonhome=chpythonpath(1:ipos(2,nwords-1))
        chpythoncom=chpythonpath(ipos(1,nwords):len_trim(chpythonpath))
        chpythonpath=''
        print*,"  Python home: ",trim(chpythonhome)
        print*,"  File seperator: ",trim(chpathsep)
        print*,"  Python command: ",trim(chpythoncom)
        print*,""
        print*,"-- Testing Python --"

        call util_test_python3(trim(chpythonhome) // chpathsep // trim(chpythoncom),kpython)

        if (kpython.eq.0) then
          print*," "
          print*,"--------------------------------------- ---------------------"
          print*,"Python3 seems not to be available..., check wave_platform.nam"
          print*,"-----------------------------------    ----------------------"
          print*," "
        endif

        call util_string_trim(chwavepath,k1,k2)
        cline=chwavepath(k1:k2)

        if (cline.eq.'auto'.or.cline.eq.'AUTO'.or.cline.eq.'') then
          chwavepath=trim(chwavehome) // chpathsep // 'bin'
          call util_string_trim(chwavepath,k1,k2)
          chwavepath=chwavepath(k1:k2)
        endif
      else
        print*,"-- wave_platform.nam not found --"
        print*,"-- Will try to guess system -- "
      endif

      print*,""
      print*,"-- Checking SHELL, PATH, and HOME ---"

      call get_environment_variable("SHELL", chshell)
      call get_environment_variable("PATH",chpath)
      call get_environment_variable("HOME",chuserhome)

      print*,""
      print*,"  - Shell is: "
      print*,trim(chshell)
      print*,""
      print*,"  - Path is: "
      print*,trim(chpath)
      print*,""
      print*,"  - Home is: "
      print*,trim(chuserhome)
      print*,""

      do i=1,len_trim(chpath)
        c1=chpath(i:i)
        if (c1.eq.'/') then
          nslash=nslash+1
        else if (c1.eq.'\') then
          nbacksl=nbacksl+1
        endif
      enddo

      if (nbacksl.gt.nslash) then
        print*,"-- Found more backslashes then slashes in PATH, guess WINDOWS system --"
        print*,""

        call get_environment_variable("OS", chos)
        chsys=chos

        if (len_trim(chos).gt.0) then
          chplatform='WINDOWS'
          call get_environment_variable("HOMEPATH",chuserhome)
        else
          print*,'Warning in wavesystemcheck: Assumed WINDOWS, but system-variable OS is not set??'
        endif

        chplatform='WINDOWS'
        chpythonhome=''
        chpythoncom='python'
        chwavepath='auto'
        chrmtree='rmdir /s /q' ! system command to remove directory tree, even if not empty
        chrmdir='rmdir' ! system command to remove empty directory
        chmvfile='move' ! system command to move  file
        chrmfile='del' ! system command to remove file
        chmkdir='mkdir' ! system command to create directory
        chcptree='xcopy /h' ! system command to copy tree
        chcpfile='copy' ! system command to copy file
        chpathsep='\'     ! path seperator      else

      else

        print*,"-- Found less backslashes then slashes in PATH, guess LINUX like system --"

        print*,"-- Will now execute 'uname > .wavesystem' --"
        print*,""

        istat = system('uname > .wavesystem')

        if (istat.eq.0) then
          open(newunit=luns,file='.wavesystem',status='old')
          read(luns,'(a)') chsys
          close(luns)
          chplatform=chsys(1:5)
          if (chsys(1:5).ne.'MINGW') chplatform='LINUX'
        endif !uname

      endif !(nbacksl.gt.nslash) then

      print*,""
      print*,"-- System is " // trim(chplatform) // " --"

      call hostnm(chhostname)
      call getcwd(chwavedir)

      chcwd=chwavedir
      chcompiler=compiler_version()

      call getlog(chuser)

      if (inam.eq.0.and.chplatform.ne.'LINUX' .and.
     &    chplatform.ne.'WINDOWS' .and.
     &    chplatform.ne.'MINGW') then
        print*
        print*,'*** Warning in wavesystemcheck: System seems not to be Windows, Linux, or MINGW ***'
        print*,'*** Provide wavesystem.nam or be aware of problems ***'
        print*
      endif !(chplatform.eq.'LINUX') then

      call get_environment_variable("WAVE",chwavehome)

      if (len_trim(chwavehome).eq.0) then
        call util_string_split_sep(chwavedir,1000,nwords,ipos,chpathsep,istat)
        chtest=chwavedir(1:ipos(2,nwords-1))
        chtest=trim(chtest) // chpathsep // 'bin' // chpathsep // 'wave.exe'
        if (iutil_fexist(trim(chtest)).ne.0) then
          chwavehome=chwavedir(1:ipos(2,nwords-1))
        else
          print*,""
          print*,"*** Warning in wavesystemcheck, could not find: " //
     &      trim(chtest)
          print*,"*** WAVE not properly installed ***"
          print*
        endif
      endif

      print*,""
      if (inam.eq.0) then

        if (chplatform.eq.'LINUX' .or. chplatform.eq.'MINGW') then
          print*,"-- Trying to find python3 path using command 'which' --"
          istat = system('which python3 > .wavesystem')
          open(newunit=luns,file='.wavesystem',status='old')
          read(luns,'(a)') chpythonpath
          close(luns)
        else
          print*,"-- Trying to find python3 path using command 'where' --"
          istat = system('where python')
          if (istat.eq.0) then
            istat = system('where python > .wavesystem')
            if (istat.eq.0) then
              open(newunit=luns,file='.wavesystem',status='old')
              read(luns,'(a)') chpythonpath
              close(luns)
            endif
          endif
        endif !(chplatform.eq.'LINUX') then

        call util_test_python3(chpythonpath,kpython)
        chtest=chpythonpath
        call util_string_split_sep(chtest,1000,nwords,ipos,chpathsep,istat)
        chpythonhome=chpathsep // chtest(ipos(1,1):ipos(2,nwords-1))
        chpythoncom=chtest(ipos(1,nwords):ipos(2,nwords))

      endif

      if (kpython.eq.0.and.inam.eq.0) then
        print*," "
        print*,"---------------------------------------------------------"
        print*,"*** Warning in wavesystemcheck ***"
        print*
        print*,"Python3 seems not to be available nor wave_platform.nam"
        print*,"---------------------------------------------------------"
        print*," "
      else
        print*,""
        print*,"-- Python command is: " // trim(chpythoncom)
      endif

      return
      end
