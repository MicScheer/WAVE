*CMZ :  4.01/00 12/03/2023  17.18.39  by  Michael Scheer
*CMZ :  4.00/14 31/12/2021  15.13.30  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_path(chin)

      use mhbook_mod
      implicit none

      integer i
      character(*) chin
      character(lenpathp_mh) chpath


      chpath=trim(adjustl(chin))
      if (chpath.eq.chdir_mh) return

      chdir_old_mh=chdir_mh

      do i=1,npath_mh
        if (chpath.eq.chpath_mh(i)) then
          chdir_mh=chpath
          kpath_mh=i
          goto 9999
        endif
      enddo

      npath_mh=npath_mh+1
      if (npath_mh.gt.npathp_mh) stop "*** Error in mh_path: Dimension npath_mh exceeded ***"
      chpath_mh(npath_mh)=chpath
      chdir_mh=chpath
      kpath_mh=npath_mh

9999  continue
      return
      end
