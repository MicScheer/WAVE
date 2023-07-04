*CMZ :  4.01/00 12/03/2023  17.18.39  by  Michael Scheer
*CMZ :  4.00/14 31/12/2021  15.13.30  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_fopen(chfin)

      use mhbook_mod
      implicit none

      integer i
      character(*) chfin
      character(2048) chfile


      chfile=trim(adjustl(chfin))

      do i=1,nfile_mh
        if (chfile.eq.chfiles_mh(i)) then
          if (luns_mh(i).eq.0) then
            kfile_mh=i
            chfiles_mh(i)=chfile
            open(newunit=luns_mh(i),file=chfile)
          endif
          goto 9999
        endif
      enddo

      nfile_mh=nfile_mh+1
      if (nfile_mh.gt.nfilep_mh) stop "*** Error in mh_fopen: Dimension nfile_mh exceeded ***"
      chfiles_mh(nfile_mh)=chfile
      open(newunit=luns_mh(nfile_mh),file=chfile)
      kfile_mh=nfile_mh

9999  continue
      return
      end
