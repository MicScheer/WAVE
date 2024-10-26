*CMZ :          02/05/2024  12.38.37  by  Michael Scheer
*CMZ :  4.01/05 20/04/2024  09.48.41  by  Michael Scheer
*CMZ :  4.00/15 14/03/2022  10.05.01  by  Michael Scheer
*CMZ :  4.00/14 10/02/2022  17.40.28  by  Michael Scheer
*CMZ :  4.00/13 16/12/2021  12.18.27  by  Michael Scheer
*CMZ :  3.02/00 01/09/2014  11.11.30  by  Michael Scheer
*CMZ :  3.01/06 23/06/2014  16.20.53  by  Michael Scheer
*CMZ :  3.01/05 13/06/2014  09.06.46  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  13.12.28  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  09.54.48  by  Michael Scheer
*-- Author :    Michael Scheer   19/01/2012
      subroutine hfm(nid,buffD)

      use clustermod
      use mhbook_mod

*KEEP,gplhint.
*KEND.

      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEEP,ntupinfo.
      include 'ntupinfo.cmn'
*KEND.

      double precision buffD(*)
      integer nid,i,k,nvar
      real buff(100),rlow(100),rhigh(100)
      character(80) tit
      character(10) chtag(100)
      character(7) cname,cname2
      character c

      if (icluster.eq.0.or.icluster.eq.9999) then
        call mh_filln(nid,buffd)
        return
      endif

      if (icluster.lt.0) then
        if (nid.eq.3601) then
          write(nscr3601,*)buffd(1:36)
        else if (nid.eq.3600) then
          write(nscr3600,*)buffd(1:2)
        else if (nid.eq.3700) then
          write(nscr3700,*)buffd(1:34)
        else if (nid.eq.4600) then
          write(nscr4600,*)buffd(1:5)
        else if (nid.eq.4700) then
          write(nscr4700,*)buffd(1:12)
        else if (nid.eq.30) then
          write(nscr30,*)buffd(1:41)
        else if (nid.eq.7777) then
          write(nscr7777,*)buffd(1:14)
        endif
      else if ( ! they are filled in wpampntup
     &    nid.ne.3601.and.
     &    nid.ne.3600.and.
     &    nid.ne.3700.and.
     &    nid.ne.4600.and.
     &    nid.ne.4700.and.
     &    nid.ne.30.and.
     &    nid.ne.7777) then
        call mh_filln(nid,buffd)
        return
      endif

      if (iroottrees.ge.0) then
      endif

      if (iroottrees.eq.0) return

      return
      end
