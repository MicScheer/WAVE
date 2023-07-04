*CMZ :  4.00/15 07/04/2022  22.17.21  by  Michael Scheer
*CMZ :  4.00/14 21/12/2021  11.30.58  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_book1(id,chtit,nx,xmin,xmax)

      use mhbook_mod
      implicit none

      integer mh_exists,id,nx,ihkind
      double precision xmin,xmax
      character(*) chtit


      if (nalloc_mh.le.0) then
        call mh_limit(1000)
      endif

      if (mh_exists(id,ihkind).ne.0) then
        print*
        print*,"*** Error in mh_book1: Already existing histogram ",id
        print*,trim(chtit)
        print*
        goto 9999 !return
      endif

      if (nx.lt.1) then
        print*,"*** Error in mh_book1: nx < 1",
     &    id, trim(chtit)
        goto 9999 !return
      endif

      if (xmin.eq.xmax) then
        xmin=xmin-0.5
        xmax=xmax+0.5
      else if (xmin.gt.xmax) then
        print*,"*** Error in mh_book1: xmin > xmax",
     &    id, trim(chtit)
        goto 9999 !return
      endif

      nhist1_mh=nhist1_mh+1
      nhist_mh=nhist_mh+1

      if (nhist_mh.gt.nalloc_mh) then
        call mh_limit(nhist_mh*2)
      endif

      histos_mh(nhist_mh)=hempty_mh

      histos_mh(nhist_mh)%id=id
      histos_mh(nhist_mh)%title=trim(chtit)
      histos_mh(nhist_mh)%nx=nx
      histos_mh(nhist_mh)%xmin=xmin
      histos_mh(nhist_mh)%xmax=xmax
      histos_mh(nhist_mh)%dx=(xmax-xmin)/dble(nx)

      allocate(histos_mh(nhist_mh)%channels(5,nx+2,1))
      histos_mh(nhist_mh)%channels=0.0d0

      if (lun_index_mh.eq.0) open(newunit=lun_index_mh,file='mh_book.lis')

      nhbooked_mh=nhbooked_mh+1
      if (nhist_mh.gt.nhistmax_mh) nhistmax_mh=nhist_mh
      write(lun_index_mh,*)nhbooked_mh,", 1D Histogramm, ",id," ,",trim(chtit)

      lastid_mh=id
      lastind_mh=nhist_mh

9999  continue
      return
      end
