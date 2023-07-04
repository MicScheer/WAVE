*CMZ :  4.00/14 21/12/2021  11.26.21  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_book2(id,chtit,nx,xmin,xmax,ny,ymin,ymax)

      use mhbook_mod
      implicit none

      integer mh_exists,id,nx,ny,ihkind
      double precision xmin,xmax,ymin,ymax
      character(*) chtit


      if (nalloc_mh.le.0) then
        call mh_limit(1000)
      endif

      if (mh_exists(id,ihkind).ne.0) then
        print*,"*** Error in mh_book2: Already existing histogram ",
     &    id, trim(chtit)
        print*,trim(chtit)
        goto 9999 !return
      endif

      if (nx.lt.1) then
        print*,"*** Error in mh_book2: nx < 1",
     &    id, trim(chtit)
        goto 9999 !return
      endif

      if (nx.le.0) then
        print*,"*** Error in mh_book2: nx <=0",id, trim(chtit)
        goto 9999 !return
      endif

      if (xmin.ge.xmax) then
        print*,"*** Error in mh_book2: xmin >= xmax",
     &    id, trim(chtit)
        goto 9999 !return
      endif

      if (ny.lt.1) then
        print*,"*** Error in mh_book2: ny < 1",
     &    id, trim(chtit)
        goto 9999 !return
      endif

      if (ny.le.0) then
        print*,"*** Error in mh_book2: ny <=0",id, trim(chtit)
        goto 9999 !return
      endif

      if (ymin.ge.ymax) then
        print*,"*** Error in mh_book2: ymin >= ymax",
     &    id, trim(chtit)
        goto 9999 !return
      endif

      nhist2_mh=nhist2_mh+1
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
      histos_mh(nhist_mh)%ny=ny
      histos_mh(nhist_mh)%ymin=ymin
      histos_mh(nhist_mh)%ymax=ymax
      histos_mh(nhist_mh)%dy=(ymax-ymin)/dble(ny)

      allocate(histos_mh(nhist_mh)%channels(6,nx+2,ny+2))
      histos_mh(nhist_mh)%channels=0.0d0

      if (lun_index_mh.eq.0) open(newunit=lun_index_mh,file='mh_book.lis')

      nhbooked_mh=nhbooked_mh+1
      if (nhist_mh.gt.nhistmax_mh) nhistmax_mh=nhist_mh
      write(lun_index_mh,*)nhbooked_mh,", 2D Histogramm, ",id," ,",trim(chtit)

      lastid_mh=id
      lastind_mh=nhist_mh

9999  continue
      return
      end
