*CMZ :  4.00/14 16/12/2021  15.05.33  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_h2out(ind)

      use mhbook_mod
      implicit none

      double precision xmin,xmax,dx,x,ymin,ymax,dy,y
      integer ind,lunhis,ihkind,nx,ix,ny,iy


      if (histos_mh(ind)%ktouched.eq.2) goto 9999

      if (ind.le.0.or.ind.gt.nhist_mh) then
        print*,"*** Error in mh_h2out: Histogram index out of range",ind
        goto 9999 !return
      endif

      if (kfile_mh.le.0.or.kfile_mh.gt.nfile_mh) then
        print*,"*** Error in mh_h2out: File index out of range",kfile_mh
        goto 9999 !return
      endif

      if (luns_mh(kfile_mh).eq.0) then
        print*,"*** Error in mh_h2out: File not open for writing",
     &    trim(chfiles_mh(kfile_mh))
        goto 9999 !return
      endif

      lunhis=luns_mh(kfile_mh)
      if (histos_mh(ind)%ny.gt.0) then
        ihkind=2
      else
        ihkind=1
      endif

      write(lunhis,*)'! ------------------------------------------------'
      write(lunhis,*)histos_mh(ind)%id,ihkind,'   ! id and kind of histogram or Ntuple'

      write(lunhis,*)len_trim(histos_mh(ind)%title),'   ! length of title'
      write(lunhis,'(a)')trim(histos_mh(ind)%title)

      nx=histos_mh(ind)%nx
      xmin=histos_mh(ind)%xmin
      xmax=histos_mh(ind)%xmax
      dx=histos_mh(ind)%dx

      ny=histos_mh(ind)%ny
      ymin=histos_mh(ind)%ymin
      ymax=histos_mh(ind)%ymax
      dy=histos_mh(ind)%dy

      write(lunhis,'(i10,2e25.12e3,a)') nx,histos_mh(ind)%xmin,histos_mh(ind)%xmax,
     &  '   ! nx, xmin, xmax'
      write(lunhis,'(i10,2e25.12e3,a)') ny,histos_mh(ind)%ymin,histos_mh(ind)%ymax,
     &  '   ! ny, ymin, ymax'
      write(lunhis,'(i10,e25.12e3,a)')
     &  histos_mh(ind)%nentries,histos_mh(ind)%hsum,
     &  '   ! number of entries and weighted sum'

      write(lunhis,'(2e25.12e3,a)')histos_mh(ind)%xmean,histos_mh(ind)%xrms,
     &  '   ! mean and rms'
      write(lunhis,'(2e25.12e3,a)')histos_mh(ind)%ymean,histos_mh(ind)%yrms,
     &  '   ! mean and rms'

      write(lunhis,*)
     &  histos_mh(ind)%nuu+histos_mh(ind)%num+histos_mh(ind)%nuo,
     &  histos_mh(ind)%nou+histos_mh(ind)%nom+histos_mh(ind)%noo,
     &  '   ! under- and overflows in x'

      write(lunhis,'(2e25.12e3,a)')
     &  histos_mh(ind)%wnuu+histos_mh(ind)%wnum+histos_mh(ind)%wnuo,
     &  histos_mh(ind)%wnou+histos_mh(ind)%wnom+histos_mh(ind)%wnoo,
     &  '   ! weighted under- and overflows in x'

      write(lunhis,*)
     &  histos_mh(ind)%nuu+histos_mh(ind)%nmu+histos_mh(ind)%nou,
     &  histos_mh(ind)%nuo+histos_mh(ind)%nmo+histos_mh(ind)%noo,
     &  '   ! under- and overflows in y'

      write(lunhis,'(2e25.12e3,a)')
     &  histos_mh(ind)%wnuu+histos_mh(ind)%wnmu+histos_mh(ind)%wnou,
     &  histos_mh(ind)%wnuo+histos_mh(ind)%wnmo+histos_mh(ind)%wnoo,
     &  '   ! weighted under- and overflows in y'

      write(lunhis,*)
     &  histos_mh(ind)%nuu,histos_mh(ind)%nou,
     &  '   ! under- and overflows in x for underflows in y'
      write(lunhis,'(2e25.12e3,a)')
     &  histos_mh(ind)%wnuu,histos_mh(ind)%wnou,
     &  '   ! weighted under- and overflows in x for underflows in y'

      write(lunhis,*)
     &  histos_mh(ind)%nuo,histos_mh(ind)%noo,
     &  '   ! under- and overflows in x for overflows in y'
      write(lunhis,'(2e25.12e3,a)')
     &  histos_mh(ind)%wnuo,histos_mh(ind)%wnoo,
     &  '   ! weighted under- and overflows in x for overflows in y'

      y=ymin-dy/2.0d0
      do iy=1,ny
        y=y+dy
        x=xmin-dx/2.0d0
        do ix=1,nx
          x=x+dx
          write(lunhis,'(2i10,3e25.12e3,a)') ix,iy,x,y,
     &      histos_mh(ind)%channels(2,ix,iy),
     &      '   ! ix, iy, x, y, content'
        enddo !nx
      enddo !ny

      histos_mh(ind)%ktouched=2

      lastid_mh=histos_mh(ind)%id
      lastind_mh=ind

9999  continue
      return
      end
