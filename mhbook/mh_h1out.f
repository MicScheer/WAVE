*CMZ :  4.00/14 16/12/2021  15.05.33  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_h1out(ind)

      use mhbook_mod
      implicit none

      double precision xmin,xmax,dx,x
      integer ind,lunhis,ihkind,nx,ix


      if (histos_mh(ind)%ktouched.eq.2) goto 9999

      if (ind.le.0.or.ind.gt.nhist_mh) then
        print*,"*** Error in mh_h1out: Histogram index out of range",ind
        goto 9999 !return
      endif

      if (kfile_mh.le.0.or.kfile_mh.gt.nfile_mh) then
        print*,"*** Error in mh_h1out: File index out of range",kfile_mh
        goto 9999 !return
      endif

      if (luns_mh(kfile_mh).eq.0) then
        print*,"*** Error in mh_h1out: File not open for writing",
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

      write(lunhis,'(i10,2e25.12e3,a)') nx,histos_mh(ind)%xmin,histos_mh(ind)%xmax,
     &  '   ! nx, xmin, xmax'
      write(lunhis,'(i10,e25.12e3,a)')
     &  histos_mh(ind)%nentries,histos_mh(ind)%hsum,'   ! number of entries and weighted sum'
      write(lunhis,'(2e25.12e3,a)')histos_mh(ind)%xmean,histos_mh(ind)%xrms,'   ! mean and rms'

      write(lunhis,*)histos_mh(ind)%nunder,histos_mh(ind)%nover,
     &  '   ! under- and overflows in x'
      write(lunhis,'(2e25.12e3,a)')
     &  histos_mh(ind)%channels(2,nx+1,1),histos_mh(ind)%channels(2,nx+2,1),
     &  '   ! weighted under- and overflows in x'

      x=xmin-dx/2.0d0

      do ix=1,nx
        x=x+dx
        write(lunhis,'(i10,2e25.12e3,a)') ix,x,histos_mh(ind)%channels(2,ix,1),
     &    '   ! ix, x, content'
      enddo !nx

      lastid_mh=histos_mh(ind)%id
      lastind_mh=ind

      histos_mh(ind)%ktouched=2

9999  continue
      return
      end
