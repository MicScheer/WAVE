*CMZ :  4.00/14 16/12/2021  15.05.33  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_fill2(id,x,y,w)

      use mhbook_mod
      implicit none

      double precision x,w,y,hn,hsumtot,hsum,hsum2,hmean,hrms
      integer mh_exists,id,ind,ihkind,ix,iy


      if (lastid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in mh_fill2: Non-existing histogram",id
          goto 9999 !return
        endif
      else
        ind=lastind_mh
      endif

      if (x.lt.histos_mh(ind)%xmin) then
        ix=histos_mh(ind)%nx+1
      else if (x.gt.histos_mh(ind)%xmax) then
        ix=histos_mh(ind)%nx+2
      else
        ix=int((x-histos_mh(ind)%xmin)/histos_mh(ind)%dx)+1
      endif

      if (y.lt.histos_mh(ind)%ymin) then
        iy=histos_mh(ind)%ny+1
      else if (y.gt.histos_mh(ind)%ymax) then
        iy=histos_mh(ind)%ny+2
      else
        iy=int((y-histos_mh(ind)%ymin)/histos_mh(ind)%dy)+1
      endif

      histos_mh(ind)%nentries=histos_mh(ind)%nentries+1

      hn=histos_mh(ind)%channels(1,ix,1)+1.0d0
      hsum=histos_mh(ind)%channels(2,ix,1)+w
      hsum2=histos_mh(ind)%channels(5,ix,1)+w*w
      hmean=hsum/hn
      hrms=sqrt(max(0.0d0,hsum2/hn-hmean**2))

      histos_mh(ind)%channels(1,ix,1)=hn
      histos_mh(ind)%channels(2,ix,1)=hsum
      histos_mh(ind)%channels(3,ix,1)=hmean
      histos_mh(ind)%channels(4,ix,1)=hrms
      histos_mh(ind)%channels(5,ix,1)=hsum2

      if (x.lt.histos_mh(ind)%xmin) then
        if (y.lt.histos_mh(ind)%ymin) then
          histos_mh(ind)%nuu=histos_mh(ind)%nuu+1
          histos_mh(ind)%wnuu=histos_mh(ind)%wnuu+w
        else if (y.gt.histos_mh(ind)%ymax) then
          histos_mh(ind)%nuo=histos_mh(ind)%nuo+1
          histos_mh(ind)%wnuo=histos_mh(ind)%wnuo+w
        else
          histos_mh(ind)%num=histos_mh(ind)%num+1
          histos_mh(ind)%wnum=histos_mh(ind)%wnum+w
        endif
      else if (x.gt.histos_mh(ind)%xmax) then
        if (y.lt.histos_mh(ind)%ymin) then
          histos_mh(ind)%nou=histos_mh(ind)%nou+1
          histos_mh(ind)%wnou=histos_mh(ind)%wnou+w
        else if (y.gt.histos_mh(ind)%ymax) then
          histos_mh(ind)%noo=histos_mh(ind)%noo+1
          histos_mh(ind)%wnoo=histos_mh(ind)%wnoo+w
        else
          histos_mh(ind)%nom=histos_mh(ind)%nom+1
          histos_mh(ind)%wnom=histos_mh(ind)%wnom+w
        endif
      else
        if (y.lt.histos_mh(ind)%ymin) then
          histos_mh(ind)%nmu=histos_mh(ind)%nmu+1
          histos_mh(ind)%wnmu=histos_mh(ind)%wnmu+w
        else if (y.gt.histos_mh(ind)%ymax) then
          histos_mh(ind)%nmo=histos_mh(ind)%nmo+1
          histos_mh(ind)%wnmo=histos_mh(ind)%wnmo+w
        endif
      endif

      if (w.lt.histos_mh(ind)%hmin) histos_mh(ind)%hmin=w
      if (w.gt.histos_mh(ind)%hmax) histos_mh(ind)%hmax=w

      hsumtot=histos_mh(ind)%hsum+w
      histos_mh(ind)%hsum=hsumtot
      histos_mh(ind)%hint=hsumtot*histos_mh(ind)%dx*histos_mh(ind)%dy
      histos_mh(ind)%hsum2=histos_mh(ind)%hsum2+w*w

      hsum=histos_mh(ind)%xsum+x*w
      hsum2=histos_mh(ind)%xsum2+x*x*w
      if (hsumtot.ne.0) then
        hmean=hsum/hsumtot
        hrms=sqrt(max(0.0d0,hsum2/hsumtot-hmean**2))
      else
        hmean=0.0d0
        hrms=0.0d0
      endif
      histos_mh(ind)%xsum=hsum
      histos_mh(ind)%xsum2=hsum2
      histos_mh(ind)%xmean=hmean
      histos_mh(ind)%xrms=hrms

      hsum=histos_mh(ind)%ysum+y*w
      hsum2=histos_mh(ind)%ysum2+y*y*w
      if (hsumtot.ne.0) then
        hmean=hsum/hsumtot
        hrms=sqrt(max(0.0d0,hsum2/hsumtot-hmean**2))
      else
        hmean=0.0d0
        hrms=0.0d0
      endif
      histos_mh(ind)%ysum=hsum
      histos_mh(ind)%ysum2=hsum2
      histos_mh(ind)%ymean=hmean
      histos_mh(ind)%yrms=hrms

      histos_mh(ind)%ktouched=1

      lastid_mh=id
      lastind_mh=ind

9999  continue
      return
      end
