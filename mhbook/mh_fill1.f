*CMZ :  4.00/14 16/12/2021  15.05.33  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_fill1(id,x,w)

      use mhbook_mod
      implicit none

      double precision x,w,hn,hsum,hsum2,hmean,hrms,hsumtot
      integer mh_exists,id,ind,ix,ihkind


      if (lastid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in mh_fill1: Non-existing histogram",id
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

      if (w.lt.histos_mh(ind)%hmin) histos_mh(ind)%hmin=w
      if (w.gt.histos_mh(ind)%hmax) histos_mh(ind)%hmax=w

      hsumtot=histos_mh(ind)%hsum+w
      histos_mh(ind)%hsum=hsumtot
      histos_mh(ind)%hint=hsumtot*histos_mh(ind)%dx
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

      histos_mh(ind)%ktouched=1

      lastid_mh=id
      lastind_mh=ind

9999  continue
      return
      end
