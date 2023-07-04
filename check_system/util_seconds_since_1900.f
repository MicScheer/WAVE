*CMZ : 00.00/19 07/09/2015  09.49.29  by  Michael Scheer
*CMZ : 00.00/06 21/09/2007  11.19.19  by  Michael Scheer
*CMZ : 00.00/02 14/08/2006  13.09.52  by  Michael Scheer
*-- Author :    Michael Scheer   14/08/2006
      subroutine util_seconds_since_1900(
     &  iyear,imon,iday,ihour,imin,dsec,ds1900,mode,istat)

      implicit none

      double precision dsec,ds1900

      integer idaytot,iyear,imon,iday,ihour,imin,mode,i,istat,idaymon(12),
     &  iday_of_year,n29

      data idaymon/31,28,31,30,31,30,31,31,30,31,30,31/

c mode >= 0:
c Input: iyear, imon, iday, ihour, imin, dsec
c Output: seconds_since_1900

c mode < 0:
c Input: seconds_since_1900
c Output: iyear, imon, iday, ihour, imin, ds1900

      istat=0

      if (mode.ge.0) then

        if (mod(iyear-1900,4).eq.0..and.iyear.ne.1900) then
          idaymon(2)=29
        endif

        if (iyear.lt.1900.or.iyear.gt.2058
     &      .or.imon.lt.1.or.imon.gt.12
     &      .or.iday.lt.1
     &      .or.iday.gt.idaymon(imon)
     &      .or.ihour.gt.24
     &      .or.
     &      imin.gt.60.or.
     &      dsec.gt.60.0d0) then
          istat=-1
          return
        endif

c Anzahl der vergangenen Schaltjahre

        n29=(iyear-1901)/4 !1900 war kein, aber 2000 war ein Schaltjahr

        idaytot=365*(iyear-1900)+iday-1+n29

        do i=1,imon-1
          idaytot=idaytot+idaymon(i)
        enddo

        ds1900=dble(((idaytot*24+ihour)*60)+imin)*60+dsec

      else !mode

        idaytot=int(ds1900/(24.0d0*3600.0d0+0.00001))
        iyear=idaytot/365

        if (mod(iyear,4).eq.0.and.iyear.ne.0) then
          idaymon(2)=29
        endif

        n29=(iyear-1)/4 !1900 war kein, aber 2000 war ein Schaltjahr
        iday_of_year=idaytot-iyear*365-n29

        iday=iday_of_year
        do imon=1,12
          iday=iday-idaymon(imon)
          if (iday.lt.0) then
            iday=iday+idaymon(imon)+1
            goto 9
          endif
        enddo

9       continue

        iyear=iyear+1900
        dsec=ds1900-dble(idaytot*24.0d0*3600.0d0)

        ihour=int(dsec/3600.0d0)
        imin=int(dsec/60.)
        imin=imin-ihour*60

        dsec=dsec-dble((ihour*60+imin)*60)

        if (ihour.eq.24) then
          ihour=0
          iday=iday+1
          if (iday.gt.idaymon(imon)) then
            iday=1
            imon=imon+1
            if (imon.gt.12) then
              imon=1
              iyear=iyear+1
            endif
          endif
        endif

      endif

      return
      end
