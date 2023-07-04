*CMZ : 00.00/06 21/09/2007  11.11.02  by  Michael Scheer
*CMZ : 00.00/02 14/08/2006  13.09.52  by  Michael Scheer
*-- Author :    Michael Scheer   14/08/2006
      subroutine util_seconds(
     &  iyear,imon,iday,ihour,imin,isec,isectot,mode,istat)

      implicit none

      integer idaytot,iyear,imon,iday,ihour,imin,isec,isectot,
     &  mode,i,istat,idaymon(12),
     &  iday_of_year,n29

      data idaymon/31,28,31,30,31,30,31,31,30,31,30,31/

c mode >= 0:
c Input: iyear, imon, iday, ihour, imin, isec
c Output: seconds since 1.1.1970

c mode < 0:
c Input: seconds since 1.1.1970
c Output: iyear, imon, iday, ihour, imin, seconds

      istat=0

      if (mode.ge.0) then

        if (mod(iyear-1900,4).eq.0..and.iyear.ne.1900) then
          idaymon(2)=29
        endif

        if (iyear.lt.1970.or.iyear.gt.2038
     &      .or.iday.lt.1.or.imon.lt.1.or.imon.gt.12.or.
     &      iday.gt.idaymon(imon).or.ihour.gt.24
     &      .or.
     &      imin.gt.60.or.
     &      isec.gt.60) then
          istat=-1
          return
        endif

c Anzahl der vergangenen Schaltjahre

        n29=(iyear-1969)/4 !2000 war ein Schaltjahr!

        idaytot=365*(iyear-1970)+iday-1+n29

        do i=1,imon-1
          idaytot=idaytot+idaymon(i)
        enddo

        isectot=((idaytot*24+ihour)*60+imin)*60+isec

      else !mode

        idaytot=isectot/(24*3600)
        iyear=idaytot/365

        if (mod(iyear-2,4).eq.0) then
          idaymon(2)=29
        endif

        n29=(iyear-1)/4 !2000 war Schaltjahr
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

        iyear=iyear+1970
        isec=isectot-idaytot*24*3600

        ihour=isec/3600.0d0
        imin=isec/60
        imin=imin-ihour*60

        isec=isec-((ihour*60+imin)*60)

      endif

9999  return
      end
