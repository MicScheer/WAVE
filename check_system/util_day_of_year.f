*CMZ : 00.00/06 19/09/2007  14.34.57  by  Michael Scheer
*CMZ : 00.00/02 14/08/2006  13.09.52  by  Michael Scheer
*-- Author :    Michael Scheer   14/08/2006
      subroutine util_day_of_year(idate,iyear,mon,iday,iday_of_year,mode,istat)

      implicit none

      integer idate,iyear,mon,iday,mode,i,iday_of_year,istat,idaymon(12),
     &  montag

      data idaymon/31,28,31,30,31,30,31,31,30,31,30,31/

c mode=0:
c Input: idate (yyyy.mm.dd)
c Output: iday_of_year

c mode=1:
c Input: iyear, mon, iday
c Output: idate (yyyymmdd), iday_of_year

c mode=2:
c Input: iday_of_year, yyyy in idate
c Output: idate (yyyymmdd), iyear, mon, iday

      istat=0

      if (mode.eq.0) then

        if (idate/10000.lt.1000) then
          istat=-1
          return
        endif

        iday=0
        iyear=idate/10000

        if (
     &      (mod(iyear,4).eq.0.and.mod(iyear,100).ne.0)
     &      .or.mod(iyear,400).eq.0
     &      ) then
          idaymon(2)=29
        endif

        mon=(idate-iyear*10000)/100
        iday=idate-iyear*10000-mon*100

        iday_of_year=iday
        do i=1,mon-1
          iday_of_year=iday_of_year+idaymon(i)
        enddo

        if (mon.lt.1.or.mon.gt.12) then
          istat=-1
          return
        endif

        if (iday.lt.1.or.iday.gt.idaymon(mon)) then
          istat=-1
          return
        endif

      else if (mode.eq.1) then

        if (iyear/1000.lt.1) then
          istat=-1
          return
        endif

        idate=iyear*10000+mon*100+iday

        iday_of_year=iday
        do i=1,mon-1
          iday_of_year=iday_of_year+idaymon(i)
        enddo

      else if (mode.eq.2) then

        if (
     &      (mod(iyear,4).eq.0.and.mod(iyear,100).ne.0)
     &      .or.mod(iyear,400).eq.0
     &      ) then
          idaymon(2)=29
        endif

        iday=0
        mon=1

        montag=0
        do i=1,12
          iday=iday+idaymon(i)
          if (iday_of_year.gt.iday) then
            mon=mon+1
            montag=montag+idaymon(i)
          else
            iday=iday_of_year-montag
            goto 9
          endif
        enddo

9       idate=iyear*10000+mon*100+iday

      else
        istat=1
      endif

9999  return
      end
