*CMZ :  3.02/03 03/11/2014  10.42.03  by  Michael Scheer
*CMZ :  2.68/05 02/10/2012  09.01.01  by  Michael Scheer
*CMZ :  2.66/20 06/07/2011  11.58.18  by  Michael Scheer
*CMZ :  2.66/12 24/06/2010  12.50.52  by  Michael Scheer
*CMZ :  2.66/08 12/03/2010  15.02.37  by  Michael Scheer
*CMZ :  2.66/07 10/03/2010  14.08.46  by  Michael Scheer
*-- Author :    Michael Scheer   24/02/2010
      subroutine ubunch(xbunch,ybunch,zbunch,ypbunch,zpbunch,gambunch,
     &  dtphase)
*KEEP,gplhint.
*KEND.

      use bunchmod

      implicit none

      double precision dtphase,xbunch,ybunch,zbunch,
     &  ypbunch,zpbunch,gambunch

      integer ical

      data ical/0/

      if (ical.eq.0) then

        xbunch=xbunch
        ybunch=ybunch
        zbunch=zbunch
        ypbunch=ypbunch
        zpbunch=zpbunch
        gambunch=gambunch
        dtphase=dtphase

        ical=1
      endif !ical

      return
      end

