*CMZ : 00.00/07 28/03/2008  09.41.32  by  Michael Scheer
*-- Author :    Michael Scheer   25/03/2008
*CMZ :          25/03/2008  12.30.00  by  Michael Scheer
*-- Author :    Michael Scheer   25/03/2008
      subroutine util_soft_edge(xin,xlen,ypeak,xpos,wedge,yout,mode,ifail)

c calcutes yout as Fermi-edge shaped block
c mode: negative: Initialization
c mode: abs(mode)=1: xlen=yint/ypeak
c mode: abs(mode)=2: ypeak=yint/xlen

      implicit none

      double precision x,xin,xlen,ypeak,xpos,wedge,yout,yint,
     &  ay1,ay2,by1,by2,explw,xlen2,xlenw,xw,rwedge

      integer mode,ifail,ical

      ifail=0
      yout=0.0d0

      x=xin-xpos

      if (ypeak.eq.0.0d0.or.xlen.eq.0.0d0) then
        return
      endif

      if (ical.eq.0.or.mode.lt.0) then
        rwedge=wedge/xlen !anders als in magseq, nachträglich eingefügt...
        xlen2=xlen/2.0d0
        xlenw=xlen*rwedge
        if (rwedge.eq.0.0d0.or.xlenw.gt.70.d0) then
          xw=0.0d0
        else
          xw=rwedge
        endif
        if (rwedge.lt.0.0d0) then
          ifail=-1
          return
        else if (xw.gt.0.0d0) then
          explw=exp(xlenw)
          if (explw.le.1.0d0) then
            ifail=-1
            return
          endif
          yint=explw*ypeak/(explw-1.0d0)*xlen
          if (abs(mode).eq.1) then
            xlen=yint/ypeak
          else if (abs(mode).eq.2) then
            if (xlen.eq.0.0d0) then
              ifail=-3
              return
            endif
            ypeak=yint/xlen
          else !mode
            ifail=-3
            return
          endif !mode
        endif !rwedge.eq.0.0d0

        ical=1

      endif !(ical.eq.0.or.mode.lt.0)

      if (xw.eq.0.0d0.and.abs(x).le.xlen2)then

          yout=ypeak

      else if (xlen*ypeak.ne.0.0d0) then

        ay1=(+x-xlen2)*rwedge
        ay2=(-x-xlen2)*rwedge

        if (ay1.gt.70.0d0) then
          by1=0.0d0
        else if (ay1.lt.-70.) then
          by1=1.0d0
        else
          by1=1.0d0/(1.0d0+dexp(ay1))
        endif

        if (ay2.gt.70.0d0) then
          by2=0.0d0
        else if (ay2.lt.-70.0d0) then
          by2=1.0d0
        else
          by2=1.0d0/(1.0d0+dexp(ay2))
        endif

        yout=ypeak*by1*by2

      endif

      return
      end
