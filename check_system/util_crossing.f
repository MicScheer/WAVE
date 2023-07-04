*CMZ : 00.00/02 05/07/2006  15.48.36  by  Michael Scheer
*-- Author :    Michael Scheer   30/06/2006
      subroutine util_crossing(ax1,ay1,ax2,ay2,bx1,by1,bx2,by2,px,py,istat)

c calculates crossing point of lines
c lines a and b are given two points, (px,py) is crossing point if it exists.
c istat=0 means point exists between the points
c istat=1 means point exists outside the points
c istat=-1 means point does not exist

      double precision ax1,ay1,ax2,ay2,bx1,by1,bx2,by2,px,py
      double precision slope1,slope2
      integer istat

      istat=-1

      if (ax1.eq.ax2.and.bx1.eq.bx2) return

      if (ax1.ne.ax2.and.bx1.ne.bx2) then

        slope1=(ay2-ay1)/(ax2-ax1)
        slope2=(by2-by1)/(bx2-bx1)

        if (slope1.eq.slope2) return

c Normal case, both slopes not infinit and different

c solve(ay1+slope1*(px-ax1)=by1+slope2*(px-bx1),px)

        px=(-(bx1*slope2-by1+ay1-ax1*slope1))/(slope1-slope2)
        py = ay1 + slope1*(px-ax1)

        if (
     &      (px.ge.ax1.and.px.le.ax2.or.px.ge.ax2.and.px.le.ax1).and.
     &      (px.ge.bx1.and.px.le.bx2.or.px.ge.bx2.and.px.le.bx1)
     &      ) then
            istat=0
          else
            istat=1
        endif

      else if (ax1.ne.ax2) then

        slope1=(ay2-ay1)/(ax2-ax1)

        px = bx1
        py = ay1 + slope1*(px-ax1)

        if (
     &      (px.ge.ax1.and.px.le.ax2.or.px.ge.ax2.and.px.le.ax1).and.
     &      (py.ge.by1.and.py.le.by2.or.px.ge.by2.and.px.le.by1)
     &   ) then
            istat=0
          else
            istat=1
        endif

      else if (bx1.ne.bx2) then

        slope2=(by2-by1)/(bx2-bx1)

        px = ax1
        py = by1 + slope2*(px-bx1)

        if (
     &      (px.ge.bx1.and.px.le.bx2.or.px.ge.bx2.and.px.le.bx1).and.
     &      (py.ge.ay1.and.py.le.ay2.or.py.ge.ay2.and.py.le.ay1)
     &      ) then
            istat=0
          else
            istat=1
        endif

      endif !slope1

      return
      end
