*CMZ : 00.00/15 24/09/2012  12.56.10  by  Michael Scheer
*-- Author :    Michael Scheer   24/09/2012
      subroutine util_fwhm(wfwhmx,wfwhmy,wfwhm,ni,ne)

c +PATCH,//UTIL/FOR
c +DECK,util_fwhm.

      implicit none

      double precision wfwhmx(ne),wfwhmy(ne),wfwhm
      double precision wmxy,wm,wp,wmxy2
      integer i,ni,ne,nm,np,nmx

      wmxy=-1.0d30

      do i=ni,ne
        if (wfwhmy(i).gt.wmxy) then
          nmx=i
          wmxy=wfwhmy(i)
        endif
      enddo

      wmxy2=wmxy/2.0d0

      nm=ni
      do i=ni,nmx
        if (wfwhmy(i).le.wmxy2) then
          nm=i
        endif
      enddo

      np=ne
      do i=ne,nmx,-1
        if (wfwhmy(i).le.wmxy2) then
          np=i
        endif
      enddo

      if (wfwhmy(nm).gt.wmxy2) then
        wm=wfwhmx(nm)
      else
        wm=wfwhmx(nm)
     &    +(wfwhmx(nm+1)-wfwhmx(nm))/(wfwhmy(nm+1)-wfwhmy(nm))
     &    *(wmxy2-wfwhmy(nm))
      endif

      if (wfwhmy(np).gt.wmxy2) then
        wm=wfwhmx(nm)
      else
        wp=wfwhmx(np-1)
     &    +(wfwhmx(np)-wfwhmx(np-1))/(wfwhmy(np)-wfwhmy(np-1))
     &    *(wmxy2-wfwhmy(np-1))
      endif

      wfwhm=abs(wp-wm)

      return
      end
