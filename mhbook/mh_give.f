*CMZ :  4.00/14 14/12/2021  18.18.26  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_give(id,title,nx,xmin,xmax,ny,ymin,ymax,ind,idum)

      use mhbook_mod
      implicit none

      double precision xmin,xmax,ymin,ymax
      integer mh_exists,id,ind,ihkind,nx,ny,idum
      character(*) title


      if (lastid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in mh_give: Non-existing histogram",id
          goto 9999 !return
        endif
      else
        ind=lastind_mh
      endif

      title=trim(histos_mh(ind)%title)

      nx=histos_mh(ind)%nx
      ny=histos_mh(ind)%ny

      xmin=histos_mh(ind)%xmin
      ymin=histos_mh(ind)%ymin
      xmax=histos_mh(ind)%xmax
      ymax=histos_mh(ind)%ymax

      idum=0

      lastid_mh=id
      lastind_mh=ind

9999  continue
      return
      end
