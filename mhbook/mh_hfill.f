*CMZ :  4.00/14 15/12/2021  11.54.34  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_hfill(id,x,y,w)

      use mhbook_mod
      implicit none

      double precision x,w,y
      integer mh_exists,id,ind,ihkind


      if (lastid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in mh_hfill: Non-existing histogram",id
          goto 9999 !return
        endif
      else
        ind=lastind_mh
        if (histos_mh(ind)%ny.gt.0) then
          ihkind=2
        else
          ihkind=1
        endif
      endif

      if (ihkind.eq.1) then
        call mh_fill1(id,x,w)
      else
        call mh_fill2(id,x,y,w)
      endif

9999  continue
      return
      end
