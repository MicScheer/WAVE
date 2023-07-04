*CMZ :  4.00/14 31/12/2021  15.23.16  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      function  hi_mh(id,i)

      use mhbook_mod
      implicit none

      double precision hi_mh
      integer mh_exists,id,ind,i,ihkind


      if (lastid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in hi_mh: Non-existing histogram",id
          goto 9999 !return
        endif
        if (ihkind.ne.1) then
          print*,"*** Error in hi_mh: Identifier does not belong to a 1d-histogram",id
          goto 9999 !return
        endif
      else
        ind=lastind_mh
      endif

      hi_mh=histos_mh(ind)%channels(2,i,1)

      lastid_mh=id
      lastind_mh=ind

9999  continue
      return
      end
