*CMZ :  4.00/14 31/12/2021  15.23.10  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      function  hij_mh(id,i,j)

      use mhbook_mod
      implicit none

      double precision hij_mh
      integer mh_exists,id,ind,i,j,ihkind


      if (lastid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in hij_mh: Non-existing histogram",id
          goto 9999 !return
        endif
        if (ihkind.ne.2) then
          print*,"*** Error in hij_mh: Identifier does not belong to a 2d-histogram",id
          goto 9999 !return
        endif
      else
        ind=lastind_mh
      endif

      hij_mh=histos_mh(ind)%channels(2,i,j)

      lastid_mh=id
      lastind_mh=ind

9999  continue
      return
      end
