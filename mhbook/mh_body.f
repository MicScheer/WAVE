*CMZ :  4.00/14 31/12/2021  11.03.05  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_body(id)

      use mhbook_mod
      implicit none

      integer mh_exists,id,ind,ihkind


      if (lastid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in mh_body: Non-existing histogram",id
          goto 9999 !return
        endif
      else
        ind=lastind_mh
      endif

      lastid_mh=id
      lastind_mh=ind

9999  continue
      return
      end
