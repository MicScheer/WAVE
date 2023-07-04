*CMZ :  4.00/14 30/12/2021  12.38.06  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      function hmax_mh(id)

      use mhbook_mod
      implicit none

      double precision hmax_mh
      integer id,ihkind,ind,mh_exists


      if (nalloc_mh.le.0.or.nhist_mh.le.0) then
        return
      endif

      if (id.eq.lastid_mh) then
        hmax_mh=histos_mh(lastind_mh)%hmax
      else
        ind=mh_exists(id,ihkind)
        lastid_mh=histos_mh(ind)%id
        lastind_mh=ind
        hmax_mh=histos_mh(ind)%hmax
      endif

      return
      end
