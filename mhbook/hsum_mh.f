*CMZ :  4.00/14 30/12/2021  12.55.40  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      function hsum_mh(id)

      use mhbook_mod
      implicit none

      double precision hsum_mh
      integer id,ihkind,ind,mh_exists


      if (nalloc_mh.le.0.or.nhist_mh.eq.0) then
        return
      endif

      if (id.eq.lastid_mh) then
        hsum_mh=histos_mh(lastind_mh)%hsum
      else
        ind=mh_exists(id,ihkind)
        lastid_mh=histos_mh(ind)%id
        lastind_mh=ind
        hsum_mh=histos_mh(ind)%hsum
      endif

      return
      end
