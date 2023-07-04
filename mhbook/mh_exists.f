*CMZ :  4.00/14 31/12/2021  10.57.48  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      function mh_exists(id,ihkind)

      use mhbook_mod
      implicit none

      integer mh_exists,id,i,ihkind


      mh_exists=0

      if (nalloc_mh.le.0.or.nhist_mh+nntup_mh.le.0) then
        goto 9999 !return
      endif

      if (id.eq.lastid_mh.and.lastind_mh.gt.0.and.lastind_mh.le.nhist_mh) then
        if (histos_mh(lastind_mh)%id.eq.id) then
          mh_exists=lastind_mh
          if (histos_mh(lastind_mh)%ny.gt.0) then
            ihkind=2
          else
            ihkind=1
          endif
          goto 9999
        endif
      else if (id.eq.lastnid_mh.and.lastnind_mh.gt.0.and.
     &    lastnind_mh.le.nntup_mh) then
        if (tups_mh(lastnind_mh)%id.eq.id) then
          mh_exists=lastnind_mh
          ihkind=3
          goto 9999
        endif
      endif

      lastid_mh=0
      lastind_mh=0
      ihkind=0

      do i=1,nhist_mh
        if (histos_mh(i)%id.eq.id) then
          mh_exists=i
          if (histos_mh(i)%ny.gt.0) then
            ihkind=2
          else
            ihkind=1
          endif
          lastid_mh=id
          lastind_mh=i
          goto 9999 !return
        endif
      enddo

      do i=1,nntup_mh
        if (tups_mh(i)%id.eq.id) then
          mh_exists=i
          ihkind=3
          lastnid_mh=id
          lastnind_mh=i
          goto 9999 !return
        endif
      enddo

9999  continue
      return
      end
