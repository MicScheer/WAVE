*CMZ :  4.00/14 20/12/2021  16.03.49  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_rout(id)

      use mhbook_mod
      implicit none

      integer mh_exists,id,ind,ihkind,i,k
      integer,save :: ical=0


      if (id.eq.0) then

        call mh_sort
        do i=1,nhist_mh
          k=indhis_mh(i)
          if (histos_mh(k)%ny.gt.0) then
            call mh_h2out(k)
          else
            call mh_h1out(k)
          endif
        enddo
        do i=1,nntup_mh
          call mh_nout(indntup_mh(i))
        enddo

      else

        if (lastid_mh.ne.id.and.lastind_mh.ne.id) then
          ind=mh_exists(id,ihkind)
          if (ind.le.0) then
            print*,"*** Error in mh_rout: Non-existing histogram",id
            goto 9999 !return
          endif
        else if (lastid_mh.eq.id) then
          ind=lastind_mh
          if (histos_mh(ind)%ny.gt.0) then
            ihkind=2
          else
            ihkind=1
          endif
        else
          ind=lastnind_mh
          ihkind=3
        endif

        if (ihkind.eq.3) then
          call mh_nout(ind)
          lastnid_mh=id
          lastnind_mh=ind
        else if (ihkind.eq.2) then
          call mh_h2out(ind)
          lastid_mh=id
          lastind_mh=ind
        else
          call mh_h1out(ind)
          lastid_mh=id
          lastind_mh=ind
        endif

      endif !id.eq.0

9999  continue
      return
      end
