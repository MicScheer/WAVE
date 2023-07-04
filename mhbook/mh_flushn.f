*CMZ :  4.00/14 19/12/2021  15.38.08  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_flushn(id)

      use mhbook_mod
      implicit none

      integer mh_exists,id,ind,ihkind,i


      if (lastnid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in mh_flushn: Non-existing histogram",id
          goto 9999 !return
        endif
      else
        ind=lastnind_mh
      endif

      if (tups_mh(ind)%lun.eq.0) then
        open(newunit=tups_mh(ind)%lun,status='scratch',form='unformatted')
      endif

      do i=1,tups_mh(ind)%neve
        write(tups_mh(ind)%lun) tups_mh(ind)%eve(1:tups_mh(ind)%nvar,i)
      enddo

      tups_mh(ind)%neve=0

      lastnid_mh=id
      lastnind_mh=ind

9999  continue
      return
      end
