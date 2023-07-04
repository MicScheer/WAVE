*CMZ :  4.00/14 19/12/2021  10.52.36  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_delete(id)

      use mhbook_mod
      implicit none

      integer mh_exists,id,ind,nx,ny,ihkind


      ind=mh_exists(id,ihkind)

      if (ihkind.le.2) then

        if (ind.eq.0) then
          print*,"*** Error in mh_delete: Non-existing histogram ",id
          goto 9999 !return
        endif

        if (nhist_mh.gt.1.and.ind.lt.nhist_mh) then
          nx=histos_mh(nhist_mh)%nx
          ny=histos_mh(nhist_mh)%ny
          deallocate(histos_mh(ind)%channels)
          if (ny.le.0) then
            allocate(histos_mh(ind)%channels(5,nx+2,1))
          else
            allocate(histos_mh(ind)%channels(5,nx+2,ny+2))
          endif
          histos_mh(ind)=histos_mh(nhist_mh)
          deallocate(histos_mh(nhist_mh)%channels)
        endif

        nhist1_mh=nhist1_mh-1
        nhist_mh=nhist_mh-1

c      else
c
c        stop "*** Error in mh_delete: Deleting of Ntuples not yet implemented ***"

      endif

      lastid_mh=0
      lastind_mh=0

9999  continue

      return
      end
