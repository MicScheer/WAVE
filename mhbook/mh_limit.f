*CMZ :  4.00/14 14/12/2021  18.01.13  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_limit(limit)

      use mhbook_mod
      implicit none

      integer*8 iutil_memory_size
      integer limit,istat,i,nx,ny


      if (memsize_mh.le.0) memsize_mh=iutil_memory_size(memstart_mh)

      if (limit.le.0) then
        if (nalloc_mh.gt.0) then
          deallocate(histos_mh,tups_mh,stat=istat)
        endif
        nalloc_mh=-1
        goto 9999 !return
      endif

      if (nalloc_mh.le.0) then

        allocate(histos_mh(limit),tups_mh(limit))

      else if (limit.gt.nalloc_mh) then

        allocate(histos_copy_mh(nalloc_mh))

        do i=1,nalloc_mh
          nx=histos_mh(i)%nx
          ny=histos_mh(i)%ny
          if (ny.gt.0) then
            allocate(histos_copy_mh(i)%channels(5,nx+2,ny+2))
          else
            allocate(histos_copy_mh(i)%channels(5,nx+2,1))
          endif
          histos_copy_mh(i)=histos_mh(i)
        enddo
        deallocate(histos_mh)
        allocate(histos_mh(limit))
        do i=1,nalloc_mh
          nx=histos_copy_mh(i)%nx
          ny=histos_copy_mh(i)%ny
          if (ny.gt.0) then
            allocate(histos_mh(i)%channels(5,nx+2,ny+2))
          else
            allocate(histos_mh(i)%channels(5,nx+2,1))
          endif
          histos_mh(i)=histos_copy_mh(i)
        enddo
        deallocate(histos_copy_mh)

        allocate(tups_copy_mh(nalloc_mh))
        do i=1,nalloc_mh
          nx=tups_mh(i)%nalloc
          allocate(tups_copy_mh(i)%eve(nvarp_mh,nx))
          tups_copy_mh(i)=tups_mh(i)
        enddo
        deallocate(tups_mh)
        allocate(tups_mh(limit))
        do i=1,nalloc_mh
          nx=tups_copy_mh(i)%nalloc
          allocate(tups_mh(i)%eve(nvarp_mh,nx))
          tups_mh(i)=tups_copy_mh(i)
        enddo
        deallocate(tups_copy_mh)

      endif

      nalloc_mh=limit

      deallocate(idhis_mh,idntup_mh,indhis_mh,indntup_mh,
     &  idhbuff_mh,idnbuff_mh,stat=istat)
      allocate(indhis_mh(nalloc_mh),indntup_mh(nalloc_mh),
     &  idhis_mh(nalloc_mh),idntup_mh(nalloc_mh),
     &  idhbuff_mh(nalloc_mh),idnbuff_mh(nalloc_mh))

9999  continue
      return
      end
