*CMZ :  4.00/14 31/12/2021  15.13.30  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_bookn(id,chtit,nvar,chvar,nalloc)

      use mhbook_mod
      implicit none

      integer mh_exists,id,nalloc,nvar,i,ihkind

      character(*) chvar(*)
      character(*) chtit


      if (nvar.gt.nvarp_mh) then
        print*
        print*,"*** Error in mh_bookn: Too many variables, limit is ",nvarp_mh
        print*,id,trim(chtit)
        print*
        goto 9999 !return
      endif

      if (nalloc_mh.le.0) then
        call mh_limit(1000)
      endif

      if (mh_exists(id,ihkind).ne.0) then
        print*
        print*,"*** Error in mh_bookn: Already existing Ntuple ",id
        print*,trim(chtit)
        print*
        goto 9999 !return
      endif

      if (nalloc.lt.1) then
        nntup_mh=nntup_mh+1
        tups_mh(nntup_mh)%id=id
        tups_mh(nntup_mh)%nalloc=-1
        tups_mh(nntup_mh)%nvar=-1
        tups_mh(nntup_mh)%title=chtit
        tups_mh(nntup_mh)%chpath=chvar(1)
        goto 9999 !return
      endif

      if (nvar.lt.1) then
        print*,"*** Error in mh_bookn: nvar < 1",
     &    id, trim(chtit)
        goto 9999 !return
      endif

      nntup_mh=nntup_mh+1

      if (nntup_mh.gt.nalloc_mh) then
        print*,"*** mh_bookn: Increasing max. number of Ntuples ***"
        call mh_limit(nntup_mh*2)
      endif

      do i=1,nvar
        tups_mh(nntup_mh)%chvar(i)=trim(adjustl(chvar(i)))
        tups_mh(nntup_mh)%varm(1,i)=1.0d30
        tups_mh(nntup_mh)%varm(2,i)=-1.0d30
      enddo

      memntuptot_mh=memntuptot_mh+nalloc*nvar

      if (memntuptot_mh.gt.memsize_mh*memsizemax_mh/100) then
        print*,"*** mh_bookn: Switching to disk mode ***"
        open(newunit=tups_mh(nntup_mh)%lun,status='scratch',form='unformatted')
        allocate(tups_mh(nntup_mh)%eve(nvarp_mh,nflushp_mh))
        tups_mh(nntup_mh)%nalloc=nflushp_mh
      else
        allocate(tups_mh(nntup_mh)%eve(nvarp_mh,nalloc))
        tups_mh(nntup_mh)%nalloc=nalloc
      endif

      tups_mh(nntup_mh)%id=id
      tups_mh(nntup_mh)%nvar=nvar
      tups_mh(nntup_mh)%title=chtit
      tups_mh(nntup_mh)%chpath=chdir_mh

      if (lun_index_mh.eq.0) open(newunit=lun_index_mh,file='mh_book.lis')
      nnbooked_mh=nnbooked_mh+1
      if (nntup_mh.gt.nntupmax_mh) nntupmax_mh=nntup_mh
      write(lun_index_mh,*)nnbooked_mh,", Ntuple, ",id," ,",trim(chtit)

      lastnid_mh=id
      lastnind_mh=nntup_mh

9999  continue
      return
      end
