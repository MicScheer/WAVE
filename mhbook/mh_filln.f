*CMZ :  4.00/15 11/03/2022  10.47.52  by  Michael Scheer
*CMZ :  4.00/14 21/12/2021  12.47.20  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_filln(id,var)

      use mhbook_mod
      implicit none

      integer mh_exists,id,ind,ihkind,i
      double precision var(*)


      if (lastnid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in mh_filln: Non-existing Ntuple",id
          goto 9999 !return
        endif
      else
        ind=lastnind_mh
      endif

      if (tups_mh(ind)%nvar.le.0) goto 9999

      if (tups_mh(ind)%neve.ge.nflushp_mh.or.
     &    tups_mh(ind)%neve.ge.tups_mh(ind)%nalloc) then
c        print*,"*** mh_filln: Flushing Ntuple",id,", ",trim(tups_mh(ind)%title)
        call mh_flushn(id)
      endif

      tups_mh(ind)%nentries=tups_mh(ind)%nentries+1
      tups_mh(ind)%neve=tups_mh(ind)%neve+1

      do i=1,tups_mh(ind)%nvar
        if (var(i).lt.tups_mh(ind)%varm(1,i)) tups_mh(ind)%varm(1,i)=var(i)
        if (var(i).gt.tups_mh(ind)%varm(2,i)) tups_mh(ind)%varm(2,i)=var(i)
        tups_mh(ind)%eve(i,tups_mh(ind)%neve)=var(i)
      enddo

      tups_mh(ind)%ktouched=1

      lastnid_mh=id
      lastnind_mh=ind

9999  continue
      return
      end
