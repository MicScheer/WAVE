*CMZ :  4.00/14 16/12/2021  15.05.33  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_copy(id,id2,title)

      use mhbook_mod
      implicit none

      integer mh_exists,id,id2,ind,ind2,ihkind
      character(*) title


      ind2=mh_exists(id2,ihkind)
      if (ind2.ne.0) then
        print*
        print*,"*** Error in mh_book1: Already existing histogram ",id2
        print*, trim(histos_mh(ind2)%title)
        print*
        goto 9999 !return
      endif

      if (lastid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in mh_copy: Non-existing histogram",id
          goto 9999 !return
        endif
      else
        ind=lastind_mh
      endif

      nhist_mh=nhist_mh+1
      ind2=nhist_mh

      allocate(histos_mh(ind2)%channels(4,histos_mh(ind)%nx+2,histos_mh(ind)%ny+2))
      histos_mh(ind2)=histos_mh(ind)

      histos_mh(ind2)%id=id2
      histos_mh(ind2)%title=trim(title)

      if (histos_mh(id)%ny.le.0) then
        nhist1_mh=nhist1_mh+1
      else
        nhist2_mh=nhist2_mh+1
      endif

      histos_mh(ind2)%ktouched=1

      lastid_mh=id2
      lastind_mh=ind2

9999  continue
      return
      end
