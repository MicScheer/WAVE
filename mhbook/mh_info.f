*CMZ :  4.00/14 15/12/2021  08.16.01  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_info(id)

      use mhbook_mod
      implicit none

      integer mh_exists,id,ind,ihkind,nx


      if (lastid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in mh_info: Non-existing histogram",id
          goto 9999 !return
        endif
      else
        ind=lastind_mh
      endif

      !call mh_update(id)

      print*,""
      print*,"id:",id
      print*,"title: ",trim(histos_mh(ind)%title)

      nx=histos_mh(ind)%nx

      if (histos_mh(ind)%ny.le.0) then
        print*,"number of entries:",histos_mh(ind)%nentries
        print*,""
        print*,"nx,underflow,overflow:",nx,histos_mh(ind)%nunder,histos_mh(ind)%nover
        print*,"underflow,overflow (weighted):",
     &    histos_mh(ind)%channels(2,nx+1,1),histos_mh(ind)%channels(2,nx+2,1)
        print*,"xmin,xmax,dx:",sngl(histos_mh(ind)%xmin),sngl(histos_mh(ind)%xmax),sngl(histos_mh(ind)%dx)
        print*,"xmean,xrms:",sngl(histos_mh(ind)%xmean),sngl(histos_mh(ind)%xrms)
        print*,"hmean,hrms:",sngl(histos_mh(ind)%hmean),sngl(histos_mh(ind)%hrms)
        print*,"hsum,hint:",sngl(histos_mh(ind)%hsum),sngl(histos_mh(ind)%hint)
      else
        print*,"Entries and under/over-flows:"
        print*,histos_mh(ind)%nou,histos_mh(ind)%nom,histos_mh(ind)%noo
        print*,histos_mh(ind)%num,histos_mh(ind)%nentries,histos_mh(ind)%nom
        print*,histos_mh(ind)%nuu,histos_mh(ind)%nmu,histos_mh(ind)%nuo
        print*,""
        print*,"Entries and under/over-flows (weighted):"
        print*,histos_mh(ind)%wnou,histos_mh(ind)%wnom,histos_mh(ind)%wnoo
        print*,histos_mh(ind)%wnum,histos_mh(ind)%hsum,histos_mh(ind)%wnom
        print*,histos_mh(ind)%wnuu,histos_mh(ind)%wnmu,histos_mh(ind)%wnuo
        print*,""
        print*,"xmin,xmax,dx:",sngl(histos_mh(ind)%xmin),sngl(histos_mh(ind)%xmax),sngl(histos_mh(ind)%dx)
        print*,"xmean,xrms:",sngl(histos_mh(ind)%xmean),sngl(histos_mh(ind)%xrms)
        print*,"ymin,ymax,dy:",sngl(histos_mh(ind)%ymin),sngl(histos_mh(ind)%ymax),sngl(histos_mh(ind)%dy)
        print*,"ymean,yrms:",sngl(histos_mh(ind)%ymean),sngl(histos_mh(ind)%yrms)
        print*,"hmean,hrms:",sngl(histos_mh(ind)%hmean),sngl(histos_mh(ind)%hrms)
        print*,"hsum,hint:",sngl(histos_mh(ind)%hsum),sngl(histos_mh(ind)%hint)
      endif

9999  continue


      return
      end
