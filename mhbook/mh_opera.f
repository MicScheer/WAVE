*CMZ :  4.00/14 16/12/2021  15.05.33  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_opera(id1,cop,id2,id3,c1,c2)

      use mhbook_mod
      implicit none

      double precision xmin1,xmax1,xmin2,xmax2,ymin1,ymax1,ymin2,ymax2,c1,c2
      integer mh_exists,ind1,ind2,ind3,id1,id2,id3,nx1,ny1,nx2,ny2,idum,ix,iy,
     &  ihkind
      character(2048) title
      character cop


      call mh_give(id2,title,nx2,xmin2,xmax2,ny2,ymin2,ymax2,ind2,idum)
      call mh_give(id1,title,nx1,xmin1,xmax1,ny1,ymin1,ymax1,ind1,idum)

      if (ind1.le.0) then
        print*,"*** Error in mh_opera: Non-existing histogram",id1
        goto 9999 !return
      endif

      if (ind2.le.0) then
        print*,"*** Error in mh_opera: Non-existing histogram",id2
        goto 9999 !return
      endif

      if (nx1.ne.nx2.or.xmin1.ne.xmin2.or.xmax1.ne.xmax2.or.
     &    ny1.ne.ny2.or.ymin1.ne.ymin2.or.ymax1.ne.ymax2) then
        write(6,*)'*** Error in mh_opera: Histograms do not match, ID1, ID2=',
     &    id1,id2
        return
      endif

      if (
     &    cop.ne.'+' .and. cop.ne.'-'.and. cop.ne.'*' .and. cop.ne.'/') then
        write(6,*)'*** Error in mh_opera: Option ',cop(1:len_trim(cop)),
     &    'not yet implemented! Only +,-,*,/ are available.'
        goto 9999
      endif

      ind3=mh_exists(id3,ihkind)
      if (ind3.le.0) then
        call mh_copy(id1,id3,histos_mh(ind1)%title)
      endif

      ind3=lastind_mh

      histos_mh(ind3)%nentries=
     &  histos_mh(ind1)%nentries+histos_mh(ind2)%nentries

      if (cop.eq.'+') then
        histos_mh(ind3)%channels=c1*histos_mh(ind1)%channels+
     &    c2*histos_mh(ind2)%channels
      else if (cop.eq.'-') then
        histos_mh(ind3)%channels=c1*histos_mh(ind1)%channels-
     &    c2*histos_mh(ind2)%channels
      else if (cop.eq.'*') then
        histos_mh(ind3)%channels=c1*histos_mh(ind1)%channels*
     &    c2*histos_mh(ind2)%channels
      else if (cop.eq.'/') then
        if (c2.eq.0.0d0) then
          histos_mh(ind3)%channels=0.0d0
        else
          do iy=1,ny2
            do ix=1,nx2
              if(histos_mh(ind2)%channels(2,ix,iy).ne.0.0d0) then
                histos_mh(ind3)%channels=c1*histos_mh(ind1)%channels/
     &            c2*histos_mh(ind2)%channels
              else
                histos_mh(ind3)%channels(:,ix,iy)=0.0d0
              endif
            enddo
          enddo
        endif
      endif

      !call mh_update(id3)

      histos_mh(ind3)%ktouched=1

9999  continue
      return
      end
