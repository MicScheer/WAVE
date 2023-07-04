*CMZ :  4.00/14 14/12/2021  18.18.26  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_type(id)

      use mhbook_mod
      implicit none

      double precision x,y,xmin,ymin,dx,dy
      integer mh_exists,id,ind,ihkind,ix,iy


      if (lastid_mh.ne.id) then
        ind=mh_exists(id,ihkind)
        if (ind.le.0) then
          print*,"*** Error in mh_type: Non-existing histogram",id
          goto 9999 !return
        endif
      else
        ind=lastind_mh
      endif

      call mh_info(id)

      if (histos_mh(ind)%ny.gt.0) then

        write(6,*) ""
        write(6,*) "  ix   iy     x            y             n           h            hmean        hrms"
        write(6,*) "------------------------------------------------------------------------------------------"

        xmin=histos_mh(ind)%xmin
        ymin=histos_mh(ind)%ymin
        dx=histos_mh(ind)%dx
        dy=histos_mh(ind)%dy

        do iy=1,histos_mh(ind)%ny
          y=ymin+(iy-1)*dy
          do ix=1,histos_mh(ind)%nx
            x=xmin+(ix-1)*dx
            write(6,'(2I5,"  ",7(1PE13.4))') ix,iy,x,y,
     &        histos_mh(ind)%channels(1:4,ix,iy)
          enddo
        enddo
        write(6,*) ""

      else

        write(6,*) ""
        write(6,*) "  ix   x           n           h           hmean       hrms"
        write(6,*) "------------------------------------------------------------------"

        do ix=1,histos_mh(ind)%nx
          x=xmin+(ix-1)*dx
          write(6,'(I5,"  ",5(1PE12.5))') ix,x,
     &      histos_mh(ind)%channels(1:4,ix,1)
        enddo

      endif !2d

      lastid_mh=id
      lastind_mh=nhist_mh

9999  continue
      return
      end
