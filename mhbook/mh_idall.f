*CMZ :  4.00/14 31/12/2021  09.58.08  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_idhall(idvect,nid)

      use mhbook_mod
      implicit none

      integer idvect(*),nid,i


      nid=nhist_mh+nntup_mh

      do i=1,nhist_mh
        idvect(i)=histos_mh(i)%id
      enddo

      do i=1,nntup_mh
        idvect(i+nhist_mh)=tups_mh(i)%id
      enddo

9999  continue
      return
      end
