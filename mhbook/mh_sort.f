*CMZ :  4.00/14 14/12/2021  21.18.27  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_sort

      use mhbook_mod
      implicit none

      integer i


      do i=1,nntup_mh
        idnbuff_mh(i)=i
        idntup_mh(i)=tups_mh(i)%id
      enddo
      call util_sort_func_integer(nntup_mh,idntup_mh,idnbuff_mh)
      do i=1,nntup_mh
        indntup_mh(i)=idnbuff_mh(nntup_mh-i+1)
      enddo

      do i=1,nhist_mh
        idhbuff_mh(i)=i
        idhis_mh(i)=histos_mh(i)%id
      enddo
      call util_sort_func_integer(nhist_mh,idhis_mh,idhbuff_mh)
      do i=1,nhist_mh
        indhis_mh(i)=idhbuff_mh(nhist_mh-i+1)
      enddo

      return
      end
