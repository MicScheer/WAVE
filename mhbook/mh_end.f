*CMZ :  4.00/14 21/12/2021  11.30.58  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_end

      use mhbook_mod
      implicit none


!9999  continue

      if (lun_index_mh.ne.0) then
        write(lun_index_mh,*)" "
        write(lun_index_mh,*)"Maximum index for storage:",
     &    nhistmax_mh+nntupmax_mh
        flush(lun_index_mh)
        close(lun_index_mh)
      endif

      lastid_mh=0
      lastind_mh=0

      return
      end
