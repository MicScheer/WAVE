*CMZ :  4.00/14 09/12/2021  16.07.44  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_list1

      use mhbook_mod
      implicit none

      integer i


      print*,""
      print*,"Index of 1d-Histograms"
      print*,""
      print*,"          Index       Id    Title      nEntries"

      do i=1,nhist_mh
        if (histos_mh(i)%ny.gt.0) cycle
        print*,i,histos_mh(i)%id,"    ",trim(histos_mh(i)%title),
     &    histos_mh(i)%nentries
      enddo
      print*,""

      lastid_mh=0
      lastind_mh=0

!9999  continue
      return
      end
