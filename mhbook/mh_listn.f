*CMZ :  4.00/14 19/12/2021  10.34.22  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_listn

      use mhbook_mod
      implicit none

      integer i


      print*,""
      print*,"Index of Ntuples"
      print*,""
      print*,"          Index       Id    Title      nEntries"

      do i=1,nntup_mh
        print*,i,tups_mh(i)%id,"    ",trim(tups_mh(i)%title),
     &    tups_mh(i)%nentries
      enddo
      print*,""

      lastid_mh=0
      lastind_mh=0

!9999  continue
      return
      end
