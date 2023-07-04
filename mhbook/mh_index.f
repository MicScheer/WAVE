*CMZ :  4.00/14 19/12/2021  10.40.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_index

      use mhbook_mod
      implicit none

      integer i,luno


      open(newunit=luno,file='wave_index.his')

      write(luno,*)""
      write(luno,*)"Index of 1d-Histograms"
      write(luno,*)""
      write(luno,*)"          Index       Id    Title      nEntries"

      do i=1,nhist_mh
        if (histos_mh(i)%ny.gt.0) cycle
        write(luno,*)i,histos_mh(i)%id,"    ",trim(histos_mh(i)%title),
     &    histos_mh(i)%nentries
      enddo

      write(luno,*)""
      write(luno,*)"Index of 2d-Histograms"
      write(luno,*)""
      write(luno,*)"          Index       Id    Title      nEntries"

      do i=1,nhist_mh
        if (histos_mh(i)%ny.le.0) cycle
        write(luno,*)i,histos_mh(i)%id,"    ",trim(histos_mh(i)%title),
     &    histos_mh(i)%nentries
      enddo
      write(luno,*)""

      write(luno,*)""
      write(luno,*)"Index of Ntuples"
      write(luno,*)""
      write(luno,*)"          Index       Id    Title      nEntries"

      do i=1,nntup_mh
        write(luno,*)i,tups_mh(i)%id,"    ",trim(tups_mh(i)%title),
     &    tups_mh(i)%nentries
      enddo
      write(luno,*)""

      close(luno)

      lastid_mh=0
      lastind_mh=0

      return
      end
