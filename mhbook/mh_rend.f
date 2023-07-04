*CMZ :  4.00/14 31/12/2021  15.13.30  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_rend(chfin)

      use mhbook_mod
      implicit none

      integer i,ind
      character(*) chfin
      character(2048) chfile


      chfile=trim(adjustl(chfin)) // ".mhb"

      do i=1,nfile_mh
        if (chfile.eq.chfiles_mh(i)) then
          if (luns_mh(i).ne.0) then
            do ind=1,nhist_mh
              if (histos_mh(ind)%ny.le.0) then
                call mh_h1out(ind)
              else
                call mh_h2out(ind)
              endif
            enddo
            do ind=1,nntup_mh
              call mh_nout(ind)
            enddo
            flush(luns_mh(i))
            close(luns_mh(i))
            goto 9999
          endif
        endif
      enddo

      print*,"*** Error in mh_rend: File " // trim(chfile) // " was not open ***"
9999  continue
      return
      end
