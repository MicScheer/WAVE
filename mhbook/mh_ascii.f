*CMZ :  4.00/14 31/12/2021  17.39.44  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_ascii(jd,ihisascii,icode,chisascii,code)

      use mhbook_mod
      implicit none

      integer mh_exists,jd,id,ind,ihkind,ihisascii,icode,
     &  jhisascii,i1dim,i2dim,intup,inoheader

      character(*) code
      character chisascii


      jhisascii=abs(ihisascii)

      i1dim=(jhisascii-jhisascii/10*10)/1
      i2dim=(jhisascii-jhisascii/100*100)/10
      intup=(jhisascii-jhisascii/1000*1000)/100

      if (ihisascii.lt.0) then
        inoheader=1
      else
        inoheader=0
      endif

      if (jd.eq.0) then
        if (i1dim.ne.0.or.i2dim.ne.0) then
          do ind=1,nhist_mh
            if (histos_mh(ind)%ny.le.0.and.i1dim.ne.0) then
              call mh_h1ascii(ind,icode,inoheader,chisascii,code)
            else if (histos_mh(ind)%ny.gt.0.and.i2dim.ne.0) then
              call mh_h2ascii(ind,icode,inoheader,chisascii,code)
            endif
          enddo
        endif !(i1dim.ne.0) then
        if (intup.ne.0) then
          do id=1,nntup_mh
            call mh_nascii(id,icode,inoheader,chisascii,code)
          enddo
        endif !(i1dim.ne.0) then
      else
        ind=mh_exists(jd,ihkind)
        if (ind.le.0) goto 9999
        if (ihkind.eq.1.and.i1dim.ne.0) then
          call mh_h1ascii(ind,icode,inoheader,chisascii,code)
        else if (ihkind.eq.2.and.i2dim.ne.0) then
          call mh_h2ascii(ind,icode,inoheader,chisascii,code)
        else if (ihkind.gt.2.and.intup.ne.0) then
          call mh_nascii(ind,icode,inoheader,chisascii,code)
        endif
      endif !jd.eq.0

9999  continue
      return
      end
