*CMZ :  4.00/14 02/01/2022  10.13.40  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_h1ascii(ind,icode,inoheader,chisascii,code)

      use mhbook_mod
      implicit none

      integer mh_exists,id,ind,ihkind,inoheader,icode,ic,lunu,ix,itit
      double precision dx,x

      character(*) code
      character(16) crun,cid
      character(80)c80
      character(128)cline,file
      character c,chisascii
      equivalence(c,ic)


      write(crun,*) abs(icode)
      id=histos_mh(ind)%id
      write(cid,*) id

      c80=trim(histos_mh(ind)%title(1:80))

      do itit=1,len_trim(c80)
        c=c80(itit:itit)

        IF (C.GE.'A'.AND.C.LE.'Z') ic=ic+32

        IF (
     &      C.GE.'A'.AND.C.LE.'Z'
     &      .OR.C.GE.'a'.AND.C.LE.'z'
     &      .OR.C.GE.'0'.AND.C.LE.'9'
     &      .OR.C.EQ.'('.or.C.EQ.')'
     &      ) THEN
          c80(itit:itit)=c
        else
          c80(itit:itit)='_'
        endif
      ENDDO

      file=trim(adjustl(c80))//'_'//trim(adjustl(cid))//'.wvh'

      open(newunit=lunu,file=file)

      if (inoheader.eq.0) then
        write(lunu,'(a)')chisascii//' 1d histogram'
        cline(1:1)=chisascii
        cline(2:2)=' '
        write(cline(3:10),'(i8)')icode
        cline=cline(1:10)//' | '//trim(adjustl(code))
        write(lunu,'(a)')trim(cline)
        write(cline(3:10),'(i8)')id
        cline=cline(1:10)//' | '//trim(c80)
        write(lunu,'(a)')trim(cline)
      endif

      dx=histos_mh(ind)%dx
      x=histos_mh(ind)%xmin-0.5*dx

      do ix=1,histos_mh(ind)%nx
        x=x+dx
        write(lunu,*) x,histos_mh(ind)%channels(2,ix,1)
      enddo   !ix

      close(lunu)

9999  continue
      return
      end
