*CMZ :  4.00/14 02/01/2022  10.25.32  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_nascii(ind,icode,inoheader,chisascii,code)

      use mhbook_mod
      implicit none

      double precision v(nvarp_mh)

      integer mh_exists,ind,ihkind,icode,inoheader,ic,itit,nvar,lunu,nevent,
     &  ntagmx,i,lunn

      character(16) crun,cid
      character(80) c80
      character(124) file,cline
      character(2048) cwrite
      character chisascii
      character c
      character(*) code
      equivalence(ic,c)

      write(crun,*) abs(icode)
      c80=trim(adjustl(tups_mh(ind)%title))

      do itit=1,len_trim(c80)
        c=c80(itit:itit)
        if (c.ge.'A'.and.c.le.'Z') ic=ic+32
        if (
     &      c.ge.'A'.and.c.le.'Z'
     &      .or.c.ge.'a'.and.c.le.'z'
     &      .or.c.ge.'0'.and.c.le.'9'
     &      .or.c.eq.'('.and.c.eq.')'
     &      ) THEN
          c80(itit:itit)=c
        else
          c80(itit:itit)='_'
        endif
      enddo

      write(cid,*)tups_mh(ind)%id
      cid=trim(adjustl(cid))

      file=trim(c80)//'_'//trim(cid)//'.wvh'

      nevent=tups_mh(ind)%nentries
      nvar=tups_mh(ind)%nvar

      if (nvar.le.5) then
        open(newunit=lunu,file=file,recl=256)
      else if (nvar.le.10) then
        open(newunit=lunu,file=file,recl=512)
      else if (nvar.le.20) then
        open(newunit=lunu,file=file,recl=1024)
      else if (nvar.le.50) then
        open(newunit=lunu,file=file,recl=2048)
      else
        open(newunit=lunu,file=file)
      endif

      ntagmx=0
      do i=1,nvar
        if (len_trim(tups_mh(ind)%chvar(i)).gt.ntagmx) ntagmx=len_trim(tups_mh(ind)%chvar(i))
      enddo

      if (inoheader.eq.0) then
        write(lunu,'(a)')chisascii//' Ntuple'
        cline(1:1)=chisascii
        cline(2:2)=' '
        write(cline(3:10),'(i8)')icode
        cline=cline(1:10)//' | '//code
        write(lunu,'(a)')trim(cline)
        write(cline(3:10),'(i8)')tups_mh(ind)%id
        cline=cline(1:10)//' | '//trim(c80)
        write(lunu,'(a)')trim(cline)
        write(cline,*)chisascii,' ',nvar,ntagmx,nevent," | number and length of parameters, number of entries"
        write(lunu,'(a)')trim(adjustl(cline))
        write(cwrite,*)chisascii,' ',(tups_mh(ind)%chvar(i)(1:ntagmx+1),i=1,nvar)
        write(lunu,'(a)')trim(adjustl(cwrite))
        write(lunu,'(a)')chisascii
      endif   !inoheader

      lunn=tups_mh(ind)%lun

      if (lunn.ne.0) then
        rewind(lunn)
        do i=1,tups_mh(ind)%nentries-tups_mh(ind)%neve
          read(lunn)v(1:nvar)
          write(cwrite,*) v(1:nvar)
          write(lunu,*)cwrite(2:len_trim(cwrite))
        enddo
      endif

      do i=1,tups_mh(ind)%neve
        write(cwrite,*)tups_mh(ind)%eve(1:nvar,i)
        write(lunu,*)cwrite(2:len_trim(cwrite))
      enddo

      close(lunu)

9999  continue
      return
      end
