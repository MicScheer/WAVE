*CMZ :  4.00/14 19/12/2021  15.38.08  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_nout(ind)

      use mhbook_mod
      implicit none

      integer ind,lunhis,nvar,ivar,ient,lunn
      double precision v(nvarp_mh)


      if (tups_mh(ind)%ktouched.eq.2.or.tups_mh(ind)%nvar.le.0) goto 9999

      if (ind.le.0.or.ind.gt.nntup_mh) then
        print*,"*** Error in mh_nout: Ntuple index out of range",ind
        goto 9999 !return
      endif

      if (kfile_mh.le.0.or.kfile_mh.gt.nfile_mh) then
        print*,"*** Error in mh_nout: File index out of range",kfile_mh
        goto 9999 !return
      endif

      if (luns_mh(kfile_mh).eq.0) then
        print*,"*** Error in mh_nout: File not open for writing",
     &    trim(chfiles_mh(kfile_mh))
        goto 9999 !return
      endif

      lunhis=luns_mh(kfile_mh)

      write(lunhis,*)'! ------------------------------------------------'
      write(lunhis,*)tups_mh(ind)%id,' 3   ! id and kind of histogram of Ntuple'

      write(lunhis,*)len_trim(tups_mh(ind)%title),'   ! length of title'
      write(lunhis,'(a)')trim(tups_mh(ind)%title)

      write(lunhis,*)
     &  len_trim(tups_mh(ind)%chpath),tups_mh(ind)%nvar,lenvarp_mh,
     &  '   ! length of pathname, number of variables, and length of variable names'

      write(lunhis,'(a)') trim(tups_mh(ind)%chpath) // ' ! path'

      nvar=tups_mh(ind)%nvar

      do ivar=1,nvar
        write(lunhis,'(a,2e25.12e3,a)')
     &    tups_mh(ind)%chvar(ivar), tups_mh(ind)%varm(:,ivar),'   ! variable, mininum and maximum'
      enddo

      write(lunhis,*)tups_mh(ind)%nentries,'   ! number of entries'

      lunn=tups_mh(ind)%lun

      if (lunn.ne.0) then
        rewind(lunn)
        do ient=1,tups_mh(ind)%nentries-tups_mh(ind)%neve
          read(lunn)v(1:nvar)
          write(lunhis,*)ient
          write(lunhis,'(5e25.12e3)') v(1:nvar)
        enddo
      endif

      do ient=1,tups_mh(ind)%neve
        write(lunhis,*) tups_mh(ind)%nentries-tups_mh(ind)%neve+ient
        write(lunhis,'(5e25.12e3)')tups_mh(ind)%eve(1:nvar,ient)
      enddo

      tups_mh(ind)%ktouched=2

      lastnid_mh=tups_mh(ind)%id
      lastnind_mh=ind

9999  continue
      return
      end
