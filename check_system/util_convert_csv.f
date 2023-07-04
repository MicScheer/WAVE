*CMZ : 00.00/16 20/08/2014  16.00.10  by  Michael Scheer
*-- Author :    Michael Scheer   15/08/2014
      subroutine util_convert_csv(lunin,lunout,filein,fileout,title,istat)

c Converts CSV-files, where each variable is followed by datas until LF.

      implicit none

c      include 'phyconnew.cmn'
*KEEP,PHYCONNEW.
      DOUBLE PRECISION ::
     &  HBAREV1=6.58211889D-16
     &, CLIGHT1=2.99792458D8
     &, EMASSKG1=9.10938188D-31
     &, EMASSE1=0.510998902D6
     &, EMASSG1=0.510998902D-3
     &, ECHARGE1=1.602176462D-19
     &, ERAD1=2.8179380D-15
     &, EPS01=8.854187817D-12
     &, HPLANCK1=6.626176D-34
     &, PI1=3.141592653589793D0
     &, rmu04pi1=1.0D-7
     &, HBAR1=1.0545715955643569D-034
     &, WTOE1=1239.8619236366035d0
     &, CQ1=3.8319386004820543D-013
     &, POL1CON1=  0.92376043070340130d0
     &, POL2CON1=4.3617244578340901D-003
     &, HALFPI1=1.5707963267948966d0
     &, DNULL1=0.0d0
     &, DONE1=1.0d0
     &, RMU01=1.2566370614359173D-006
     &, ALPHA1=7.2973525357634249D-003
     &, GAUSSN1=0.39894228040143270d0
     &, CGAM1=8.8462690127388661D-005
     &, TWOPI1=6.2831853071795862D0
     &  , SQRTTWOPI1=2.5066282746310002D0
     &  ,RADGRA1=57.295779513082323d0
     &  ,grarad1=1.7453292519943295D-002

      character:: clf,ccr,cquote='"'
      integer :: ilf=10, icr=13
      equivalence (ilf,clf)
      equivalence (icr,ccr)
*KEND.

      integer :: lunin,lunout,istat,i,mchar,npos,nvarmaxp,ivar,ival,
     &  nvar=1,
     &  nval=0,
     &  nchar=1

      parameter (nvarmaxp=100)

      double precision buff(nvarmaxp)
      integer luns(nvarmaxp),lenmax

      character c1
      character(*) filein,fileout,title
      character(64):: fscratch(nvarmaxp),cname(nvarmaxp),cform
      character(128):: cline
      character, dimension (:), allocatable :: cbuff

      open(unit=lunin,file=filein(1:len_trim(filein)),status='old',
     &  access="stream")

      mchar=0
      write(cline,*)'fcsv.dum.',nvar
      do i=1,len_trim(cline)
        if(cline(i:i).ne.' ') then
          mchar=mchar+1
          fscratch(1)(mchar:mchar)=cline(i:i)
        endif
      enddo

      open(unit=lunout,file=fscratch(1)(1:mchar),status='unknown',
     &  access="stream")

      nchar=1
      npos=0
 1    read(lunin,pos=nchar,end=9)c1

      npos=npos+1

      if (c1.eq.',') then
        if (nvar.eq.1) nval=nval+1
        write(lunout,pos=npos) clf
      else

        write(lunout,pos=npos) c1

        if (c1.eq.clf) then

          close(lunout)
          nvar=nvar+1
          if (nvar.ge.nvarmaxp) then
            stop '*** Error in util_convert_csv: Too many variables!'
          endif
          write(cline,*)'fcsv.dum.',nvar
          mchar=0
          npos=0
          do i=1,len_trim(cline)
            if(cline(i:i).ne.' ') then
              mchar=mchar+1
              fscratch(nvar)(mchar:mchar)=cline(i:i)
            endif
          enddo

          open(unit=lunout,file=fscratch(nvar)(1:mchar),status='unknown',
     &      access="stream")

        endif
      endif

      nchar=nchar+1

      goto 1
 9    continue

      close(lunin)
      close(lunout)

      nvar=nvar-1

      open(unit=lunout,file=fileout(1:len_trim(fileout)),status='unknown',
     &  recl=512)

      do ivar=1,nvar
        luns(ivar)=10000+ivar
        open(unit=luns(ivar),
     &    file=fscratch(ivar)(1:len_trim(fscratch(ivar))),
     &    status='old')
      enddo

      lenmax=0
      do ivar=1,nvar
        read(luns(ivar),'(a)')cname(ivar)
        lenmax=lenmax+len_trim(cname(ivar))
      enddo

      lenmax=lenmax+nvar-1
      allocate (cbuff(lenmax))

      nchar=0
      do ivar=1,nvar
        do i=1,len_trim(cname(ivar))
          nchar=nchar+1
          cbuff(nchar:nchar)=cname(ivar)(i:i)
        enddo
        nchar=nchar+1
        cbuff(nchar:nchar)=':'
      enddo

      write(cform,*)'(',lenmax+20+len_trim(title)+len_trim(filein),'a)'
      write(lunout,cform) '*r* ncre("',
     &  title(1:len_trim(title)),'","',
     &  filein(1:len_trim(filein)),
     &  '","',
     &  cbuff(1:lenmax),
     &  '");'

      write(cform,*)'(',20+len_trim(title)+len_trim(filein),'a)'
      write(lunout,cform) '*r* nread(',
     &  title(1:len_trim(title)),',"',
     &  fileout(1:len_trim(fileout)),'");'

      write(cform,*)'(',20+len_trim(title),'a)'
      write(lunout,cform) '*r* nprint(',
     &  title(1:len_trim(title)),');'

      do ival=1,nval
        do ivar=1,nvar
          read(luns(ivar),*)buff(ivar)
        enddo
        write(lunout,*)buff(1:nvar)
      enddo

      do ivar=1,nvar+1
        cline='rm '//fscratch(ivar)(1:len_trim(fscratch(ivar)))
        call system(cline(1:len_trim(cline)))
      enddo

      deallocate(cbuff)

      return
      end
