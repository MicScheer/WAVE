*CMZ :  3.01/02 23/08/2013  12.58.34  by  Michael Scheer
*-- Author :    Michael Scheer   23/08/2013
      program bmap_to_rec

c +PATCH,//WAVE/MAIN
c +DECK,bmap_to_rec.

c Read bmap_to_rec.bmap and writes bmap_to_rec.par for additional magnets
c for REC input to WAVE.

      implicit none

      double precision x,y,z,bx,by,bz,xlentot,xlen,xhgt,width,xmin,xmax,
     &  by1,by2,bx1,bx2

      integer last,nread,istat,ieof

      character(256)filein,fileout
      character(2048)cline
      character c1

      data filein/'bmap_to_rec.bmap'/
      data fileout/'bmap_to_rec.par'/

      data xhgt/0.0001/ !m
      data width/0.01/ !m

      open(unit=20,file=filein,status='old')

      call util_string_lastcharacter(fileout,last,c1,istat)
      cline='mv '//fileout(1:last)//' '//fileout(1:last)//'.bck'//' 2>/dev/null'
      call util_string_lastcharacter(cline,last,c1,istat)
      call system(cline(1:last))
      open(unit=21,file=fileout,status='unknown', recl=256)

      xmin=1.0e30
      xmax=-1.0e30
      nread=0
1     call util_skip_comment_end(20,ieof)
      if (ieof.ne.0) goto 9
      read(20,*)x,y,z,bx,by,bz
      nread=nread+1
      if (x.lt.xmin) xmin=x
      if (x.gt.xmax) xmax=x
      goto 1
9     rewind(20)

      xlentot=xmax-xmin
      xlen=xlentot/(nread-1)

      write(21,*)2*nread-2,'   !Berechnet mit bmap_to_rec.exe'
      nread=0
11    call util_skip_comment_end(20,ieof)
      if (ieof.ne.0) goto 99
      bx1=bx2
      by1=by2
      read(20,*)x,y,z,bx2,by2,bz
      nread=nread+1
      if (nread.eq.1) goto 11
      x=x-xlen/2.0d0
      by=(by2+by1)/2.0d0
      bx=(bx2+bx1)/2.0d0
      write(21,*)xlen*1000.0d0,', ',xhgt*1000.0d0,' ,', width*1000.0d0,
     &  '   ! xlen, xhgt, width ============='
      if (bx.ge.0.0d0) then
        write(21,*)' 90.0, 0.0', 0.60d0*abs(bx)/sqrt(bx**2+by**2),'   ! theta, phi, bc'
      else
        write(21,*)'270.0, 0.0', 0.60d0*abs(bx)/sqrt(bx**2+by**2),'   ! theta, phi, bc'
      endif
      write(21,*)x*1000.0d0,', ',y*1000.0d0,', ',z*1000.0d0,
     &  '   ! dx0, dy0, dz0'
      write(21,*)xlen*1000.0d0,', ',xhgt*1000.0d0,' ,', width*1000.0d0,
     &  '   ! xlen, xhgt, width ============='
      if (by.ge.0.0d0) then
        write(21,*)'  0.0, 0.0', 0.60d0*abs(by)/sqrt(bx**2+by**2),'   ! theta, phi, bc'
      else
        write(21,*)'180.0, 0.0', 0.60d0*abs(by)/sqrt(bx**2+by**2),'   ! theta, phi, bc'
      endif
      write(21,*)x*1000.0d0,', ',y*1000.0d0,', ',z*1000.0d0,
     &  '   ! dx0, dy0, dz0'
      goto 11

99    close(21)
      close(20)

      stop
      end
