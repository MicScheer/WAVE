*CMZ :  2.61/02 09/03/2007  18.29.49  by  Michael Scheer
*-- Author :    Michael Scheer   09/03/2007
      program addtherot_main
c +PATCH,//WAVE/MAIN
c +DECK,addtherot_main.

      implicit none

      double precision radius1,width1,xcen1,ycen1,zcen1,therot1,bcr1
      double precision radius2,width2,xcen2,ycen2,zcen2,therot2,bcr2

      integer lun1,lun2,lunsum,iread1,iread2,i

      character(256) file1,file2,filesum

      data lun1,lun2,lunsum/21,22,23/

      data file1/'therot.add1'/
      data file2/'therot.add2'/
      data filesum/'therot.sum'/

      open(unit=lun1,file=file1,status='old')
      open(unit=lun2,file=file2,status='old')

      call system('mv -f therot.sum.bck therot.sum.bck.tmp')
      call system('mv -f therot.sum therot.sum.bck')

      open(unit=lunsum,file=filesum,status='new')

      read(lun1,*)iread1
      read(lun2,*)iread2

      if (iread2.lt.iread1) then
        i=lun2
        lun2=lun1
        lun1=i
        i=iread2
        iread2=iread1
        iread1=iread2
      endif

      rewind(lun1)
      rewind(lun2)

      read(lun1,*)iread1
      read(lun2,*)iread2

      do i=1,iread2

        read(lun2,*)radius2,width2,xcen2,ycen2,zcen2,therot2,bcr2

        if (i.le.iread1) then

          read(lun1,*)radius1,width1,xcen1,ycen1,zcen1,therot1,bcr1

          if (
     &        abs(radius2-radius1)+
     &        abs(width2-width1)+
     &        abs(xcen2-xcen1)+
     &        abs(ycen2-ycen1)+
     &        abs(zcen2-zcen1)+
     &        abs(therot2-therot1)+
     &        abs(bcr2-bcr1).gt.1.0d-6
     &        ) then

            close(lunsum)
            close(lun1)
            close(lun2)

            call system('mv -f therot.sum.tmp therot.sum')
            call system('mv -f therot.sum.bck.tmp therot.sum.bck')

            print*,'*** Error in addtherot_main: Magnets do not match'
            print*,'    or rotation angles not zero'
            print*,'*** Giving up leaving files untouched.'

            stop

          endif

          bcr2=bcr2+bcr1

        endif

        write(lunsum,*)radius2,width2,xcen2,ycen2,zcen2,therot2,bcr2

      enddo

      close(lunsum)
      close(lun1)
      close(lun2)


      stop
      end
