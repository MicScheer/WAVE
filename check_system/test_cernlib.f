*CMZ :  3.01/06 20/06/2014  10.33.52  by  Michael Scheer
*-- Author :    Michael Scheer   20/06/2014
      program test_cernlib

c +PATCH,//WAVE/MAIN
c +DECK,test_cernlib.

      implicit none

      real h(1000000),tup(2)
      integer i,istat,icycle

      character(4) chtags(2)

      common/pawc/h

      call hlimit(1000000)

      call hropen(20,'test','test.his','NQ',1024,istat)
      chtags(1)='x1'
      chtags(2)='x2'
      tup(1)=1.
      tup(2)=2.
      call hbookn(1,'Ntest',2,'//test',1024,chtags)
      do i=1,10000
        call hfn(1,tup)
      enddo
      call hrout(0,icycle,' ')
      call hrend('test')
      close(20)

      stop
      end
