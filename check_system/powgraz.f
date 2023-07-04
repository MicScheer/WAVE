*CMZ :  4.00/11 15/06/2021  11.16.36  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine powgraz

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEND.
      use sourcef90
      use observf90

      implicit none
c      include 'uservar.cmn'

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,observ.
      include 'observ.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEND.

      double precision sou(3),x,y,z,dx,dy,dz,dxs,dys,dzs

      integer :: iobsv,isour,lun,lunob,k

      do isour=1,nsource
        do iobsv=1,nobsv
          k=isour+(iobsv-1)*nsource
          dxs=obsv(1,iobsv)-schwingercen(1,iobsv,isour)
          dys=obsv(2,iobsv)-schwingercen(2,iobsv,isour)
          dzs=obsv(3,iobsv)-schwingercen(3,iobsv,isour)
          specpowtgraz(iobsv)=specpowtgraz(iobsv)+specpow(k)*
     &      sqrt((dys**2+dzs**2)/(dxs**2+dys**2+dzs**2))
        enddo !nobsv
      enddo ! nsource

      return
      end
