*CMZ :  4.00/11 17/05/2021  14.09.13  by  Michael Scheer
*-- Author :    Michael Scheer   17/05/2021
      subroutine fbtab_ini(x,irbtab,irbtabzy,irbtabxyz)

      use fbtabzymod

      implicit none
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEND.

      DOUBLE PRECISION x,y,z,bx,by,bz,ax,ay,az,
     &  XA(NBTABP),BYA(NBTABP),Y2A(NBTABP),
     &  XAz(NBTABP),BzA(NBTABP),z2A(NBTABP)

      integer irbtab,irbtabzy,irbtabxyz

      integer :: ical=0

      COMMON/BTABC/XA,BYA,Y2A
      COMMON/BTABCz/XAz,BzA,z2A

      save ical

      if (ical.ne.0) return

      if (nfbtabc.eq.0) then
        if(irbtabxyz.ne.0) then
          call btabxyz(x,0.0d0,0.0d0,bx,by,bz,ax,ay,az)
        else if(irbtabzy.ne.0) then
          call btabzy(x,0.0d0,0.0d0,bx,by,bz,ax,ay,az)
        else if(irbtab.ne.0) then
          call btab(x,0.0d0,0.0d0,bx,by,bz,ax,ay,az)
        else
          irbtab=1
          call btab(x,0.0d0,0.0d0,bx,by,bz,ax,ay,az)
          irbtab=0
        endif
      endif

      nxbyfbt=nfbtabc
      allocate(xbyfbt(nxbyfbt),byfbt(nxbyfbt))
      xbyfbt(1:nxbyfbt)=xa(1:nxbyfbt)
      byfbt(1:nxbyfbt)=bya(1:nxbyfbt)

      if (nfbtabcz.eq.0) then
        call btabz(x,0.0d0,0.0d0,bx,by,bz,ax,ay,az,0.0d0,0.0d0)
      endif

      nxbzfbt=nfbtabc
      allocate(xbzfbt(nxbzfbt),bzfbt(nxbzfbt))
      xbzfbt(1:nxbzfbt)=xaz(1:nxbzfbt)
      bzfbt(1:nxbzfbt)=bza(1:nxbzfbt)

      ical=1

      return
      end
