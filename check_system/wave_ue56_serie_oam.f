*CMZ :  2.67/02 20/04/2012  09.13.59  by  Michael Scheer
*-- Author :    Michael Scheer   20/04/2012
      program wave_ue56_serie

c +PATCH,//WAVE/MAIN
c +DECK,wave_ue56_serie_oam.
c *CMZ :          20/04/2012  09.13.24  by  Michael Scheer
c *-- Author :    Michael Scheer   20/04/2012
      implicit none

      real u1harm1, u1b01, bm,buff(10000),bmo
      integer iold,im,nm

      open(unit=20,file='wave_ue56_serie_oam.out',status='old',readonly)

      nm=0
      bmo=-9999.
1     read(20,*,end=9) u1harm1, u1b01, bm
      if (bm.ne.bmo) then
        iold=0
        do im=1,nm
          if (bm.eq.buff(im)) then
            iold=1
            goto 2
          endif
        enddo
2       if (iold.eq.0) then
          nm=nm+1
          buff(nm)=bm
        endif
      endif
      bmo=bm
      goto 1

9     close(20)

      call util_sort_single(nm,buff)

      open(unit=21,file='wave_ue56_serie_oam.bm')
      write(21,*)nm
      do im=1,nm
        write(21,*)buff(im)
      enddo
      close(21)

      stop
      end
