*CMZ :  4.00/08 06/08/2020  15.47.58  by  Michael Scheer
*-- Author :    Michael Scheer   06/08/2020
      subroutine undumag_bpolyeder_range_para(ith,nmaxth,nmag,xin,yin,zin,
     &  bxout,byout,bzout,ifail)

      implicit none

      double precision xin,yin,zin,bxout,byout,bzout

      integer, parameter :: nthreadp=1000 !see also module undumagf90m

      integer ith,nmaxth,nmag,ifail,ifirst(nthreadp),
     &  nmagrange,magi(nthreadp),mage(nthreadp)

      integer :: ical=0

      save

      nmagrange=nmag/nmaxth
      ifail=0

      if (ical.eq.0) then
        ifirst=0
        ical=1
      endif

      if (ifirst(ith).ne.1) then
        if (ith.lt.nmaxth) then
          magi(ith)=nmagrange*(ith-1)+1
          mage(ith)=magi(ith)+nmagrange-1
        else
          magi(ith)=nmagrange*(ith-1)+1
          mage(ith)=nmag
        endif
        print*,ith,magi(ith),mage(ith),nmag,ifirst(ith)
        call undumag_bpolyeder_range(magi(ith),mage(ith),xin,yin,zin,
     &    bxout,byout,bzout,ifail)
        ifirst(ith)=1
      endif

      return
      end
