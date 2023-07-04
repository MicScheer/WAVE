*CMZ :  4.00/09 13/08/2020  13.17.34  by  Michael Scheer
*CMZ :  4.00/08 07/08/2020  11.10.43  by  Michael Scheer
*CMZ :  4.00/07 06/08/2020  10.51.12  by  Michael Scheer
*CMZ :  1.25/00 16/03/2018  14.11.34  by  Michael Scheer
*CMZ :  1.23/03 19/09/2017  19.25.01  by  Michael Scheer
*CMZ :  1.23/02 30/08/2017  13.27.12  by  Michael Scheer
*CMZ :  1.22/02 31/07/2017  10.32.51  by  Michael Scheer
*CMZ :  1.22/01 20/07/2017  14.46.06  by  Michael Scheer
*CMZ :  1.22/00 05/07/2017  09.55.55  by  Michael Scheer
*CMZ :  1.20/03 29/06/2017  09.17.17  by  Michael Scheer
*CMZ :  1.20/01 22/06/2017  13.26.26  by  Michael Scheer
*CMZ :  1.20/00 22/06/2017  11.26.04  by  Michael Scheer
*CMZ :  1.15/11 24/04/2017  16.58.30  by  Michael Scheer
*CMZ :  1.15/10 12/04/2017  14.53.10  by  Michael Scheer
*CMZ :  1.15/04 03/04/2017  12.30.25  by  Michael Scheer
*CMZ :  1.15/03 03/04/2017  10.59.22  by  Michael Scheer
*CMZ :  1.15/02 02/04/2017  07.35.42  by  Michael Scheer
*CMZ :  1.15/01 28/03/2017  13.53.21  by  Michael Scheer
*CMZ :  1.13/01 08/03/2017  16.31.38  by  Michael Scheer
*CMZ :  1.11/03 16/01/2017  12.22.22  by  Michael Scheer
*CMZ :  1.10/02 24/11/2016  09.47.59  by  Michael Scheer
*CMZ :  1.10/01 18/11/2016  15.02.58  by  Michael Scheer
*CMZ :  1.07/00 23/09/2016  09.19.06  by  Michael Scheer
*CMZ :  1.04/01 14/09/2016  15.10.51  by  Michael Scheer
*CMZ :  1.00/00 19/08/2016  18.27.23  by  Michael Scheer
*CMZ :  0.00/13 28/07/2016  16.09.29  by  Michael Scheer
*CMZ :  0.00/09 06/07/2016  08.42.18  by  Michael Scheer
*CMZ :  0.00/06 16/06/2016  14.14.37  by  Michael Scheer
*CMZ :  0.00/04 13/05/2016  13.18.24  by  Michael Scheer
*CMZ :  0.00/02 29/04/2016  09.17.13  by  Michael Scheer
*CMZ :  0.00/01 25/04/2016  16.03.15  by  Michael Scheer
*CMZ :  0.00/00 20/04/2016  12.41.34  by  Michael Scheer
*CMZ :  1.17/14 13/04/2016  09.46.51  by  Michael Scheer
*CMZ :  1.17/11 05/04/2016  13.27.16  by  Michael Scheer
*CMZ :  1.17/08 04/04/2016  08.57.43  by  Michael Scheer
*CMZ :  1.17/07 04/04/2016  08.31.31  by  Michael Scheer
*CMZ :  1.17/06 01/04/2016  13.53.25  by  Michael Scheer
*CMZ :  1.17/05 27/03/2016  10.43.50  by  Michael Scheer
*CMZ :  1.17/03 21/03/2016  18.38.48  by  Michael Scheer
*-- Author :    Michael Scheer   02/12/2003
      subroutine undumag_bpolyeder_omp(xin,yin,zin,bxout,byout,bzout,ifail)

      use omp_lib

      use bpolyederf90m
      use undumagf90m

      implicit none
*KEEP,seqdebug.
      include 'seqdebug.cmn'
*KEND.

      double precision xin,yin,zin,bxout,byout,bzout,bo(3,nthreadp),
     &  bxm,bym,bzm,bxp,byp,bzp

      integer nmaxth,ith,ifail,ical,kfail(nthreadp),ic,imag,ifailp,ifailm,
     &  ifailin,kinsidelocal(nthreadp)

      save

      data ical/0/

      if (magmag.lt.0) then
        return
      endif !magmag.le.0

      if (ical.eq.0) then

        nmaxth=1
        ith=1
        nmaxth=nthreads
        nmaxth=OMP_GET_MAX_THREADS()
        if (nthreads.gt.0) nmaxth=min(nmaxth,nthreads,nthreadp)
        bo=0.0d0

      endif

      kfail=0
      ifailin=ifail
      kinsidelocal=0

!$OMP PARALLEL NUM_THREADS(nmaxth) DEFAULT(PRIVATE)
!$OMP& SHARED(bo,kfail,kinsidelocal)
!$OMP& FIRSTPRIVATE(nmaxth,nmag,xin,yin,zin,magmag)

      ith=OMP_GET_THREAD_NUM()+1
      bo(1:3,ith)=0.0d0

!$OMP DO

      do imag=1,nmag

        call undumag_bpolyeder_range(imag,imag,xin,yin,zin,bxout,
     &    byout,bzout,ifail)

        if (ifail.ne.0) kfail(ith)=ifail
        if (kinside.ne.0) kinsidelocal(ith)=kinside

        bo(1,ith)=bo(1,ith)+bxout
        bo(2,ith)=bo(2,ith)+byout
        bo(3,ith)=bo(3,ith)+bzout

      enddo

!$OMP END DO
!$OMP END PARALLEL

      bxout=0.0d0
      byout=0.0d0
      bzout=0.0d0

      ifail=0
      kinside=0

      do ic=1,nmaxth
        ifail=ifail+kfail(ic)
        if (kinsidelocal(ic).ne.0) kinside=kinsidelocal(ic)
        bxout=bxout+bo(1,ic)
        byout=byout+bo(2,ic)
        bzout=bzout+bo(3,ic)
      enddo

      if (ifail.ne.0) goto 7799

      goto 7979

7799  continue

      if (corrtiny.eq.0.0) then
        ifail=-4
        goto 7979
      endif

      kwarncom=1
      ifail=-3

C Not working for unknown reasons for OMP ?? still true 20.4.2017??

c      print*,"*** calling undumag_bpolyeder_corr for ",xin,yin,zin
      call undumag_bpolyeder_corr(xin-corrtiny,yin-corrtiny,zin-corrtiny,
     &  bxm,bym,bzm,ifailm)

      call undumag_bpolyeder_corr(xin+corrtiny,yin+corrtiny,zin+corrtiny,
     &  bxp,byp,bzp,ifailp)

      if (ifailm.eq.0.and.ifailp.eq.0) then
        bxout=(bxm+bxp)/2.0d0
        byout=(bym+byp)/2.0d0
        bzout=(bzm+bzp)/2.0d0
        ifail=-1
        kwarncom=1
      else if (ifailm.eq.0) then
        bxout=bxm
        byout=bym
        bzout=bzm
        kwarncom=2
        ifail=2
      else if (ifailp.eq.0) then
        bxout=bxp
        byout=byp
        bzout=bzp
        kwarncom=2
        ifail=2
      else
        if (ifailin.ge.0) then
          print*, "*** Warning in undumag_bpolyeder_omp: Could not recover for x,y,z:",
     &      sngl(xin*1000.),sngl(yin*1000.0),sngl(zin*1000.)
          print*,"Differences in Bx,By,Bz, abs. and. rel.:",
     &      sngl(abs(bxp-bxm)),sngl(abs(byp-bym)),sngl(abs(bzp-bzm))
        endif

        bxout=(bxm+bxp)/2.0d0
        byout=(bym+byp)/2.0d0
        bzout=(bzm+bzp)/2.0d0

        if (ifailin.ge.0) then
          print '(6e15.4)',
     &      abs((bxp-bxm)/(bxout+1.0d-15)),abs((byp-bym)/(byout+1.0d-15)),
     &      abs((bzp-bzm)/(bzout+1.0d-15))
        endif

        kwarncom=3
        ifail=3
      endif

7979  continue

      ical=1

      return
      end
