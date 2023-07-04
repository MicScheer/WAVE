*CMZ :  4.00/08 06/08/2020  16.13.52  by  Michael Scheer
*CMZ :  3.07/01 21/03/2019  10.29.57  by  Michael Scheer
*CMZ :  3.05/02 15/05/2018  17.04.11  by  Michael Scheer
*CMZ :  3.05/01 04/05/2018  16.03.55  by  Michael Scheer
*-- Author :    Michael Scheer   04/05/2018
      subroutine omp_ini(lungfo,mthreads,komp)

      use omp_lib
      use ompmod
      use bpolyederf90m

      implicit none

      integer lungfo,mthreads,komp

      lungfo_bpoly=lungfo

      iomp=0
      if (komp.eq.0.and.mthreads.ne.0) then
        write(6,*)' '
        write(6,*)' *** Warning in omp_ini: Mthreads is not zero, but this version of WAVE does not support OpenMP ***'
        write(6,*)' '
        write(lungfo,*)' '
        write(lungfo,*)' *** Warning in omp_ini: Mthreads is not zero, but this version of WAVE does not support OpenMP ***'
        write(lungfo,*)' '
        mthreads=0
      endif
      if (mthreads.ne.0) then
        mmaxthreads=OMP_GET_MAX_THREADS()
        iomp=1
        write(lungfo,*)' '
        write(lungfo,*)' '
        if (mthreads.lt.0.or.mthreads.gt.mmaxthreads) mthreads=mmaxthreads
        write(lungfo,*)'       Max number of threads available, number of used ones:'
        write(lungfo,*)'       ',mmaxthreads, mthreads
      endif
      nthreads=mthreads

      return
      end
