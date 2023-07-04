*CMZ :  4.00/14 11/12/2021  15.51.50  by  Michael Scheer
*-- Author :    Michael Scheer   11/12/2021
      function iutil_memory_size(istart8)

      implicit none

      integer*8 :: iutil_memory_size,istart8,isize=1e6
      integer :: istat=0,ide=0

      double precision, dimension(:), allocatable :: a

      isize=max(isize,istart8)
      deallocate(a,stat=ide)

      do while (istat.eq.0)
        allocate(a(isize),stat=istat)
        isize=2*isize
        deallocate(a,stat=ide)
        !print*,isize, istat,ide
      enddo

      iutil_memory_size=isize/2

      return
      end
