*CMZ :          02/04/2024  12.36.09  by  Michael Scheer
*-- Author :    Michael Scheer   02/04/2024
*KEEP,WIGNERMOD.
      module wignermod

      double precision, dimension(:,:,:,:,:,:), allocatable :: wigr,wigi
      double precision, dimension(:), allocatable :: wigthez,wigthey,wigz,wigy
      double precision :: dzprop,dyprop,dzwig,dywig,pinwwig,pinhwig

      integer nywig,nzwig,nzthewig,nythewig,nzfringe,nyfringe

      end module wignermod
