*CMZ :  3.05/11 11/01/2018  11.54.22  by  Michael Scheer
*-- Author :    Michael Scheer   10/01/2018
      subroutine util_invert_matrix_2x2(a,ainv,ifail)

c +PATCH,//UTIL/UTIL
c +DECK,util_invert_matrix_2x2.

      ! Calcutate

      implicit none

      double precision a(2,2),ainv(2,2),uni(2,2)
      integer :: ifail

      data uni/1.0d0,0.0d0,0.0d0,1.0d0/

      call util_solve_matrix_2x2(a,uni,ainv,ifail)

      return
      end
