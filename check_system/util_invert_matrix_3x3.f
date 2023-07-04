*CMZ :          28/11/2021  11.53.56  by  Michael Scheer
*-- Author :    Michael Scheer   10/01/2018
      subroutine util_invert_matrix_3x3(a,ainv,ifail)

c +PATCH,//UTIL/UTIL
c +DECK,util_invert_matrix_3x3.

      ! Calcutate

      implicit none

      double precision a(3,3),ainv(3,3),uni(3,3)
      integer :: ifail

      data uni/1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0/

      call util_solve_matrix_3x3(a,uni,ainv,ifail)

      return
      end
