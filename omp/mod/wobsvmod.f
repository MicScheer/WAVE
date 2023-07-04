*CMZ :  4.01/00 13/03/2023  13.32.35  by  Michael Scheer
*CMZ :  3.07/00 08/03/2019  18.46.38  by  Michael Scheer
*-- Author :    Michael Scheer   05/03/2019
      module wobsvmod

      double precision, dimension (:,:), allocatable :: obsv_th

      double precision, dimension (:), allocatable ::
     &  x_th,f_th,spec_th,
     &  wobsv1_th,wobsv2_th,wobsv3_th,wobsv4_th,wobsv5_th,wobsv6_th,wobsv7_th

!$OMP THREADPRIVATE(wobsv1_th,wobsv2_th,wobsv3_th,wobsv4_th,wobsv5_th,wobsv6_th,wobsv7_th,
!$OMP&       spec_th,x_th,f_th)

      end module wobsvmod
