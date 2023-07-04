*CMZ :  4.01/00 13/03/2023  13.32.35  by  Michael Scheer
*CMZ :  3.07/00 07/03/2019  13.44.20  by  Michael Scheer
*-- Author :    Michael Scheer   07/03/2019
      module circpinmod

      double precision, dimension (:), allocatable ::
     &  fphir_th, sumycirc, speccirc

!$OMP THREADPRIVATE(fphir_th, sumycirc, speccirc)

      end module circpinmod
