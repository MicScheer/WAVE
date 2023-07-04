*CMZ :  3.05/23 21/11/2018  13.03.11  by  Michael Scheer
*-- Author :    Michael Scheer   21/11/2018
      function util_atan2(y,x)

      implicit none

      double precision x,y,util_atan2

      util_atan2=atan2(y,x)
      if (util_atan2.lt.0.0d0) util_atan2=util_atan2+6.2831853071795862d0

      return
      end
