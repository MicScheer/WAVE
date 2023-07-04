*CMZ :  4.00/16 09/09/2022  17.29.00  by  Michael Scheer
*CMZ :  2.48/04 12/03/2004  15.40.31  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.14/02 20/04/2000  14.28.52  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.36.30  by  Michael Scheer
*-- Author :    Michael Scheer   08/02/2000
      DOUBLE PRECISION FUNCTION DCOSD(X)

      DOUBLE PRECISION X,GR
      PARAMETER (GR=0.0174532925199D0)

      DCOSD=COS(X*GR)

      RETURN
      END
