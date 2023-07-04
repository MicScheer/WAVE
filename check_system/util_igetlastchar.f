*CMZ :          08/03/2018  14.14.52  by  Michael Scheer
*CMZ : 00.00/15 29/08/2012  12.02.35  by  Michael Scheer
*CMZ : 00.00/06 07/03/2007  16.58.44  by  Michael Scheer
*CMZ : 00.00/05 07/03/2007  12.39.29  by  Michael Scheer
*CMZ : 00.00/02 13/09/2002  13.37.13  by  Michael Scheer
*CMZ :  1.19/07 22/08/2002  15.44.21  by  Michael Scheer
*CMZ :  1.00/00 22/11/2001  16.14.10  by  Michael Scheer
*CMZ :  0.01/02 20/11/2001  15.36.03  by  Michael Scheer
*CMZ :  0.01/01 19/11/2001  15.17.00  by  Michael Scheer
*CMZ :  0.01/00 16/11/2001  17.38.53  by  Michael Scheer
*-- Author :    Michael Scheer   09/11/2001

      FUNCTION UTIL_IGETLASTCHAR(N1,N,CLINE,C)

      IMPLICIT NONE

      INTEGER N1,I,N,UTIL_IGETLASTCHAR,IC,ILEN,NN
      CHARACTER(*) CLINE
      CHARACTER C,C1

      EQUIVALENCE (IC,C1)
        ic=0

        ic=len_trim(cline)
      util_igetlastchar=ic
        c=cline(ic:ic)

      RETURN
      END
