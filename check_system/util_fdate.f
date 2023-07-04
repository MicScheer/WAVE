*CMZ :          05/09/2020  08.46.01  by  Michael Scheer
*CMZ : 00.00/02 26/01/2004  09.37.26  by  Michael Scheer
*-- Author :    Michael Scheer   27/11/98

      SUBROUTINE util_FDATE(TIMEDATE)

      CHARACTER(24) TIMEDATE
      CHARACTER(9) CDATE
      CHARACTER(8) CTIME

      CALL TIME (CTIME)
      CALL DATE(CDATE)

      TIMEDATE='xxx '//CDATE(4:6)//' '//CTIME//' 19'//CDATE(8:9)

      RETURN
      END
