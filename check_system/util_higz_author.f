*CMZ : 00.00/01 23/06/95  20.48.36  by  Michael Scheer
*-- Author :    Michael Scheer   23/06/95

      SUBROUTINE UTIL_HIGZ_AUTHOR

      CALL HPLZON(2,2,2,'S')
      CALL HPLFRA(0.,20.,0.,20.,'AB')
      CALL IGSET('CHHE',0.75)
      CALL IGSET('TXFP',-40.)
      CALL ITX(18.7,21.75,'M. Scheer')

      RETURN
      END
