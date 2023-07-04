*CMZ : 00.00/01 20/07/95  11.21.28  by  Michael Scheer
*-- Author :    Michael Scheer   23/06/95

      SUBROUTINE UTIL_HIGZ_STYLE

C--- ROUTINE TO STYLE HIGZ PLOTS

      CALL HPLSET('YGTI',1.0)
      CALL HPLSET('GSIZ',0.35)

      CALL HPLSET('BWID',2.)
      CALL HPLSET('HWID',2.)
      CALL HPLSET('PWID',2.)

      CALL HPLSET('VFON',-10.)
      CALL HPLSET('LFON',-10.)
      CALL HPLSET('GFON',-60.)
      CALL HPLSET('CFON',-40.)
      CALL HPLSET('TFON',-40.)


      CALL IGSET('TXFP',-40.)

      RETURN
      END
