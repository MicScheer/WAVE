*CMZ :          22/01/2018  16.20.00  by  Michael Scheer
*-- Author :    Michael Scheer   22/01/2018
      SUBROUTINE util_time_comment(LUN,COMMENT)

C To determine date and time and write it to logical unit LUN

      IMPLICIT NONE
      INTEGER LUN
      CHARACTER(64) COMMENT

      call util_zeit_kommentar(lun,comment)

      RETURN
      END
