      SUBROUTINE READ2 ( NUNIT, IREC, X, NSIZE )
      DIMENSION X(NSIZE)
      READ (NUNIT,REC=IREC) X
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine READ2>>>>>>>>>>>>>>>>>>>>>>>>>>
      END
