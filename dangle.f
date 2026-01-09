C -------------------------------------------------------------
        FUNCTION DANGLE( Z , Y )
C   
C       CALCULATE ATAN( Z/Y ) IN RANGE OF 0 TO 2*PI
C       HANSEONG LEE
C-------------------------------------------------------------
        PI = ACOS(-1.0)
        IF( Z .EQ. 0.0 .AND. Y .EQ. 0.0) THEN
           DANGLE = 0.0
        ELSE
           DANGLE = ATAN2( Z , Y )
        ENDIF

        IF( DANGLE.LT.0.0 ) DANGLE = 2.00 * PI + DANGLE
        RETURN
        END

