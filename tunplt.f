C ======================================================================
      SUBROUTINE TUNPLT
C
C     Tunnel geometry with circular cross section
C
C ======================================================================
      INCLUDE 'PUFCAV.INC'

      DELTAK = TWOPI / NBLADE
      DTHETA = DELTAK / REAL(MTUNEL)

      MTUNEL1 = MTUNEL + 1
      NAXT1 = NAXT + 1

      IF(NBLADE .GT. 1) THEN
         NBLDD = NBLADE/2
         IF(NBLDD*2 .NE. NBLADE) NBLDD = NBLDD + 1
         NBLDD2 = NBLADE - NBLDD
      ELSE
         NBLDD = 1
         NBLDD2 = 0
      ENDIF

c ==================================================
      OPEN(57,FILE='tunnel.plt',STATUS='UNKNOWN')
c ==================================================

      WRITE (57,*) ' VARIABLES =X, Y, Z'

      WRITE (57,*) 
     %  'ZONE T="FIRST", I=', Naxt1, ', J = ', 
     %   NBLDD*MTUNEL1 , ', K = 1 '
            
      DO KK = 1 , NBLDD
         DO J=1, MTUNEL1
            DO I = 1 , NAXT1
               XTDUM = XTUN(I,J)
               RRR = SQRT( YTUN(I,J)**2 + ZTUN(I,J)**2 )
               THETA = DANGLE(ZTUN(I,J), YTUN(I,J))
               THETA2 = THETA + DELTAK*(KK-1)
               YTDUM = RRR * COS(THETA2)
               ZTDUM = RRR * SIN(THETA2)
               WRITE(57,*) XTDUM,YTDUM,ZTDUM
            ENDDO
         ENDDO
      ENDDO

      IF(NBLDD2 .NE. 0) THEN
         WRITE (57,*) 
     %   'ZONE T="2ND", I=', Naxt1, ', J=', MTUNEL1*NBLDD2 ,', K = 1 '
            
         DO KK = NBLDD+1, NBLADE
            DO J=1, MTUNEL1
               DO I = 1 , NAXT1
                  XTDUM = XTUN(I,J)
                  RRR = SQRT( YTUN(I,J)**2 + ZTUN(I,J)**2 )
                  THETA = DANGLE(ZTUN(I,J), YTUN(I,J))
                  THETA2 = THETA + DELTAK*(KK-1)
                  YTDUM = RRR * COS(THETA2)
                  ZTDUM = RRR * SIN(THETA2)
                  WRITE(57,*) XTDUM,YTDUM,ZTDUM
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      CLOSE(57)
      
      RETURN
      END


C ======================================================================
      SUBROUTINE TUNPOT
C
C     Potential on Tunnel
C
C ======================================================================
      INCLUDE 'PUFCAV.INC'

      DELTAK = TWOPI / NBLADE
      DTHETA = DELTAK / REAL(MTUNEL)

c ==================================================
      OPEN(57,FILE='tunpot.plt',STATUS='UNKNOWN')
c ==================================================

      WRITE (57,*) 'VARIABLES="X","Y","Z","POT"'
            
      DO KK = 1 , NBLADE-1
      WRITE (57,*) 'ZONE T="',KK,'", I=', Naxt, ', J = ', MTUNEL
         DO J=1, MTUNEL
            DO I = 1 , NAXT
               IJ = INDEXTN(I,J)
               XTDUM = XCT(IJ,1)
               RRR = SQRT( XCT(IJ,2)**2 + XCT(IJ,3)**2 )
               THETA = DANGLE(XCT(IJ,3), XCT(IJ,2))
               THETA2 = THETA + DELTAK*(KK-1)
               YTDUM = RRR * COS(THETA2)
               ZTDUM = RRR * SIN(THETA2)
               WRITE(57,*) XTDUM,YTDUM,ZTDUM,POT(IJ)
            ENDDO
         ENDDO
      ENDDO

      WRITE (57,*) 
     %  'ZONE T="SECOND", I=', Naxt, ', J = ', MTUNEL 
            
      DO J=1, MTUNEL
         DO I = 1 , NAXT
            IJ = INDEXTN(I,J)
            XTDUM = XCT(IJ,1)
            RRR = SQRT( XCT(IJ,2)**2 + XCT(IJ,3)**2 )
            THETA = DANGLE(XCT(IJ,3), XCT(IJ,2))
            THETA2 = THETA + DELTAK*(NBLADE-1)
            YTDUM = RRR * COS(THETA2)
            ZTDUM = RRR * SIN(THETA2)
            WRITE(57,*) XTDUM,YTDUM,ZTDUM,POT(IJ)
         ENDDO
      ENDDO

      CLOSE(57)
      
      RETURN
      END

