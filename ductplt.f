C ======================================================================
      SUBROUTINE DUCTPLT
C
C     Duct geometry 
C
C ======================================================================
      INCLUDE 'PUFCAV.INC'

      DELTAK = TWOPI / NBLADE
      DTHETA = DELTAK / REAL(MDUCT)

c ==================================================
      OPEN(57,FILE='duct.plt',STATUS='UNKNOWN')
c ==================================================

      NB22 = NBLADE/2
      NBREM = NBLADE - NB22

      WRITE (57,*) ' VARIABLES =X, Y, Z'
      WRITE (57,*) 
     %  'ZONE T="FIRST", I=', NDUCTP, ', J = ', 
     %  NB22*MDUCTP , ', K = 1 '
            
      DO KK = 1 , NB22
         DO J=1, MDUCTP
            DO I = 1 , NDUCTP
               XDDUM = XD(I,J)
               RRR = SQRT( YD(I,J)**2 + ZD(I,J)**2 )
               THETA = DANGLE(ZD(I,J), YD(I,J))
               THETA2 = THETA + DELTAK*(KK-1)
               YDDUM = RRR * COS(THETA2)
               ZDDUM = RRR * SIN(THETA2)
               WRITE(57,*) XDDUM,YDDUM,ZDDUM
            ENDDO
         ENDDO
      ENDDO

      WRITE (57,*) 
     %  'ZONE T="SECOND", I=', NDUCTP, ', J = ', NBREM*MDUCTP ,', K=1'
            
      DO KK = NB22+1,NBLADE
         DO J=1, MDUCTP
            DO I = 1 , NDUCTP 
               XDDUM = XD(I,J)
               RRR = SQRT( YD(I,J)**2 + ZD(I,J)**2 )
               THETA = DANGLE(ZD(I,J), YD(I,J))
               THETA2 = THETA + DELTAK*(KK-1)
               YDDUM = RRR * COS(THETA2)
               ZDDUM = RRR * SIN(THETA2)
               WRITE(57,*) XDDUM,YDDUM,ZDDUM
            ENDDO
         ENDDO
      ENDDO
   
      WRITE (57,*) 
     %  'ZONE T="DUCTWAKE-1", I=', NDWKP, ', J = ', 
     %  NB22*MDUCTP , ', K = 1 '
            
      DO KK = 1 , NB22
         DO J=1, MDUCTP
            DO I = 1 , NDWKP
               XDDUM = XDW(I,J)
               RRR = SQRT( YDW(I,J)**2 + ZDW(I,J)**2 )
               THETA = DANGLE(ZDW(I,J), YDW(I,J))
               THETA2 = THETA + DELTAK*(KK-1)
               YDDUM = RRR * COS(THETA2)
               ZDDUM = RRR * SIN(THETA2)
               WRITE(57,*) XDDUM,YDDUM,ZDDUM
            ENDDO
         ENDDO
      ENDDO

      WRITE (57,*) 
     %  'ZONE T="DUCTWAKE-2", I=', NDWKP, ', J = ',NBREM*MDUCTP ,',K=1'
            
      DO KK = NB22+1,NBLADE
         DO J=1, MDUCTP
            DO I = 1 , NDWKP
               XDDUM = XDW(I,J)
               RRR = SQRT( YDW(I,J)**2 + ZDW(I,J)**2 )
               THETA = DANGLE(ZDW(I,J), YDW(I,J))
               THETA2 = THETA + DELTAK*(KK-1)
               YDDUM = RRR * COS(THETA2)
               ZDDUM = RRR * SIN(THETA2)
               WRITE(57,*) XDDUM,YDDUM,ZDDUM
            ENDDO
         ENDDO
      ENDDO

      WRITE (57,*) 
     %  'ZONE T="DUCTDISK-1", I=', MDUCTP, ', J = ', 
     %  NB22*2 , ', K = 1 '
            
      DO KK = 1 , NB22
         DO J=1, 2
            DO I = 1 , MDUCTP
               XDDUM = XDWDK(I,J)
               RRR = SQRT( YDWDK(I,J)**2 + ZDWDK(I,J)**2 )
               THETA = DANGLE(ZDWDK(I,J), YDWDK(I,J))
               THETA2 = THETA + DELTAK*(KK-1)
               YDDUM = RRR * COS(THETA2)
               ZDDUM = RRR * SIN(THETA2)
               WRITE(57,*) XDDUM,YDDUM,ZDDUM
            ENDDO
         ENDDO
      ENDDO

      WRITE (57,*) 
     %  'ZONE T="DUCTDISK-2", I=', MDUCTP, ',J= ', NBREM*2 ,',K=1 '
            
      DO KK = NB22+1, NBLADE
         DO J=1, 2
            DO I = 1 , MDUCTP
               XDDUM = XDWDK(I,J)
               RRR = SQRT( YDWDK(I,J)**2 + ZDWDK(I,J)**2 )
               THETA = DANGLE(ZDWDK(I,J), YDWDK(I,J))
               THETA2 = THETA + DELTAK*(KK-1)
               YDDUM = RRR * COS(THETA2)
               ZDDUM = RRR * SIN(THETA2)
               WRITE(57,*) XDDUM,YDDUM,ZDDUM
            ENDDO
         ENDDO
      ENDDO

      CLOSE(57)
      
      RETURN
      END

C ======================================================================
      SUBROUTINE DUCTPLT_CHECK
C
C     Duct geometry 
C
C ======================================================================
      INCLUDE 'PUFCAV.INC'

      DELTAK = TWOPI / NBLADE
      DTHETA = DELTAK / REAL(MDUCT)

      NB22 = NBLADE/2
      NBREM = NBLADE - NB22
      
      WRITE (3367,*) 
     %  'ZONE T="FIRST", I=', NDUCTP, ', J = ', 
     %  NB22*MDUCTP , ', K = 1 '
      WRITE(3367,*) 'SOLUTIONTIME= ',icavt
c      WRITE(808+icavt-1,*) 'SOLUTIONTIME= ',icavt
            
      DO KK = 1 , NB22
         DO J=1, MDUCTP
            DO I = 1 , NDUCTP
               XDDUM = XD(I,J)
               RRR = SQRT( YD(I,J)**2 + ZD(I,J)**2 )
               THETA = DANGLE(ZD(I,J), YD(I,J))
               THETA2 = THETA + DELTAK*(KK-1)
               YDDUM = RRR * COS(THETA2)
               ZDDUM = RRR * SIN(THETA2)
               WRITE(3367,*) XDDUM,YDDUM,ZDDUM
            ENDDO
         ENDDO
      ENDDO

      WRITE (3367,*) 
     %  'ZONE T="SECOND", I=', NDUCTP, ', J = ', NBREM*MDUCTP ,', K=1'
      WRITE(3367,*) 'SOLUTIONTIME= ',icavt
            
      DO KK = NB22+1,NBLADE
         DO J=1, MDUCTP
            DO I = 1 , NDUCTP 
               XDDUM = XD(I,J)
               RRR = SQRT( YD(I,J)**2 + ZD(I,J)**2 )
               THETA = DANGLE(ZD(I,J), YD(I,J))
               THETA2 = THETA + DELTAK*(KK-1)
               YDDUM = RRR * COS(THETA2)
               ZDDUM = RRR * SIN(THETA2)
               WRITE(3367,*) XDDUM,YDDUM,ZDDUM
            ENDDO
         ENDDO
      ENDDO

      WRITE (3367,*) 
     %  'ZONE T="DUCTWAKE-1", I=', NDWKP, ', J = ', 
     %  NB22*MDUCTP , ', K = 1 '
      WRITE(3367,*) 'SOLUTIONTIME= ',icavt
            
      DO KK = 1 , NB22
         DO J=1, MDUCTP
            DO I = 1 , NDWKP
               XDDUM = XDW(I,J)
               RRR = SQRT( YDW(I,J)**2 + ZDW(I,J)**2 )
               THETA = DANGLE(ZDW(I,J), YDW(I,J))
               THETA2 = THETA + DELTAK*(KK-1)
               YDDUM = RRR * COS(THETA2)
               ZDDUM = RRR * SIN(THETA2)
               WRITE(3367,*) XDDUM,YDDUM,ZDDUM
            ENDDO
         ENDDO
      ENDDO

      WRITE (3367,*) 
     %  'ZONE T="DUCTWAKE-2", I=', NDWKP, ', J = ',NBREM*MDUCTP ,',K=1'
      WRITE(3367,*) 'SOLUTIONTIME= ',icavt
c      WRITE(808+icavt-1,*) 'SOLUTIONTIME= ',icavt
            
      DO KK = NB22+1,NBLADE
         DO J=1, MDUCTP
            DO I = 1 , NDWKP
               XDDUM = XDW(I,J)
               RRR = SQRT( YDW(I,J)**2 + ZDW(I,J)**2 )
               THETA = DANGLE(ZDW(I,J), YDW(I,J))
               THETA2 = THETA + DELTAK*(KK-1)
               YDDUM = RRR * COS(THETA2)
               ZDDUM = RRR * SIN(THETA2)
               WRITE(3367,*) XDDUM,YDDUM,ZDDUM
            ENDDO
         ENDDO
      ENDDO

c      WRITE (3367,*) 
c     %  'ZONE T="DUCTDISK-1", I=', MDUCTP, ', J = ', 
c     %  NB22*2 , ', K = 1 '
c      WRITE(3367,*) 'SOLUTIONTIME= ',icavt
c            
c      DO KK = 1 , NB22
c         DO J=1, 2
c            DO I = 1 , MDUCTP
c               XDDUM = XDWDK(I,J)
c               RRR = SQRT( YDWDK(I,J)**2 + ZDWDK(I,J)**2 )
c               THETA = DANGLE(ZDWDK(I,J), YDWDK(I,J))
c               THETA2 = THETA + DELTAK*(KK-1)
c               YDDUM = RRR * COS(THETA2)
c               ZDDUM = RRR * SIN(THETA2)
c               WRITE(3367,*) XDDUM,YDDUM,ZDDUM
c            ENDDO
c         ENDDO
c      ENDDO
c
c      WRITE (3367,*) 
c     %  'ZONE T="DUCTDISK-2", I=', MDUCTP, ',J= ', NBREM*2 ,',K=1 '
c      WRITE(3367,*) 'SOLUTIONTIME= ',icavt
c            
c      DO KK = NB22+1, NBLADE
c         DO J=1, 2
c            DO I = 1 , MDUCTP
c               XDDUM = XDWDK(I,J)
c               RRR = SQRT( YDWDK(I,J)**2 + ZDWDK(I,J)**2 )
c               THETA = DANGLE(ZDWDK(I,J), YDWDK(I,J))
c               THETA2 = THETA + DELTAK*(KK-1)
c               YDDUM = RRR * COS(THETA2)
c               ZDDUM = RRR * SIN(THETA2)
c               WRITE(3367,*) XDDUM,YDDUM,ZDDUM
c            ENDDO
c         ENDDO
c      ENDDO
      
      RETURN
      END



