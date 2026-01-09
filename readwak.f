      SUBROUTINE READWAK
************************************************************************
*     READWAK: read and extrapolate data from *.wak                    *
*                                                                      *
*     Date     Comment or Revision                                     *
*     -------- -------------------                                     *
*     JY090799 Moved this part of the subroutine from conpt.f to here. *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

C-----------------------------------------------------------------------
C     Read wake harmonics from a wake input file at zero time step
C      notice that I = 1--axial, 2--radial, 3--tangential
C-----------------------------------------------------------------------
  401 CONTINUE
      OPEN(UNIT=1,FILE=WKFILE,STATUS='OLD',ERR=402)
      GOTO 403
  402 CONTINUE
C.....Enter inflow wake filename........................................
      WRITE(*,*)' ****WAKE FILE DOES NOT EXIST, TRY AGAIN!****'
      WRITE(*,*) ' PROPCAV> ENTER INFLOW WAKE FILENAME: '
      READ(*,*) WKFILE
      WRITE(*,*) ' '
      GOTO 401
  403 CONTINUE

      READ(1,'(A)') ADUMMY
      READ(1,*) NWKCOE,NHARM(1),NHARM(2),NHARM(3)

C-----------------------------------------------------------------------
C     NWKCOE was limited to 9 because there was a dimensioning error 
C     for variable XRW in PUFCAVB.INC.  I've already fixed it so that
C     NWKCOE is limited to 15.                                  JY090799
C-----------------------------------------------------------------------
CJY      IF(NWKCOE.GT.9)THEN
      IF(NWKCOE.GT.15)THEN
         WRITE(*,5000) NWKCOE
 5000    FORMAT('You have specified',I3,' input radii in wake file.')
         WRITE(*,5020)
CJY 5020    FORMAT('This is too many!  Try again.')
 5020    FORMAT('This is too many! It must be less than 15. Try again.')
         STOP
      ELSEIF(NWKCOE.LT.5)THEN
         WRITE(*,5000) NWKCOE
         WRITE(*,5040)
 5040    FORMAT('This is not enough!  Try again.')
         STOP
      ENDIF
      READ(1,'(A)') ADUMMY
      READ(1,*) (XRW(M),M=1,NWKCOE)
      DO 46 I=1,3
        DO 45 K=1,2
          READ(1,'(A)') ADUMMY
          DO 42 J=1,NHARM(I)      
            READ(1,*) (WAKHAR(M,J,K,I), M=1,NWKCOE)        
            CALL UGLYDK(NWKCOE,1,1,XRW,WAKHAR(1,J,K,I),0.0,0.0,
     *                  XWCUB(1,J,K,I))

C-----------------------------------------------------------------------
C     Extrapolate the mean axial, radial, and tangential velocities
C     from the zeroth harmonic of the cosine components.        JY090799
C-----------------------------------------------------------------------
            IF(K.EQ.1.AND.J.EQ.1) THEN
               IF(I.EQ.1) THEN
                  CALL EVALDK(NWKCOE,NX,XRW,XR,XVA,XWCUB(1,J,K,I))
               ELSE IF(I.EQ.2) THEN
                  CALL EVALDK(NWKCOE,NX,XRW,XR,XVR,XWCUB(1,J,K,I))
               ELSE
                  CALL EVALDK(NWKCOE,NX,XRW,XR,XVT,XWCUB(1,J,K,I))
               END IF
            END IF

   42     CONTINUE
   45   CONTINUE
   46 CONTINUE 
      READ(1,'(A)') ADUMMY
      CLOSE(1)                                                      

      RETURN
      END

