C ======================================================================
      SUBROUTINE TUNNELGEO2
C
C     Tunnel geometry with circular cross section
C
C ======================================================================
      INCLUDE 'PUFCAV.INC'
      
      DIMENSION XHUB(30),YHUB(30),ZHUB(30)
      DIMENSIOn XBB(100,2),YBB(100,2),ZBB(100,2)

      PI = ACOS(-1.)
      TWOPI = 2.0 * PI
      
      IF(IHUB .NE. 0) THEN
         Mtunel = MHBT
         Mtunel1 = Mtunel+1
         DO M = 1 , MTUNEL1
            XHUB(M) = XH(1,M)
            YHUB(M) = YH(1,M)
            ZHUB(M) = ZH(1,M)
         ENDDO
         PPDHUB = XPI(1)
      ENDIF
      
      PPD = XPI(NX)
      
      DO N = 1 , NHP
         I1 = NHP - N + 1
         I2 = NHP + N - 1
         IF(TUNRAD .EQ. 1.0) THEN
            XBB(N,1) = XB(I2,MRP)
            XBB(N,2) = XB(I1,MRP)
            YBB(N,1) = YB(I2,MRP)
            YBB(N,2) = YB(I1,MRP)
            ZBB(N,1) = ZB(I2,MRP)
            ZBB(N,2) = ZB(I1,MRP)
         ELSE
            XBB(N,1) = 0.5 * (XB(I1,MRP)+XB(I2,MRP))
            YBB(N,1) = 0.5 * (YB(I1,MRP)+YB(I2,MRP))
            ZBB(N,1) = 0.5 * (ZB(I1,MRP)+ZB(I2,MRP))
            XBB(N,2) = XBB(N,1)
            YBB(N,2) = YBB(N,1)
            ZBB(N,2) = ZBB(N,1)
         ENDIF
      ENDDO

      NAXM = NHP-1
      NAXM1 = NHP
      PITCH2 = 0.5 * PI - ATAN(PPD/PI)
      
      NAXa1 = NAXa + 1
      Naxf1 = Naxf + 1
      NAXMD1 = NAXF + NAXM + 1
      NAX =NAXf+NAXm+NAXa
      NAX1=NAX+1
      Mtunel1 = Mtunel + 1
      Nside1 = Nside + 1

C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)
C      NAXT = NAX + 2*NSIDE
C      NAXT1 = NAX1 + 2*NSIDE
      NAXT = NAX + NSIDE
      NAXT1 = NAX1 + NSIDE
C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)

      DELTAK = TWOPI / NBLADE
      DTHETA = DELTAK / REAL(MTUNEL)

c ===============================
c     construct TUNNEL WALL
c ===============================

      write(*,*) ' Mid Part : '

C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)
C      I11 = NSIDE + NAXF
      I11 = NAXF
C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)

      DO I = 1 , NAXM1
         THETA11 = ATAN2(ZBB(I,1),YBB(I,1))
         THETA22 = ATAN2(ZBB(I,2),YBB(I,2)) + DELTAK
         DTHETA = (THETA22 - THETA11) / REAL(MTUNEL)
         DDXX = (XBB(I,2)-XBB(I,1))/REAL(MTUNEL) 
         DO J = 1, MTUNEL1
            THETA = THETA11 + DTHETA*REAL(J-1)
            XTUN(I11+I,J) = XBB(I,1) + DDXX * REAL(J-1)
            YTUN(I11+I,J) = TUNRAD * COS(THETA)
            ZTUN(I11+I,J) = TUNRAD * SIN(THETA)
         ENDDO
      ENDDO

      write(*,*) ' FORWARD PART : '

      DTN = 0.5 * PI / REAL(NAXF)

C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)
C      I11 = NSIDE + NAXF1
      I11 = NAXF1
C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)

      DO J=1,MTUNEL1
         THETA = DANGLE(ZTUN(I11,J),YTUN(I11,J))

         IF(IHUB .NE. 0) THEN
            DELX = -XHUB(1) + XTUN(I11,J)
         ELSE
            DELX = -XFOR + XTUN(I11,J)
         ENDIF

         DO I = 1 , NAXF
            XDI = DELX* (1.-COS(DTN*(I)))
            XTUN(I11-I,J) = XTUN(I11,J) - XDI
            
            DTHETA = XDI * TAN(PITCH2)/TUNRAD
            THETA2 = THETA - DTHETA
            
            YTUN(I11-I,J) = TUNRAD * COS(THETA2)
            ZTUN(I11-I,J) = TUNRAD * SIN(THETA2)
         ENDDO
      ENDDO

      write(*,*) ' AFT PART : '

      DTN = 0.5 * PI / REAL(NAXA)

C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)
C      I11 = NSIDE + NAXF + NAXM1
      I11 = NAXF + NAXM1
C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)

      DO J=1,MTUNEL1
         THETA = DANGLE(ZTUN(I11,J),YTUN(I11,J))
         DELX = XAFT - XTUN(I11,J)
         DO I = 1, NAXA
            XDI = DELX * (1.-COS(DTN*I))
            XTUN(I11+I,J) = XTUN(I11,J) + XDI
            DTHETA = XDI * TAN(PITCH2) / TUNRAD
            THETA2 = THETA + DTHETA
            YTUN(I11+I,J) = TUNRAD * COS(THETA2)
            ZTUN(I11+I,J) = TUNRAD * SIN(THETA2)
         ENDDO
      ENDDO


c =============================
c     Left Lid (upstream)
c =============================

C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)
C      write(*,*) ' LEFT LID  : '
      write(*,*) ' NO LEFT LID  : '

C      IF(IHUB .NE. 0) THEN
C         RRR1 = SQRT(YHUB(1)**2 + ZHUB(1)**2)
C      ELSE
C         RRR1 = 0.0
C      ENDIF

C      DO J = 1 , Mtunel1
C         THETA = DANGLE(ZTUN(NSIDE1,J),YTUN(NSIDE1,J))   
C         DR = (TUNRAD-RRR1) / REAL(NSIDE)
C         DO I=NSIDE1,1,-1
C            RRR = TUNRAD - DR * REAL(NSIDE1-I)
C            XTUN(I,J) = XTUN(NSIDE1,J)
C            YTUN(I,J) = RRR * COS(THETA)
C            ZTUN(I,J) = RRR * SIN(THETA)
C         ENDDO
C      ENDDO
C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)

c ===============================     
c     right lid (downstream)
c ===============================
     
      write(*,*) ' RIGHT LID  : '

C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)
C      I11 = NSIDE+NAX+1
      I11 = NAX+1
C --- S.H.CHANG 03/31/2010  (REMOVING THE FRONT LID)

      DO J=1,Mtunel1
         THETA = DANGLE(ZTUN(I11,J),YTUN(I11,J))
         DR = TUNRAD / REAL(NSIDE)
         DO I=1, NSIDE
            RRR = TUNRAD - DR * REAL(I)
            XTUN(I11+I,J) = XTUN(I11,J)
            YTUN(I11+I,J) = RRR * COS(THETA)
            ZTUN(I11+I,J) = RRR * SIN(THETA)
         ENDDO
      ENDDO

      DO I = 1,NAXT
         WRITE(789,*) XTUN(I,1),YTUN(I,1),ZTUN(I,1)
         WRITE(790,*) I,XTUN(I,MTUNEL)
      END DO 

      RETURN
      END
