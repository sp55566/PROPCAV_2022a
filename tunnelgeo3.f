C ======================================================================
      SUBROUTINE TUNNELGEO3
C
C     Tunnel geometry with circular cross section
C
C ======================================================================
      INCLUDE 'PUFCAV.INC'
      
      DIMENSION XHUB(100),YHUB(100),ZHUB(100)
      DIMENSIOn XBB(100,2),YBB(100,2),ZBB(100,2)
      DIMENSION THE1(50), THE2(50)

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
         
         PPD = XPI(NX)
         PPDHUB = XPI(1)

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
C         NAXM = NHP-1
C         NAXM1 = NHP
          NAXM = (NHP-1)/2
          NAXM1 = NAXM+1
         PITCH2 = 0.5 * PI - ATAN(PPD/PI)
      ENDIF

      NAXa1 = NAXa + 1
      Naxf1 = Naxf + 1
      NAXMD1 = NAXF + NAXM + 1
      NAX =NAXf+NAXm+NAXa
      NAX1=NAX+1
      Mtunel1 = Mtunel + 1
      Nside1 = Nside + 1
      
      NAXT = NAX + 2 * NSIDE
      NAXT1 = NAX1 + 2 * NSIDE

      DELTAK = TWOPI / NBLADE
      DTHETA = DELTAK / REAL(MTUNEL)

c ===============================
c     construct TUNNEL WALL
c ===============================

      write(*,*) ' Mid Part : '

      I11 = NSIDE + NAXF

      DO I = 1 , NAXM1
         JJ = 2 * I - 1
         THETA11 = ATAN2(ZBB(JJ,1),YBB(JJ,1))
         THETA22 = ATAN2(ZBB(JJ,2),YBB(JJ,2)) + DELTAK
         DTHETA = (THETA22 - THETA11) / REAL(MTUNEL)
         DDXX = (XBB(JJ,2)-XBB(JJ,1))/REAL(MTUNEL) 
         DO J = 1, MTUNEL1
            THETA = THETA11 + DTHETA*REAL(J-1)
            XTUN(I11+I,J) = XBB(I,1) + DDXX * REAL(J-1)
            YTUN(I11+I,J) = TUNRAD * COS(THETA)
            ZTUN(I11+I,J) = TUNRAD * SIN(THETA)
         ENDDO
      ENDDO

      write(*,*) ' FORWARD PART : '

      DTN = 0.5 * PI / REAL(NAXF)

      I11 = NSIDE + NAXF1
      DO J=1,MTUNEL1
         THETA = DANGLE(ZTUN(I11,J),YTUN(I11,J))
         DELX = -XHUB(1) + XTUN(I11,J)
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

C      XBI = XBB(NHP)

      I11 = NSIDE + NAXF + NAXM1
      DO J=1,MTUNEL1
         THETA = DANGLE(ZTUN(I11,J),YTUN(I11,J))
         DELX = XAFT - XTUN(I11,J)
         DO I = 1, NAXA
            XDI = DELX * (1.-COS(DTN*I))
C            XTUN(I11+I,J) = XBI + XDI
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
     
      write(*,*) ' LEFT LID  : '

      DO M = 1 , MTUNEL1
         THE1(M) = DANGLE(ZTUN(NSIDE1,M),YTUN(NSIDE1,M))   
         THE2(M) = DANGLE(ZHUB(M),YHUB(M))
      ENDDO

      DFER1 = THE1(MTUNEL1) - THE1(1)
      DFER2 = THE2(MTUNEL1) - THE2(1)

      IF(PPDHUB .GE. 0.0) THEN
         IF(DFER1 .LT. 0.0) THEN
            DO M = MTUNEL, 1 , -1
               IF(THE1(M+1) .LT. THE1(M)) THEN
                  THE1(M) = THE1(M) - 2. * PI
               ENDIF
            ENDDO
         ENDIF
      ELSE
         IF(DFER1 .GT. 0.0) THEN
            DO M = MTUNEL, 1 , -1
               IF(THE1(M+1) .GT. THE1(M)) THEN
                  THE1(M) = 2. * PI - THE1(M)
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      
      IF(PPD .GE. 0.0) THEN
         IF(DFER2 .LT. 0.0) THEN
            DO M = MTUNEL, 1 , -1
               IF(THE2(M+1) .LT. THE2(M)) THEN
                  THE2(M) = THE2(M) - 2. * PI
               ENDIF
            ENDDO
         ENDIF
      ELSE
         IF(DFER2 .GT. 0.0) THEN
            DO M = MTUNEL, 1 , -1
               IF(THE2(M+1) .GT. THE2(M)) THEN
                  THE2(M) = 2. * PI - THE2(M)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      DFER1 = THE1(1) - THE2(1)

      IF(DFER1 .GT. 1.5*PI) THEN
         DO M = 1 , MTUNEL1
            THE2(M) = 2.*PI + THE2(M) 
         ENDDO
      ELSEIF(DFER1 .LT. -1.5*PI) THEN
         DO M = 1, MTUNEL1
            THE2(M) = THE2(M) - 2.*PI
         ENDDO
      ENDIF


      DO J = 1 , Mtunel1
         DTHETA12 = (THE1(J) - THE2(J))/REAL(NSIDE)
         RRR1 = SQRT(YHUB(J)**2 + ZHUB(J)**2) 
         DR = (TUNRAD-RRR1) / REAL(NSIDE)
         DO I=NSIDE1,1,-1

            RRR = TUNRAD - DR * REAL(NSIDE1-I)

            THETA12 = THE1(J) - DTHETA12 * (NSIDE1 - I)
            XTUN(I,J) = XTUN(NSIDE1,J)
            YTUN(I,J) = RRR * COS(THETA12)
            ZTUN(I,J) = RRR * SIN(THETA12)
         ENDDO
      ENDDO

c ===============================     
c     right lid (downstream)
c ===============================
     
      write(*,*) ' RIGHT LID  : '

      I11 = NSIDE+NAX+1
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

      RETURN
      END

