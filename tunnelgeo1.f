C ======================================================================
      SUBROUTINE TUNNELGEO1
C
C     Tunnel geometry with circular cross section
C
C ======================================================================
      INCLUDE 'PUFCAV.INC'
      
      DIMENSION XHUB(100),YHUB(100),ZHUB(100)

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

      DELAXa=0.0
      DELAXf=0.0

      if(NAXa.gt.0) DELAXa=(XAFT-Xmida)/FLOAT(NAXa)
      if(NAXf.gt.0) DELAXf=(-XHUB(1)-Xmidf)/FLOAT(NAXf)
      if(Naxm.gt.0) DELAXm= (Xmidf+Xmida) / REAL(NAXM)

      DELTAK = TWOPI / NBLADE
      DTHETA = DELTAK / REAL(MTUNEL)

c ===============================
c     construct TUNNEL WALL
c ===============================

c =============================
c     Left Lid (upstream)
c =============================
     
      write(*,*) ' LEFT LID  : '

      DO J=1,Mtunel1
         
         IF(IHUB .EQ. 1) THEN
            THETA = DANGLE(ZHUB(J),YHUB(J))
            R1 = SQRT(YHUB(J)**2+ZHUB(J)**2)
            DR = (TUNRAD-R1) / REAL(NSIDE)
         ELSE
            THETA = DTHETA * REAL(J-1)
            DR = TUNRAD / REAL(NSIDE)
         ENDIF
      
         DO I=1,NSIDE1
            IF(IHUB .EQ. 0) THEN
               RRR = DR * REAL(I-1)
            ELSE
               RRR = R1 + DR*REAL(I-1)
            ENDIF

            XTUN(I,J) = XHUB(J)
            YTUN(I,J) = RRR * COS(THETA)
            ZTUN(I,J) = RRR * SIN(THETA)
         ENDDO
      ENDDO

      write(*,*) ' MID PART : '
  
      I1 = NSIDE

      DO J=1,MTUNEL1
         XBI=XTUN(NSIDE1,J) - DELAXF
         THETA = DANGLE(ZTUN(NSIDE1,J),YTUN(NSIDE1,J))
         DO I = 1 , NAX1
            IF(I .LE. NAXF1) THEN
               XBI = XBI + DELAXF
            ELSEIF(I .GT. NAXMD1) THEN
               XBI = XBI + DELAXA
            ELSE
               XBI = XBI + DELAXM
            ENDIF

            XTUN(NSIDE+I,J) = XBI
            YTUN(NSIDE+I,J) = TUNRAD * COS(THETA)
            ZTUN(NSIDE+I,J) = TUNRAD * SIN(THETA)
         ENDDO
      ENDDO
            
c ===============================     
c     right lid (downstream)
c ===============================
     
      write(*,*) ' RIGHT LID  : '

      I1 = NAX1 + NSIDE

      DO J=1,Mtunel1
         
         THETA = DANGLE(ZTUN(I1,J),YTUN(I1,J))        
         DR = TUNRAD / REAL(NSIDE)
         
         DO I=1, NSIDE1
            RRR = DR * REAL(NSIDE-I+1)
            
            XTUN(I1+I-1,J) = XTUN(I1,J)
            YTUN(I1+I-1,J) = RRR * COS(THETA)
            ZTUN(I1+I-1,J) = RRR * SIN(THETA)
         ENDDO
      ENDDO
      
      RETURN
      END

