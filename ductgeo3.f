C ====================================================================
      SUBROUTINE DUCTGEO3(IDPT) 
C
C     For IDGEO=1 Generate duct panel using duct pitch
C                 Read Pitch angle for duct panel  
C
C     IDPT = 1 : No. of panels along circumferential direction = MHBT
C            0 : No. of panels along circumferential direction = 1
C                This option is used for DUCTGEO2 routine
C ====================================================================
      INCLUDE 'PUFCAV.INC'

      DIMENSION XID1(101),ETAD1(101),DCUBICX1(400),DCUBICR1(400)
      DIMENSION XID2(101),ETAD2(101),DCUBICX2(400),DCUBICR2(400)
      DIMENSION SSX1(101),SSX2(101)

      NN1 = NDDAT/2 + 1
   
      DO N = 1 , NN1
         
         N1 = NN1 - N + 1
         N2 = NN1 + N - 1
         
         XID1(N) = XID(N1)
         ETAD1(N) = ETAD(N1)
         
         XID2(N) = XID(N2)
         ETAD2(N) = ETAD(N2)
         
         IF(N .EQ. 1) THEN
            SSX1(1) = 0.0
            SSX2(1) = 0.0
         ELSE
            
            DDX1 = SQRT((XID1(N)-XID1(N-1))**2
     %           +(ETAD1(N)-ETAD1(N-1))**2)
            DDX2 = SQRT((XID2(N)-XID2(N-1))**2
     %           +(ETAD2(N)-ETAD2(N-1))**2)
            
            SSX1(N) = SSX1(N-1) + DDX1
            SSX2(N) = SSX2(N-1) + DDX2
         ENDIF
         
      ENDDO
      
      CALL UGLYDK(NN1,1,1,SSX1,XID1,ZERO,ZERO,DCUBICX1)
      CALL UGLYDK(NN1,1,1,SSX1,ETAD1,ZERO,ZERO,DCUBICR1)
      CALL UGLYDK(NN1,1,1,SSX2,XID2,ZERO,ZERO,DCUBICX2)
      CALL UGLYDK(NN1,1,1,SSX2,ETAD2,ZERO,ZERO,DCUBICR2)
      
C      NDUCTH = NDFWD + NDAFT
      NDUCTH = NDUCT / 2
      NDUCTHP = NDUCTH + 1
      NDUCTP = NDUCT + 1
C      NDUCT = NDUCTP - 1

C      MDUCT = MHBT

      MMTMP = MDUCT
C     IF(IDPT .EQ. 0) MDUCT=1

      MDUCTP = MDUCT + 1

      DELTAK = TWOPI / NBLADE
      DTHETA = DELTAK / REAL(MMTMP)

      DTN = PI / REAL(NDUCTH)
      
      DELX1 = SSX1(NN1)
      DELX2 = SSX2(NN1)
      
      DO M = 1 , MDUCTP
         THETA = DTHETA*REAL(M-1)            
         XD(NDUCTHP,M) = XID1(1)
         YD(NDUCTHP,M) = ETAD1(1) * COS(THETA)
         ZD(NDUCTHP,M) = ETAD1(1) * SIN(THETA)
      ENDDO
      
      DO N = 1 , NDUCTH
         XDI = 0.5*DELX1 * (1. - COS(DTN*(N)))
         
         CALL EVALDKs(NN1,1,SSX1,XDI,XDD1,DCUBICX1)
         CALL EVALDKs(NN1,1,SSX1,XDI,RRR1,DCUBICR1)
         
         XDI = 0.5*DELX2 * (1. - COS(DTN*(N)))
         
         CALL EVALDKs(NN1,1,SSX2,XDI,XDD2,DCUBICX2)
         CALL EVALDKs(NN1,1,SSX2,XDI,RRR2,DCUBICR2)
         
         DO M = 1 , MDUCTP
            THETA = DTHETA*REAL(M-1)

            XD(NDUCTHP-N,M) = XDD1
            YD(NDUCTHP-N,M) = RRR1 * COS(THETA)
            ZD(NDUCTHP-N,M) = RRR1 * SIN(THETA)

            XD(NDUCTHP+N,M) = XDD2
            YD(NDUCTHP+N,M) = RRR2 * COS(THETA)
            ZD(NDUCTHP+N,M) = RRR2 * SIN(THETA)
         ENDDO
      ENDDO

      RETURN
      END
