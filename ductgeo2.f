C ===============================================================
      SUBROUTINE DUCTGEO2
C     
C     For IDGEO=1 Generate duct panel using duct pitch
C                 Read Pitch angle for duct panel      
C ===============================================================
      INCLUDE 'PUFCAV.INC'

      DIMENSION XID1(101),ETAD1(101)
C     DIMENSION DCUBICX1(400),DCUBICR1(400)
      DIMENSION XID2(101),ETAD2(101)
C     DCUBICX2(400),DCUBICR2(400)
C      DIMENSION SSX1(101),SSX2(101)

      CALL DUCTGEO3(0)
      
      DO N=1,NDUCT
         XG(N,1,1)=XD(N,2)
         XG(N,1,2)=YD(N,2)
         XG(N,1,3)=ZD(N,2)
         XG(N,2,1)=XD(N,1)
         XG(N,2,2)=YD(N,1)
         XG(N,2,3)=ZD(N,1)
         XG(N,3,1)=XD(N+1,1)
         XG(N,3,2)=YD(N+1,1)
         XG(N,3,3)=ZD(N+1,1)
         XG(N,4,1)=XD(N+1,2)
         XG(N,4,2)=YD(N+1,2)
         XG(N,4,3)=ZD(N+1,2)
      ENDDO

      CALL GEOM3D(NDUCT,XG,CHRLEPS,IER)

      DO I = 1 , NDUCT
         XXCP(I) = XCT(I,1)
         RRCP(I) = SQRT(XCT(I,2)**2+XCT(I,3)**2)
      ENDDO

      PITCH2 = DUCTPT

      NDUCTH = NDUCT/2
C      NDUCTH = NDFWD + NDAFT
      NDUCTHP = NDUCTH + 1
      NDUCTP = NDUCT + 1
C      NDUCT = NDUCTP - 1

C      MDUCT = MHBT
      MDUCTP = MDUCT + 1

      DELTAK = TWOPI / NBLADE
      DTHETA = DELTAK / REAL(MDUCT)

      DO M = 1 , MDUCTP
         DO N = 1 , NDUCTP
            XD(N,M) = XD(N,1)
         ENDDO
      ENDDO

C -- Find Min. Value

      XXM = XD(1,1)
      DO I = 1 , NDUCTP
         IF(XD(I,1) .LE. XXM) THEN
            XXM = XD(I,1)
            NN1 = I
         ENDIF
      ENDDO
            
      NN2 = NDUCTP - NN1 + 1

      DO N = 1 , NN1
         N1 = NN1 - N + 1         
         XID1(N) = XD(N1,1)
         ETAD1(N) = SQRT( YD(N1,1)**2+ZD(N1,1)**2)
      ENDDO

      DO N = 1 , NN2
         N2 = NN1 + N - 1         
         XID2(N) = XD(N2,1)
         ETAD2(N) = SQRT( YD(N2,1)**2+ZD(N2,1)**2)
      ENDDO
      
      DO M = 1 , MDUCTP
            THETA = DTHETA * REAL(M-1)
            YD(NN1,M) = ETAD1(1) * COS(THETA)
            ZD(NN1,M) = ETAD1(1) * SIN(THETA)
      ENDDO

      DO M = 1 , MDUCTP
         THETA = DTHETA * REAL(M-1)
         DO N = 2 , NN1
            N1 = NN1 - N + 1
            DXX = XID1(N) - XID1(N-1)
            DTHETA1 = DXX * TAN(PITCH2)
            THETA = THETA + DTHETA1
            YD(N1,M) = ETAD1(N) * COS(THETA)
            ZD(N1,M) = ETAD1(N) * SIN(THETA)
         ENDDO

         THETA = DTHETA * REAL(M-1)
         DO N = 2 , NN2
            N1 = NN1 + N - 1
            DXX = XID2(N) - XID2(N-1)
            DTHETA1 = DXX * TAN(PITCH2)
            THETA = THETA + DTHETA1
            YD(N1,M) = ETAD2(N) * COS(THETA)
            ZD(N1,M) = ETAD2(N) * SIN(THETA)
         ENDDO

      ENDDO

      RETURN
      END
