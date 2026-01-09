      SUBROUTINE RTABLE_DUCT(RULT,RHULT,RTBL,XTBL,XULT,RW,MM,
     1                  XID1,ETAD1,NN1)
C***********************************************************************
C     RTABLE: Radii TABLE of the transition wake geometry
C       --- Generate a table of radii for contracted transition wake
C

      DIMENSION RTBL(11,MM),XTBL(11),RW(MM)
      DIMENSION XID1(101),ETAD1(101),DCUBICX1(400)

      ZERO=0.0
      MM1=MM-1

      DO M=1,MM
         RTBL(1,M)=RW(M)
      ENDDO

      CALL UGLYDK(NN1,1,1,XID1,ETAD1,ZERO,ZERO,DCUBICX1)

      XULTUN=XULT
      XULTUN=AMAX1(XULTUN,0.3)
      XTBL(1)=ZERO

      DO N=2,11

         RTBL(N,1) = RHULT

         XTBL(N)=0.1*(N-1)*XULT

         IF(XTBL(N) .LE. XID1(NN1)) THEN
            CALL EVALDKs(NN1,1,XID1,XTBL(N),RRR,DCUBICX1)
            IF(RRR .LE. RW(MM)) THEN
               RTBL(N,MM) = RRR
            ELSE
               RTBL(N,MM) = RW(MM)
            ENDIF         
         ELSE
            IF(ETAD1(NN1) .LE. RW(MM)) THEN
               RTBL(N,MM) = ETAD1(NN1)
            ELSE
               RTBL(N,M) = RW(M)
            ENDIF
         ENDIF

         DO M=2,MM1
            D1=SQRT(RTBL(N,M-1)**2+RTBL(1,M)**2-RTBL(1,M-1)**2)
     *        -RTBL(N,M-1)
            IF(M.GT.MM/2+2) THEN
               D1=AMIN1( D1,(RTBL(N,M-1)/RTBL(N,MM))**2
     *           *(RTBL(N,MM)-RTBL(N,M-1)) )
            END IF 
            RTBL(N,M)=RTBL(N,M-1)+D1
         ENDDO
      ENDDO

      RETURN
C))))))))))))))))))))) End of subroutine RTABLE ((((((((((((((((((((((((
      END
