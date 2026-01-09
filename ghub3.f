      SUBROUTINE GHUB3
C***********************************************************************
C     GHUB: Geometry of the HUB
C      --- Generate a hub geometry
C
      INCLUDE 'PUFCAV.INC'
      DIMENSION TH(NHPZ,MHPZ)

C-----------------------------------------------------------------------
C     Parameters for the hub geometry 
C-----------------------------------------------------------------------
      NHBB=NH
      DO 10 N=1,NSW(1)
         IF(XW(N,1) .GT. XHBFD) THEN
            NHBDT=N-1
            NHBU1=NHBU+1
            NHBUB=NHBU+NH
            NHBX=NHBUB+NHBDT
            MHBT1=MHBT+1
            GO TO 20
         END IF
10    CONTINUE

      NHBDT=NSW(1)
      NHBU1=NHBU+1
      NHBUB=NHBU+NH
      NHBX=NHBUB+NHBDT
      MHBT1=MHBT+1

20    CONTINUE

C-----------------------------------------------------------------------
C     Hub geometry in the near region of the blade
c       (blade section as a bi-sector)
C-----------------------------------------------------------------------
      PHIA=ATAN( PITCH(1)/(PI*RHUB) )
      PSI1=PHIA
      PSI2=ATAN2( XB(NHP+1,1)-XB(NHP,1),(THR(NHP+1)-THR(NHP))*RHUB )
      PSI0=HALF*(PI-PSI1-PSI2)
      DEL=RHUB*DELK*SIN(PHIA+PSI0)/(MHBT*SIN(PHIA))
      PSI=PSI0
      DO 30 N=1,NHP
         XH(NHBU+N,1)=XB(NHBB+N,1)
         TH(NHBU+N,1)=THR(NHBB+N)
         YH(NHBU+N,1)=YB(NHBB+N,1)
         ZH(NHBU+N,1)=ZB(NHBB+N,1)
         IF(N.LE.NHBB/2) THEN
            PSI=HALF*PSI
         ELSE
            PSI=ZERO
         END IF
         XH(NHBU+N,2)=XH(NHBU+N,1)-DEL*SIN(PSI)
         TH(NHBU+N,2)=TH(NHBU+N,1)+DEL*COS(PSI)/RHUB
         YH(NHBU+N,2)=RHUB*COS(TH(NHBU+N,2))
         ZH(NHBU+N,2)=RHUB*SIN(TH(NHBU+N,2))
30    CONTINUE
      DO 40 N=1,NHP
         XH(NHBU+N,MHBT1)=XB(NHP+1-N,1)
         TH(NHBU+N,MHBT1)=THR(NHP+1-N) + DELK
         YH(NHBU+N,MHBT1)=RHUB*COS(TH(NHBU+N,MHBT1))
         ZH(NHBU+N,MHBT1)=RHUB*SIN(TH(NHBU+N,MHBT1))
40    CONTINUE
C-----------------------------------------------------------------------
C     Hub geometry at upstream of the blade
C       --- half cosine spacing axially
C       --- constant spacing circumferentially
C-----------------------------------------------------------------------
      DTN=HALFPI/NHBU
      DO 50 N=1,NHBU
         XH(N,1)=XB(NHP,1) - XHBU*( 1.0-SIN(DTN*(N-1)) )
         TH(N,1)=THR(NHP)-PI/PITCH(1)*(XB(NHP,1)-XH(N,1))
         YH(N,1)=RHUB*COS(TH(N,1))
         ZH(N,1)=RHUB*SIN(TH(N,1))
         XH(N,2)=XH(1,1) + (XH(NHBU1,2)-XH(1,1))* SIN( DTN*(N-1) )
         TH(N,2)=TH(NHBU1,2)-PI/PITCH(1)*(XH(NHBU1,2)-XH(N,2))
         YH(N,2)=RHUB*COS(TH(N,2))
         ZH(N,2)=RHUB*SIN(TH(N,2))
         XH(N,MHBT1)=XH(N,1)
         TH(N,MHBT1)=TH(N,1)+DELK
         YH(N,MHBT1)=RHUB*COS(TH(N,MHBT1))
         ZH(N,MHBT1)=RHUB*SIN(TH(N,MHBT1))
50    CONTINUE
      DO 70 N=1,NHBUB
         DTH=TH(N,MHBT1)-TH(N,2)
         DO 60 M=3,MHBT
            TH(N,M)=TH(N,2)+ DTH/(MHBT-1)*(M-2)
            XH(N,M)=XH(N,2)+(XH(N,MHBT1)-XH(N,2))*(TH(N,M)-TH(N,2))/DTH
            YH(N,M)=RHUB*COS(TH(N,M))
            ZH(N,M)=RHUB*SIN(TH(N,M))
60       CONTINUE
70    CONTINUE
C-----------------------------------------------------------------------
C     Hub geometry at downstream of the blade
C-----------------------------------------------------------------------
      DELT1=TH(NHBUB+1,2)-TH(NHBUB+1,1)
      DO 90 N=1,NHBDT+1
         NN=NHBUB+N
         XH(NN,1)=XW(N,1)
         THW=ATAN2(ZW(N,1),YW(N,1))
         TH(NN,1)=THW
         RR=SQRT(YW(N,1)**2+ZW(N,1)**2)
         IF(IREADW.EQ.1) RR=RHUB !Wake panel input method. S.N.KIM 2018.
         YH(NN,1)=YW(N,1)
         ZH(NN,1)=ZW(N,1)
         IF(IREADW.EQ.1) YH(NN,1)=RR*COS(TH(NN,1)) !Wake panel input method. S.N.KIM 2018.
         IF(IREADW.EQ.1) ZH(NN,1)=RR*SIN(TH(NN,1)) !Wake panel input method. S.N.KIM 2018.
         FAC= (XH(NN,1)-XH(NHBUB+1,1)) / (XW(NHBDT,1)-XH(NHBUB+1,1))
         XH(NN,2)=XH(NHBUB+1,2)+ (XW(NHBDT,1)-XH(NHBUB+1,2)) * FAC
         TH(NN,2)=TH(NN,1)+DELT1
         YH(NN,2)=RR*COS(TH(NN,2))
         ZH(NN,2)=RR*SIN(TH(NN,2))
         XH(NN,MHBT1)=XH(NN,1)
         TH(NN,MHBT1)=TH(NN,1) + DELK
         YH(NN,MHBT1)=RR*COS(TH(NN,MHBT1))
         ZH(NN,MHBT1)=RR*SIN(TH(NN,MHBT1))
         DELT=TH(NN,MHBT1)-TH(NN,2)
         DO 80 M=3,MHBT
            TH(NN,M)=TH(NN,2) + DELT/(MHBT-1)*(M-2)
            XH(NN,M)=XH(NN,2) + (XH(NN,MHBT1)-XH(NN,2))
     *              *(TH(NN,M)-TH(NN,2))/DELT
            YH(NN,M)=RR*COS(TH(NN,M))
            ZH(NN,M)=RR*SIN(TH(NN,M))
80       CONTINUE
90    CONTINUE

C-----------------------------------------------------------------------
C     Hub at far upstream and far downstream of the blade is each 
C     replaced by a dipole disk
C-----------------------------------------------------------------------
      NHBDK=MHBT
      DTH=DELK/NHBDK
      DO 110 N=1,NHBDK+1

C.......Upstream disk
         XHDK(N,1)=XH(1,N)
         YHDK(N,1)=ZERO
         ZHDK(N,1)=ZERO
         XHDK(N,2)=XH(1,N)
         YHDK(N,2)=RHUB*COS(DTH*(N-1))
         ZHDK(N,2)=RHUB*SIN(DTH*(N-1))

C.......Downstream disk
         XHDKDT(N,1)=XH(NHBX+1,N)
         YHDKDT(N,1)=ZERO
         ZHDKDT(N,1)=ZERO
         XHDKDT(N,2)=XH(NHBX+1,N)
         YHDKDT(N,2)=RHUB*COS(DTH*(N-1))
         ZHDKDT(N,2)=RHUB*SIN(DTH*(N-1))

110   CONTINUE

      RETURN
C)))))))))))))))))))))) End of subroutine GHUB3(((((((((((((((((((((((((
      END
