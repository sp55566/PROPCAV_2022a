C***************************** MIT PSF10 ********************************
C
C                PROPELLER STEADY FLOW ANALYSIS PROGRAM
C                        PANEL METHOD SOLUTION         
C
C
C       Copyright (c) Massachusetts Institute of Technology 1997
C
C                     Release Date: 1 January 1997
C
C***********************************************************************
      SUBROUTINE GHUB4
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
 10   CONTINUE

      NHBDT=NSW(1)
      NHBU1=NHBU+1
      NHBUB=NHBU+NH
      NHBX=NHBUB+NHBDT
      MHBT1=MHBT+1

 20   CONTINUE

C-----------------------------------------------------------------------
C     Hub geometry in the near region of the blade
c       (blade section as a bi-sector)
C-----------------------------------------------------------------------
      PHI=ATAN( PITCH(1)/(PI*RHUB) )
      PSI1=PHI
      PSI2=ATAN2( XB(NHP+1,1)-XB(NHP,1),(THR(NHP+1)-THR(NHP))*RHUB )
      PSI0=HALF*(PI-PSI1-PSI2)
      DEL=RHUB*DELK*SIN(PHI+PSI0)/(MHBT*SIN(PHI))
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
 30   CONTINUE
      DO 40 N=1,NHP
         XH(NHBU+N,MHBT1)=XB(NHP+1-N,1)
         TH(NHBU+N,MHBT1)=THR(NHP+1-N) + DELK
         YH(NHBU+N,MHBT1)=RHUB*COS(TH(NHBU+N,MHBT1))
         ZH(NHBU+N,MHBT1)=RHUB*SIN(TH(NHBU+N,MHBT1))
 40   CONTINUE
C-----------------------------------------------------------------------
C     Hub geometry at upstream of the blade
C       --- half cosine spacing axially
C       --- constant spacing circumferentially
C-----------------------------------------------------------------------
      DTN=HALFPI/NHBU
      DO 50 N=1,NHBU
         XH(N,1)=XB(NHP,1) - XHBU*( 1.0-SIN(DTN*(N-1)) )
         TH(N,1)=THR(NHP)-PI/PITCH(1)*(XB(NHP,1)-XH(N,1))
       X1 = XH(N,1) - XH(1,1) 
       X2 = XHBU 
       RAT = X1/X2 
       RNOS = RNOSE2(RAT, RHUB)
C       write(*,*) 'RNOS =', n, rnos
       IF(RNOS .GT. RHUB) THEN
         RR = RHUB
       ELSEIF(RNOS .LE. 0.) THEN
         RR = 0.0
       ELSE
         RR = RNOS
       ENDIF
         YH(N,1)=RR*COS(TH(N,1))
         ZH(N,1)=RR*SIN(TH(N,1))
         XH(N,2)=XH(1,1) + (XH(NHBU1,2)-XH(1,1))* SIN( DTN*(N-1) )
         TH(N,2)=TH(NHBU1,2)-PI/PITCH(1)*(XH(NHBU1,2)-XH(N,2))
         YH(N,2)=RR*COS(TH(N,2))
         ZH(N,2)=RR*SIN(TH(N,2))
         XH(N,MHBT1)=XH(N,1)
         TH(N,MHBT1)=TH(N,1)+DELK
         YH(N,MHBT1)=RR*COS(TH(N,MHBT1))
         ZH(N,MHBT1)=RR*SIN(TH(N,MHBT1))
         DTH=TH(N,MHBT1)-TH(N,2)
         DO 70 M=3,MHBT
            TH(N,M)=TH(N,2)+ DTH/(MHBT-1)*(M-2)
            XH(N,M)=XH(N,2)+(XH(N,MHBT1)-XH(N,2))*(TH(N,M)-TH(N,2))/DTH
            YH(N,M)=RR*COS(TH(N,M))
            ZH(N,M)=RR*SIN(TH(N,M))
 70      CONTINUE
 50   CONTINUE

      DO N = NHBU+1, NHBUB
         DTH=TH(N,MHBT1)-TH(N,2)
         DO M=3,MHBT
            TH(N,M)=TH(N,2)+ DTH/(MHBT-1)*(M-2)
            XH(N,M)=XH(N,2)+(XH(N,MHBT1)-XH(N,2))*(TH(N,M)-TH(N,2))/DTH
            YH(N,M)=RHUB*COS(TH(N,M))
            ZH(N,M)=RHUB*SIN(TH(N,M))
       ENDDO
      ENDDO
      
C-----------------------------------------------------------------------
C     Hub geometry at downstream of the blade
C-----------------------------------------------------------------------
      DELT1=TH(NHBUB+1,2)-TH(NHBUB+1,1)
      DO 80 N=1,NHBDT+1
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
         DO 90 M=3,MHBT
            TH(NN,M)=TH(NN,2) + DELT/(MHBT-1)*(M-2)
            XH(NN,M)=XH(NN,2) + (XH(NN,MHBT1)-XH(NN,2))
     *              *(TH(NN,M)-TH(NN,2))/DELT
            YH(NN,M)=RR*COS(TH(NN,M))
            ZH(NN,M)=RR*SIN(TH(NN,M))
 90      CONTINUE
 80   CONTINUE

C-----------------------------------------------------------------------
C     Hub at far downstream of the blade is replaced by a dipole disk
C-----------------------------------------------------------------------
      NHBDK=MHBT
      DTH=DELK/NHBDK
      DO 110 N=1,NHBDK+1

C.......Downstream disk
         XHDKDT(N,1)=XH(NHBX+1,N)
         YHDKDT(N,1)=ZERO
         ZHDKDT(N,1)=ZERO
         XHDKDT(N,2)=XH(NHBX+1,N)
         YHDKDT(N,2)=RHUB*COS(DTH*(N-1))
         ZHDKDT(N,2)=RHUB*SIN(DTH*(N-1))

 110  CONTINUE

      RETURN
C)))))))))))))))))))))) End of subroutine GHUB4 ((((((((((((((((((((((((
      END






