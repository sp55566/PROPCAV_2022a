C------------------------------------------------------------------
      SUBROUTINE XFLAPPT
C------------------------------------------------------------------
C
C     X-FLAP TO BE DETERMINED

      use m_INPGEO2
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'INPGEO.INC'
CSH----getting XFLAP------------------------------------------------
            XFLAP=CHFLAP*(-XB((NCP+1)/2,1)+XB(1,1))+XB((NCP+1)/2,1)
            WRITE(*,*)
            WRITE(*,*) XFLAP,"X CO-ORD OF FLAP PIVOT"
CSH---------------------------------------------------------------   

      RETURN
      END

C------------------------------------------------------------------
      SUBROUTINE GFLAP
C------------------------------------------------------------------
C
C     FLAP TO BE ROTATED ABOUT FLAP PIVOT

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'INPGEO.INC'
!     PARAMETER(NCOFBT=4*NBZ)
      DIMENSION XFB(NCP,2),XFT(NCP,2)
      DIMENSION XIN(NCP),ZIN(NCP),XOUT(NCP),ZOUT(NCP)
      DIMENSION XCOFB(4*NBZ),ZCOFB(4*NBZ)
      DIMENSION XCOFT(4*NBZ),ZCOFT(4*NBZ)
      DIMENSION SFB1(NCP),SFB2(NCP),SFT1(NCP),SFT2(NCP)
      DIMENSION XBF1(NCP,MRP),ZBF1(NCP,MRP),SPFL(NCP)
CSH----getting AND RORATING FLAP CORD-----------------------------------
      integer NCOFBT
      NCOFBT = 4*NBZ

      COSFL=COS(ANGFLAP*PI/180.)
      SINFL=SIN(ANGFLAP*PI/180.)

      DO M=1,MRP

C     BOTTOM
         IFB=0
      DO I=1,(NCP+1)/2
         IF(XB(I,M).GE.XFLAP) THEN
            XFMP=XB(I,M)-XFLAP
            ZFMP=ZB(I,M)
            XFB(I,1)=XFMP*COSFL-ZFMP*SINFL+XFLAP
            XFB(I,2)=XFMP*SINFL+ZFMP*COSFL
            IFB=IFB+1
C            WRITE(*,*)XFB(I,1),XFB(I,2),I
         ELSE
            XFB(I,1)=XB(I,M)
            XFB(I,2)=ZB(I,M)
         ENDIF  
C     TOP
      ENDDO  
        IFT=0
      DO I=(NCP+1)/2,NCP
         I1=I-(NCP+1)/2+1
         IF(XB(I,M).GE.XFLAP) THEN
            XFMP=XB(I,M)-XFLAP
            ZFMP=ZB(I,M)
            XFT(I1,1)=XFMP*COSFL-ZFMP*SINFL+XFLAP
            XFT(I1,2)=XFMP*SINFL+ZFMP*COSFL
C            WRITE(*,*)XFT(I1,1),XFT(I1,2),I1
         ELSE
            IFT=IFT+1
            XFT(I1,1)=XB(I,M)
            XFT(I1,2)=ZB(I,M)            
          ENDIF  
          
      ENDDO  
CSH----arranging  FLAP CORD----------------------------------------------
C      WRITE(*,*)IFT,IFB
C     BOTTOM
       NRB=0
      DO I=1,IFB
         IF(XFB(I,1).GE.XFLAP)THEN
            NRB=NRB+1
         ENDIF  
      ENDDO   

      DO I=IFB+1,(NCP+1)/2
         NRB=NRB+1
         XFB(NRB,1)=XFB(I,1)
         XFB(NRB,2)=XFB(I,2)
      ENDDO   
C      WRITE(*,*)NRB,(NCP+1)/2
C      WRITE(*,*)
C     TOP
       NRT=0
      DO I=1,IFT
         NRT=NRT+1
         XFT(NRT,1)=XFT(I,1)
         XFT(NRT,2)=XFT(I,2)
      ENDDO  

      DO I=IFT+1,(NCP+1)/2
         IF(XFT(I,1).GE.XFLAP)THEN
            NRT=NRT+1
         XFT(NRT,1)=XFT(I,1)
         XFT(NRT,2)=XFT(I,2)
         ENDIF   
      ENDDO   
C      WRITE(*,*)NRB,NRT

CSH----interpolating top and bot  FLAP CORD---------------------------------
C     BOT
        SFB1(1)=0.0
        DO N=2,NRB
           DSX=XFB(N,1)-XFB(N-1,1)
           DSZ=XFB(N,2)-XFB(N-1,2)
           DSXZ=SQRT(DSX*DSX+DSZ*DSZ)
           SFB1(N)=SFB1(N-1)+DSXZ
        ENDDO   
C     UNIFORM
C        SFB2(1)=0.0
C        DSFB=SFB1(NRB)/REAL((NCP-1)/2)
C        DO N=2,(NCP+1)/2
C           SFB2(N)=SFB2(N-1)+DSFB
C        ENDDO  

C     FULL COSINE
        DO N=1,(NCP+1)/2
           SPFL(N)=0.5*(1.-COS(PI*REAL(N)/REAL((NCP-1)/2+2)))
        ENDDO   

        DO N=1,(NCP+1)/2
           SFB2(N)=SPFL(N)*SFB1(NRB)
        ENDDO  

        DO N=1,NRB
           XIN(N)=XFB(N,1)
           ZIN(N)=XFB(N,2)
        ENDDO   

        CALL UGLYDK(NRB,1,1,SFB1,XIN,0,0,XCOFB)
        CALL UGLYDK(NRB,1,1,SFB1,ZIN,0,0,ZCOFB)

        CALL EVALDK(NRB,(NCP+1)/2,SFB1,SFB2,XOUT,XCOFB)
        CALL EVALDK(NRB,(NCP+1)/2,SFB1,SFB2,ZOUT,ZCOFB)

        DO N=1,(NCP+1)/2
           XBF1(N,M)=XOUT(N)
           ZBF1(N,M)=ZOUT(N)
        ENDDO
   
C     TOT
        SFT1(1)=0.0
        DO N=2,NRT
           DSX=XFT(N,1)-XFT(N-1,1)
           DSZ=XFT(N,2)-XFT(N-1,2)
           DSXZ=SQRT(DSX*DSX+DSZ*DSZ)
           SFT1(N)=SFT1(N-1)+DSXZ
        ENDDO  

C     UNIFORM
C        SFT2(1)=0.0
C        DSFT=SFT1(NRT)/REAL((NCP-1)/2)
C        DO N=2,(NCP+1)/2
C           SFT2(N)=SFT2(N-1)+DSFT
C        ENDDO  

C     FULL COSINE
        DO N=1,(NCP+1)/2
           SFT2(N)=SPFL(N)*SFT1(NRT)
        ENDDO  

        DO N=1,NRT
           XIN(N)=XFT(N,1)
           ZIN(N)=XFT(N,2)
        ENDDO   

        CALL UGLYDK(NRT,1,1,SFT1,XIN,0,0,XCOFT)
        CALL UGLYDK(NRT,1,1,SFT1,ZIN,0,0,ZCOFT)

        CALL EVALDK(NRT,(NCP+1)/2,SFT1,SFT2,XOUT,XCOFT)
        CALL EVALDK(NRT,(NCP+1)/2,SFT1,SFT2,ZOUT,ZCOFT)

        DO N=1,(NCP+1)/2-1
           N1=(NCP+1)/2+N
           XBF1(N1,M)=XOUT(N+1)
           ZBF1(N1,M)=ZOUT(N+1)
        ENDDO   
CSH------------------------------------------------------------------   
      ENDDO
CSH----writing  FOIL-FLAP CORD------------------------------------------
        OPEN(487,FILE='foilflap.plt',STATUS='UNKNOWN')
        WRITE(487,*)'ZONE T="FLAP"','I=',NCP,'J=',MRP
      DO M=1,MRP
         DO N=1,NCP
            WRITE(487,*)XBF1(N,M),YB(N,M),ZBF1(N,M)
            XB(N,M)=XBF1(N,M)
            ZB(N,M)=ZBF1(N,M)
         ENDDO
      ENDDO   
      CLOSE(487)
CSH------------------------------------------------------------------   

 
CSH------------------------------------------------------------------   
C      STOP
      RETURN
      END
