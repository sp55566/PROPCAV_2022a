      SUBROUTINE GWAKESCS
************************************************************************
*     GWAKES: Geometry of the WAKE Subpanels                           *
*     ---  Generate geometries of the wake subpanels                   *
*          wake                                                        *
*                                                                      *
*     Date     Comment or Revision                                     *
*     -------- -------------------                                     *
*     JY081499 This file was modified to generate subpanels for        *
*              supercavity panels.                                     *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'


      MN=MRP

C-----------------------------------------------------------------------
C     Generate the wake subpanels
C     Hong Sun modified for the first 4  wake panels
C-----------------------------------------------------------------------

C.....N1SUB is the number of subpanels the first wake panel will be 
C.....devided into.  
      N1SUB=INT(DTPROP)

C.....NSUB is the total number of macro panels that will be subdivided.
C.....Note: If ISP=1, NSUB must equal to NWMINFW. (JY112499)
      NSUB=NWMINFW
      
C.....NWSUB1 is the number of panels each macro panel (except for the
C.....first panel) will be divided into.  This should always be 1.
      NWSUB1=1

C.....NWSUB is the total number of subpanels.
      IF(IVISC.EQ.1) THEN
C      IF(IVISC.EQ.1.OR.IREADW.EQ.1) THEN
         NWSUB=N1SUB*4+(NSUB-4)*NWSUB1    
      ELSE  
         NWSUB=N1SUB+(NSUB-1)*NWSUB1
      END IF
 

C....Next IF statement added by JY072200
      IF(NWSUB.GE.NWZ) THEN
         WRITE(*,'('' M = '',I5,'' NWSUB = '',I5,'' NWZ = '',I5)') 
     *        M,NWSUB,NWZ
         WRITE(*,'('' Decrease XUWDK in .geo, OR '')')
         WRITE(*,'('' Increase XWZ --> (NWSUB+1) in PARAM.INC! '')')
         STOP
      END IF

C-----------------------------------------------------------------------
C     Generate the transition wake geometry
C-----------------------------------------------------------------------
      DO 10 M=1,MN

c........Initial coordinates of the transition wake.....................
         XWS(1,M)=XB(1,M) !XW(1,M)
         YWS(1,M)=YB(1,M) !YW(1,M)
         ZWS(1,M)=ZB(1,M) !ZW(1,M)

         DO 20 N=1,NSUB
 
            IF(IVISC.EQ.1)THEN
C            IF(IVISC.EQ.1.OR.IREADW.EQ.1)THEN
              IF(N.LE.4) THEN

                 DWX=(XW(N+1,M)-XW(N,M))/FLOAT(N1SUB)
                 DWY=(YW(N+1,M)-YW(N,M))/FLOAT(N1SUB)
                 DWZ=(ZW(N+1,M)-ZW(N,M))/FLOAT(N1SUB)
                 DO 35 L=1,N1SUB
                    NIDX=(N-1)*N1SUB+L+1
                    XWS(NIDX,M)=XWS(NIDX-1,M)+DWX
                    YWS(NIDX,M)=YWS(NIDX-1,M)+DWY
                    ZWS(NIDX,M)=ZWS(NIDX-1,M)+DWZ
 35              CONTINUE
c
                 XWS(4*N1SUB+1,M)=XW(5,M)
                 YWS(4*N1SUB+1,M)=YW(5,M)
                 ZWS(4*N1SUB+1,M)=ZW(5,M)

              ELSE
                 DWX=(XW(N+1,M)-XW(N,M))/FLOAT(NWSUB1)
                 DWY=(YW(N+1,M)-YW(N,M))/FLOAT(NWSUB1)
                 DWZ=(ZW(N+1,M)-ZW(N,M))/FLOAT(NWSUB1)
                 DO 45 L=1,NWSUB1
                    NIDX=N1SUB*4+(N-5)*NWSUB1+L+1
                    XWS(NIDX,M)=XWS(NIDX-1,M)+DWX
                    YWS(NIDX,M)=YWS(NIDX-1,M)+DWY
                    ZWS(NIDX,M)=ZWS(NIDX-1,M)+DWZ
 45              CONTINUE
              END IF    !N
          ELSE                            
              IF(N.EQ.1) THEN

                 DWX=(XW(2,M)-XW(1,M))/FLOAT(N1SUB)
                 DWY=(YW(2,M)-YW(1,M))/FLOAT(N1SUB)
                 DWZ=(ZW(2,M)-ZW(1,M))/FLOAT(N1SUB)
                 DO 30 L=1,N1SUB-1
                    XWS(L+1,M)=XWS(L,M)+DWX
                    YWS(L+1,M)=YWS(L,M)+DWY
                    ZWS(L+1,M)=ZWS(L,M)+DWZ
 30              CONTINUE

                 XWS(N1SUB+1,M)=XW(2,M)
                 YWS(N1SUB+1,M)=YW(2,M)
                 ZWS(N1SUB+1,M)=ZW(2,M)

              ELSE
                 DWX=(XW(N+1,M)-XW(N,M))/FLOAT(NWSUB1)
                 DWY=(YW(N+1,M)-YW(N,M))/FLOAT(NWSUB1)
                 DWZ=(ZW(N+1,M)-ZW(N,M))/FLOAT(NWSUB1)
                 DO 40 L=1,NWSUB1
                    NIDX=N1SUB+(N-2)*NWSUB1+L+1
                    XWS(NIDX,M)=XWS(NIDX-1,M)+DWX
                    YWS(NIDX,M)=YWS(NIDX-1,M)+DWY
                    ZWS(NIDX,M)=ZWS(NIDX-1,M)+DWZ
 40              CONTINUE
              END IF
            END IF       !IVISC               
     
 20      CONTINUE
 10   CONTINUE
      
      RETURN
C))))))))))))))))))))) End of subroutine GWAKE (((((((((((((((((((((((((
      END



