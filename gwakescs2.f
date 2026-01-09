C ============================================================
      SUBROUTINE GWAKESCS2(KK)
C ============================================================
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      MN=MR+1
C/s S.N.KIM - I multiplied 3 to NDLTAT to split the first wake panel
C             into denser one. This is to avoid stability issue when
C             the cavity length falls around blade T.E. or the first wake
C             panel. Otherwise, stupid 'Ping Pong' loop may occur. 
      N1SUB=3*NDLTAT
      NSUB=NWMIN
      NWSUB1=1

      NWSUB=N1SUB+(NSUB-1)*NWSUB1

      DO 10 M=1,MN
         XWS(1,M)=XWW(1,M,KK)
         YWS(1,M)=YWW(1,M,KK)
         ZWS(1,M)=ZWW(1,M,KK)
         IF(KK.EQ.1) THEN
            XWS(1,M)=XB(1,M)
            YWS(1,M)=YB(1,M)
            ZWS(1,M)=ZB(1,M)             
         ENDIF
         DO 20 N=1,NSUB
            IF(N.EQ.1) THEN
               DWX=(XWW(2,M,KK)-XWW(1,M,KK))/FLOAT(N1SUB)
               DWY=(YWW(2,M,KK)-YWW(1,M,KK))/FLOAT(N1SUB)
               DWZ=(ZWW(2,M,KK)-ZWW(1,M,KK))/FLOAT(N1SUB)
               DO 30 L=1,N1SUB-1
                  XWS(L+1,M)=XWS(L,M)+DWX
                  YWS(L+1,M)=YWS(L,M)+DWY
                  ZWS(L+1,M)=ZWS(L,M)+DWZ
 30            CONTINUE
               XWS(N1SUB+1,M)=XWW(2,M,KK)
               YWS(N1SUB+1,M)=YWW(2,M,KK)
               ZWS(N1SUB+1,M)=ZWW(2,M,KK)
            ELSE
               DWX=(XWW(N+1,M,KK)-XWW(N,M,KK))/FLOAT(NWSUB1)
               DWY=(YWW(N+1,M,KK)-YWW(N,M,KK))/FLOAT(NWSUB1)
               DWZ=(ZWW(N+1,M,KK)-ZWW(N,M,KK))/FLOAT(NWSUB1)
               DO 40 L=1,NWSUB1
                  NIDX=N1SUB+(N-2)*NWSUB1+L+1
                  XWS(NIDX,M)=XWS(NIDX-1,M)+DWX
                  YWS(NIDX,M)=YWS(NIDX-1,M)+DWY
                  ZWS(NIDX,M)=ZWS(NIDX-1,M)+DWZ
 40            CONTINUE
            END IF
 20      CONTINUE
 10   CONTINUE
      
      RETURN
      END



