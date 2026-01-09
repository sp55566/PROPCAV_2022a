        SUBROUTINE GEO2DC

C*********************************************************************
C  This subroutine calculates the node coordinates (x,y) of the 2-D  *
C      blade strips used in the viscous boundary layer analysis.     *
C                                                                    *
C                   Created by Hong Sun, 2005                        *
C                                                                    *
C                       CAVITATING CASE                              * 
C                                                                    *  
C*********************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C     XIVIS, ETAVIS : N=1,NC+1+NSW(M)  M=1,MR  

        DO 30 M = 1, MRP
         
         IF(RZ(M).EQ.0) THEN
           RZ(M)=1.0
         ENDIF

       TANP = PITCH(M)/PI/RZ(M)
       COSP = 1.0/SQRT(1.+TANP**2)
         SINP = TANP*COSP

C         write(83,*) 'ZONE T= "M= ', M, ' " '
         NN = 0 
         DO 20 N = 2, NSW(M)
         TH = ATAN2(ZW(N,M),YW(N,M))
         THM1 = ATAN2(ZW(N-1,M),YW(N-1,M))

c          make sure TH is continuous
         IF(TH.LT.0.AND.THM1.GT.0) NN = NN+1
         IF(TH.GT.0.AND.THM1.LT.0) NN = NN+1
         IF(TH.LT.0.) THEN
             TH = TH + 2.*PI*NN
           ELSE IF(TH.GE.0.AND.NN.GT.0) THEN
             TH = TH + 2.*PI*(NN-1)
           ENDIF

         R0 = SQRT(YW(N,M)**2+ZW(N,M)**2)
         X1 = RAKE(M)*2.
         TH1 = SKEW(M)*RAD
C           write(83,*) N, TH, NN, THM1

         XI1 = (XW(N,M)-X1)*SINP + COSP*(R0*(TH-TH1))
         ETA1 = -(XW(N,M)-X1)*COSP + SINP*(R0*(TH-TH1)) 
  
         XIVW(N-1,M) = SQRT(XI1**2+ETA1**2)           
           IF(ICON.EQ.5)  XIVW(N-1,M) = XW(N,M)-X1
         ETAVW(N-1,M) = ETAV(NCP,M)
 
 20     CONTINUE

         DO N = 1, NSUB
           IF(N.EQ.1) THEN
             DWXI = (XIVW(N,M)-XIV(1,M))/FLOAT(N1SUB)
             DWETA = (ETAVW(N,M)-ETAV(1,M))/FLOAT(N1SUB)
             XIVWs(1,M) = XIV(1,M) + DWXI
           ETAVWs(1,M) = ETAV(1,M) + DWETA 
             DO L = 2, N1SUB
                NIDX=(N-1)*N1SUB+L
              XIVWs(NIDX,M) = XIVWs(NIDX-1,M) + DWXI
              ETAVWs(NIDX,M) = ETAVWs(NIDX-1,M) + DWETA 
             ENDDO
           ELSE IF(N.LE.4) THEN
             DWXI = (XIVW(N,M)-XIVW(N-1,M))/FLOAT(N1SUB)
             DWETA = (ETAVW(N,M)-ETAVW(N-1,M))/FLOAT(N1SUB)
             DO L = 1, N1SUB
                NIDX=(N-1)*N1SUB+L
              XIVWs(NIDX,M) = XIVWs(NIDX-1,M) + DWXI
              ETAVWs(NIDX,M) = ETAVWs(NIDX-1,M) + DWETA 
             ENDDO
           ELSE
             DWXI = (XIVW(N,M)-XIVW(N-1,M))/FLOAT(NWSUB1)
             DWETA = (ETAVW(N,M)-ETAVW(N-1,M))/FLOAT(NWSUB1)
             DO L = 1, NWSUB1
                NIDX=4*N1SUB+(N-5)*NWSUB1+L
              XIVWs(NIDX,M) = XIVWs(NIDX-1,M) + DWXI
              ETAVWs(NIDX,M) = ETAVWs(NIDX-1,M) + DWETA 
             ENDDO
           END IF
        ENDDO     

 30    CONTINUE



C       Find the middle section on the blade       
      DO 40 M=1,MR
           write(840,*) 'ZONE T= "M= ', M, ' " '

         DO 5 N=1,NCP
            DX1 = XIV(N,M+1)-XIV(N,M)
            DY1 = ETAV(N,M+1)-ETAV(N,M)
            TH = ATAN2(DY1,DX1)
            RL = SQRT(DX1**2+DY1**2)
            XIVIS(N,M) = XIV(N,M)+RL/2.*COS(TH)
              ETAVIS(N,M) = ETAV(N,M)+RL/2.*SIN(TH)
             write(840,100) N,XIVIS(N,M),ETAVIS(N,M)
 5         CONTINUE

         DO 15 N=1,NTRA
            DX1 = XIVWs(N,M+1)-XIVWs(N,M)
            DY1 = ETAVWs(N,M+1)-ETAVWs(N,M)
            TH = ATAN2(DY1,DX1)
            RL = SQRT(DX1**2+DY1**2)
            XIVIS(N+NCP,M) = XIVWs(N,M)+RL/2.*COS(TH)
              ETAVIS(N+NCP,M) = ETAVWs(N,M)+RL/2.*SIN(TH)
             write(840,100) N+NCP,XIVIS(N+NCP,M),ETAVIS(N+NCP,M)
 15         CONTINUE

         DO 25 N=NSUB+1,NSW(M)-2
            DX1 = XIVW(N,M+1)-XIVW(N,M)
            DY1 = ETAVW(N,M+1)-ETAVW(N,M)
            TH = ATAN2(DY1,DX1)
            RL = SQRT(DX1**2+DY1**2)
                N1 = N+NCP+NTRA-NSUB
            XIVIS(N1,M) = XIVW(N,M)+RL/2.*COS(TH)
              ETAVIS(N1,M) = ETAVW(N,M)+RL/2.*SIN(TH)
             write(840,100) N1,XIVIS(N1,M),ETAVIS(N1,M)
 25         CONTINUE


40      CONTINUE

 100    format(1x,I5,2(3x,f10.5))
 101      format(4(f10.5,3x))

        RETURN 
        END
