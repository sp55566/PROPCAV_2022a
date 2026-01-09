        SUBROUTINE GEO2D_D

C*********************************************************************
C  This subroutine calculates the node coordinates (x,y) of the 2-D  *
C      blade strips used in the viscous boundary layer analysis.     *
C                                                                    *
C                                                                    *
C                            DUCT  CASE                              *
C                                                                    *       
C                  by Hong Sun          06/01/2006                   *
C*********************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C     XIVIS, ETAVIS : N=1,NC+1+NSW(M)  M=1,MR  
        
        DO 40 M = 1, MDUCTP
           DO 10 N = 1, NDUCTP
              XIVD(N,M) = XD(N,M) -XLED
              ETAVD(N,M) = SQRT(YD(N,M)**2+ZD(N,M)**2) - YLED
 10        CONTINUE
           DO 20 IN = 1, NWSUB+1
              XIVDWs(IN,M) = XWSD(IN,M) - XLED
              ETAVDWs(IN,M) = SQRT(YWSD(IN,M)**2+ZWSD(IN,M)**2) - YLED
 20        CONTINUE           
           DO 30 IN = 1, NDWK+1
              XIVDW(IN,M) = XDW(IN,M) - XLED
              ETAVDW(IN,M) = SQRT(YDW(IN,M)**2+ZDW(IN,M)**2) - YLED
 30        CONTINUE 
 40     CONTINUE 
 
C       Find the middle section on the blade       
      DO 80 M=1,MDUCT
           write(840,*) 'ZONE T= "M= ', M, ' " '

         DO 50 N=1,NDUCTP
            DX1 = XIVD(N,M+1)-XIVD(N,M)
            DY1 = ETAVD(N,M+1)-ETAVD(N,M)
            TH = ATAN2(DY1,DX1)
            RL = SQRT(DX1**2+DY1**2)
            XIVISD(N,M) = XIVD(N,M)+RL/2.*COS(TH)
              ETAVISD(N,M) = ETAVD(N,M)+RL/2.*SIN(TH)
                write(840,100) N,XIVISD(N,M),ETAVISD(N,M)
 50        CONTINUE

         DO 60 N=1,NWSUB
                N1 = N+1
            DX1 = XIVDWs(N1,M+1)-XIVDWs(N1,M)
            DY1 = ETAVDWs(N1,M+1)-ETAVDWs(N1,M)
            TH = ATAN2(DY1,DX1)
            RL = SQRT(DX1**2+DY1**2)
            XIVISD(N+NDUCTP,M) = XIVDWs(N1,M)+RL/2.*COS(TH)
              ETAVISD(N+NDUCTP,M) = ETAVDWs(N1,M)+RL/2.*SIN(TH)
                write(840,100) N+NDUCTP,XIVISD(N+NDUCTP,M),
     &                         ETAVISD(N+NDUCTP,M)
 60        CONTINUE

         DO 70 N=NSUB+1,NDWK
                N1 = N+1
            DX1 = XIVDW(N1,M+1)-XIVDW(N1,M)
            DY1 = ETAVDW(N1,M+1)-ETAVDW(N1,M)
            TH = ATAN2(DY1,DX1)
            RL = SQRT(DX1**2+DY1**2)
                N2 = N+NDUCTP+NWSUB-NSUB
            XIVISD(N2,M) = XIVDW(N1,M)+RL/2.*COS(TH)
              ETAVISD(N2,M) = ETAVDW(N1,M)+RL/2.*SIN(TH)
             write(840,100) N2,XIVISD(N2,M),ETAVISD(N2,M)
 70        CONTINUE
 80      CONTINUE

 100    format(1x,I5,2(3x,f10.5))
 101      format(4(f10.5,3x))
        close(840)

        RETURN 
        END
