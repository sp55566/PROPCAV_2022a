       SUBROUTINE VCAVPLT(IPR)
************************************************************************
*                                                                      *
*  Subroutine Vertical CAVity PLoT plots out vertical cavities at each *
*  strip so that they may viewed with Tecplot                          *
*                           ----------        ----------               *
*                          |          |      |          |              *
*                          |  CAVOUT  |> --->|  VCAVPLT |              *
*                          |          |      |          |-->           *
*                           -----^----        ----------   |           *
*                                |                         |           *
*                                 --<---------<------------v           *
*                                                                      *
*  Author: Julie Young 080698                                          *
*  Date:      Revision/comments                                        *
*  -----      ---------------                                          *
*  JY080698   This routine currently will not plot supercavitating     *
*             panels.  However, it can handle face cavitation. This    *
*             routine is created because there was an index error in   *
*             the old VCAVPLT plus it did not plot the point where the *
*             panel splits.  All those problems are fixed in this new  *
*             version.                                                 *
*  JY090198   Modified routine so supercavity (non-split & split)      *
*             panels are also plotted.                                 *
*  JY021499   Modified subroutine to allow cavity to grow on both the  *
*             back and face of the foil.                               *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
*  JY060700   Modified subroutine to allow simultaneous face and back  *
*             supercavities.                                           *
*                                                                      *
************************************************************************
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       INCLUDE 'PUFCAVC.INC'
       DIMENSION CAVMRX(3,NBPZ,2),CAVMRXB(3,NBPZ,2)

       IF(IPR.EQ.1.AND.IDXREV.EQ.1) THEN
          CLOSE(700)
          OPEN(700,FILE='cav.plt',STATUS='UNKNOWN')
          WRITE(700,*) 'TITLE="Cavity planform : 3-D"' 
          WRITE(700,*) 'VARIABLES="x", "y", "z"'
       END IF

       ISR=1
       IF(IFACE.EQ.2) ISR=2

C-----------------------------------------------------------------------
C     Plot blade outline:
C     -----------------------------------------
C     NELMTS0 = no. of elements on the blade
C     NPTS0 = no. of nodes on the blade
C-----------------------------------------------------------------------
       NELMTS0=2*(NC+MR)
       NPTS0=2*(NCP+MRP)

C-----------------------------------------------------------------------
C     Plot the free surface for ISP=1
C-----------------------------------------------------------------------
       IF(ISP.EQ.1) THEN
          NELMTS0=NELMTS0+1
          NPTS0=NPTS0+4
       END IF

C-----------------------------------------------------------------------
C     NELMTS = no. of elements on the cavity
C-----------------------------------------------------------------------
       NFLAG=0
       NELMTS=0
       DO 70 M=1,MR
          IF(JCV(M,1)+JCV(M,2).GE.1) THEN
             NELMTS=NELMTS+JCV(M,1)+NSPP(M,1)+
     *            JCV(M,2)+NSPP(M,2)+NWC(M,1)+NWC(M,2)+
     *            NSPS(M,1)+NSPS(M,2)
             NFLAG=1
          END IF
 70    CONTINUE

C-----------------------------------------------------------------------
C     NPTS = no. of nodes on the cavity
C----------------------------------------------------------------------
       IF(IPR.EQ.1) THEN
          NPTS=0
          DO M=1,MR      

C...........Count no. of nodes on the blade       
             DO II=1,ISR
                IF((IFACE.EQ.0).OR.(II.EQ.1.AND.IFACE.EQ.2)) THEN      
                   IDR=1
                ELSE
                   IDR=2
                END IF
                JCAV=JCV(M,IDR)
                IF(JCAV.GE.1) THEN
                   DO N=1,JCAV+1
                      NPTS=NPTS+1
                   END DO
                   IF(NSPP(M,IDR).NE.0) NPTS=NPTS+1
                END IF
             END DO

C...........count no. of nodes on the supercavitating wake.
             IF(NNWC(M).GE.1) THEN
                DO N=2,NNWC(M)+NSPS(M,NWDIR(M))+1
                   NPTS=NPTS+1
                   IF(NWC(M,1).GT.0.AND.NWC(M,2).GT.0) NPTS=NPTS+1
                END DO
             END IF
          END DO

 5040     FORMAT(1X,'ZONE T="TIME= ',F4.0,'", N=',I8,', E=',I8,
     *         ', F=FEPOINT, ET=QUADRILATERAL')
 5060     FORMAT(3(1X,E16.8))
 5080     FORMAT(4(1X,I8))
          IF(NPTS.EQ.0) THEN
             WRITE(700,5040) TT(IDXREV),NPTS0,NELMTS0
          ELSE
             WRITE(700,5040) TT(IDXREV),NPTS0+2*NPTS,NELMTS0+NELMTS
          END IF
       END IF

       IF(IPR.EQ.1) THEN

          DT1=-TSTEP

C........plot blade root
          M=1
          DO N=1,NCP
             XX0 = XB(N,M)
             YY0 = YB(N,M)
             IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
             YY1=YY0*COS(DT1)-ZB(N,M)*SIN(DT1)
             ZZ1=YY0*SIN(DT1)+ZB(N,M)*COS(DT1)
             IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YY1)
             WRITE(700,5060) XX0,YY1,ZZ1
          END DO

C........plot blade TE
          N=1
          DO M=1,MRP
             XX0 = XB(N,M)
             YY0 = YB(N,M)
             IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
             YY1=YY0*COS(DT1)-ZB(N,M)*SIN(DT1)
             ZZ1=YY0*SIN(DT1)+ZB(N,M)*COS(DT1)
             IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YY1)
             WRITE(700,5060) XX0,YY1,ZZ1
          END DO

C........plot blade tip
          M=MRP
          DO N=1,NCP
             XX0 = XB(N,M)
             YY0 = YB(N,M)
             IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
             YY1=YY0*COS(DT1)-ZB(N,M)*SIN(DT1)
             ZZ1=YY0*SIN(DT1)+ZB(N,M)*COS(DT1)
             IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YY1)
             WRITE(700,5060) XX0,YY1,ZZ1
          END DO

C........plot blade LE
          N=NHP
          DO M=MRP,1,-1
             XX0 = XB(N,M)
             YY0 = YB(N,M)
             IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
             YY1=YY0*COS(DT1)-ZB(N,M)*SIN(DT1)
             ZZ1=YY0*SIN(DT1)+ZB(N,M)*COS(DT1)
             IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YY1)
             WRITE(700,5060) XX0,YY1,ZZ1
          END DO

C........plot free surface if ISP=1
          IF(ISP.EQ.1) THEN
             WRITE(700,5060) -1.,YFS,-1.2
             WRITE(700,5060)  1.,YFS,-1.2
             WRITE(700,5060)  1.,YFS, 1.2
             WRITE(700,5060) -1.,YFS, 1.2
          END IF

          DO 10 M=1,MR
             
             DO 20 II=1,ISR
                IF((IFACE.EQ.0).OR.(II.EQ.1.AND.IFACE.EQ.2)) THEN      
                   IDR=1
                   K=1
                   ISF=0
                ELSE
                   IDR=2
                   K=-1
                   ISF=1
                END IF
                
                M0M=M0(M,IDR)
                M0M1=M0M-1
                JCAV=JCV(M,IDR)

C-----------------------------------------------------------------------
C     Plot cavities on blade (LE -> TE)
C-----------------------------------------------------------------------
                IF(JCAV.GE.1) THEN
                   DO 30 N=1,JCAV+1
                      N1=M0M+K*(N-1)
                      DO 40 KK=1,3
                         CAVMRX(KK,N,IDR)=CONMRX(KK,N1,M)+
     *                        HT(N,M,IDR)*VECMRX(KK,N1,M)
                         CAVMRXB(KK,N,IDR)=CONMRX(KK,N1,M)
 40                   CONTINUE
 30                CONTINUE
                   
                   IF(NSPP(M,IDR).NE.0) THEN
                      N=JCAV+2
                      NN1=N1
                      NN2=N1+K*1
                      DXB=CONMRX(1,NN2,M)-CONMRX(1,NN1,M)
                      DYB=CONMRX(2,NN2,M)-CONMRX(2,NN1,M)
                      DZB=CONMRX(3,NN2,M)-CONMRX(3,NN1,M)
                      CAVMRX(1,N,IDR)=CONMRX(1,NN1,M)+DXB*FLP(M,IDR)
                      CAVMRX(2,N,IDR)=CONMRX(2,NN1,M)+DYB*FLP(M,IDR)
                      CAVMRX(3,N,IDR)=CONMRX(3,NN1,M)+DZB*FLP(M,IDR)
                      DO 50 KK=1,3
                         CAVMRXB(KK,N,IDR)=CAVMRX(KK,N,IDR)
 50                   CONTINUE
                   END IF
                END IF
                
 20          CONTINUE
             
             IF(NFLAG.EQ.1) THEN
                
C-----------------------------------------------------------------------
C      Plot cavity on blade
C-----------------------------------------------------------------------
                DO 90 II=1,ISR
                   IF((IFACE.EQ.0).OR.(II.EQ.1.AND.IFACE.EQ.2)) THEN   
                      IDR1=1
                   ELSE
                      IDR1=2
                   END IF
                   IF(JCV(M,IDR1).GE.1) THEN
                      DO 100 N=1,JCV(M,IDR1)+NSPP(M,IDR1)+1
                         XX0 = CAVMRX(1,N,IDR1)
                         YY0 = CAVMRX(2,N,IDR1)
                         ZZ0 = CAVMRX(3,N,IDR1)
                         IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
                         CALL ROTATE(YY0,ZZ0,YY1,ZZ1,DT1)
                         IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YY1)
                         WRITE(700,5060) XX0,YY1,ZZ1

                         XX0 = CAVMRXB(1,N,IDR1)
                         YY0 = CAVMRXB(2,N,IDR1)
                         ZZ0 = CAVMRXB(3,N,IDR1)
                         IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
                         CALL ROTATE(YY0,ZZ0,YY1,ZZ1,DT1)
                         IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YY1)
                         WRITE(700,5060) XX0,YY1,ZZ1
 100                  CONTINUE
                   END IF
                   
C-----------------------------------------------------------------------
C      Plot supercavity on wake
C-----------------------------------------------------------------------
                   IF(NNWC(M).GT.0) THEN
                      JW=(MR-M)*NTRA+1
                      IDR=NWDIR(M)
                      IDX=NHP-NLEP(M,IDXREV,IDR)                  
                      IF(NWC(M,1).GT.0.AND.NWC(M,2).GT.0) THEN
                         CALL CALHTW(1,M,HTETOP)
                         CALL CALHTW(2,M,HTEBOT)
                         HTE=HTETOP+HTEBOT
                         TP=HTETOP/HTE
                         TB=HTEBOT/HTE
                         WKCAM1=CONMRX(1,NCP,M)
                         WKCAM2=CONMRX(2,NCP,M)
                         WKCAM3=CONMRX(3,NCP,M)
                      ELSE IF(NWC(M,1).GT.0.AND.NWC(M,2).EQ.0) THEN
                         CALL CALHTW(1,M,HTE)
                         HTEHF=HTE/2.
                         WKCAM1=CONMRX(1,NCP,M)-HTEHF*VELW(JW,1)
                         WKCAM2=CONMRX(2,NCP,M)-HTEHF*VELW(JW,2)
                         WKCAM3=CONMRX(3,NCP,M)-HTEHF*VELW(JW,3)
                         TP=0.5
                         TB=0.5
                      ELSE IF(NWC(M,2).GT.0.AND.NWC(M,1).EQ.0) THEN
                         CALL CALHTW(2,M,HTE)
                         HTEHF=HTE/2.
                         WKCAM1=CONMRX(1,1,M)+HTEHF*VELW(JW,1)
                         WKCAM2=CONMRX(2,1,M)+HTEHF*VELW(JW,2)
                         WKCAM3=CONMRX(3,1,M)+HTEHF*VELW(JW,3)
                         TP=0.5
                         TB=0.5
                      END IF
                     
                      DO N=2,NNWC(M)+NSPS(M,NWDIR(M))+1
                         N1=(MR-M)*NTRA+N-1
                     
                         IF(N.LE.NNWC(M)+1) THEN
                        
                            WKCAM1=-DZW(N-1,M)*ULW(N1,1)+WKCAM1
                            WKCAM2=-DZW(N-1,M)*ULW(N1,2)+WKCAM2
                            WKCAM3=-DZW(N-1,M)*ULW(N1,3)+WKCAM3
                            
                            CAVTOP1=WKCAM1-
     *                           HT(IDX-1+N,M,IDR)*TP*VELW(N1,1)
                            CAVTOP2=WKCAM2-
     *                           HT(IDX-1+N,M,IDR)*TP*VELW(N1,2)
                            CAVTOP3=WKCAM3-
     *                           HT(IDX-1+N,M,IDR)*TP*VELW(N1,3)
                        
                            CAVBOT1=WKCAM1+
     *                           HT(IDX-1+N,M,IDR)*TB*VELW(N1,1)
                            CAVBOT2=WKCAM2+
     *                           HT(IDX-1+N,M,IDR)*TB*VELW(N1,2)
                            CAVBOT3=WKCAM3+
     *                           HT(IDX-1+N,M,IDR)*TB*VELW(N1,3)
                        
                         ELSE 

                            WKCAM1=-DZW(N-1,M)*FLS(M,IDR)*ULW(N1+1,1)
     *                           +WKCAM1
                            WKCAM2=-DZW(N-1,M)*FLS(M,IDR)*ULW(N1+1,2)
     *                           +WKCAM2
                            WKCAM3=-DZW(N-1,M)*FLS(M,IDR)*ULW(N1+1,3)
     *                           +WKCAM3
                     
                            CAVTOP1=WKCAM1
                            CAVTOP2=WKCAM2
                            CAVTOP3=WKCAM3

                            CAVBOT1=WKCAM1
                            CAVBOT2=WKCAM2
                            CAVBOT3=WKCAM3
                            
                         END IF

                         IF(NWC(M,1).GT.0.AND.NWC(M,2).GT.0) THEN
                            IF(IDR1.EQ.1) THEN
                               IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,
     %                              CAVTOP1,CAVTOP2)
                               CALL ROTATE(CAVTOP2,CAVTOP3,YY1,ZZ1,DT1)
                               IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,
     %                              CAVTOP1,YY1)
                               WRITE(700,5060) CAVTOP1,YY1,ZZ1
                               IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,
     %                              WKCAM1,WKCAM2)
                               CALL ROTATE(WKCAM2,WKCAM3,YY1,ZZ1,DT1)
                               IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,
     %                              WKCAM1,YY1)
                               WRITE(700,5060) WKCAM1,YY1,ZZ1
                            ELSE IF(IDR1.EQ.2) THEN
                               IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,
     %                              CAVBOT1,CAVBOT2)
                               CALL ROTATE(CAVBOT2,CAVBOT3,YY1,ZZ1,DT1)
                               IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,
     %                              CAVBOT1,YY1)
                               WRITE(700,5060) CAVBOT1,YY1,ZZ1
                               IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,
     %                              WKCAM1,WKCAM2)
                               CALL ROTATE(WKCAM2,WKCAM3,YY1,ZZ1,DT1)
                               IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,
     %                              WKCAM1,YY1)
                               WRITE(700,5060) WKCAM1,YY1,ZZ1
                            END IF
                         ELSE IF(NWC(M,1).GT.0.AND.NWC(M,2).EQ.0.
     *                           AND.IDR1.EQ.1) THEN
                            
                            IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,
     %                           CAVTOP1,CAVTOP2)
                            CALL ROTATE(CAVTOP2,CAVTOP3,YY1,ZZ1,DT1)
                            IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,
     %                           CAVTOP1,YY1)
                            WRITE(700,5060) CAVTOP1,YY1,ZZ1
                            IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,
     %                           CAVBOT1,CAVBOT2)
                            CALL ROTATE(CAVBOT2,CAVBOT3,YY1,ZZ1,DT1)
                            IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,
     %                           CAVBOT1,YY1)
                            WRITE(700,5060) CAVBOT1,YY1,ZZ1
                         ELSE IF(NWC(M,2).GT.0.AND.NWC(M,1).EQ.0.
     *                           AND.IDR1.EQ.2) THEN
                            IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,
     %                           CAVBOT1,CAVBOT2)
                            CALL ROTATE(CAVBOT2,CAVBOT3,YY1,ZZ1,DT1)
                            IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,
     %                           CAVBOT1,YY1)
                            WRITE(700,5060) CAVBOT1,YY1,ZZ1
                            IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,
     %                           CAVTOP1,CAVTOP2)
                            CALL ROTATE(CAVTOP2,CAVTOP3,YY1,ZZ1,DT1)
                            IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,
     %                           CAVTOP1,YY1)
                            WRITE(700,5060) CAVTOP1,YY1,ZZ1
                         END IF
                      END DO

                   END IF

 90             CONTINUE
             END IF
 10       CONTINUE

C-----------------------------------------------------------------------
C     Establishing connectivity
C-----------------------------------------------------------------------
C........plot blade root
          DO N=1,NC
             WRITE(700,5080) N, N+1, N+1, N
          END DO

C........plot blade TE
          IMM = NCP
          DO M=1,MR
             WRITE(700,5080) IMM+M,IMM+M+1,IMM+M+1,IMM+M
          END DO

C........plot blade tip
          IMM = IMM + MRP
          DO N=1,NC
             WRITE(700,5080) IMM+N, IMM+N+1, IMM+N+1, IMM+N
          END DO

C........plot blade LE
          IMM = IMM+NCP
          DO M=1,MR
             WRITE(700,5080) IMM+M,IMM+M+1,IMM+M+1,IMM+M
          END DO

C........plot free surface
          IF(ISP.EQ.1) THEN
             IMM=IMM+MRP
             WRITE(700,5080) IMM+1,IMM+2,IMM+3,IMM+4
          END IF
          
          IF(NFLAG.EQ.1) THEN
             ICNT=NPTS0+1
             DO 130 M=1,MR
                NCAV1=JCV(M,1)+NSPP(M,1)
                NCAV2=JCV(M,2)+NSPP(M,2)  
                IF(NNWC(M).GT.0) THEN
                   NCAVW=NNWC(M)+NSPS(M,NWDIR(M))
                ELSE
                   NCAVW=0
                END IF

                IF(NWC(M,1).GT.0) NCAV1=NCAV1+NCAVW
                IF(NWC(M,2).GT.0) NCAV2=NCAV2+NCAVW

                IF((NCAV1+NCAV2).GE.1) THEN
                   DO 140 N=1,NCAV1
                      WRITE(700,5080) ICNT,ICNT+1,ICNT+3,ICNT+2
                      ICNT=ICNT+2
 140               CONTINUE
                   
                   IF(NCAV2.GT.0.AND.NCAV1.GT.0) ICNT=ICNT+2
                   DO 150 N=1,NCAV2
                      WRITE(700,5080) ICNT,ICNT+1,ICNT+3,ICNT+2
                      ICNT=ICNT+2
 150               CONTINUE

                   ICNT=ICNT+2
                END IF
 130         CONTINUE
             
          END IF
       END IF

C-----------------------------------------------------------------------
C     If IFILE=1, close file 700 at the end of every revolution.
C-----------------------------------------------------------------------
       IF(IFILE.EQ.1.AND.IDXREV.EQ.NTPREV) CLOSE(700)
       IF(IFILE.EQ.1.AND.ISTEADY.EQ.0) CLOSE(700)

       IF(IPR.EQ.1.AND.IDXREV.EQ.NTPREV) CLOSE(700)

       RETURN
       END

                

       SUBROUTINE ROTATE(YIN,ZIN,YB1,ZB1,DT1)

       YB1=YIN*COS(DT1)-ZIN*SIN(DT1)
       ZB1=YIN*SIN(DT1)+ZIN*COS(DT1)

       RETURN
       END
