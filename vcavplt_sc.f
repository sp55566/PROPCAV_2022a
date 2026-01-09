       SUBROUTINE VCAVPLT_SC(IPR)
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
*  Author: Julie Young 112001                                          *
*                                                                      *
*  Date:      Revision/comments                                        *
*  -----      ---------------                                          *
************************************************************************
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       INCLUDE 'PUFCAVC.INC'
       DIMENSION CAVMRX(3,NBPZ,2),CAVMRXB(3,NBPZ,2)

       IF(IPR.EQ.1.AND.IDXREV.EQ.1) THEN
          CLOSE(700)
          OPEN(700,FILE='cav.plt',STATUS='UNKNOWN')
          WRITE(700,*) 'TITLE="Cavity Planform : 3-D"'
          WRITE(700,*) 'VARIABLES="x", "y", "z"'
       END IF

       NCOLDP=NCOLD+1

C-----------------------------------------------------------------------
C     Plot blade outline:
C     -----------------------------------------
C     NELMTS0 = no. of elements on the blade
C     NPTS0 = no. of nodes on the blade
C-----------------------------------------------------------------------
       NELMTS0=2*(NCOLD+MR)
       NPTS0=2*(NCOLDP+MRP)+MRP
       
C-----------------------------------------------------------------------
C     Plot the free surface for ISP=1
C-----------------------------------------------------------------------
       IF(ISP.EQ.1) THEN
          NELMTS0=NELMTS0+1
          NPTS0=NPTS0+4
       END IF

C-----------------------------------------------------------------------
C     NELMTS = no. of elements on the cavity
C     NPTSS = no. of nodes on the cavity
C-----------------------------------------------------------------------
       NELMTS=0
       NPTS=0

       IF(ISP.EQ.1.AND.NBW.EQ.0) GO TO 1000

       DO 70 M=1,MR
          DO IDR=1,2
             IF(ISP.EQ.0) THEN
                IF(JCV(M,IDR).GT.0) THEN
                   N1=JCV(M,IDR)
                   IF(N1.GT.0) THEN
                      NELMTS=NELMTS+N1
                      NPTS=NPTS+(N1+1)*2 
                   END IF
                END IF
             ELSE
                IF(ISUB1(M,IDR,IDXREV).GT.0) THEN
                   IF(IDR.EQ.1) THEN
                      NI1=JCV(M,1)
                      IF(NI1.GT.0) THEN
                         NELMTS=NELMTS+NI1
                         NPTS=NPTS+(NI1+1)*2
                      END IF
                   END IF
                END IF
             END IF
          END DO

          IF(NNWC(M).GT.0) THEN          
             NDUM=NNWC(M)+NSPS(M,NWDIR(M))+1
             NELMTS=NELMTS+NDUM
             NPTS=NPTS+(NDUM+1)*2
          END IF
 70    CONTINUE

 1000  CONTINUE

 5040  FORMAT(1X,'ZONE T="TIME= ',F4.0,'", N=',I8,', E=',I8,
     *      ', F=FEPOINT, ET=QUADRILATERAL')
 5060  FORMAT(3(1X,E16.8))
 5080  FORMAT(4(1X,I8))
       IF(NPTS.EQ.0) THEN
          WRITE(700,5040) TT(IDXREV),NPTS0,NELMTS0
       ELSE
          WRITE(700,5040) TT(IDXREV),NPTS0+NPTS,NELMTS0+NELMTS
       END IF

       IF(IPR.EQ.1) THEN

          DT1=-TSTEP

C........plot blade root
          M=1
          DO N=N0(2),N0(1)
             XX0 = XB(N,M)
             YY0 = YB(N,M)
             IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
             YY1=YY0*COS(DT1)-ZB(N,M)*SIN(DT1)
             ZZ1=YY0*SIN(DT1)+ZB(N,M)*COS(DT1)
             IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YY1)
             WRITE(700,5060) XX0,YY1,ZZ1
          END DO

C........plot blade TE
          N=N0(1)
          DO M=1,MRP
             XX0 = XB(N,M)
             YY0 = YB(N,M)
             IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
             YY1=YY0*COS(DT1)-ZB(N,M)*SIN(DT1)
             ZZ1=YY0*SIN(DT1)+ZB(N,M)*COS(DT1)
             IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YY1)
             WRITE(700,5060) XX0,YY1,ZZ1
          END DO
          N=N0(2)
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
          DO N=N0(2),N0(1)
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

          IF(ISP.EQ.1.AND.NBW.EQ.0) GO TO 2000
C-----------------------------------------------------------------------
C     Plot cavities 
C-----------------------------------------------------------------------
          DO M=1,MR
             DO IDR=1,2
                IF(ISP.EQ.0) THEN
                   IF(JCV(M,IDR).EQ.0) THEN
                      GO TO 1500
                   ELSE
                      NSTART=NLEP(M,IDXREV,IDR)+1
                      NEND=NSTART+JCV(M,IDR)
                   END IF
                ELSE
                   IF(ISUB1(M,IDR,IDXREV).EQ.0) GO TO 1500
                   IF(JCV(M,IDR).EQ.0) GO TO 1500
                   IF(IDR.EQ.1) THEN
                      NSTART=NLEP(M,IDXREV,IDR)+1
                      NEND=NSTART+JCV(M,IDR)
                   END IF
                END IF

                DO N=NSTART,NEND
                   IF(IDR.EQ.1) THEN
                      N1=NH+N
                   ELSE
                      N1=NHP+1-N
                   END IF

                   DO KK=1,3
                      CAVMRX(KK,N,IDR)=CONMRX(KK,N1,M)+
     *                     HTP(N,M,IDR)*VECMRX(KK,N1,M)
                      CAVMRXB(KK,N,IDR)=CONMRX(KK,N1,M)
                   END DO

                   XX0 = CAVMRX(1,N,IDR)
                   YY0 = CAVMRX(2,N,IDR)
                   ZZ0 = CAVMRX(3,N,IDR)
                   IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
                   CALL ROTATE(YY0,ZZ0,YY1,ZZ1,DT1)
                   IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YY1)
                   WRITE(700,5060) XX0,YY1,ZZ1

                   XX0 = CAVMRXB(1,N,IDR)
                   YY0 = CAVMRXB(2,N,IDR)
                   ZZ0 = CAVMRXB(3,N,IDR)
                   IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,XX0,YY0)
                   CALL ROTATE(YY0,ZZ0,YY1,ZZ1,DT1)
                   IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,XX0,YY1)
                   WRITE(700,5060) XX0,YY1,ZZ1
                END DO

 1500           CONTINUE

             END DO

             if (NWDIR(M).eq.0) then
!              write(*,*) 'vcavplt_sc.f:229',M
!              stop
               NWDIR(M) = 1
             end if
             IDR=NWDIR(M)
             NSC=NNWC(M)+NSPS(M,IDR)

             IF(NSC.GT.0) THEN
C..............Blade T.E.
                NT1=N0(1)-NH
                NB1=NHP-N0(2)+1
                CAVTOP1=CONMRX(1,N0(1),M)+
     *               HTP(NT1,M,1)*VECMRX(1,N0(1),M)
                CAVTOP2=CONMRX(2,N0(1),M)+
     *               HTP(NT1,M,1)*VECMRX(2,N0(1),M)
                CAVTOP3=CONMRX(3,N0(1),M)+
     *               HTP(NT1,M,1)*VECMRX(3,N0(1),M)
                CAVBOT1=CONMRX(1,N0(2),M)+
     *               HTP(NB1,M,2)*VECMRX(1,N0(2),M)
                CAVBOT2=CONMRX(2,N0(2),M)+
     *               HTP(NB1,M,2)*VECMRX(2,N0(2),M)
                CAVBOT3=CONMRX(3,N0(2),M)+
     *               HTP(NB1,M,2)*VECMRX(3,N0(2),M)

                IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,CAVTOP1,CAVTOP2)
                CALL ROTATE(CAVTOP2,CAVTOP3,YY1,ZZ1,DT1)
                IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,CAVTOP1,YY1)
                WRITE(700,5060) CAVTOP1,YY1,ZZ1

                IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,CAVBOT1,CAVBOT2)
                CALL ROTATE(CAVBOT2,CAVBOT3,YY1,ZZ1,DT1)
                IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,CAVBOT1,YY1)
                WRITE(700,5060) CAVBOT1,YY1,ZZ1

C..............SR T.E.
                CAVTOP1=CONMRX(1,NCP,M)+
     *               HTP(NHP,M,1)*VECMRX(1,NCP,M)
                CAVTOP2=CONMRX(2,NCP,M)+
     *               HTP(NHP,M,1)*VECMRX(2,NCP,M)
                CAVTOP3=CONMRX(3,NCP,M)+
     *               HTP(NHP,M,1)*VECMRX(3,NCP,M)
                CAVBOT1=CONMRX(1,1,M)+
     *               HTP(NHP,M,2)*VECMRX(1,1,M)
                CAVBOT2=CONMRX(2,1,M)+
     *               HTP(NHP,M,2)*VECMRX(2,1,M)
                CAVBOT3=CONMRX(3,1,M)+
     *               HTP(NHP,M,2)*VECMRX(3,1,M)
                WKCAM1=(CAVTOP1+CAVBOT1)/2.
                WKCAM2=(CAVTOP2+CAVBOT2)/2.
                WKCAM3=(CAVTOP3+CAVBOT3)/2.
                HTEHF=HTWP(1,M)/2.
                CAVTOP1=WKCAM1+HTEHF*WAKVEC(1,1,M)
                CAVTOP2=WKCAM2+HTEHF*WAKVEC(2,1,M)
                CAVTOP3=WKCAM3+HTEHF*WAKVEC(3,1,M)
                CAVBOT1=WKCAM1-HTEHF*WAKVEC(1,1,M)
                CAVBOT2=WKCAM2-HTEHF*WAKVEC(2,1,M)
                CAVBOT3=WKCAM3-HTEHF*WAKVEC(3,1,M)

                IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,CAVTOP1,CAVTOP2)
                CALL ROTATE(CAVTOP2,CAVTOP3,YY1,ZZ1,DT1)
                IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,CAVTOP1,YY1)
                WRITE(700,5060) CAVTOP1,YY1,ZZ1

                IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,CAVBOT1,CAVBOT2)
                CALL ROTATE(CAVBOT2,CAVBOT3,YY1,ZZ1,DT1)
                IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,CAVBOT1,YY1)
                WRITE(700,5060) CAVBOT1,YY1,ZZ1

                DO N=1,NSC
                   N1=(MR-M)*NTRA+N
                   IF(N.EQ.NNWC(M)+1) THEN
                      WKCAM1=WKCAM1-DZW(N,M)*ULW(N1,1)*FLS(M,IDR)
                      WKCAM2=WKCAM2-DZW(N,M)*ULW(N1,2)*FLS(M,IDR)
                      WKCAM3=WKCAM3-DZW(N,M)*ULW(N1,3)*FLS(M,IDR)
                      HTEHF=0.
                   ELSE
                      WKCAM1=WKCAM1-DZW(N,M)*ULW(N1,1)
                      WKCAM2=WKCAM2-DZW(N,M)*ULW(N1,2)
                      WKCAM3=WKCAM3-DZW(N,M)*ULW(N1,3)
                      HTEHF=HTWP(N+1,M)/2.
                   END IF
                   CAVTOP1=WKCAM1+HTEHF*WAKVEC(1,N+1,M)
                   CAVTOP2=WKCAM2+HTEHF*WAKVEC(2,N+1,M)
                   CAVTOP3=WKCAM3+HTEHF*WAKVEC(3,N+1,M)
                   CAVBOT1=WKCAM1-HTEHF*WAKVEC(1,N+1,M)
                   CAVBOT2=WKCAM2-HTEHF*WAKVEC(2,N+1,M)
                   CAVBOT3=WKCAM3-HTEHF*WAKVEC(3,N+1,M)

                 IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,CAVTOP1,CAVTOP2)
                 CALL ROTATE(CAVTOP2,CAVTOP3,YY1,ZZ1,DT1)
                 IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,CAVTOP1,YY1)
                 WRITE(700,5060) CAVTOP1,YY1,ZZ1

                 IF(ISP .NE. 0) CALL ROTATE2(-1,SPANGLE,CAVBOT1,CAVBOT2)
                 CALL ROTATE(CAVBOT2,CAVBOT3,YY1,ZZ1,DT1)
                 IF(ISP .NE. 0) CALL ROTATE2(0,SPANGLE,CAVBOT1,YY1)
                 WRITE(700,5060) CAVBOT1,YY1,ZZ1
                END DO
             END IF
          END DO

 2000     CONTINUE

C-----------------------------------------------------------------------
C     Establishing connectivity
C-----------------------------------------------------------------------
C........plot blade root
          DO N=1,NCOLD
             WRITE(700,5080) N, N+1, N+1, N
          END DO

C........plot blade TE
          IMM = NCOLDP
          DO M=1,MR
             WRITE(700,5080) IMM+M,IMM+M+1,
     *            IMM+M+1+MRP,IMM+M+MRP
          END DO

C........plot blade tip
          IMM = IMM + MRP*2
          DO N=1,NCOLD
             WRITE(700,5080) IMM+N, IMM+N+1, IMM+N+1, IMM+N
          END DO

C........plot blade LE
          IMM = IMM+NCOLDP
          DO M=1,MR
             WRITE(700,5080) IMM+M,IMM+M+1,IMM+M+1,IMM+M
          END DO

C........plot free surface
          IF(ISP.EQ.1) THEN
             IMM=IMM+MRP
             WRITE(700,5080) IMM+1,IMM+2,IMM+3,IMM+4
          END IF

          IF(ISP.EQ.1.AND.NBW.EQ.0)  GO TO 3000

C........Plot cavity
          ICNT=NPTS0+1
          DO M=1,MR
             DO IDR=1,2
                IF(ISP.EQ.0) THEN
                   N00=JCV(M,IDR)
                ELSE    
                   N00=0
                   IF(ISUB1(M,IDR,IDXREV).GT.0) THEN
                      N00=JCV(M,IDR)
                   END IF
                END IF

                DO N=1,N00
                   WRITE(700,5080) ICNT,ICNT+1,ICNT+3,ICNT+2
                   ICNT=ICNT+2
                END DO
                IF(N00.GT.0) ICNT=ICNT+2
             END DO
             
             IF(NNWC(M).GT.0) THEN
                DO N=1,NNWC(M)+NSPS(M,NWDIR(M))+1
                   WRITE(700,5080) ICNT,ICNT+1,ICNT+3,ICNT+2
                   ICNT=ICNT+2
                END DO
                ICNT=ICNT+2
             END IF
          END DO

 3000     CONTINUE
       END IF

C-----------------------------------------------------------------------
C     If IFILE=1, close file 700 at the end of every revolution.
C-----------------------------------------------------------------------
       IF(IFILE.EQ.1.AND.IDXREV.EQ.NTPREV) CLOSE(700)
       IF(IFILE.EQ.1.AND.ISTEADY.EQ.0) CLOSE(700)

       IF(IPR.EQ.1.AND.IDXREV.EQ.NTPREV) CLOSE(700)
       
       RETURN
       END
                

