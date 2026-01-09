      SUBROUTINE CAVPLT2D
************************************************************************
*                                                                      *
*      This subroutine plots the 2-D cavity shape of a particular      *
*      strip on the projected X-Z plane.  It's mainly use for checking *
*      the convergence of the cavity shape.                            *
*                                                                      *
*      Date     Revision or Note                                       *
*      -------  ----------------                                       *
*      JY092898 Subroutine created.                                    *
*      JY030599 Modified subroutine to allow cavity to grow on both the*
*               back and face of the foil.                             *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
*      JY060200 Modified subroutine to plot the correct form of the    *
*               supercavity.                                           *
*                                                                      *
************************************************************************
      
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

       CHARACTER*30 FN2D

C-----------------------------------------------------------------------
C      This is the strip number that you wish to plot
C-----------------------------------------------------------------------
       M=5

C-----------------------------------------------------------------------
C      This is the step in each revolution that you wish to check 
C      convergence (e.g. NANG=6 implies at a blade angle of 30 degrees
C-----------------------------------------------------------------------
       IF(ISTEADY.EQ.0) THEN
          NANG=1
       ELSE
          NANG=NTPREV
       END IF
       
C-----------------------------------------------------------------------
C      Open output file, *.2d
C-----------------------------------------------------------------------
       IF(NTSTEP.EQ.1) THEN
          CALL CHRLEN(FN,LENCH)
          FN2D=FN(1:LENCH)//'.2d'
          OPEN(730,FILE=FN2D,STATUS='UNKNOWN')

C-----------------------------------------------------------------------
C      Plot the 2-D blade
C-----------------------------------------------------------------------
          write(730,*) 'Variables = "x","z"'
          WRITE(730,5000) M
 5000     FORMAT(1X,'ZONE T="FOIL AT STRIP #',I2,'"')
          
          DO 10 N=1,NCP
             WRITE(730,*) CONMRX(1,N,M),CONMRX(3,N,M)
 10       CONTINUE
       END IF

       IF(IDXREV.EQ.NANG) THEN

 5019  FORMAT(1X,'ZONE T="BOTH: T=',I3,'"')
 5020  FORMAT(1X,'ZONE T="BACK: T=',I3,'"')
 5021  FORMAT(1X,'ZONE T="FACE: T=',I3,'"')

C-----------------------------------------------------------------------
C      Plot back cavity if applicable (LE -> TE)
C-----------------------------------------------------------------------
       IF(JCV(M,1).GT.0) THEN

          IF(NWC(M,1).GT.0.AND.NWC(M,2).GT.0) THEN
             WRITE(730,5019) NTSTEP
          ELSE
             WRITE(730,5020) NTSTEP      
          END IF

          M0M=M0(M,1)
          M0M1=M0M-1
          JCAV=JCV(M,1)

          DO N=1,JCAV+1
             N1=M0M+(N-1)
             WRITE(730,*) CONMRX(1,N1,M)+HT(N,M,1)*VECMRX(1,N1,M),
     *            CONMRX(3,N1,M)+HT(N,M,1)*VECMRX(3,N1,M)
          END DO

          IF(NSPP(M,1).NE.0) THEN
             N=JCAV+2
             NN1=N1
             NN2=N1+1
             DXB=CONMRX(1,NN2,M)-CONMRX(1,NN1,M)
             DZB=CONMRX(3,NN2,M)-CONMRX(3,NN1,M)
             WRITE(730,*) CONMRX(1,NN1,M)+DXB*FLP(M,1),
     *            CONMRX(3,NN1,M)+DZB*FLP(M,1)
          END IF
       END IF

C-----------------------------------------------------------------------
C      Plot supercavity on wake
C-----------------------------------------------------------------------
       IF(NNWC(M).GT.0) THEN
          IF(NWC(M,1).EQ.0) THEN
             WRITE(730,5021) NTSTEP
             WRITE(730,*) CONMRX(1,NCP,M),CONMRX(3,NCP,M)
          END IF

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
             WKCAM3=CONMRX(3,NCP,M)
          ELSE IF(NWC(M,1).GT.0.AND.NWC(M,2).EQ.0) THEN
             CALL CALHTW(1,M,HTE)
             HTEHF=HTE/2.
             WKCAM1=CONMRX(1,NCP,M)-HTEHF*VELW(JW,1)
             WKCAM3=CONMRX(3,NCP,M)-HTEHF*VELW(JW,3)
             TP=0.5
             TB=0.5
          ELSE IF(NWC(M,2).GT.0.AND.NWC(M,1).EQ.0) THEN
             CALL CALHTW(2,M,HTE)
             HTEHF=HTE/2.
             WKCAM1=CONMRX(1,1,M)+HTEHF*VELW(JW,1)
             WKCAM3=CONMRX(3,1,M)+HTEHF*VELW(JW,3)
             TP=0.5
             TB=0.5
          END IF

          DO N=2,NNWC(M)+1
             N1=(MR-M)*NTRA+N-1
             WKCAM1=-DZW(N-1,M)*ULW(N1,1)+WKCAM1
             WKCAM3=-DZW(N-1,M)*ULW(N1,3)+WKCAM3
             CAVTOP1=WKCAM1-HT(IDX-1+N,M,IDR)*TP*VELW(N1,1)
             CAVTOP3=WKCAM3-HT(IDX-1+N,M,IDR)*TP*VELW(N1,3)
             WRITE(730,*) CAVTOP1,CAVTOP3  
          END DO

          IF(NSPS(M,IDR).NE.0) THEN
             N=NNWC(M)+NSPS(M,IDR)+1
             WKCAM1=-DZW(N-1,M)*FLS(M,IDR)*ULW(N1+1,1)+WKCAM1
             WKCAM3=-DZW(N-1,M)*FLS(M,IDR)*ULW(N1+1,3)+WKCAM3
             WRITE(730,*) WKCAM1,WKCAM3
             WKCAM1=WKCAM1+DZW(N-1,M)*FLS(M,IDR)*ULW(N1+1,1)
             WKCAM3=WKCAM3+DZW(N-1,M)*FLS(M,IDR)*ULW(N1+1,3)
          END IF

          DO N=NNWC(M)+1,2,-1
             N1=(MR-M)*NTRA+N-1
             CAVBOT1=WKCAM1+HT(IDX-1+N,M,IDR)*TB*VELW(N1,1)
             CAVBOT3=WKCAM3+HT(IDX-1+N,M,IDR)*TB*VELW(N1,3)
             WRITE(730,*) CAVBOT1,CAVBOT3      
             WKCAM1=WKCAM1+DZW(N-1,M)*ULW(N1,1)
             WKCAM3=WKCAM3+DZW(N-1,M)*ULW(N1,3)
          END DO

          IF(NWC(M,2).EQ.0) THEN
             WRITE(730,*) CONMRX(1,1,M),CONMRX(3,1,M)
          END IF
       END IF
       
C-----------------------------------------------------------------------
C      Plot face cavity if applicable (TE -> LE)
C-----------------------------------------------------------------------
       IF(JCV(M,2).GT.0) THEN
          IF(NWC(M,2).EQ.0) WRITE(730,5021) NTSTEP

          M0M=M0(M,2)
          M0M1=M0M-1
          JCAV=JCV(M,2)

          IF(NSPP(M,2).NE.0) THEN
             N=JCAV+2
             NN1=M0M-JCAV
             NN2=NN1-1
             DXB=CONMRX(1,NN2,M)-CONMRX(1,NN1,M)
             DZB=CONMRX(3,NN2,M)-CONMRX(3,NN1,M)
             WRITE(730,*) CONMRX(1,NN1,M)+DXB*FLP(M,2),
     *            CONMRX(3,NN1,M)+DZB*FLP(M,2)
          END IF

          DO N=JCAV+1,1,-1
             N1=M0M-(N-1)
             WRITE(730,*) CONMRX(1,N1,M)+HT(N,M,2)*VECMRX(1,N1,M),
     *            CONMRX(3,N1,M)+HT(N,M,2)*VECMRX(3,N1,M)
          END DO
       END IF

       END IF

       IF(NTSTEP.EQ.NCTIME) CLOSE(730)

       RETURN
       END
