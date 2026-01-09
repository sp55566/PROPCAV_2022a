      SUBROUTINE CAVPLT2D_SC
************************************************************************
*                                                                      *
*      This subroutine plots the 2-D cavity shape of a particular      *
*      strip on the projected X-Z plane.  It's mainly use for checking *
*      the convergence of the cavity shape.                            *
*                                                                      *
*      Date     Revision or Note                                       *
*      -------  ----------------                                       *
************************************************************************
      
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

       DIMENSION CAVTOP(2,NSCZP,2)
       CHARACTER*30 FN2D

C-----------------------------------------------------------------------
C      This is the strip number that you wish to plot
C-----------------------------------------------------------------------
       M=MR/2
c       M=1

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
c       IF(NTSTEP.EQ.1.AND.ITER.EQ.1) THEN
          CALL CHRLEN(FN,LENCH)
          FN2D=FN(1:LENCH)//'.2d'
          OPEN(730,FILE=FN2D,STATUS='UNKNOWN')
          
          write(730,*) 'Variables = "x","z"'
          WRITE(730,5000) M
 5000     FORMAT(1X,'ZONE T="FOIL AT STRIP #',I2,'"')
          
          DO 10 N=N0(2),N0(1)
             WRITE(730,*) CONMRX(1,N,M),CONMRX(3,N,M)
 10       CONTINUE
          N=N0(2)
          WRITE(730,*) CONMRX(1,N,M),CONMRX(3,N,M)
       END IF

       IF(IDXREV.EQ.NANG) THEN
 5020     FORMAT(1X,'ZONE T="T=',I3,'"')
          WRITE(730,5020) NTSTEP
c          WRITE(730,5020) ITER

          WKCAM1=0.
          WKCAM3=0.
          DO IDR=1,2
             DO N=1,NHP
                IF(IDR.EQ.1) THEN
                   N1=NH+N
                ELSE
                   N1=NHP+1-N
                END IF
                CAVTOP(1,N,IDR)=CONMRX(1,N1,M)+
     *               HTP(N,M,IDR)*VECMRX(1,N1,M)
                CAVTOP(2,N,IDR)=CONMRX(3,N1,M)+
     *               HTP(N,M,IDR)*VECMRX(3,N1,M)
                IF(N.EQ.NHP) THEN
                   WKCAM1=WKCAM1+CONMRX(1,N1,M)+
     *                  HTP(N,M,IDR)*VECMRX(1,N1,M)
                   WKCAM3=WKCAM3+CONMRX(3,N1,M)+
     *                  HTP(N,M,IDR)*VECMRX(3,N1,M)
                END IF
             END DO
          END DO
          
          WKCAM1=WKCAM1/2.
          WKCAM3=WKCAM3/2.
          IDR=NWDIR(M)
          IF(NNWC(M).GT.0) THEN
             DO N=1,NNWC(M)+NSPS(M,IDR)                
                N1=(MR-M)*NTRA+N
                IF(N.EQ.NWC(M,IDR)+1) THEN
                   WKCAM1=WKCAM1-DZW(N,M)*ULW(N1,1)*FLS(M,IDR)
                   WKCAM3=WKCAM3-DZW(N,M)*ULW(N1,3)*FLS(M,IDR)
                   HTEHF=0.
                ELSE
                   WKCAM1=WKCAM1-DZW(N,M)*ULW(N1,1)
                   WKCAM3=WKCAM3-DZW(N,M)*ULW(N1,3)
                   HTEHF=HTWP(N+1,M)/2.
                END IF
                N00=NHP+N
                CAVTOP(1,N00,1)=WKCAM1+HTEHF*WAKVEC(1,N+1,M)
                CAVTOP(2,N00,1)=WKCAM3+HTEHF*WAKVEC(3,N+1,M)
                CAVTOP(1,N00,2)=WKCAM1-HTEHF*WAKVEC(1,N+1,M)
                CAVTOP(2,N00,2)=WKCAM3-HTEHF*WAKVEC(3,N+1,M)
             END DO
          END IF

          N1=NHP+NNWC(M)+NSPS(M,IDR) 
          DO N=1,N1
             WRITE(730,*) CAVTOP(1,N,1),CAVTOP(2,N,1)
          END DO
          DO N=N1,1,-1
             WRITE(730,*) CAVTOP(1,N,2),CAVTOP(2,N,2)
          END DO
       END IF

       RETURN
       END
