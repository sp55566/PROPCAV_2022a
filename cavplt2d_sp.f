       SUBROUTINE CAVPLT2D_SP
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
*                                                                      *
************************************************************************
      
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

       CHARACTER*30 FN2D1,FN2D2

       IF(NREV.EQ.NTREV) THEN

C-----------------------------------------------------------------------
C      Open output file, *.2d
C-----------------------------------------------------------------------
       IF(IDXREV.EQ.1) THEN
          CALL CHRLEN(FN,LENCH)
          FN2D1=FN(1:LENCH)//'.zy'
          FN2D2=FN(1:LENCH)//'.zx'
          OPEN(730,FILE=FN2D1,STATUS='UNKNOWN')
          WRITE(730,*) 'VARIABLES="Z","Y"'
          OPEN(731,FILE=FN2D2,STATUS='UNKNOWN')
          WRITE(731,*) 'VARIABLES="Z","X"'
       END IF

C-----------------------------------------------------------------------
C      Plot the cavity up to the T.E. or to where it's cut off by 
C      the blade surface.
C-----------------------------------------------------------------------
 5000  FORMAT(1X,'ZONE T="T=',F7.1,',M=',I2,'"')

C       DO 10 M=1,MR
       M=MR/2
       DT1=-TSTEP
       SINT=SIN(DT1)
       COST=COS(DT1)

       IF(NPERM(M).GT.0) THEN

          WRITE(730,5000) TT(IDXREV),M
          WRITE(731,5000) TT(IDXREV),M

C........Plot the blade section.
          IF(ISC.EQ.0) THEN
             NW1=IW(1,M,IDXREV)
             NW2=IW(2,M,IDXREV)+1
             NC1=IC(1,M,IDXREV)
             NC2=IC(2,M,IDXREV)+1
          ELSE
             NW1=MAX(IW(1,M,IDXREV),N0(2))
             NW2=IW(2,M,IDXREV)+1
             NC1=IC(1,M,IDXREV)
             NC2=MIN(IC(2,M,IDXREV)+1,N0(1))
          END IF

          IF(ISUB1(M,1,IDXREV).GT.0) THEN
             N=NC1
             XX0=CONMRX(1,N,M)
             YY1=CONMRX(2,N,M)*COST-CONMRX(3,N,M)*SINT
             ZZ1=CONMRX(2,N,M)*SINT+CONMRX(3,N,M)*COST
             WRITE(730,*) ZZ1,YY1
             WRITE(731,*) ZZ1,XX0
          END IF
          IF(ISUB1(M,2,IDXREV).GT.0) THEN
             DO N=NW2,NW1,-1
                XX0=CONMRX(1,N,M)
                YY1=CONMRX(2,N,M)*COST-CONMRX(3,N,M)*SINT
                ZZ1=CONMRX(2,N,M)*SINT+CONMRX(3,N,M)*COST
                WRITE(730,*) ZZ1,YY1
                WRITE(731,*) ZZ1,XX0
             END DO
          END IF
          IF(ISUB1(M,1,IDXREV).GT.0) THEN
             DO N=NC2,NC1,-1
                XX0=CONMRX(1,N,M)
                YY1=CONMRX(2,N,M)*COST-CONMRX(3,N,M)*SINT
                ZZ1=CONMRX(2,N,M)*SINT+CONMRX(3,N,M)*COST
                WRITE(730,*) ZZ1,YY1
                WRITE(731,*) ZZ1,XX0
             END DO
          END IF

C........Plot the cavity if it exists.
          IF(ISUB1(M,1,IDXREV).GT.0) THEN
             JCAV=JCV(M,1)
             IF(JCAV.EQ.0) GO TO 10
          
             M0M=M0(M,1)
             M0M1=M0M-1

             DO N=NC1,M0M
                XX0=CONMRX(1,N,M)
                YY1=CONMRX(2,N,M)*COST-CONMRX(3,N,M)*SINT
                ZZ1=CONMRX(2,N,M)*SINT+CONMRX(3,N,M)*COST
                WRITE(730,*) ZZ1,YY1
                WRITE(731,*) ZZ1,XX0
             END DO

             DO N=1,JCAV+1
                N1=M0M1+N
                XX0=CONMRX(1,N1,M)+HT(N,M,1)*VECMRX(1,N1,M)
                YY0=CONMRX(2,N1,M)+HT(N,M,1)*VECMRX(2,N1,M)
                ZZ0=CONMRX(3,N1,M)+HT(N,M,1)*VECMRX(3,N1,M)
                YY1=YY0*COST-ZZ0*SINT
                ZZ1=YY0*SINT+ZZ0*COST
                WRITE(730,*) ZZ1,YY1
                WRITE(731,*) ZZ1,XX0
             END DO 

 10          CONTINUE
          END IF
       END IF

       IF(NTSTEP.EQ.NCTIME) THEN
          CLOSE(730)
          CLOSE(731)
       END IF

       END IF

       RETURN
       END
