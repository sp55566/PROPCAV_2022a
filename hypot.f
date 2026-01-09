C**************************** MIT PSF10 ********************************
C
C                PROPELLER STEADY FLOW ANALYSIS PROGRAM
C                        PANEL METHOD SOLUTION         
C
C                             Version 1.0
C
C       Copyright (c) Massachusetts Institute of Technology 1997
C
C                     Release Date: 1 January 1997
C
C***********************************************************************
      SUBROUTINE HYPOT(X1R,X2R,X3R,X4R,XCCR,FDR,FSR)
C**********************************************************************
C     DOUBLE PRECISION
C     Compute the potential due to sources and dipoles based on 
C      Morino's formula
C     -- x1, x2, x3 and x4 should be input in the clock wise direction
C     -- This subroutine CAN NOT calculate the influence functions of 
C        a TRIANGULAR panel
C     
C     10-04-88 C.Y.HSIN @MHL
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL X1R,X2R,X3R,X4R,XCCR,FDR,FSR 
      DIMENSION X1(3),X2(3),X3(3),X4(3),P(3),PC(3),P1(3),P2(3),P3(3),
     *          XCC(3),RC(3),XI0(4),ETA0(4),FDP(4),FSP(4),  
     *          R(3),A1(3),A2(3),U(3),XNC(3),XN(3),RXA1(3),RXA2(3),
     *          NNEG(4)
      DIMENSION X1R(3),X2R(3),X3R(3),X4R(3),XCCR(3) 
      DATA XI0 /-1.0D0,-1.0D0,1.0D0,1.0D0/
      DATA ETA0/-1.0D0,1.0D0,1.0D0,-1.0D0/
      DATA ZERO,QUAD,HALF,TWO,FOUR/0.0D0,0.25D0,0.50D0,
     *                                 2.0D0,4.0D0/
C
C.....Transfer to double precision
      DO 10 I=1,3
         X1(I)=DBLE(X1R(I))         
         X2(I)=DBLE(X2R(I))         
         X3(I)=DBLE(X3R(I))         
         X4(I)=DBLE(X4R(I))         
         XCC(I)=DBLE(XCCR(I))         
 10   CONTINUE
C      PI=3.14159265D0
      PI = DACOS(-1.D0)
      TWOPI=TWO*PI 
      FOURPI=FOUR*PI 
      PI2=HALF*PI
      DO 20 I=1,3
         PC(I)=QUAD*( X1(I)+X2(I)+X3(I)+X4(I) )
         P1(I)=QUAD*( X3(I)+X4(I)-X2(I)-X1(I) )           
         P2(I)=QUAD*( X3(I)-X4(I)+X2(I)-X1(I) )           
         P3(I)=QUAD*( X3(I)-X4(I)-X2(I)+X1(I) )           
 20   CONTINUE  
C
C.....normal vector at the center point
C
      XI=ZERO
      ETA=ZERO
      DO 30 I=1,3
         P(I)=PC(I)+XI*P1(I)+ETA*P2(I)+XI*ETA*P3(I)
         RC(I)=P(I)-XCC(I)            
         A1(I)=P1(I)+ETA*P3(I)
         A2(I)=P2(I)+XI*P3(I)
 30   CONTINUE
      CALL EXPROD(A1,A2,U)
      UL=DSQRT( U(1)*U(1)+U(2)*U(2)+U(3)*U(3) )
      RCL=DSQRT( RC(1)*RC(1)+RC(2)*RC(2)+RC(3)*RC(3) )
      DO 40 I=1,3 
         XNC(I)=U(I)/UL
 40   CONTINUE        
C
C.....Compute Is, Id at four corner points       
      DO 50 J=1,4
         XI=XI0(J)           
         ETA=ETA0(J)
         DO 60 I=1,3
            P(I)=PC(I)+XI*P1(I)+ETA*P2(I)+XI*ETA*P3(I)
            R(I)=P(I)-XCC(I)
            A1(I)=P1(I)+ETA*P3(I)
            A2(I)=P2(I)+XI*P3(I)
 60      CONTINUE
         CALL EXPROD(A1,A2,U)
         UL=DSQRT( U(1)*U(1)+U(2)*U(2)+U(3)*U(3) )
         DO 70 I=1,3 
            XN(I)=U(I)/UL
 70      CONTINUE
         CALL EXPROD(R,A1,RXA1)
         CALL EXPROD(R,A2,RXA2)
         CALL ENPROD(RXA1,RXA2,RA1RA2)
         CALL ENPROD(R,U,RU)
         CALL ENPROD(R,A1,RA1)
         CALL ENPROD(R,A2,RA2)
         CALL ENPROD(RXA1,XNC,RXA1NC)
         CALL ENPROD(RXA2,XNC,RXA2NC)
         CALL ENPROD(RC,XNC,RNC)
         RL=DSQRT( R(1)*R(1)+R(2)*R(2)+R(3)*R(3) )
         RXA1L=DSQRT( RXA1(1)*RXA1(1)+RXA1(2)*RXA1(2)+RXA1(3)*RXA1(3) )        
         RXA2L=DSQRT( RXA2(1)*RXA2(1)+RXA2(2)*RXA2(2)+RXA2(3)*RXA2(3) )        
         A1L=DSQRT( A1(1)*A1(1)+A1(2)*A1(2)+A1(3)*A1(3) )        
         A2L=DSQRT( A2(1)*A2(1)+A2(2)*A2(2)+A2(3)*A2(3) )        
         DEN=RL*RU
         FDP(J)=ATAN3( DEN,RA1RA2 ) 
C........Determine the quardrants
         IF(FDP(J).GE.ZERO) THEN
            NNEG(J)=0
         ELSE
            NNEG(J)=1
         END IF
         FS1=RNC*FDP(J)
         FS2=RXA1NC*DASINH(RA1/RXA1L)/A1L
         FS3=RXA2NC*DASINH(RA2/RXA2L)/A2L
         FSP(J)=-FS2+FS3
 50   CONTINUE
C
C     Compute potential of source and dipole
      FD=FDP(1)-FDP(2)+FDP(3)-FDP(4)
      MNEG=NNEG(1)+NNEG(2)+NNEG(3)+NNEG(4)
      IF(MNEG.EQ.1) THEN
         IF(NNEG(1).EQ.1.OR.NNEG(3).EQ.1) THEN
            FD=FD+TWOPI 
         END IF
      ELSE IF(MNEG.EQ.3) THEN
         IF(NNEG(1).EQ.0.OR.NNEG(3).EQ.0) THEN
            FD=FD-TWOPI 
         END IF
      END IF
      IF(FD.LT.-6.28318531) THEN
         FD=FOURPI+FD
      END IF
      FS1=FSP(1)-FSP(2)+FSP(3)-FSP(4)
      FS=-FD*RNC+FSP(1)-FSP(2)+FSP(3)-FSP(4)
      FDR=SNGL(FD)
      FSR=SNGL(FS)
      RETURN
C))))))))))))))))))))))) End of subroutine HYPOT  (((((((((((((((((((((
      END

