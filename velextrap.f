       SUBROUTINE VELEXTRAP
************************************************************************
*                                                                      *
*      This subroutine extrapolates velocities for M>MRTIP.            *
*                                                                      *
*      Author: Julie Young                                             *
*                                                                      *
*      Date        Revision/comments                                   *
*      --------    ---------------                                     *
*      JY012001    Subroutine created.                                 *
*               
C  Variable RZP(M) changed to HRZP(N,M)
************************************************************************

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       DIMENSION XX(3),YY1(3),C1(3)

       MR1=MRTIP-2
       MR2=MRTIP-1
       MR3=MRTIP
       IF(ISC.EQ.0) THEN
          NF=1
          NL=NC
       ELSE
          NF=N0(2)
          NL=N0(1)-1
       END IF

       DO N=NF,NL

C........Ux in global coordinates
          XX(1)=HRZP(N,MR1)
          YY1(1)=UXTOT(N,MR1)
          XX(2)=HRZP(N,MR2)
          YY1(2)=UXTOT(N,MR2)
          XX(3)=HRZP(N,MR3)
          YY1(3)=UXTOT(N,MR3)
          CALL QUADCOEF(XX,YY1,C1)
          DO M=MRTIP+1,MR
             UXTOT(N,M)=C1(1)+C1(2)*HRZP(N,M)+C1(3)*HRZP(N,M)**2.
          END DO          

C........Uy in global coordinates
          XX(1)=HRZP(N,MR1)
          YY1(1)=UYTOT(N,MR1)
          XX(2)=HRZP(N,MR2)
          YY1(2)=UYTOT(N,MR2)
          XX(3)=HRZP(N,MR3)
          YY1(3)=UYTOT(N,MR3)
          CALL QUADCOEF(XX,YY1,C1)
          DO M=MRTIP+1,MR
             UYTOT(N,M)=C1(1)+C1(2)*HRZP(N,M)+C1(3)*HRZP(N,M)**2.
          END DO

C........Uz in global coordinates
          XX(1)=HRZP(N,MR1)
          YY1(1)=UZTOT(N,MR1)
          XX(2)=HRZP(N,MR2)
          YY1(2)=UZTOT(N,MR2)
          XX(3)=HRZP(N,MR3)
          YY1(3)=UZTOT(N,MR3)
          CALL QUADCOEF(XX,YY1,C1)
          DO M=MRTIP+1,MR
             UZTOT(N,M)=C1(1)+C1(2)*HRZP(N,M)+C1(3)*HRZP(N,M)**2.
          END DO

       END DO

       RETURN
       END


       SUBROUTINE QUADCOEF(X,Y,C)
       DIMENSION X(3),Y(3),C(3),A(3,3),B(3)
       
       DO I=1,3
          A(I,1)=1.0
       END DO
       A(1,2)=X(1)
       A(1,3)=X(1)**2.
       A(2,2)=X(2)
       A(2,3)=X(2)**2.
       A(3,2)=X(3)
       A(3,3)=X(3)**2.
       B(1)=Y(1)
       B(2)=Y(2)
       B(3)=Y(3)

       CALL MTXSOL(3,A,B,C)

       RETURN
       END


       SUBROUTINE MTXSOL(N,A,F,X)
************************************************************************
*      This is a simple matrix solver using Doolittle LU decomposition *
*      and forward & backward substitution with no pivoting.           *
*                                                                      *
*      Input:  N      =    No. of equations                            *
*              A(N,N) =    LHS of system of equations                  *
*              F(N,N) =    RHS of system of equations                  *
*      Output: X(N)   =    solution vector                             *
*                                                                      *
*      Author: Julie Young                                             *
************************************************************************
       DIMENSION A(N,N),F(N),X(N),AL(N,N),U(N,N),Z(N)
C
C      LU decomposition using Doolittle's method
C
       DO K=1,N
          AL(K,K)=1.
          DO I=K,N
             SUM=0.
             DO L=1,K-1
                SUM=SUM+AL(K,L)*U(L,I)
             END DO
             U(K,I)=A(K,I)-SUM
          END DO
          DO I=K+1,N
             SUM=0.
             DO L=1,K-1
                SUM=SUM+AL(I,L)*U(L,K)
             END DO
             AL(I,K)=(A(I,K)-SUM)/U(K,K)
          END DO
       END DO
C
C      Forward substitution
C
       Z(1)=F(1)
       DO I=2,N
          SUM=0.
          DO K=1,I-1
             SUM=SUM+AL(I,K)*Z(K)
          END DO
          Z(I)=F(I)-SUM
       END DO
C
C      Backward substitution
C
       X(N)=Z(N)/U(N,N)
       DO I=N-1,1,-1
          SUM=0.
          DO K=I+1,N
             SUM=SUM+U(I,K)*X(K)
          END DO
          X(I)=(Z(I)-SUM)/U(I,I)
       END DO

       RETURN
       END
