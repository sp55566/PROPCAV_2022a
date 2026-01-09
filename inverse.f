C     THE FOLLOWING SUBROUTINES ARE COPIED FROM CAV2DBL BY HONG SUN
C
C =====================================================
      SUBROUTINE INVRSE(A,AINV,N,NDIM)
C =====================================================
C
C-----FINDS THE INVERSE OF SQUARE MATRICES
C-----FINDS THE INVERSE OF THE LU DECOMPOSITION OF THE INPUT MATRIX
C-----CALL SUBROUTINES FACTOR, AND SUBST
c-----s.a.kinnas -- 12/5/99 -- change dim. to 3000 (>60X36=2160)
      DIMENSION A(NDIM,NDIM),W(N,N),B(10000),
     1          AINV(NDIM,NDIM),IPIVOT(10000),D(10000)
c-----s.a.kinnas -- 12/5/99 -- change limit of do loop to 3000
      DO 10 I=1,10000
           IPIVOT(I)=0
           D(I)=0.
           B(I)=0.
10    CONTINUE
c------12/5/99 ------------------------------------------------
      CALL FACTORV(A,W,IPIVOT,D,N,NDIM,IFLAG)
      IF (IFLAG.EQ.2) THEN
           WRITE(*,100)
100        FORMAT (' THE MATRIX IS SINGULAR')
           STOP
      ELSE
           CONTINUE
      ENDIF
      DO 1 I=1,N
           B(I)=0.
1     CONTINUE
      DO 3 J=1,N
           B(J)=1.
           CALL SUBSTV(W,B,D,IPIVOT,N)
           DO 2 I=1,N
                AINV(I,J)=D(I)
2          CONTINUE
      B(J)=0.
3     CONTINUE
      RETURN
      END
C
C =====================================================
      SUBROUTINE FACTORV(A,W,IPIVOT,D,N,NDIM,IFLAG)
C =====================================================
C
C-----FACTORS SQUARE MATRICES. SUBROUTINE SUBST USED FOR
C-----BACK SUBSTITUITION TO SOLVE SYSTEM OF EQUATIONS
      DIMENSION A(NDIM,NDIM),W(N,N),IPIVOT(*),D(*)
      IFLAG=1
      DO 10 I=1,N
      IPIVOT(I)=I
      ROWMAX=0.
      DO 9 J=1,N
      W(I,J)=A(I,J)
9     ROWMAX=AMAX1(ROWMAX,ABS(W(I,J)))
      IF(ROWMAX.EQ.0.) GO TO 999
10    D(I)=ROWMAX
      NM1=N-1
      IF(NM1.EQ.0) RETURN
      DO 20 K=1,NM1
      J=K
      KP1=K+1
      IP=IPIVOT(K)
      COLMAX=ABS(W(IP,K))/D(IP)
      DO 11 I=KP1,N
      IP=IPIVOT(I)
      AWIKOV=ABS(W(IP,K))/D(IP)
      IF(AWIKOV.LE.COLMAX) GO TO 11
      COLMAX=AWIKOV
      J=I
11    CONTINUE
      IF(COLMAX.EQ.0.) GO TO 999
      IPK=IPIVOT(J)
      IPIVOT(J)=IPIVOT(K)
      IPIVOT(K)=IPK
      DO 20 I=KP1,N
      IP=IPIVOT(I)
      W(IP,K)=W(IP,K)/W(IPK,K)
      RATIO=-W(IP,K)
      DO 20 J=KP1,N
20    W(IP,J)=RATIO*W(IPK,J)+W(IP,J)
      IF(W(IP,N).EQ.0.) GO TO 999
      RETURN
999   IFLAG=2
      RETURN
      END


C
C =====================================================
      SUBROUTINE SUBSTV(W,B,X,IPIVOT,N)
C =====================================================
C
C-----SUBROUTINE TO PERFORM BACK SUBSTITUTION AFTER CALLING
C-----SUBROUTINE FACTOR
      DIMENSION W(N,N),B(*),X(*),IPIVOT(*)
      IF(N.GT.1) GO TO 10
      X(1)=B(1)/W(1,1)
      RETURN
10    IP=IPIVOT(1)
      X(1)=B(IP)
      DO 15 K=2,N
      IP=IPIVOT(K)
      KM1=K-1
      SUM=0.
      DO 14 J=1,KM1
14    SUM=W(IP,J)*X(J)+SUM
15    X(K)=B(IP)-SUM
      X(N)=X(N)/W(IP,N)
      K=N
      DO 20 NP1MK=2,N
      KP1=K
      K=K-1
      IP=IPIVOT(K)
      SUM=0.0
      DO 19 J=KP1,N
19    SUM=W(IP,J)*X(J)+SUM
20    X(K)=(X(K)-SUM)/W(IP,K)
      RETURN
      END


C*********************************************************************


      subroutine matinv(a,n,np,indx,y)
c-----------------------------------------------------------
c   a - original matrix to invert
c   n - dimension actually contained in a
c   np - dimension of a declared in dimension statement (physical n)
c   indx - output vector which records the row permutation effected
c          by partial pivoting
c   y - inverse of a
c  note:  original a matrix is destroyed
c-----------------------------------------------------------------

      dimension a(np,np),y(np,np),indx(np)
      do 12 i=1,n
            do 11 j=1,n
                  y(i,j)=0.0
 11            continue
            y(i,i)=1.0
 12      continue
      call sludcmp(a,n,np,indx,d)
      do 13 j=1,n
            call slubksb(a,n,np,indx,y(1,j))
 13      continue
      return
      end




