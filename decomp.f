      SUBROUTINE DECOMP(NDIM,N,A,COND,IPVT,WORK)
************************************************************************
*     This subroutine........                                          *
*                                                                      *
*     Date     Revision or Comment                                     *
*     -------- -------------------                                     *
*     CM102897 A,B,T are doubleprecision.  See header in blic2.        *
*     CM121597 Changed COND, ZNORM, and YNORM to double precision.     *
*              The code was crashing during the calculation of COND.   *
*              It will crash eventually anyway with such a high        *
*              conditional number.                                     *
*     CM121697 Put a control statement in the subroutine asking if the *
*              user wants to continue if the conditional number is high*
*                                                                      *
************************************************************************
!Allen Du 01/06/2018 add the option of coupling with cavopt3d 
      use m_cavopt

      DOUBLEPRECISION A, ANORM,T,WORK,ZNORM,YNORM,COND

      DIMENSION A(NDIM,N)
      DIMENSION WORK(N)
      INTEGER IPVT(N)
      IPVT(N)=1

      
      IF(N.EQ.1) GO TO 80

      NM1=N-1
      ANORM=0.0

      DO 10 J=1,N
         T=0.0

         DO 5 I=1,N
            T=T+DABS(A(I,J))
C           write(*,*)A(I,J)
 5       CONTINUE

         IF(T.GT.ANORM) ANORM=T

 10   CONTINUE

      DO 35 K=1,NM1
         KP1=K+1
         M=K

         DO 15 I=KP1,N

            IF(DABS(A(I,K)).GT.DABS(A(M,K))) M=I

 15      CONTINUE

         IPVT(K)=M
         IF(M.NE.K) IPVT(N)=-IPVT(N)

         T=A(M,K)
         A(M,K)=A(K,K)
         A(K,K)=T

         IF(T.EQ.0.0) GO TO 35

         DO 20 I=KP1,N
            A(I,K)=-A(I,K)/T
 20      CONTINUE

         DO 30 J=KP1,N
            T=A(M,J)
            A(M,J)=A(K,J)
            A(K,J)=T

            IF(T.EQ.0.0) GO TO 30

            DO 25 I=KP1,N
               A(I,J)=A(I,J)+A(I,K)*T
 25         CONTINUE

 30      CONTINUE

 35   CONTINUE


      DO 50 K=1,N
         T=0.0

         IF(K.EQ.1) GO TO 45

         KM1=K-1

         DO 40 I=1,KM1
            T=T+A(I,K)*WORK(I)
 40      CONTINUE

 45      EK=1.0

         IF(T.LT.0.0) EK=-1.0

         IF(A(K,K).EQ.0.0) GO TO 90

         WORK(K)=-(DBLE(EK)+T)/A(K,K)

 50   CONTINUE

      DO 60 KB=1,NM1
         K=N-KB
         T=0.0
         KP1=K+1

         DO 55 I=KP1,N
            T=T+A(I,K)*WORK(K)
 55      CONTINUE

         WORK(K)=T
         M=IPVT(K)
         
         IF(M.EQ.K) GO TO 60

         T=WORK(M)
         WORK(M)=WORK(K)
         WORK(K)=T
 60   CONTINUE

      YNORM=0D0

      DO 65 I=1,N
         YNORM=YNORM+DABS(WORK(I))
 65   CONTINUE

      CALL DSOLVE(NDIM,N,A,WORK,IPVT)

      ZNORM=0D0

      DO 70 I=1,N
         ZNORM=ZNORM+DABS(WORK(I))
 70   CONTINUE

      COND=ANORM*ZNORM/YNORM

      IF(COND.LT.1.0) COND=1.0

      RETURN

 80   COND=1.0

      IF(A(1,1).NE.0.0) RETURN

 90   COND=1.0D+32

C-----------------------------------------------------------------------
C     Ask user if he/she wants to continue with high cond number      CM
C-----------------------------------------------------------------------
      WRITE(*,*) 'NOTICE: High conditional number in matrix solver.'
!Allen Du 01/06/2018 add the option of coupling with cavopt3d 
      open(70401,file="cavopt_parameter.dat")
      write(70401,*)"1"
      close(70401)
CJY      WRITE(*,*) 'Trying another matrix solver.'
CJY      IBAD = 1

      RETURN

      END
