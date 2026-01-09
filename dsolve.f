      SUBROUTINE DSOLVE(NDIM,N,A,B,IPVT)
************************************************************************
*     Subroutine Dsolve.  I don't know what this subroutine does!      *
*     Maybe I will tell you when I find out :).  CM                    *
*                                                                      *
*     Date          Revision or Comment                                *
*     --------      ---------------------                              *
*     CM100397      Header Created                                     *
*     CM122897      A,B,T are doubleprecision.  See header in blic2.   *
*                                                                      *
************************************************************************

      DOUBLEPRECISION A,B,T

      INTEGER IPVT(N)
      DIMENSION A(NDIM,N),B(N)

      IF(N.EQ.1) GO TO 50

      NM1=N-1

      DO 20 K=1,NM1
         KP1=K+1
         M=IPVT(K)
         T=B(M)
         B(M)=B(K)
         B(K)=T

         DO 10 I=KP1,N
            B(I)=B(I)+A(I,K)*T
 10      CONTINUE

 20   CONTINUE

       DO 40 KB=1,NM1
         KM1=N-KB
         K=KM1+1
!s---YE TIAN --- 06/11/2013---
!        if (abs(A(K,K)).lt.1d-8) then
!          write(*,*) 'dsolve.f:40'
!          write(*,*) 'KB,KM1,K',KB, KM1, K, A(K,K)
!          write(*,*) 'NDIM,N', NDIM, N
!        end if
!e---YE TIAN --- 06/11/2013---
         B(K)=B(K)/A(K,K)
         T=-B(K)

         DO 30 I=1,KM1
            B(I)=B(I)+A(I,K)*T
 30      CONTINUE

 40   CONTINUE

 50   B(1)=B(1)/A(1,1)
      RETURN
C>>>>>>>>>>>>>>>>>>>>>>>>>>END OF DSOLVE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      END
