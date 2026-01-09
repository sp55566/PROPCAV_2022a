       SUBROUTINE CALHTW(IDR,M,HTW)
************************************************************************
*      This subroutine calculates the cavity height at the blade       *
*      trailing edge with respect to the normal of the first wake      *
*      panel.                                                          *
*                                                                      *
*      Date      Revision or Note                                      *
*      --------  ----------------                                      *
*      JY060200  Subroutine created.                                   *
*      JY071100  Subroutine modified.  Now it uses up to 6 points for  *
*                interpolation.                                        *
*                                                                      *
************************************************************************
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       DIMENSION XX(NBHZP),YY(NBHZP),NN1(NBHZP),JB1(NBHZ),SS1(NBHZ),
     *      HH1(NBHZP),CT1(NBHZP),ST1(NBHZP),CB1(NBHZ),SB1(NBHZ)

       IF(ISC.EQ.0) THEN
          JCAV1=JCV(M,IDR)+1
       ELSE
          IF(ISP.EQ.0) THEN
             IF(JCV(M,IDR).EQ.0) THEN
                JCAV1=NSR2P
             ELSE
                JCAV1=NHP-NLEP(M,IDXREV,IDR)
             END IF
          ELSE
             IF(IDR.EQ.1) THEN
                JCAV1=MAX(0,IC(2,M,IDXREV)-N0(1)+1)
                IF(JCV(M,1).GT.0) JCAV1=JCAV1+JCV(M,1)
             ELSE
                JCAV1=MAX(0,N0(2)-IW(1,M,IDXREV))
                IF(JCV(M,2).GT.0) JCAV1=JCAV1+JCV(M,2)
             END IF
          END IF
       END IF

       NT1=MIN(10,JCAV1)
       JW=(MR-M)*NTRA+1

       IF(IDR.EQ.1) THEN
          IC1=NCP
          DO N=NT1,1,-1
             NN1(N)=IC1
             IF(N.NE.NT1) THEN
                JB1(N)=INDEXB(IC1,M)
                SS1(N)=DZ(IC1,M)
             END IF
             IC1=IC1-1

          END DO
       ELSE IF(IDR.EQ.2) THEN
          IC1=1
          DO N=NT1,1,-1
             NN1(N)=IC1
             IF(N.NE.NT1) THEN
                JB1(N)=INDEXB(IC1-1,M)
                SS1(N)=DZ(IC1-1,M)
             END IF
             IC1=IC1+1
          END DO
       END IF

       IF(JCAV1.LE.2.OR.HT(JCAV1,M,IDR).LE.ZERO) THEN
          DOT=0.
          DO KK=1,3
             DOT=DOT+VEL(JB1(NT1-1),KK)*VELW(JW,KK)
          END DO
          HTW=HT(JCAV1,M,IDR)*ABS(DOT)
       ELSE
          DO N=NT1,1,-1
             HH1(N)=HT(JCAV1-NT1+N,M,IDR)
             CT1(N)=ABS(VECMRX(1,NN1(N),M)*ULW(JW,1)+
     *            VECMRX(2,NN1(N),M)*ULW(JW,2)+
     *            VECMRX(3,NN1(N),M)*ULW(JW,3))
             ST1(N)=ABS(VECMRX(1,NN1(N),M)*VELW(JW,1)+
     *            VECMRX(2,NN1(N),M)*VELW(JW,2)+
     *            VECMRX(3,NN1(N),M)*VELW(JW,3))
             IF(N.NE.NT1) THEN
                CB1(N)=ABS(UL(JB1(N),1)*ULW(JW,1)+
     *               UL(JB1(N),2)*ULW(JW,2)+
     *               UL(JB1(N),3)*ULW(JW,3))
                SB1(N)=ABS(UL(JB1(N),1)*VELW(JW,1)+
     *               UL(JB1(N),2)*VELW(JW,2)+
     *               UL(JB1(N),3)*VELW(JW,3))
             END IF
          END DO

          DO N=1,NT1
             XX(N)=HH1(N)*CT1(N)
             DO NN=1,N-1
                XX(N)=XX(N)+SS1(NN)*CB1(NN)
             END DO
             YY(N)=HH1(N)*ST1(N)
             DO NN=N,NT1-1
                YY(N)=YY(N)+SS1(NN)*SB1(NN)
             END DO
          END DO

          XTWW=0.
          DO N=1,NT1-1
             XTWW=XTWW+SS1(N)*CB1(N)
          END DO

          CALL QUADFIT(NT1,XX,YY,A0,A1,A2)
          HTW=A0+A1*XTWW+A2*XTWW**2.

       END IF

       RETURN
       END



       SUBROUTINE QUADFIT(N,X,Y,A0,A1,A2)

       DIMENSION X(N),Y(N),A(3,3),B(3),INDX(3)
       
       NORDER=2

       SUMX=0.
       SUMY=0.
       SUMX2=0.
       SUMX3=0.
       SUMX4=0.
       SUMXY=0.
       SUMX2Y=0.
       DO I=1,N
          SUMX=SUMX+X(I)
          SUMX2=SUMX2+X(I)**2.
          SUMX3=SUMX3+X(I)**3.
          SUMX4=SUMX4+X(I)**4.
          SUMY=SUMY+Y(I)
          SUMXY=SUMXY+X(I)*Y(I)
          SUMX2Y=SUMX2Y+X(I)**2.*Y(I)
       END DO       
       A(1,1)=FLOAT(N)
       A(1,2)=SUMX
       A(1,3)=SUMX2
       A(2,1)=SUMX
       A(2,2)=SUMX2
       A(2,3)=SUMX3
       A(3,1)=SUMX2
       A(3,2)=SUMX3
       A(3,3)=SUMX4
       B(1)=SUMY
       B(2)=SUMXY
       B(3)=SUMX2Y

       CALL SLUDCMP(A,3,3,INDX,D)
       CALL SLUBKSB(A,3,3,INDX,B)

       A0=B(1)
       A1=B(2)
       A2=B(3)

       RETURN
       END



      SUBROUTINE sludcmp(a,n,np,indx,d)
      PARAMETER (NMAX=500,TINY=1.0e-20)
      DIMENSION a(np,np), vv(nmax)
      dimension indx(n)

      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n

          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)

          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END

      SUBROUTINE slubksb(a,n,np,indx,b)
      DIMENSION indx(n)
      DIMENSION a(np,np),b(n)

      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END



