       SUBROUTINE MVLS1D_2(NNODE,X,Z,XG,B)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER(NMAX=5000)
       DIMENSION X(NMAX),Z(NMAX),P(NMAX,3),A(3,3),B(3),INDX(3)       

       XM=0D0
       DO N=1,NNODE
          XM=DMAX1(XM,DABS(X(N)-XG))
       END DO

       DO N=1,NNODE
          P(N,1)=1D0
          P(N,2)=X(N)-XG
          P(N,3)=P(N,2)**2D0
       END DO

       C=6D0
       DO I=1,3
          DO J=1,3
             A(I,J)=0D0
          END DO
          B(I)=0D0
       END DO
       
       DO N=1,NNODE
          W=DEXP((-C**2D0)*(((X(N)-XG)/XM)**2D0))
          DO I=1,3
             DO J=1,3
                A(I,J)=A(I,J)+P(N,I)*P(N,J)*W
             END DO
             B(I)=B(I)+P(N,I)*W*Z(N)
          END DO
       END DO

       CALL LUDCMP(A,3,3,INDX,D)
       CALL LUBKSB(A,3,3,INDX,B)

       RETURN
       END


       SUBROUTINE MVLS1D_3(NNODE,X,Z,XG,B)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER(NMAX=5000)
       DIMENSION X(NMAX),Z(NMAX),P(NMAX,4),A(4,4),B(4),INDX(4)       

       XM=0D0
       DO N=1,NNODE
          XM=DMAX1(XM,DABS(X(N)-XG))
       END DO

       DO N=1,NNODE
          P(N,1)=1D0
          P(N,2)=X(N)-XG
          P(N,3)=P(N,2)**2D0
          P(N,4)=P(N,2)**3D0
       END DO

       C=6D0
       DO I=1,4
          DO J=1,4
             A(I,J)=0D0
          END DO
          B(I)=0D0
       END DO
       
       DO N=1,NNODE
          W=DEXP((-C**2D0)*(((X(N)-XG)/XM)**2D0))
          DO I=1,4
             DO J=1,4
                A(I,J)=A(I,J)+P(N,I)*P(N,J)*W
             END DO
             B(I)=B(I)+P(N,I)*W*Z(N)
          END DO
       END DO

       CALL LUDCMP(A,4,4,INDX,D)
       CALL LUBKSB(A,4,4,INDX,B)

       RETURN
       END


      SUBROUTINE ludcmp(a,n,np,indx,d)
      double precision d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      double precision aamax,dum,sum,vv(NMAX)
      dimension indx(n)

      d=1d0
      do 12 i=1,n
        aamax=0d0
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.0d0) pause 'singular matrix in ludcmp'
        vv(i)=1d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0d0
        do 16 i=j,n

          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
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
        if(a(j,j).eq.0d0) a(j,j)=TINY
        if(j.ne.n)then
          dum=1d0/a(j,j)

          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END

      SUBROUTINE lubksb(a,n,np,indx,b)
      dimension indx(n)
      double precision a(np,np),b(n),sum

      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0d0) then
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

