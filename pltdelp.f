       SUBROUTINE PLTDELP
************************************************************************
*     Subroutine created to plot change of DELP with time in wakesheet.*
*                                                                      *
*     Date              Comments                                       *
*     ----------------  --------------------------------------------   *
*     JY071700          subroutine created.                            *
************************************************************************

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       CHARACTER*30 FNWPLT

       IF(NREV.EQ.NTREV) THEN
       
C-----------------------------------------------------------------------
C      Open output file, *.delp
C-----------------------------------------------------------------------
       IF(IDXREV.EQ.1) THEN
          CALL CHRLEN(FN,LENCH)
          FNWPLT=FN(1:LENCH)//'.delp'
          OPEN(631,FILE=FNWPLT,STATUS='UNKNOWN')
          WRITE(631,*) 'VARIABLES="x","y","z","delp"'
       END IF

C-----------------------------------------------------------------------
C      Read DELP
C-----------------------------------------------------------------------
       IREC=NTPOS(1)-1
       IF(IREC.LE.0)IREC=IREC+TWOPI/DELTAT

       NREAD=NWMIN*MR
       CALL READ2(46,IREC,TEMP5,NREAD)

       DO L=NWMINFW,1,-1
          DO M=MR,1,-1
             IDX0=MR*(L-1)+M
             IDX1=MR*L+M
             IF(L.EQ.NWMINFW) THEN
                TEMP5(IDX1)=ZERO
             ELSE
                TEMP5(IDX1)=TEMP5(IDX0)
             END IF
          END DO
       END DO
       DO M=MR,1,-1
          TEMP5(M)=DELP(M)
       END DO

C-----------------------------------------------------------------------
C      Generate contour plot of DELP
C-----------------------------------------------------------------------
       DT1=-DELTAT*(IDXREV-1)

       DO M=1,MR+1
          DO N=1,NWMIN+1
             YW1(N,M)=YW(N,M)*COS(DT1)-ZW(N,M)*SIN(DT1)
             ZW1(N,M)=YW(N,M)*SIN(DT1)+ZW(N,M)*COS(DT1)
          END DO
       END DO

       NEL=0
       NPTS=0
       DO M=1,MR
          IF(ISP.EQ.1) THEN
             IF(MSW(M,IDXREV).GT.0) THEN
                NEL=NEL+MSW(M,IDXREV)
                NPTS=NPTS+(MSW(M,IDXREV)+1)*2
             END IF
          ELSE
             NEL=NEL+NWMIN
             NPTS=NPTS+(NWMIN+1)*2
          END IF
       END DO

 1000  FORMAT(1x,'ZONE T="Time=',F4.0,'", N=',I5,', E=',I5,
     *      ' F=FEPOINT, ET=QUADRILATERAL')
 1010  FORMAT(4(1X,I4))
 1015  FORMAT(4(1X,F12.8))

       IF(NEL.GT.0) THEN
          IF(ISP.EQ.1) THEN
             WRITE(631,1000) TT(IDXREV),NPTS+4,NEL+1

             WRITE(631,1015) -1.,YFS,-1.2,0.
             WRITE(631,1015)  1.,YFS,-1.2,0.
             WRITE(631,1015)  1.,YFS, 1.2,0.             
             WRITE(631,1015) -1.,YFS, 1.2,0.
             
             DO M=1,MR
             IF(MSW(M,IDXREV).GT.0) THEN
                DO L=1,MSW(M,IDXREV)+1
                   IDX=MR*(L-1)+M
                   WRITE(631,1015) XW(L,M),YW1(L,M),ZW1(L,M),TEMP5(IDX)
                   WRITE(631,1015) XW(L,M+1),YW1(L,M+1),ZW1(L,M+1),
     *                  TEMP5(IDX)
                END DO
             END IF
             END DO

             WRITE(631,1010) 1,2,3,4
             IEL1=4
             IEL=0
             DO M=1,MR
                IF(MSW(M,IDXREV).GT.0) THEN
                   DO L=1,MSW(M,IDXREV)
                      IEL=IEL+1
                      WRITE(631,1010) 2*IEL-1+IEL1,2*IEL+IEL1,
     *                     2*IEL+2+IEL1,2*IEL+1+IEL1
                   END DO
                   IEL=IEL+1
                END IF
             END DO
          ELSE
             WRITE(631,1000) TT(IDXREV),NPTS,NEL

             DO M=1,MR
                DO L=1,NWMIN+1
                   IDX=MR*(L-1)+M
                   WRITE(631,1015) XW(L,M),YW1(L,M),ZW1(L,M),TEMP5(IDX)
                   WRITE(631,1015) XW(L,M+1),YW1(L,M+1),ZW1(L,M+1),
     *                  TEMP5(IDX)
                END DO
             END DO

             IEL1=0
             IEL=0
             DO M=1,MR
                DO L=1,NWMIN
                   IEL=IEL+1
                   WRITE(631,1010) 2*IEL-1+IEL1,2*IEL+IEL1,
     *                  2*IEL+2+IEL1,2*IEL+1+IEL1
                END DO
                IEL=IEL+1
             END DO
          END IF

       END IF

       IF(IDXREV.EQ.NTPREV) CLOSE(631)

       END IF

       RETURN
       END
