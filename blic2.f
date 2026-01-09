
      SUBROUTINE BLIC2(SIG,RHSM,TOL,NORD,NBLOCK,NPERB,
     * ITMAX,IPASS,NNR,NBLK,ITMX)
************************************************************************
*     This subroutine solves the greens's formula matrix in a BLock    *
*     ITerative method.                                                *
*                                                                      *
*     NNR  = ORDER OF MATRIX [A]  (# of elements in the panel code)    *
*     NNB  = ORDER OF THE SUBMATRIX  (Matrix Block)  IN MATRIX [A]     *
*     NBLK = TOTAL NUMBER OF BLOCKS IN MATRIX [A]                      *
*     ITMX = TOTAL NUMBER OF ITERATIONS + 1  (Same as ITMAX+1 )        *
*                                                                      *
*     Replace the dimension of AD,F,IPVT,Q11,Q12,P,RHSD and WORK by    *
*     'NNB' if : NNB > ITMX                                            *
*                                                                      *
*                                                                      *
*     Date     Comment or Revision                                     *
*     -------- -------------------                                     *
*     CM102897 A lot of variables were converted to Double precision.  *
*     This is in an attempt to keep the code from crashing in          *
*     the matrix solver so much.  The code crashes less with           *
*     these changes, but still has problems because some of            *
*     variables that are read into the matrix solver are               *
*     originally single precision. Changes were also made in           *
*     decomp and dsolve.  It may be possible to increase               *
*     matrix stability by having more of the inputs to this            *
*     subroutine calculated as double precision.                       *
*     CM121597 Changed COND to double precision                        *
*                                                                      *
************************************************************************
      USE MEMSOL
      INCLUDE 'PARAM.INC'
      DOUBLEPRECISION AD, Q11, Q12, WORK, RHSDBL,COND


!     PARAMETER(NNR=NTZ,ITMX=101)
!     PARAMETER(NBLK=NBLKMAX)
C     COMMON /MEMSOL/ AINF(NNR,NNR)

CVV   INTENT::IN
      DIMENSION RHSM(NNR),SIG(NNR)
      DIMENSION NPERB(NBLK)

C      DIMENSION ATEM(NNR),DSIG(NNR),RES(NNR),RESNEW(NNR),              
C     * SIGNEW(NNR),F(NNB),
C     * P(NNB,NNB),RHSD(NNB),
C     * URS1(2*NNR*ITMX),URS2(2*NNR*ITMX)
C      DIMENSION IDSTRT(NBLK),IDEND(NBLK),ILSTRT(NBLK),NSL(NBLK),
C     * NSU(NBLK),NPERSQ(NBLK),NPERB(NBLK),IPVT1(ITMX),
C     * IPVT(NNB,NBLK)
C      DIMENSION AD(NNB*NNB,NBLK),Q11(NNB,NNB),Q12(NNB),WORK(NNB),
C     * RHSDBL(NNB)

CVV
      ALLOCATABLE :: ATEM(:),DSIG(:),RES(:),RESNEW(:),SIGNEW(:)
      ALLOCATABLE :: F(:),P(:,:),RHSD(:)
      ALLOCATABLE :: URS1(:),URS2(:)
      ALLOCATABLE :: IDSTRT(:),IDEND(:),ILSTRT(:),NSL(:)
      ALLOCATABLE :: NSU(:),NPERSQ(:),IPVT1(:)
      ALLOCATABLE :: IPVT(:,:)
      ALLOCATABLE :: AD(:,:),Q11(:,:),Q12(:),WORK(:),RHSDBL(:)

      SAVE IDSTRT
C      SAVE ILSTRT
C      SAVE IDEND 
C      SAVE NSL
C      SAVE NSU
C      SAVE ATEM,DSIG,RES,RESNEW,SIGNEW
C      SAVE F,P,RHSD
C      SAVE URS1,URS2
C      SAVE NPERSQ,IPVT1,IPVT
C      SAVE AD,Q11,Q12,WORK
C      SAVE RHSDBL  

!     write(*,*) 'ITMX=',ITMX

      NNB=0
      DO NB=1,NBLOCK
       IF(NPERB(NB).GT.NNB) NNB=NPERB(NB)
      ENDDO

      MAXDIM=MAX0(NNB,ITMX)
      NPERMX=NNB
      IF(.NOT.ALLOCATED(IDSTRT)) THEN 
       ALLOCATE(IDSTRT(NBLK),ILSTRT(NBLK))
       ALLOCATE(IDEND(NBLK))
       ALLOCATE(NSL(NBLK))
       ALLOCATE(NSU(NBLK))
       ALLOCATE(ATEM(NNR),DSIG(NNR),RES(NNR),RESNEW(NNR),SIGNEW(NNR))
       ALLOCATE(F(NNB),P(NNB,NNB),RHSD(NNB))
       ALLOCATE(URS1(2*NNR*ITMX),URS2(2*NNR*ITMX))
       ALLOCATE(NPERSQ(NBLK),IPVT1(ITMX))
       ALLOCATE(IPVT(NNB,NBLK))
       ALLOCATE(AD(NNB*NNB,NBLK),Q11(NNB,NNB),Q12(NNB),WORK(NNB))
       ALLOCATE(RHSDBL(NNB))
      ENDIF

CVV
     
CP****
C      IF(IPASS.GT.1) THEN 
C      GO TO 65
C      ENDIF
CP****
      IDSTRT(1)=1
      IDEND(1)=IDSTRT(1)+NPERB(1)-1
      ILSTRT(1)=IDEND(1)+1
      NSL(1)=NORD-IDEND(1)
      NSU(1)=IDSTRT(1)-1
      NPERSQ(1)=NPERB(1)*NPERB(1)
      DO 20 NB=2,NBLOCK
       IDSTRT(NB)=IDSTRT(NB-1)+NPERB(NB-1)
       IDEND(NB)=IDSTRT(NB)+NPERB(NB)-1
       ILSTRT(NB)=IDEND(NB)+1
       NSL(NB)=NORD-IDEND(NB)
       NSU(NB)=IDSTRT(NB)-1
       NPERSQ(NB)=NPERB(NB)*NPERB(NB)
 20   CONTINUE
 9010 FORMAT(2X,'BLOCK TOO BIG IN SOLVER, NB and NPER : ',2I5)

C-----------------------------------------------------------------------
C     Diagonal blocks
C-----------------------------------------------------------------------
      DO 50 NB=1,NBLOCK
       ID0=IDSTRT(NB)-1
       DO 40 JCOL=1,NPERB(NB)
        KST=(JCOL-1)*NPERB(NB)
        KSS=IDSTRT(NB)+JCOL-1
        DO 30 I=1,NPERB(NB)
C---------Convert AD to double precision from ALHS CM102897-------------
         AD(KST+I,NB)=DBLE(ALHS(ID0+I,KSS))
 30     CONTINUE
 40    CONTINUE
 50   CONTINUE
C     
      DO 60 NB=1,NBLOCK
       CALL DECOMP(NPERB(NB),NPERB(NB),AD(1,NB),COND,
     *  IPVT(1,NB),WORK)
C      WRITE(*,9020) NB,NPERB(NB),COND
 60   CONTINUE
 9020 FORMAT(2X,'NB=',I4,3X,'NPER=',I4,3X,'COND.NO.=',1PD11.4)

C-----------------------------------------------------------------------
C     Preprocess of the solutions (by the first guess of RES)
C-----------------------------------------------------------------------
 65   CONTINUE

      IT=0
      SUM=0.0
      DO 70 N=1,NORD
       DSIG(N)=0.0
       SIG(N)=0.0
       DUM=RHSM(N)
       RES(N)=DUM
       SUM=SUM+DUM*DUM
 70   CONTINUE
      P(1,1)=SUM
c     
c.....URS1?
      URS1(1)=IT
      DO 80 J=1,NORD
       URS1(1+J)=RES(J)
 80   CONTINUE
      DO 90 J=1,NORD
       URS1(NORD+1+J)=DSIG(J)
 90   CONTINUE
      
C-----------------------------------------------------------------------
C     Starting point of the iteration
C-----------------------------------------------------------------------
 1000 CONTINUE

C-----------------------------------------------------------------------
C     Lower blocks
C-----------------------------------------------------------------------
      DO 170 NB=1,NBLOCK
       DO 110 IROW=1,NPERB(NB)
        RHSD(IROW)=RES(IDSTRT(NB)-1+IROW)
C---------Need to have rhsd as double precision before it passes to-----
C---------dsolve CM102897-----------------------------------------------
        RHSDBL(IROW)=DBLE(RHSD(IROW))
 110   CONTINUE
       CALL DSOLVE(NPERB(NB),NPERB(NB),AD(1,NB),RHSDBL,IPVT(1,NB))
       DO 120 IROW=1,NPERB(NB)
c     write(*,*) irow, nb,rhsdbl(irow),idstrt(nb)-1+irow
        DSIG(IDSTRT(NB)-1+IROW)=SNGL(RHSDBL(IROW))
C        write(*,*)irow,IDSTRT(NB)-1+IROW,size(dsig)
 120   CONTINUE
       DO 160 JM=IDSTRT(NB),IDEND(NB)
        DO 150 J=1,NSL(NB)
         RES(ILSTRT(NB)-1+J)=RES(ILSTRT(NB)-1+J)
     *    -ALHS(ILSTRT(NB)-1+J,JM)*DSIG(JM)
 150    CONTINUE
 160   CONTINUE
 170  CONTINUE
      DO 180 N=1,NORD
       RES(N)=0.0
 180  CONTINUE

C-----------------------------------------------------------------------
C     UPPER blocks
C-----------------------------------------------------------------------
      DO 230 NB=2,NBLOCK
       DO 220 JM=IDSTRT(NB),IDEND(NB)
        DO 210 IROW=1,NSU(NB)
         RES(IROW)=RES(IROW)-ALHS(IROW,JM)*DSIG(JM)
 210    CONTINUE
 220   CONTINUE
 230  CONTINUE
C     
C.....The next lines may net be needed
      DO 240 IROW=IDSTRT(NBLOCK),IDEND(NBLOCK)
       RES(IROW)=0.0
 240  CONTINUE

C-----------------------------------------------------------------------
C     Establish matrix P
C-----------------------------------------------------------------------
      IT=IT+1
      SUM=0.0
      DO 250 N=1,NORD
       SIGNEW(N)=SIG(N)+DSIG(N)
       DUM=RES(N)
       RESNEW(N)=RES(N)
       SUM=SUM+RES(N)*RES(N)
 250  CONTINUE
      P(IT+1,IT+1)=SUM
      ITPREV=URS1(1)

      I1=0
      DO 280 IRES=1,ITPREV+1
       SUM=0.0
       DO 260 N=1,NORD
        SUM=SUM+URS1(I1+1+N)*RES(N)
 260   CONTINUE
       I1=I1+NORD
       P(IT+1,IRES)=SUM
       P(IRES,IT+1)=SUM
 280  CONTINUE
C-----------------------------------------------------------------------
C     Establish matrix Q
C-----------------------------------------------------------------------

      DO 310 I=1,IT
       DO 300 J=1,IT
C---------Calculate Q11 as double precision CM102897--------------------
        Q11(I,J)=DBLE(P(I,J)-P(IT+1,I)-P(IT+1,J)+P(IT+1,IT+1))
300   CONTINUE
C---------Q12 is also going to be double precision CM102897-------------
       Q12(I)=DBLE(P(IT+1,I)-P(IT+1,IT+1))
       Q12(I)=-Q12(I)
 310  CONTINUE
C-----------------------------------------------------------------------
C     Establish matrix F
C-----------------------------------------------------------------------
      CALL DECOMP(NNB,IT,Q11,COND,IPVT1,WORK)
C      WRITE(777,9030) IT,COND
      CALL DSOLVE(NNB,IT,Q11,Q12,IPVT1)
      SUM=0.0
      DO 330 I=1,IT
C---------F is single precision  CM102897-------------------------------
       F(I)=SNGL(Q12(I))
       SUM=SUM+F(I)
 330  CONTINUE
      F(IT+1)=1.0-SUM
 9030 FORMAT(2X,'IT=',I4,3X,'Q11 COND.NO.=',1PD11.4)
C-----------------------------------------------------------------------
C     Update solutions I
C-----------------------------------------------------------------------
      DO 340 N=1,NORD
       SIG(N)=0.0
       RES(N)=0.0
 340  CONTINUE
C    
      I1=0
      I2=0
      ITPREV=URS1(1)
      URS2(1)=IT
      I1=I1+1
      I2=I2+1
      DO 400 ITC=1,ITPREV+1
       FMULT=F(ITC)
       DO 350 J=1,NORD
        ATEM(J)=URS1(I1+J)
 350   CONTINUE
       I1=I1+NORD
       DO 360 J=1,NORD
        URS2(I2+J)=ATEM(J)
 360   CONTINUE
       I2=I2+NORD
       DO 370 I=1,NORD
        RES(I)=RES(I)+FMULT*ATEM(I)
 370   CONTINUE
 400  CONTINUE
C    
      RESERR=0.0
      FMULT=F(IT+1)
      DO 410 I=1,NORD
       RES(I)=RES(I)+F(IT+1)*RESNEW(I)
       RESERR=AMAX1(RESERR,ABS(RES(I)))
 410  CONTINUE
      DO 420 J=1,NORD
       URS2(I2+J)=RES(J)
 420  CONTINUE
      I2=I2+NORD
!     WRITE(*,9040) IT,RESERR
 9040 FORMAT(2X,'IT=',I4,3X,'MAX RESIDUAL=',1PE11.4)
C-----------------------------------------------------------------------
C     Update solutions II
C-----------------------------------------------------------------------
      DO 500 ITC=1,ITPREV+1
       DO 450 J=1,NORD
        ATEM(J)=URS1(I1+J)
 450   CONTINUE
       I1=I1+NORD
       DO 460 J=1,NORD
        URS2(I2+J)=ATEM(J)
 460   CONTINUE
       I2=I2+NORD
       DO 470 I=1,NORD
        SIG(I)=SIG(I)+F(ITC)*ATEM(I)
 470   CONTINUE
 500  CONTINUE

C     
      DO 510 I=1,NORD
       SIG(I)=SIG(I)+F(IT+1)*SIGNEW(I)
 510  CONTINUE
      DO 520 J=1,NORD
       URS2(I2+J)=SIG(J)
 520  CONTINUE
      I2=I2+NORD
C     
      IF(RESERR.LE.TOL) THEN
       DEALLOCATE (ATEM,DSIG,RES,RESNEW,SIGNEW)
       DEALLOCATE (F,P,RHSD)
       DEALLOCATE (URS1,URS2)
       DEALLOCATE (IDSTRT,IDEND,ILSTRT,NSL)
       DEALLOCATE (NSU,NPERSQ,IPVT1)
       DEALLOCATE (IPVT)
       DEALLOCATE (AD,Q11,Q12,WORK,RHSDBL)
       RETURN
      END IF
C     
      IF(IT.GE.ITMAX) THEN
       WRITE(*,9050) ITMAX,RESERR
       DEALLOCATE (ATEM,DSIG,RES,RESNEW,SIGNEW)
       DEALLOCATE (F,P,RHSD)
       DEALLOCATE (URS1,URS2)
       DEALLOCATE (IDSTRT,IDEND,ILSTRT,NSL)
       DEALLOCATE (NSU,NPERSQ,IPVT1)
       DEALLOCATE (IPVT)
       DEALLOCATE (AD,Q11,Q12,WORK,RHSDBL)
       RETURN
      END IF
 9050 FORMAT(/2X,'NO CONVERGENCE IN BLIC2 AFTER',I4,
     * ' ITERATIONS',/2X,'MAX RESIDUAL=',1PE11.4)
C-----------------------------------------------------------------------
C     Final clean
C-----------------------------------------------------------------------
      SUM=0.0
      DO 530 ITC=1,IT+1
       SUM=SUM+F(ITC)*P(ITC,IT+1)
 530  CONTINUE
      DO 540 ITC=1,IT+1
       P(IT+1,ITC)=SUM
       P(ITC,IT+1)=SUM
 540  CONTINUE
      IMAX=MAX0(I1,I2)
      DO 550 J=1,IMAX
       TEMP=URS1(J)
       URS1(J)=URS2(J)
       URS2(J)=TEMP
 550  CONTINUE
      DO 600 N=1,NORD
       DSIG(N)=0.0
 600  CONTINUE
      
      GO TO 1000

      END
