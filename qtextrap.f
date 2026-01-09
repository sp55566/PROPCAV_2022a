      SUBROUTINE QTEXTRAP(M,IDR)
************************************************************************
*   left split-panel source uses a polynomial extrapolation            *
*   of the form a(sl-s)^2+b(sl-s)+c/sqrt(sl-s)+d, dropping one term at *
*   at time (starting with the highest order) as the number of cavity  *
*   panels drops below 4 (actually, the last to go is the square-root  *
*   singular term)                                                     *
*   moved here from delr.f on 06-23-92 NF                              *
*   moved here from PSFCAV on 06-26-92 NF                              *
*                                                                      *
*   DATE             CORRECTIONS/COMMENTS                              *
*   -------------    ------------------------------------------------  *
*   11/27/98 JY      Made several corrections relating to face         *
*                    cavitation (see CfaceJY).                         *
*   02/06/99 JY      modified subroutine to allow cavity to grow on    *
*                    both the back and face of the foil.               *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
************************************************************************
      INCLUDE 'PUFCAV.INC'

      M0M=M0(M,IDR)
      M0M1=M0M-1
      JCAV=JCV(M,IDR)
      
      IF(IDR.EQ.1) THEN
         K=1
         ISF=0
      ELSE IF(IDR.EQ.2) THEN
         K=-1
         ISF=1
      END IF
      AK1=FLOAT(K)
      
      N1=M0M1+K*JCAV+ISF
      NN1=N1+K*1

C      IF(JCAV.GE.4)THEN
C.......4-point extrapolation...........................................
C         S3=AK1*(SPZ(N1-K*2,M)-SPZ(N1-K*3,M))
C         S2=AK1*(SPZ(N1-K*1,M)-SPZ(N1-K*3,M))
C         S1=AK1*(SPZ(N1,M)-SPZ(N1-K*3,M))
C         S0=AK1*(SPZ(NN1,M)-SPZ(N1-K*3,M))
C         SL=S0+.25*DZL(M,IDR)
C         IF(M.LE.MRTIP) THEN
C            CALL QEXT1(S0,S1,S2,S3,SL,QT(M,1,IDR),QT(M,2,IDR),
C     *           QT(M,3,IDR),QT(M,4,IDR))
C         ELSE
C            CALL CUBEXT(S0,S1,S2,S3,QT(M,1,IDR),QT(M,2,IDR),
C     *           QT(M,3,IDR),QT(M,4,IDR))
C         END IF         
C      ELSEIF(JCAV.EQ.3)THEN
C.......3-point extrapolation...........................................
C         S2=AK1*(SPZ(N1-K*1,M)-SPZ(N1-K*2,M))
C         S1=AK1*(SPZ(N1,M)-SPZ(N1-K*2,M))
C         S0=AK1*(SPZ(NN1,M)-SPZ(N1-K*2,M))
C         SL=S0+.25*DZL(M,IDR)
C         IF(M.LE.MRTIP) THEN
C            CALL QEXT2(S0,S1,S2,SL,QT(M,1,IDR),QT(M,2,IDR),
C     *           QT(M,3,IDR))
C         ELSE
C            CALL QUADEXT(S0,S1,S2,QT(M,1,IDR),QT(M,2,IDR),
C     *           QT(M,3,IDR))
C         END IF
C         QT(M,4,IDR)=ZERO
C      ELSEIF(JCAV.EQ.2)THEN
      IF(JCAV.GE.2)THEN
C.......2-point extrapolation...........................................
         S1=AK1*(SPZ(N1,M)-SPZ(N1-K*1,M))
         S0=AK1*(SPZ(NN1,M)-SPZ(N1-K*1,M))
         SL=S0+.25*DZL(M,IDR)
         IF(M.LE.MRTIP) THEN
            CALL QEXT3(S0,S1,SL,QT(M,1,IDR),QT(M,2,IDR))
         ELSE
            QT(M,1,IDR)=S0/S1
            QT(M,2,IDR)=ONE-S0/S1
         END IF
         QT(M,3,IDR)=ZERO
         QT(M,4,IDR)=ZERO
      ELSEIF(JCAV.EQ.1)THEN
C.......1-point extrapolation...........................................
         QT(M,1,IDR)=ONE
         QT(M,2,IDR)=ZERO
         QT(M,3,IDR)=ZERO
         QT(M,4,IDR)=ZERO
      ENDIF

      RETURN
C<<<<<<<<<<<<<<<<<<<<<<<<<end of subroutine QTEXTRAP>>>>>>>>>>>>>>>>>>>>
      END



