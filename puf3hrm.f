      SUBROUTINE PUF3HRM(NBLADE,NSR,IWET,XKTV,INSTEP)
c-----------------------------------------------------------------------
c       ccc
ccc    ccccccc
ccc  cccccccccccc
cccc ccccccccccc               +++++++
ccccccccccc                +++++++++++++++                         CC C
cccccccccccc             ++++++++++++++++++++                 CC CCCCCC
cc ccc cccc             ++++++++++++++++++++++             CCCCCCCCCCCC
**  *  * *    __________________________________________  CCCCCCCCCCC**
**  *  *      I- Massachusetts Institute of Technology-I  CCCC* CCC*CCC
**  **        ------------------------------------------    CCCCC *CCC*
*  **           III    III    III    III    III    III        CCCCC CCC
* **            III    III    III    III    III    III            * **
**              III    III    III    III    III    III             * **
*               III    III    III    III    III    III               **
*               III    III    III    III    III    III               **
*               III    III    III    III    III    III               **
*             ==========================================             **
*            ============================================            **
**          ==============================================          ***
*
C---------------------------- MIT PUF3-HRM ----------------------------
C
C                      HARMONIC ANALYSIS OF SIX
C                        COMPONENTS OF FORCE
C
C                  PROGRAMMER  CHANG-SUP LEE  (1978)
C
C           COPYRIGHT (C) MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C                              (DEC 1986)
C
C     DATE                   COMMENTS
C     -----------------      --------------------------------------
C     JY071399               This subroutine was modified for PROPCAV
C                            input.
C
C----------------------------------------------------------------------

      use m_P3DAT1
      INCLUDE 'PARAM.INC'
c!     parameter (n4=30)
c!     real,allocatable :: xktv(:,:)
      integer INSTEP
      DIMENSION XKTV(INSTEP,6)
c!     COMMON/P3DAT1/F(NSTEP,6),C(N4,NSTEP),S(N4,NSTEP),
c!    *     AO(6),A(N4,6),B(N4,6),
c!    *     FT(NSTEP,6),AMP(N4,6),PH(N4,6),HAMP(N4),HPH(N4),VAMP(N4),
c!    *     VPH(N4)
      COMMON/I3DAT1/NH

      if (.NOT.allocated(F)) then
        allocate(F(INSTEP,6),C(N4,INSTEP),S(N4,INSTEP),FT(INSTEP,6))
      end if

      DO K=1,NSR
         DO J=1,6
            F(K,J)=XKTV(K,J)
         END DO
      END DO

      IF(IWET.EQ.1) THEN
         OPEN(16,FILE='harmny.wet',STATUS='UNKNOWN')
      ELSE IF(IWET.EQ.0) THEN
         OPEN(16,FILE='harmny.cav',STATUS='UNKNOWN')
      END IF

      NH=NSR/4
      NSR2=NSR/2
      NBLH=NH/NBLADE
      STEP=6.283185/NSR
      DO 1 N=1,NH
      DO 1 K=1,NSR
      T=N*STEP*(K-1.0)
      C(N,K)=COS(T)
 1    S(N,K)=SIN(T)
      DO 3 J=1,6
      DO 3 N=1,NH
      AO(J)=0.0
      A(N,J)=0.0
      B(N,J)=0.0
      DO 5 K=1,NSR
      AO(J)=AO(J)+F(K,J)/NSR
      A(N,J)=A(N,J)+F(K,J)*C(N,K)/NSR2
      B(N,J)=B(N,J)+F(K,J)*S(N,K)/NSR2
 5    CONTINUE
      IF(N.EQ.NSR2) A(N,J)=A(N,J)/2.0
      AMP(N,J)=SQRT(A(N,J)*A(N,J)+B(N,J)*B(N,J))
      PH(N,J)=ATAN2(A(N,J),B(N,J))*57.296*(-1)+90.
      IF(PH(N,J).LT.0.0) PH(N,J)=PH(N,J)+360.
 3    CONTINUE
      WRITE(16,101)
 101  FORMAT(////,32X,'*** PROPCAV ***'//
     9           20X,'HARMONIC ANALYSIS OF SIX COMPONENTS OF FORCE'
     9/20X,44(1H-)//25X,'*** ACTING ON THE KEY BLADE ***'
     9//15X,'FX',13X,'FY',13X,'FZ',13X,'MX',13X,'MY',13X,
     9'MZ'//2X,'N-TH',6(4X,'AMP',5X,'PH',1X))
      N=0
      WRITE(16,103) N,(AO(J),J=1,6)
 103  FORMAT(/I6,8(F9.5,6X))
      DO 6 N=1,NH
 6    WRITE(16,102) N,(AMP(N,J),PH(N,J),J=1,6)
 102  FORMAT(I6,8(F9.5,F6.1))
CJY      IF(NSR.EQ.60) GOTO 7
CJY      DO 2 M=1,NBLH
CJY      N=M*NBLADE
CJY      DO 2 K=1,NSR
CJY      T=N*0.10472*(K-1.0)
CJY      C(N,K)=COS(T)
CJY 2    S(N,K)=SIN(T)
CJY 7    CONTINUE
      DO 10 J=1,4,3
      DO 10 K=1,NSR
      FT(K,J)=0.0
      DO 11 M=1,1
      N=M*NBLADE
 11   FT(K,J)=FT(K,J)+A(N,J)*C(N,K)+B(N,J)*S(N,K)
 10   FT(K,J)=FT(K,J)*NBLADE
      AO(1)=NBLADE*AO(1)
      AO(4)=AO(4)*NBLADE
      DO 122 N=1,NBLH
      PH(N,4)=PH(N*NBLADE,4)
      PH(N,1)=PH(N*NBLADE,1)
      AMP(N,1)=NBLADE*AMP(N*NBLADE,1)
 122  AMP(N,4)=NBLADE*AMP(N*NBLADE,4)
      WRITE(16,123)
 123  FORMAT(//25X,'*** ACTING ON PROPELLER SHAFT ***'//27X,
     9'**** THRUST HARMONIC ANALYSIS ****'/31X,
     9'N-BLH      AMP      PH')
      N=0
      WRITE(16,124) N,AO(1)
 124  FORMAT(33X,I2,F11.5,F8.1)
      DO 126 N=1,NBLH
 126  WRITE(16,124) N,AMP(N,1),PH(N,1)
      WRITE(16,127)
 127  FORMAT(//27X,'**** AXIAL TORQUE HARMONIC ANALYSIS ****'/31X,
     9'N-BLH      AMP      PH')
      N=0
      WRITE(16,124) N,AO(4)
      DO 128 N=1,NBLH
 128  WRITE(16,124) N,AMP(N,4),PH(N,4)
      DO 14 J=1,4,3
      DO 14 K=1,NSR
      FT(K,J)=0.0
      DO 15 M=1,NBLH
      N=M*NBLADE
 15   FT(K,J)=FT(K,J)+A(N,J)*C(N,K)+B(N,J)*S(N,K)
 14   FT(K,J)=FT(K,J)*NBLADE
      DO 16 J=1,4,3
      STEADY=AO(J)*NBLADE
      DO 16 K=1,NSR
  16   FT(K,J)=FT(K,J)+STEADY
C
C     ------ HORIZONTAL FORCE AND MOMENT CALCULATION ------
      WRITE(16,125)
 125  FORMAT(////27X,'**** HORIZONTAL FORCE HARMONIC ANALYSIS ****'/
     9 31X,'N-BLH      AMP      PH')
      CALL HORZFM(NBLADE,NSR,3,2,NBLH,STEADY)
      DO 35 K=1,NSR
 35   FT(K,1)=FT(K,2)+FT(K,3)-STEADY
      WRITE(16,134)
 134  FORMAT(////27X,'**** HORIZONTAL TORQUE HARMONIC ANALYSIS ****'/
     9 31X,'N-BLH      AMP      PH')
      CALL HORZFM(NBLADE,NSR,6,5,NBLH,STEADY)
      DO 36 K=1,NSR
 36   FT(K,1)=FT(K,6)+FT(K,5)-STEADY
C
C     ------ VERTICAL FORCE AND MOMENT CALCULATION -----
      WRITE(16,141)
 141  FORMAT(////27X,'**** VERTICAL FORCE HARMONIC ANALYSIS ****'/
     9 31X,'N-BLH      AMP      PH')
      CALL VERTFM(NBLADE,NSR,3,2,NBLH,STEADY)
      DO 37 K=1,NSR
 37   FT(K,1)=FT(K,3)+FT(K,2)-STEADY
      WRITE(16,150)
 150  FORMAT(////27X,'**** VERTICAL TORQUE HARMONIC ANALYSIS ****'/
     9 31X,'N-BLH      AMP      PH')
      CALL VERTFM(NBLADE,NSR,6,5,NBLH,STEADY)
      DO 38 K=1,NSR
 38   FT(K,1)=FT(K,6)+FT(K,5)-STEADY
      CLOSE(16)
      RETURN
      END


      SUBROUTINE HORZFM(NBLADE,NSR,JT,JR,NBLHIN,STEADY)
      use m_P3DAT1
      INCLUDE 'PARAM.INC'
c!     parameter (n4=30)
c!     COMMON/P3DAT1/F(NSTEP,6),C(N4,NSTEP),S(N4,NSTEP),
c!    *     AO(6),A(N4,6),B(N4,6),
c!    *     FT(NSTEP,6),AMP(N4,6),PH(N4,6),HAMP(N4),HPH(N4),VAMP(N4),
c!    *     VPH(N4)
      COMMON/I3DAT1/NH

      DO 23 K=1,NSR
      FT(K,JT)=A(1,JT)
      DO  20 M=1,NBLHIN
      N=M*NBLADE
      IF(NBLADE.LE.1) GOTO 21
      FT(K,JT)=FT(K,JT)+A(N-1,JT)*C(N,K)+B(N-1,JT)*S(N,K)
 21   IF(N.GE.NH) GOTO 20
      FT(K,JT)=FT(K,JT)+A(N+1,JT)*C(N,K)+B(N+1,JT)*S(N,K)
 20   CONTINUE
 23   FT(K,JT)=FT(K,JT)*NBLADE/2.0
      DO 33 K=1,NSR
      FT(K,JR)=B(1,JR)
      DO 30 M=1,NBLHIN
      N=M*NBLADE
      IF(NBLADE.LE.1) GOTO 31
      FT(K,JR)=FT(K,JR)+A(N-1,JR)*S(N,K)-B(N-1,JR)*C(N,K)
 31   IF(N.GE.NH) GOTO 30
      FT(K,JR)=FT(K,JR)-A(N+1,JR)*S(N,K)+B(N+1,JR)*C(N,K)
 30   CONTINUE
 33   FT(K,JR)=-FT(K,JR)*NBLADE/2.0
      STEADY=(A(1,JT)-B(1,JR))*NBLADE/2.0
      BUG=0.0
      CAT=0.0
      DO 40 M=1,NBLHIN
      N=M*NBLADE
      IF(NBLADE.LE.1) GOTO 41
      BUG=A(N-1,JT)+B(N-1,JR)
      CAT=B(N-1,JT)-A(N-1,JR)
 41   IF(N.GE.NH) GOTO 42
      BUG=BUG+A(N+1,JT)-B(N+1,JR)
      CAT=CAT+B(N+1,JT)+A(N+1,JR)
 42   HAMP(M)=SQRT(BUG**2+CAT**2)*NBLADE/2.0
      HPH(M)=ATAN2(BUG,CAT)*57.296*(-1)+90.
 40   IF(HPH(M).LT.0.0) HPH(M)=HPH(M)+360.
      M=0
      WRITE(16,101) M,STEADY
      DO 50 M=1,NBLHIN
 50   WRITE(16,101) M,HAMP(M),HPH(M)
 101  FORMAT(33X,I2,F11.5,F8.1)
      RETURN
      END


      SUBROUTINE VERTFM(NBLADE,NSR,JT,JR,NBLHIN,STEADY)
      use m_P3DAT1
      INCLUDE 'PARAM.INC'
c!     parameter (n4=30)
c!     COMMON/P3DAT1/F(NSTEP,6),C(N4,NSTEP),S(N4,NSTEP),
c!    *     AO(6),A(N4,6),B(N4,6),
c!    *     FT(NSTEP,6),AMP(N4,6),PH(N4,6),HAMP(N4),HPH(N4),VAMP(N4),
c!    *     VPH(N4)
      COMMON/I3DAT1/NH

      DO 43 K=1,NSR
      FT(K,JT)=B(1,JT)
      DO 40 M=1,NBLHIN
      N=M*NBLADE
      IF(NBLADE.LE.2) GOTO 41
      FT(K,JT)=FT(K,JT)+A(N-1,JT)*S(N,K)-B(N-1,JT)*C(N,K)
 41   IF(N.GE.NH) GOTO 40
      FT(K,JT)=FT(K,JT)-A(N+1,JT)*S(N,K)+B(N+1,JT)*C(N,K)
 40   CONTINUE
 43   FT(K,JT)=FT(K,JT)*NBLADE/2.0
      DO 53 K=1,NSR
      FT(K,JR)=A(1,JR)
      DO 50 M=1,NBLHIN
      N=M*NBLADE
      IF(NBLADE.LE.2) GOTO 51
      FT(K,JR)=FT(K,JR)+A(N-1,JR)*C(N,K)+B(N-1,JR)*S(N,K)
 51   IF(N.GE.NH) GOTO 50
      FT(K,JR)=FT(K,JR)+A(N+1,JR)*C(N,K)+B(N+1,JR)*S(N,K)
 50   CONTINUE
 53   FT(K,JR)=FT(K,JR)*NBLADE/2.0
      STEADY=(B(1,JT)+A(1,JR))*NBLADE/2.0
      BUG=0.0
      CAT=0.0
      DO 60 M=1,NBLHIN
      N=M*NBLADE
      IF(NBLADE.LE.1) GOTO 61
      BUG=-B(N-1,JT)+A(N-1,JR)
      CAT=A(N-1,JT)+B(N-1,JR)
 61   IF(N.GE.NH) GOTO 62
      BUG=BUG+B(N+1,JT)+A(N+1,JR)
      CAT=CAT-A(N+1,JT)+B(N+1,JR)
 62   VAMP(M)=SQRT(BUG**2+CAT**2)*NBLADE/2.0
      VPH(M)=ATAN2(BUG,CAT)*57.296*(-1)+90.
 60   IF(VPH(M).LT.0.0) VPH(M)=VPH(M)+360.
      M=0
      WRITE(16,101) M,STEADY
      DO 70 M=1,NBLHIN
 70   WRITE(16,101) M,VAMP(M),VPH(M)
 101  FORMAT(33X,I2,F11.5,F8.1)
      RETURN
      END








