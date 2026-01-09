      SUBROUTINE DPDVEXTRAP(M,IDR)
************************************************************************
*   cubic (or lower) extrapolation of dpdv to the left split panel     *
*   moved here from delr.f on 06-23-92 NF                              *
*   moved here from PSFCAV on 06-26-92 NF                              *
*   02-19-97 NF Changed VO[XYZ] to VO[XYZ]1                            *
*                                                                      *
*   DATE             CORRECTIONS/COMMENTS                              *
*   -------------    ------------------------------------------------  *
*   11-27-98 JY      Made several corrections relating to face         *
*                    cavitation (see CfaceJY).                         *
*   02-06-99 JY      modified subroutine to allow cavity to grow on    *
*                    both the back and face of the foil.               *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
************************************************************************
      INCLUDE 'PUFCAV.INC'

      Q1=ZERO
      Q2=ZERO
      Q3=ZERO
      Q4=ZERO

      IF(IDR.EQ.1) THEN
         K=1
         ISF=0
      ELSE IF(IDR.EQ.2) THEN
         K=-1
         ISF=1
      END IF
      AK1=FLOAT(K)

      M0M=M0(M,IDR)
      M0M1=M0M-1
      JCAV=JCV(M,IDR)
      LCAV=M0M+K*(JCAV+1)
      N1=M0M1+K*JCAV+ISF
      NN1=N1+K*1

      DP1=DPHIDV(JCAV,M,IDR)
      IF(JCAV-1.GT.0) DP2=DPHIDV(JCAV-1,M,IDR)
      IF(JCAV-2.GT.0) DP3=DPHIDV(JCAV-2,M,IDR)
      IF(JCAV-3.GT.0) DP4=DPHIDV(JCAV-3,M,IDR)

c      IF(JCAV.GE.4)THEN
C.......4-point extrapolation...........................................
c         S3=AK1*(SPZ(N1-K*2,M)-SPZ(N1-K*3,M))
c         S2=AK1*(SPZ(N1-K*1,M)-SPZ(N1-K*3,M))
c         S1=AK1*(SPZ(N1,M)-SPZ(N1-K*3,M))
c         S0=AK1*(SPZ(NN1,M)-SPZ(N1-K*3,M))
c         CALL CUBEXT(S0,S1,S2,S3,Q1,Q2,Q3,Q4)
c      ELSEIF(JCAV.EQ.3)THEN
C.......3-point extrapolation...........................................
c         S2=AK1*(SPZ(N1-K*1,M)-SPZ(N1-K*2,M))
c         S1=AK1*(SPZ(N1,M)-SPZ(N1-K*2,M))
c         S0=AK1*(SPZ(NN1,M)-SPZ(N1-K*2,M))
c         CALL QUADEXT(S0,S1,S2,Q1,Q2,Q3)

c         Q4=ZERO
c      ELSEIF(JCAV.EQ.2)THEN
      IF(JCAV.GE.2)THEN
C.......2-point extrapolation...........................................
         S1=AK1*(SPZ(N1,M)-SPZ(N1-K*1,M))
         S0=AK1*(SPZ(NN1,M)-SPZ(N1-K*1,M))
         SL=S0+.25*DZL(M,IDR)
         Q1=S0/S1
         Q2=ONE-S0/S1

         Q3=ZERO
         Q4=ZERO
      ELSEIF(JCAV.EQ.1)THEN
C.......1-point extrapolation...........................................
         Q1=ONE
         
         Q2=ZERO
         Q3=ZERO
         Q4=ZERO
      ENDIF
      DPDVSP(M,IDR)=DP1*Q1+DP2*Q2+DP3*Q3+DP4*Q4
      
      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine DPDVEXTRAP>>>>>>>>>>>>>>>>>>>>>>
      END




