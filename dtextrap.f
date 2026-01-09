      SUBROUTINE DTEXTRAP(M,IDR)
************************************************************************
*   right split-panel dipole uses a cubic, quadratic or linear         *
*   extrapolation to define phi_r                                      *
*   moved here from delr.f on 06-23-92 NF                              *
*   moved here from PSFCAV on 06-26-92 NF                              *
*                                                                      *
*   02-06-99 JY      modified subroutine to allow cavity to grow on    *
*                    both the back and face of the foil.               *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
************************************************************************
      INCLUDE 'PUFCAV.INC'

      IF(IDR.EQ.2)THEN
         N1=M0(M,IDR)-JCV(M,IDR)-2
         IF(ISC.EQ.1) N1=N1-N0(2)+1
         LCAV=M0(M,IDR)-JCV(M,IDR)-2
         I=-1
         SGN=-1.0
      ELSE IF(IDR.EQ.1) THEN
         LCAV=LCV(M,IDR)
         IF(ISC.EQ.0) THEN
            N1=NCP-LCAV
         ELSE
            N1=N0(1)-LCAV
         END IF
         I=1
         SGN=1.0
      ENDIF

C      IF(N1.GE.4)THEN
C.......4-point cubic extrapolation
C         S3=(SPZ(LCAV+(I*3),M)-SPZ(LCAV+(I*2),M))*SGN
C         S2=(SPZ(LCAV+(I*3),M)-SPZ(LCAV+(I*1),M))*SGN
C         S1=(SPZ(LCAV+(I*3),M)-SPZ(LCAV,M))*SGN
C         S0=S1+HALF*(DS(LCAV,M)+DZR(M,IDR))
C         CALL CUBEXT(S0,S1,S2,S3,DT(M,1,IDR),DT(M,2,IDR),
C     *        DT(M,3,IDR),DT(M,4,IDR))
C      ELSEIF(N1.EQ.3)THEN
C.......3-point quadratic extrapolation
C         S2=(SPZ(LCAV+(I*2),M)-SPZ(LCAV+(I*1),M))*SGN
C         S1=(SPZ(LCAV+(I*2),M)-SPZ(LCAV,M))*SGN
C         S0=S1+HALF*(DS(LCAV,M)+DZR(M,IDR))
C         CALL QUADEXT(S0,S1,S2,DT(M,1,IDR),DT(M,2,IDR),
C     *        DT(M,3,IDR))
C         
C         DT(M,4,IDR)=ZERO
C      ELSEIF(N1.EQ.2)THEN
      IF(N1.GE.2) THEN
C.......2-point linear extrapolation
         S1=(SPZ(LCAV+(I*1),M)-SPZ(LCAV,M))*SGN
         S0=S1+HALF*(DS(LCAV,M)+DZR(M,IDR))
         DT(M,1,IDR)=S0/S1
         DT(M,2,IDR)=ONE-S0/S1

         DT(M,3,IDR)=ZERO
         DT(M,4,IDR)=ZERO
      ELSE
         DT(M,1,IDR)=ONE

         DT(M,2,IDR)=ZERO
         DT(M,3,IDR)=ZERO
         DT(M,4,IDR)=ZERO
      ENDIF

      RETURN
C<<<<<<<<<<<<<<<<<<<<<<<<<end of subroutine DTEXTRAP>>>>>>>>>>>>>>>>>>>>
      END

