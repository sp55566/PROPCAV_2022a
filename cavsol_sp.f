      SUBROUTINE CAVSOL_SP
************************************************************************
*     Notes:  The different between this subroutine and cavsol.f is    *
*             that we only consider the fully submerged panels.        *
*                                                                      *
*     Author: Julie Young  11-24-99                                    *
*                                                                      *
*     Date of last Revision         Revision                           *
*     ---------------------         --------                           *
*     JY071700  Modified subroutine to allow for detachment for ISP=1. *
************************************************************************
      USE MEMSOL
      USE CVRHS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

C      COMMON/MEMSOL/ALHS(NTZ,NTZ)
C      COMMON/CVRHS/RHS(NTZ)

C.....Note: If ITMAX>100, You MUST set ITMX=ITMAX+1 in BLIC2.F(JY051800)
      ITMAX=100
      IPASS=1
      NBLOCK=0

      DO 10 M=1,MR
         IF(NPERM(MR-M+1).GT.0) THEN
            NBLOCK=NBLOCK+1
            NPERB(NBLOCK)=NPERM(MR-M+1)
         END IF
 10   CONTINUE

      NORD=NBW

      IF(IHUB.NE.0) THEN
         DO N=1,NHBX
            IF(ISUBH(N,IDXREV).GT.0) THEN
               NBLOCK=NBLOCK+1
               NPERB(NBLOCK)=ISUBH(N,IDXREV)
               NORD=NORD+ISUBH(N,IDXREV)
            END IF
         END DO
      END IF

C.....Tolerance of the matrix solution is set to be 0.000005
      TOL=5.0E-06

C-----------------------------------------------------------------------
C     Call block iterative solver BLIC2
C-----------------------------------------------------------------------
      CALL BLIC2(SOL,RHS,TOL,NORD,NBLOCK,NPERB,ITMAX,IPASS,
     &           NTZ,NBLKMAX,101)

C-----------------------------------------------------------------------
C     dissect the results
C-----------------------------------------------------------------------
      IFRST=0

      DO 20 M=MR,1,-1

         IF(ISC.EQ.0) THEN
            NFRST=IW(1,M,IDXREV)
            NLAST=IC(2,M,IDXREV)
         ELSE
            NFRST=MAX(N0(2),IW(1,M,IDXREV))
            NLAST=MIN(IC(2,M,IDXREV),N0(1)-1)
         END IF

C.......If all panels on face side are dry.
         IF(ICB(M,2,IDXREV).EQ.0) THEN
            DO N=1,NH
               J=INDEXB(N,M)
               POT(J)=ZERO
               DPDNC(J)=ZERO
            END DO

C.......submerged panels on the face side
         ELSE

C..........dry panels on face side near TE
            IF(IW(1,M,IDXREV).GT.1) THEN
               DO N=1,IW(1,M,IDXREV)-1
                  J=INDEXB(N,M)
                  POT(J)=ZERO
                  DPDNC(J)=ZERO
               END DO

               IF(IW(1,M,IDXREV).GE.N0(2)) PHI0T(M,2)=ZERO

            END IF

C..........dry panels on face side near LE
            IF(IW(2,M,IDXREV).LT.NH) THEN
               DO N=IW(2,M,IDXREV)+1,NH
                  J=INDEXB(N,M)
                  POT(J)=ZERO
                  DPDNC(J)=ZERO
               END DO
            END IF

            IF(ISC.EQ.0) THEN
C.............wettted panels on the face
               DO 30 N=IW(1,M,IDXREV),IW(2,M,IDXREV)
                  J=INDEXB(N,M)
                  IFRST=IFRST+1
                  POT(J)=SOL(IFRST)
 30            CONTINUE
            ELSE

               IF(IW(1,M,IDXREV).LT.N0(2)) THEN

C................cavitating panels on face SR region
                  DO N=IW(1,M,IDXREV),N0(2)-1
                     J=INDEXB(N,M)
                     IFRST=IFRST+1
                     DPDNC(J)=SOL(IFRST)
                  END DO
               END IF

C.............wetted panel on face side
               DO N=NFRST,IW(2,M,IDXREV)
                  J=INDEXB(N,M)
                  IFRST=IFRST+1
                  POT(J)=SOL(IFRST)
               END DO

C.............only if wetted or partially cavitating
               IWLEF2=0
               IF(SOP(M,2).EQ.ZERO.AND.IW(2,M,IDXREV).GE.N0(2).
     *              AND.IW(1,M,IDXREV).LT.N0(2)) THEN     
                  IWLEF2=IW(2,M,IDXREV)-N0(2)+1

                  IF(IWLEF2.GT.0) THEN
                     N1=N0(2)
                     J1=INDEXB(N1,M)
                     J2=INDEXB(N1+1,M)
                     J3=INDEXB(N1+2,M)
                     PHI0T(M,2)=CT2(M,1,2)*POT(J1)+
     *                    CT2(M,2,2)*POT(J2)+CT2(M,3,2)*POT(J3)
     *                    +CT2(M,4,2)*QCSR(1,M,2)
                  ELSE
                     PHI0T(M,2)=ZERO
                  END IF

C................potential, POT=PHI0T+PHI2, beneath SR
C................surface on face
                  DO N=IW(1,M,IDXREV),N0(2)-1                 
                     N123 = N0(2)-N
                     J=INDEXB(N,M)
                     POT(J)=PHI0T(M,2)+PHI2(N123,M,2)
                  END DO
               END IF

            END IF

         END IF

C.......If all panels on back side are dry.
         IF(ICB(M,1,IDXREV).EQ.0) THEN
            DO N=NHP,NC
               J=INDEXB(N,M)
               POT(J)=ZERO
               DPDNC(J)=ZERO
            END DO

C.......submerged panels on the back side
         ELSE

C..........dry panels on back side near LE
            IF(IC(1,M,IDXREV).GT.NHP) THEN
               DO N=NHP,IC(1,M,IDXREV)-1
                  J=INDEXB(N,M)
                  POT(J)=ZERO
                  DPDNC(J)=ZERO
               END DO
            END IF

C..........dry panels on back side near TE
            IF(IC(2,M,IDXREV).LT.NC) THEN
               DO N=IC(2,M,IDXREV)+1,NC
                  J=INDEXB(N,M)
                  POT(J)=ZERO
                  DPDNC(J)=ZERO
               END DO

               IF(IC(2,M,IDXREV).LT.N0(1)) PHI0T(M,1)=ZERO
            END IF

            JCAVB=JCV(M,1)
            
            IF(JCAVB.LE.0) THEN 

C.............Wetted panels on the back
               DO N=IC(1,M,IDXREV),NLAST
                  J=INDEXB(N,M)
                  IFRST=IFRST+1
                  POT(J)=SOL(IFRST)
               END DO

            ELSE

               M0M=M0(M,1)
               M0M1=M0M-1

C.............Wetted panels on the back before the cavity
               DO 35 N=IC(1,M,IDXREV),M0M1
                  J=INDEXB(N,M)
                  IFRST=IFRST+1
                  POT(J)=SOL(IFRST)
 35            CONTINUE

C.............Cavitating panels on the back
               DO 40 N=M0M,NLAST
                  J=INDEXB(N,M)
                  IFRST=IFRST+1
                  DPDNC(J)=SOL(IFRST)
 40            CONTINUE

C.............calculate phi0 for cavity on the back
               IF(IC(1,M,IDXREV)-1.EQ.IW(2,M,IDXREV)) THEN
                  IWLEB=NLEP(M,IDXREV,1)+
     *                 IW(2,M,IDXREV)-IW(1,M,IDXREV)+1
               ELSE
                  IWLEB=NLEP(M,IDXREV,1)-(IC(1,M,IDXREV)-NHP)
               END IF
               N1=M0M1
               
               IF(IWLEB.GT.0) THEN
                  J1=INDEXB(N1,M)
                  J2=INDEXB(N1-1,M)
                  J3=INDEXB(N1-2,M)
                  PHI0(M,1)=CT(M,1,1)*POT(J1)+CT(M,2,1)*POT(J2)+
     *                 CT(M,3,1)*POT(J3)+CT(M,4,1)*QC(1,M,1)
               ELSE
                  PHI0(M,1)=ZERO
               END IF

C.............potential, POT=PHI0+PHI1, beneath cavity 
C.............surface on back
               ICAV=0
               DO 50 N=M0M,NLAST
                  ICAV=ICAV+1
                  J=INDEXB(N,M)
                  POT(J)=PHI0(M,1)+PHI1(ICAV,M,1)
 50            CONTINUE

               IF(ISC.EQ.1.AND.IC(2,M,IDXREV).GE.N0(1)) THEN

                  IF(SOP(M,1).EQ.ONE) THEN
C...................if supercavitating
                     PHI0T(M,1)=PSI0T(M,1)+PHI0(M,1)
                  ELSE

C...................only if wetted or partially cavitating
                     IWLEB2=0
                     IF(IC(1,M,IDXREV).LT.N0(1).AND.
     *                    IC(2,M,IDXREV).GE.N0(1)) IWLEB2=
     *                    N0(1)-IC(1,M,IDXREV)
                     
                     IF(IWLEB2.GT.0) THEN
                        N1=N0(1)-1
                        J1=INDEXB(N1,M)
                        J2=INDEXB(N1-1,M)
                        J3=INDEXB(N1-2,M)
                        PHI0T(M,1)=CT2(M,1,1)*POT(J1)+
     *                       CT2(M,2,1)*POT(J2)+CT2(M,3,1)*POT(J3)
     *                       +CT2(M,4,1)*QCSR(1,M,1)
                     ELSE
                        PHI0T(M,1)=ZERO
                     END IF
                  END IF

C................source and potential (POT=PHI0+PHI1) beneath SR
C................surface on back
                  DO N=N0(1),IC(2,M,IDXREV)
                     IFRST=IFRST+1
                     N123 = N-N0(1)+1
                     J=INDEXB(N,M)
                     DPDNC(J)=SOL(IFRST)
                     POT(J)=PHI0T(M,1)+PHI2(N123,M,1)
                  END DO
               END IF

            END IF

         END IF
         
C.......submerged, cavitating wake panels...............................
         DO 60 N=1,NNWC(M)
            IFRST=IFRST+1
            J=(MR-M)*NTRA+N

C..........source strength of supercavitating wake panels
            SORW(J)=SOL(IFRST)

C..........potential, POT=PHI0+PHI1, beneath supercavity surface on wake
            IF(ISC.EQ.0) THEN
               ICAV=ICAV+1
               POTW(J)=PHI0(M,1)+PHI1(ICAV,M,1)
            ELSE
               POTW(J)=PHI0T(M,1)+PHI2(NSR2+N,M,1)
            END IF
 60      CONTINUE

         IF(NNWC(M).LT.NTRA) THEN
            DO N=NNWC(M)+1,NTRA
               J=(MR-M)*NTRA+N
               POTW(J)=ZERO
               SORW(J)=ZERO
            END DO
         END IF

 20   CONTINUE

C....Obtain potential for the hub from the solution
      IF(IHUB.NE.0) THEN
         J1=0
         DO NN=1,NHBX               
            DO MM=1,MHBT
               J=INDEXH(NN,MM)
               IF(ISUBM(J,IDXREV).EQ.1) THEN
                  J1=J1+1
                  POT(J)=SOL(NBW+J1)
               END IF
            END DO
         END DO
      END IF

      RETURN
      END
