       SUBROUTINE CAVSC1_SP(IEQN,I)
************************************************************************
*                                                                      *
*  SET up the linear system of equations for solving the unsteady      *
*  cavitating propeller problem for ISC=1. Originally was setup for the*
*  wing problem. CAVSET computes the entire LHS (which changes         *
*  with every timestep because cavity size changes) and the part of    *
*  the RHS which is associated with the cavity on the key blade. The   *
*  remainder of the RHS, namely the influence of the other blades and  *
*  the wakes of all the blades, is computed in CAVRHS.                 *
*                                                                      *
*  This subroutine is called inside CAVSETSC_SP.  This sets up the eqns*
*  for equation IEQN for the given column I.  This part is the         *
*  for Green's formula on the blade, the hub, the bulb, and the        *
*  tip vortex.                                                         *
*                                                                      *
*     ALHS(I,J) = array which contains the left hand side of the matrix*
*                equation                                              *
*     RHS(I) = array which contains the right hand side of the matrix  *
*                equation                                              *
*                                                                      *
************************************************************************
      USE MEMSOL 
      USE CVRHS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C      COMMON/MEMSOL/ALHS(NTZ,NTZ)
C      COMMON/CVRHS/RHS(NTZ)
C      DIMENSION WW(NPANZ,NPAWZ)

CVV
      ALLOCATABLE :: WW(:,:)  
CVV

CVV
      ALLOCATE(WW(NPANZ,NPAWZ))
CVV

      IF(NTSTEP.EQ.1.AND.NBW.NE.0) THEN
         REWIND 81
         IF(IMG.EQ.1) REWIND 181
         DO L=1,NWMINFW
            DO M=MR,1,-1
               J=IDXWAK(L,M)
               CALL READ1(81,TEMP1,NPANEL)
               IF(IMG.EQ.1) CALL READ1(181,TEMP2,NPANEL)
               DO I1=1,NPANEL
                  IF(L.LE.MSW(M,IDXREV)) THEN
                     IF(IMG.EQ.0) THEN
                        WW(I1,J)=TEMP1(I1)
                     ELSE
                        WW(I1,J)=TEMP1(I1)-TEMP2(I1)
                     END IF
                  ELSE
                     WW(I1,J)=ZERO
                  END IF
               END DO
            END DO
         END DO
      END IF

C....IFRST1 marks the beginning of the unknowns on the foil strip.
C....within a column of the matrix................................
      IFRST1=0

      RHS(IEQN)=ZERO

      IF(NBW.EQ.0) GO TO 35

      DO 30 M=MR,1,-1
         IFRST=INDEXB(0,M)
         JSUBM=0

         JCAVB=JCV(M,1)
         M0M=M0(M,1)
         M0M1=M0M-1

         FCTF=ZERO
         FCTB=ZERO
         IF(ISC.EQ.1) THEN
            FCTF2=ZERO
            FCTB2=ZERO
         END IF

         IF(ISC.EQ.0) THEN
            NFRST=IW(1,M,IDXREV)
            NLAST=IC(2,M,IDXREV)
         ELSE
            NFRST=MAX(N0(2),IW(1,M,IDXREV))
            NLAST=MIN(IC(2,M,IDXREV),N0(1)-1)
         END IF

C.......Submerged panels on the face side 
         IF(ICB(M,2,IDXREV).EQ.1) THEN

            IF(ISC.EQ.0) THEN
C.............wettted panels on the face
               DO 40 N=IW(1,M,IDXREV),IW(2,M,IDXREV)
                  JSUBM=JSUBM+1
                  ALHS(IEQN,IFRST1+JSUBM)=AA(I,IFRST+N)
 40            CONTINUE
            ELSE

               IF(IW(1,M,IDXREV).LT.N0(2)) THEN

C................cavitating panels on face SR region
                  DO N=IW(1,M,IDXREV),N0(2)-1
                     JSUBM=JSUBM+1
                     ALHS(IEQN,IFRST1+JSUBM)=-BB(I,IFRST+N)
                     FCTF2=FCTF2+AA(I,IFRST+N)
                  END DO

C................term associated with Kutta condition
                  IF(IW(1,M,IDXREV).EQ.1) THEN
                     IF(NTSTEP.EQ.1) THEN
                        DO L=1,MSW(M,IDXREV)
                           FCTF2=FCTF2-WW(I,IDXWAK(L,M))
                        END DO
                     ELSE
                        FCTF2=FCTF2-(HALF*W(I,M)+WSUBIF(I,M))
                     END IF
                  END IF
               END IF

C.............wetted panels on the face side
               DO N=NFRST,IW(2,M,IDXREV)
                  JSUBM=JSUBM+1
                  ALHS(IEQN,IFRST1+JSUBM)=AA(I,IFRST+N)
               END DO

C.............only if wetted or partially cavitating
               IWLEF2=0
               IF(IW(2,M,IDXREV).GE.N0(2).
     *              AND.IW(1,M,IDXREV).LT.N0(2)) THEN     
                  IWLEF2=IW(2,M,IDXREV)-N0(2)+1

                  N1=IFRST1+N0(2)-IW(1,M,IDXREV)+1
                  N2=N1+1
                  N3=N1+2
                  
                  IF(IWLEF2.GE.1) ALHS(IEQN,N1)=
     *                 ALHS(IEQN,N1)+FCTF2*CT2(M,1,2)
                  IF(IWLEF2.GE.2) ALHS(IEQN,N2)=
     *                 ALHS(IEQN,N2)+FCTF2*CT2(M,2,2)
                  IF(IWLEF2.GE.3) ALHS(IEQN,N3)=
     *                 ALHS(IEQN,N3)+FCTF2*CT2(M,3,2)
                  
               END IF

            END IF

         END IF

C.......submerged panels on the back side
         IF(ICB(M,1,IDXREV).EQ.1) THEN

            IF(ISC.EQ.1.AND.IC(2,M,IDXREV).GE.N0(1)) THEN

C............cavitating panels on back SR region
               NW1=IW(2,M,IDXREV)-IW(1,M,IDXREV)+1
               NJ=NW1+N0(1)-IC(1,M,IDXREV)
               DO N=N0(1),IC(2,M,IDXREV)                  
                  NJ=NJ+1
                  ALHS(IEQN,IFRST1+NJ)=-BB(I,IFRST+N)
                  FCTB2=FCTB2+AA(I,IFRST+N)
               END DO

C.............term associated with Kutta condition
               IF(IC(2,M,IDXREV).EQ.NC) THEN
                  IF(NTSTEP.EQ.1) THEN
                     DO L=1,MSW(M,IDXREV)
                        FCTB2=FCTB2+WW(I,IDXWAK(L,M))
                     END DO
                  ELSE
                     FCTB2=FCTB2+(HALF*W(I,M)+WSUBIF(I,M))
                  END IF
               END IF
            END IF

            IF(JCAVB.LE.0) THEN               

C.............If all the back blade panels are wetted
               DO N=IC(1,M,IDXREV),NLAST
                  JSUBM=JSUBM+1
                  ALHS(IEQN,IFRST1+JSUBM)=AA(I,IFRST+N)
               END DO

            ELSE

C.............Wetted panels before cavity on the back side
               DO 45 N=IC(1,M,IDXREV),M0M1
                  JSUBM=JSUBM+1
                  ALHS(IEQN,IFRST1+JSUBM)=AA(I,IFRST+N)
 45            CONTINUE

C.............Cavitating panels on the back side
               DO 50 N=M0M,NLAST
                  JSUBM=JSUBM+1
                  ALHS(IEQN,IFRST1+JSUBM)=-BB(I,IFRST+N)
                  FCTB=FCTB+AA(I,IFRST+N)
 50            CONTINUE

C.............Keep straight with the index for SR panels
               IF(ISC.EQ.1.AND.IC(2,M,IDXREV).GE.N0(1)) THEN
                  JSUBM=JSUBM+IC(2,M,IDXREV)-N0(1)+1
               END IF

C.............term associated with Kutta condition
               IF(ISC.EQ.0.AND.IC(2,M,IDXREV).EQ.NC.
     *              AND.SOP(M,1).EQ.ONE) THEN
                 IF(NTSTEP.EQ.1) THEN
                     DO L=1,MSW(M,IDXREV)
                        FCTB=FCTB+WW(I,IDXWAK(L,M))
                     END DO
                  ELSE
                     FCTB=FCTB+(HALF*W(I,M)+WSUBIF(I,M))
                  END IF
               END IF

               IF(ISC.EQ.1.AND.SOP(M,1).EQ.ONE) FCTB=FCTB+FCTB2

C.............terms associated with phi0 extrapolation 
C.............of the cavity
               IF(IC(1,M,IDXREV)-1.EQ.IW(2,M,IDXREV)) THEN
                  IWLEB=NLEP(M,IDXREV,1)+
     *                 IW(2,M,IDXREV)-IW(1,M,IDXREV)+1
                  N1=IFRST1+IWLEB
               ELSE
                  IWLEB=NLEP(M,IDXREV,1)-(IC(1,M,IDXREV)-NHP)
                  N1=IFRST1+IW(2,M,IDXREV)-IW(1,M,IDXREV)+1+IWLEB
               END IF

               IF(IWLEB.GE.1) ALHS(IEQN,N1)=
     *              ALHS(IEQN,N1)+FCTB*CT(M,1,1)
               IF(IWLEB.GE.2) ALHS(IEQN,N1-1)=
     *              ALHS(IEQN,N1-1)+FCTB*CT(M,2,1)
               IF(IWLEB.GE.3) ALHS(IEQN,N1-2)=
     *              ALHS(IEQN,N1-2)+FCTB*CT(M,3,1)
               
            END IF

C..........only if wetted or partially cavitating
            IWLEB2=0
            IF(SOP(M,1).EQ.ZERO.AND.ISC.EQ.1.AND.
     *           IC(1,M,IDXREV).LT.N0(1).AND.
     *           IC(2,M,IDXREV).GE.N0(1)) THEN
               IWLEB2=N0(1)-IC(1,M,IDXREV)

               N1=IFRST1+IW(2,M,IDXREV)-IW(1,M,IDXREV)+1+IWLEB2
               N2=N1-1
               N3=N1-2                  
               IF(IWLEB2.GE.1) ALHS(IEQN,N1)=
     *              ALHS(IEQN,N1)+FCTB2*CT2(M,1,1)
               IF(IWLEB2.GE.2) ALHS(IEQN,N2)=
     *              ALHS(IEQN,N2)+FCTB2*CT2(M,2,1)
               IF(IWLEB2.GE.3) ALHS(IEQN,N3)=
     *              ALHS(IEQN,N3)+FCTB2*CT2(M,3,1)
            END IF

         END IF

C.......term associated with Kutta condition
         IF(ISC.EQ.0) THEN
            N1=IFRST1+1
            N2=IFRST1+IW(2,M,IDXREV)-IW(1,M,IDXREV)+
     *           IC(2,M,IDXREV)-IC(1,M,IDXREV)+2
            IF(NTSTEP.EQ.1) THEN
               DUM1=ZERO
               DO L=1,MSW(M,IDXREV)
                  DUM1=DUM1+WW(I,IDXWAK(L,M))
               END DO
            ELSE
               DUM1=HALF*W(I,M)+WSUBIF(I,M)
            END IF
            
            IF(ICB(M,2,IDXREV).EQ.1.AND.IW(1,M,IDXREV).EQ.1) THEN
               ALHS(IEQN,N1)=ALHS(IEQN,N1)-DUM1
            END IF
            IF(ICB(M,1,IDXREV).EQ.1.AND.IC(2,M,IDXREV).EQ.NC.
     *           AND.SOP(M,1).EQ.ZERO) THEN
               ALHS(IEQN,N2)=ALHS(IEQN,N2)+DUM1
            END IF
         END IF

C.......submerged, supercavitating wake panels
         DO 60 N=1,NNWC(M)
            J=(MR-M)*NTRA+N
            JSUBM=JSUBM+1
            ALHS(IEQN,IFRST1+JSUBM)=-C(I,J)
 60      CONTINUE

C.......compute the right-hand-side due to the key blade
         SUM2=ZERO
         SUM3=ZERO
         TERM1=ZERO
         TERM2=ZERO
         TERM3=ZERO
         
C.......submerged panels on the face side
         IF(ICB(M,2,IDXREV).EQ.1) THEN

            IF(ISC.EQ.0) THEN
C.............wettted panels on the face
               DO 70 N=IW(1,M,IDXREV),IW(2,M,IDXREV)
                  N1=IFRST+N
                  SUM3=SUM3+BB(I,N1)*DPDNC(N1)
 70            CONTINUE
            ELSE

               IF(IW(1,M,IDXREV).LT.N0(2)) THEN

C................cavitating panels on face SR region
                  DO N=IW(1,M,IDXREV),N0(2)-1
                     N123=N0(2)-N
                     SUM2=SUM2+AA(I,IFRST+N)*PHI2(N123,M,2)
                  END DO

C................only if wetted or partially cavitating
                  TERM3=TERM3-CT2(M,4,2)*QCSR(1,M,2)*FCTF2
               END IF

C.............wettted panels on the face
               DO N=NFRST,IW(2,M,IDXREV)
                  N1=IFRST+N
                  SUM3=SUM3+BB(I,N1)*DPDNC(N1)
               END DO

            END IF
         END IF

C.......submerged panels on the back side
         IF(ICB(M,1,IDXREV).EQ.1) THEN
            
            IF(ISC.EQ.1.AND.IC(2,M,IDXREV).GE.N0(1)) THEN

C.............cavitating panels on back SR region
               DO N=N0(1),IC(2,M,IDXREV)
                  N123=N-N0(1)+1
                  SUM2=SUM2+AA(I,IFRST+N)*PHI2(N123,M,1)
               END DO

               IF(SOP(M,1).EQ.ZERO) THEN
C................only if wetted or partially cavitating
                  TERM3=TERM3-CT2(M,4,1)*QCSR(1,M,1)*FCTB2
               ELSE
C................only if supercavitating
                  TERM3=TERM3-PSI0T(M,1)*FCTB2
               END IF

            END IF

            IF(JCAVB.LE.0) THEN               

C.............If all the back blade panels are wetted
               DO N=IC(1,M,IDXREV),NLAST
                  N1=IFRST+N
                  SUM3=SUM3+BB(I,N1)*DPDNC(N1)
               END DO

            ELSE

C.............Wetted panels before cavity on the back side
               DO 75 N=IC(1,M,IDXREV),M0M1
                  N1=IFRST+N
                  SUM3=SUM3+BB(I,N1)*DPDNC(N1)
 75            CONTINUE
            
C.............Cavitating panels on the back side
               ICAV=0
               DO 80 N=M0M,NLAST
                  ICAV=ICAV+1
                  SUM2=SUM2+AA(I,IFRST+N)*PHI1(ICAV,M,1)
 80            CONTINUE

C.............term associated with phi0 extrapolation
               IF(IWLEB.GT.0) TERM3=TERM3-
     *              CT(M,4,1)*QC(1,M,1)*FCTB

            END IF

         END IF

C.......sum the sources and dipoles on the non-split panels
         TERM1=-SUM2+SUM3
                  
C.......term associated with Kutta condition...................
         IF(NTSTEP.EQ.1) THEN
            DUM1=ZERO
            DO L=1,MSW(M,IDXREV)
               DUM1=DUM1+WW(I,IDXWAK(L,M))
            END DO
         ELSE
            DUM1=HALF*W(I,M)+WSUBIF(I,M)
         END IF
         
         IF(ISC.EQ.0) THEN
            IF(ICB(M,1,IDXREV).EQ.1.AND.IC(2,M,IDXREV).EQ.NC.
     *           AND.SOP(M,1).EQ.ONE) THEN
               TERM2=TERM2-DUM1*PHI1(ICAV,M,1)
            END IF
         ELSE
            IF(ICB(M,2,IDXREV).EQ.1.AND.IW(1,M,IDXREV).EQ.1) THEN
               TERM2=TERM2+DUM1*PHI2(NSR2,M,2)
            END IF
            IF(ICB(M,1,IDXREV).EQ.1.AND.IC(2,M,IDXREV).EQ.NC) THEN
               TERM2=TERM2-DUM1*PHI2(NSR2,M,1)
            END IF
         END IF
                  
         RHS(IEQN)=RHS(IEQN)+TERM1+TERM2+TERM3

         IFRST1=IFRST1+NPERM(M)
         
 30   CONTINUE

 35   IF(IHUB.NE.0)THEN
         J1=0

C.......cal. terms to LHS that are associated with the hub...........
         DO N=1,NHBX
            DO M=1,MHBT
               J=INDEXH(N,M)

C.............Submerged panels on the hub
               IF(ISUBM(J,IDXREV).EQ.1) THEN
                  J1=J1+1
                  ALHS(IEQN,NBW+J1)=AA(I,J)
               END IF

            END DO
                  
C..........cal. terms to RHS that are associated with the hub.....
            SUM1=ZERO
            DO M=1,MHBT
               J=INDEXH(N,M)

C.............Submerged panels on the hub
               IF(ISUBM(J,IDXREV).EQ.1) THEN
                  SUM1=SUM1+BB(I,J)*DPDNC(J)
               END IF

            END DO
            RHS(IEQN)=RHS(IEQN)+SUM1
         END DO
      ENDIF

CVV
      DEALLOCATE(WW)
CVV
      RETURN
      END
