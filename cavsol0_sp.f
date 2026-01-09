       SUBROUTINE CAVSOL0_SP
************************************************************************
*      This subroutine solves for the submerged hub panels (if any)    *
*      for the case of ISP=1 with NBW=0.                               *
*                                                                      *
*      Author: Julie Young  01-19-02                                   *
*                                                                      *
*      Date of last Revision         Revision                          *
*      ---------------------         --------                          *
************************************************************************
       USE MEMSOL
       USE CVRHS
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

C       COMMON/MEMSOL/ALHS(NTZ,NTZ)
C       COMMON/CVRHS/RHS(NTZ)
C       DIMENSION TMP1(NPANZ),TMP2(NPANZ)

CVV
      ALLOCATABLE :: TMP1(:),TMP2(:)
CVV


CVV
      ALLOCATE(TMP1(NPANZ),TMP2(NPANZ))
CVV

C.....Set all the source and dipole strengths in the blade and wake to 
C.....be zero since NBW=0 (all these panels are dry).
       DO M=1,MR
          DO N=1,NC
             J=INDEXB(N,M)
             POT(J)=ZERO
             DPDNC(J)=ZERO
             CPB(N,M)=ZERO
          END DO
          DO N=1,NTRA
             J=(MR-M)*NTRA+N
             POTW(J)=ZERO
             SORW(J)=ZERO
          END DO
          DO IDR=1,2
             DO N=1,NHP+NTRA
                HT(N,M,IDR)=ZERO
             END DO
             DELTA(M,IDR)=ZERO
          END DO
          IF(ISC.EQ.1) DELTAT1(M)=ZERO
       END DO

C.....Solve for submerged hub panels
       IF(IHUB.NE.0) THEN

C........Count the number of submerged hub panels
          NORD=0
          NBLOCK=0 

C........Solve for the submerged hub panels associated with all the blades
          DO KK=1,NBLADE
             IREC=NTPOS(KK)

C...........Count the total number of submerged panels
             DO N=1,NHBX
                IF(ISUBH(N,IREC).GT.0) THEN
                   NBLOCK=NBLOCK+1
                   NPERB(NBLOCK)=ISUBH(N,IREC)
                   NORD=NORD+ISUBH(N,IREC)
                END IF
             END DO
          END DO

          IF(NORD.GT.0) THEN

             REWIND 41
             IF(IMG.EQ.1) REWIND 141
             REWIND 42
             IF(IMG.EQ.1) REWIND 142
             
             IEQN=0

C...........source influence
             DO KK=1,NBLADE
                IREC=NTPOS(KK)
                CALL READ2(47,IREC,STRGTH(1,KK),NPANEL)
             END DO

C-----------------------------------------------------------------------
C     Set up system of equations for the submerged hub panels.
C-----------------------------------------------------------------------
             DO I=1,NPANEL
                DO KK=1,NBLADE

C.................Det. angular position of hub panels for all blades
                   IREC=NTPOS(KK)

C.................Read dipole (41) & source (42) influence coef.              
                   CALL READ1(41,TEMP1,NPANEL)
                   CALL READ1(42,TEMP2,NPANEL)
       
                   IF(IMG.EQ.1) THEN
                      CALL READ1(141,TMP1,NPANEL)
                      CALL READ1(142,TMP2,NPANEL)
                   END IF
                   IF(I.GT.NPANB.AND.ISUBM(I,IREC).EQ.1) THEN
                      IEQN=IEQN+1
                      J1=0
                      RHS(IEQN)=ZERO
                      
                      DO J=1,NPANEL
                         IF(J.GT.NPANB.AND.ISUBM(J,IREC).EQ.1) THEN
                            J1=J1+1

                            IF(IMG.EQ.1) THEN
                               ALHS(IEQN,J1)=TEMP1(J)-TMP1(J)
                               RHS(IEQN)=RHS(IEQN)+
     *                              (TEMP2(J)-TMP2(J))*STRGTH(J,KK)
                            ELSE
                               ALHS(IEQN,J1)=TEMP1(J)
                               RHS(IEQN)=RHS(IEQN)+
     *                              TEMP2(J)*STRGTH(J,KK)
                            END IF
                         END IF
                      END DO

                   END IF

                END DO

             END DO

       
C-----------------------------------------------------------------------
C     Solve the system of equations for submerged hub panels.
C-----------------------------------------------------------------------
C...........Note: If ITMAX>100, You MUST set ITMX=ITMAX+1 in BLIC2.F
             ITMAX=100
             IPASS=1

C...........Tolerance of the matrix solution is set to be 0.000005
             TOL=5.0E-06

             CALL BLIC2(SOL,RHS,TOL,NORD,NBLOCK,NPERB,ITMAX,IPASS,
     &                  NTZ,NBLKMAX,101)

C...........Obtain potential for the hub from the solution
             J1=0
             DO I1=1,NPANH
                I=NPANB+I1
                DO KK=1,NBLADE
                   IREC=NTPOS(KK)
                   IF(ISUBM(I,IREC).EQ.1) THEN
                      J1=J1+1
                      POT(I)=SOL(J1)
                   ELSE

C....................Set all the dry hub panel singularity 
C....................strengths to zero
                      POT(I)=ZERO
                      DPDNC(I)=ZERO
                   END IF
                END DO
             END DO
          
          ELSE
             
C...........Set all the dry hub panel singularity strengths to zero
             DO NN=1,NHBX               
                DO MM=1,MHBT
                   J=INDEXH(NN,MM)
                   POT(J)=ZERO
                   DPDNC(J)=ZERO
                END DO
             END DO
          END IF

       END IF
CVV
      DEALLOCATE(TMP1,TMP2)
CVV
       RETURN
       END
