      SUBROUTINE CAVSET_SC
************************************************************************
*                                                                      *
*  SET up the linear system of equations for solving the unsteady      *
*  cavitating propeller problem. Originally was setupb for the wing    *
*  problem. CAVSET computes the entire left-hand-side (which changes   *
*  with every timestep because cavity size changes) and the part of    *
*  the RHS which is associated with the cavity on the key blade. The   *
*  remainder of the RHS, namely the influence of the other blades and  *
*  the wakes of all the blades, is computed in CAVRHS.                 *
*                                                                      *
*  Notes:  The l.h. side of the matrix should contain ((NC)xMR)^2      *
*          elements.  The r.h. side will contain (NCxMR).  Plus any    *
*          wake panels and minus any split panels. The lhs and rhs are *
*          passed into  blic2 so that a solution may be found.  CM     *
*                                                                      *
*  Author: Neal Fine   10-07-91                                        *
*                                                                      *
*     ALHS(I,J) = array which contains the left hand side of the matrix*
*                equation                                              *
*     RHS(I) = array which contains the right hand side of the matrix  *
*                equation                                              *
*  Date of last revision                      Revision                 *
* -----------------------                  --------------              *
*  JY091901               Copied from cavset.f.  Modified for ISC=1.   *
*                                                                      *
************************************************************************
      USE MEMSOL
      USE CVRHS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C      COMMON/MEMSOL/ALHS(NTZ,NTZ)
C      COMMON/CVRHS/RHS(NTZ)

      DIMENSION ISPE(2)

C-----------------------------------------------------------------------
C     Green's formula on the blade
C-----------------------------------------------------------------------
      IEQN=0

      DO 10 MM=MR,1,-1
         ISPE(1)=999
         ISPE(2)=999

         NSPSUM=0
         DO 20 IDR=1,2
            IF(IDR.EQ.1) THEN
               IF(NSPP(MM,IDR).NE.0) ISPE(IDR)=LCV(MM,IDR)-1
            ELSE IF(IDR.EQ.2) THEN
               IF(NSPP(MM,IDR).NE.0) ISPE(IDR)=M0(MM,IDR)-1-JCV(MM,IDR)
            END IF
            NSPSUM=NSPSUM+NSPP(MM,IDR)
 20      CONTINUE

         DO 30 NN=1,NC-NSPSUM

C..........set indices for equation number and influence coefficients...
C..........IEQN is the equation number (row index)......................
C..........I1 is the index of the control point where GREEN is satisfied
C..........IFRST1 marks the beginning of the unknowns on the foil strip.
C..........within a column of the matrix................................
            IEQN=IEQN+1
            I=INDEXB(NN,MM)
            I1=I

            IF(NSPSUM.EQ.1) THEN
               IF(NN.GE.ISPE(1).OR.NN.GE.ISPE(2)) I1=I1+1
            ELSE IF(NSPSUM.EQ.2) THEN
               IF(NN.GE.ISPE(2)) I1=I1+1
               IF(NN+1.GE.ISPE(1)) I1=I1+1
            END IF

            CALL CAVSC1(IEQN,I1)
            
 30      CONTINUE
         
         IEQN=IEQN+NNWC(MM)
         
 10   CONTINUE
      
      IF(IHUB.EQ.0) GOTO 270
C-----------------------------------------------------------------------
C     Green's formula on the hub
C-----------------------------------------------------------------------
      IEQN=NBW
      DO 280 NN=1,NHBX
         DO 290 MM=1,MHBT

C..........set indices for equation number and influence coefficients...
C..........IEQN is the equation number (row index)......................
            IEQN=IEQN+1

C..........I1 is the index of the control point where GREEN is .........
C..........satisfied....................................................
            I1=INDEXH(NN,MM)

            CALL CAVSC1(IEQN,I1)

 290     CONTINUE
 280  CONTINUE
 270  CONTINUE

cC-----------------------------------------------------------------------
cC     Green's formula on the bulb and on the tip vortex cavity
cC-----------------------------------------------------------------------
c       IF(IAN.NE.2) GO TO 1200
c
c       DO NN=1,NTHX
c          DO MM=1,MCVT
c
cC..........set indices for equation number and influence coefficients...
cC..........IEQN is the equation number (row index)......................
c             IEQN=IEQN+1
c             
cC..........I1 is the index of the control point where GREEN is .........
cC..........satisfied....................................................
c             I1=INDEXT(NN,MM)
c
c             CALL CAVSC1(IEQN,I1)
c            
c          END DO
c       END DO
c
c       DO NN=1,NCVX
c          DO MM=1,MCVT
c
cC..........set indices for equation number and influence coefficients...
cC..........IEQN is the equation number (row index)......................
c             IEQN=IEQN+1
c
cC..........I1 is the index of the control point where GREEN is .........
cC..........satisfied....................................................
c             I1=INDEXC(NN,MM)
c
c             CALL CAVSC1(IEQN,I1)
c
c          END DO
c       END DO
c
c 1200  CONTINUE

C-----------------------------------------------------------------------
C     Green's formula in the wake
C-----------------------------------------------------------------------
       IEQN=0
       DO 520 MM=MR,1,-1
          IEQN=IEQN+NC-NSPP(MM,2)-NSPP(MM,1)

          DO 530 NN=1,NNWC(MM)

C...........I1 is the index of the control point where GREEN is .......
C...........satisfied..................................................
             I1=(MR-MM)*NTRA+NN

             IEQN=IEQN+1
             
             CALL CAVSC2(MM,NN,IEQN,I1)

 530     CONTINUE
 520  CONTINUE

C-----------------------------------------------------------------------
C     Put the unsteady wake and the other blade influence in RHS
C-----------------------------------------------------------------------
      CALL CAVRHS

      RETURN
C<<<<<<<<<<<<<<<<<<<<end of Subroutine cavset>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END




