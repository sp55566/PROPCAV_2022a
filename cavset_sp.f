      SUBROUTINE CAVSET_SP
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
*  JY010702               Copied from cavset_sc.f.  Modified for ISP=1.*
*                                                                      *
************************************************************************
      USE MEMSOL
      USE CVRHS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C      COMMON/MEMSOL/ALHS(NTZ,NTZ)
C      COMMON/CVRHS/RHS(NTZ)

      IEQN=0
      IF(NBW.EQ.0) GO TO 15

C-----------------------------------------------------------------------
C     Green's formula on the blade
C-----------------------------------------------------------------------
      DO 10 MM=MR,1,-1
         DO 20 NN=1,NC
            I=INDEXB(NN,MM)

C..........only solve for submerged panels..............................
            IF(ISUBM(I,IDXREV).EQ.1) THEN
               IEQN=IEQN+1

               CALL CAVSC1_SP(IEQN,I)

            END IF

 20      CONTINUE

         IEQN=IEQN+NNWC(MM)

 10   CONTINUE

 15   IF(IHUB.EQ.0) GOTO 30
C-----------------------------------------------------------------------
C     Green's formula on the hub
C-----------------------------------------------------------------------
      IEQN=NBW
      DO 40 NN=1,NHBX
         DO 50 MM=1,MHBT
             I=INDEXH(NN,MM)

C..........only solve for submerged panels..............................
             IF(ISUBM(I,IDXREV).EQ.1) THEN
                IEQN=IEQN+1
                
                CALL CAVSC1_SP(IEQN,I)

             END IF

 50      CONTINUE
 40   CONTINUE

 30   CONTINUE

      IF(NBW.EQ.0) GO TO 95
C-----------------------------------------------------------------------
C     Green's formula in the wake
C-----------------------------------------------------------------------
      IEQN=0
      DO 90 MM=MR,1,-1
         IEQN=IEQN+NPERM(MM)-NNWC(MM)
         DO 100 NN=1,NNWC(MM)
            I=(MR-MM)*NTRA+NN
            
            IEQN=IEQN+1

            CALL CAVSC2_SP(MM,NN,IEQN,I)

 100     CONTINUE
 90   CONTINUE

C-----------------------------------------------------------------------
C     Put the unsteady wake and the other blade influence in RHS
C-----------------------------------------------------------------------
 95   IF(NTSTEP.GT.1) CALL CAVRHS_SP

      RETURN
      END
