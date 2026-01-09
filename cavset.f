      SUBROUTINE CAVSET
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
*  Author: Neal Fine   10-07-91                                        *
*                                                                      *
*     ALHS(I,J) = array which contains the left hand side of the matrix*
*                equation                                              *
*     RHS(I) = array which contains the right hand side of the matrix  *
*                equation                                              *
*  Date of last revision                      Revision                 *
* -----------------------                  --------------              *
*  JYHL102401                  Subroutine undergone major revision.    *
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
      ISR=1
      IF(IFACE.EQ.2) ISR=2

      DO 10 MM=MR,1,-1
         ISPE(1)=999
         ISPE(2)=999

         NSPSUM=0
         DO 20 II=1,ISR
            IF((IFACE.EQ.0).OR.(II.EQ.1.AND.IFACE.EQ.2)) THEN
               IDR=1
            ELSE
               IDR=2
            END IF
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

            CALL CAVSETSUB(IEQN,I1)

 30      CONTINUE
     
         IEQN=IEQN+NNWC(MM)
         
 10   CONTINUE

      IF(IHUB.EQ.0) GOTO 170
C-----------------------------------------------------------------------
C     Green's formula on the hub
C-----------------------------------------------------------------------
      IEQN=NBW
      DO 180 NN=1,NHBX
         DO 190 MM=1,MHBT

C..........set indices for equation number and influence coefficients...
C..........IEQN is the equation number (row index)......................
            IEQN=IEQN+1

C..........I1 is the index of the control point where GREEN is .........
C..........satisfied....................................................
            I1=INDEXH(NN,MM)

            CALL CAVSETSUB(IEQN,I1)

 190      CONTINUE
 180   CONTINUE

 170   CONTINUE


      IF(IDUCT.EQ.0) GOTO 200
C-----------------------------------------------------------------------
C     Green's formula on the duct
C-----------------------------------------------------------------------

      DO MM=1,MDUCT
         DO NN=1,NDUCT

C..........set indices for equation number and influence coefficients...
C..........IEQN is the equation number (row index)......................
            IEQN=IEQN+1

C..........I1 is the index of the control point where GREEN is .........
C..........satisfied....................................................
            I1=INDEXD(NN,MM)

            CALL CAVSETSUB(IEQN,I1)

         ENDDO
      ENDDo

 200   CONTINUE

      IF(ITUN .EQ. 0) GOTO 270
C-----------------------------------------------------------------------
C     Green's formula on the Tunnel
C-----------------------------------------------------------------------

      DO 280 NN=1,NAXT
         DO 290 MM=1,MTUNEL

            IEQN=IEQN+1

            I1=INDEXTN(NN,MM)

            CALL CAVSETSUB(IEQN,I1)

 290      CONTINUE
 280   CONTINUE

 270   CONTINUE


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
c             CALL CAVSETSUB(IEQN,I1)
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
c             CALL CAVSETSUB(IEQN,I1)
c
c          END DO
c       END DO
c
c 1200  CONTINUE

C-----------------------------------------------------------------------
C     Green's formula in the wake
C-----------------------------------------------------------------------
       IEQN=0
       DO MM=MR,1,-1

C........set indices for equation number and influence coefficients...
C........IEQN is the equation number (row index)......................
          IEQN=IEQN+NC-NSPP(MM,2)-NSPP(MM,1)

          DO NN=1,NNWC(MM)

C...........I1 is the index of the control point where GREEN is .......
C...........satisfied..................................................
             I1=(MR-MM)*NTRA+NN

             IEQN=IEQN+1
             
             CALL CAVSETSUBSC(MM,NN,IEQN,I1)

          END DO
       END DO

C-----------------------------------------------------------------------
C     Put the unsteady wake and the other blade influence in RHS
C-----------------------------------------------------------------------
       CALL CAVRHS

       RETURN
C<<<<<<<<<<<<<<<<<<<<end of Subroutine cavset>>>>>>>>>>>>>>>>>>>>>>>>>>>
       END




