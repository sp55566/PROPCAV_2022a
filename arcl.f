      SUBROUTINE ARCL
************************************************************************
*     define the linear arclength at the panel boundaries and cp's     *
*     03-22-92 NF                                                      *
*      ----------           ----------                                 *
*     |          |         |          |                                *
*     |  ITER1   |-------->|   ARCL   |                                *
*     |          |         |          |-->                             *
*      -----^----           ----------   |                             *
*           |                            v                             *
*           <-----------------------------                             *
*     Date       Revision                                              *
*     --------   -----------                                           *
*     041597 CM  Added the matrix, ARCLNG.  It is needed for plotting  *
*     out the heights of cavities w.r.t. arclength                     *
*     091098 JY  Corrected the calculation of ARCLNG.  Now the ARLNG   *
*     calculation will also work for face cavitation.                  *
*     020699 JY  modified subroutine to allow cavity to grow on both   *
*     the back and face of the foil.                                   *
*     IDR = 1  (back)                                                  *
*     IDR = 2  (face)                                                  *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      DO 10 M=1,MR
       SZ(1,M)=ZERO
       DO 20 N=1,NC
        DU=HALF*(XB(N+1,M)+XB(N+1,M+1)-XB(N,M)-XB(N,M+1))
        DV=HALF*(YB(N+1,M)+YB(N+1,M+1)-YB(N,M)-YB(N,M+1))
        DW=HALF*(ZB(N+1,M)+ZB(N+1,M+1)-ZB(N,M)-ZB(N,M+1))
        DZ(N,M)=SQRT(DU*DU+DV*DV+DW*DW)
        SZ(N+1,M)=SZ(N,M)+DS(N,M)
        SPZ(N,M)=SZ(N,M)+HALF*DS(N,M)
 20    CONTINUE

       DO 30 N=1,NTRA
        N1=NC+N
        DU=HALF*(XWS(N+1,M)+XWS(N+1,M+1)-XWS(N,M)-XWS(N,M+1))
        DV=HALF*(YWS(N+1,M)+YWS(N+1,M+1)-YWS(N,M)-YWS(N,M+1))
        DW=HALF*(ZWS(N+1,M)+ZWS(N+1,M+1)-ZWS(N,M)-ZWS(N,M+1))
        DZW(N,M)=SQRT(DU*DU+DV*DV+DW*DW)
        SZ(N1+1,M)=SZ(N1,M)+DZW(N,M)
        SPZ(N1,M)=SZ(N1,M)+HALF*DZW(N,M)            
 30    CONTINUE
       DO 40 IDR=1,2
        ARCLNG(1,M,IDR)=ZERO
        IF(IDR.EQ.1) THEN
         K=1
         ISF=0
        ELSE IF(IDR.EQ.2) THEN
         K=-1
         ISF=1
        END IF
        DO 50 N=2, NC/2+1
         INDX=NC/2+K*(N-1)+ISF
         ARCLNG(N,M,IDR)=ARCLNG(N-1,M,IDR)+DS(INDX,M)
 50     CONTINUE
        DO 60 N=2,NTRA+1
         ARCLNG(NC/2+N,M,IDR)=ARCLNG(NC/2+N-1,M,IDR)+DZW(N-1,M)
 60     CONTINUE
 40    CONTINUE

 10   CONTINUE

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine ARCL>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END
