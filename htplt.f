      SUBROUTINE HTPLT
************************************************************************
*                                                                      *
*  Subroutine HeighT PLoT plots out the cavity so that it may be       *
*  viewed with Tecplot                                                 *
*                           ----------        ----------               *
*                          |          |      |          |              *
*                          |    HT    |> --->|  CAVPLT  |              *
*                          |          |      |          |-->           *
*                           -----^----        ----------   |           *
*                                |                         |           *
*                                 --<---------<------------v           *
*  Date      Comment or Revision                                       *
*  --------  -------------------                                       *
*  CM050598  Added value of Delta as last height.                      *
*  JY091198  Subroutine modified by J.Young.                           *
*  JY030999  Modified subroutine to allow cavity to grow on both the   *
*            back and face of the foil.                                *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
C      SAVE

      CHARACTER*29 FNHT

      IF(IDXREV.EQ.1)THEN
         CALL CHRLEN(FN,LENCH)
         FNHT=FN(1:LENCH)//'.ht'
         OPEN(611, FILE=FNHT, STATUS='UNKNOWN')
         WRITE(611,5000)
 5000    FORMAT(1x,'TITLE="Heights"')
         WRITE(611,5500)
 5500    FORMAT(1x,'VARIABLES="X/C", "Height"')
      ENDIF

      ISR=1            
      IF(IFACE.EQ.2) ISR=2
 6000 FORMAT(1X,'ZONE T="BACK: strip',I2,' at time',I3,'"')
 6500 FORMAT(1X,'ZONE T="FACE: strip',I2,' at time',I3,'"')

      DO 10 M=1, MR           
         CH=HALF*(CHORD(M)+CHORD(M+1))
         DO 20 II=1,ISR
            IF((IFACE.EQ.0).OR.(II.EQ.1.AND.IFACE.EQ.2)) THEN        
               IDR=1
               K=1
            ELSE
               IDR=2
               K=-1
            END IF

            N1=JCV(M,IDR)+NSPP(M,IDR)+1
            IF(NWC(M,IDR).GT.0.AND.IDR.EQ.NWDIR(M)) 
     *           N1=N1+NWC(M,IDR)+NSPS(M,IDR)
            
C..........Add IF statement so that the cavity heights are only plotted
C..........if cavity do exist in that strip. (JY012300)
            IF(N1.GT.1) THEN

            IF(IDR.EQ.1) THEN
               WRITE(611,6000) M,NTSTEP
            ELSE IF(IDR.EQ.2) THEN
               WRITE(611,6500) M,NTSTEP
            END IF

            ISPLIT=0
            DNOM=ARCLNG(NHP,M,IDR)
            DO 30 N=1,N1
               IF(NSPP(M,IDR).EQ.1.AND.N.EQ.JCV(M,IDR)+NSPP(M,IDR)+1) 
     *              ISPLIT=1
               IF(NSPS(M,IDR).EQ.1.AND.N.EQ.N1) ISPLIT=2
               INDX=NLEP(M,IDXREV,IDR)+N
               IF(ISPLIT.EQ.0) THEN
                  WRITE(611,*) ARCLNG(INDX,M,IDR)/DNOM,HT(N,M,IDR)
               ELSE IF(ISPLIT.EQ.1) THEN
                  WRITE(611,*) (ARCLNG(INDX-1,M,IDR)+
     *                 (ARCLNG(INDX,M,IDR)-ARCLNG(INDX-1,M,IDR))*
     *                 FLP(M,IDR))/DNOM,DELTA(M,IDR)*CH/HALF
               ELSE IF(ISPLIT.EQ.2) THEN
                  WRITE(611,*) (ARCLNG(INDX-1,M,IDR)+
     *                 (ARCLNG(INDX,M,IDR)-ARCLNG(INDX-1,M,IDR))*
     *                 FLS(M,IDR))/DNOM,DELTA(M,IDR)*CH/HALF
               END IF
 30         CONTINUE

            END IF
C..........End of changes (JY012300)....................................

 20      CONTINUE
 10   CONTINUE

C....If IFILE=1, close file 611 at the end of every revolution.
      IF(IFILE.EQ.1.AND.IDXREV.EQ.NTPREV) CLOSE(611)
      IF(IFILE.EQ.1.AND.ISTEADY.EQ.0) CLOSE(611)

      RETURN

C<<<<<<<<<<<<<<<<<<<<<end of subroutine CAVPLT>>>>>>>>>>>>>>>>>>>>>>>>>>
      END





