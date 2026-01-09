      SUBROUTINE WAKEPLT2
************************************************************************
*                                                                      *
*  Date     Revision or Comment                                        *
*  -------- -------------------                                        *
*  JY081499 This subroutine was copied from wakeplt.f so we can see    *
*           the wake geometry where the supercavity panels are placed. *
*  JY072800 This subroutine was modified to calculate WAKMRX & WAKVEC. *
*                                                                      *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      DIMENSION WAKMRX(3),DUM1(3),DUM2(3),DUM3(3)

C-----------------------------------------------------------------------
C     Plot wake subpanels where the supercavity panels are placed.
C-----------------------------------------------------------------------
      IF(NTSTEP.EQ.1) OPEN(800,FILE="wake2.plt",STATUS='UNKNOWN')

      IF(NTSTEP.EQ.1) WRITE(800,*) 'VARIABLES="x","y","z"'

      IF(NTRA.GT.0) THEN
         NNODE=(MR+1)*(NTRA+1)
         NELE=MR*NTRA

         WRITE(800,1000) NTSTEP,NNODE,NELE
 1000    FORMAT(1x,'ZONE T="TIMESTEP=',I3,'" N=',I5,' E=',I5,
     *        ' F=FEPOINT ET=QUADRILATERAL')

         DO M=1,MR+1
            DO N123=1,NTRA+1
               WRITE(800,*) XWS(N123,M),YWS(N123,M),ZWS(N123,M)
            END DO
         END DO

         DO M=1,MR
            DO N123=1,NTRA
               WRITE(800,*) (M-1)*(NTRA+1)+N123,(M-1)*(NTRA+1)+N123+1,
     *              M*(NTRA+1)+N123+1,M*(NTRA+1)+N123
            END DO
         END DO
      END IF

      IF(NTSTEP.EQ.NCTIME) CLOSE(800)

C-----------------------------------------------------------------------
C     Calculate and plot WAKMRX & WAKVEC.  
C-----------------------------------------------------------------------

      Do M=1,MR
         DO N=1,NTRA+1
            WAKMRX(1)=(XWS(N,M)+XWS(N,M+1))/2
            WAKMRX(2)=(YWS(N,M)+YWS(N,M+1))/2
            WAKMRX(3)=(ZWS(N,M)+ZWS(N,M+1))/2

C....Same for ICON=8 as well (JY110100)

            IF((ICON.EQ.5).OR.(ICON.EQ.6).OR.(ICON.EQ.8)) THEN
               WAKVEC(1,N,M)=0.
               WAKVEC(2,N,M)=0.
               WAKVEC(3,N,M)=1.
            ELSE
               IF(N.EQ.1) THEN
                  DUM1(1)=HALF*(XWS(2,M)+XWS(2,M+1)-
     *                 XWS(1,M)-XWS(1,M+1))
                  DUM1(2)=HALF*(YWS(2,M)+YWS(2,M+1)-
     *                 YWS(1,M)-YWS(1,M+1))
                  DUM1(3)=HALF*(ZWS(2,M)+ZWS(2,M+1)-
     *                 ZWS(1,M)-ZWS(1,M+1))
               ELSE IF(N.EQ.NTRA+1) THEN
                  DUM1(1)=HALF*(XWS(N,M)+XWS(N,M+1)-
     *                 XWS(N-1,M)-XWS(N-1,M+1))
                  DUM1(2)=HALF*(YWS(N,M)+YWS(N,M+1)-
     *                 YWS(N-1,M)-YWS(N-1,M+1))
                  DUM1(3)=HALF*(ZWS(N,M)+ZWS(N,M+1)-
     *                 ZWS(N-1,M)-ZWS(N-1,M+1))
               ELSE
                  NN1=N+1
                  NN2=N-1
                  DUM1(1)=HALF*(XWS(NN1,M)+XWS(NN1,M+1)-
     *                 XWS(NN2,M)-XWS(NN2,M+1))
                  DUM1(2)=HALF*(YWS(NN1,M)+YWS(NN1,M+1)-
     *                 YWS(NN2,M)-YWS(NN2,M+1))
                  DUM1(3)=HALF*(ZWS(NN1,M)+ZWS(NN1,M+1)-
     *                 ZWS(NN2,M)-ZWS(NN2,M+1))
               END IF

               DUM2(1)=XWS(N,M+1)-XWS(N,M)
               DUM2(2)=YWS(N,M+1)-YWS(N,M)
               DUM2(3)=ZWS(N,M+1)-ZWS(N,M)

               CALL EXPRO(DUM1,DUM2,DUM3)
            
               VECLEN=ZERO
               DO KK=1,3
                  VECLEN=VECLEN+DUM3(KK)**2.
               END DO
               VECLEN=SQRT(VECLEN)

               DO KK=1,3
                  WAKVEC(KK,N,M)=DUM3(KK)/VECLEN
               END DO
            END IF
         END DO
      END DO

      RETURN
      END

      SUBROUTINE WAKEPLT3(KK)
*************************************************************************
*                                                                       *
*  Date        Revision or Comment
*  *
*  ----------- -------------------
*  *
*  SK.Dec.2018 This subroutine was copied from wakeplt.f so we can see
*  *
*                                                                       *
*************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      DIMENSION WAKMRX(3),DUM1(3),DUM2(3),DUM3(3)

C-----------------------------------------------------------------------
C     Plot wake subpanels where the supercavity panels are placed.
C-----------------------------------------------------------------------
      IF(NTSTEP.EQ.1.AND.KK.EQ.NBLADE) THEN
        OPEN(369,FILE="wake_uns.plt",STATUS='UNKNOWN')
        WRITE(369,*) 'VARIABLES="x","y","z"'
      ENDIF

      NNODE=(MR+1)*(NTRA+1)
      NELE=MR*NTRA

      WRITE(369,1000) NTSTEP,NNODE,NELE
 1000    FORMAT(1x,'ZONE T="TIMESTEP=',I3,'" N=',I5,' E=',I5,
     *        ' F=FEPOINT ET=QUADRILATERAL')

      DO M=1,MR+1
        DO N123=1,NTRA+1
           WRITE(369,*) XWS(N123,M),YWS(N123,M),ZWS(N123,M)
        END DO
      END DO

      DO M=1,MR
        DO N123=1,NTRA
          WRITE(369,*) (M-1)*(NTRA+1)+N123,(M-1)*(NTRA+1)+N123+1,
     *              M*(NTRA+1)+N123+1,M*(NTRA+1)+N123
        END DO
      END DO

      IF(NTSTEP.EQ.NCTIME.AND.KK.EQ.1) CLOSE(369)

      RETURN
      END SUBROUTINE
