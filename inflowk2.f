      SUBROUTINE INFLOWK2
************************************************************************
*     INFLOW2: INFLOW velocities  at the wake control points           *
*                                                                      *
*     Date     Comment or Revision                                     *
*     -------- -------------------                                     *
*     10/13/99   HSLEE                                           *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      DIMENSION U(3),NH1(3)

      IF(NTSTEP.EQ.0) THEN
        DO 10 IH=1,3
          NH1(IH)=1
   10   CONTINUE
      ELSE
        DO 15 IH=1,3
          NH1(IH)=NHARM(IH)
   15   CONTINUE
      END IF

      DO 40 J=1,NPWAKE
        RCP=SQRT(XCTW(J,2)**2+XCTW(J,3)**2)
        THP=ATAN2(XCTW(J,3),XCTW(J,2))
        DO 30 IH=1,3
          U(IH)=0.0
          DO 20 JH=1,NH1(IH)
            CALL EVALDKs(NWKCOE,1,XRW,RCP,COEC,XWCUB(1,JH,1,IH))
            CALL EVALDKs(NWKCOE,1,XRW,RCP,COES,XWCUB(1,JH,2,IH))
            T=FLOAT(JH-1)*(THP-TSTEP)
            ST=SIN(T)
            CT1=COS(T)
            U(IH)=U(IH)+COEC*CT1+COES*ST
   20     CONTINUE
   30   CONTINUE
        VOXW(J)=U(1)
        VORW=U(2)
        VOTW=U(3)
        COCP=XCTW(J,2)/RCP
        SICP=XCTW(J,3)/RCP
        WROVS=PI*RCP/ADVCO
        VOYW(J)=VORW*COCP-(WROVS+VOTW)*SICP
        VOZW(J)=VORW*SICP+(WROVS+VOTW)*COCP
   40 CONTINUE


      RETURN
      END
