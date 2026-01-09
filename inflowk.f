      SUBROUTINE INFLOWK
************************************************************************
*     INFLOW: INFLOW velocities  at the wake control points            *
*     08-04-92 NF                                                      *
*                                                                      *
*     Date     Comment or Revision                                     *
*     -------- -------------------                                     *
*     CM010598 Changes made by H.S. Lee noted below for ICON = 5       *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

      DIMENSION U(3),NH1(3)

      IF(NTSTEP.EQ.0.OR.ISTEADY.EQ.0) THEN
        DO 10 IH=1,3
          NH1(IH)=1               
   10   CONTINUE
      ELSE
         DO 15 IH=1,3
            NH1(IH)=NHARM(IH)               
 15      CONTINUE
      END IF

C-----------------------------------------------------------------------
C        Next 7 lines added by HL for hydrofoil case           CM010598
C-----------------------------------------------------------------------
C....Same for ICON=8 as well (JY110100)
      if((icon.eq.5).or.(icon.eq.6).or.(icon.eq.8)) then
       do 100 j = 1 , npwaks
            voxw(j) = cos(alpha*rad)
            voyw(j) = 0.0
            vozw(j) = sin(alpha*rad)
 100     continue

CSHREE------------------------------------/04/24/2002-------

      else

      DO 40 J=1,NPWAKS

        RCP=SQRT(XCPW(J,2,1)**2+XCPW(J,3,1)**2)
        THP=ATAN2(XCPW(J,3,1),XCPW(J,2,1))
C.......Compute inflow velocities based on the wake harmonics...........
        DO 30 IH=1,3
          U(IH)=0.0
          DO 20 JH=1,NH1(IH)
C...........Spline wake harmonics to control radii (assume same.........
C...........number of panels for each blade)............................
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
        COCP=XCPW(J,2,1)/RCP
        SICP=XCPW(J,3,1)/RCP
        WROVS=PI*RCP/ADVCO
        VOYW(J)=VORW*COCP-(WROVS+VOTW)*SICP
        VOZW(J)=VORW*SICP+(WROVS+VOTW)*COCP
   40 CONTINUE

      endif

      RETURN
C>>>>>>>>>>>>>>>>>>>>>>End of subroutine INFLOWK>>>>>>>>>>>>>>>>>>>>>>>>
      END
