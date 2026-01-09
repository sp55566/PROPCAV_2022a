      SUBROUTINE CAVHT
************************************************************************
*                                                                      *
*  Determine the cavity height such that delta is everywhere positive  *
*  2-8-91  Neal Fine                                                   *
*                                                                      *
*  Date of last revision                 Revision                      *
*  ---------------------             ----------------                  *
*  03-19-91 NF                   -search for correct l(r) by finding   *
*                                 minimum delta                        *
*  05-14-91 NF                   -include dhdv term                    *
*  07-16-91 NF                   -do not search for correct l(r)       *
*  07-24-91 NF                   -split-panel model: last source is    *
*                                 the split-panel source               *
*  10-25-91 NF                   -added cavity height in the wake for  *
*                                 supercavitation.                     *
*  05-11-92 NF      -removed supercavitation and revised for unsteady  *
*                    propeller solution                                *
*  07-28-92 NF      -added supercavitation                             *
*  CM011998      Working on a symmetric height planform                *
*  CM030398      Added some face cavitation changes.                   *
*  JY082898      Removed fort file 54 and replace it with variable HTP *
*  JY090998      Corrected calculation of dh/dt term.                  *
*  JY021499      Modified subroutine to allow cavity to grow on both   *
*                the back and face of the foil.                        *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
*  JY060100      Modified calculation of supercavity height at blade   *
*                trailing edge.                                        *
*  JY060800      Change HTP in wake to HTWP.  This is necessary because*
*                the definition of normal is different on the blade and*
*                on the wake.                                          *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

C-----------------------------------------------------------------------
C     compute the cavity height
C-----------------------------------------------------------------------

C....DTN is the time step...............................................
      DTN=DELTAT*ADVCO/PI

      ISR=1
      IF(IFACE.EQ.2) ISR=2

C....T5 controls the on and off of the time term..................
      T5=ZERO
      IF((NREV.GT.2).AND.(ISTEADY.NE.0).AND.
     *     (ISTEADY.NE.1)) T5=ONE

      IF(IAN .EQ. 2 .AND. ISCAV .EQ. 2 
     *     .AND. NREV .EQ. NNREV) T5 = ZERO

      DO 10 M=1,MR

C.......compute cavity height on blade..................................
         DO 20 II=1,ISR
            IF((IFACE.EQ.0).OR.(II.EQ.1.AND.IFACE.EQ.2)) THEN
               IDR=1
               K=1
               ISF=0
            ELSE
               IDR=2
               K=-1
               ISF=1
            END IF

            SGN=FLOAT(K)
            M0M=M0(M,IDR)
            M0M1=M0M-1
            JCAV=JCV(M,IDR)
            LCAV=M0M+K*(JCAV+NSPP(M,IDR))-ISF
            ISP1=LCAV-K*1

            CH=HALF*(CHORD(M)+CHORD(M+1))

C..........initialize cavity height and cavity T.E. thickness...........
            DO 30 N=1,NHP+NTRA
               HT(N,M,IDR)=ZERO
 30         CONTINUE
            DELTA(M,IDR)=ZERO

            IF(JCAV.EQ.0) GOTO 20

C..........calculate cavity height on the foil..........................
            DO 40 N=1,JCAV+NSPP(M,IDR)
               N1=M0M1+ISF+K*N
               J=INDEXB(N1,M)
               IF(M.NE.1) J1=INDEXB(N1,M-1)
               T1=COSPHI(J)
c               IF(IDR.EQ.1) THEN
                 T2=-SINPHI(J)
c               ELSE
c                 T2=SINPHI(J)
c               ENDIF
               T12=T1*T1
               T3=DPHIDS(N,M,IDR)-T2*DPHIDV(N,M,IDR)
               T4=DPHIDV(N,M,IDR)-T2*DPHIDS(N,M,IDR)

               NP1=NLEP(M,IDXREV,IDR)+N
               IF(T5.NE.ZERO) THEN
                  HOLD=T5*HALF*(HTP(NP1+1,M,IDR)+HTP(NP1,M,IDR))/DTN
               ELSE
                  HOLD=ZERO
               END IF

               DPHIDN=DPDN(J)-DPDNC(J)+HOLD
               DSI=DS(N1,M)

C.............if N is a plit panel, then do the following...............
               IF(N.EQ.(JCAV+1))THEN
                  DSI=DZL(M,IDR)
                  T3=DPHIDS(N,M,IDR)-T2*DPDVSP(M,IDR)
                  T4=DPDVSP(M,IDR)-T2*DPHIDS(N,M,IDR)
                  
                  HTENM1=HTP(NP1,M,IDR)+(HTP(NP1+1,M,IDR)-
     *                 HTP(NP1,M,IDR))*DZL(M,IDR)/DZ(ISP1,M)
                  IF(T5.NE.ZERO) THEN
                     HOLD=T5*HALF*(HTENM1+HTP(NP1,M,IDR))/DTN
                  ELSE
                     HOLD=ZERO
                  END IF
                  DPHIDN=QSPR(M,IDR)-QSPL(M,IDR)+HOLD
               ENDIF

               IF(M.EQ.1)THEN

C................If M=1, assume the crossflow term is zero..............
C................use central difference in s-dir, and backward..........
C................difference in time.....................................
                  DV1=HALF*DELV(J)
                  IF(T5.NE.ZERO) THEN
                    ANUM=T12*DPHIDN*DSI+HT(N,M,IDR)*(T3-
     *                   T5*HALF*T12*DSI/DTN-HALF*T4*DSI/DV1)
                    DNOM=T3+T5*HALF*T12*DSI/DTN
     *                  +HALF*DSI*T4/DV1
                  ELSE
                    ANUM=T12*DPHIDN*DSI+HT(N,M,IDR)*(T3-HALF*T4*DSI/DV1)
                    DNOM=T3+HALF*DSI*T4/DV1
                  END IF

               ELSE

C................Use central difference in s-dir, backward difference...
C................in v-dir, and backward difference in time..............
                  DV1=HALF*(DELV(J)+DELV(J1))
                  IF(T5.NE.ZERO) THEN
                     DUM1=T4/DV1+T5*T12/DTN
                  ELSE
                     DUM1=T4/DV1
                  END IF
                  HBARM1=HALF*(HT(N+1,M-1,IDR)+HT(N,M-1,IDR))

                  ANUM=T12*DPHIDN*DSI+HT(N,M,IDR)*(T3-HALF*DSI*DUM1)
     *                 +T4*DSI/DV1*HBARM1
                  DNOM=T3+HALF*DSI*DUM1
               ENDIF

               HT(N+1,M,IDR)=ANUM/DNOM

 40         CONTINUE

C..........openness at the cavity trailing edge.........................
            NTOT=JCAV+NSPP(M,IDR)
            DELTA(M,IDR)=HALF*HT(NTOT+1,M,IDR)/CH 

 20      CONTINUE
 10   CONTINUE

C....compute supercavity height on wake.................................
      DO 70 M=1,MR

         IF(NNWC(M).EQ.0) GO TO 70
          
         CH=HALF*(CHORD(M)+CHORD(M+1))
         IF(NWDIR(M).EQ.1) THEN
            IDR=1
         ELSE IF(NWDIR(M).EQ.2) THEN
            IDR=2
         END IF

         JCAV=JCV(M,IDR)
         JW=(MR-M)*NTRA+1


         IF(NWC(M,1).GT.0.AND.NWC(M,2).GT.0) THEN
            CALL CALHTW(1,M,HTWTOP)
            CALL CALHTW(2,M,HTWBOT)
            HTW1(M)=HTWTOP+HTWBOT
         ELSE

            CALL CALHTW(IDR,M,HTW1(M))

         END IF

         N0J=JCAV+1
         NSC=NNWC(M)
         NSP=NSPS(M,IDR)

         DO 80 N=1,NSC+NSP
            N1=N0J+N
            J=(MR-M)*NTRA+N
            DSI=DZW(N,M)

            IF(T5.NE.ZERO) THEN
               HOLD=T5*HALF*(HTWP(N,M)+HTWP(N+1,M))/DTN
            ELSE
               HOLD=ZERO
            END IF

            IF(N.EQ.1) THEN
               HN=HTW1(M)
            ELSE
               HN=HT(N1-1,M,IDR)
            END IF
            IF(N.EQ.NSC+1)THEN
               DSI=DZL(M,IDR)               
               HTENM1=HTWP(N,M)+(HTWP(N+1,M)-
     *              HTWP(N,M))*DZL(M,IDR)/DZW(NSC+1,M)
               IF(T5.NE.ZERO) THEN
                  HOLD=T5*HALF*(HTENM1+HTWP(N,M))/DTN
               ELSE
                  HOLD=ZERO
               END IF
            ENDIF

            T1=DSI/ABS(DPHIDS(JCAV+N,M,IDR))
            IF(T5.NE.ZERO) THEN
               T2=T5*HALF/DTN*T1
            ELSE
               T2=ZERO
            END IF
            DPHIDN=-SORW(J)+HOLD 
            ANUM=HN*(ONE-T2)+DPHIDN*T1
            DNOM=ONE+T2
            HT(N1,M,IDR)=ANUM/DNOM

 80      CONTINUE

C.......openness at the cavity trailing edge............................
         NTOT=JCAV+NSC+NSP
         DELTA(M,IDR)=HALF*HT(NTOT+1,M,IDR)/CH

 70   CONTINUE

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine CAVHT>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END



