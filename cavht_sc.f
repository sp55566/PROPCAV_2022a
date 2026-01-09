      SUBROUTINE CAVHT_SC
************************************************************************
*                                                                      *
*  Determine the cavity height such that delta is everywhere positive  *
*  2-8-91  Neal Fine                                                   *
*                                                                      *
*  Date of last revision                 Revision                      *
*  ---------------------             ----------------                  *
*  JY091901               Copied from cavht.f.  Modified for ISC=1.    *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

C-----------------------------------------------------------------------
C     compute the cavity height
C-----------------------------------------------------------------------

C....DTN is the time step...............................................
      DTN=DELTAT*ADVCO/PI

C....T5 controls the on and off of the time term..................
      T5=ZERO

      DO 10 M=1,MR

C.......compute cavity height on blade..................................
         DO 20 IDR=1,2
            IF(IDR.EQ.1) THEN
               K=1
               ISF=0
            ELSE
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

            IF(JCAV.EQ.0) GOTO 25

C..........calculate cavity height on the foil..........................
            DO 40 N=1,JCAV+NSPP(M,IDR)
               N1=M0M1+ISF+K*N
               J=INDEXB(N1,M)
               IF(M.NE.1) J1=INDEXB(N1,M-1)
               T1=COSPHI(J)
               T2=-SINPHI(J)
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
                  IF(T5.NE.ZERO) THEN
                     ANUM=T12*DPHIDN*DSI+HT(N,M,IDR)*(T3-
     *                    T5*HALF*T12*DSI/DTN)
                     DNOM=T3+T5*HALF*T12*DSI/DTN
                  ELSE
                     ANUM=T12*DPHIDN*DSI+HT(N,M,IDR)*(T3)
                     DNOM=T3
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

            NC1=JCAV+NSPP(M,IDR)+1

C..........openness at the cavity trailing edge.........................
            DELTA(M,IDR)=HALF*HT(NC1,M,IDR)/CH 

C..........If there is a supercavity, then set DELTA(M,IDR)=ZERO so
C..........the length on the blade won't change.
            IF(SOP(M,IDR).EQ.ONE) THEN
               IF(HT(NC1,M,IDR).LT.ZERO) THEN
                  CAVLTE=ARCLNG(NHP-NSR2,M,IDR)/ARCLNG(NHP,M,IDR)
                  CAVL(M,IDR)=CAVLTE
               ELSE
                  DELTA(M,IDR)=ZERO
               END IF
            ELSE

C.............setting the cavity height equal to zero at blade panels
C.............between the cavity and the SR.
               N1=JCAV+NSPP(M,IDR)+1
               DO N=LCAV,N0(IDR)-1+ISF,K
                  N1=N1+1
                  HT(N1,M,IDR)=ZERO
               END DO
               NC1=N1
            END IF

 25         IF(JCAV.EQ.0) NC1=1

C..........computing cavity height on the SR.
            N1=N0(IDR)-1+ISF
            DO N=1,NSR2
               NC1=NC1+1
               N1=N1+K*1
               J=INDEXB(N1,M)
               IF(M.NE.1) J1=INDEXB(N1,M-1)

               T1=COSPHI(J)
               T2=-SINPHI(J)
               T12=T1*T1
               T3=DPHIDS2(N,M,IDR)-T2*DPHIDV2(N,M,IDR)
               T4=DPHIDV2(N,M,IDR)-T2*DPHIDS2(N,M,IDR)
            
               NP1=NHOLD+N

               IF(T5.NE.ZERO) THEN
                  HOLD=T5*HALF*(HTP(NP1+1,M,IDR)+HTP(NP1,M,IDR))/DTN
               ELSE
                  HOLD=ZERO
               END IF

               DPHIDN=DPDN(J)-DPDNC(J)+HOLD
               DSI=DS(N1,M)

               IF(M.EQ.1)THEN

C................If M=1, assume the crossflow term is zero..............
C................use central difference in s-dir, and backward..........
C................difference in time.....................................
                  IF(T5.NE.ZERO) THEN
                     ANUM=T12*DPHIDN*DSI+HT(NC1-1,M,IDR)*(T3-
     *                    T5*HALF*T12*DSI/DTN)
                     DNOM=T3+T5*HALF*T12*DSI/DTN
                  ELSE
                     ANUM=T12*DPHIDN*DSI+HT(NC1-1,M,IDR)*(T3)
                     DNOM=T3
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
                  HBARM1=HALF*(HT(NC1,M-1,IDR)+HT(NC1-1,M-1,IDR))

                  ANUM=T12*DPHIDN*DSI+HT(NC1-1,M,IDR)*(T3-HALF*DSI*DUM1)
     *                 +T4*DSI/DV1*HBARM1
                  DNOM=T3+HALF*DSI*DUM1
               ENDIF

               HT(NC1,M,IDR)=ANUM/DNOM

            END DO

 20      CONTINUE
 10   CONTINUE

C....compute supercavity height on wake.................................
      DO 70 M=1,MR

         CH=HALF*(CHORD(M)+CHORD(M+1))
         IF(NWDIR(M).EQ.1) THEN
            IDR=1
         ELSE IF(NWDIR(M).EQ.2) THEN
            IDR=2
         END IF

         CALL CALHTW(1,M,HTWTOP)
         CALL CALHTW(2,M,HTWBOT)
         HTW1(M)=HTWTOP+HTWBOT
         
         IF(NNWC(M).EQ.0) THEN
            DELTAT1(M)=HALF*HTW1(M)/CH
            GO TO 70
         END IF

         IF(JCV(M,IDR).EQ.0) THEN
            N00=NSR2P
         ELSE
            N00=NHP-NLEP(M,IDXREV,IDR)
         END IF

         NSC=NNWC(M)
         NSP=NSPS(M,IDR)

         DO 80 N=1,NSC+NSP
            N1=N00+N
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
               DSI=FLS(M,IDR)*DZW(NSC+1,M)
               HTENM1=HTWP(N,M)+(HTWP(N+1,M)-
     *              HTWP(N,M))*FLS(M,IDR)
               IF(T5.NE.ZERO) THEN
                  HOLD=T5*HALF*(HTENM1+HTWP(N,M))/DTN
               ELSE
                  HOLD=ZERO
               END IF
            ENDIF

            T1=DSI/ABS(DPHIDS2(NSR2+N,M,IDR))
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
         DELTAT1(M)=HALF*HT(N1,M,IDR)/CH

 70   CONTINUE

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine CAVHT>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END



