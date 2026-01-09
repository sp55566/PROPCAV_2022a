      SUBROUTINE CAVHT_SP
************************************************************************
*                                                                      *
*  Determine the cavity height such that delta is everywhere positive. *
*                                                                      *
*     Notes:  The different between this subroutine and cavht.f is     *
*             (1) dh/dv is ignored in the cavity height calculation on *
*                 the blade.                                           *
*             (2) dh/dt is also ignored in this subroutine.            *
*             (3) split panels are ignored in this subroutine.         *
*                                                                      *
*  Author: Julie Young  01-16-00                                       *
*                                                                      *
*  Date of last revision                 Revision                      *
*  ---------------------             ----------------                  *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

C-----------------------------------------------------------------------
C     compute the cavity height
C-----------------------------------------------------------------------
      DO 10 M=1,MR

C.......compute cavity height on blade.
         IDR=1
         K=1
         ISF=0

C.......initialize cavity height and cavity T.E. thickness.
         DO 30 N=1,NHP+NTRA
            HT(N,M,IDR)=ZERO
 30      CONTINUE
         DELTA(M,IDR)=ZERO

         JCAV=JCV(M,IDR)
         IF(JCAV.EQ.0) GO TO 25

         M0M=M0(M,IDR)
         M0M1=M0M-1

         CH=HALF*(CHORD(M)+CHORD(M+1))

C.......calculate cavity height on the foil.
         DO 40 N=1,JCAV
            N1=M0M1+ISF+K*N
            J=INDEXB(N1,M)
            
            T1=COSPHI(J)
            T2=-SINPHI(J)
            T12=T1*T1
            T3=DPHIDS(N,M,IDR)-T2*DPHIDV(N,M,IDR)

            DPHIDN=DPDN(J)-DPDNC(J)
            DSI=DS(N1,M)

C..........Use central difference in s-dir
            ANUM=T12*DPHIDN*DSI+HT(N,M,IDR)*T3
            DNOM=T3
            
            HT(N+1,M,IDR)=ANUM/DNOM

 40      CONTINUE

         NC1=JCAV+1

C.......openness at the cavity trailing edge.
         DELTA(M,IDR)=HALF*HT(NC1,M,IDR)/CH 

C.......If there is a supercavity, then set DELTA(M,IDR)=ZERO so
C.......the length on the blade won't change.
         IF(ISC.EQ.1.AND.SOP(M,IDR).EQ.ONE.AND.
     *        HT(NC1,M,IDR).GE.ZERO) DELTA(M,IDR)=ZERO
            
 25      CONTINUE
         
         IF(ISC.EQ.1) THEN
            DO IDR=1,2
               IF(IDR.EQ.1) THEN
                  K=1
                  ISF=0
               ELSE
                  K=-1
                  ISF=1
               END IF
               JCAV=JCV(M,IDR)

               IF(JCAV.EQ.0) THEN
                  NC1=1
               ELSE
                  NC1=JCAV+1
               END IF

C.............computing cavity height on the SR.
               N1=N0(IDR)-1+ISF

               NL1=0
               DO N=1,NSR2
                  IF(IDR.EQ.2) THEN
                     N11=N0(2)-N
                  ELSE
                     N11=N0(1)-1+N
                  END IF
                  IF(ISUBM(INDEXB(N11,M),IDXREV).EQ.1) THEN
                     NL1=NL1+1
                  END IF
               END DO

               DO N=1,NL1
                  NC1=NC1+1
                  N1=N1+K*1
                  J=INDEXB(N1,M)
                  
                  T1=COSPHI(J)
                  T2=-SINPHI(J)
                  T12=T1*T1
                  T3=DPHIDS2(N,M,IDR)-T2*DPHIDV2(N,M,IDR)
                  
                  DPHIDN=DPDN(J)-DPDNC(J)
                  DSI=DS(N1,M)

                  ANUM=T12*DPHIDN*DSI+HT(NC1-1,M,IDR)*T3
                  DNOM=T3
                  
                  HT(NC1,M,IDR)=ANUM/DNOM
               END DO

            END DO

         END IF

 10   CONTINUE

C....compute supercavity height on wake.
      DO 70 M=1,MR
         IF(NNWC(M).EQ.0) GO TO 70

         CH=HALF*(CHORD(M)+CHORD(M+1))
         IDR=1

         IF(ISC.EQ.0) THEN
            CALL CALHTW(IDR,M,HTW1(M))
         ELSE
            CALL CALHTW(1,M,HTWTOP)
            CALL CALHTW(2,M,HTWBOT)
            HTW1(M)=HTWTOP+HTWBOT
         END IF

         JCAV=JCV(M,IDR)
         IF(ISC.EQ.0) THEN
            N00=JCAV+1
         ELSE
            IF(JCAV.EQ.0) THEN
               N00=NSR2P
            ELSE
               N00=NHP-NLEP(M,IDXREV,IDR)
            END IF
         END IF

         NSC=NNWC(M)

         DO 80 N=1,NSC
            N1=N00+N
            J=(MR-M)*NTRA+N
            DSI=DZW(N,M)

            IF(N.EQ.1) THEN
               HN=HTW1(M)
            ELSE
               HN=HT(N1-1,M,IDR)
            END IF

            IF(ISC.EQ.0) THEN
               T1=DSI/ABS(DPHIDS(JCAV+N,M,IDR))
            ELSE
               T1=DSI/ABS(DPHIDS2(NSR2+N,M,IDR))
            END IF

            DPHIDN=-SORW(J)

            ANUM=HN+DPHIDN*T1
            DNOM=ONE

            HT(N1,M,IDR)=ANUM/DNOM

 80      CONTINUE

C.......openness at the cavity trailing edge.
         IF(ISC.EQ.0) THEN
            DELTA(M,IDR)=HALF*HT(N1,M,IDR)/CH
         ELSE
            DELTAT1(M)=HALF*HT(N1,M,IDR)/CH
         END IF

 70   CONTINUE

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine CAVHT_SP>>>>>>>>>>>>>>>>>>>>>>>>
      END



