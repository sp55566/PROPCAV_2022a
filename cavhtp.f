       SUBROUTINE CAVHTP
 
       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'

C-----------------------------------------------------------------------
C    Save the cavity height to HTP.  Note HTP is defined from blade LE,
C    NOT cavity LE.                                             
C-----------------------------------------------------------------------
       ISR=1
       IF(IFACE.EQ.2) ISR=2
       
       DO 165 M=1,MR
          
          DO 180 I=1,ISR
             IF((IFACE.EQ.0).OR.(I.EQ.1.AND.IFACE.EQ.2)) THEN
                IDR=1
             ELSE
                IDR=2
             END IF      
             
             IF(JCV(M,IDR).EQ.0) THEN
                IF(ISC.EQ.0) THEN
                   DO N=1,NHP
                      HTP(N,M,IDR)=ZERO
                   END DO
                ELSE
                   N1=NHOLD+1
                   DO N=1,N1
                      HTP(N,M,IDR)=ZERO
                   END DO
                   DO N=1,NSR2
                      N1=N1+1
                      HTP(N1,M,IDR)=HT(N+1,M,IDR)
                   END DO
                END IF
             ELSE
                N1=1
                DO N=1,NHP
                   IF(N.GT.NLEP(M,IDXREV,IDR)+1) THEN
                      N1=N1+1
                      IF(NSPP(M,IDR).EQ.1.AND.
     *                     N1.EQ.JCV(M,IDR)+2) THEN
                         HTP(N,M,IDR)=ZERO
                      ELSE
                         HTP(N,M,IDR)=HT(N1,M,IDR)
                      END IF
                   ELSE
                      HTP(N,M,IDR)=ZERO
                   END IF
                END DO
             END IF
 180      CONTINUE 
          
          IF(NNWC(M).EQ.0) THEN
             DO N=1,NTRA+1
                HTWP(N,M)=ZERO
             END DO
          ELSE
             IDR=NWDIR(M)
             HTWP(1,M)=HTW1(M)
             DO N=1,NNWC(M)
                N2=JCV(M,IDR)+N+1
                
                IF(ISC.EQ.1) THEN
                   N2=NSR2+N+1
                   IF(JCV(M,IDR).GT.0) N2=N2+
     *                  NHOLD-NLEP(M,IDXREV,IDR)
                END IF
                HTWP(N+1,M)=HT(N2,M,IDR)
             END DO
             DO N=NNWC(M)+1,NTRA
                HTWP(N+1,M)=ZERO
             END DO
          END IF
          
 165   CONTINUE
       
       RETURN
       END
