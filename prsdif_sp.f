      SUBROUTINE PRSDIF_SP
************************************************************************
*     PRSDIF: PReSsure calculations for ISP=1.                         *
************************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

!     PARAMETER(NMM=NBZ+2)
!     DIMENSION DUMX(NMM),DUMY(NMM),DUMC(4*NMM-4)
      DIMENSION DUMX(NBZ+2),DUMY(NBZ+2),DUMC(4*(NBZ+2)-4)
      integer NMM
      NMM = NBZ+2

C....DTN is the nondimensionalized delta(t).............................
      DTN=DELTAT*ADVCO/PI

      DO M=1,MR

         IF(ICB(M,2,IDXREV).EQ.1) THEN         
            IF(ISC.EQ.0) THEN
               NWET=IW(2,M,IDXREV)-IW(1,M,IDXREV)+1
               NFRST=IW(1,M,IDXREV)
            ELSE
               NFRST=MAX(N0(2),IW(1,M,IDXREV))
               NWET=IW(2,M,IDXREV)-NFRST+1
            END IF

            NLAST=IW(2,M,IDXREV)

            DS1=SZ(NFRST,M)
            DS2=SZ(NLAST+1,M)
            
            IF(NWET.GE.5) THEN
               N1=NFRST+2
               N2=NLAST-2
               CALL CALCPB(N1,N2,M)

               DUMX(1)=ZERO
               DUMY(1)=ZERO
               IC1=1
               DO N=N1,N2
                  IC1=IC1+1
                  DUMX(IC1)=SPZ(N,M)-DS1
                  DUMY(IC1)=CPB(N,M)
               END DO
               IC1=IC1+1
               DUMX(IC1)=DS2-DS1
               DUMY(IC1)=ZERO

               IF(IC1.LT.5) THEN
                  CALL QUADFIT(IC1,DUMX,DUMY,A0,A1,A2)
                  DO N=NFRST,NLAST
                     XX1=SPZ(N,M)-DS1
                     CPB(N,M)=A0+A1*XX1+A2*XX1**2.
                  END DO
               ELSE
                  CALL UGLYDK(IC1,1,1,DUMX,DUMY,0.,0.,DUMC)
                  DO N=NFRST,NLAST
                     XX1=SPZ(N,M)-DS1
                     CALL EVALDKs(IC1,1,DUMX,XX1,YY1,DUMC)
                     CPB(N,M)=YY1
                  END DO
               END IF               
            ELSE IF(NFRST.GT.2.AND.NWET.GT.2) THEN
               CALL CALCPB(NFRST,NLAST,M)

               N1=NFRST-2
               N2=NLAST+2

               DS1=SZ(N1,M)
               DS2=SZ(N2+1,M)
               
               DO N=1,3
                  DUMX(N)=SZ(N1+N-1,M)-DS1
                  DUMY(N)=ZERO
               END DO
               IC1=3
               DO N=NFRST,NLAST
                  IC1=IC1+1
                  DUMX(IC1)=SPZ(N,M)-DS1
                  DUMY(IC1)=CPB(N,M)
               END DO
               DO N=1,3
                  IC1=IC1+1
                  DUMX(IC1)=SZ(NLAST+N,M)-DS1
                  DUMY(IC1)=ZERO
               END DO
               CALL QUADFIT(IC1,DUMX,DUMY,A0,A1,A2)
               DO N=NFRST,NLAST
                  XX1=SPZ(N,M)-DS1
                  CPB(N,M)=A0+A1*XX1+A2*XX1**2.
               END DO
            ELSE
               DO N=NFRST,NLAST
                  CPB(N,M)=ZERO
               END DO
            END IF
         
            DO N=NFRST,NLAST
               IF(-CPB(N,M).GT.ZERO) CPB(N,M)=ZERO
            END DO

            DO N=1,NFRST-1
               CPB(N,M)=ZERO
            END DO

            DO N=NLAST+1,NH
               CPB(N,M)=ZERO
            END DO

         ELSE
            DO N=1,NH
               CPB(N,M)=ZERO
            END DO
         END IF

         IF(ICB(M,1,IDXREV).EQ.1) THEN

            NWET=MAX(0,NLEP(M,IDXREV,1)-(IC(1,M,IDXREV)-NHP))
            NFRST=IC(1,M,IDXREV)

            IF(JCV(M,1).EQ.0) THEN
               IF(ISC.EQ.0) THEN
                  NLAST=IC(2,M,IDXREV)
               ELSE
                  NLAST=MIN(IC(2,M,IDXREV),N0(1)-1)
               END IF
            ELSE
               NLAST=M0(M,1)-1
            END IF

            DS1=SZ(NFRST,M)
            DS2=SZ(NLAST+1,M)
            
            IF(NWET.GE.5) THEN
               N1=NFRST+2
               N2=NLAST-2
               CALL CALCPB(N1,N2,M)
               
               DUMX(1)=ZERO
               DUMY(1)=ZERO
               IC1=1
               DO N=N1,N2
                  IC1=IC1+1
                  DUMX(IC1)=SPZ(N,M)-DS1
                  DUMY(IC1)=CPB(N,M)
               END DO
               IC1=IC1+1
               DUMX(IC1)=DS2-DS1
               DUMY(IC1)=ZERO
               
               IF(IC1.LT.5) THEN
                  CALL QUADFIT(IC1,DUMX,DUMY,A0,A1,A2)
                  DO N=NFRST,NLAST
                     XX1=SPZ(N,M)-DS1
                     CPB(N,M)=A0+A1*XX1+A2*XX1**2.
                  END DO
               ELSE
                  CALL UGLYDK(IC1,1,1,DUMX,DUMY,0.,0.,DUMC)
                  DO N=NFRST,NLAST
                     XX1=SPZ(N,M)-DS1
                     CALL EVALDKs(IC1,1,DUMX,XX1,YY1,DUMC)
                     CPB(N,M)=YY1
                  END DO
               END IF      
               
            ELSE IF(NWET.GT.2.AND.NLAST.LT.NC-1) THEN
               CALL CALCPB(NFRST,NLAST,M)
               
               N1=NFRST-2
               N2=NLAST+2
               
               DS1=SZ(N1,M)
               DS2=SZ(N2+1,M)
               
               DO N=1,3
                  DUMX(N)=SZ(N1+N-1,M)-DS1
                  DUMY(N)=ZERO
               END DO
               IC1=3
               DO N=NFRST,NLAST
                  IC1=IC1+1
                  DUMX(IC1)=SPZ(N,M)-DS1
                  DUMY(IC1)=CPB(N,M)
               END DO
               DO N=1,3
                  IC1=IC1+1
                  DUMX(IC1)=SZ(NLAST+N,M)-DS1
                  DUMY(IC1)=ZERO
               END DO
               
               CALL QUADFIT(IC1,DUMX,DUMY,A0,A1,A2)
               DO N=NFRST,NLAST
                  XX1=SPZ(N,M)-DS1
                  CPB(N,M)=A0+A1*XX1+A2*XX1**2.
                  IF(-CPB(N,M).GT.ZERO) CPB(N,M)=ZERO
               END DO
               
            ELSE 

               DO N=NFRST,NLAST
                  CPB(N,M)=ZERO
               END DO

            END IF


            DO N=NFRST,NLAST
               IF(-CPB(N,M).GT.ZERO) CPB(N,M)=ZERO
            END DO

            DO N=NHP,NFRST-1
               CPB(N,M)=ZERO
            END DO

            DO N=NLAST+1,NC
               CPB(N,M)=ZERO
            END DO

         ELSE
            DO N=NHP,NC
               CPB(N,M)=ZERO
            END DO
         END IF
            
         DO N=1,NC
            L=INDEXB(N,M)

            IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1) THEN
               DPDT=0.0
            ELSE
               NMIN=2
               IF(NBLADE.GT.1) NMIN=3
               IF(NREV.GT.NMIN) THEN
                  DPDT=DPDTPRE(IDXREV,L)
               ELSE
                  DPDT=0.
               END IF
            END IF
            
            IF(ISUBM(L,IDXREV).EQ.1) THEN
               QC2=-CPB(N,M)+VINFSB(N,M)-2.0*DPDT
            ELSE
               QC2=ZERO
            END IF

            IF(QC2.LT.ZERO) THEN
               VTOTS(L)=ZERO
            ELSE
               VTOTS(L)=QC2
            END IF

         END DO

         DELCP(M)=CPB(1,M)-CPB(NC,M)
      END DO

C-----------------------------------------------------------------------
C     Compute the potential differences at the trailing edge
C-----------------------------------------------------------------------
      DO M=1,MR
         L1=INDEXB(1,M)
         L2=INDEXB(NC,M)
         IF(ISUBM(L1,IDXREV).EQ.1.AND.
     *        ISUBM(L2,IDXREV).EQ.1) THEN
            DELP(M)=POT(L2)-POT(L1)
         ELSE
            DELP(M)=ZERO
         END IF
      END DO

      RETURN
C<<<<<<<<<<<<<<<<<<<<<<End of subroutine PRSDIF_SP>>>>>>>>>>>>>>>>>>>>>>
      END


      SUBROUTINE CALCPB(N1,N2,M)
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      
      DO N=N1,N2
         L0= INDEXB(N,M)
         L2= INDEXB(N-1,M)
         L1= INDEXB(N+1,M)
         X1=HALF*(DELU(L1)+DELU(L0))
         X2=HALF*(DELU(L2)+DELU(L0))
         DPDUB(N,M)=DPDXQ(X1,X2,POT(L1),POT(L0),POT(L2))

         L=INDEXB(N,M)
         UXIT=VXIB(N,M)+DPDUB(N,M)
         VTOTS(L)=UXIT**2.
         
         IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1) THEN
            DPDT=0.0
         ELSE
            NMIN=2
            IF(NBLADE.GT.1) NMIN=3
            IF(NREV.GT.NMIN) THEN
               DPDT=DPDTPRE(IDXREV,L)
            ELSE
               DPDT=0.
            END IF
         END IF
         
         CPB(N,M)=VINFSB(N,M)-VTOTS(L)-2.0*DPDT
         CPBN(N,M) = CPB(N,M)
      END DO
      
      RETURN
      END


