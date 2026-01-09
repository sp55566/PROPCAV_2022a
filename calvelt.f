C --------------------------------------------
      subroutine calvelt
C --------------------------------------------
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'

      mcvtm = mcvt - 1 
        DO M = 1 , mcvt
          DO N = 2 , nthx-1

            L0 = INDEXT(N,M)
            L1 = INDEXT(N-1,M)
            L2 = INDEXT(N+1,M)

            X1 = HALF * (DELU(L1) + DELU(L0))
            X2 = HALF * (DELU(L2) + DELU(L0))

            DPDUTH(N,M) = DPDXQ(X1,X2,Pot(l1),Pot(l0),Pot(l2))
        ENDDO

          DPDUTH(nthx,M)=DPDUTH(nthx-1,M)+
     %                  (DPDUTH(nthx-1,M)-DPDUTH(nthx-2,M))*X2/X1
          L0= INDEXT(2,M)
          L1= INDEXT(3,M)
          L2= INDEXT(1,M)
          DPDUTH(1,M)=DPDUTH(2,M)+(DPDUTH(2,M)-DPDUTH(3,M))
     %             *(DELU(L2)+DELU(L0))/(DELU(L1)+DELU(L0))
      ENDDO

        DO N = 1 , nthx
         DO M=2,mcvt - 1

            L0= INDEXT(N,M)
            L1= INDEXT(N,M+1)
            L2= INDEXT(N,M-1)

            X1=0.5*(DELV(L1)+DELV(L0))
            X2=0.5*(DELV(L2)+DELV(L0))
            DPDVTH(N,M)=DPDXQ(X1,X2,Pot(l1),Pot(l0),Pot(l2))
       ENDDO
         DPDVTH(N,mcvt)=DPDVTH(N,mcvt-1)
     %                   +(DPDVTH(N,mcvt-1)-DPDVTH(N,mcvt-2))*X1/X2
         L0= INDEXT(N,2)
         L1= INDEXT(N,1)
         L2= INDEXT(N,3)
         DPDVTH(N,1)=DPDVTH(N,2)+(DPDVTH(N,2)-DPDVTH(N,3))
     *             *(DELV(L1)+DELV(L0))/(DELV(L2)+DELV(L0))
      ENDDO


       THETAT=-FLOAT(NTSTEP-1)*DELTAT
       ADVCO2=ADVCO*ADVCO
       PI2=PI*PI

       DO M = 1 , mcvt
          DO N = 1 , nthx
            L = INDEXT(N,M)
            RCP=SQRT(XCT(L,2)*XCT(L,2)+XCT(L,3)*XCT(L,3))
            THETA=THETAT+ACOS(XCT(L,2)/RCP)
            RCP2=RCP*RCP
            IF((ISTEADY.EQ.0).OR.(ISTEADY.EQ.3))THEN
               TWOGY = 0.
CSH---------------shree----------07/08/2002--------
          IF((ITERGB.EQ.1).OR.(IFLAP.EQ.1))THEN
             ZTEMP = RCP*COS(THETA)
             TWOGY=-ZTEMP*2./(FNRUDDER)
                     ENDIF
CSH---------------shree----------07/08/2002--------
            ELSE
               TWOGY=RCP*COS(THETA)/(FROUDE*ADVCO2)
            ENDIF

C..........time dericavtive of potential................................
            IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.IDXREV.EQ.0) THEN
               DPDT=0.0
            ELSE
               IF(NREV.GT.2) THEN
                  DPDT=DPDTPRE(IDXREV,L)
               ELSE
                  DPDT=0.
               END IF
            END IF

           UXI = DPDUTH(N,M)
           UETA = (DPDVTH(N,M) - DPDUTH(N,M)*SINPHI(L))/COSPHI(L)

         UXIT = VXITH(N,M) + UXI
         UETAT = VETATH(N,M) + UETA
         
         UXTHTOT(N,M)=UXIT*DIR(L,1,1)+UETAT*DIR(L,2,1)
         UYTHTOT(N,M)=UXIT*DIR(L,1,2)+UETAT*DIR(L,2,2)
         UZTHTOT(N,M)=UXIT*DIR(L,1,3)+UETAT*DIR(L,2,3)  
           VTOTS(L) = UXIT**2 + UETAT**2

           CPTHN(N,M) = VINFSTH(N,M) - VTOTS(L) - 2.0*DPDT
           CPTH(N,M) = CPTHN(N,M) - TWOGY

           VTTH(N,M) = SQRT(VTOTS(L))
        ENDDO
      ENDDO

      RETURN
      END


