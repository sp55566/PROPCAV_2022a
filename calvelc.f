C --------------------------------------------
      subroutine calvelc
C --------------------------------------------
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'

      dth = twopi / mcvt

      mcvtm = mcvt - 1 
        DO M = 1 , mcvt
          DO N = 2 , ncvx-1

            L0 = INDEXC(N,M)
            L1 = INDEXC(N-1,M)
            L2 = INDEXC(N+1,M)

            X1 = HALF * (DELU(L1) + DELU(L0))
            X2 = HALF * (DELU(L2) + DELU(L0))

            DPDUC(N,M) = DPDXQ(X1,X2,Pot(l1),Pot(l0),Pot(l2))
        ENDDO

          DPDUC(ncvx,M)=DPDUC(ncvx-1,M)+
     %                  (DPDUC(ncvx-1,M)-DPDUC(ncvx-2,M))*X2/X1
          L0= INDEXC(2,M)
          L1= INDEXC(3,M)
          L2= INDEXC(1,M)
          DPDUC(1,M)=DPDUC(2,M)+(DPDUC(2,M)-DPDUC(3,M))
     %             *(DELU(L2)+DELU(L0))/(DELU(L1)+DELU(L0))
      ENDDO


        DO N = 1 , ncvx
         DO M=2,mcvt - 1

            L0= INDEXC(N,M)
            L1= INDEXC(N,M+1)
            L2= INDEXC(N,M-1)

            X1=0.5*(DELV(L1)+DELV(L0))
            X2=0.5*(DELV(L2)+DELV(L0))
          x1 = rc1 * dth
          x2 = rc1 * dth
            DPDVC(N,M)=DPDXQ(X1,X2,Pot(l1),Pot(l0),Pot(l2))
       ENDDO
         DPDVC(N,mcvt)=DPDVC(N,mcvt-1)
     %                   +(DPDVC(N,mcvt-1)-DPDVC(N,mcvt-2))*X1/X2
         L0= INDEXC(N,2)
         L1= INDEXC(N,1)
         L2= INDEXC(N,3)
       x1 = rc1 * dth
       x2 = rc1 * dth
         DPDVC(N,1)=DPDVC(N,2)+(DPDVC(N,2)-DPDVC(N,3))*x1/x2
c         DPDVC(N,1)=DPDVC(N,2)+(DPDVC(N,2)-DPDVC(N,3))
c     *             *(DELV(L1)+DELV(L0))/(DELV(L2)+DELV(L0))
      ENDDO

       THETAT=-FLOAT(NTSTEP-1)*DELTAT
       ADVCO2=ADVCO*ADVCO
       PI2=PI*PI

       DO M = 1 , mcvt
          DO N = 1 , ncvx
          L = INDEXC(N,M)
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
           UXI = DPDUC(N,M)
           UETA = (DPDVC(N,M) - DPDUC(N,M)*SINPHI(L))/COSPHI(L)

         UXIT = VXIC(N,M) + UXI
         UETAT = VETAC(N,M) + UETA
         UXCHTOT(N,M)=UXIT*DIR(L,1,1)+UETAT*DIR(L,2,1)
         UYCHTOT(N,M)=UXIT*DIR(L,1,2)+UETAT*DIR(L,2,2)
         UZCHTOT(N,M)=UXIT*DIR(L,1,3)+UETAT*DIR(L,2,3)

           VTOTS(L) = UXIT**2 + UETAT**2

           CPCN(N,M) = VINFSC(N,M) - VTOTS(L) - 2.0*DPDT
           CPC(N,M) = CPCN(N,M) - TWOGY

           VTC(N,M) = SQRT(VTOTS(L))
        ENDDO
      ENDDO

      RETURN
      END


