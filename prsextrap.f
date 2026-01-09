       SUBROUTINE PRSEXTRAP
************************************************************************
*                                                                      *
*      This subroutine re-calculates the pressures for M>MRTIP using   *
*      the extrapolated velocities.                                    *
*                                                                      *
*      Author: Julie Young                                             *
*                                                                      *
*      Date        Revision/comments                                   *
*      --------    ---------------                                     *
*      JY100901    Copied from prsextrap.f.  Modified for ISC=1.       *
*                                                                      *
************************************************************************

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       DIMENSION ICAV(NBZ),NCAV(NBZ),NDIR(NBZ)

       THETAT=-FLOAT(NTSTEP-1)*DELTAT
       ADVCO2=ADVCO*ADVCO

       ISR=1
       IF(IFACE.EQ.2) ISR=2

       IF(ISC.EQ.0) THEN
          NF=1
          NL=NC
       ELSE
          NF=N0(2)
          NL=N0(1)-1
       END IF

       DO M=MRTIP+1,MR

C.......Flag with cavitating panels (including the split panel w/ ICAV..
         DO N=NF,NL
            ICAV(N)=0
         END DO

         IF(IWET.EQ.0) THEN
            DO II=1,ISR
               IF((IFACE.EQ.0).OR.(II.EQ.1.AND.IFACE.EQ.2)) THEN
                  IDR=1
                  K=1
                  ISF=0
               ELSE
                  IDR=2
                  K=-1
                  ISF=1
               END IF
               M0M1=M0(M,IDR)-1
               DO N=1,JCV(M,IDR)+NSPP(M,IDR)
                  N1=M0M1+K*N+ISF
                  ICAV(N1)=1
                  NCAV(N1)=N
                  NDIR(N1)=IDR
               END DO
            END DO
         END IF

         DO N=NF,NL
            L=INDEXB(N,M)

C.........graviational term.............................................
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

C --- Use DPOTDT earlier in UNS case ---- Yiran Su 09/26/2017
C === DpDt should be enabled long enough after the transition to fully coupling
C === The is to ensure it does not create strong fake potential field
            IF ((IUNS.EQ.1).AND.(NTSTEP_KEY.GT.30)) THEN
              DPDT=DPDTPRE(IDXREV,L)
            END IF


C..........calculate the magnitude of the total velocity base on the
C..........extrapolated velocities in global coordinates
            VTOTS(L)=UXTOT(N,M)*UXTOT(N,M)+UYTOT(N,M)*UYTOT(N,M)
     *           +UZTOT(N,M)*UZTOT(N,M)

C..........calculate the total velocities in local panel coordinates
            UXUT=UXTOT(N,M)*UL(L,1)+UYTOT(N,M)*UL(L,2)+
     *           UZTOT(N,M)*UL(L,3)
            UXVT=UXTOT(N,M)*VL(L,1)+UYTOT(N,M)*VL(L,2)+
     *           UZTOT(N,M)*VL(L,3)

C..........calculate the inflow velocities in local panel coordinates
            UXU=VOX1(L)*UL(L,1)+VOY1(L)*UL(L,2)+VOZ1(L)*UL(L,3)
            UXV=VOX1(L)*VL(L,1)+VOY1(L)*VL(L,2)+VOZ1(L)*VL(L,3)

C..........calculate the perturbation velocities in local panel
C..........coordinates
            DPDUB(N,M)=UXUT-UXU
            DPDVB(N,M)=UXVT-UXV

C..........calculate perturbation and total velocity in the w-direction
C..........of the local panel coordinates
            UXETA=(DPDVB(N,M)-DPDUB(N,M)*SINPHI(L)) / COSPHI(L)
            UXETAT=VETAB(N,M)+UXETA

C..........If ITAN=1, the crossflow will be zeroed. (JY053000)
            IF(ISC.EQ.0) THEN
               ITAN=0
            ELSE
               ITAN=1
            END IF

            IF(IWET.EQ.0.AND.ICON.NE.5.AND.ICON.NE.6) THEN
c               IF(ISP.EQ.1.OR.IFACE.EQ.2) THEN
               IF(ISP.EQ.1) THEN
                  ITAN=1
               END IF
               IF(ICON.EQ.7) ITAN=1
            END IF

            IF(ITAN.EQ.1) THEN
               UXETAT=ZERO
               IF(N.LE.NC/2) THEN
                  DPDVB(N,M)=ABS(UXUT)*SINPHI(L)-UXU
               ELSE
                  DPDVB(N,M)=-ABS(UXUT)*SINPHI(L)-UXU
               END IF
               VTOTS(L)=UXUT**2.
            END IF

C..........make sure the cavitating velocities are consistant with DBC
            IF(ICAV(N).EQ.1) THEN
               UTS=DPHIDS(NCAV(N),M,NDIR(N))
               UTV=DPHIDV(NCAV(N),M,NDIR(N))
               IF(NDIR(N).EQ.2) THEN
                  UTW=(UTV-UTS*SINPHI(L))/COSPHI(L)
               ELSE IF(NDIR(N).EQ.1) THEN
                  UTW=(UTV+UTS*SINPHI(L))/COSPHI(L)
               END IF

               VTOTS(L)=UTS**2.+UTW**2.
               UXUT=DPHIDS(NCAV(N),M,NDIR(N))
               IF(NDIR(N).EQ.1) UXUT=-UXUT
               UXETAT=UTW
            END IF

C..........re-calculate the total velocities in global coordinate if
C..........they have been changed
            IF(ITAN.EQ.1.OR.ICAV(N).EQ.1) THEN
               UXTOT(N,M)=UXUT*DIR(L,1,1)+UXETAT*DIR(L,2,1)
               UYTOT(N,M)=UXUT*DIR(L,1,2)+UXETAT*DIR(L,2,2)
               UZTOT(N,M)=UXUT*DIR(L,1,3)+UXETAT*DIR(L,2,3)
            END IF

C..........calculate the new pressure
            CPB(N,M)=VINFSB(N,M)-VTOTS(L)-2.0*DPDT-TWOGY

C..........IF cavitating, the max -Cp should be sigma on the
C..........extrapolated portion.
            IF(IWET.EQ.0) THEN
               IF(-CPB(N,M)*ADVCO2.GT.SIGMA)
     *              CPB(N,M)=-SIGMA/ADVCO2
            END IF

         END DO
      END DO

      RETURN
      END


