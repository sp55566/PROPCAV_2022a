       SUBROUTINE CAVSC1(IEQN,I1)
************************************************************************
*                                                                      *
*  SET up the linear system of equations for solving the unsteady      *
*  cavitating propeller problem for ISC=1. Originally was setup for the*
*  wing problem. CAVSET computes the entire LHS (which changes         *
*  with every timestep because cavity size changes) and the part of    *
*  the RHS which is associated with the cavity on the key blade. The   *
*  remainder of the RHS, namely the influence of the other blades and  *
*  the wakes of all the blades, is computed in CAVRHS.                 *
*                                                                      *
*  This subroutine is called inside CAVSETSC.  This sets up the eqns   *
*  for equation IEQN for the given column I1.  This part is the        *
*  for Green's formula on the blade, the hub, the bulb, and the        *
*  tip vortex.                                                         *
*                                                                      *
*     ALHS(I,J) = array which contains the left hand side of the matrix*
*                equation                                              *
*     RHS(I) = array which contains the right hand side of the matrix  *
*                equation                                              *
*                                                                      *
************************************************************************
      USE MEMSOL
      USE CVRHS 
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'

C      COMMON/MEMSOL/ALHS(NTZ,NTZ)
C      COMMON/CVRHS/RHS(NTZ)

C....IFRST1 marks the beginning of the unknowns on the foil strip.
C....within a column of the matrix................................
      IFRST1=0

      RHS(IEQN)=ZERO
            
      DO 40 M=MR,1,-1

C.......set indices for inner summations..........................
         M0MB=M0(M,1)
         JCAVB=JCV(M,1)
         LCAVB=M0MB+JCAVB+NSPP(M,1)
         ISPB=LCAVB-1
               
         M0MF=M0(M,2)
         JCAVF=JCV(M,2)
         LCAVF=M0MF-JCAVF-NSPP(M,2)-1
         ISPF=LCAVF+1

         IFRST=INDEXB(0,M)

C.......compute the left-hand-side due to the blade...............
         SUMF=0.
         SUMB=0.
         FCTF=0.
         FCTB=0.

         SUMF2=0.
         SUMB2=0.
         FCTF2=0.
         FCTB2=0.
               
C.......cavitating panels on face SR region
         DO N=1,N0(2)-1
            ALHS(IEQN,IFRST1+N)=-BB(I1,IFRST+N)
            SUMF2=SUMF2+AA(I1,IFRST+N)
         END DO

C.......terms associated with PHI0 extrapolation of the
C.......face SR region.
         IF(ISTEADY.EQ.0) THEN
            FCTF2=SUMF2-W(I1,M)
         ELSE
            FCTF2=SUMF2-(HALF*W(I1,M)+WSUBIF(I1,M))
         END IF
               
         IF(JCAVF.GE.1) THEN

C..........wetted panels on face near T.E.........................
            N123=LCAVF
            DO 50 N=N0(2),N123
               ALHS(IEQN,IFRST1+N)=AA(I1,IFRST+N)
 50         CONTINUE

C..........cavitating panels on face..............................
            DO 70 N=1,JCAVF
               N123 = M0MF-N-NSPP(M,2)
               N456 = M0MF-N
               ALHS(IEQN,IFRST1+N123)=-BB(I1,IFRST+N456)
               SUMF=SUMF+AA(I1,IFRST+N456)
 70         CONTINUE

C..........fully wetted panels on face near L.E. and fully wetted.
C..........panels on back before cavity L.E.......................
            IF(JCAVB.EQ.0) THEN
               NLAST=N0(1)-1
            ELSE
               NLAST=M0MB-1
            END IF

            DO 80 N=M0MF,NLAST
               N1=N-NSPP(M,2)
               ALHS(IEQN,IFRST1+N1)=AA(I1,IFRST+N)
 80         CONTINUE

C..........terms associated with phi0 extrapolation on face.......
            N1=IFRST1+M0MF-NSPP(M,2)
            N2=N1+1
            N3=N1+2

            FCTF=SUMF
            IF(NSPP(M,2).NE.0) FCTF=FCTF+
     *           AA(I1,IFRST+ISPF)*FLP(M,2)
            IF(SOP(M,2).EQ.ONE) FCTF=FCTF+FCTF2
                  
            ALHS(IEQN,N1)=ALHS(IEQN,N1)+FCTF*CT(M,1,2)
            ALHS(IEQN,N2)=ALHS(IEQN,N2)+FCTF*CT(M,2,2)
            ALHS(IEQN,N3)=ALHS(IEQN,N3)+FCTF*CT(M,3,2)
            
            IF(NSPP(M,2).EQ.0) GO TO 93

C..........terms associated with dipole extrapolation right of....
C..........split panel on face....................................
            N1=IFRST1+ISPF-1
            N2=N1-1
            N3=N2-1
            N4=N3-1
            FCT=AA(I1,IFRST+ISPF)*FRP(M,2)
            IF((ISPF-1).GE.N0(2)) ALHS(IEQN,N1)=ALHS(IEQN,N1)+
     *           DT(M,1,2)*FCT
            IF((ISPF-2).GE.N0(2)) ALHS(IEQN,N2)=ALHS(IEQN,N2)+
     *           DT(M,2,2)*FCT
            IF((ISPF-3).GE.N0(2)) ALHS(IEQN,N3)=ALHS(IEQN,N3)+
     *           DT(M,3,2)*FCT
            IF((ISPF-4).GE.N0(2)) ALHS(IEQN,N4)=ALHS(IEQN,N4)+
     *           DT(M,4,2)*FCT

C..........terms associated with source extrapolatin left of .....
C..........split panel on face....................................
            N1=IFRST1+ISPF
            N2=N1+1
            N3=N2+1
            N4=N3+1
            FCT=BB(I1,IFRST+ISPF)*FLP(M,2)
            ALHS(IEQN,N1)=ALHS(IEQN,N1)-QT(M,1,2)*FCT
            ALHS(IEQN,N2)=ALHS(IEQN,N2)-QT(M,2,2)*FCT
            ALHS(IEQN,N3)=ALHS(IEQN,N3)-QT(M,3,2)*FCT
            ALHS(IEQN,N4)=ALHS(IEQN,N4)-QT(M,4,2)*FCT

 93         CONTINUE

C.......define LHS when face is fully wetted......................
         ELSE

C..........fully wetted panels on face and fully wetted panels....
C..........before L.E. of cavity on back..........................
            IF(JCAVB.EQ.0) THEN
               NLAST=N0(1)-1
            ELSE
               NLAST=M0MB-1
            END IF
                 
            DO 90 N=N0(2),NLAST
               ALHS(IEQN,IFRST1+N)=AA(I1,IFRST+N)
 90         CONTINUE
         END IF

C.......only if wetted or partially cavitating
         IF(SOP(M,2).EQ.ZERO) THEN
            N1=IFRST1+N0(2)
            N2=N1+1
            N3=N1+2
            ALHS(IEQN,N1)=ALHS(IEQN,N1)+FCTF2*CT2(M,1,2)
            ALHS(IEQN,N2)=ALHS(IEQN,N2)+FCTF2*CT2(M,2,2)
            ALHS(IEQN,N3)=ALHS(IEQN,N3)+FCTF2*CT2(M,3,2)
         END IF

C.......cavitating panels on back SR region
         DO N=N0(1),NC
            N123=N-NSPP(M,2)-NSPP(M,1)
            ALHS(IEQN,IFRST1+N123)=-BB(I1,IFRST+N)
            SUMB2=SUMB2+AA(I1,IFRST+N)
         END DO

C.......terms associated with PHI0 extrapolation of the
C.......back SR region.
         IF(ISTEADY.EQ.0) THEN
            FCTB2=SUMB2+W(I1,M)
         ELSE
            FCTB2=SUMB2+(HALF*W(I1,M)+WSUBIF(I1,M))
         END IF
               
C.......define LHS when back is cavitating........................
         IF(JCAVB.GE.1) THEN

C..........cavitating panels on back..............................
            DO 110 N=1,JCAVB
               N123 = M0MB-1+N-NSPP(M,2)
               N456 = M0MB-1+N 
               ALHS(IEQN,IFRST1+N123)=-BB(I1,IFRST+N456)
               SUMB=SUMB+AA(I1,IFRST+N456)
 110        CONTINUE

C..........fully wetted panels on back near T.E...................
            DO 120 N=LCAVB,N0(1)-1
               N1=N-NSPP(M,2)-NSPP(M,1)
               ALHS(IEQN,IFRST1+N1)=AA(I1,IFRST+N)
 120        CONTINUE

C..........terms associated with phi0 extrapolation on back.......
            N1=IFRST1+M0MB-1-NSPP(M,2)
            N2=N1-1
            N3=N1-2                 
            FCTB=SUMB
            IF(NSPP(M,1).NE.0) FCTB=FCTB+
     *           AA(I1,IFRST+ISPB)*FLP(M,1)
            IF(SOP(M,1).EQ.ONE) FCTB=FCTB+FCTB2

            ALHS(IEQN,N1)=ALHS(IEQN,N1)+FCTB*CT(M,1,1)
            ALHS(IEQN,N2)=ALHS(IEQN,N2)+FCTB*CT(M,2,1)
            ALHS(IEQN,N3)=ALHS(IEQN,N3)+FCTB*CT(M,3,1)

            IF(NSPP(M,1).EQ.0) GO TO 125

C..........terms associated with dipole extrapolation right of....
C..........split panel on back....................................
            N1=IFRST1+ISPB+1-NSPP(M,1)-NSPP(M,2)
            N2=N1+1
            N3=N2+1
            N4=N3+1
            FCT=AA(I1,IFRST+ISPB)*FRP(M,1)
            IF((ISPB+1).LT.N0(1)) ALHS(IEQN,N1)=ALHS(IEQN,N1)+
     *           DT(M,1,1)*FCT
            IF((ISPB+2).LT.N0(1)) ALHS(IEQN,N2)=ALHS(IEQN,N2)+
     *           DT(M,2,1)*FCT
            IF((ISPB+3).LT.N0(1)) ALHS(IEQN,N3)=ALHS(IEQN,N3)+
     *           DT(M,3,1)*FCT
            IF((ISPB+4).LT.N0(1)) ALHS(IEQN,N4)=ALHS(IEQN,N4)+
     *           DT(M,4,1)*FCT

C..........terms associated with source extrapolatin left of .....
C..........split panel on back....................................
            N1=IFRST1+ISPB-1-NSPP(M,2)
            N2=N1-1
            N3=N2-1
            N4=N3-1
            FCT=BB(I1,IFRST+ISPB)*FLP(M,1)

            ALHS(IEQN,N1)=ALHS(IEQN,N1)-QT(M,1,1)*FCT
            ALHS(IEQN,N2)=ALHS(IEQN,N2)-QT(M,2,1)*FCT
            ALHS(IEQN,N3)=ALHS(IEQN,N3)-QT(M,3,1)*FCT
            ALHS(IEQN,N4)=ALHS(IEQN,N4)-QT(M,4,1)*FCT
            
 125        CONTINUE

         END IF

C.......only if wetted or partially cavitating
         IF(SOP(M,1).EQ.ZERO) THEN
            N1=IFRST1+N0(1)-1-NSPP(M,2)-NSPP(M,1)
            N2=N1-1
            N3=N1-2

            ALHS(IEQN,N1)=ALHS(IEQN,N1)+FCTB2*CT2(M,1,1)
            ALHS(IEQN,N2)=ALHS(IEQN,N2)+FCTB2*CT2(M,2,1)
            ALHS(IEQN,N3)=ALHS(IEQN,N3)+FCTB2*CT2(M,3,1)

         END IF

C.......terms associated with supercavitating wake panels.........
         IF(NNWC(M).EQ.0) GO TO 150

         DO 140 N=1,NNWC(M)
            J=(MR-M)*NTRA+N
            J1=IFRST1+NC-NSPP(M,2)-NSPP(M,1)+N
            ALHS(IEQN,J1)=-C(I1,J)
 140     CONTINUE

         IF(NSPS(M,NWDIR(M)).EQ.0) GOTO 150

C.......terms associated with source extrapolation of left split..
C.......panel on supercavitating wake.............................
         J2=J1-1
         J3=J1-2
         J4=J1-3
         J=(MR-M)*NTRA+NNWC(M)+1
         FCT=C(I1,J)*FLS(M,NWDIR(M))
         ALHS(IEQN,J1)=ALHS(IEQN,J1)-FCT*QW(M,1,NWDIR(M))
         ALHS(IEQN,J2)=ALHS(IEQN,J2)-FCT*QW(M,2,NWDIR(M))
         ALHS(IEQN,J3)=ALHS(IEQN,J3)-FCT*QW(M,3,NWDIR(M))
         ALHS(IEQN,J4)=ALHS(IEQN,J4)-FCT*QW(M,4,NWDIR(M))
 150     CONTINUE

C.......compute the right-hand-side due to the blade..............
         SUM2=ZERO
         SUM3=ZERO
         TERM1=ZERO
         TERM2=ZERO
         TERM3=ZERO
         TERM4=ZERO
         TERM5=ZERO

C.......cavitating panels on face SR region
         DO N=1,N0(2)-1
            N123=N0(2)-N
            SUM2=SUM2+AA(I1,IFRST+N)*PHI2(N123,M,2)
         END DO

         IF(SOP(M,2).EQ.ZERO) THEN
C.......only if wetted or partially cavitating
            TERM3=TERM3-CT2(M,4,2)*QCSR(1,M,2)*FCTF2
         ELSE
C.......only if supercavitating
            TERM3=TERM3-PSI0T(M,2)*FCTF2
         END IF

C.......define RHS when face is cavitating........................
         IF(JCAVF.GE.1) THEN

C..........wetted panels on face near T.E.........................
            N123=LCAVF
            DO 160 N=N0(2),N123
               N1=IFRST+N
               SUM3=SUM3+BB(I1,N1)*DPDNC(N1)
 160        CONTINUE                  

C..........cavitating panels on face..............................
            DO 170 N=1,JCAVF
               SUM2=SUM2+AA(I1,IFRST+M0MF-N)*PHI1(N,M,2)          
 170        CONTINUE

C..........fully wetted panels on the face near L.E. and fully....
C..........wetted panels on back before cavity L.E................
            IF(JCAVB.EQ.0) THEN
               NLAST=N0(1)-1
            ELSE
               NLAST=M0MB-1
            END IF

            DO 180 N=M0MF,NLAST
               N1=IFRST+N
               SUM3=SUM3+BB(I1,N1)*DPDNC(N1)
 180        CONTINUE

C..........define terms associated with phi0 extrapolation........
            TERM3=TERM3-CT(M,4,2)*QC(1,M,2)*FCTF

C..........define terms associated with dipole extrapolation......
C..........left of split panel and source extrapolation right of..
C..........split panel on face....................................
            IF(NSPP(M,2).NE.0) THEN
               TERM4=TERM4-AA(I1,IFRST+ISPF)*DLISP(M,2)*FLP(M,2)+
     *              BB(I1,IFRST+ISPF)*QSPR(M,2)*FRP(M,2)
            END IF

C.......define RHS when face is fully wetted......................
         ELSE
                  
C..........fully wetted panels on the face and fully wetted ......
C..........panels on back before cavity LE........................
            IF(JCAVB.EQ.0) THEN
               NLAST=N0(1)-1
            ELSE
               NLAST=M0MB-1
            END IF

            DO 190 N=N0(2),NLAST
               N1=IFRST+N
               SUM3=SUM3+BB(I1,N1)*DPDNC(N1)
 190        CONTINUE
         END IF

C.......cavitating panels on back SR region
         DO N=N0(1),NC
            N123=N-N0(1)+1
            SUM2=SUM2+AA(I1,IFRST+N)*PHI2(N123,M,1)
         END DO

         IF(SOP(M,1).EQ.ZERO) THEN
C..........only if wetted or partially cavitating
            TERM3=TERM3-CT2(M,4,1)*QCSR(1,M,1)*FCTB2
         ELSE
C..........only if supercavitating
            TERM3=TERM3-PSI0T(M,1)*FCTB2
         END IF

C.......define RHS when back is cavitating........................
         IF(JCAVB.GE.1) THEN

C..........cavitating panels on back..............................
            DO 200 N=1,JCAVB
               SUM2=SUM2+AA(I1,IFRST+M0MB-1+N)*PHI1(N,M,1)  
 200        CONTINUE

C..........fully wetted panels on back near T.E...................
            DO 210 N=LCAVB,N0(1)-1
               N1=IFRST+N
               SUM3=SUM3+BB(I1,N1)*DPDNC(N1)
 210        CONTINUE

C..........define terms associated with phi0 extrapolation........
            TERM3=TERM3-CT(M,4,1)*QC(1,M,1)*FCTB

C..........define terms associated with dipole extrapolation......
C..........left of split panel and source extrapolation right of..
C..........split panel on back....................................
            IF(NSPP(M,1).NE.0) THEN
               TERM4=TERM4-AA(I1,IFRST+ISPB)*DLISP(M,1)*FLP(M,1)+
     *              BB(I1,IFRST+ISPB)*QSPR(M,1)*FRP(M,1)
            END IF

         END IF

C.......sum the sources and dipoles on the non-split panels.......
         TERM1=-SUM2+SUM3

C.......define term associateds with kutta condition on foil......
         IF(ISTEADY.EQ.0) THEN
            DUM1=W(I1,M)
         ELSE
            DUM1=HALF*W(I1,M)+WSUBIF(I1,M)
         END IF

         TERM2=TERM2+DUM1*(PHI2(NSR2,M,2)-PHI2(NSR2,M,1))

C.......define term associated with trailing sinks due to.........
C.......collapsed supercavities (JY060800)
         IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NREV.LE.2) THEN
            TERM5=ZERO
         ELSE
C.HS071007..IF(NNWC(M).GT.0.AND.NSPS(M,NWDIR(M)).EQ.1) THEN 
            IF(NNWC(M).GT.0) THEN
              IF(NSPS(M,NWDIR(M)).EQ.1) THEN
                J=(MR-M)*NTRA+NNWC(M)+1
                TERM5=TERM5+C(I1,J)*QSSR(M)*FRS(M,NWDIR(M))
              ENDIF
            END IF
C........Corrected by Hong (HS071007)

            IF(NNWC(M).GT.0) THEN
               N11=NNWC(M)+NSPS(M,NWDIR(M))+1
            ELSE
               N11=1
            END IF
            DO 230 N=N11,NTRA
               J=(MR-M)*NTRA+N
               TERM5=TERM5+C(I1,J)*SORW(J)
 230        CONTINUE
         END IF

C.......sum all terms to form RHS.................................
         RHS(IEQN)=RHS(IEQN)+TERM1+TERM2+TERM3+TERM4+TERM5

         IFRST1=IFRST1+NC-NSPP(M,2)-NSPP(M,1)+NNWC(M)

 40   CONTINUE

      IF(IHUB.NE.0)THEN
C.......cal. terms to LHS that are associated with the hub...........
         DO 240 N=1,NHBX
            DO 250 M=1,MHBT
               J=INDEXH(N,M)
               J1=J-NPANB
               ALHS(IEQN,NBW+J1)=AA(I1,J)
 250        CONTINUE
                  
C..........cal. terms to RHS that are associated with the hub.....
            SUM1=ZERO
            DO 260 M=1,MHBT
               J=INDEXH(N,M)
               SUM1=SUM1+BB(I1,J)*DPDNC(J)
 260        CONTINUE
            RHS(IEQN)=RHS(IEQN)+SUM1
 240     CONTINUE
      ENDIF

c      IF(IAN.EQ.2) THEN
c         ICCOUNT=NBW+NPANH
c         
cC.......cal. terms to LHS that are associated with the 
cC..........bulb on the blade tip
c         DO N=1,NTHX
c            DO M=1,MCVT
c               J=INDEXT(N,M)
c               J1=J - NPANB - NPANH
c               ALHS(IEQN,ICCOUNT+J1)=AA(I1,J)
c            ENDDO
c            
cC..........cal. terms to RHS that are associated with the 
cC..........bulb on the blade tip
c            SUM1=ZERO
c            DO  M=1,MCVT
c               J=INDEXT(N,M)
c               SUM1=SUM1+BB(I1,J)*DPDNC(J)
c            ENDDO
c            RHS(IEQN)=RHS(IEQN)+SUM1
c         ENDDO
c
cC.......cal. terms to LHS that are associated with the 
cC.......tip vortex
c         DO N=1,NCVX
c            DO M=1,MCVT
c               J=INDEXC(N,M)
c               J1=J - NPANB - NPANH
c               ALHS(IEQN,ICCOUNT+J1)=AA(I1,J)
c            ENDDO
c            
cC..........cal. terms to RHS that are associated with the 
cC..........tip vortex
c            SUM1=ZERO
c            DO  M=1,MCVT
c               J=INDEXC(N,M)
c               SUM1=SUM1+BB(I1,J)*DPDNC(J)
c            ENDDO
c            RHS(IEQN)=RHS(IEQN)+SUM1
c         ENDDO
c         
c      END IF

      RETURN
      END
