       SUBROUTINE CAVSETSUB(IEQN,I1)
************************************************************************
*                                                                      *
*  SET up the linear system of equations for solving the unsteady      *
*  cavitating propeller problem. Originally was setup for the wing     *
*  problem. CAVSET computes the entire left-hand-side (which changes   *
*  with every timestep because cavity size changes) and the part of    *
*  the RHS which is associated with the cavity on the key blade. The   *
*  remainder of the RHS, namely the influence of the other blades and  *
*  the wakes of all the blades, is computed in CAVRHS.                 *
*                                                                      *
*  This subroutine is called inside CAVSET.  This sets up the equations*
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
            
      DO M=MR,1,-1

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

         IF(IFACE.EQ.1.OR.IFACE.EQ.2) THEN

C..........define LHS when face is cavitating........................
            IF(JCAVF.GE.1) THEN

C.............wetted panels on face near T.E.........................
               N123=LCAVF
               DO 310 N=1,N123
                  ALHS(IEQN,IFRST1+N)=AA(I1,IFRST+N)
 310           CONTINUE

C.............cavitating panels on face..............................
               DO 330 N=1,JCAVF
                  N123 = M0MF-N-NSPP(M,2)
                  N456 = M0MF-N
                  ALHS(IEQN,IFRST1+N123)=-BB(I1,IFRST+N456)
                  SUMF=SUMF+AA(I1,IFRST+N456)
 330           CONTINUE

C.............fully wetted panels on face near L.E. and fully wetted.
C.............panels on back before cavity L.E.......................
               IF(IFACE.EQ.1) THEN
                  NLAST=NC
               ELSE IF(IFACE.EQ.2) THEN
                  IF(JCAVB.EQ.0) THEN
                     NLAST=NC
                  ELSE
                     NLAST=M0MB-1
                  END IF
               END IF
               
               DO 340 N=M0MF,NLAST
                  N1=N-NSPP(M,2)
                  ALHS(IEQN,IFRST1+N1)=AA(I1,IFRST+N)
 340           CONTINUE

C.............terms associated with phi0 extrapolation on face.......
               N1=IFRST1+M0MF-NSPP(M,2)
               N2=N1+1
               N3=N1+2

               IF(ISTEADY.EQ.0) THEN
                  FCTF=SUMF-W(I1,M)*SOP(M,2)
               ELSE
                  FCTF=SUMF-(HALF*W(I1,M)+WSUBIF(I1,M))*SOP(M,2)
               END IF
               
               IF(NSPP(M,2).NE.0) THEN
                  FCTF=FCTF+AA(I1,IFRST+ISPF)*FLP(M,2)
               END IF
               
               ALHS(IEQN,N1)=ALHS(IEQN,N1)+FCTF*CT(M,1,2)
               ALHS(IEQN,N2)=ALHS(IEQN,N2)+FCTF*CT(M,2,2)
               ALHS(IEQN,N3)=ALHS(IEQN,N3)+FCTF*CT(M,3,2)

               IF(NSPP(M,2).EQ.0) GO TO 345

C.............terms associated with dipole extrapolation right of....
C.............split panel on face....................................
               N1=IFRST1+ISPF-1
               N2=N1-1
               N3=N2-1
               N4=N3-1
               FCT=AA(I1,IFRST+ISPF)*FRP(M,2)
               IF((ISPF-1).GT.0) ALHS(IEQN,N1)=ALHS(IEQN,N1)+
     *              DT(M,1,2)*FCT
               IF((ISPF-2).GT.0) ALHS(IEQN,N2)=ALHS(IEQN,N2)+
     *              DT(M,2,2)*FCT
               IF((ISPF-3).GT.0) ALHS(IEQN,N3)=ALHS(IEQN,N3)+
     *              DT(M,3,2)*FCT
               IF((ISPF-4).GT.0) ALHS(IEQN,N4)=ALHS(IEQN,N4)+
     *              DT(M,4,2)*FCT

C.............terms associated with source extrapolatin left of .....
C.............split panel on face....................................
               N1=IFRST1+ISPF
               N2=N1+1
               N3=N2+1
               N4=N3+1
               FCT=BB(I1,IFRST+ISPF)*FLP(M,2)
               ALHS(IEQN,N1)=ALHS(IEQN,N1)-QT(M,1,2)*FCT
               ALHS(IEQN,N2)=ALHS(IEQN,N2)-QT(M,2,2)*FCT
               ALHS(IEQN,N3)=ALHS(IEQN,N3)-QT(M,3,2)*FCT
               ALHS(IEQN,N4)=ALHS(IEQN,N4)-QT(M,4,2)*FCT

 345           CONTINUE

C..........define LHS when face is fully wetted......................
            ELSE

C.............fully wetted panels on face and fully wetted panels....
C.............before L.E. of cavity on back..........................
               IF(IFACE.EQ.1) THEN
                  NLAST=NC
               ELSE IF(IFACE.EQ.2) THEN
                  IF(JCAVB.EQ.0) THEN
                     NLAST=NC
                  ELSE
                     NLAST=M0MB-1
                  END IF
               END IF
               
               DO 350 N=1,NLAST
                  ALHS(IEQN,IFRST1+N)=AA(I1,IFRST+N)
 350           CONTINUE
            END IF
            
         END IF
       
         IF(IFACE.EQ.0.OR.IFACE.EQ.2) THEN
            
C..........define LHS when back is cavitating........................
            IF(JCAVB.GE.1) THEN

C.............wetted panels on face..................................
               IF(IFACE.EQ.0) THEN
                  DO 355 N=1,M0MB-1
                     ALHS(IEQN,IFRST1+N)=AA(I1,IFRST+N)
 355              CONTINUE
               END IF

C.............cavitating panels on back..............................
               DO 370 N=1,JCAVB
                  N123 = M0MB-1+N-NSPP(M,2)
                  N456 = M0MB-1+N 
                  ALHS(IEQN,IFRST1+N123)=-BB(I1,IFRST+N456)
                  SUMB=SUMB+AA(I1,IFRST+N456)
 370           CONTINUE

C.............fully wetted panels on back near T.E...................
               DO 380 N=LCAVB,NC
                  N1=N-NSPP(M,2)-NSPP(M,1)
                  ALHS(IEQN,IFRST1+N1)=AA(I1,IFRST+N)
 380           CONTINUE

C.............terms associated with phi0 extrapolation on back.......
               N1=IFRST1+M0MB-1-NSPP(M,2)
               N2=N1-1
               N3=N1-2

               IF(ISTEADY.EQ.0) THEN
                  FCTB=SUMB+W(I1,M)*SOP(M,1)
               ELSE
                  FCTB=SUMB+(HALF*W(I1,M)+WSUBIF(I1,M))*SOP(M,1)
               END IF
               
               IF(NSPP(M,1).NE.0) THEN
                  FCTB=FCTB+AA(I1,IFRST+ISPB)*FLP(M,1)
               END IF
               
               ALHS(IEQN,N1)=ALHS(IEQN,N1)+FCTB*CT(M,1,1)
               ALHS(IEQN,N2)=ALHS(IEQN,N2)+FCTB*CT(M,2,1)
               ALHS(IEQN,N3)=ALHS(IEQN,N3)+FCTB*CT(M,3,1)
               
               IF(NSPP(M,1).EQ.0) GO TO 385

C.............terms associated with dipole extrapolation right of....
C.............split panel on back....................................
               N1=IFRST1+ISPB+1-NSPP(M,1)-NSPP(M,2)
               N2=N1+1
               N3=N2+1
               N4=N3+1
               FCT=AA(I1,IFRST+ISPB)*FRP(M,1)
               IF((ISPB+1).LE.NC) ALHS(IEQN,N1)=ALHS(IEQN,N1)+
     *              DT(M,1,1)*FCT
               IF((ISPB+2).LE.NC) ALHS(IEQN,N2)=ALHS(IEQN,N2)+
     *              DT(M,2,1)*FCT
               IF((ISPB+3).LE.NC) ALHS(IEQN,N3)=ALHS(IEQN,N3)+
     *              DT(M,3,1)*FCT
               IF((ISPB+4).LE.NC) ALHS(IEQN,N4)=ALHS(IEQN,N4)+
     *              DT(M,4,1)*FCT

C.............terms associated with source extrapolatin left of .....
C.............split panel on back....................................
               N1=IFRST1+ISPB-1-NSPP(M,2)
               N2=N1-1
               N3=N2-1
               N4=N3-1
               FCT=BB(I1,IFRST+ISPB)*FLP(M,1)
               ALHS(IEQN,N1)=ALHS(IEQN,N1)-QT(M,1,1)*FCT
               ALHS(IEQN,N2)=ALHS(IEQN,N2)-QT(M,2,1)*FCT
               ALHS(IEQN,N3)=ALHS(IEQN,N3)-QT(M,3,1)*FCT
               ALHS(IEQN,N4)=ALHS(IEQN,N4)-QT(M,4,1)*FCT
               
 385           CONTINUE
               
C..........define LHS when back is fully wetted......................
            ELSE
               
C.............fully wetted panels on back............................
               IF(IFACE.EQ.0) THEN
                  DO 390 N=1,NC
                     N1=N-NSPP(M,2)
                     ALHS(IEQN,IFRST1+N1)=AA(I1,IFRST+N)
 390              CONTINUE
               END IF
            END IF
            
         END IF
         
C.......terms associated with kutta conditon on foil..............
         N1=IFRST1+1
         N2=IFRST1+NC-NSPP(M,2)-NSPP(M,1)
         
         IF(ISTEADY.EQ.0) THEN
            ALHS(IEQN,N1)=ALHS(IEQN,N1)-(ONE-SOP(M,2))*W(I1,M)
            ALHS(IEQN,N2)=ALHS(IEQN,N2)+(ONE-SOP(M,1))*W(I1,M)
         ELSE
            ALHS(IEQN,N1)=ALHS(IEQN,N1)-(ONE-SOP(M,2))*
     *           (HALF*W(I1,M)+WSUBIF(I1,M))
            ALHS(IEQN,N2)=ALHS(IEQN,N2)+(ONE-SOP(M,1))*
     *           (HALF*W(I1,M)+WSUBIF(I1,M))
         END IF
         
C.......terms associated with supercavitating wake panels.........
         IF(NNWC(M).EQ.0) GO TO 400
         
         DO 410 N=1,NNWC(M)
            J=(MR-M)*NTRA+N
            J1=IFRST1+NC-NSPP(M,2)-NSPP(M,1)+N
            ALHS(IEQN,J1)=-C(I1,J)
 410     CONTINUE

         IF(NSPS(M,NWDIR(M)).EQ.0) GOTO 400
         
C.......terms associated with source extrapolation of left split..
C.......panel on supercavitating wake.............................
         J2=J1-1
         J3=J1-2
         J4=J1-3
         J=(MR-M)*NTRA+NNWC(M)+1
         ALHS(IEQN,J1)=ALHS(IEQN,J1)-C(I1,J)*QW(M,1,NWDIR(M))*
     *        FLS(M,NWDIR(M))
         ALHS(IEQN,J2)=ALHS(IEQN,J2)-C(I1,J)*QW(M,2,NWDIR(M))*
     *        FLS(M,NWDIR(M))
         ALHS(IEQN,J3)=ALHS(IEQN,J3)-C(I1,J)*QW(M,3,NWDIR(M))*
     *        FLS(M,NWDIR(M))
         ALHS(IEQN,J4)=ALHS(IEQN,J4)-C(I1,J)*QW(M,4,NWDIR(M))*
     *        FLS(M,NWDIR(M))
 400     CONTINUE

C.......compute the right-hand-side due to the blade..............
         SUM2=ZERO
         SUM3=ZERO
         TERM1=ZERO
         TERM2=ZERO
         TERM3=ZERO
         TERM4=ZERO
         TERM5=ZERO
         
         IF(IFACE.EQ.1.OR.IFACE.EQ.2) THEN
            
C..........define RHS when face is cavitating........................
            IF(JCAVF.GE.1) THEN
               
C.............wetted panels on face near T.E.........................
               N123=LCAVF
               DO 420 N=1,N123
                  N1=IFRST+N
                  SUM3=SUM3+BB(I1,N1)*DPDNC(N1)
 420           CONTINUE                  

C.............cavitating panels on face..............................
               DO 430 N=1,JCAVF
                  SUM2=SUM2+AA(I1,IFRST+M0MF-N)*PHI1(N,M,2)          
 430           CONTINUE
               
C.............fully wetted panels on the face near L.E. and fully....
C.............wetted panels on back before cavity L.E................
               IF(IFACE.EQ.1) THEN
                  NLAST=NC
               ELSE IF(IFACE.EQ.2) THEN
                  IF(JCAVB.EQ.0) THEN
                     NLAST=NC
                  ELSE
                     NLAST=M0MB-1
                  END IF
               END IF
               
               DO 440 N=M0MF,NLAST
                  N1=IFRST+N
                  SUM3=SUM3+BB(I1,N1)*DPDNC(N1)
 440           CONTINUE
               
C.............define terms associated with phi0 extrapolation........
               TERM3=TERM3-CT(M,4,2)*QC(1,M,2)*FCTF

C.............define terms associated with dipole extrapolation......
C.............left of split panel and source extrapolation right of..
C.............split panel on face....................................
               IF(NSPP(M,2).NE.0) THEN                  
                  TERM4=TERM4-AA(I1,IFRST+ISPF)*DLISP(M,2)*FLP(M,2)+
     *                 BB(I1,IFRST+ISPF)*QSPR(M,2)*FRP(M,2)
               END IF

C..........define RHS when face is fully wetted......................
            ELSE
               
C.............fully wetted panels on the face and fully wetted ......
C.............panels on back before cavity LE........................
               IF(IFACE.EQ.1) THEN
                  NLAST=NC
               ELSE IF(IFACE.EQ.2) THEN
                  IF(JCAVB.EQ.0) THEN
                     NLAST=NC
                  ELSE
                     NLAST=M0MB-1
                  END IF
               END IF

               DO 450 N=1,NLAST
                  N1=IFRST+N
                  SUM3=SUM3+BB(I1,N1)*DPDNC(N1)
 450           CONTINUE
            END IF

         END IF
         
         IF(IFACE.EQ.0.OR.IFACE.EQ.2) THEN
            
C..........define RHS when back is cavitating........................
            IF(JCAVB.GE.1) THEN
               
C.............wetted panels on face..................................
               IF(IFACE.EQ.0) THEN
                  DO 455 N=1,M0MB-1
                     N1=IFRST+N
                     SUM3=SUM3+BB(I1,N1)*DPDNC(N1)
 455              CONTINUE
               END IF

C.............cavitating panels on back..............................
               DO 460 N=1,JCAVB
                  SUM2=SUM2+AA(I1,IFRST+M0MB-1+N)*PHI1(N,M,1)        
 460           CONTINUE

C.............fully wetted panels on back near T.E...................
               DO 470 N=LCAVB,NC
                  N1=IFRST+N
                  SUM3=SUM3+BB(I1,N1)*DPDNC(N1)
 470           CONTINUE

C.............define terms associated with phi0 extrapolation........
               TERM3=TERM3-CT(M,4,1)*QC(1,M,1)*FCTB
               
C.............define terms associated with dipole extrapolation......
C.............left of split panel and source extrapolation right of..
C.............split panel on back....................................
               IF(NSPP(M,1).NE.0) THEN
                  TERM4=TERM4-AA(I1,IFRST+ISPB)*DLISP(M,1)*FLP(M,1)+
     *                 BB(I1,IFRST+ISPB)*QSPR(M,1)*FRP(M,1)
               END IF

C..........define RHS when back is fully wetted......................
            ELSE

C.............fully wetted panels on back............................
               IF(IFACE.EQ.0) THEN
                  DO 480 N=1,NC
                     N1=IFRST+N
                     SUM3=SUM3+BB(I1,N1)*DPDNC(N1)
 480              CONTINUE
               END IF
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
         
         IF(SOP(M,2).NE.ZERO) TERM2=TERM2+DUM1*SOP(M,2)*
     *        PHI1(JCAVF,M,2)
         IF(SOP(M,1).NE.ZERO) TERM2=TERM2-DUM1*SOP(M,1)*
     *        PHI1(JCAVB,M,1)
         
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
              END IF
            END IF
C........Corrected by Hong (HS071007)    
            
            IF(NNWC(M).GT.0) THEN
               N11=NNWC(M)+NSPS(M,NWDIR(M))+1
            ELSE
               N11=1
            END IF
            DO 485 N=N11,NTRA
               J=(MR-M)*NTRA+N
               TERM5=TERM5+C(I1,J)*SORW(J)
 485        CONTINUE
         END IF
         
C.......sum all terms to form RHS.................................
         RHS(IEQN)=RHS(IEQN)+TERM1+TERM2+TERM3+TERM4+TERM5
         
         IFRST1=IFRST1+NC-NSPP(M,2)-NSPP(M,1)+NNWC(M)
         
      END DO
      
      IF(IHUB.NE.0) THEN
C.......cal. terms to LHS that are associated with the hub...........
         DO 490 N=1,NHBX
            DO 500 M=1,MHBT
               J=INDEXH(N,M)
               J1=J-NPANB
               ALHS(IEQN,NBW+J1)=AA(I1,J)
 500        CONTINUE
         
C..........cal. terms to RHS that are associated with the hub........
            SUM1=ZERO
            DO 510 M=1,MHBT
               J=INDEXH(N,M)
               SUM1=SUM1+BB(I1,J)*DPDNC(J)
 510        CONTINUE
            RHS(IEQN)=RHS(IEQN)+SUM1
 490     CONTINUE
      END IF
      
      IF(IDUCT .NE. 0) THEN
         ICCOUNT = NBW + NPANH
C.......cal. terms to LHS that are associated with the Duct...........
         DO M=1,MDUCT
            DO N=1,NDUCT
               J=INDEXD(N,M)
               J1=J-NPANB-NPANH
               ALHS(IEQN,ICCOUNT+J1)=AA(I1,J)
            ENDDO
         
C..........cal. terms to RHS that are associated with the Duct........

            SUM1=ZERO
            DO N=1,NDUCT
               J=INDEXD(N,M)
               SUM1=SUM1+BB(I1,J)*DPDNC(J)
            ENDDO
            RHS(IEQN)=RHS(IEQN)+SUM1

            N1 = ICCOUNT + (M-1)*NDUCT + 1
            N2 = ICCOUNT + M*NDUCT

            IF(ISTEADY .EQ. 0) THEN
               ALHS(IEQN,N1) = ALHS(IEQN,N1) - WD(I1,M)
               ALHS(IEQN,N2) = ALHS(IEQN,N2) + WD(I1,M)
            ELSE
               ALHS(IEQN,N1)=ALHS(IEQN,N1)-
     *              (HALF*WD(I1,M)+WSUBIFD(I1,M))
               ALHS(IEQN,N2)=ALHS(IEQN,N2)+
     *              (HALF*WD(I1,M)+WSUBIFD(I1,M)) 
            ENDIF
         ENDDO
      END IF


      IF(ITUN .NE. 0) THEN
         ICCOUNT = NBW + NPANH +NPAND
C.......cal. terms to LHS that are associated with the tunnel...........
         DO N=1,NAXT
            DO M=1,MTUNEL
               J=INDEXTN(N,M)
               J1=J-NPANB-NPANH - NPAND
               ALHS(IEQN,ICCOUNT+J1)=AA(I1,J)
            ENDDO
         
C..........cal. terms to RHS that are associated with the TUNNEL........
            SUM1=ZERO
            DO M=1,MTUNEL
               J=INDEXTN(N,M)
               SUM1=SUM1+BB(I1,J)*DPDNC(J)
            ENDDO
            RHS(IEQN)=RHS(IEQN)+SUM1
         ENDDO
      END IF


c      IF(IAN.EQ.2) THEN
c         ICCOUNT=NBW+NPANH+NPAND+NPANTN
c         
cC.......cal. terms to LHS that are associated with the 
cC..........bulb on the blade tip
c         DO N=1,NTHX
c            DO M=1,MCVT
c               J=INDEXT(N,M)
c               J1=J - NPANB - NPANH - NPAND - NPANTN
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
c
c         ICCOUNT=NBW+NPANH+NPAND+NPANTN+NPANT
c
c         DO N=1,NCVX
c            DO M=1,MCVT
c               J=INDEXC(N,M)
c               J1=J - NPANB - NPANH - NPAND - NPANTN - NPANT
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
