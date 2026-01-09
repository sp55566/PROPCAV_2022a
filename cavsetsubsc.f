       SUBROUTINE CAVSETSUBSC(MM,NN,IEQN,I1)
************************************************************************
*                                                                      *
*  SET up the linear system of equations for solving the unsteady      *
*  cavitating propeller problem. Originally was setupb for the wing    *
*  problem. CAVSET computes the entire left-hand-side (which changes   *
*  with every timestep because cavity size changes) and the part of    *
*  the RHS which is associated with the cavity on the key blade. The   *
*  remainder of the RHS, namely the influence of the other blades and  *
*  the wakes of all the blades, is computed in CAVRHS.                 *
*                                                                      *
*  This subroutine is called inside CAVSET.  This sets up the equations*
*  for equation IEQN for the given column I1 = function (NN,MM).       *
*  This part is for Green's formula on the super-cavitating wake.      *
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
      
C       COMMON/MEMSOL/ALHS(NTZ,NTZ)
C       COMMON/CVRHS/RHS(NTZ)

       FOURPI=4.0*PI

C....IFRST1 marks the beginning of the unknowns on the foil strip.
C....within a column of the matrix................................
       IFRST1=0
      
       RHS(IEQN)=ZERO
       
       DO 540 M=MR,1,-1

C........Define WK1 as either WK or WKFACE depends if it's face 
C........or back supercavitation. (JY060100)
          IF(NWDIR(MM).EQ.1) THEN
             WK1=WK(I1,M)
          ELSE IF(NWDIR(MM).EQ.2) THEN
             WK1=WKFACE(I1,M,1)
          END IF

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

C..............wetted panels on face near T.E.........................
                N123=LCAVF
                DO 550 N=1,N123
                   ALHS(IEQN,IFRST1+N)=D(I1,IFRST+N)
 550            CONTINUE

C..............cavitating panels on face..............................
                DO 560 N=1,JCAVF
                   N123=M0MF-N
                   SUMF=SUMF+D(I1,IFRST+N123)
 560            CONTINUE

                DO 570 N=1,JCAVF
                   N123 = M0MF-N-NSPP(M,2)
                   N456 = M0MF-N
                   ALHS(IEQN,IFRST1+N123)=-E(I1,IFRST+N456)
 570            CONTINUE

C..............fully wetted panels on face near L.E. and fully wetted.
C..............panels on back before cavity L.E.......................
                IF(IFACE.EQ.1) THEN
                   NLAST=NC
                ELSE IF(IFACE.EQ.2) THEN
                   IF(JCAVB.EQ.0) THEN
                      NLAST=NC
                   ELSE
                      NLAST=M0MB-1
                   END IF
                END IF

                DO 580 N=M0MF,NLAST
                   N1=N-NSPP(M,2)
                   ALHS(IEQN,IFRST1+N1)=D(I1,IFRST+N)
 580            CONTINUE

C..............terms associated with phi0 extrapolation on face.......
                N1=IFRST1+M0MF-NSPP(M,2)
                N2=N1+1
                N3=N1+2

                IF(ISTEADY.EQ.0) THEN
                   IF(NWDIR(MM).EQ.1) THEN
                      FCTF=SUMF-WK(I1,M)*SOP(M,2)
                   ELSE IF(NWDIR(MM).EQ.2) THEN
                      FCTF=SUMF-WK2(I1,M)*SOP(M,2)
                   END IF
                ELSE
                   FCTF=SUMF-HALF*WK1*SOP(M,2)
                END IF

                IF(NSPP(M,2).NE.0) THEN
                   FCTF=FCTF+D(I1,IFRST+ISPF)*FLP(M,2)
                END IF

                IF(M.EQ.MM.AND.NWDIR(MM).EQ.2) FCTF=FCTF+FOURPI

                ALHS(IEQN,N1)=ALHS(IEQN,N1)+FCTF*CT(M,1,2)
                ALHS(IEQN,N2)=ALHS(IEQN,N2)+FCTF*CT(M,2,2)
                ALHS(IEQN,N3)=ALHS(IEQN,N3)+FCTF*CT(M,3,2)

                IF(NSPP(M,2).EQ.0) GO TO 585

C..............terms associated with dipole extrapolation right of....
C..............split panel on face....................................
                N1=IFRST1+ISPF-1
                N2=N1-1
                N3=N2-1
                N4=N3-1
                FCT=D(I1,IFRST+ISPF)*FRP(M,2)
                IF((ISPF-1).GT.0) ALHS(IEQN,N1)=ALHS(IEQN,N1)+
     *               DT(M,1,2)*FCT
                IF((ISPF-2).GT.0) ALHS(IEQN,N2)=ALHS(IEQN,N2)+
     *               DT(M,2,2)*FCT
                IF((ISPF-3).GT.0) ALHS(IEQN,N3)=ALHS(IEQN,N3)+
     *               DT(M,3,2)*FCT
                IF((ISPF-4).GT.0) ALHS(IEQN,N4)=ALHS(IEQN,N4)+
     *               DT(M,4,2)*FCT

C..............terms associated with source extrapolatin left of .....
C..............split panel on face....................................
                N1=IFRST1+ISPF
                N2=N1+1
                N3=N2+1
                N4=N3+1
                FCT=E(I1,IFRST+ISPF)*FLP(M,2)
                ALHS(IEQN,N1)=ALHS(IEQN,N1)-QT(M,1,2)*FCT
                ALHS(IEQN,N2)=ALHS(IEQN,N2)-QT(M,2,2)*FCT
                ALHS(IEQN,N3)=ALHS(IEQN,N3)-QT(M,3,2)*FCT
                ALHS(IEQN,N4)=ALHS(IEQN,N4)-QT(M,4,2)*FCT

 585            CONTINUE
                  
C...........define LHS when face is fully wetted......................
             ELSE
            
C..............fully wetted panels on face and fully wetted panels....
C..............before L.E. of cavity on back..........................
                IF(IFACE.EQ.1) THEN
                   NLAST=NC
                ELSE IF(IFACE.EQ.2) THEN
                   IF(JCAVB.EQ.0) THEN
                      NLAST=NC
                   ELSE
                      NLAST=M0MB-1
                   END IF
                END IF
                
                DO 590 N=1,NLAST
                   ALHS(IEQN,IFRST1+N)=D(I1,IFRST+N)
 590            CONTINUE
             END IF

          END IF

          IF(IFACE.EQ.0.OR.IFACE.EQ.2) THEN

C...........define LHS when back is cavitating........................
             IF(JCAVB.GE.1) THEN

C..............wetted panels on face..................................
                IF(IFACE.EQ.0) THEN
                   DO 595 N=1,M0MB-1
                      ALHS(IEQN,IFRST1+N)=D(I1,IFRST+N)
 595               CONTINUE
                END IF

C..............cavitating panels on back..............................
                DO 600 N=1,JCAVB
                   N123 = M0MB-1+N
                   SUMB=SUMB+D(I1,IFRST+N123)
 600            CONTINUE

                DO 610 N=1,JCAVB
                   N123 = M0MB-1+N-NSPP(M,2)
                   N456 = M0MB-1+N 
                   ALHS(IEQN,IFRST1+N123)=-E(I1,IFRST+N456)
 610            CONTINUE

C..............fully wetted panels on back near T.E...................
                DO 620 N=LCAVB,NC
                   N1=N-NSPP(M,2)-NSPP(M,1)
                   ALHS(IEQN,IFRST1+N1)=D(I1,IFRST+N)
 620            CONTINUE

C..............terms associated with phi0 extrapolation on back.......
                N1=IFRST1+M0MB-1-NSPP(M,2)
                N2=N1-1
                N3=N1-2

                IF(ISTEADY.EQ.0) THEN
                   IF(NWDIR(MM).EQ.1) THEN
                      FCTB=SUMB+WK(I1,M)*SOP(M,1)
                   ELSE IF(NWDIR(MM).EQ.2) THEN
                      FCTB=SUMB+WK2(I1,M)*SOP(M,1)
                   END IF
                ELSE
                   FCTB=SUMB+HALF*WK1*SOP(M,1)
                END IF

                IF(NSPP(M,1).NE.0) THEN
                   FCTB=FCTB+D(I1,IFRST+ISPB)*FLP(M,1)
                END IF

                IF(M.EQ.MM.AND.NWDIR(MM).EQ.1) FCTB=FCTB+FOURPI

                ALHS(IEQN,N1)=ALHS(IEQN,N1)+FCTB*CT(M,1,1)
                ALHS(IEQN,N2)=ALHS(IEQN,N2)+FCTB*CT(M,2,1)
                ALHS(IEQN,N3)=ALHS(IEQN,N3)+FCTB*CT(M,3,1)

                IF(NSPP(M,1).EQ.0) GO TO 625

C..............terms associated with dipole extrapolation right of....
C..............split panel on back....................................
                N1=IFRST1+ISPB+1-NSPP(M,1)-NSPP(M,2)
                N2=N1+1
                N3=N2+1
                N4=N3+1
                FCT=D(I1,IFRST+ISPB)*FRP(M,1)
                IF((ISPB+1).LE.NC) ALHS(IEQN,N1)=ALHS(IEQN,N1)+
     *               DT(M,1,1)*FCT
                IF((ISPB+2).LE.NC) ALHS(IEQN,N2)=ALHS(IEQN,N2)+
     *               DT(M,2,1)*FCT
                IF((ISPB+3).LE.NC) ALHS(IEQN,N3)=ALHS(IEQN,N3)+
     *               DT(M,3,1)*FCT
                IF((ISPB+4).LE.NC) ALHS(IEQN,N4)=ALHS(IEQN,N4)+
     *               DT(M,4,1)*FCT

C..............terms associated with source extrapolatin left of .....
C..............split panel on back....................................
                N1=IFRST1+ISPB-1-NSPP(M,2)
                N2=N1-1
                N3=N2-1
                N4=N3-1
                FCT=E(I1,IFRST+ISPB)*FLP(M,1)
                ALHS(IEQN,N1)=ALHS(IEQN,N1)-QT(M,1,1)*FCT
                ALHS(IEQN,N2)=ALHS(IEQN,N2)-QT(M,2,1)*FCT
                ALHS(IEQN,N3)=ALHS(IEQN,N3)-QT(M,3,1)*FCT
                ALHS(IEQN,N4)=ALHS(IEQN,N4)-QT(M,4,1)*FCT

 625            CONTINUE

C...........define LHS when back is fully wetted......................
             ELSE

C..............fully wetted panels on back............................
                IF(IFACE.EQ.0) THEN
                   DO 630 N=1,NC
                      N1=N-NSPP(M,2)
                      ALHS(IEQN,IFRST1+N1)=D(I1,IFRST+N)
 630               CONTINUE
                END IF
             END IF
                  
          END IF

C........terms associated with kutta conditon on foil..............
          N1=IFRST1+1
          N2=IFRST1+NC-NSPP(M,2)-NSPP(M,1)

          IF(ISTEADY.EQ.0) THEN
             IF(NWDIR(MM).EQ.1) THEN
                ALHS(IEQN,N1)=ALHS(IEQN,N1)-
     *               (ONE-SOP(M,2))*WK(I1,M)
                ALHS(IEQN,N2)=ALHS(IEQN,N2)+
     *               (ONE-SOP(M,1))*WK(I1,M)
             ELSE IF(NWDIR(MM).EQ.2) THEN
                ALHS(IEQN,N1)=ALHS(IEQN,N1)-
     *               (ONE-SOP(M,2))*WK2(I1,M)
                ALHS(IEQN,N2)=ALHS(IEQN,N2)+
     *               (ONE-SOP(M,1))*WK2(I1,M)
             END IF
          ELSE
             ALHS(IEQN,N1)=ALHS(IEQN,N1)-HALF*(ONE-SOP(M,2))*WK1
             ALHS(IEQN,N2)=ALHS(IEQN,N2)+HALF*(ONE-SOP(M,1))*WK1
          END IF

C........terms associated with supercavitating wake panels.........
          IF(NNWC(M).EQ.0) GO TO 650

          DO 640 N=1,NNWC(M)
             J=(MR-M)*NTRA+N
             J1=IFRST1+NC-NSPP(M,2)-NSPP(M,1)+N
             ALHS(IEQN,J1)=-F(I1,J)
 640      CONTINUE

          IF(NSPS(M,NWDIR(M)).EQ.0) GOTO 650

C........terms associated with source extrapolation of left split..
C........panel on supercavitating wake.............................
          J2=J1-1
          J3=J1-2
          J4=J1-3
          J=(MR-M)*NTRA+NNWC(M)+1
          ALHS(IEQN,J1)=ALHS(IEQN,J1)-F(I1,J)*QW(M,1,NWDIR(M))*
     *         FLS(M,NWDIR(M))
          ALHS(IEQN,J2)=ALHS(IEQN,J2)-F(I1,J)*QW(M,2,NWDIR(M))*
     *         FLS(M,NWDIR(M))
          ALHS(IEQN,J3)=ALHS(IEQN,J3)-F(I1,J)*QW(M,3,NWDIR(M))*
     *         FLS(M,NWDIR(M))
          ALHS(IEQN,J4)=ALHS(IEQN,J4)-F(I1,J)*QW(M,4,NWDIR(M))*
     *         FLS(M,NWDIR(M))

 650      CONTINUE

C........compute the right-hand-side due to the blade..............
          SUM2=ZERO
          SUM3=ZERO
          TERM1=ZERO
          TERM2=ZERO
          TERM3=ZERO
          TERM4=ZERO
          TERM5=ZERO

          IF(IFACE.EQ.1.OR.IFACE.EQ.2) THEN

C...........define RHS when face is cavitating........................
             IF(JCAVF.GE.1) THEN

C..............wetted panels on face near T.E.........................
                N123=LCAVF
                DO 660 N=1,N123
                   N1=IFRST+N
                   SUM3=SUM3+E(I1,N1)*DPDNC(N1)
 660            CONTINUE                  

C..............cavitating panels on face..............................
                DO 670 N=1,JCAVF
                   SUM2=SUM2+D(I1,IFRST+M0MF-N)*PHI1(N,M,2)           
 670            CONTINUE

C..............fully wetted panels on the face near L.E. and fully....
C..............wetted panels on back before cavity L.E................
                IF(IFACE.EQ.1) THEN
                   NLAST=NC
                ELSE IF(IFACE.EQ.2) THEN
                   IF(JCAVB.EQ.0) THEN
                      NLAST=NC
                   ELSE
                      NLAST=M0MB-1
                   END IF
                END IF

                DO 680 N=M0MF,NLAST
                   N1=IFRST+N
                   SUM3=SUM3+E(I1,N1)*DPDNC(N1)
 680            CONTINUE

C...............define terms associated with phi0 extrapolation........
                TERM3=TERM3-CT(M,4,2)*QC(1,M,2)*FCTF

C..............define terms associated with dipole extrapolation......
C..............left of split panel and source extrapolation right of..
C..............split panel on face....................................
                IF(NSPP(M,2).NE.0) THEN 
                   TERM4=TERM4-D(I1,IFRST+ISPF)*DLISP(M,2)*FLP(M,2)+
     *                  E(I1,IFRST+ISPF)*QSPR(M,2)*FRP(M,2)
                END IF

C...........define RHS when face is fully wetted......................
             ELSE
                  
C..............fully wetted panels on the face and fully wetted ......
C..............panels on back before cavity LE........................
                IF(IFACE.EQ.1) THEN
                   NLAST=NC
                ELSE IF(IFACE.EQ.2) THEN                   
                   IF(JCAVB.EQ.0) THEN
                      NLAST=NC
                   ELSE
                      NLAST=M0MB-1
                   END IF
                END IF

                DO 690 N=1,NLAST
                   N1=IFRST+N
                   SUM3=SUM3+E(I1,N1)*DPDNC(N1)
 690            CONTINUE
             END IF
             
          END IF

          IF(IFACE.EQ.0.OR.IFACE.EQ.2) THEN

C...........define RHS when back is cavitating........................
             IF(JCAVB.GE.1) THEN

C..............wetted panels on face..................................
                IF(IFACE.EQ.0) THEN
                   DO 695 N=1,M0MB-1
                      N1=IFRST+N
                      SUM3=SUM3+E(I1,N1)*DPDNC(N1)
 695               CONTINUE
                END IF

C..............cavitating panels on back..............................
                DO 700 N=1,JCAVB
                   SUM2=SUM2+D(I1,IFRST+M0MB-1+N)*PHI1(N,M,1)    
 700            CONTINUE

C..............fully wetted panels on back near T.E...................
                DO 710 N=LCAVB,NC
                   N1=IFRST+N
                   SUM3=SUM3+E(I1,N1)*DPDNC(N1)
 710            CONTINUE

C..............define terms associated with phi0 extrapolation........
                TERM3=TERM3-CT(M,4,1)*QC(1,M,1)*FCTB

C..............define terms associated with dipole extrapolation......
C..............left of split panel and source extrapolation right of..
C..............split panel on back....................................
                IF(NSPP(M,1).NE.0) THEN
                   TERM4=TERM4-D(I1,IFRST+ISPB)*DLISP(M,1)*FLP(M,1)+
     *                  E(I1,IFRST+ISPB)*QSPR(M,1)*FRP(M,1)
                END IF

C...........define RHS when back is fully wetted......................
             ELSE

C..............fully wetted panels on back............................
                IF(IFACE.EQ.0) THEN
                   DO 720 N=1,NC
                      N1=IFRST+N
                      SUM3=SUM3+E(I1,N1)*DPDNC(N1)
 720               CONTINUE
                END IF
             END IF

          END IF

C........sum the sources and dipoles on the non-split panels.......
          TERM1=-SUM2+SUM3

C........define term associateds with kutta condition on foil......
          IF(ISTEADY.EQ.0) THEN
             IF(NWDIR(MM).EQ.1) THEN
                DUM1=WK(I1,M)
             ELSE IF(NWDIR(MM).EQ.2) THEN
                DUM1=WK2(I1,M)
             END IF
          ELSE
             DUM1=HALF*WK1
          END IF

          IF(SOP(M,2).NE.ZERO) TERM2=TERM2+DUM1*SOP(M,2)*
     *         PHI1(JCAVF,M,2)
          IF(SOP(M,1).NE.ZERO) TERM2=TERM2-DUM1*SOP(M,1)*
     *         PHI1(JCAVB,M,1)

C........define term associated with trailing sinks due to.........
C........collapsed supercavities (JY060800)
          IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NREV.LE.2) THEN
             TERM5=ZERO
          ELSE
C.HS071007..IF(NNWC(M).GT.0.AND.NSPS(M,NWDIR(M)).EQ.1) THEN 
            IF(NNWC(M).GT.0) THEN
              IF(NSPS(M,NWDIR(M)).EQ.1) THEN
                J=(MR-M)*NTRA+NNWC(M)+1
                TERM5=TERM5+F(I1,J)*QSSR(M)*FRS(M,NWDIR(M))
              END IF
            END IF
C........Corrected by Hong (HS071007)   
             
             IF(NNWC(M).GT.0) THEN
                N11=NNWC(M)+NSPS(M,NWDIR(M))+1
             ELSE
                N11=1
             END IF
             DO 725 N=N11,NTRA
                J=(MR-M)*NTRA+N
                TERM5=TERM5+F(I1,J)*SORW(J)
 725         CONTINUE
          END IF

C........sum all terms to form RHS.................................
          RHS(IEQN)=RHS(IEQN)+TERM1+TERM2+TERM3+TERM4+TERM5

          IFRST1=IFRST1+NC-NSPP(M,2)-NSPP(M,1)+NNWC(M)

 540   CONTINUE


       RHS(IEQN)=RHS(IEQN)-
     *      FOURPI*PHI1(JCV(MM,NWDIR(MM))+NN,MM,NWDIR(MM))

       IF(IHUB.NE.0)THEN
      
C........cal. terms to LHS that are associated with the hub...........
          DO 730 N=1,NHBX
             DO 740 M=1,MHBT
                J=INDEXH(N,M)
                J1=J-NPANB
                ALHS(IEQN,NBW+J1)=D(I1,J)
               
 740         CONTINUE

C...........cal. terms to RHS that are associated with the hub.....
             SUM1=ZERO
             DO 750 M=1,MHBT
                J=INDEXH(N,M)
                SUM1=SUM1+E(I1,J)*DPDNC(J)
 750         CONTINUE
             RHS(IEQN)=RHS(IEQN)+SUM1
 730      CONTINUE
       ENDIF
      
      IF(IDUCT .NE. 0) THEN
         ICCOUNT=NBW+NPANH
         
C.......cal. terms to LHS that are associated with the duct

         DO N=1,NDUCT
            DO M=1,MDUCT
               J=INDEXD(N,M)
               J1=J - NPANB - NPANH
               ALHS(IEQN,ICCOUNT+J1)=D(I1,J)
            ENDDO
            
C..........cal. terms to RHS that are associated with the duct

            SUM1=ZERO
            DO  M=1,MDUCT
               J=INDEXD(N,M)
               SUM1=SUM1+E(I1,J)*DPDNC(J)
            ENDDO
            RHS(IEQN)=RHS(IEQN)+SUM1
         ENDDO


         DO M = 1, MDUCT
            N1 = ICCOUNT + (M-1)*NDUCT + 1
            N2 = ICCOUNT + M*NDUCT 

            IF(ISTEADY .EQ. 0) THEN
               ALHS(IEQN,N1)=ALHS(IEQN,N1)-WKD(I1,M)
               ALHS(IEQN,N2)=ALHS(IEQN,N2)+WKD(I1,M)
            ELSE
               ALHS(IEQN,N1)=ALHS(IEQN,N1)-HALF*WKD(I1,M)
               ALHS(IEQN,N2)=ALHS(IEQN,N2)+HALF*WKD(I1,M)
            END IF
         ENDDO
      ENDIF
     
      IF(ITUN .NE. 0) THEN
         ICCOUNT=NBW+NPANH+NPAND
         
C.......cal. terms to LHS that are associated with the duct
         
         DO N=1,NAXT
            DO M=1,MTUNEL
               J=INDEXTN(N,M)
               J1=J - NPANB - NPANH -NPAND
               ALHS(IEQN,ICCOUNT+J1)=D(I1,J)
            ENDDO
            
C..........cal. terms to RHS that are associated with the duct
            
            SUM1=ZERO
            DO  M=1,MTUNEL
               J=INDEXTN(N,M)
               SUM1=SUM1+E(I1,J)*DPDNC(J)
            ENDDO
            RHS(IEQN)=RHS(IEQN)+SUM1
         ENDDO         
      ENDIF

c      IF(IAN.EQ.2) THEN
c         ICCOUNT=NBW+NPANH+NPAND+NPANTN
c               
cC.......cal. terms to LHS that are associated with the 
cC..........bulb on the blade tip
c         DO N=1,NTHX
c            DO M=1,MCVT
c               J=INDEXT(N,M)
c               J1=J - NPANB - NPANH -NPAND - NPANTN
c               ALHS(IEQN,ICCOUNT+J1)=D(I1,J)
c            ENDDO
c            
cC..........cal. terms to RHS that are associated with the 
cC..........bulb on the blade tip
c            SUM1=ZERO
c            DO  M=1,MCVT
c               J=INDEXT(N,M)
c               SUM1=SUM1+E(I1,J)*DPDNC(J)
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
c               J1=J - NPANB - NPANH -NPAND - NPANTN - NPANT
c               ALHS(IEQN,ICCOUNT+J1)=D(I1,J)
c            ENDDO
c            
cC..........cal. terms to RHS that are associated with the 
cC..........tip vortex
c            SUM1=ZERO
c            DO  M=1,MCVT
c               J=INDEXC(N,M)
c               SUM1=SUM1+E(I1,J)*DPDNC(J)
c            ENDDO
c            RHS(IEQN)=RHS(IEQN)+SUM1
c         ENDDO
c         
c      END IF

      RETURN
      END
