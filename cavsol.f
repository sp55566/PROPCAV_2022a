      SUBROUTINE CAVSOL
************************************************************************
*                                                                      *
*  Author: Neal Fine  June 20, 1991                                    *
*                                                                      *
*  Date of last Revision                     Revision                  *
*  ---------------------                 ----------------              *
*  JYHL102401                  Subroutine undergone major revision.    *
*                                                                      *
************************************************************************
      USE MEMSOL
      USE CVRHS
      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'

C      COMMON/MEMSOL/ALHS(NTZ,NTZ)
C      COMMON/CVRHS/RHS(NTZ)

C.....Note: If ITMAX>100, You MUST set ITMX=ITMAX+1 in BLIC2.F(JY051800)
      ITMAX=100
      IPASS=1
      NORD=0

      NBLOCK = MR
      IF(IHUB .NE. 0) NBLOCK = NBLOCK + NHBX
      IF(IDUCT .NE. 0) NBLOCK = NBLOCK + MDUCT
      IF(ITUN .NE. 0) NBLOCK = NBLOCK + NAXT
c      IF(IAN .EQ. 2) NBLOCK = NBLOCK + NTHX + NCVX
      DO 10 M=1,MR
         NORD=NORD+NC+NNWC(M)-NSPP(M,2)-NSPP(M,1)
 10   CONTINUE

      IF(IHUB .NE. 0) NORD=NORD+NPANH
      IF(IDUCT .NE. 0) NORD=NORD+NPAND
      IF(ITUN .NE. 0) NORD = NORD + NPANTN
c      IF(IAN.EQ.2) NORD=NORD+NPANC+NPANT

      DO 20 M=1,MR
         NPERB(M)=NC-NSPP(MR-M+1,2)-NSPP(MR-M+1,1)+NNWC(MR-M+1)
 20   CONTINUE

      NNN = MR

      IF(IHUB .NE. 0) THEN
         DO N=1,NHBX
            NPERB(NNN+N)=MHBT
         ENDDO
         NNN = NNN + NHBX
      ENDIF

      IF(IDUCT .NE. 0) THEN
         DO N = 1, MDUCT
            NPERB(NNN+N) = NDUCT
         ENDDO
         NNN = NNN + MDUCT
      ENDIF

      IF(ITUN .NE. 0) THEN
         DO N = 1, NAXT
            NPERB(NNN+N) = MTUNEL
         ENDDO
         NNN = NNN + NAXT
      ENDIF

c      IF(IAN.EQ.2) THEN
c         DO N=1,NTHX
c            NPERB(NNN+N)=MCVT
c         END DO
c         NNN = NNN + NTHX 
c
c         DO N=1,NCVX
c            NPERB(NNN+N)=MCVT
c         END DO
c
c         NNN=NNN+NCVX
c      END IF

C.....Tolerance of the matrix solution is set to be 0.000005
      TOL=5.0E-06

C-----------------------------------------------------------------------
C     Call block iterative solver BLIC2
C-----------------------------------------------------------------------
      CALL BLIC2(SOL,RHS,TOL,NORD,NBLOCK,NPERB,ITMAX,IPASS,
     &           NTZ,NBLKMAX,101)

C-----------------------------------------------------------------------
C     dissect the results
C-----------------------------------------------------------------------
      IFRST=0

      DO 40 M=MR,1,-1

C.......indices for back cavity.........................................
         M0MB=M0(M,1)
         JCAVB=JCV(M,1)
         LCAVB=M0MB+JCAVB+NSPP(M,1)
         ISPB=LCAVB-1
               
C.......indices for face cavity.........................................
         M0MF=M0(M,2)
         JCAVF=JCV(M,2)
         LCAVF=M0MF-JCAVF-NSPP(M,2)-1
         ISPF=LCAVF+1

         IF(IFACE.EQ.1.OR.IFACE.EQ.2) THEN

C.......dipole and source strengths if face is cavitating...............
            IF(JCAVF.GE.1) THEN

C..........dipole strengths on wetted panels on face near T.E...........
               DO 50 N=1,LCAVF
                  J=INDEXB(N,M)
                  POT(J)=SOL(IFRST+N)
 50            CONTINUE
            
C..........source strengths beneath cavitating panels on face...........
               DO 60 N=1,JCAVF
                  N123 = M0MF-N
                  N456 = M0MF-N-NSPP(M,2)
                  J=INDEXB(N123,M)
                  DPDNC(J)=SOL(IFRST+N456)
 60            CONTINUE

C..........dipole strengths on wetted panels on face near L.E. and......
C..........wetted panels on back before cavity L.E......................
               IF(IFACE.EQ.1) THEN
                  NLAST=NC
               ELSE IF(IFACE.EQ.2) THEN
                  IF(JCAVB.EQ.0) THEN
                     NLAST=NC
                  ELSE
                     NLAST=M0MB-1
                  END IF
               END IF

               DO 70 N=M0MF,NLAST
                  N1=N-NSPP(M,2)
                  J=INDEXB(N,M)
                  POT(J)=SOL(IFRST+N1)
 70            CONTINUE

C..........calculate phi0 for cavity on the face........................
               N1=M0MF
               J1=INDEXB(N1,M)
               J2=INDEXB(N1+1,M)
               J3=INDEXB(N1+2,M)

               PHI0(M,2)=CT(M,1,2)*POT(J1)+CT(M,2,2)*POT(J2)+
     *              CT(M,3,2)*POT(J3)+CT(M,4,2)*QC(1,M,2)

C..........potential, POT=PHI0+PHI1, beneath cavity surface on face.....
               DO 80 N=1,JCAVF
                  N123 = M0MF-N
                  J=INDEXB(N123,M)
                  POT(J)=PHI0(M,2)+PHI1(N,M,2)
 80            CONTINUE

C..........terms associtated with split panel on face...................
               IF(NSPP(M,2).EQ.1)THEN

C.............dipole strength at left and right portion of split panel..
C.............on face...................................................
                  D1=0.
                  D2=0.
                  D3=0.
                  D4=0.

                  N00=ISPF-1
                  IF(N00.GT.0) D1=POT(INDEXB(N00,M))
                  IF(N00-1.GT.0) D2=POT(INDEXB(N00-1,M))
                  IF(N00-2.GT.0) D3=POT(INDEXB(N00-2,M))
                  IF(N00-3.GT.0) D4=POT(INDEXB(N00-3,M))
                  
                  PHIR(M,2)=DT(M,1,2)*D1+DT(M,2,2)*D2+DT(M,3,2)*D3+
     *                 DT(M,4,2)*D4
                  PHIL(M,2)=PHI0(M,2)+DLISP(M,2)

C.............dipole strength of split panel on face....................
                  POT(INDEXB(N00+1,M))=PHIL(M,2)*FLP(M,2)+
     *                 PHIR(M,2)*FRP(M,2)

C.............source strength at left side of split panel on face.......
                  Q1=0.
                  Q2=0.
                  Q3=0.
                  Q4=0.
            
                  N00=ISPF+1
                  IF(JCAVF.GE.1) Q1=DPDNC(INDEXB(N00,M))
                  IF(JCAVF.GE.2) Q2=DPDNC(INDEXB(N00+1,M))
                  IF(JCAVF.GE.3) Q3=DPDNC(INDEXB(N00+2,M))
                  IF(JCAVF.GE.4) Q4=DPDNC(INDEXB(N00+3,M))
                  
                  QSPL(M,2)=QT(M,1,2)*Q1+QT(M,2,2)*Q2+QT(M,3,2)*Q3+
     *                 QT(M,4,2)*Q4
                  
C.............panel number of split panel on face.......................
                  NISPF=INDEXB(ISPF,M)
                  
C.............source strength of split panel on face....................
                  DPDNC(NISPF)=QSPL(M,2)*FLP(M,2)+QSPR(M,2)*FRP(M,2)

               ENDIF

C.......dipole and source strengths if face is wetted...................
            ELSE

C..........dipole strengths on wetted panels on face and wetted panels..
C..........before cavity L.E. on back...................................
               IF(IFACE.EQ.1) THEN
                  NLAST=NC
               ELSE IF(IFACE.EQ.2) THEN
                  IF(JCAVB.EQ.0) THEN
                     NLAST=NC
                  ELSE
                     NLAST=M0MB-1
                  END IF
               END IF
                  
               DO 90 N=1,NLAST
                  J=INDEXB(N,M)
                  POT(J)=SOL(IFRST+N)
 90            CONTINUE
            
            END IF

         END IF

         IF(IFACE.EQ.0.OR.IFACE.EQ.2) THEN

C.......dipole and source strengths if back is cavitating...............
            IF(JCAVB.GE.1) THEN

C..........wetted panels on face........................................
               IF(IFACE.EQ.0) THEN
                  DO 95 N=1,M0MB-1
                     J=INDEXB(N,M)
                     POT(J)=SOL(IFRST+N)
 95               CONTINUE
               END IF

C..........source strengths beneath cavitating panels on back...........
               DO 100 N=1,JCAVB
                  N123 = M0MB-1+N
                  N456 = M0MB-1+N-NSPP(M,2)
                  J=INDEXB(N123,M)
                  DPDNC(J)=SOL(IFRST+N456)
 100           CONTINUE
            
C..........dipole strengths on wetted panels on back near T.E...........
               DO 110 N=LCAVB,NC
                  N1=N-NSPP(M,2)-NSPP(M,1)
                  J=INDEXB(N,M)
                  POT(J)=SOL(IFRST+N1)
 110           CONTINUE

C..........calculate phi0 for cavity on the back........................
               N1=M0MB-1
               J1=INDEXB(N1,M)
               J2=INDEXB(N1-1,M)
               J3=INDEXB(N1-2,M)

               PHI0(M,1)=CT(M,1,1)*POT(J1)+CT(M,2,1)*POT(J2)+
     *              CT(M,3,1)*POT(J3)+CT(M,4,1)*QC(1,M,1)

C..........potential, POT=PHI0+PHI1, beneath cavity surface on back.....
               DO 120 N=1,JCAVB
                  N123=M0MB-1+N
                  J=INDEXB(N123,M)
                  POT(J)=PHI0(M,1)+PHI1(N,M,1)
 120           CONTINUE
               
C..........terms associtated with split panel on back...................
               IF(NSPP(M,1).EQ.1)THEN

C.............dipole strength at left and right portion of split panel..
C.............on back...................................................
                  D1=0.
                  D2=0.
                  D3=0.
                  D4=0.

                  N00=LCAVB
                  IF(N00.LE.NC) D1=POT(INDEXB(N00,M))
                  IF(N00+1.LE.NC) D2=POT(INDEXB(N00+1,M))
                  IF(N00+2.LE.NC) D3=POT(INDEXB(N00+2,M))
                  IF(N00+3.LE.NC) D4=POT(INDEXB(N00+3,M))
                  
                  PHIR(M,1)=DT(M,1,1)*D1+DT(M,2,1)*D2+DT(M,3,1)*D3+
     *                 DT(M,4,1)*D4
                  PHIL(M,1)=PHI0(M,1)+DLISP(M,1)

C.............dipole strength of split panel on back....................
                  POT(INDEXB(N00-1,M))=PHIL(M,1)*FLP(M,1)+
     *                 PHIR(M,1)*FRP(M,1)

C.............source strength at left side of split panel on back.......
                  Q1=0.
                  Q2=0.
                  Q3=0.
                  Q4=0.
                  
                  N00=ISPB-1
                  IF(JCAVB.GE.1) Q1=DPDNC(INDEXB(N00,M))
                  IF(JCAVB.GE.2) Q2=DPDNC(INDEXB(N00-1,M))
                  IF(JCAVB.GE.3) Q3=DPDNC(INDEXB(N00-2,M))
                  IF(JCAVB.GE.4) Q4=DPDNC(INDEXB(N00-3,M))
               
                  QSPL(M,1)=QT(M,1,1)*Q1+QT(M,2,1)*Q2+QT(M,3,1)*Q3+
     *                 QT(M,4,1)*Q4

C.............panel number of split panel on back.......................
                  NISPB=INDEXB(ISPB,M)

C.............source strength of split panel on back....................
                  DPDNC(NISPB)=QSPL(M,1)*FLP(M,1)+QSPR(M,1)*FRP(M,1)

               END IF

C.......dipole and source strengths if back is wetted...................
            ELSE

C..........dipole strengths on wetted panels on back....................
               IF(IFACE.EQ.0) THEN
                  DO 130 N=1,NC
                     N1=N-NSPP(M,2)
                     J=INDEXB(N,M)
                     POT(J)=SOL(IFRST+N1)
 130              CONTINUE
               END IF
            
            END IF

         END IF

C.......potential, POT=PHI0+PHI1, beneath supercavity surface on wake...
         DO 140 N=1,NNWC(M)
            J=(MR-M)*NTRA+N
            POTW(J)=PHI0(M,NWDIR(M))+PHI1(JCV(M,NWDIR(M))+N,M,NWDIR(M))
 140     CONTINUE

C.......source strength of supercavitating wake panels..................

         DO 150 N=1,NNWC(M)
            J=(MR-M)*NTRA+N
            SORW(J)=SOL(IFRST+NC-NSPP(M,2)-NSPP(M,1)+N)
 150     CONTINUE

C.......source strength at split panel on wake..........................

C.HS071007..IF(NNWC(M).GT.0.AND.NSPS(M,NWDIR(M)).EQ.1) THEN 
         IF(NNWC(M).GT.0) THEN
           IF(NSPS(M,NWDIR(M)).EQ.1) THEN

            Q1=0.
            Q2=0.
            Q3=0.
            Q4=0.
            
            J=(MR-M)*NTRA+NNWC(M)
            IF(NNWC(M).GE.1) Q1=SORW(J)
            IF(NNWC(M).GE.2) Q2=SORW(J-1)
            IF(NNWC(M).GE.3) Q3=SORW(J-2)
            IF(NNWC(M).GE.4) Q4=SORW(J-3)
               
            QSSL=QW(M,1,NWDIR(M))*Q1+QW(M,2,NWDIR(M))*Q2+
     *           QW(M,3,NWDIR(M))*Q3+QW(M,4,NWDIR(M))*Q4
            
             IF(ISTEADY.EQ.0.OR.ISTEADY.EQ.1.OR.NREV.LE.2) THEN
               SORW(J+1)=QSSL*FLS(M,NWDIR(M))
             ELSE
               SORW(J+1)=QSSL*FLS(M,NWDIR(M))+QSSR(M)*FRS(M,NWDIR(M))
             END IF
           END IF
         ENDIF

         IFRST=IFRST+NC-NSPP(M,2)-NSPP(M,1)+NNWC(M)

 40   CONTINUE

C....Obtain potential for the hub from the solution
      IF(IHUB.NE.0) THEN
         ICC = NBW         
         DO NN=1,NHBX               
            DO MM=1,MHBT
               J=INDEXH(NN,MM)
               III = ICC + J - NPANB
               POT(J)=SOL(III)
            END DO
         END DO
      END IF

      IF(IDUCT.NE.0) THEN
         ICC = NBW + NPANH         
         DO MM=1,MDUCT              
            DO NN=1,NDUCT
               J=INDEXD(NN,MM)
               III = ICC + J - NPANB - NPANH
               POT(J)=SOL(III)
            END DO
         END DO
      END IF

      IF(ITUN .NE. 0) THEN
         ICC = NBW + NPANH + NPAND         
         DO NN=1,NAXT               
            DO MM=1,MTUNEL
               J=INDEXTN(NN,MM)
               J1 = J - NPANB - NPANH - NPAND
               POT(J)=SOL(ICC + J1)
            END DO
         END DO
      END IF

cC....Obtain potential for the bulb & tip vortex from the solution
c      IF(IAN.EQ.2) THEN
c         ICC = NBW + NPANH + NPAND + NPANTN
c         DO NN=1,NTHX
c            DO MM=1,MCVT
c               J=INDEXT(NN,MM)
c               J1 = J - NPANB - NPANH -NPAND - NPANTN
c               POT(J)=SOL(ICC + J1)
c            END DO
c         END DO
c
c         ICC=ICC + NPANT
c         DO NN=1,NCVX
c            DO MM=1,MCVT
c               J=INDEXC(NN,MM)
c               J1 = J -  NPANB - NPANH -NPAND - NPANTN - NPANT
c               POT(J)=SOL(ICC + J1)
c            END DO
c         END DO
c      END IF

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine CAVSOL>>>>>>>>>>>>>>>>>>>>>>>>>>
      END
