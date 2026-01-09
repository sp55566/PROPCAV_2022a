      SUBROUTINE CAVSOL_SC
************************************************************************
*                                                                      *
*  Author: Neal Fine  June 20, 1991                                    *
*                                                                      *
*  Date of last Revision                     Revision                  *
*  ---------------------                 ----------------              *
*  JY091801               Copied from cavsol.f.  Modified for ISC=1.   *
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

      IF(IAN.EQ.2) THEN
         NBLOCK=MR+NHBX+NTHX+NCVX
      ELSE
         NBLOCK=MR+NHBX
      END IF

      DO 10 M=1,MR
         NORD=NORD+NC+NNWC(M)-NSPP(M,2)-NSPP(M,1)
 10   CONTINUE

      if(ihub .ne. 0) NORD=NORD+NPANH
      IF(IAN.EQ.2) NORD=NORD+NPANC+NPANT

      DO 20 M=1,MR
         NPERB(M)=NC-NSPP(MR-M+1,2)-NSPP(MR-M+1,1)+NNWC(MR-M+1)
 20   CONTINUE

      if(ihub .ne. 0) then
         DO 30 M=1,NHBX
            NPERB(MR+M)=MHBT
 30      CONTINUE
      endif

c      IF(IAN.EQ.2) THEN
c         NNN=MR+NHBX
c         DO N=1,NTHX
c            NPERB(NNN+N)=MCVT
c         END DO
c
c         NNN=NNN+NTHX
c         DO N=1,NCVX
c            NPERB(NNN+N)=MCVT
c         END DO
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

C.......cavitating panels on face SR region
         DO N=1,N0(2)-1         
            J=INDEXB(N,M)
            DPDNC(J)=SOL(IFRST+N)
         END DO

C.......dipole and source strengths if face is cavitating...............
         IF(JCAVF.GE.1) THEN

C..........dipole strengths on wetted panels on face near T.E...........
            DO 50 N=N0(2),LCAVF
               J=INDEXB(N,M)
               POT(J)=SOL(IFRST+N)
 50         CONTINUE
            
C..........source strengths beneath cavitating panels on face...........
            DO 60 N=1,JCAVF
               N123 = M0MF-N
               N456 = M0MF-N-NSPP(M,2)
               J=INDEXB(N123,M)
               DPDNC(J)=SOL(IFRST+N456)
 60         CONTINUE

C..........dipole strengths on wetted panels on face near L.E. and......
C..........wetted panels on back before cavity L.E......................
            IF(JCAVB.EQ.0) THEN
               NLAST=N0(1)-1
            ELSE
               NLAST=M0MB-1
            END IF

            DO 70 N=M0MF,NLAST
               N1=N-NSPP(M,2)
               J=INDEXB(N,M)
               POT(J)=SOL(IFRST+N1)
 70         CONTINUE

C..........calculate phi0 for cavity on the face........................
            N1=M0MF
            J1=INDEXB(N1,M)
            J2=INDEXB(N1+1,M)
            J3=INDEXB(N1+2,M)

            PHI0(M,2)=CT(M,1,2)*POT(J1)+CT(M,2,2)*POT(J2)+
     *           CT(M,3,2)*POT(J3)+CT(M,4,2)*QC(1,M,2)

C..........potential, POT=PHI0+PHI1, beneath cavity surface on face.....
            DO 80 N=1,JCAVF
               N123 = M0MF-N
               J=INDEXB(N123,M)
               POT(J)=PHI0(M,2)+PHI1(N,M,2)
 80         CONTINUE

C..........if supercavitating
            IF(SOP(M,2).EQ.ONE) THEN
               PHI0T(M,2)=PSI0T(M,2)+PHI0(M,2)

C.............potential, POT=PHI0+PHI1, beneath cavity surface on face.....
               DO N=1,N0(2)-1                 
                  N123 = N0(2)-N
                  J=INDEXB(N,M)
                  POT(J)=PHI0T(M,2)+PHI2(N123,M,2)
               END DO
            END IF

C..........terms associtated with split panel on face...................
            IF(NSPP(M,2).EQ.1)THEN

C.............panel number of split panel on face.......................
               NISPF=INDEXB(ISPF,M)

C.............dipole strength at left and right portion of split panel..
C.............on face...................................................
               D1=0.
               D2=0.
               D3=0.
               D4=0.

               N0M1=N0(2)-1
               N00=ISPF-1
               IF(N00.GT.N0M1) D1=POT(INDEXB(N00,M))
               IF(N00-1.GT.N0M1) D2=POT(INDEXB(N00-1,M))
               IF(N00-2.GT.N0M1) D3=POT(INDEXB(N00-2,M))
               IF(N00-3.GT.N0M1) D4=POT(INDEXB(N00-3,M))
            
               PHIR(M,2)=DT(M,1,2)*D1+DT(M,2,2)*D2+DT(M,3,2)*D3+
     *              DT(M,4,2)*D4
               PHIL(M,2)=PHI0(M,2)+DLISP(M,2)

C.............dipole strength of split panel on face....................
               POT(NISPF)=PHIL(M,2)*FLP(M,2)+PHIR(M,2)*FRP(M,2)

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
     *              QT(M,4,2)*Q4

C.............source strength of split panel on face....................
               DPDNC(NISPF)=QSPL(M,2)*FLP(M,2)+QSPR(M,2)*FRP(M,2)

            ENDIF

C.......dipole and source strengths if face is wetted...................
         ELSE

C..........dipole strengths on wetted panels on face and wetted panels..
C..........before cavity L.E. on back...................................
            IF(JCAVB.EQ.0) THEN
               NLAST=N0(1)-1
            ELSE
               NLAST=M0MB-1
            END IF

            DO 90 N=N0(2),NLAST
               J=INDEXB(N,M)
               POT(J)=SOL(IFRST+N)
 90         CONTINUE
            
         END IF

C.......only if wetted or partially cavitating
         IF(SOP(M,2).EQ.ZERO) THEN
            N1=N0(2)
            J1=INDEXB(N1,M)
            J2=INDEXB(N1+1,M)
            J3=INDEXB(N1+2,M)
            PHI0T(M,2)=CT2(M,1,2)*POT(J1)+CT2(M,2,2)*POT(J2)+
     *           CT2(M,3,2)*POT(J3)+CT2(M,4,2)*QCSR(1,M,2)

C..........potential, POT=PHI0T+PHI2, beneath SR surface on face
            DO N=1,N0(2)-1                 
               N123 = N0(2)-N
               J=INDEXB(N,M)
               POT(J)=PHI0T(M,2)+PHI2(N123,M,2)
            END DO
         END IF

C.......cavitating panels on back SR region
         DO N=N0(1),NC
            N123=N-NSPP(M,2)-NSPP(M,1)
            J=INDEXB(N,M)
            DPDNC(J)=SOL(IFRST+N123)
         END DO

C.......dipole and source strengths if back is cavitating...............
         IF(JCAVB.GE.1) THEN

C..........source strengths beneath cavitating panels on back...........
            DO 100 N=1,JCAVB
               N123 = M0MB-1+N
               N456 = M0MB-1+N-NSPP(M,2)
               J=INDEXB(N123,M)
               DPDNC(J)=SOL(IFRST+N456)
 100        CONTINUE

C..........dipole strengths on wetted panels on back near T.E...........
            DO 110 N=LCAVB,N0(1)-1
               N1=N-NSPP(M,2)-NSPP(M,1)
               J=INDEXB(N,M)
               POT(J)=SOL(IFRST+N1)
 110        CONTINUE

C..........calculate phi0 for cavity on the back........................
            N1=M0MB-1
            J1=INDEXB(N1,M)
            J2=INDEXB(N1-1,M)
            J3=INDEXB(N1-2,M)

            PHI0(M,1)=CT(M,1,1)*POT(J1)+CT(M,2,1)*POT(J2)+
     *           CT(M,3,1)*POT(J3)+CT(M,4,1)*QC(1,M,1)

C..........potential, POT=PHI0+PHI1, beneath cavity surface on back.....
            DO 120 N=1,JCAVB
               N123=M0MB-1+N
               J=INDEXB(N123,M)
               POT(J)=PHI0(M,1)+PHI1(N,M,1)
 120        CONTINUE

C..........if supercavitating
            IF(SOP(M,1).EQ.ONE) THEN
               PHI0T(M,1)=PSI0T(M,1)+PHI0(M,1)

C.............potential, POT=PHI0+PHI1, beneath SR on back
               DO N=N0(1),NC
                  N123 = N-N0(1)+1
                  J=INDEXB(N,M)
                  POT(J)=PHI0T(M,1)+PHI2(N123,M,1)
               END DO
            END IF

C..........terms associtated with split panel on back...................
            IF(NSPP(M,1).EQ.1)THEN
C.............panel number of split panel on back.......................
               NISPB=INDEXB(ISPB,M)

C.............dipole strength at left and right portion of split panel..
C.............on back...................................................
               D1=0.
               D2=0.
               D3=0.
               D4=0.

               N0M1=N0(1)-1
               N00=LCAVB
               IF(N00.LE.N0M1) D1=POT(INDEXB(N00,M))
               IF(N00+1.LE.N0M1) D2=POT(INDEXB(N00+1,M))
               IF(N00+2.LE.N0M1) D3=POT(INDEXB(N00+2,M))
               IF(N00+3.LE.N0M1) D4=POT(INDEXB(N00+3,M))
            
               PHIR(M,1)=DT(M,1,1)*D1+DT(M,2,1)*D2+DT(M,3,1)*D3+
     *              DT(M,4,1)*D4
               PHIL(M,1)=PHI0(M,1)+DLISP(M,1)

C.............dipole strength of split panel on back....................
               POT(NISPB)=PHIL(M,1)*FLP(M,1)+
     *              PHIR(M,1)*FRP(M,1)

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
     *              QT(M,4,1)*Q4

C.............source strength of split panel on back....................
               DPDNC(NISPB)=QSPL(M,1)*FLP(M,1)+QSPR(M,1)*FRP(M,1)

            END IF

         END IF

C.......only if wetted or partially cavitating
         IF(SOP(M,1).EQ.ZERO) THEN
            N1=N0(1)-1
            J1=INDEXB(N1,M)
            J2=INDEXB(N1-1,M)
            J3=INDEXB(N1-2,M)
            PHI0T(M,1)=CT2(M,1,1)*POT(J1)+CT2(M,2,1)*POT(J2)+
     *           CT2(M,3,1)*POT(J3)+CT2(M,4,1)*QCSR(1,M,1)

C..........potential, POT=PHI0+PHI1, beneath cavity surface on face.....
            DO N=N0(1),NC
               N123 = N-N0(1)+1
               J=INDEXB(N,M)
               POT(J)=PHI0T(M,1)+PHI2(N123,M,1)
            END DO
         END IF

         DO 140 N=1,NNWC(M)
            J=(MR-M)*NTRA+N

C..........potential, POT=PHI0+PHI1, beneath supercavity surface on wake
            POTW(J)=PHI0T(M,NWDIR(M))+PHI2(NSR2+N,M,NWDIR(M))

C..........source strength of supercavitating wake panels
            SORW(J)=SOL(IFRST+NC-NSPP(M,2)-NSPP(M,1)+N)
 140     CONTINUE

C.......source strength at split panel on wake..........................
C.HS071007..IF(NNWC(M).GT.0.AND.NSPS(M,NWDIR(M)).EQ.1)THEN
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
         ICC=NBW         
         DO NN=1,NHBX               
            DO MM=1,MHBT
               J=INDEXH(NN,MM)               
               POT(J)=SOL(ICC+J-NPANB)
            END DO
         END DO
      END IF

cC....Obtain potential for the bulb & tip vortex from the solution
c      IF(IAN.EQ.2) THEN
c         ICC=NBW
c         DO NN=1,NTHX
c            DO MM=1,MCVT
c               J=INDEXT(NN,MM)
c               POT(J)=SOL(ICC+J-NPANB)
c            END DO
c         END DO
c         
c         DO NN=1,NCVX
c            DO MM=1,MCVT
c               J=INDEXC(NN,MM)
c               POT(J)=SOL(ICC+J-NPANB)
c            END DO
c         END DO
c      END IF

      RETURN
C<<<<<<<<<<<<<<<<<<<<<end of subroutine CAVSOL>>>>>>>>>>>>>>>>>>>>>>>>>>
      END
