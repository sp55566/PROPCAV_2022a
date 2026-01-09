      SUBROUTINE CAVBL2D_1

C**********************************************************************
C     Coupled with 2D Boundary Layer Analysis Code XFOIL               *
C     (Cavitating Case)                               *   
C     *
C     By Hong Sun                   April, 2006                    *
C**********************************************************************

      INCLUDE 'PUFCAV.INC'
      INCLUDE 'PUFCAVB.INC'
      INCLUDE 'PUFCAVC.INC'
      INCLUDE 'PUFBL.INC' 

      DIMENSION XN(NZVIS),YN(NZVIS),XC(NZVIS),YC(NZVIS),UINV(NZVIS)
     &  ,    XP(NZVIS),YP(NZVIS),DELS(NZVIS)
     &  ,    FW(NZVIS),ATEMP(NZVIS,NZVIS),BTEMP(NZVIS,NZVIS)
     &  ,    CPBLO(NBZ,MBZ)
      DIMENSION RV1(NZVIS),RRTEMP(NZVIS)
C     DIMENSION DELUE(NZVIS,MBZ)
      DIMENSION UEDGE(NZVIS), DSTARS(NZVIS), DSTARP(NZVIS)
     &  ,    THETAS(NZVIS), THETAP(NZVIS)
      REAL*8    UEDGE, DSTARS, DSTARP, THETAS, THETAP
      REAL*8    CDVIS
      DIMENSION UXIV(NZVIS),UETAV(NZVIS),UEDGV(NZVIS)
      INTEGER LVISCON 

      CALL OPFILEVIS
      
      CALL CLEAR(PPOTW, NSCWZ)
      CALL CLEAR(PPOTWS, NSCWZ)

C*************LOOP 1000 solves BL variables at each strip**************
      
      DO 1000 M = 1, MR

        CHDM = (CHORD(M)+CHORD(M+1))/2.0  
C     Reynolds number for each strip 
        IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
          REC(M) = REYD/2.0     ! based on D
        ELSE
          REC(M) = REYD*CHDM/ADVCO*SQRT(ADVCO**2+PI**2*RZP(M)**2)
        ENDIF      
C,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
C     Calculate the inflow angle of attack (pitch angle -inflow angle)
        
        CALL EVALDKs(NX,1,XR,RZP(M),VAR0,VACUB)  
        WROVS=PI*RZP(M)/ADVCO
        AOA = ATAN(0.5*(PITCH(M)+PITCH(M+1))/PI/RZP(M))
     &    - ATAN(VAR0/WROVS)
        ALPHA =AOA/RAD
        
        NTW(M) = NC+NTRA+NWMIN-NSUB 

C,,,,,,,,,,,,,,,,,,Calculate Geometric Information  ,,,,,,,,,,,,,,,,,,,

        DO 10 N = 1, NTW(M)+1
          XC(N) = XIVIS(N,M)
          YC(N) = ETAVIS(N,M)
 10     CONTINUE 
        
        DO 20 N = 1, NTW(M)
          XP(N) = (XC(N+1) + XC(N))/2.
          YP(N) = (YC(N+1) + YC(N))/2.0  
          DX1 = XC(N+1) - XC(N)
          DY1 = YC(N+1) - YC(N) 
          DS1 = SQRT(DX1*DX1 + DY1*DY1) 
          XN(N) = -1.0*DY1/DS1
          YN(N) =  DX1/DS1  
          DELTZ(N,M) = DS1             
 20     CONTINUE

        DO 30 N = 2, NTW(M)
          CALL ARCLEN( SS1,SS2,XC(N-1),YC(N-1),XC(N),YC(N)
     &      ,XC(N+1),YC(N+1),DELTZ(N-1,M),DELTZ(N,M) )
          IF(N.EQ.2) THEN                                                  
            DELTS(1,M)=SS1
          ELSE                                                  
            DELTS(N-1,M)=HALF*(SS1+SPRE)  
          END IF                                                
          SPRE=SS2                                                         
 30     CONTINUE
        DELTS(NTW(M),M)=SS2 
        
        SC(1,M) = 0.0
        DO 40 N = 1, NTW(M)       
          SC(N+1,M) = SC(N,M) + DELTS(N,M)
          SP(N,M) = SC(N,M) + HALF*DELTS(N,M)
          DELS(N) = DELTS(N,M)       
 40     CONTINUE
        
        
C'''''''''''''Calculate Matrices for each strip ''''''''''''''''''''''
C     write(*,*) 'Calculate Matrices for each strip'

        CALL SETUP1(NC,NC,NC,BTEMP,ATEMP,XC,YC,XP,YP,XN,YN,
     &    DELS,NZVIS)
        
C     FW is the wake influence matrix
        DO N = 1, NC
          IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
            FW(N) = - ALPHA*RAD - ATAN2( YP(N) - YC(1), XC(1) - XP(N) ) 
          ELSE
            FW(N) = - ATAN2( YP(N)-YC(1), XC(1)-XP(N) ) - AOA
          ENDIF
        ENDDO

        DO IM = 1, NC
          DO IN = 1, NC
            IF(IM.EQ.1) THEN
              AV(IN,IM) = ATEMP(IN,IM) - FW(IN)
            ELSE IF(IM.EQ.NC) THEN
              AV(IN,IM) = ATEMP(IN,IM) + FW(IN)              
            ELSE
              AV(IN,IM) = ATEMP(IN,IM)
            ENDIF
          END DO
        END DO                  ! AV

        CALL INVRSE(AV,AVINV,NC,NZVIS) !AVINV

        NWv = NTW(M)-NC

        CALL SETUP2(NWv,NC,NC,CS,AWV,XC,YC,XP,YP,XN,YN,
     &    DELS,NZVIS)           !AWV

        CALL SETUP1(NC,NTW(M),NC,BETA,ATEMP,XC,YC,XP,YP,XN,YN,
     &    DELS,NZVIS)           !BETA

        CALL MAT2DOT(AVINV,BETA,CA,NC,NZVIS,NC,NZVIS,NTW(M),NZVIS)
        
        I1=NC/4
        I2=3*I1      
        NSTAG=I1
        CPMIN=1.0
        DO I=I1,I2
          CP = -CPBN(I,M)*ADVCO**2.
C     CP = UINV(I)**2-1.
          IF(CP.LT.CPMIN)THEN
            NSTAG = I
            CPMIN = CP
          ENDIF
        ENDDO
C     CPM = UINV(NSTAG-1)**2-1.
C     CPP = UINV(NSTAG+1)**2-1.
C     IF(CPM.LT.CPP) NSTAG=NSTAG-1
C     write(*,*)'m=','nstag=',nstag,'cpmin=',cpmin
C     NSTAG=NC/2+1      
        CPM = -CPBN(NSTAG-1,M)*ADVCO**2.
        CPP = -CPBN(NSTAG+1,M)*ADVCO**2.
        IF(CPM.LT.CPP) NSTAG=NSTAG-1
        
        DO 140 I=2,NC
          DO 130 J=1, NTW(M)
            DH(I,J)=(CA(I,J)-CA(I-1,J))/(SP(I,M)-SP(I-1,M))
            IF(I.LT.NSTAG)THEN
              DH(I,J)=-DH(I,J) 
            ENDIF
 130      CONTINUE
 140    CONTINUE

        CALL SETUP2(NWv,NTW(M),NC,CSIGMA,ATEMP,XC,YC,XP,YP,XN,YN,
     &    DELS,NZVIS)           !CSIGMA      

        CALL MAT2DOT(AWV,AVINV,CB,NWv,NWZ,NC,NBZ,NC,NBZ)

        CALL MAT2DOT(CB,BETA,DH1,NWv,NWZ,NC,NBZ,NTW(M),NZVIS)

        DO 160 I=1,NWv
          DO 150 J=1,NTW(M)
            DH1(I,J)=(DH1(I,J)+CSIGMA(I,J))/2./PI
 150      CONTINUE
 160    CONTINUE
        
        DO 180 I=2,NWv
          DO 170 J=1,NTW(M)
            DH(NC+I,J)=(DH1(I,J)-DH1(I-1,J))/(SP(NC+I,M)-SP(NC+I-1,M))
 170      CONTINUE
 180    CONTINUE

        DO 190 J=1,NTW(M)
          DH(1,J)   = DH(2,J)+(DH(2,J)-DH(3,J))*(SC(1,M)-SC(2,M))
     &      / (SC(2,M)-SC(3,M))
          DH(NC+1,J)= DH(NC,J)+(DH(NC,J)-DH(NC-1,J))
     &      * (SC(NC+1,M)-SC(NC,M))/(SC(NC,M)-SC(NC-1,M))
          IF(J.EQ.1)   DH(1,J)=-DH(1,J)
          IF(J.EQ.NC)  DH(NC+1,J)=-DH(NC+1,J)
          IF(J.LT.NSTAG) DH(NC+1,J)=-DH(1,J)
          IF(J.GE.NSTAG) DH(1,J)=-DH(NC+1,J)
 190    CONTINUE
c==   %==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%=
c-----------------check symmetry of H matrix-------------------------
C     CALL MAT2DOT(AVINV,AV,BTEMP,NC,NZVIS,NC,NZVIS,NC,NZVIS)
C     IF(M.EQ.1) THEN
C     WRITE(126,*) 'I,J,AV,AVINV,ATEMP,CB,BTEMP'
C     DO I = 1, NC
C     DO J = 1, NC
C     WRITE(126,'(2I6,5F10.4)')I,J,AV(I,J),AVINV(I,J),
C     &       ATEMP(I,J),CB(I,J),BTEMP(I,J)
C     ENDDO 
C     ENDDO
C     
C     WRITE(127,*) 'I,J,DH1(I,J)'
C     DO I = 1, NWv
C     DO J = 1, NTW(M)
C     IF(J.LE.NC) THEN
C     WRITE(127,'(2I6,4F10.4)') I,J,DH1(I,J),CSIGMA(I,J)
C     &                                  ,AWV(I,J)
C     ELSE
C     WRITE(127,'(2I6,4F10.4)') I,J,DH1(I,J),CSIGMA(I,J),0.
C     ENDIF
C     ENDDO 
C     ENDDO
C     
C     WRITE(128,*) 'I,J,CA(I,J),beta(i,j)'
C     DO I = 1, NC
C     DO J = 1, NTW(M)
C     WRITE(128,'(2I6,2F10.4)') I,J,CA(I,J),beta(i,j)
C     ENDDO 
C     ENDDO
C     ENDIF

c-------------------check symmetry of H matrix------------------------
        
        DO 210 I=1,NTW(M)-1
          DO 200 J=2,NTW(M)-1
            DD(I,J) = DH(I,J-1)/(SC(J,M)-SC(J-1,M)) 
     &        - DH(I,J)/(SC(J+1,M)-SC(J,M))
            IF(I.LT.NSTAG.AND.J.LE.NC+1) THEN
              DD(I,J) = -DD(I,J)
            ENDIF
            IF(J.GT.NC+1) DD(1,J)=-DD(1,J) 
 200      CONTINUE
 210    CONTINUE

        DO 220 I=1,NTW(M)-1
          IF(I.GT.NC+1)THEN
            DD(I,1)=DD(I,NC+1)
          ELSE
            DD(I,1)=DD(NC+2-I,NC+1)
          ENDIF
 220    CONTINUE

C'''''''''''''''''Calculate UINV  for each strip '''''''''''''''''''''
        write(*,*)'Calculate UINV  for each strip' 

C.......CALCULATE PERTUBATION POTENTIAL IN THE WAKE

        JM = M
        CALL VISPOT2(JM)

        DO 50 N = 2, NC
          N1 = N-1
          L  = INDEXB(N,M)
          L1 = INDEXB(N1,M)
          UXI = DPDUB(N,M)+VXIB(N,M)
          UXIM = DPDUB(N1,M)+VXIB(N1,M)
C     VXI1 = VXIB(N,M)
C     VXI2 = VXIB(N1,M)
C     VXIA = (VXI1+VXI2)/2.0 
C     UINV(N)= (POT(L)-POT(L1))/(SP(N,M)-SP(N-1,M))-VXIA
          UINV(N)= -(UXI+UXIM)/2.0
          VINFBA = SQRT((VINFSB(N,M)+VINFSB(N1,M))/2.0)
          UINV(N) = UINV(N)/VINFBA
 50     CONTINUE

        write(138,*) 'Zone T="M= ', M, ' " '
        DO 60 IN = 2,NTRA
          L = NTRA*(MR-M)+IN
          L1 = NTRA*(MR-M)+IN-1
          VXI1 = ABS(VXIWs(IN,M))
          VXI2 = ABS(VXIWs(IN-1,M))
          VXIWA = (VXI1+VXI2)/2.0    
          UINV(NC+IN) = (PPOTWs(L)-PPOTWs(L1)) / 
     &      (SP(NC+IN,M)-SP(NC+IN-1,M))+VXIWA
          VINFWsA = SQRT((VINFSWs(IN,M)+VINFSWs(IN-1,M))/2.0)
          UINV(NC+IN) = UINV(NC+IN)/VINFWsA
          write(138,*) NC+IN,PPOTWs(L1) 
 60     CONTINUE
        
        write(139,*) 'Zone T="M= ', M, ' " '
        DO 70 IN = NSUB+2, NWMIN
          L = INDEXW(IN,M) 
          L1 = INDEXW(IN-1,M)
          VXI1 = ABS(VXIW(IN,M))
          VXI2 = ABS(VXIW(IN-1,M))
          VXIWA = (VXI1+VXI2)/2.0     
          UINV(NC+IN+NTRA-NSUB) = (PPOTW(L)-PPOTW(L1)) / 
     &      (SP(NC+IN+NTRA-NSUB,M)-SP(NC+IN-1+NTRA-NSUB,M))+VXIWA
          VINFWA = SQRT((VINFSW(IN,M)+VINFSW(IN-1,M))/2.0)
          UINV(NC+IN+NTRA-NSUB) = UINV(NC+IN+NTRA-NSUB)/VINFWA
          write(139,*) NC+IN+NTRA-NSUB,PPOTW(L1) 
 70     CONTINUE

        IF(NWMIN.GT.NSUB) THEN 
          L2 = NTRA*(MR-M+1)
          L3 = INDEXW(NSUB+1,M)
          VXI1 = ABS(VXIW(NSUB+1,M))
          VXI2 = ABS(VXIW(NSUB,M))
          VXIWA = (VXI1+VXI2)/2.0  
          UINV(NC+NTRA+1) = (PPOTW(L3)-PPOTWs(L2)) / 
     &      (SP(NC+NTRA+1,M)-SP(NC+NTRA,M))+VXIWA
          VINFWA = SQRT((VINFSW(NSUB+1,M)+VINFSW(NSUB,M))/2.0)
          UINV(NC+NTRA+1) =  UINV(NC+NTRA+1)/VINFWA 
        ENDIF
        
c     UINV(1) = UINV(2) + (SC(1,M)-SC(2,M))*
c     &              (UINV(2)-UINV(3))/(SC(2,M)-SC(3,M))
c     UINV(NC+1) = - UINV(1)
        UINV(1) = -UINV(NC+2)
        UINV(NC+1) = UINV(NC+2)      

C*********************PUT TO VISCAL CONVENTIONS   *****************
C     NON-DIMENSIONALIZATION OF THE MATRICIES BY CHORD/PROPELLER RADIUS
C**********************************************************************
        
C     NON-DIMENSIONALIZATION 
c     CHDVIS = 2.*XCVIS(1)
        DO I = 1, NTW(M)
c     XCVIS(I)= XCVIS(I)/CHDVIS + 0.5
c     YCVIS(I)= YCVIS(I)/CHDVIS
c     SVL(I)= SVL(I)/CHDVIS
          XC(I)= XC(I)/CHDM/2. + 0.5
          YC(I)= YC(I)/CHDM/2.
          SC(I,M)= SC(I,M)/CHDM/2.        

C     DO J=1,NTW(M)
C     DD(I,J)= DD(I,J)*CHDVIS
C     ENDDO
        ENDDO

        WRITE(*,*) 'PUT TO VISCAL SIGN CONVENTIONS '

C     PUT TO VISCAL SIGN CONVENTIONS  
        DO 230 I=1,NTW(M)-1
          IF(I.LE.NC+1) THEN
            UVL(I)= UINV(NC+2-I)
            SVL(I)= SC(NC+1,M)-SC(NC+2-I,M)
            XNVIS(I)= XN(NC+2-I)
            YNVIS(I)= YN(NC+2-I)    
            XCVIS(I)= XC(NC+2-I)
            YCVIS(I)= YC(NC+2-I)
          ELSE
            UVL(I)= -UINV(I)
            SVL(I)= SC(I,M)
            XNVIS(I)= XN(I)
            YNVIS(I)= YN(I)    
            XCVIS(I)= XC(I)
            YCVIS(I)= YC(I)
          ENDIF
 230    CONTINUE


C     CALCULATE DVL 
        DO 250 I = 1,NTW(M)-1
          DO 240 J = 1,NTW(M)-1
            IF(I.LE.NC+1.AND.J.LE.NC+1)THEN
              DVL(I,J)=DD(NC+2-I,NC+2-J)
            ELSE IF(I.LE.NC+1.AND.J.GT.NC+1)THEN
              DVL(I,J)=DD(NC+2-I,J)
            ELSE IF(I.GT.NC+1.AND.J.LE.NC+1)THEN
              DVL(I,J)=DD(I,NC+2-J)
            ELSE IF(I.GT.NC+1.AND.J.GT.NC+1)THEN
              DVL(I,J)=DD(I,J)
            ENDIF
 240      CONTINUE
 250    CONTINUE

C     PUT TO VISCAL SIGN CONVENTION 
        DO 270 I=1,NTW(M)
          DO 260 J=1,NTW(M)
            IF(I.NE.NC+1) DVL(I,J)=-DVL(I,J)
            IF(I.GT.NC+1.AND.J.LE.NC) DVL(I,J)=-DVL(I,J)
            IF(J.GT.NC+1.AND.I.LE.NSTAG) DVL(I,J)=-DVL(I,J)
            IF(J.GT.NC+1.AND.I.EQ.NC+1) DVL(I,J)=-DVL(I,J)
 260      CONTINUE
 270    CONTINUE

        XCVIS(NC+2) = XCVIS(NC+1)
        YCVIS(NC+2) = YCVIS(NC+1)
        SVL(NC+2) = SVL(NC+1)
        UVL(NC+2) = UVL(NC+1)

        DO 280 J=1,NTW(M)
          DVL(NC+2,J)= DVL(NC+1,J)
 280    CONTINUE
        DVL(1,NC+1)=-DVL(1,NC+1)
        DVL(NC+1,NC+1)=-DVL(NC+1,NC+1)
        DVL(NC+2,NC+1)=-DVL(NC+2,NC+1)

        DO 290 J=NSTAG+1,NC
          DVL(1,J)=-DVL(1,J)
          DVL(NC+1,J)=-DVL(NC+1,J)
          DVL(NC+2,J)=-DVL(NC+2,J)
 290    CONTINUE

C-------------------------Code Test -----------------------------
        IREAD = 0               !!!
        if (IREAD.EQ.1) then
          NWv = 30
          open(123, file='121.dat')
          read(123,*)
          do I = 1, NC+NWv-1
            read(123,*) IA,UVL(I),SVL(I),XNVIS(I),YNVIS(I)
     &        ,          XC(I),YC(I)
          enddo 
          close(123)

          open(124, file='122.dat')
          DO I = 1, NC+NWv
            READ(124,293) (DVL(I,J),J=1,NC+NWv)
          ENDDO 
          close(124)
        endif
C---------------------------------------------------------------

        WRITE(121,*)'Zone T= "r/R=',RZP(M),', AOA=', AOA/RAD, 
     &    ', RE=', REC(M), '"'
        DO I = 1, NC+NWv-1
          WRITE(121,292) I,UVL(I),SVL(I),XNVIS(I),YNVIS(I)
     &      ,          XCVIS(I),YCVIS(I)
        ENDDO
 292    FORMAT(1X,I3,7F12.6)   

C     WRITE(122,*) 'Zone T= "M = ', M, '"'
C     WRITE(122,*) 'RE = ', REC(M), 'r/R = ', RZP(M)
C     DO I = 1, NC+NWv
C     WRITE(122,293) (DVL(I,J),J=1,NC+NWv)
C     ENDDO
 293    FORMAT(1X,7F10.4)
        
c==   %==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%==%=
C     IF(M.EQ.1) THEN
C     WRITE(125,*) 'DH(I,J),DD(I,J),DVL(I,J),DVL(J,I)'
C     DO I = 1, NTW(M)
C     DO J = 1, NTW(M)
C     WRITE(125,221) DH(I,J),DD(I,J),DVL(I,J),DVL(J,I)
C     ENDDO 
C     ENDDO
C     221    FORMAT(1X,4F10.4)
C     ENDIF

c     GOTO 1000
C     VISCAL ANALYSIS AT EACH STRIP
C**********************************************************************
        WRITE(10,*) 'Zone T= "r/R = ', RZP(M), '"'
        WRITE(11,*) 'Zone T= "r/R = ', RZP(M), '"'
        WRITE(12,*) 'Zone T= "r/R = ', RZP(M), '"'
c     WRITE(14,*) 'Zone T= "r/R = ', RZP(M), '"'
c     WRITE(20,*) 'Zone T= "r/R = ', RZP(M), '"'
c     WRITE(97,*) 'Zone T= "r/R = ', RZP(M), '"'   

C     WRITE(*,*) 'BEFORE CALL VISCAL ROUTINE'
        NZ = NZVIS
        WRITE(*,*) 'M=',M,'r/R=',RZP(M), 'NSTAG=',NSTAG
        
C     ILEP,ITEP: Use XFoil Index Numbering (anti-clockwise)
        IILEP = M0(M,1)
        IITEP = LCV(M,1)
        write(*,*) 'IILEP=', IILEP, 'IITEP=', IITEP
        
        IF(IILEP.NE.IILEP) THEN
          IF(IILEP.LE.NC+1) THEN
            ILEP = NC+2-IILEP
          ELSE 
            ILEP = IILEP
          ENDIF
          IF(IITEP.LE.NC+1) THEN
            ITEP = NC+2-IITEP  
          ELSE
            ITEP = IITEP
          ENDIF
          ICV=2
          
          IF(NNWC(M).GT.0) THEN
            ICV=3
            ITEP = NC+1
          ENDIF
        ENDIF             

        XTRANS = XTRANLS(M)
        XTRANP = XTRANLP(M)     
        WRITE(*,*) 'XTRANS =', XTRANS, 'XTRANP =', XTRANP
        WRITE(*,*) 
        
        write(1236,*) 'strip ', M
        CALL VISCAL(NC+1,NWv-2,SVL,XCVIS,YCVIS,XTRANS,XTRANP,REC(M),
     &    ALPHA,UVL,DVL,NZ,UEDGE,DSTARS,DSTARP,THETAS,THETAP,ICV,
     &    ILEP,ITEP,ICLE,ICTE,XNVIS,YNVIS,RVCRIT,MAXIT,EPS1,CDVIS,
     &    LVISCON, BETA_BL_INP) 
        
        LVFLAG(M) = LVISCON      

        write(*,*) 'ICLE=',ICLE,'ICTE=',ICTE,'LVFLAG=',LVFLAG(M)

        CDV(M)= SNGL(CDVIS)
        IF(LVFLAG(M).EQ.1)  CDV(M)=XCDF
        WRITE(22,'(I5,2F12.8)') M,  RZP(M), CDV(M)

C**********************************************************************
C     CALCULATE 3-D VISCOUS PRESSURE DISTRIBUTION AT PANEL CENTROID
C     BY USING VISCOUS EDGE VELOCITY

C     
C     PUT BACK TO PROPELLER COORDINATE SYSTEM

        DO I = 1,NC+1
          UEDGV(I)=SNGL(UEDGE(NC+2-I))
        ENDDO

C     FIND DU CORRECTION AT BLADE T.E.(REFER TO HUFFORD, 1992)

        L1 = INDEXB(1,M)
        LNC = INDEXB(NC,M)
        UXIL = ( ABS(UEDGV(2))+ABS(UEDGV(1)) )/2.
        UXIU = ( ABS(UEDGV(NC+1))+ABS(UEDGV(NC)) )/2.
        UXIL = UXIL*SQRT(VINFSB(1,M)) - VXIB(1,M)
        UXIU = UXIU*SQRT(VINFSB(NC,M)) - VXIB(NC,M)
        ULOWER = ABS(VXIB(1,M)+UXIL)
        UUPPER = ABS(VXIB(NC,M)+UXIU)
        UBLTE = 0.5*(UUPPER+ULOWER)
        UETAL = (DPDVB(1,M)-DPDUB(1,M)*SINPHI(L1))/COSPHI(L1)
        UETAU = (DPDVB(NC,M)-DPDUB(NC,M)*SINPHI(LNC))/COSPHI(LNC)
        WLOWER = VETAB(1,M)+UETAL
        WUPPER = VETAB(NC,M)+UETAU
        DU = (WLOWER**2-WUPPER**2)/(4.*UBLTE)

        DO I =1, NC
          L = INDEXB(I,M)
          UXIV(I) = ( ABS(UEDGV(I+1))+ABS(UEDGV(I)) )/2.
          UXIV(I) = UXIV(I)*SQRT(VINFSB(I,M))-VXIB(I,M)
          UETAV(I) = (DPDVB(I,M)-DPDUB(I,M)*SINPHI(L))/COSPHI(L)
          IF(I.EQ.1) THEN
            UXITV = UBLTE - DU 
          ELSE IF(I.EQ.NC) THEN
            UXITV = UBLTE + DU           
          ELSE
            UXITV = VXIB(I,M)+UXIV(I)
          END IF
          UETATV = VETAB(I,M)+UETAV(I)        
          
          VTOTS_V(L) = UXITV**2 + UETATV**2
          CPBL(I,M)=CPBN(I,M)+VTOTS(L)-VTOTS_V(L)
          CPBLO(I,M) = CPB(I,M)+VTOTS(L)-VTOTS_V(L)           

          UXTOT_V(I,M)=UXITV*DIR(L,1,1)+UETATV*DIR(L,2,1)
          UYTOT_V(I,M)=UXITV*DIR(L,1,2)+UETATV*DIR(L,2,2)
          UZTOT_V(I,M)=UXITV*DIR(L,1,3)+UETATV*DIR(L,2,3)
          
          IF(LVFLAG(M).EQ.1) THEN   
            CPBL(I,M)=CPBN(I,M)
            CPBLO(I,M) = CPB(I,M)
            VTOTS_V(L) = VTOTS(L)           
            UXTOT_V(I,M)=UXTOT(I,M)
            UYTOT_V(I,M)=UYTOT(I,M)
            UZTOT_V(I,M)=UZTOT(I,M)
          ENDIF         
        ENDDO

        WRITE(20,*)'ZONE T="M=', M,' T=', ITSTEP,'"'    
C.......Writing wetted pressure for pressure side.
        DO NN=1,NC/2
c     WRITE(20,*) SBP(NC/2-NN+1),-CPBL(NN,M)*ADVCO**2.
          WRITE(20,*) SBP(NC/2-NN+1),-CPBLO(NN,M)*ADVCO**2.         
        ENDDO
C.......Writing wettedpressure for suction side.
        DO NN=1,NC/2
c     WRITE(20,*) SBP(NN),-CPBL(NC/2+NN,M)*ADVCO**2.
          WRITE(20,*) SBP(NN),-CPBLO(NC/2+NN,M)*ADVCO**2.         
        ENDDO
        

C**********************************************************************

C     
C     COMPUTE MASS DEFECT ALONG EACH STRIP
C     
        DO 310  I = 1, NC/2
          RV1(I) = SNGL(UEDGE(I)*DSTARS(I))*(CHDM*2.0)      
 310    CONTINUE
        DO 315  I = NC/2+1, NC+NWv-1      
          RV1(I) = SNGL(UEDGE(I)*DSTARP(I))*(CHDM*2.0)       
 315    CONTINUE

C     PUT BACK TO PROPCAV CONVENTION
        DO 320  I = 1, NC
          RRTEMP(I) = RV1(NC+2-I)
 320    CONTINUE
        DO 325  I = 1, NC
          RV1(I) = RRTEMP(I) 
 325    CONTINUE
        
c     if(M.EQ.5) then
c     do i = 1, nc+nwv+1
c     write(131,*) I, RV1(I)
c     enddo
c     endif 
C     
C     
C     COMPUTE BLOWING SOURCE STRGTH (SIGMA) ALONG EACH STRIP
C     
        DO 330  I = 1, NC
          L = INDEXB(I,M)
          BSRCB(L,1) = -(RV1(I+1)-RV1(I))/(SC(I+1,M)-SC(I,M))
          IF(I.GT.NC/2) THEN
            BSRCB(L,1) = - BSRCB(L,1)
          ENDIF

c     if(M.EQ.5) then     
c     write(131,*) I, BSRCB(L,1)   
c     endif 
 330    CONTINUE

        DO 340 I = 1, NTRA
          L = NTRA*(MR-M)+I
          BSRCWS(L,1)=(RV1(NC+I+1)-RV1(NC+I))/(SC(NC+I+1,M)-SC(NC+I,M))

c     if(M.EQ.5) then     
c     write(131,*) I, BSRCWS(L,1)   
c     endif
 340    CONTINUE

        DO 350 I = NSUB+1, NWMIN
          L = INDEXW(I,M) 
          BSRCW(L,1)=(RV1(NC+I+NTRA-NSUB+1)-RV1(NC+I+NTRA-NSUB)) 
     &      /(SC(NC+I+NTRA-NSUB+1,M)-SC(NC+I+NTRA-NSUB,M))

c     if(M.EQ.5) then     
c     write(131,*) I, BSRCW(L,1)   
c     endif
 350    CONTINUE 


 1000 CONTINUE                  !!     LOOP 1000 ENDED

      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
c     CLOSE(13)
c     CLOSE(14)
c     CLOSE(20)
      CLOSE(22)
c     CLOSE(97)

      RETURN
      END



      
