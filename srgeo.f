       SUBROUTINE SRGEO(XIOLD,ETAOLD,XIMOLD,ETAMOLD)

       use m_BLADEM
       INCLUDE 'PUFCAV.INC'

!      COMMON/BLADEM/XBM(NBHZP,MBPZ),YBM(NBHZP,MBPZ),ZBM(NBHZP,MBPZ)
       DIMENSION XIM(NBHZP),ETAM(NBHZP)
       DIMENSION XIOLD(NBPZ,MBPZ),ETAOLD(NBPZ,MBPZ),XIMOLD(NBHZP,MBPZ),
     *      ETAMOLD(NBHZP,MBPZ),XXSR(NZSRP),YYSR(NZSRP),SBTMP(NBHZP),
     *      XMSR(NZSR2P),YMSR(NZSR2P)
       DIMENSION XXB(NBPZ),YYB(NBPZ),XXT(NBPZ),YYT(NBPZ),
     *      XXM(NBPZ),YYM(NBPZ),CUBTMPB(NBHZP*4-4),CUBTMPT(NBHZP*4-4),
     *      CUBTMPM(NBHZP*4-4)

C.....Read separated region information
       WRITE(*,*) 
       WRITE(*,*) ' PROPCAV> Separated Region (XSR, NSR2):'
       READ(*,*) XSR,NSR2

C.....Number of separated region panels.
       NSR2P=NSR2+1
       NSR=NSR2*2
       NSRP=NSR+1
       NTSR=NSR*MR

C.....Check parameters
       write(99,*) ' --------------------------------------'
       write(99,*) '          Separated Region Part        '
       write(99,*) ' --------------------------------------'
       IERSR=0
       IF(NSR2.GT.NZSR2) IERSR=1
       WRITE(99,110) NZSR2,NSR2,IERSR
 110   format(' Max NSR2   = ',i3,'    Input NSR2   = ',i3,
     %      '    Err = ',i1) 

       IF(IERSR.EQ.1) THEN
          write(*,*) ' ------------------------------------'
          write(*,*) '  You Have ERROR in INPUT DATA FILE! '
          write(*,*) '  Check <ERR.LOG> File to check error'
          write(*,*) ' ------------------------------------'
          stop
       END IF

C.....Original number of panels per each section on the blade.
       NCOLD=NC
       NHOLD=NCOLD/2

C.....Determine the end point of the SR on the wake surface.
       CALL GWAKE1

       DO M=1,MRP

C........Change the wake panel to local coordinates.
          IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
             DO N=1,NHP
                XIMOLD(N,M)=(XIOLD(NH+N,M)+XIOLD(NHP+1-N,M))/2.
                ETAMOLD(N,M)=(ETAOLD(NH+N,M)+ETAOLD(NHP+1-N,M))/2.
             END DO
             DO N=1,5
                XXB(N)=XW(N,M)-TWO*RAKE(M)
                YYB(N)=ZW(N,M)-SKEW(M)*TWO
             END DO
             XIW=XIMOLD(NHP,M)+XSR*CHORD(M)*TWO
             ETAW=YYB(5)
          ELSE
             TANP=PITCH(M)/PI/RZ(M)
             COSP=1.0/SQRT(1+TANP*TANP)
             SINP=TANP*COSP
             DO N=1,5
                TTMP=ATAN(ZW(N,M)/YW(N,M))
                ATMP=(TTMP-SKEW(M)*RAD)*RZ(M)
                BTMP=XW(N,M)-RAKE(M)*TWO
                XXB(N)=COSP*ATMP+SINP*BTMP
                YYB(N)=SINP*ATMP-COSP*BTMP             
             END DO

C...........Determine the end point of the SR in XI and ETA coord.
             CALL UGLYDK(5,1,1,XXB,YYB,0.,0.,CUBTMPB)
             XIW=XIMOLD(NHP,M)+XSR*CHORD(M)*TWO
             CALL EVALDKs(5,1,XXB,XIW,ETAW,CUBTMPB)
          END IF

C........Determine the camber line, and upper & lower surface
C........of the SR.
          DO N=1,NHP
             XXB(N)=XIOLD(NHP+1-N,M)
             YYB(N)=ETAOLD(NHP+1-N,M)
             XXT(N)=XIOLD(NH+N,M)
             YYT(N)=ETAOLD(NH+N,M)
             XXM(N)=XIMOLD(N,M)
             YYM(N)=ETAMOLD(N,M)
          END DO

          NT1=NHP+1
          XXB(NT1)=XIW
          YYB(NT1)=ETAW
          XXT(NT1)=XXB(NT1)
          YYT(NT1)=YYB(NT1)
          XXM(NT1)=XXB(NT1)
          YYM(NT1)=YYB(NT1)

          CALL UGLYDK(NT1,1,0,XXB,YYB,0.,0.,CUBTMPB)
          CALL UGLYDK(NT1,1,0,XXT,YYT,0.,0.,CUBTMPT)
          CALL UGLYDK(NT1,1,0,XXM,YYM,0.,0.,CUBTMPM)

          XXSR(1)=XIOLD(NCP,M)
          YYSR(1)=ETAOLD(NCP,M)
          XXSR(NSRP)=XIOLD(1,M)
          YYSR(NSRP)=ETAOLD(1,M)
          XXSR(NSR2P)=XXB(NT1)
          YYSR(NSR2P)=YYB(NT1)
          XMSR(1)=XIMOLD(NHP,M)
          YMSR(1)=ETAMOLD(NHP,M)
          XMSR(NSR2P)=XXM(NT1)
          YMSR(NSR2P)=YYM(NT1)

          DXSRT=(XIW-XIOLD(NCP,M))/FLOAT(NSR2)
          DXSRB=(XIW-XIOLD(1,M))/FLOAT(NSR2)
          DXSRM=(XIW-XIMOLD(NHP,M))/FLOAT(NSR2)
          DO N=2,NSR2
             NTOP=N
             NBOT=NSRP+1-N

C...........Bottom SR surface.
             XXX=XIOLD(1,M)+DXSRB*FLOAT(N-1)
             CALL EVALDKs(NT1,1,XXB,XXX,YYY,CUBTMPB)
             XXSR(NBOT)=XXX
             YYSR(NBOT)=YYY

C...........Top SR surface.
             XXX=XIOLD(NCP,M)+DXSRT*FLOAT(N-1)
             CALL EVALDKs(NT1,1,XXT,XXX,YYY,CUBTMPT)
             XXSR(NTOP)=XXX
             YYSR(NTOP)=YYY

             IF(YYSR(NTOP).LT.YYSR(NBOT)) THEN
                WRITE(*,*) 'Neg. thickness at initial SR region:',N,M
                WRITE(*,*) 'Try different XSR or NSR'
                WRITE(*,*) 'Stop Program'
                STOP
             END IF

C...........Mid SR surface.
             XXX=XIMOLD(NHP,M)+DXSRM*FLOAT(N-1)
             CALL EVALDKs(NT1,1,XXM,XXX,YYY,CUBTMPM)
             XMSR(NTOP)=XXX
             YMSR(NTOP)=YYY
          END DO

C........Re-calculate number of panels in chordwise direction
          NC=NCOLD+NSR
          NH=NC/2
          NCP=NC+1
          NHP=NH+1
       
          IF(NC.GT.NBZ) THEN
             WRITE(*,*) 'Currently, NBZ=',NBZ
             WRITE(*,*) 'Need to increase NBZ to ', NC
             WRITE(*,*) 'Please make changes and recompile program!'
             STOP
          END IF

C........Calculate new coordinates of the SR.      
          DO N=1,NSR2
             XI(N)=XXSR(NSR2P+N-1)
             ETA(N)=YYSR(NSR2P+N-1)
             XI(NCP+1-N)=XXSR(NSR2P+1-N)
             ETA(NCP+1-N)=YYSR(NSR2P+1-N)

             XIM(NHOLD+1+N)=XMSR(N+1)
             ETAM(NHOLD+1+N)=YMSR(N+1)
          END DO          

          DO N=1,NCOLD+1
             XI(NSR2+N)=XIOLD(N,M)
             ETA(NSR2+N)=ETAOLD(N,M)
             
             IF(N.LE.NHOLD+1) THEN
                XIM(N)=XIMOLD(N,M)
                ETAM(N)=ETAMOLD(N,M)
             END IF
          END DO

          IF(M.EQ.1) THEN
             DO N=1,NHOLD
                SBTMP(N)=SB(N)
             END DO
             DO N=1,NSR2
                DSB=XIM(NHOLD+N+1)-XIM(NHOLD+N)
                SBTMP(NHOLD+N)=SBTMP(NHOLD+N-1)+DSB
             END DO
          END IF

C........Plot blade sectsions (JY071201)
          IF(M.EQ.1) THEN
             OPEN(56,FILE='bldsec-sr.plt',STATUS='UNKNOWN')
             WRITE(56,*) 'VARIABLES="c/R","y/R"'
          END IF
          WRITE(56,*) 'ZONE T="M=',M,'"'
          DO N=1,NCP
             WRITE(56,*) XI(N)+CHORD(M),ETA(N)+RZ(M)
          END DO
          IF(M.EQ.MRP) CLOSE(56)

C........Generate new X,Y,Z coordinates of the blade with the SR.
          DO N=1,NCP
             IF(ICON.EQ.5.OR.ICON.EQ.6.OR.ICON.EQ.8) THEN
                XB(N,M)=RAKE(M) * TWO + XI(N)  
                YB(N,M)=RZ(M)
                ZB(N,M)=ETA(N) + SKEW(M) * TWO
             ELSE
                DX= XI(N)*SINP-ETA(N)*COSP
                XB(N,M)=RAKE(M)*TWO +DX
                THETA=SKEW(M)*RAD+(XI(N)*COSP+ETA(N)*SINP)/RZ(M)
                IF(M.EQ.1) THEN
                   THR(N)=THETA
                END IF
                IF(N.EQ.1) THEN
                   THT(M)=THETA
                END IF
                YB(N,M)=RZ(M)*COS(THETA)
                ZB(N,M)=RZ(M)*SIN(THETA)
                
                IF(N.LE.NHP) THEN
                   DX= XIM(N)*SINP-ETAM(N)*COSP
                   XBM(N,M)=RAKE(M)*TWO +DX
                   THETA=SKEW(M)*RAD+ (XIM(N)*COSP+ETAM(N)*SINP)/RZ(M)
                   YBM(N,M)=RZ(M)*COS(THETA)
                   ZBM(N,M)=RZ(M)*SIN(THETA)
                END IF
             END IF
          END DO

C........Restore the original index (JY070401)
          IF(M.LT.MRP) THEN
             NC=NCOLD
             NH=NC/2
             NCP=NC+1
             NHP=NH+1
          END IF

       END DO

C.....Define new chordwise spacings
       DO N=1,NH
          SB(N)=SBTMP(N)
       END DO

C.....Parameters relating to 
       N0(2)=NSR2+1
       N0(1)=NCP-NSR2

C.....Redefine total number of panels on the blade
       NPANB=NC*MR
    
C.....Calculate area and normal vector of trailing edge panel of
C.....original blade. (JY070901)
       N1=N0(2)
       N2=N0(1)

       IF(ISP .NE. 0) THEN
          DO M = 1 , MR+1
             CALL ROTATE2(0,SPANGLE,XB(N1,M),YB(N1,M))
             CALL ROTATE2(0,SPANGLE,XB(N2,M),YB(N2,M))
          ENDDO
       ENDIF
       
       DO M=1,MR
          XGW(M,1,1)=XB(N1,M)
          XGW(M,2,1)=XB(N2,M)
          XGW(M,3,1)=XB(N2,M+1)
          XGW(M,4,1)=XB(N1,M+1)
          XGW(M,1,2)=YB(N1,M)
          XGW(M,2,2)=YB(N2,M)
          XGW(M,3,2)=YB(N2,M+1)
          XGW(M,4,2)=YB(N1,M+1)
          XGW(M,1,3)=ZB(N1,M)
          XGW(M,2,3)=ZB(N2,M)
          XGW(M,3,3)=ZB(N2,M+1)
          XGW(M,4,3)=ZB(N1,M+1)
       END DO

       CALL GEO3DW(MR,XGW,CHRLEWS,IER)
       IF(IER.EQ.0) THEN
          WRITE(*,*) 'UNACCEPTABLE PANEL IN XXTE'
          STOP
       END IF

       DO M=1,MR
          SSTE(M)=SSW(M,1)
          DO KK=1,6
             VELTE(M,KK)=VELW(M,KK)
          END DO
       END DO

       IF(ISP .NE. 0) THEN
          DO M = 1 , MR+1
             CALL ROTATE2(-1,SPANGLE,XB(N1,M),YB(N1,M))
             CALL ROTATE2(-1,SPANGLE,XB(N2,M),YB(N2,M))
          ENDDO
       ENDIF


      
       RETURN
       END



