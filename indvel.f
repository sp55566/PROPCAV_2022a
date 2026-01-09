************************************************************************
      SUBROUTINE INDVEL
************************************************************************
*     INDVEL: Induced velocity on the wake surface                     *
************************************************************************

        INCLUDE 'PUFCAV.INC'
        INCLUDE 'PUFCAVB.INC'

C      DIMENSION VX(NPAWZ),VY(NPAWZ),VZ(NPAWZ)
      DIMENSION XPP(3)

CVV
      ALLOCATABLE :: VX(:),VY(:),VZ(:)
CVV

      CALL INFLOWK2

C --- Calculate Image of Control points to evaluate induced velocities
C      due to all blades, However, the induced velocities
C       at the key wake control points are directly evaluated using the
C       other blade wakes panels.
C       The effects of tip vortex and Tip bulb of other blades 
C       are ignored. 


CVV
      ALLOCATE(VX(NPAWZ),VY(NPAWZ),VZ(NPAWZ))
CVV

      CALL CLEAR(VX,NPAWZ)
      CALL CLEAR(VY,NPAWZ)
      CALL CLEAR(VZ,NPAWZ)

C -- Cal. induced velocity at the wake C.P. due to the blade panels

      IMR1 = 1
      
      DO I = 1 , NPWAKE
         
         DO J = 1 , NPANB

            DO K = 1 , 4
             XV(K) = XVP(J,K)
             YV(K) = YVP(J,K)
             SIDE(K) = SID(J,K)
            ENDDO
            
            DO K = 1 , 15
             S(K) = SS(J,K)
            ENDDO

            DO KK = 1 , NBLADE
             XPP(1) = XCTWO(I,1,KK) + DELTAM * DIRWO(I,3,1,KK) 
             XPP(2) = XCTWO(I,2,KK) + DELTAM * DIRWO(I,3,2,KK) 
             XPP(3) = XCTWO(I,3,KK) + DELTAM * DIRWO(I,3,3,KK)

             XLOC = 0.0
             YLOC = 0.0
             ZLOC = 0.0
                
             DO K = 1 , 3
                XLOC = XLOC+(XPP(K)-XCT(J,K))*DIR(J,1,K)
                YLOC = YLOC+(XPP(K)-XCT(J,K))*DIR(J,2,K)
                ZLOC = ZLOC+(XPP(K)-XCT(J,K))*DIR(J,3,K)
             ENDDO

             IMR = IMR1
             CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),
     %                     FS,FD,FSX,FSY,FDX,FDY,FDZ,1,IMR)
             
             VDX=FDX*DIR(J,1,1)+FDY*DIR(J,2,1)+FDZ*DIR(J,3,1)
             VDY=FDX*DIR(J,1,2)+FDY*DIR(J,2,2)+FDZ*DIR(J,3,2)
             VDZ=FDX*DIR(J,1,3)+FDY*DIR(J,2,3)+FDZ*DIR(J,3,3)
             VSX=FSX*DIR(J,1,1)+FSY*DIR(J,2,1)-FD*DIR(J,3,1)
             VSY=FSX*DIR(J,1,2)+FSY*DIR(J,2,2)-FD*DIR(J,3,2)
             VSZ=FSX*DIR(J,1,3)+FSY*DIR(J,2,3)-FD*DIR(J,3,3)

             THE=-DELK*FLOAT(KK-1)

             VDY1=VDY*COS(THE)-VDZ*SIN(THE)
             VDZ1=VDY*SIN(THE)+VDZ*COS(THE)
             VSY1=VSY*COS(THE)-VSZ*SIN(THE)
             VSZ1=VSY*SIN(THE)+VSZ*COS(THE)
             VDY=VDY1
             VDZ=VDZ1
             VSY=VSY1
             VSZ=VSZ1
             VX(I)=VX(I)-(POT(J)*VDX-DPDN(J)*VSX)/4.0/PI
             VY(I)=VY(I)-(POT(J)*VDY-DPDN(J)*VSY)/4.0/PI
             VZ(I)=VZ(I)-(POT(J)*VDZ-DPDN(J)*VSZ)/4.0/PI
            ENDDO
         ENDDO
      ENDDO

      IF(IHUB .NE. 0) THEN

         DO I = 1 , NPWAKE
            
            DO J = NPANB+1 , NPANB+NPANH

             DO K = 1 , 4
                XV(K) = XVP(J,K)
                YV(K) = YVP(J,K)
                SIDE(K) = SID(J,K)
             ENDDO
             
             DO K = 1 , 15
                S(K) = SS(J,K)
             ENDDO

             DO KK = 1 , NBLADE
                XPP(1) = XCTWO(I,1,KK) + DELTAM * DIRWO(I,3,1,KK) 
                XPP(2) = XCTWO(I,2,KK) + DELTAM * DIRWO(I,3,2,KK) 
                XPP(3) = XCTWO(I,3,KK) + DELTAM * DIRWO(I,3,3,KK)

                XLOC = 0.0
                YLOC = 0.0
                ZLOC = 0.0
                
                DO K = 1 , 3
                   XLOC = XLOC+(XPP(K)-XCT(J,K))*DIR(J,1,K)
                   YLOC = YLOC+(XPP(K)-XCT(J,K))*DIR(J,2,K)
                   ZLOC = ZLOC+(XPP(K)-XCT(J,K))*DIR(J,3,K)
                ENDDO
                
                IMR = IMR1
                CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),
     %                   FS,FD,FSX,FSY,FDX,FDY,FDZ,1,IMR)
                
                VDX=FDX*DIR(J,1,1)+FDY*DIR(J,2,1)+FDZ*DIR(J,3,1)
                VDY=FDX*DIR(J,1,2)+FDY*DIR(J,2,2)+FDZ*DIR(J,3,2)
                VDZ=FDX*DIR(J,1,3)+FDY*DIR(J,2,3)+FDZ*DIR(J,3,3)
                VSX=FSX*DIR(J,1,1)+FSY*DIR(J,2,1)-FD*DIR(J,3,1)
                VSY=FSX*DIR(J,1,2)+FSY*DIR(J,2,2)-FD*DIR(J,3,2)
                VSZ=FSX*DIR(J,1,3)+FSY*DIR(J,2,3)-FD*DIR(J,3,3)
                
                THE=-DELK*FLOAT(KK-1)
                
                VDY1=VDY*COS(THE)-VDZ*SIN(THE)
                VDZ1=VDY*SIN(THE)+VDZ*COS(THE)
                VSY1=VSY*COS(THE)-VSZ*SIN(THE)
                VSZ1=VSY*SIN(THE)+VSZ*COS(THE)
                VDY=VDY1
                VDZ=VDZ1
                VSY=VSY1
                VSZ=VSZ1
                VX(I)=VX(I)-(POT(J)*VDX-DPDN(J)*VSX)/4.0/PI
                VY(I)=VY(I)-(POT(J)*VDY-DPDN(J)*VSY)/4.0/PI
                VZ(I)=VZ(I)-(POT(J)*VDZ-DPDN(J)*VSZ)/4.0/PI
             ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF(IAN .EQ. 2) THEN        

         DO I = 1 , NPWAKE   

            DO J = NPANB+NPANH+1 , NPANEL

             DO K = 1 , 4
                XV(K) = XVP(J,K)
                YV(K) = YVP(J,K)
                SIDE(K) = SID(J,K)
             ENDDO
             
             DO K = 1 , 15
                S(K) = SS(J,K)
             ENDDO

             XPP(1) = XCTWO(I,1,1) + DELTAM * DIRWO(I,3,1,1) 
             XPP(2) = XCTWO(I,2,1) + DELTAM * DIRWO(I,3,2,1) 
             XPP(3) = XCTWO(I,3,1) + DELTAM * DIRWO(I,3,3,1)

             XLOC = 0.0
             YLOC = 0.0
             ZLOC = 0.0
                
             DO K = 1 , 3
                XLOC = XLOC+(XPP(K)-XCT(J,K))*DIR(J,1,K)
                YLOC = YLOC+(XPP(K)-XCT(J,K))*DIR(J,2,K)
                ZLOC = ZLOC+(XPP(K)-XCT(J,K))*DIR(J,3,K)
             ENDDO
                
             IMR = IMR1
             CALL RPAN(XLOC,YLOC,ZLOC,CHRLEPS(J),
     %                   FS,FD,FSX,FSY,FDX,FDY,FDZ,1,IMR)

             VDX=FDX*DIR(J,1,1)+FDY*DIR(J,2,1)+FDZ*DIR(J,3,1)
             VDY=FDX*DIR(J,1,2)+FDY*DIR(J,2,2)+FDZ*DIR(J,3,2)
             VDZ=FDX*DIR(J,1,3)+FDY*DIR(J,2,3)+FDZ*DIR(J,3,3)
             VSX=FSX*DIR(J,1,1)+FSY*DIR(J,2,1)-FD*DIR(J,3,1)
             VSY=FSX*DIR(J,1,2)+FSY*DIR(J,2,2)-FD*DIR(J,3,2)
             VSZ=FSX*DIR(J,1,3)+FSY*DIR(J,2,3)-FD*DIR(J,3,3)
             VX(I)=VX(I)-(POT(J)*VDX-DPDN(J)*VSX)/4.0/PI
             VY(I)=VY(I)-(POT(J)*VDY-DPDN(J)*VSY)/4.0/PI
             VZ(I)=VZ(I)-(POT(J)*VDZ-DPDN(J)*VSZ)/4.0/PI
            ENDDO
         ENDDO
      ENDIF
            
C -- Cal. induced velcoity at the wake C.P. due to the wake panels

      NSAPW2 = NWPANEL *  MR

      DO I = 1 , NPWAKE

         XPP(1) = XCTWO(I,1,1) + DELTAM * DIRWO(I,3,1,1) 
         XPP(2) = XCTWO(I,2,1) + DELTAM * DIRWO(I,3,2,1) 
         XPP(3) = XCTWO(I,3,1) + DELTAM * DIRWO(I,3,3,1) 
         
         DO KK = 1 , NBLADE

            IREC = NTPOS(KK)
            CALL READ2(46,IREC,TEMP5,NSAPW2)
            
            DO J = 1 , NPWAKE
             DO K = 1 , 4
                XV(K) = XVPWO(J,K,KK)
                YV(K) = YVPWO(J,K,KK)
                SIDE(K) = SIDWO(J,K,KK)
             ENDDO
         
             DO K = 1 , 15
                S(K) = SSWO(J,K,KK)
             ENDDO
             
             XLOC = 0.0
             YLOC = 0.0
             ZLOC = 0.0
                
             DO K = 1 , 3
                XLOC = XLOC + (XPP(K)-XCTWO(J,K,KK)) 
     %                     * DIRWO(J,1,K,KK)
                YLOC = YLOC + (XPP(K)-XCTWO(J,K,KK)) 
     %                     * DIRWO(J,2,K,KK)
                ZLOC = ZLOC + (XPP(K)-XCTWO(J,K,KK))
     %                        * DIRWO(J,3,K,KK)
                ENDDO

                CRTL = CHRLEWSO(J,KK)

              IMR = IMR1
              CALL RPAN(XLOC,YLOC,ZLOC,CRTL,FS,FD,FSX,FSY,
     %                     FDX,FDY,FDZ,1,IMR)

C              IF(I .EQ. J .and. KK .eq. 1) THEN
C                 FDX = 0.0
C                 FDY = 0.0
C                 FDZ = 0.0
C              ENDIF
                  VDX=FDX*DIRWO(J,1,1,KK)+FDY*DIRWO(J,2,1,KK)
     %                        +FDZ*DIRWO(J,3,1,KK)
                  VDY=FDX*DIRWO(J,1,2,KK)+FDY*DIRWO(J,2,2,KK)
     %                        +FDZ*DIRWO(J,3,2,KK)
                  VDZ=FDX*DIRWO(J,1,3,KK)+FDY*DIRWO(J,2,3,KK)
     %                        +FDZ*DIRWO(J,3,3,KK)

              VX(I)=VX(I)-(TEMP5(J)*VDX)/4.0/PI
              VY(I)=VY(I)-(TEMP5(J)*VDY)/4.0/PI
              VZ(I)=VZ(I)-(TEMP5(J)*VDZ)/4.0/PI
             ENDDO
          ENDDO
       ENDDO

C --- Total velocity on the wake surface

       DO I = 1 , NPWAKE
          XIND(I) = VX(I) + VOXW(I) 
          YIND(I) = VY(I) + VOYW(I)
          ZIND(I) = VZ(I) + VOZW(I)
       ENDDO

       DEALLOCATE(VX,VY,VZ)
       RETURN
       END



C ===================================
C     Calculate Wake Pitch
C
C     P/D (r) = Pi*r*Va / Vt
C
      Subroutine Wakepitch
C==================================
      include 'PUFCAV.INC'
      include 'PUFCAVB.INC'
      DIMENSION PPCUB(1000),RRCUB(1000)
      dimension rrr(150),ppdd(150),xxxx(150)
      dimension pv1(30),pv2(30),rv1(30),rv2(30)

      PI = ACOS(-1.)
      
      DO J = 1, MR
         DO I = 1 , NWPANEL
            I1 = INDEXW2(I,J)
            RRR(I) = SQRT (XCTWo(I1,2,1)**2+XCTWo(I1,3,1)**2)
            COSP = XCTWo(I1,2,1)/RRR(i)
            SINP = XCTWo(I1,3,1)/RRR(i)
            UAAA = XIND(I1)
            UTTT = -SINP*YIND(I1) + COSP*ZIND(I1)
            PPDD(I) = PI * RRR(i) * UAAA / UTTT
            XXXX(I) = Xctwo(I1,1,1)
         ENDDO

       CALL UGLYDK(NWPANEL,1,1,XXXX,PPDD,0.0,0.0,PPCUB)
       CALL UGLYDK(NWPANEL,1,1,XXXX,RRR,0.0,0.0,RRCUB)

         CALL EVALDKs(NWPANEL,1,XXXX,0.328,PV1(J),PPCUB)
         CALL EVALDKs(NWPANEL,1,XXXX,0.328,RV1(J),RRCUB)

         CALL EVALDKs(NWPANEL,1,XXXX,0.95,PV2(J),PPCUB)
         CALL EVALDKs(NWPANEL,1,XXXX,0.95,RV2(J),RRCUB)       
      ENDDO

C      write(777,*) 'zone T="x=0.328"'
C      DO I = 1 , MR
C         WRITE(777,*) RV1(I),PV1(I)
C      ENDDO

C      write(777,*) 'zone T="x=0.95"'
C      DO I = 1 , MR
C         WRITE(777,*) RV2(I),PV2(I)
C      ENDDO

      RETURN
      END
