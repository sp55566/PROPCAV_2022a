C***************************** MIT PSF10 ********************************
C
C                PROPELLER STEADY FLOW ANALYSIS PROGRAM
C                        PANEL METHOD SOLUTION         
C
C
C       Copyright (c) Massachusetts Institute of Technology 1997
C
C                     Release Date: 1 January 1997
C
C***********************************************************************
      SUBROUTINE GHUB6
C***********************************************************************
C     GHUB: Geometry of the POD
C     --- Generate a pod geometry (wake aligned)
C     --- Scheme A (INTER-BLADE REGION NON-UNIFORM MESH)

      use m_NACAIO
      INCLUDE 'PUFCAV.INC'
!     COMMON /NACAIO/ XKEY(NBPZ),YKEY(NBPZ),ZKEY(NBPZ)
!    #     ,XNEAR(NBPZ),YNEAR(NBPZ),CROSS(3,NBPZ)
      DIMENSION XFC(100),YFC(100),ZFC(100),YTFC(100)
      DIMENSION XBK(100),YBK(100),ZBK(100),YTBK(100)
      DIMENSION TH(NHPZ,MHPZ)
      DIMENSION XXX(NBPZ,MHPZ+1),ZZZ(NBPZ,MHPZ+1)

C-----------------------------------------------------------------------

      MHBT1 = MHBT+1
      NH1   = NH/2
      PI2  = 0.5*PI

C-----------------------------------------------------------------------
C     PREPARE DATA FOR GENERATING HUB PANELS IN INTER-BLADE REGION
C     RIGHT-HANDED PROPELLER 
C-----------------------------------------------------------------------

      
      NNN1 = 1
      NNN2 = NHBU + 1
      NNN3 = NHBU + NHP
      NNN4 = NNN3+NHBDT
      MHBT1 = MHBT + 1
      NHBB=NH

C*****BLADE : FACE CURVE      
      DO N=1,NHP
         XFC(N)  = XB(NHP+1-N,1)
         YFC(N)  = YB(NHP+1-N,1)
         ZFC(N)  = ZB(NHP+1-N,1)
         YTFC(N) = ATAN2(ZFC(N),YFC(N))+DELK
      ENDDO

C*****BLADE : BACK CURVE

      DO N=1,NHP
         XBK(N)  = XB(NHP-1+N,1)
         YBK(N)  = YB(NHP-1+N,1)
         ZBK(N)  = ZB(NHP-1+N,1)
         YTBK(N) = ATAN2(ZBK(N),YBK(N))
      ENDDO
      
      DO N=1,NCP
         CROSS(1,N) = XB(N,1)
         CROSS(2,N) = YB(N,1)
         CROSS(3,N) = ZB(N,1)
      ENDDO

      CALL INTERBLADE(XXX,ZZZ,NBLADE,NH,MHBT,HPCH1,HPCH2)
     
      ICC = 0
      DO N=NNN2,NNN3
         ICC = ICC + 1
         DO M=1,MHBT+1
            XH(N,M) = XXX(ICC,M)
            ZH(N,M) = ZZZ(ICC,M)
         ENDDO
      ENDDO

C
C-- Calcualte Y,Z component from NNN2 to NNN3
C
      DO M=1,MHBT+1
         DO N=NNN2,NNN3
            CALL EVALDKs(NHIN,1,XHIN,XH(N,M),YTEMP,XYHINCUB)            
            ZTEMP = ZH(N,M)
            YH(N,M) = YTEMP * COS( ZTEMP )
            ZH(N,M) = YTEMP * SIN( ZTEMP )
         ENDDO
      ENDDO

      DO N = NNN2,NNN3
         DO M = 1, MHBT1
            TH(N,M) = ATAN2(ZH(N,M),YH(N,M))
         ENDDO
      ENDDO

C-----------------------------------------------------------------------            
C     UPSTREAM HUB      
C-----------------------------------------------------------------------
                  
      X2  = XFC(1)
      X1  = (XFC(2)+XBK(2))/2.0
      TH2 = YTFC(1)
      TH1 = (YTFC(2)+YTBK(2)+DELK)/2.0
      HC  = (X2-X1)/(TH2-TH1)
                  
      DPSI = PI2/(NHBU)

      DO M=1,MHBT1
         DL = XHIN(1)-XH(NHBU+1,M)
         DO N=1,NHBU
            XH(NHBU+1-N,M)=XH(NHBU+1,M)+DL*(1.0-COS(DPSI*REAL(N)))
            XTMP = XH(NHBU+1-N,M)
            DX = XH(NHBU+1-N,M)-XH(NHBU+2-N,M)
            DTH = DX/HC
            CALL EVALDKs(NHIN,1,XHIN,XTMP,RTMP,XYHINCUB)
            TH(NHBU+1-N,M) = TH(NHBU+2-N,M)+DTH
            YH(NHBU+1-N,M) = RTMP*COS(TH(NHBU+1-N,M))
            ZH(NHBU+1-N,M) = RTMP*SIN(TH(NHBU+1-N,M))
         ENDDO
      ENDDO
             
C-----------------------------------------------------------------------            
C     WAKE-ALIGNED HUB REGION       
C-----------------------------------------------------------------------
 1111 CONTINUE

      DO NN = 1,NSW(1)
         IF(XW(NN,1).GE.XHIN(NHIN)) EXIT
      ENDDO

      NHBDT = NN-1
      NHBX  = NHBU+NH+NHBDT
      NN1   = NHBU+NH      
      NHBUB = NN1
      NHBUB1 = NHBUB + 1
      NHBX1 = NHBX + 1
      MHBT1 = MHBT + 1
 
      DO M=1,MHBT1
         XH(NN1+NHBDT+1,M) = XHIN(NHIN)
         YH(NN1+NHBDT+1,M) = 0.0
         ZH(NN1+NHBDT+1,M) = 0.0
      ENDDO

      DO N=2,NHBDT
         TH1 = ATAN2(REAL(ZW(N,1)),REAL(YW(N,1)))
         IF(TH1 .LT. TH(NN1,1)) TH1 = TH1 + TWOPI
         DO M=1,MHBT1,MHBT
            XH(NN1+N,M) = XW(N,1)
            TH(NN1+N,M) = TH1+DELK*REAL(M-1)/REAL(MHBT)
            XTMP = XW(N,1)
            CALL EVALDKs(NHIN,1,XHIN,XTMP,RTMP,XYHINCUB)
            YH(NN1+N,M) = RTMP*COS(TH(NN1+N,M))
            ZH(NN1+N,M) = RTMP*SIN(TH(NN1+N,M))
         ENDDO
      ENDDO

      DELDK = DELK / REAL(MHBT)
      DL1 = XH(NHBX1,1) - XH(NHBUB1,1)      
      
      DO M = 2, MHBT
         DL2 = XH(NHBX1,1) - XH(NHBUB1,M)
         RATIO = DL2/DL1


         DO N = 1, NHBDT-1

            NN = NHBUB1+N
            DDX =  XH(NN,1) - XH(NHBUB1,1)

            DTH1 = TH(NN-1,1)
            DTH2 = TH(NN,1)
            DTH = DTH2 - DTH1

            XH(NN,M) = XH(NHBUB1,M) + DDX*RATIO
            XTMP = XH(NN,M)

            CALL EVALDKs(NHIN,1,XHIN,XTMP,RTMP,XYHINCUB)

            TH(NN,M) = TH(NN-1,M) + RATIO*DTH
            YH(NN,M)=RTMP*COS(TH(NN,M))
            ZH(NN,M)=RTMP*SIN(TH(NN,M))
         ENDDO

         XH(NN+1,M) = XH(NN+1,1)
         YH(NN+1,M) = YH(NN+1,1)
         ZH(NN+1,M) = ZH(NN+1,1)

      ENDDO

      NHBDK=MHBT
      DTH=DELK/NHBDK
      DO N=1,NHBDK+1
         XHDK(N,1)=XH(1,N)
         YHDK(N,1)=ZERO
         ZHDK(N,1)=ZERO
         XHDK(N,2)=XH(1,N)
         YHDK(N,2)=YHIN(1)*COS(DTH*(N-1))
         ZHDK(N,2)=YHIN(1)*SIN(DTH*(N-1))
      ENDDO

      RETURN
      END
