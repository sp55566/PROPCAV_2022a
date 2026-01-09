C ============================================
C
C  This routine generates hub panels for the Thick blade section 
C  and/or low pitch propeller.
C
C  Programmed by Hanseong Lee
C
C ------------------------------------
      SUBROUTINE GHUB7
C ------------------------------------
      use m_NACAIO
      INCLUDE 'PUFCAV.INC'
!     COMMON /NACAIO/ XKEY(NBPZ),YKEY(NBPZ),ZKEY(NBPZ)
!    %     ,XNEAR(NBPZ),YNEAR(NBPZ),CROSS(3,NBPZ)
      DIMENSION XXX(NBPZ,MHPZ+1),ZZZ(NBPZ,MHPZ+1)
      DIMENSION TH(NHPZ,MHPZ)
      DIMENSION XFC(NBHZP),YFC(NBHZP),ZFC(NBHZP),YTFC(NBHZP)
      DIMENSION XBK(NBHZP),YBK(NBHZP),ZBK(NBHZP),YTBK(NBHZP)
      
C-----------------------------------------------------------------------
C     PREPARE DATA FOR GENERATING HUB PANELS IN INTER-BLADE REGION
C     RIGHT-HANDED PROPELLER 
C-----------------------------------------------------------------------

C*****BLADE : FACE CURVE  

      DO N=1,NHP
         XFC(N)  = XB(NHP+1-N,1)
         YFC(N)  = YB(NHP+1-N,1)
         ZFC(N)  = ZB(NHP+1-N,1)
         YTFC(N) = DANGLE(ZFC(N),YFC(N))+DELK
      ENDDO

C*****BLADE : BACK CURVE

      DO N=1,NHP
         XBK(N)  = XB(NHP-1+N,1)
         YBK(N)  = YB(NHP-1+N,1)
         ZBK(N)  = ZB(NHP-1+N,1)
         YTBK(N) = DANGLE(ZBK(N),YBK(N))
      ENDDO

      NHBDT=NSW(1)
      DO N=1,NSW(1)-1
         IF(XW(N,1) .LT. XHBFD .AND. XW(N+1,1) .GT. XHBFD)  NHBDT=N
      ENDDO

      NHBU1=NHBU+1
      NHBUB=NHBU+NH
      NHBUB1=NHBUB+1
      NHBX=NHBUB+NHBDT
      NHBX1 = NHBX+1
      NHBDT1 = NHBDT+1
      MHBT1=MHBT+1

C
C-- From Blade LE to TE
C       Transverse grid line and circumferential grid points are
C       generated
C
        DO N=1,NCP
           CROSS(1,N) = XB(N,1)
           CROSS(2,N) = YB(N,1)
           CROSS(3,N) = ZB(N,1)
        ENDDO

        CALL INTERBLADE(XXX,ZZZ,NBLADE,NH,MHBT,HPCH1,HPCH2)

        ICC = 0
        DO N = NHBU1, NHBUB1
           ICC = ICC + 1
           DO M = 1, MHBT+1
              XH(N,M) = XXX(ICC,M)
              ZH(N,M) = ZZZ(ICC,M)
           ENDDO
        ENDDO

C
        DO M = 1, MHBT+1
           DO N = NHBU1, NHBUB1
              YTEMP = RHUB
              ZTEMP = ZH(N,M)
              YH(N,M) = YTEMP * COS( ZTEMP )
              ZH(N,M) = YTEMP * SIN( ZTEMP )
           ENDDO
        ENDDO

        DO N = NHBU1, NHBUB1
           DO M = 1, MHBT1
              TH(N,M) = DANGLE(ZH(N,M),YH(N,m))
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
                  
      DPSI = 0.5*PI/REAL(NHBU)

      DO M=1,MHBT1
         DL = -XHBU-XH(NHBU+1,M)
         DO N=1,NHBU
            XH(NHBU+1-N,M)=XH(NHBU+1,M)+DL*(1.0-COS(DPSI*REAL(N)))
            XTMP = XH(NHBU+1-N,M)
            DX = XH(NHBU+1-N,M)-XH(NHBU+2-N,M)
            DTH = DX/HC
            TH(NHBU+1-N,M) = TH(NHBU+2-N,M)+DTH
            YH(NHBU+1-N,M) = RHUB*COS(TH(NHBU+1-N,M))
            ZH(NHBU+1-N,M) = RHUB*SIN(TH(NHBU+1-N,M))
         ENDDO
      ENDDO

C-----------------------------------------------------------------------
C     Hub geometry at downstream of the blade
C-----------------------------------------------------------------------

C -- When M=1, and M=MHBT1

        DO N = 1, NHBDT1
 
           NN = NHBUB + N

           IF(N .EQ. NHBDT1) THEN
              XH(NN,1) = XHBFD
              DDX1 = XH(NN,1) - XH(NN-1,1)
              DDX2 = XW(N,1) - XW(N-1,1)
              TH2 = DANGLE(ZW(N,1),YW(N,1))
              TH1 = DANGLE(ZW(N-1,1),YW(N-1,1))
              TH(NN,1) = TH(NN-1,1) + (TH2-TH1)*DDX1/DDX2
           ELSE
              XH(NN,1)=XW(N,1)
              TH(NN,1) = DANGLE(ZW(N,1),YW(N,1))
           ENDIF
           
           AAXX =  (XHBFD-XW(N,1))/XHBT
           RNOS=RNOSE( AAXX,RHUB )

           IF(RNOS.LT.RHULT) THEN
              RR=RNOS
              YH(NN,1)=RR*COS(TH(NN,1))
              ZH(NN,1)=RR*SIN(TH(NN,1))
           ELSE
              RR=SQRT(YW(N,1)**2+ZW(N,1)**2)
              YH(NN,1)=YW(N,1)
              ZH(NN,1)=ZW(N,1)
           END IF

           XH(NN,MHBT1)=XW(N,1)
           TH(NN,MHBT1)=TH(NN,1) + DELK
           YH(NN,MHBT1)=RR*COS(TH(NN,MHBT1))
           ZH(NN,MHBT1)=RR*SIN(TH(NN,MHBT1))
        ENDDO

        DELDK = DELK / REAL(MHBT)
        DL1 = XH(NHBX1,1) - XH(NHBUB1,1)

        DIFF = TH(NHBX1,1) - TH(NHBUB1,1)

        IF(DIFF .LE. 0.0) DELDK = -DELDK

        DO M = 2, MHBT
           DL2 = XH(NHBX1,1) - XH(NHBUB1,M)
           RATIO = DL2/DL1
           
           DDT1 = TH(NHBX1,1)+DELDK*(M-1)
           DDT = DDT1 -TH(NHBUB1,M)

           DTH = DDT/DL2
    
           DO N = 1, NHBDT
              NN = NHBUB1+N
              DDX =  XH(NN,1) - XH(NHBUB1,1)
              XH(NN,M) = XH(NHBUB1,M) + DDX*RATIO

              RNOS=RNOSE( (XHBFD-XH(NN,M))/XHBT,RHUB )

              TH(NN,M) = TH(NHBUB1,M) + DDX*RATIO*DTH

              RR=RNOS
              YH(NN,M)=RR*COS(TH(NN,M))
              ZH(NN,M)=RR*SIN(TH(NN,M))
           ENDDO
        ENDDO
        
      DO M = 1, MHBT1
         XH(NHBX1,M)=XH(NHBX1,1)
         YH(NHBX1,M)=ZERO
         ZH(NHBX1,M)=ZERO
      ENDDO

      NHBDK=MHBT
      NHBV=0
      DTH=DELK/NHBDK

      DO N=1,NHBDK+1
         XHDK(N,1)=XH(1,N)
         YHDK(N,1)=ZERO
         ZHDK(N,1)=ZERO
         XHDK(N,1)=XH(1,N)
         YHDK(N,1)=YH(1,N)
         ZHDK(N,1)=ZH(1,N)
      ENDDO

      RETURN
      END
