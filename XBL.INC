C
C
      PARAMETER (NCOM=55)
      REAL*8 MDI
      REAL*8 M1,M1_U1,M2,M2_U2
      REAL*8 MA 
      REAL*8 COM1(NCOM), COM2(NCOM)

      LOGICAL SIMI,TRAN,TURB,WAKE
      LOGICAL TRFORC,TRFREE
c---- wbrewer -- 4/10/95 --- add dtcav to common statement ---
c-----dianr --
      COMMON/CAV/ITYPCV,ICV,JPAN,IJUMP,ICLE,ICTE
      COMMON/CAV1/DTCAV,CAVL,CFTOTX,CFTOTY
c     COMMON/CAV2/NT
       
c---------------------------------
C
      COMMON/VAR1/ X1,  U1,  T1,  D1,  S1, AMPL1, U1_UEI, DW1
     &        , H1, H1_T1, H1_D1
     &        , M1, M1_U1
     &        , R1, R1_U1
     &        , V1, V1_U1, V1_RE
     &        , HK1, HK1_U1, HK1_T1, HK1_D1
     &        , HS1, HS1_U1, HS1_T1, HS1_D1,         HS1_RE
     &        , HC1, HC1_U1, HC1_T1, HC1_D1
     &        , RT1, RT1_U1, RT1_T1,                 RT1_RE
     &        , CF1, CF1_U1, CF1_T1, CF1_D1,         CF1_RE
     &        , DI1, DI1_U1, DI1_T1, DI1_D1, DI1_S1, DI1_RE
     &        , CQ1, CQ1_U1, CQ1_T1, CQ1_D1,         CQ1_RE
     &        , DE1, DE1_U1, DE1_T1, DE1_D1
      COMMON/VAR2/ X2,  U2,  T2,  D2,  S2, AMPL2, U2_UEI, DW2
     &        , H2, H2_T2, H2_D2
     &        , M2, M2_U2
     &        , R2, R2_U2
     &        , V2, V2_U2, V2_RE
     &        , HK2, HK2_U2, HK2_T2, HK2_D2
     &        , HS2, HS2_U2, HS2_T2, HS2_D2,         HS2_RE
     &        , HC2, HC2_U2, HC2_T2, HC2_D2
     &        , RT2, RT2_U2, RT2_T2,                 RT2_RE
     &        , CF2, CF2_U2, CF2_T2, CF2_D2,         CF2_RE
     &        , DI2, DI2_U2, DI2_T2, DI2_D2, DI2_S2, DI2_RE
     &        , CQ2, CQ2_U2, CQ2_T2, CQ2_D2,         CQ2_RE
     &        , DE2, DE2_U2, DE2_T2, DE2_D2
      EQUIVALENCE (X1,COM1(1)), (X2,COM2(1))
C
      COMMON/VARA/ CFA, CFA_RE
     &        , CFA_U1, CFA_T1, CFA_D1
     &        , CFA_U2, CFA_T2, CFA_D2
C
      COMMON/SAV/ C1SAV(NCOM), C2SAV(NCOM)
C
      COMMON/VAR/ XSI, AMI, CTI, UEI, THI, DSI, DSWAKI, DWTE
     &          , QINFBL, TKBL
     &          , RSTBL, HSTINV, REYBL, GAMBL, GM1BL, HVRAT
     &          , BULE, XIFORC, AMCRIT
     &          , SIMI,TRAN,TURB,WAKE
     &          , TRFORC,TRFREE,MDI
C
      COMMON/XTRAN/ XT, XT_A1
     &   , XT_X1, XT_T1, XT_D1, XT_U1
     &   , XT_X2, XT_T2, XT_D2, XT_U2
C
      COMMON/SYS/ VS1(4,5),VS2(4,5),VSREZ(4),VSR(4)
C

