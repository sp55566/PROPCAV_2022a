C**************************** MIT PSF10 ********************************
C
C                PROPELLER STEADY FLOW ANALYSIS PROGRAM
C                        PANEL METHOD SOLUTION         
C
C                             Version 1.0
C
C       Copyright (c) Massachusetts Institute of Technology 1997
C
C                     Release Date: 1 January 1997
C
C***********************************************************************

C-----------------------------------------------------------------------
C     Next common blocks are for geometries used in wake alignment 
C-----------------------------------------------------------------------
C      COMMON /INT/  NX,NBLADE,ADVCO,RULT,RHULT,DCD,XULT,RHUB
C      COMMON /WAKE/ XR(15),XRSQ(15),XPI(15),XRAKE(15),XSKEW(15),
C     *              XCHD(15),XCI(15),XTI(15),XVA(15),XVR(15),XVT(15),
C     *              XUA(15),XUAU(15),XUT(15),XUTU(15)
      COMMON /dat1/  ADVCO,RULT,RHULT,DCD,XULT,RHUB,DTPROP
      COMMON /dat2/  NX,NBLADE
      COMMON /GINP/ XR(nxmax),XRSQ(nxmax),XPI(nxmax),XRAKE(nxmax)
     *   ,          XSKEW(nxmax),XCHD(nxmax),XCI(nxmax),XTI(nxmax)
     *   ,          XVA(nxmax),XVR(nxmax),XVT(nxmax)
     *   ,          XUA(nxmax),XUAU(nxmax),XUT(nxmax),XUTU(nxmax)



