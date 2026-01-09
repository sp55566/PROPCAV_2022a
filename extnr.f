      SUBROUTINE EXTNR(M,IDR)
************************************************************************
*                                                                      *
*   Subroutine NREX uses a two-point linear extrapolation of delta(l)  *
*   to find a new cavity length l(r). The extrapolation is the same    *
*   one that Newton-Raphson uses.                                      *
*                                                                      *
*   07-15-91  Neal Fine                                                *
*                                                                      *
*   Date of last revision                       Revision               *
*   ---------------------                    --------------            *
*    07-24-91  NF                -removed the dependence on panelling  *
*                                 because of the split-panel method    *
*    CM081297       The algorithm was modified so that the length      *
*                   change is not permitted at r/R > 0.95              *
*    CM011598       Added icon.ne.5 to first if statement in the sub-  *
*                   routine, essentially voiding the above comment when*
*                   icon is 5                                          *
*    JY020699       modified subroutine to allow cavity to grow        *
*                   on both the back and face of the foil.             *
*                            IDR = 1  (back)                           *
*                            IDR = 2  (face)                           *
*    JY060500       Modified subroutine to make it more efficient.     *
*                                                                      *
************************************************************************
      INCLUDE 'PUFCAV.INC'
      COMMON/TI/NTSTEP,NREV,NTREV,NTPREV,IDXREV,ITSTEP

      AK=3.0
      IF(NWC(M,1).GT.0.OR.NWC(M,2).GT.0) AK=1.5

      AMAXMIN=0.05
      IF(ABS(DELTA(M,IDR)).GT.0.1)AMAXMIN=0.1         

      DCMAX=AMIN1(AMAXMIN,AK*ABS(DELTA(M,IDR)))

      IF(ABS(DELTA(M,IDR)).LE.DTOL) THEN
         DCAVL=SIGN(999.0,DELTA(M,IDR))
      ELSE         
         DELTAL=CAVL(M,IDR)-CAVLP(M,IDR)
         DELTAD=DELTA(M,IDR)-DELTAP(M,IDR)
         
         IF(DELTAL*DELTAD.GE.ZERO)THEN
            DCAVL=SIGN(999.0,DELTA(M,IDR))
         ELSE
            SLAB=DELTAL/DELTAD
            DCAVL=-SLAB*DELTA(M,IDR)
         END IF
      END IF 

      IF(ABS(DCAVL).GE.DCMAX)DCAVL=SIGN(DCMAX,DCAVL)

      CAVL(M,IDR)=CAVL(M,IDR)+DCAVL
      
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<<end of subroutine EXTNR>>>>>>>>>>>>>>>>>>>>>>>>>
      END







