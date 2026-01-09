      SUBROUTINE EXTNR_SC(M)
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

      AK=1.5

      AMAXMIN=0.05
      IF(ABS(DELTAT1(M)).GT.0.1)AMAXMIN=0.1         

      DCMAX=AMIN1(AMAXMIN,AK*ABS(DELTAT1(M)))

      IF(ABS(DELTAT1(M)).LE.DTOL) THEN
         DCAVL=SIGN(999.0,DELTAT1(M))
      ELSE    
         IF(DELTATP(M).EQ.ZERO.AND.ABS(DELTAT1(M)).GT.ZERO) THEN
            DCAVL=SIGN(0.01,DELTAT1(M))
         ELSE     
            DELTAL=CAVLT(M)-CAVLTP(M)
            DELTAD=DELTAT1(M)-DELTATP(M)
         
            IF(DELTAL*DELTAD.GE.ZERO)THEN
               DCAVL=SIGN(999.0,DELTAT1(M))
            ELSE
               SLAB=DELTAL/DELTAD
               DCAVL=-SLAB*DELTAT1(M)
            END IF
         END IF
      END IF 

      IF(ABS(DCAVL).GE.DCMAX)DCAVL=SIGN(DCMAX,DCAVL)

      CAVLT(M)=CAVLT(M)+DCAVL
      
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<<end of subroutine EXTNR>>>>>>>>>>>>>>>>>>>>>>>>>
      END







