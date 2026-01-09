      SUBROUTINE EXTRAP(S0,S1,S2,C1,C2,C3,C4)
************************************************************************
*                                                                      *
*  Subroutine EXTRAP computes the constants C1,...,C4 used in the      *
*  phi0 EXTRAPolation.                                                 *
*                                                                      *
*  Author  Neal E. Fine  March 1990                                    *
*                                                                      *
*   Date of last revision                    Revision                  *
*   ---------------------              -------------------             *
*    1100  April 9,1990               -included this header            *
*                                     -beautified the code             *
*    1300  Feb 6, 1990                -remodeled for wing code         *
*                                                                      *
************************************************************************
      S02=S0*S0
      S12=S1*S1
      S22=S2*S2
      R2=S1*(S12-3.0*S02)/(S2*(S2*S2-3.0*S02))
      R1=S2*R2*(S2-2.0*S0)-S12+2.0*S0*S1
      T1=(S2-2.0*S0)/(S22-3.0*S02)
      T2=1.0/(S2*S22-3.0*S02*S2)
      C1=S02*(1.0-2.0*S0*T1)/R1
      C2=-S02*(2.0*S0*R1*T2-2.0*S0*T1*R2+R2)/R1
      C3=1.0+S02*(2.0*R1*T2*S0+2.0*T1*S0*(1.0-R2)+R2-1.0)/R1
      C4=(S1-S2*R2)/R1*(2.0*S02*S0*T1-S02)+2.0*S02*S0*S2*T2+S0
      RETURN
      END
