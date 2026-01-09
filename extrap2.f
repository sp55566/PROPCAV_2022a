       SUBROUTINE EXTRAP2(S0,S1,C1,C2,C3)
************************************************************************
*                                                                      *
*  Subroutine EXTRAP2 computes the constants C1,...,C3 used in the     *
*  phi0 EXTRAPolation.                                                 *
*                                                                      *
*  Author  Julie Young  Feb. 17, 1999                                  *
*                                                                      *
************************************************************************
       FAC=S1*S1-2.*S0*S1
       C2=(S1*S1-2.*S0*S1+S0*S0)/FAC
       C1=-S0*S0/FAC
       C3=(S1*S1*S0-S0*S0*S1)/FAC
       RETURN
       END
