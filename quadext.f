      SUBROUTINE QUADEXT(S0,S1,S2,A1,A2,A3)
************************************************************************
*                                                                      *
*   QUADratic EXTrapolation constants  for extrapolation of any        *
*   function of f(s)                                                   *
*                                                                      *
*   *-----*-----*-----*                                                *
*   s0    s1    s2    s=0                                              *
*                 <---s                                                *
*   f(s0)=a1*f(s1)+a2*f(s2)+a3*f(0)                                    *
*                                                                      *
************************************************************************
c      SAVE
      A1=(S0*S0-S2*S0)/(S1*S1-S1*S2)
      A2=(S1*S0-S0*S0)/(S2*S1-S2*S2)
      A3=1.0+(S0*S0-S0*S2-S0*S1)/(S1*S2)
      RETURN
      END
