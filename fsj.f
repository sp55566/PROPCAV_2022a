      FUNCTION FSJ(S)
************************************************************************
*  function FSJ is the transition function of the cavity termination   *
*  model                                                               *
*  Author: Neal Fine    11-26-90                                       *
************************************************************************
      COMMON/STSL/ST,SL,RLAMDA1,AFJ
      AFJ=0.3
      RNU=1
      FSJ=0.0
      IF((S.GT.ST).AND.(S.LE.SL))FSJ=AFJ*((S-ST)/(SL-ST))**RNU
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<<end of function FSJ>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      END
