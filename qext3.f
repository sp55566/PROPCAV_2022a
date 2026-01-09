      SUBROUTINE QEXT3(S0,S1,SL,QT1,QT2)
************************************************************************
*  2-POINT EXTRAPOLATION                                               *
*  Extrapolate the source strength beneath the cavity to the left split*
*  panel, using the following assumed behavior of dpdn:                *
*        dpdn(s)=a/dsqrt(sl-s)+B                                       *
*  where a,b may be written:                                           *
*        a=a1*dpdn(s1)+a2*dpdn(s2)                                     *
*        b=b1*dpdn(s1)+b2*dpdn(s2)                                     *
*  thus qt1=a1/sqrt(sl-s0)+b1                                         *
*  and similarly for qt2                                               *
*  06-24-92 NF                                                         *
************************************************************************
      SQSL=SQRT(SL)
      T1=1.-SQRT(SL-S1)/SQSL
      Q1=SQRT(SL-S1)/T1
      A1=Q1
      A2=-Q1
      B1=-Q1/SQSL
      B2=1.+Q1/SQSL
      QT1=A1/SQRT((SL-S0))+B1
      QT2=A2/SQRT((SL-S0))+B2
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<<<<end of subroutine QEXT3>>>>>>>>>>>>>>>>>>>>>>>
      END
