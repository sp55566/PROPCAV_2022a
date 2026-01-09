      SUBROUTINE QEXT2(S0,S1,S2,SL,QT1,QT2,QT3)
************************************************************************
*  3-POINT EXTRAPOLATION                                               *
*  Extrapolate the source strength beneath the cavity to the left split*
*  panel, using the following assumed behavior of dpdn:                *
*        dpdn(s)=a(sl-s)+b/sqrt(sl-s)+c                                *
*  where a,b,c may be written:                                         *
*        a=a1*dpdn(s1)+a2*dpdn(s2)+a3*dpdn(s3)                         *
*        b=b1*dpdn(s1)+b2*dpdn(s2)+b3*dpdn(s3)                         *
*        c=c1*dpdn(s1)+c2*dpdn(s2)+c3*dpdn(s3)                         *
*  thus qt1=a1*(sl-s0)+b1/sqrt(sl-s0)+c1                               *
*  and similarly for qt2,qt3,qt4.                                      *
*  06-24-92 NF                                                         *
************************************************************************
      Q1=1./SQRT(SL-S2)-1./SQRT(SL)
      Q2=1./SQRT(SL-S1)-1./SQRT(SL)
      T1=S1*Q1-S2*Q2
      A1=S2*Q1/T1
      A2=-1./S2-S1*Q1/T1
      A3=-1./S2-(S2-S1)*Q1/T1
      B1=-S2/T1
      B2=S1/T1
      B3=(S2-S1)/T1
      C1=-SL*S2*Q1/T1+S2/(SQRT(SL)*T1)
      C2=SL/S2+SL*S1*Q1/T1-S1/(SQRT(SL)*T1)
      C3=1.+SL/S2+SL*(S2-S1)*Q1/T1-(S2-S1)/(SQRT(SL)*T1)
      QT1=A1*(SL-S0)+B1/SQRT((SL-S0))+C1
      QT2=A2*(SL-S0)+B2/SQRT((SL-S0))+C2
      QT3=A3*(SL-S0)+B3/SQRT((SL-S0))+C3
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<<<<end of subroutine QEXT2>>>>>>>>>>>>>>>>>>>>>>>
      END
