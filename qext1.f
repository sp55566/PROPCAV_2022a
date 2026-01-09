      SUBROUTINE QEXT1(S0,S1,S2,S3,SL,QT1,QT2,QT3,QT4)
************************************************************************
*  4-POINT EXTRAPOLATION                                               *
*  Extrapolate the source strength beneath the cavity to the left split*
*  panel, using the following assumed behavior of dpdn:                *
*        dpdn(s)=a(sl-s)**2+b(sl-s)+c/sqrt(sl-s)+d                     *
*  where a,b,c,d may be written:                                       *
*        a=a1*dpdn(s1)+a2*dpdn(s2)+a3*dpdn(s3)+a4*dpdn(s4)             *
*        b=b1*dpdn(s1)+b2*dpdn(s2)+b3*dpdn(s3)+b4*dpdn(s4)             *
*        c=c1*dpdn(s1)+c2*dpdn(s2)+c3*dpdn(s3)+c4*dpdn(s4)             *
*        d=d1*dpdn(s1)+d2*dpdn(s2)+d3*dpdn(s3)+d4*dpdn(s4)             *
*  thus qt1=a1(sl-s0)**2+b1*(sl-s0)+c1/sqrt(sl-s0)+d1                  *
*  and similarly for qt2,qt3,qt4.                                      *
*  06-17-92 NF                                                         *
************************************************************************
      T1=S3*S2*(S3-S2)
      T2=S1*S2*(S2-S1)
      T3=S3*(S3-2.*SL)
      DSL=SQRT(SL)
      Q1=S2/SQRT(SL-S3)-S3/SQRT(SL-S2)+(S3-S2)/DSL
      Q2=S1/SQRT(SL-S2)-S2/SQRT(SL-S1)+(S2-S1)/DSL
      Q3=1./SQRT(SL-S3)-1./DSL
      CD=Q2*T1-Q1*T2
      CN=(S1-S2)*T1+(S3-S2)*T2
      A1=S2*Q1/CD
      A2=-S3/T1-(S1*T1+S3*T2)*Q1/(T1*CD)
      A3=(1.+T2*Q1/CD)*S2/T1
      A4=(S3-S2+Q1*CN/CD)/T1
      B1=S2/S3*(Q1*T3-Q3*T1)/CD
      B2=-T3/T1-(S1*T1/S3+T2)*(Q1*T3/T1-Q3)/CD
      B3=S2/S3*(T3/T1+(T2*T3*Q1/T1-T2*Q3)/CD)-1./S3
      B4=1./S3-(S2-S3)*T3/(S3*T1)+(CN*Q1*T3/(S3*T1)-CN*Q3/S3)/CD
      C1=-S2*T1/CD
      C2=(S1*T1+S3*T2)/CD
      C3=-S2*T2/CD
      C4=-CN/CD
      D1=(-SL*SL*Q1-SL/S3*(Q1*T3-T1*Q3)+T1/DSL)*S2/CD
      D2=SL*SL*((S1*T1+S3*T2)*Q1/(CD*T1)+S3/T1)+
     *      SL*((S1*T1/S3+T2)*(Q1*T3/T1-Q3)/CD+T3/T1)-
     *      (S1*T1+S3*T2)/(DSL*CD)
      D3=-SL*SL*(S2/T1+S2*T2*Q1/(T1*CD))-SL*(-1./S3+S2*T3/(S3*T1)
     *   +S2*T2*Q1*T3/(S3*T1*CD)-S2*T2*Q3/(S3*CD))+S2*T2/(DSL*CD)
      D4=1.-SL*SL*((S3-S2)/T1+CN*Q1/(CD*T1))
     *     -SL*(1./S3+(S3-S2)*T3/(S3*T1)+CN*Q1*T3/(S3*CD*T1)
     *     -CN*Q3/(S3*CD))+CN/(DSL*CD)
      QT1=A1*(SL-S0)*(SL-S0)+B1*(SL-S0)+C1/SQRT((SL-S0))+D1
      QT2=A2*(SL-S0)*(SL-S0)+B2*(SL-S0)+C2/SQRT((SL-S0))+D2
      QT3=A3*(SL-S0)*(SL-S0)+B3*(SL-S0)+C3/SQRT((SL-S0))+D3
      QT4=A4*(SL-S0)*(SL-S0)+B4*(SL-S0)+C4/SQRT((SL-S0))+D4
c      WRITE(*,*) 'Qext Called'
      RETURN
C<<<<<<<<<<<<<<<<<<<<<<<<<end of subroutine QEXT1>>>>>>>>>>>>>>>>>>>>>>>
      END

























