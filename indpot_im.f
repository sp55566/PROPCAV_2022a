       SUBROUTINE INDPOT_IM
************************************************************************
*      This subroutine calculates the induced potential due to the     *
*      image panels and add them to the original influence coef's.     *
*                                                                      *
*      Date      Comment or Revision                                   *
*      --------  -------------------                                   *
*      JY011100  Subroutine created.                                   *
*                                                                      *
************************************************************************

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       INCLUDE 'PUFCAVC.INC'

       CHARACTER*50 FN141,FN142,FN210,FN211,FN212,FN213,FNWK
       CHARACTER*2 BLAIDX(10)
       DATA BLAIDX /'01','02','03','04','05','06','07','08','09','10'/

       IF(NTSTEP.EQ.1) THEN
          CALL CHRLEN(FNSCR,LENCH)

C........Blade image influence coefficients
          FN141=FNSCR(1:LENCH)//'S141.DAT'
          OPEN(141,FILE=FN141,STATUS='UNKNOWN',FORM='UNFORMATTED')
          FN142=FNSCR(1:LENCH)//'S142.DAT'
          OPEN(142,FILE=FN142,STATUS='UNKNOWN',FORM='UNFORMATTED')

C........Supercavity image influence coefficients
          FN210=FNSCR(1:LENCH)//'S210.DAT'
          OPEN(210,FILE=FN210,STATUS='UNKNOWN',FORM='UNFORMATTED')
          FN211=FNSCR(1:LENCH)//'S211.DAT'
          OPEN(211,FILE=FN211,STATUS='UNKNOWN',FORM='UNFORMATTED')
          FN212=FNSCR(1:LENCH)//'S212.DAT'
          OPEN(212,FILE=FN212,STATUS='UNKNOWN',FORM='UNFORMATTED')
          FN213=FNSCR(1:LENCH)//'S213.DAT'
          OPEN(213,FILE=FN213,STATUS='UNKNOWN',FORM='UNFORMATTED')

C........Wake image influence coefficients
          DO K=1,NBLADE
             IO=180+K
             FNWK=FNSCR(1:LENCH)//'SW3'//BLAIDX(K)//'.DAT'
             OPEN(IO,FILE=FNWK,STATUS='UNKNOWN',FORM='UNFORMATTED')
             IO=190+K
             FNWK=FNSCR(1:LENCH)//'SW4'//BLAIDX(K)//'.DAT'
             OPEN(IO,FILE=FNWK,STATUS='UNKNOWN',FORM='UNFORMATTED')
          END DO

C........Initialize the source and dipole strengths with zero
          DO I=1,NPANEL
             TEMP1(I)=ZERO
          END DO
          DO I=1,NPWAKS
             TEMP4(I)=ZERO
          END DO
          NREAD=NWMIN*MR
          DO I=1,NREAD
             TEMP5(I)=ZERO
          END DO
          NREC = 360 / NDLTAT
          DO N=1,NREC
             CALL WRITE2(45,N,TEMP1,NPANEL)
             CALL WRITE2(47,N,TEMP1,NPANEL)
             CALL WRITE2(48,N,TEMP4,NPWAKS)
             CALL WRITE2(46,N,TEMP5,NREAD)
          END DO
       END IF

       DT1=-DELTAT*(IDXREV-1)

C-----------------------------------------------------------------------
C      Calculate control points of the image panels on the blades.
C-----------------------------------------------------------------------
       DO J=1,NPANEL

C -- Rotating control points to local coord without inclination

          CALL ROTATE2(-1,SPANGLE,XCTP(J,1,1),XCTP(J,2,1))

C........rotating control points of the blade to global coord.

          XGI=XCTP(J,1,1)
          YGI=XCTP(J,2,1)*COS(DT1)-XCTP(J,3,1)*SIN(DT1)
          ZGI=XCTP(J,2,1)*SIN(DT1)+XCTP(J,3,1)*COS(DT1)

C........determine location of image in global coord.
                
          YGI=2.*YFS-YGI

C........rotating control points of image to local coord.

          XCTP_IM(J,1,1)=XGI
          XCTP_IM(J,2,1)=YGI*COS(DT1)+ZGI*SIN(DT1)
          XCTP_IM(J,3,1)=-YGI*SIN(DT1)+ZGI*COS(DT1)

C........cal. control points of images from other blades
          RCP=SQRT(XCTP_IM(J,2,1)**2.+XCTP_IM(J,3,1)**2.)
          THP=ATAN2(XCTP_IM(J,3,1),XCTP_IM(J,2,1))
          
          DO KK=2,NBLADE
             XCTP_IM(J,1,KK)=XCTP_IM(J,1,1)
             XCTP_IM(J,2,KK)=RCP*COS(THP+DELK*(KK-1))
             XCTP_IM(J,3,KK)=RCP*SIN(THP+DELK*(KK-1))
          END DO

C -- Change to local coordinate with inclination

          CALL ROTATE2(0,SPANGLE,XCTP(J,1,1),XCTP(J,2,1))

          DO KK = 1 , NBLADE
             YGI = XCTP_IM(J,2,KK) - 2.*YFS
             CALL ROTATE2(0,-SPANGLE,XCTP_IM(J,1,KK),YGI)
             XCTP_IM(J,2,KK) = YGI + 2. * YFS
          ENDDO

       END DO

       CALL INFCOF_IM

       DO I=1,NPANEL
          DO M=1,MR
             wsubif(i,m)=zero
          END DO
       END DO

       CALL INFWAK_IM

C-----------------------------------------------------------------------
C      Calculate control points of the image panels on the wakes.
C-----------------------------------------------------------------------
       DO J=1,NPWAKS

          CALL ROTATE2(-1,SPANGLE,XCPW(J,1,1),XCPW(J,2,1))

C........rotating control points of the blade to global coord.

          XGI=XCPW(J,1,1)
          YGI=XCPW(J,2,1)*COS(DT1)-XCPW(J,3,1)*SIN(DT1)
          ZGI=XCPW(J,2,1)*SIN(DT1)+XCPW(J,3,1)*COS(DT1)

C........determine location of image in global coord.                
          YGI=2.*YFS-YGI
          
C........rotating control points of image to local coord.

          XCPW_IM(J,1,1)=XGI
          XCPW_IM(J,2,1)=YGI*COS(DT1)+ZGI*SIN(DT1)
          XCPW_IM(J,3,1)=-YGI*SIN(DT1)+ZGI*COS(DT1)
          
C........cal. control points of images from other wakes
          RCP=SQRT(XCPW_IM(J,2,1)**2.+XCPW_IM(J,3,1)**2.)
          THP=ATAN2(XCPW_IM(J,3,1),XCPW_IM(J,2,1))

          DO KK=2,NBLADE
             XCPW_IM(J,1,KK)=XCPW_IM(J,1,1)
             XCPW_IM(J,2,KK)=RCP*COS(THP+DELK*(KK-1))
             XCPW_IM(J,3,KK)=RCP*SIN(THP+DELK*(KK-1))
          END DO

C -- Change to local coordinate with inclination

          CALL ROTATE2(0,SPANGLE,XCPW(J,1,1),XCPW(J,2,1))

          DO KK = 1 , NBLADE
             YGI = XCPW_IM(J,2,KK) - 2.*YFS
             CALL ROTATE2(0,-SPANGLE,XCPW_IM(J,1,KK),YGI)
             XCPW_IM(J,2,KK) = YGI + 2. * YFS
          ENDDO

       END DO
       
       CALL SCWKINF_IM

C-----------------------------------------------------------------------
C     read in influence coefficients for the key blade (and add on the
C     effect of image panels)
C-----------------------------------------------------------------------
C.....set W to be the infl. func. of key blade's first wake panel on..
C.....the key blade...................................................
       REWIND 81
       REWIND 181
       DO 20 M=MR,1,-1
          CALL READ1(81,TEMP1,NPANEL)
          CALL READ1(181,TEMP2,NPANEL)
          DO 10 I=1,NPANEL
             W(I,M)=TEMP1(I)-TEMP2(I)
 10       CONTINUE
 20    CONTINUE

C.....set WK to be the infl. func. of key blade's first wake panel on.
C.....the key blade wake..............................................
       REWIND 91
       REWIND 191
       DO 40 M=MR,1,-1
          CALL READ1(91,TEMP4,NPWAKS)
          CALL READ1(191,TEMP5,NPWAKS)
          DO 30 I=1,NPWAKS
             WK(I,M)=TEMP4(I)-TEMP5(I)
 30       CONTINUE
          DO 35 JJ=2,NSUB
             CALL READ1(91,TEMP4,NPWAKS)
             CALL READ1(191,TEMP5,NPWAKS)
 35       CONTINUE
 40    CONTINUE

C.....key blade influence coefficients AA (dipole) and BB (source)....
       REWIND 41
       REWIND 42
       REWIND 141
       REWIND 142
       DO 60 J=1,NPANEL
          CALL READ1(41,AA(1,J),NPANEL)
          CALL READ1(42,BB(1,J),NPANEL)
          CALL READ1(141,TEMP1,NPANEL)
          CALL READ1(142,TEMP2,NPANEL)
          DO I=1,NPANEL
             AA(I,J)=AA(I,J)-TEMP1(I)
             BB(I,J)=BB(I,J)-TEMP2(I)
          END DO
          DO 50 KK=2,NBLADE
             CALL READ1(41,TEMP1,NPANEL)
             CALL READ1(42,TEMP2,NPANEL)
             CALL READ1(141,TEMP1,NPANEL)
             CALL READ1(142,TEMP2,NPANEL)
 50       CONTINUE
 60    CONTINUE

C.....the next four files contain supercavitation ic's (all blades)...
       REWIND 110
       REWIND 111
       REWIND 112
       REWIND 113
       REWIND 210
       REWIND 211
       REWIND 212
       REWIND 213

       DO 80 J=1,NPWAKS
          CALL READ1(110,C(1,J),NPANEL)
          CALL READ1(113,F(1,J),NPWAKS)
          CALL READ1(210,TEMP1,NPANEL)
          CALL READ1(213,TEMP4,NPWAKS)
          DO I=1,NPANEL
             C(I,J)=C(I,J)-TEMP1(I)
          END DO
          DO I=1,NPWAKS
             F(I,J)=F(I,J)-TEMP4(I)
          END DO
          DO 70 KK=2,NBLADE
             CALL READ1(110,TEMP1,NPANEL)
             CALL READ1(113,TEMP4,NPWAKS)
             CALL READ1(210,TEMP1,NPANEL)
             CALL READ1(213,TEMP4,NPWAKS)
 70       CONTINUE
 80    CONTINUE
       DO 100 J=1,NPANEL
          CALL READ1(111,D(1,J),NPWAKS)
          CALL READ1(112,E(1,J),NPWAKS)
          CALL READ1(211,TEMP4,NPWAKS)
          CALL READ1(212,TEMP5,NPWAKS)
          DO I=1,NPWAKS
             D(I,J)=D(I,J)-TEMP4(I)
             E(I,J)=E(I,J)-TEMP5(I)
          END DO
          DO 90 KK=2,NBLADE
             CALL READ1(111,TEMP4,NPWAKS)
             CALL READ1(112,TEMP5,NPWAKS)
             CALL READ1(211,TEMP4,NPWAKS)
             CALL READ1(212,TEMP5,NPWAKS)
 90       CONTINUE
 100   CONTINUE


C =======================================================
C -- Plot to check image point
C =======================================================

       open(787,file='check_im.plt',status='unknown')
       write(787,6992)

       DO KK = 1, NBLADE
          WRITE(787,6993) kk,nc,mr
          do m = 1, mr
             do n = 1, nc
                L = indexb(n,m)
                write(787,*) xctp(l,1,kk),xctp(l,2,kk),xctp(l,3,kk)
             enddo
          enddo          
       enddo
       
       if(ihub .ne. 0) then
          DO KK = 1, NBLADE
             WRITE(787,6994) kk,nhbx,mhbt
             do m = 1, mhbt
                do n = 1, nhbx
                   L = indexh(n,m)
                   write(787,*) xctp(l,1,kk),xctp(l,2,kk),xctp(l,3,kk)
                enddo
             enddo         
          enddo
       endif

       DO KK = 1 , NBLADE
          WRITE(787,6997) kk,ntra,mr
          DO M=1,MR
             DO N=1,NTRA
                J=(MR-M)*NTRA+N
                write(787,*) xcpw(j,1,kk),xcpw(j,2,kk),xcpw(j,3,kk)
             enddo
          enddo
       enddo


       DO KK = 1, NBLADE
          WRITE(787,6995) kk,nc,mr
          do m = 1, mr
             do n = 1, nc
                L = indexb(n,m)
          write(787,*) xctp_im(l,1,kk),xctp_im(l,2,kk),xctp_im(l,3,kk)
             enddo
          enddo          
       enddo

       if(ihub .ne. 0) then
          DO KK = 1, NBLADE
             WRITE(787,6996) kk,nhbx,mhbt
             do m = 1, mhbt
                do n = 1, nhbx
                   L = indexh(n,m)
            write(787,*) xctp_im(l,1,kk),xctp_im(l,2,kk),xctp_im(l,3,kk)
                enddo
             enddo         
          enddo
       endif

       DO KK = 1 , NBLADE
          WRITE(787,6998) kk,ntra,mr
          DO M=1,MR
             DO N=1,NTRA
                J=(MR-M)*NTRA+N
           write(787,*) xcpw_im(j,1,kk),xcpw_im(j,2,kk),xcpw_im(j,3,kk)
             enddo
          enddo
       enddo

       close(787)

 6992 format('variables=x,y,z')
 6993 FORMAT('ZONE T="P=',i4,'",
     %            I=',I4,',J = ',I4,',F=POINT')
 6994 FORMAT('ZONE T="H=',i4,'",
     %            I=',I4,',J = ',I4,',F=POINT')

 6995 FORMAT('ZONE T="Pi=',i4,'",
     %            I=',I4,',J = ',I4,',F=POINT')
 6996 FORMAT('ZONE T="Hi=',i4,'",
     %            I=',I4,',J = ',I4,',F=POINT')

 6997 FORMAT('ZONE T="W=',i4,'",
     %            I=',I4,',J = ',I4,',F=POINT')
 6998 FORMAT('ZONE T="Wi=',i4,'",
     %            I=',I4,',J = ',I4,',F=POINT')

      go to 1111

C-----------------------------------------------------------------------
C     Check controls points of image panels.
C-----------------------------------------------------------------------
       ISW=0
       IF(ISW.EQ.0) GO TO 1111

 1000  FORMAT(1X,'ZONE T="K=',I1,' T=',I2,'", N=',I5,', E=',I5,',
     *      F=FEPOINT, ET=QUADRILATERAL')
 2000  FORMAT(3(1X,F10.6))
 3000  FORMAT(4(1X,I5))
       DO KK=1,NBLADE
          WRITE(999,1000) KK,IDXREV,MR*NC,(MR-1)*(NC-1)
          WRITE(899,1000) KK,IDXREV,MR*NC,(MR-1)*(NC-1)
       DO M=1,MR
             DO N=1,NC
                J=INDEXB(N,M)

                YGI=XCTP_IM(J,2,KK)*COS(DT1)-XCTP_IM(J,3,KK)*SIN(DT1)
                ZGI=XCTP_IM(J,2,KK)*SIN(DT1)+XCTP_IM(J,3,KK)*COS(DT1)
                WRITE(999,2000) XCTP_IM(J,1,KK),YGI,ZGI

                YGI=XCTP(J,2,KK)*COS(DT1)-XCTP(J,3,KK)*SIN(DT1)
                ZGI=XCTP(J,2,KK)*SIN(DT1)+XCTP(J,3,KK)*COS(DT1)
                WRITE(899,2000) XCTP(J,1,KK),YGI,ZGI

             END DO
          END DO          
          DO M=1,MR-1
             DO N=1,NC-1
                WRITE(999,3000) (M-1)*(NC)+N,(M-1)*(NC)+N+1,
     *               M*(NC)+N+1,M*(NC)+N
                WRITE(899,3000) (M-1)*(NC)+N,(M-1)*(NC)+N+1,
     *               M*(NC)+N+1,M*(NC)+N
             END DO
          END DO
          WRITE(998,1000) KK,IDXREV,MR*NTRA,(MR-1)*(NTRA-1)
          WRITE(898,1000) KK,IDXREV,MR*NTRA,(MR-1)*(NTRA-1)
          DO M=1,MR
             DO N=1,NTRA
                J=(MR-M)*NTRA+N

                YGI=XCPW_IM(J,2,KK)*COS(DT1)-XCPW_IM(J,3,KK)*SIN(DT1)
                ZGI=XCPW_IM(J,2,KK)*SIN(DT1)+XCPW_IM(J,3,KK)*COS(DT1)
                WRITE(998,2000) XCPW_IM(J,1,KK),YGI,ZGI

                YGI=XCPW(J,2,KK)*COS(DT1)-XCPW(J,3,KK)*SIN(DT1)
                ZGI=XCPW(J,2,KK)*SIN(DT1)+XCPW(J,3,KK)*COS(DT1)
                WRITE(898,2000) XCPW(J,1,KK),YGI,ZGI

             END DO
          END DO
          DO M=1,MR-1
             DO N=1,NTRA-1
                WRITE(998,3000) (M-1)*(NTRA)+N,(M-1)*(NTRA)+N+1,
     *               M*(NTRA)+N+1,M*(NTRA)+N
                WRITE(898,3000) (M-1)*(NTRA)+N,(M-1)*(NTRA)+N+1,
     *               M*(NTRA)+N+1,M*(NTRA)+N
             END DO
          END DO
       END DO

 1111  CONTINUE
C-----------------------------------------------------------------------

       RETURN
       END

