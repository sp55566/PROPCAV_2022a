       SUBROUTINE VECPLT
************************************************************************
*                                                                      *
*  Subroutine WETted VECtor plots out the wetted total velocity        *
*  vectors on the blade for the last wetted revolution.  The data      *
*  are arranged for TECPLOT format.                                    *
*                                                                      *
*                           ----------        ----------               *
*                          |          |      |          |              *
*                          |  PROPCAV |> --->|  WETVEC  |              *
*                          |          |      |          |-->           *
*                           -----^----        ----------   |           *
*                                |                         |           *
*                                 --<---------<------------v           *
*                                                                      *
*  Author: Julie Young 100598                                          *
*  Date:      Revision/comments                                        *
*  -----      ---------------                                          *
*  JY092099   Modified routine to also plot vectors in wake (subpanel) *
*             region.                                                  *
*  JY080200   Modified subroutine to plot vectors on the face and back *
*             side of the blades separately.  Also, vectors on wake    *
*             will no longer be plotted.                               *
*                                                                      *
************************************************************************

       INCLUDE 'PUFCAV.INC'
       INCLUDE 'PUFCAVB.INC'
       DIMENSION QTMP(NBZ,MBZ)

       CHARACTER*30 FNVEC       

 5000  FORMAT(1X,'TITLE="wetted vector plot"')
 5020  FORMAT(1X,'VARIABLES="X","Y","Z","U","V","W","QC"')
 5040  FORMAT(1x,'ZONE T="Time=',F8.1,'", I=',I5,', J=,',I5,
     *          'K=1, DATAPACKING=BLOCK,',
     *          'VARLOCATION=(4=CELLCENTERED,5=CELLCENTERED,',
     *                       '6=CELLCENTERED,7=CELLCENTERED)')

       ZOFF = 1.0

       IF(IDXREV.EQ.1.AND.IWET.EQ.1) THEN

          CALL CHRLEN(FN,LENCH)
          FNVEC=FN(1:LENCH)//'.wvec'
          OPEN(65,FILE=FNVEC,STATUS='UNKNOWN')

          WRITE(65,5000)
          WRITE(65,5020)
           
       END IF

C-----------------------------------------------------------------------
C      Plot velocity vectors on the blade
C-----------------------------------------------------------------------

       WRITE(65,5040) TT(IDXREV),NCP,MRP

       DO M = 1, MR
          DO N = 1, NC
             AA1 = UXTOT(N,M)**2 + UYTOT(N,M)**2+UZTOT(N,M)**2
             QTMP(N,M) = SQRT(AA1)
          ENDDO
       ENDDO
             
       WRITE(65,*) ((XB(N,M),N=1,NCP),M=1,MRP) 
       WRITE(65,*) ((YB(N,M),N=1,NCP),M=1,MRP) 
       WRITE(65,*) ((ZB(N,M),N=1,NCP),M=1,MRP) 
       WRITE(65,*) ((UXTOT(N,M),N=1,NC),M=1,MR)
       WRITE(65,*) ((UYTOT(N,M),N=1,NC),M=1,MR)
       WRITE(65,*) ((UZTOT(N,M),N=1,NC),M=1,MR)
       WRITE(65,*) ((QTMP(N,M),N=1,NC),M=1,MR)

       IF(IWET.EQ.1) THEN
          IF(NTSTEP.EQ.NTIME) CLOSE(65)
       ELSE IF(IWET.EQ.0) THEN
          IF(IDXREV.EQ.NTPREV) CLOSE(65)
          IF(ISTEADY.EQ.0) CLOSE(65)
       END IF

       RETURN
       END

