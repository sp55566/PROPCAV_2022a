
      SUBROUTINE OPFILEVIS_D
C*********************************************************************    
C        Subroutine to open file for viscous wetted results           *
C           For Duct Cease         Created by Hong Sun               *
C********************************************************************* 

        OPEN(UNIT=10, FILE='blthk_p.dat', STATUS='UNKNOWN')
        OPEN(UNIT=11, FILE='blthk_s.dat', STATUS='UNKNOWN')
        OPEN(UNIT=12, FILE='p_distn.dat', STATUS='UNKNOWN')
        OPEN(UNIT=13, FILE='CPBLD.dat', STATUS='UNKNOWN')
      
        OPEN(UNIT=20, FILE='bl.wprs', STATUS='UNKNOWN')
        OPEN(UNIT=22, FILE='Cdvis.dat', STATUS='UNKNOWN')

        OPEN(UNIT=121, FILE='bl_chk.dat', STATUS='UNKNOWN')
            
c        OPEN(UNIT=14, FILE='dfoil.dat', STATUS='UNKNOWN')
c      open(unit=97,file='suction_bl.dat',status='unknown')

        write(10,*) 'variables = "x","dstarp","thetap","hh","Ue","Uinv"'
        write(11,*) 'variables = "x","dstars","thetas","hh","Ue","Uinv"'
        write(12,*)'TITLE="Pressure Distribution"'
        write(12,*)'VARIABLES="x/c","-Cpvis","-Cpinv","H"'
      
        write(13,*)'TITLE="Pressure Distribution on Duct"'
        write(13,*)'VARIABLES="x/c","-Cpvis","-Cpinv"'
      
        write(20,*)'TITLE="Viscous Pressure on Duct (2D)"'
        write(20,*)'VARIABLES = "s/s_max","-Cp"'

c        WRITE(14,*)'TITLE="Displacemented Foil Geo"'
c        write(14,*)' VARIABLES="x","y","x0","y0"'
c        write(97,*) 'variables = "x","dstars","thetas","hh","Ue","Uinv"'

       WRITE(121,*) 'variables = "I","UVL","SVL","XNVIS","YNVIS",
     &"XC","YC" '

       RETURN
       END
