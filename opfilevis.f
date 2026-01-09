
      SUBROUTINE OPFILEVIS
C*********************************************************************    
C        Subroutine to open file for viscous wetted results           *
C                Created by Hong Sun                                 * 
C********************************************************************* 

        OPEN(UNIT=10, FILE='blthk_p.dat', STATUS='UNKNOWN')
        OPEN(UNIT=11, FILE='blthk_s.dat', STATUS='UNKNOWN')

        OPEN(UNIT=12, FILE='p_distn.dat', STATUS='UNKNOWN')
        OPEN(UNIT=20, FILE='bl.wprs', STATUS='UNKNOWN')
!Allen Du 12/22/2017 output the pressure at the control points
        OPEN(UNIT=7041, FILE='control_bl.wprs', STATUS='UNKNOWN')
c XM YU 08/2012 write out cf 
        OPEN(UNIT=2011, FILE='cfskin.dat', STATUS='UNKNOWN')
c XM YU 08/2012
        OPEN(UNIT=22, FILE='Cdvis.dat', STATUS='UNKNOWN')

        OPEN(UNIT=121, FILE='bl_chk.dat', STATUS='UNKNOWN')

C        OPEN(UNIT=14, FILE='dfoil.dat', STATUS='UNKNOWN')
C      open(unit=97,file='suction_bl.dat',status='unknown')

        write(10,*) 'variables = "x","dstarp","thetap","hh","Ue","Uinv"'
        write(11,*) 'variables = "x","dstars","thetas","hh","Ue","Uinv"'
        write(12,*)'TITLE="Pressure Distribution"'
        write(12,*)'VARIABLES="x/c","-Cpvis","-Cpinv"'

        write(20,*)'TITLE="Viscous Pressure on Blade (2D)"'
        write(20,*)'VARIABLES = "Xd/C", "-Cp"'
!Allen Du 12/22/2017 output the pressure at the control points
        write(7041,*)'TITLE="Viscous Pressure on Blade (2D)"'
        write(7041,*)'VARIABLES = "x", "-cp"'

         write(22,*)'TITLE= "Section Drag Coefficient"'
         write(22,*)'VARIABLEs = "M", "r/R", "Cd"'

c        WRITE(14,*)'TITLE="Displacemented Foil Geo"'
c        write(14,*)' VARIABLES="x","y","x0","y0"'
C        write(97,*) 'variables = "x","dstars","thetas","hh","Ue","Uinv"'

       WRITE(121,*) 'variables = "I","UVL","SVL","XNVIS","YNVIS",
     &"XC","YC" '
!Allen Du 12/22/2017 add a file to check if the viscous run converged
       OPEN(UNIT=704, FILE='bl.convergence', STATUS='UNKNOWN')

       RETURN
       END
