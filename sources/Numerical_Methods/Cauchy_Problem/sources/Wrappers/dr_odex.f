C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR ODEX ON VAN DER POL'S EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
compile odex
cfeh dr_odex odex
        subroutine dr_odex
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (NDGL=2,KM=9,NRDENS=2,
     &     LWORK=NDGL*(KM+5)+5*KM+20+(2*KM*(KM+2)+5)*NRDENS,
     &     LIWORK=2*KM+21+NRDENS)
        DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK)
        DIMENSION  RPAR(1), RTOL(1), ATOL(1), IPAR(1)
        EXTERNAL FVPOL_odex,SOLOUT_odex
C --- DIMENSION OF THE SYSTEM
        N=2
        RPAR=5.
C --- OUTPUT ROUTINE AND DENSE OUTPUT IS USED DURING INTEGRATION
        IOUT=2
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=2.0D0
        Y(2)=0.0D0
C --- ENDPOINT OF INTEGRATION
        XEND=30.0D0
C --- REQUIRED (RELATIVE) TOLERANCE
        TOL=1.0D-9
        ITOL=0
        RTOL=TOL
        ATOL=TOL
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,10
        IWORK(I)=0
  10    WORK(I)=0.D0  
        H=0.01D0   
C --- IF DENSE OUTPUT IS REQUIRED
        IWORK(8)=NRDENS
        IWORK(21)=1
        IWORK(22)=2
        
        open(8, file='./results/vanderpol_odex.plt') 
C --- CALL OF THE SUBROUTINE DOPRI5  
        CALL ODEX(N,FVPOL_odex,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  SOLOUT_odex,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
        WRITE (6,99) X,Y(1),Y(2)
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10)
C --- PRINT STATISTICS
        WRITE (6,90) TOL
 90     FORMAT('       tol=',D8.2)
        WRITE (6,91) (IWORK(J),J=17,20)
 91     FORMAT(' fcn=',I5,' step=',I4,' accpt=',I4,' rejct=',I3)
        !STOP
        END
C
C
        SUBROUTINE SOLOUT_odex (NR,XOLD,X,Y,N,CON,NCON,ICOMP,ND,
     &                     RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
C --- BY USING "CONTD5", THE CONTINUOUS COLLOCATION SOLUTION
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Y(N),CON(NCON),ICOMP(ND)
        COMMON /INTERN/XOUT
        
        write(8,*) y(1), y(2) 
        
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),Y(2),NR-1
           XOUT=X+0.1D0
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              SOL1=CONTEX(1,XOUT,CON,NCON,ICOMP,ND)
              SOL2=CONTEX(2,XOUT,CON,NCON,ICOMP,ND)
              WRITE (6,99) XOUT,SOL1,SOL2,NR-1
              XOUT=XOUT+0.1D0
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10,'    NSTEP =',I4)
        RETURN
        END
C
        SUBROUTINE FVPOL_odex(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF VAN DER POL'S EQUATION
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION Y(N),F(N)
        EPS=RPAR
        F(1)=Y(2)
        F(2)=eps * (1 - y(1)**2) * y(2)  - y(1)
        RETURN
        END 

