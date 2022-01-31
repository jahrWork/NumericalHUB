C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR DOPRI5 ON ARENSTORF ORBIT
C * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine dr_dopri5
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NDGL=4,NRDENS=2)
        PARAMETER (LWORK=8*NDGL+5*NRDENS+21,LIWORK=NRDENS+21)
        DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK),RPAR(2),IPAR(2)
        DIMENSION RTOL(NDGL),ATOL(NDGL)
        EXTERNAL FAREN,SOLOUT
        external dopri5
C --- DIMENSION OF THE SYSTEM
        N=NDGL
C --- OUTPUT ROUTINE (AND DENSE OUTPUT) IS USED DURING INTEGRATION
        IOUT=2
C --- INITIAL VALUES AND ENDPOINT OF INTEGRATION
        RPAR(1)=0.012277471D0
        RPAR(2)=1.D0-RPAR(1)
        X=0.0D0
        Y(1)=0.994D0
        Y(2)=0.0D0
        Y(3)=0.0D0
        Y(4)=-2.00158510637908252240537862224D0
        XEND=17.0652165601579625588917206249D0
C --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL=0
        RTOL=1.0D-7
        ATOL=RTOL
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.D0  
C --- DENSE OUTPUT IS USED FOR THE TWO POSITION COORDINATES 1 AND 2
        IWORK(5)=NRDENS
        IWORK(21)=1
        IWORK(22)=2
        open(8, file='./results/orbit.plt') 
        
C --- CALL OF THE SUBROUTINE DOPRI5   
        CALL DOPRI5(N,FAREN,X,Y,XEND,
     &              RTOL,ATOL,ITOL,
     &              SOLOUT,IOUT,
     &              WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
        WRITE (6,99) Y(1),Y(2)
 99     FORMAT(1X,'X = XEND     Y =',2E18.10)
C --- PRINT STATISTICS
        WRITE (6,91) RTOL,(IWORK(J),J=17,20)
 91     FORMAT('     tol=',D8.2,'   fcn=',I5,' step=',I4,
     &             ' accpt=',I4,' rejct=',I3)
       
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS BY USING "CONTD5"
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CON(5*ND),ICOMP(ND),RPAR(2)
        COMMON /INTERN/XOUT
        
        write(8,*) x, y(1), y(2) 
        
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),Y(2),NR-1
           XOUT=X+2.0D0
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              WRITE (6,99) XOUT,CONTD5(1,XOUT,CON,ICOMP,ND),
     &                     CONTD5(2,XOUT,CON,ICOMP,ND),NR-1
              XOUT=XOUT+2.0D0
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F6.2,'    Y =',2E18.10,'    NSTEP =',I4)
        RETURN
        END
C
        SUBROUTINE FAREN(N,X,Y,F,RPAR,IPAR)
C --- ARENSTORF ORBIT
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N),RPAR(2)
        AMU=RPAR(1)
        AMUP=RPAR(2)
        F(1)=Y(3)
        F(2)=Y(4)
        R1=(Y(1)+AMU)**2+Y(2)**2
        R1=R1*SQRT(R1)
        R2=(Y(1)-AMUP)**2+Y(2)**2
        R2=R2*SQRT(R2)
        F(3)=Y(1)+2*Y(4)-AMUP*(Y(1)+AMU)/R1-AMU*(Y(1)-AMUP)/R2
        F(4)=Y(2)-2*Y(3)-AMUP*Y(2)/R1-AMU*Y(2)/R2
        RETURN
        END 
