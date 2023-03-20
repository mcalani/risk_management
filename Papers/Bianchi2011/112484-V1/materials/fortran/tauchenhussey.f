C.HR MARKOV.FOR
C--------------------------------------------------------------------
C MAIN PROGRAM - MARKOV.FOR
C
C
C
C REFFERENCE:
C
C    TAUCHEN, GEORGE, 1987, "QUADRATURE-BASED METHODS FOR OBTAINING
C    APPROXIMATE SOLUTIONS TO NONLINEAR ASSET PRICING MODELS,"
C    DUKE UNIVERSITY.
C
C     Most of the routines below have beeen written by:
C    DR. A. RONALD GALLANT
C     DEPARTMANT OF STATISTICS
C     NORTH CAROLINA STATE UNIVERSITY

C THIS PROGRAM CALCULATES A DISCRETE APPROXIMATION TO A CONTINUOUS 
C GAUSSIAN VECTOR AUTREGRESSION
C
C Y(T) = B + A(1)*Y(T-1) + A(2)*Y(T-2) + ... + A(NLAG)*Y(T-NLAG) + E(T)
C
C WHERE
C
C    Y(T) IS AN (NVAR X 1) VECTOR
C    B IS AN (NVAR X 1) VECTOR
C    A(I) ARE (NVAR X NVAR) MATRICES
C    E(T) IS AN (NVAR X 1) VECTOR OF I.I.D. MULTIVARIATE NORMAL 
C       RANDOM VARIBLES WITH MEAN VECTOR ZERO AND VARIANCE SIGMA 
C 
C THREE FILES USED FOR I/O ARE DEFINED WITH OPEN STATEMENTS 
C 
C     OPEN(3,FILE='CON') 
C     OPEN(10,FILE='PARMS') 
C     OPEN(11,FILE='GHQUAD.DAT') 
C 
C UNIT 11 CONTAINS DISCRETE VALUES AND WEIGHTS FOR N-POINT 
C    GAUSSIAN QUADRATURE RULES UP TO 20 POINTS.  
C
C UNIT 10 CONTAINS THE FOLLOWING INPUT PARAMETER VALUES 
C 
C    NVAR     - NUMBER OF VARIABLES IN THE VAR 
C               FORMAT (I5) 
C    NLAG     - NUMBER OF LAGS IN THE VAR 
C               FORMAT (I5) 
C    NVAL     - NUMBER OF DISCRETE POINTS FOR EACH VARIABLE 
C               NVAR ROWS OF FORMAT (I5) 
C    B        - VECTOR OF INTERCEPTS IN THE VAR 
C               NVAR ROWS OF FORMAT (F10.0) 
C    AUTO     - HORIZONTALLY CONCATENATED MATRICES A(I) 
C               [ A(1) | A(2) | ... | A(NLAG) ]
C               THIS MATRIX IS READ COLUMNWISE
C               NVAR*NVAR*NLAG ROWS OF FORMAT (F10.0)
C    SIGMA    - VARIANCE MATRIX OF E(T) READ COLUMNWISE
C               NVAR*NVAR ROWS OF FORMAT (F10.0)
C
C UNIT 3 CONTAINS ALL OF THE OUTPUT.  THE VARIABLES NS AND 
C    NSTATE ARE DEFINED AS FOLLOWS
C
C    NS       - NUMBER OF DISCRETE VALUES THAT Y(T) CAN TAKE ON;          
C               EQUAL TO THE PRODUCT OF THE ELEMENTS OF THE MATRIX 
C               NVAL.  
C 
C    NSTATE   - NUMBER OF STATES IN THE SYSTEM; EQUAL TO NS**NLAG.
C               THE STATE IS DEFINED BY THE DISCRETE VALUES OF
C               Y(T-1), Y(T-2), ...,Y(T-NLAG); SINCE EACH LAG CAN
C               TAKE ON NS DIFFERENT VALUES, THERE ARE NS**NLAG
C               STATES IN THE SYSTEM.
C
C  EXAMPLE:  IF THERE ARE TWO VARIABLES (Y1 AND Y2) AND TWO
C            LAGS IN THE VAR AND TWO DISCRETE POINTS ARE USED FOR 
C            EACH VARIABLE, THE STATES ARE ARRANGED AS FOLLOWS:
C
C                                   DISCRETE VALUES
C                      
C            STATE NO.        TIME T-1            TIME T-2
C            ---------     ---------------     ---------------
C                           
C                1         Y1(1)     Y2(1)     Y1(1)     Y2(1)
C                2         Y1(2)     Y2(1)     Y1(1)     Y2(1)
C                3         Y1(1)     Y2(2)     Y1(1)     Y2(1)
C                4         Y1(2)     Y2(2)     Y1(1)     Y2(1)
C                5         Y1(1)     Y2(1)     Y1(2)     Y2(1)
C                6         Y1(2)     Y2(1)     Y1(2)     Y2(1)
C                7         Y1(1)     Y2(2)     Y1(2)     Y2(1)
C                8         Y1(2)     Y2(2)     Y1(2)     Y2(1)
C                9         Y1(1)     Y2(1)     Y1(1)     Y2(2)
C               10         Y1(2)     Y2(1)     Y1(1)     Y2(2)
C               11         Y1(1)     Y2(2)     Y1(1)     Y2(2)
C               12         Y1(2)     Y2(2)     Y1(1)     Y2(2)
C               13         Y1(1)     Y2(1)     Y1(2)     Y2(2)
C               14         Y1(2)     Y2(1)     Y1(2)     Y2(2)
C               15         Y1(1)     Y2(2)     Y1(2)     Y2(2)
C               16         Y1(2)     Y2(2)     Y1(2)     Y2(2)
C
C THREE OUTPUT MATRICES ARE RELEVANT
C
C    YMAT     - (NS X NVAR) MATRIX OF DISCRETE VARIABLE VALUES
C
C    PISTA    - (NSTATE X 1) VECTOR THAT IS THE STATIONARY
C               PROBABILITY DISTRIBUTION OF DISCRETE STATES
C
C    PIMAT    - (NSTATE X NS) MATRIX OF TRANSITION PROBABILITIES.  
C               IF NLAG>1, THIS MATRIX IS NOT SQUARE BECAUSE FROM
C               ANY GIVEN STATE, ONLY NS OTHER STATES ARE REACHABLE
C               IN THE NEXT PERIOD.
C
C ALSO AVAILABLE IN THE SUBROUTINE MCHAIN IS THE (NSTATE X NS)
C MATRIX IREACH WHICH CONTAINS THE REACHABLE STATE NUMBERS.  THUS
C ROW I OF IREACH CONTAINS THE NUMBERS OF THE NS STATES THAT ARE 
C REACHABLE IN THE NEXT PERIOD WHEN THE STATE NUMBER IN THE CURRENT
C PERIOD IS I.
C--------------------------------------------------------------------


       subroutine tauchenhussey(pimat,ymat2)

       use global


      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NVAR NLAG
       REAL*8  YMAT2(NVALN*NVALN,NVAR)
      INTEGER*4 NVAL(ZZNVAR),IWORK(ZZIW),IWORK1(ZZNST),IWORK2(ZZIW2)
      INTEGER*4 IWORK3(ZZIW3)
      REAL*8 THETA(ZZTHET),PIMAT(NVALN*NVALN,NVALN*NVALN)
      REAL*8 GHMAT(NRQUAD,3),WORK(ZZW)
      REAL*8 WORK1(ZZNST),WORK2(ZZNST)

       open(unit=10,file='paramvar.txt')         ! read output of  var

      OPEN(UNIT=3,FILE='COVfinal.TXT')
      OPEN(UNIT=11,FILE='GHQUAD.DAT')                 ! read files

      DO 5 I=1,NVAR

      NVAL(I)=NVALN

5     CONTINUE



      IB     = 1
      IAUTO  = IB               + NVAR
      ISIG   = IAUTO            + (NVAR**2)*NLAG
      ISIGNV = ISIG+NVAR**2
      IDETSG = ISIGNV+NVAR**2



      HH= NVAR+(NVAR**2)*(NLAG+1)


      DO 10 I=1,NVAR+(NVAR**2)*(NLAG+1)


      READ(10,*) THETA(I)

10    CONTINUE

      CALL HEADER('/B//_')
      CALL DGMPNT(THETA(IB),NVAR,1)
      CALL HEADER('/AUTO//_')
      CALL DGMPNT(THETA(IAUTO),NVAR,NVAR*NLAG)
      CALL HEADER('/SIG//_')
      CALL DGMPNT(THETA(ISIG),NVAR,NVAR)

C CALCULATE THE NUMBER OF STATES

      NS=1
      DO 40 J=1,NVAR
40    NS=NS*NVAL(J)

C ERROR TRAPS
C
      IF (NVAR .GT. ZZNVAR) THEN
         WRITE(3,1081) NVAR
         STOP
      ENDIF

      IF (NLAG .GT. ZZNLAG) THEN
         WRITE(3,1082) NLAG
         STOP
      ENDIF

      IF (NS .GT. ZZNS) THEN
         WRITE(3,1083) NS
      print*,ns, 'increase zzns to',ns
         STOP
      ENDIF

      IF ( (4*(NVAR**2)+4*NVAR) .GT. ZZW) THEN
         IFIX=(4*(NVAR**2)+4*NVAR)-ZZW
         WRITE(3,1084) IFIX
         STOP
      ENDIF
      IF (ZZW .GT. 8000) THEN
         WRITE(3,1085)
         STOP
      ENDIF
C
C READ QUADRATURE DATA POINTS
C

      DO 20 I=1,NRQUAD
   20 READ(11,1101) (GHMAT(I,J), J=1,3)
C
C COMPUTE THE INVERSE OF SIGMA
C
      DO 30 I=1,NVAR**2
   30 THETA(ISIGNV+I-1)=THETA(ISIG+I-1)
      CALL DSWEEP(THETA(ISIGNV),NVAR,1.D-13,IER)
      CALL HEADER('/SIGINV//_')
      CALL DGMPNT(THETA(ISIGNV),NVAR,NVAR)
      IF (IER .GT. 0) THEN
      print*,'error tauchen transpose'
         WRITE(3,1090) IER
         STOP
      ENDIF

C     print*,'jhere'

C POINTERS TO WORK FOR CALL TO MCHAIN
C
      MR1=1
      MR2=MR1+NS*NVAR
      MR3=MR2+NS
      MR4=MR3+NVAR
      MR5=MR4+NS*NVAR
C
C POINTERS TO IWORK FOR CALL TO MCHAIN
C
      MI1=1
      MI2=MI1+NLAG
      MI3=MI2+NS

      CALL MCHAIN(THETA(IAUTO),THETA(IB),THETA(ISIG),NLAG,NVAL,
     &  NS,GHMAT,NRQUAD,NVAR,WORK(MR1),WORK(MR2),IWORK3,
     &  IWORK(MI1),WORK(MR3),IWORK(MI2),ymat2,IWORK1,
     &  PIMAT,WORK1,THETA(ISIGNV),THETA(IDETSG),WORK2,WORK(MR5),
     &  IWORK2,IWORK(MI3))
     

 1001 FORMAT(I5)
 1002 FORMAT(F10.0)
 1101 FORMAT(F10.0,2F20.0)
 1081 FORMAT(1X,'ERROR:  INCREASE ZZNVAR TO ',I5)
 1082 FORMAT(1X,'ERROR:  INCREASE ZZNLAG TO ',I5)
 1083 FORMAT(1X,'ERROR:  INCREASE ZZNS TO ',I5)
 1084 FORMAT(1X,'ERROR:  INCREASE ZZW BY ',I5)
 1085 FORMAT(1X,'ERROR:  DECREASE ZZNS OR COMPILE UNDER HUGE MODULE')
 1090 FORMAT(1X,'ERROR:  SIGMA IS SINGULAR:  IER = ',I5)

      close(unit=10)
      close(unit=11)
      close(unit=3)


      END     subroutine
C.HR SUBROUTINE MCHAIN 
C@
C--------------------------------------------------------------------
C SUBROUTINE MCHAIN                                                  
C--------------------------------------------------------------------



            SUBROUTINE MCHAIN(AUTO,B,SIGMA,NLAG,NVAL,
     &      NS,GHMAT,NRGH,NVAR,   ZMAT,WVEC,IREACH,
     &      IS,EXVAL,IR,YMAT,ISEQA,
     &      PIMAT,PISTA,SIGINV,
     &      DETSIG,WORK1,WORK,
     &      IWORK1,IWORK)

      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 NVAL(NVAR),IS(NLAG),IR(NS),IREACH(NS**NLAG,NS)
      INTEGER*4 ISEQA(NS**NLAG),IWORK1(1),IWORK(1)
      REAL*8 AUTO(NVAR,NVAR*NLAG),B(NVAR),SIGMA(NVAR,NVAR)
      REAL*8 GHMAT(NRGH,3),ZMAT(NS,NVAR),WVEC(NS),EXVAL(NVAR)
      REAL*8 YMAT(NS,NVAR),PIMAT(NS**NLAG,NS),PISTA(NS**NLAG)
      REAL*8 SIGINV(NVAR,NVAR),WORK1(1),WORK(1)
C--------------------------------------------------------------------
C SUBROUTINE TO CALCULATE DISCRETE MARKOV CHAIN THAT APPROXIMATES VAR
C--------------------------------------------------------------------
C SET QUADRATURE ARRAYS                                               
C--------------------------------------------------------------------
      CALL QUAD(NVAL,NVAR,NS,ZMAT,WVEC,GHMAT,NRGH)
      CALL HEADER('/ZMAT//_')
      CALL DGMPNT(ZMAT,NS,NVAR)
      CALL HEADER('/WVEC//_')
      CALL DGMPNT(WVEC,NS,1)
C--------------------------------------------------------------------
C SET ADDITIONAL CONSTANTS AND ARRAYS OF PARAMETERS                  
C--------------------------------------------------------------------
      NSTATE = NS**NLAG
      CALL SPACER(2)
      WRITE(3,100) NVAR,NLAG,NS,NSTATE
C--------------------------------------------------------------------
C FORM THE REACHABLE STATES AND PUT THEM IN IREACH                    
C--------------------------------------------------------------------
      CALL SPACER(1)
      WRITE(3,110)
      DO 20 J=1,NSTATE
         CALL DECODE(J,IS,1,NS,NLAG)
         CALL REACHR(IS,IR,IWORK1,NS,NLAG)
         DO 10 K=1,NS
   10    IREACH(J,K)= IR(K)
         ISEQA(J)=J
   20 CONTINUE
      CALL SPACER(1)
      WRITE(3,120)
C--------------------------------------------------------------------
C CALCULATE THE MEAN OF THE PROCESS                                  
C--------------------------------------------------------------------
      IR1=1
      IR2=IR1+NVAR**2
      IR3=IR2+NVAR**2
      CALL ARMEAN(AUTO,B,EXVAL,WORK(IR1),WORK(IR2),WORK(IR3),NVAR,NLAG)
      CALL HEADER('/EXVAL//_')
      CALL DGMPNT(EXVAL,1,NVAR)          
C--------------------------------------------------------------------
C TRANSFORM TO THE Y'S AND ADD IN THE MEANS                          
C--------------------------------------------------------------------
      CALL HEADER('/SIGMA//_')
      CALL DGMPNT(SIGMA,NVAR,NVAR)
      CALL ZTOY(ZMAT,YMAT,WORK(IR1),EXVAL,SIGMA,NVAR,DETSIG,NS,
     &          WORK(IR2))
      CALL HEADER('/YMAT//_')
      CALL DGMPNT(YMAT,NS,NVAR)

C--------------------------------------------------------------------
C GET TRANSITION PROBABILITIES                                       
C--------------------------------------------------------------------


      CALL SPACER(2)
      WRITE(3,130)
      IR2=IR1+NVAR*NLAG
      IR3=IR2+NS
      CALL PISUB(ISEQA,NSTATE,YMAT,PIMAT,WORK(IR1),IWORK1,
     &           WORK(IR2),NS,NVAR,NLAG,WVEC,EXVAL,SIGINV,DETSIG,AUTO,
     &           B,WORK(IR3))
      CALL SPACER(1)
      WRITE(3,140)
      CALL HEADER('/PIMAT//_')
      CALL DGMPNT(PIMAT,NSTATE,NS)

C--------------------------------------------------------------------
C COMPUTE THE CHAIN'S STATIONARY DISTRIBUTION                        
C--------------------------------------------------------------------
      CALL SPACER(1)
      WRITE(3,150)
      CALL SPACER(1)
      II1=1
      II2=II1+NLAG
      CALL MARKST(PIMAT,PISTA,WORK1,IWORK(II1),IWORK1,IWORK(II2),
     &            NSTATE,NS,NLAG)
      CALL SPACER(1)
      WRITE(3,160)
      CALL HEADER('/PISTA//_')
      CALL DGMPNT(PISTA,NSTATE,1)

      RETURN
  100 FORMAT(1X,'NVAR =',I5,'   NLAG =',I5,'   NS =',I5,
     &          '   NSTATE =',I5)      
  110 FORMAT(1X,'FORMING THE REACHABLE STATES')
  120 FORMAT(1X,'REACHABLE STATES FORMED')
  130 FORMAT(1X,'STARTING TO COMPUTE TRANSITION PROBABILITIES') 
  140 FORMAT(1X,'TRANSITION PROBABILITIES COMPUTED')
  150 FORMAT(1X,'STARTING TO COMPUTE STATIONARY DISTRIBUTION')
  160 FORMAT(1X,'STATIONARY DISTRIBUTION COMPUTED')
      END  
C.HR SUBROUTINE PHIMAT 
C.@
      SUBROUTINE PHIMAT(XMAT,MU,SIGINV,FZMAT,XCTR,XSIG,DETSIG,N,NVAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XMAT(N,NVAR),MU(NVAR),SIGINV(NVAR,NVAR)
      REAL*8 XCTR(N,NVAR),XSIG(N,NVAR),FZMAT(N)
      DATA PI /3.141592653589798/
C--------------------------------------------------------------------
C SUBROUTINE FOR ROW-WISE MULTIVARIATE NORMAL DISTRIBUTION                  
C
C    NOTES: THE ELEMENTS OF RETURNED VECTOR ARE THE MULTIVARIATE     
C           NORMAL DENSITY N(MU,SIGMA) EVALUATED AT EACH ROW OF XMAT 
C--------------------------------------------------------------------
      DO 10 I=1,N
      DO 10 J=1,NVAR
   10 XCTR(I,J)=XMAT(I,J)-MU(J)
      CALL DGMPRD(XCTR,SIGINV,XSIG,N,NVAR,NVAR)
      DO 30 I=1,N
         FZMAT(I)=0
         DO 20 J=1,NVAR
   20    FZMAT(I)=FZMAT(I)+XSIG(I,J)*XCTR(I,J)
         FZMAT(I)=EXP(-0.5*FZMAT(I))/((2*PI*DETSIG)**.5)
   30 CONTINUE
      RETURN
      END
C.HR SUBROUTINE CODE 
C@
C--------------------------------------------------------------------
C SUBROUTINE CODE                                                    
C--------------------------------------------------------------------
      SUBROUTINE CODE(IMAT,ICODE,NR,NC,NBASE)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 IMAT(NR,NC),ICODE(NR)
C--------------------------------------------------------------------
C SUBROUTINE TO CODE EACH ROW OF IMAT AS                             
C                                                                    
C        X = N1 + (N2-1)*NBASE + ... + (NL-1)*NBASE(L-1)               
C                                                                    
C WHERE (N1,N2,...,NL) IS A ROW OF IMAT AND L = COLS(IMAT).          
C--------------------------------------------------------------------
C CHECK FOR ERRORS                                                   
C--------------------------------------------------------------------
      IF (NBASE .LT. 1) THEN
        WRITE(3,100) NBASE
        STOP
      ENDIF
      DO 20 I=1,NR
       DO 10 J=1,NC
        IF (IMAT(I,J) .GT. NBASE .OR. IMAT(I,J) .LE. 0) THEN
         WRITE(3,110) 
         STOP
        ENDIF
   10  CONTINUE
   20 CONTINUE
C--------------------------------------------------------------------
C START MAIN PART OF SUBROUTINE                                      
C--------------------------------------------------------------------
      DO 40 I=1,NR
         ICODE(I)=IMAT(I,1)
         DO 30 J=2,NC
   30    ICODE(I)=ICODE(I)+(IMAT(I,J)-1)*(NBASE**(J-1))
   40 CONTINUE
      RETURN
  100 FORMAT(1X,'*** CODE ERROR: BAD NBASE',I5)
  110 FORMAT(1X,'*** CODE ERROR: AN ELEMENT OF THE INPUT ARRAY EITHER',
     &       1X,'  EXCEEDS THE NBASE OR IS LESS THAN OR EQUAL TO ZERO')
      END
C.HR SUBROUTINE DECODE 
C@
      SUBROUTINE DECODE(IX,ISTATE,NR,NBASE,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 IX(NR),ISTATE(NR,L)
C--------------------------------------------------------------------
C SUBROUTINE TO DECODE THE ROWS OF IX WHERE A TYPICAL ROW IS OF THE 
C FORM
C                                                                    
C       IX(K,.) =  N1 + (N2-1)*NBASE**1 + ... (NL-1)*NBASE**(L-1)        
C--------------------------------------------------------------------
C CHECK FOR ERRORS                                                   
C--------------------------------------------------------------------
      IF (NBASE .LE. 0 .OR. L .LT. 1) THEN
         WRITE(3,100) NBASE,L
         STOP
      ENDIF
      DO 10 I=1,NR
         IF (IX(I) .LE. 0) THEN
            WRITE(3,110)
         ENDIF  
   10 CONTINUE
C--------------------------------------------------------------------
C START SUBROUTINE                                                   
C--------------------------------------------------------------------
      DO 20 I=1,NR
      DO 20 J=1,L
   20 ISTATE(I,J)=INT(MOD(AINT(REAL(IX(I)-1)/REAL(NBASE**(J-1))),
     &                          REAL(NBASE)))+1
      RETURN
  100 FORMAT(1X,'*** DECODE ERROR BAD NBASE = ',I5,
     &          '  OR BAD LENGTH = ',I5)
  110 FORMAT(1X,'*** DECODE WARNING: AT LEAST ONE ELEMENT OF INPUT',
     &       1X,'    VECTOR IS LESS THAN OR EQUAL TO ZERO')
      END
C.HR SUBROUTINE REACHR 
C@
C--------------------------------------------------------------------
C              SUBROUTINE REACHR                                     
C--------------------------------------------------------------------
      SUBROUTINE REACHR(IVEC,IY,IX,NBASE,NLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 IVEC(NLAG),IX(NBASE,NLAG),IY(NBASE)
      DO 20 I=1,NBASE
         IX(I,1)=I
         DO 10 J=2,NLAG
   10    IX(I,J)=IVEC(J-1)
   20 CONTINUE
      CALL CODE(IX,IY,NBASE,NLAG,NBASE)
      RETURN
      END
C.HR SUBROUTINE PISUB
C.@
C--------------------------------------------------------------------
C SUBROUTINE FOR SELECTED ROWS OF PROBABILITY TRANSITION MATRIX      
C--------------------------------------------------------------------
      SUBROUTINE PISUB(JVEC,NRJ,YMAT,PIOUT,X,IST,PIOUTT,NS,NVAR,NLAG,
     &           WVEC,EXVAL,SIGINV,DETSIG,AUTO,B,WORK)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 JVEC(NRJ),IST(NRJ,NLAG)
      REAL*8 YMAT(NS,NVAR),PIOUT(NRJ,NS),X(NVAR*NLAG)
      REAL*8 PIOUTT(1,NS),WVEC(NVAR**NLAG),EXVAL(NVAR)
      REAL*8 SIGINV(NVAR,NVAR),AUTO(NVAR,NVAR*NLAG),B(NVAR),WORK(1)
      I1=1
      I2=I1+NS
      I3=I2+NS
      I4=I3+NVAR*NLAG
      I5=I4+NVAR
      CALL DECODE(JVEC,IST,NRJ,NS,NLAG)
      DO 40 I=1,NRJ
         DO 10 L=1,NLAG
         DO 10 K=1,NVAR
   10    X((L-1)*NVAR+K)=YMAT(IST(I,L),K)
         CALL PIFN(X,1,PIOUTT,YMAT,NS,NVAR,WVEC,EXVAL,SIGINV,
     &             DETSIG,AUTO,B,WORK(I1),WORK(I2),WORK(I3),
     &             WORK(I4),NLAG,WORK(I5))
         VSUM=0
         DO 20 J=1,NS
   20    VSUM=VSUM+PIOUTT(1,J)
         DO 30 J=1,NS
   30    PIOUT(I,J)=PIOUTT(1,J)/VSUM
   40 CONTINUE 
      RETURN
      END
C.HR SUBROUTINE ARMEAN 
C@
      SUBROUTINE ARMEAN(AUTO,B,MEAN,SUMAUT,SUMINV,WORKSP,NVAR,NLAG)
      REAL*8 AUTO(NVAR,NLAG*NVAR),B(NVAR),MEAN(NVAR)
      REAL*8 SUMAUT(NVAR,NVAR),SUMINV(NVAR,NVAR)
      REAL*8 WORKSP(2*NVAR**2+4*NVAR)
C--------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE MEAN OF AN AUTOREGRESSIVE PROCESS      
C--------------------------------------------------------------------
C START MAIN PART OF SUBROUTINE                                      
C--------------------------------------------------------------------
      DO 30 L=1,NLAG
         DO 20 I=1,NVAR
            DO 10 J=1,NVAR
               IF (L .EQ. 1) SUMAUT(I,J)=0
               SUMAUT(I,J)= SUMAUT(I,J)+AUTO(I,(J+(L-1)*NVAR))
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      DO 50 I=1,NVAR
         DO 40 J=1,NVAR
            SUMAUT(I,J)=-SUMAUT(I,J)
            IF (I .EQ. J) SUMAUT(I,J)=1+SUMAUT(I,J)
   40    CONTINUE
   50 CONTINUE
      CALL DMPINV(SUMAUT,NVAR,NVAR,SUMINV,IRANK,WORKSP)
      CALL DGMPRD(SUMINV,B,MEAN,NVAR,NVAR,1)
      RETURN
      END
C.HR SUBROUTINE QUAD 
C@
C--------------------------------------------------------------------
C SUBROUTINE QUAD                                                    
C--------------------------------------------------------------------
      SUBROUTINE QUAD(NVAL,NVAR,NS,ZMAT,WVEC,QUADAT,NRQUAD)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 NVAL(NVAR)
      REAL*8 QUADAT(NRQUAD,3),ZMAT(NS,NVAR),WVEC(NS)
C--------------------------------------------------------------------
C SUBROUTINE TO PUT ALL POSSIBLE COMBINATIONS OF QUADRATURE ABSCISSA 
C    AND CORRESPONDING PRODUCTS OF QUADRATURE WEIGHTS IN             
C    RETURNED MATRIX                                                 
C--------------------------------------------------------------------
C START MAIN SUBROUTINE                                              
C--------------------------------------------------------------------
      NPCUM=1
      nout=3
      DO 10 I=1,NS
   10 WVEC(I)=1
      DO 50 I=1,NVAR
       NP=NVAL(I)
       DO 40 J=1,NP
        Z=QUADAT(((NP-1)*NP)/2+J,2)
        W=QUADAT(((NP-1)*NP)/2+J,3)
        DO 30 K=1,(NS/NP)/NPCUM
         DO 20 L=1,NPCUM
          ZMAT((K-1)*NPCUM*NP+(J-1)*NPCUM+L,I)=Z
          WVEC((K-1)*NPCUM*NP+(J-1)*NPCUM+L)=WVEC((K-1)*NPCUM*NP+
     &       (J-1)*NPCUM+L)*W
   20    CONTINUE
   30   CONTINUE
   40  CONTINUE
       NPCUM=NPCUM*NP
   50 CONTINUE
      RETURN
      END 
C.HR SUBROUTINE CHOL
C@
      SUBROUTINE CHOL(S,M,XX,CSTAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 S(M,M),XX((M*(M+1))/2),CSTAT(M,M)
C
C     SUBROUTINE TO FIND THE CHOLESKY FACTORIZATION OF A MATRIX S
C
C     SUBROUTINE CALLED:  DMFSD
C
      DO 10 J=1,M
      DO 10 I=1,J
   10 XX(I+(J*(J-1))/2)=S(I,J)
      CALL DMFSD(XX,M,1.D-13,IER)
      IF (IER.LT.0) THEN
      WRITE(NOUT,100) IER
      STOP
      ENDIF
      DO 30 J=1,M
      DO 20 I=1,M
      CSTAT(I,J)=0
      IF (I .LE. J) CSTAT(I,J)=XX(I+(J*(J-1))/2)
   20 CONTINUE
   30 CONTINUE
      RETURN
  100 FORMAT(1X,'ERROR, CALL TO DMFSD, IER =',I5)
      END
C.HR SUBROUTINE ZTOY 
C@
      SUBROUTINE ZTOY(ZMAT,YMAT,S,EXVAL,SIGMA,NVAR,DETSIG,NS,WORKSP)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ZMAT(NS,NVAR),YMAT(NS,NVAR),S(NVAR,NVAR),EXVAL(NVAR)
      REAL*8 SIGMA(NVAR,NVAR),WORKSP(1)
C--------------------------------------------------------------------
C SUBROUTINE TO TRANSFORM THE Z'S TO THE Y'S                         
C--------------------------------------------------------------------
C START MAIN PART OF SUBROUTINE                                      
C--------------------------------------------------------------------
      CALL CHOL(SIGMA,NVAR,WORKSP,S)
      CALL DGMPRD(ZMAT,S,YMAT,NS,NVAR,NVAR)
      DO 20 I=1,NS
         DO 10 J=1,NVAR
            YMAT(I,J)=YMAT(I,J)+EXVAL(J)
   10    CONTINUE
   20 CONTINUE
      DETSIG=1
      DO 30 J=1,NVAR
         DETSIG=DETSIG*S(J,J)
   30 CONTINUE
      RETURN    
      END
C.HR SUBROUTINE PIFN 
C.@
C--------------------------------------------------------------------
C SUBROUTINE FOR TRANSITION PROBABILITES AS A FUNCTION OF X                
C--------------------------------------------------------------------
      SUBROUTINE PIFN(XMAT,NRXMAT,PIVAL,YMAT,NS,NVAR,WVEC,EXVAL,
     &           SIGINV,DETSIG,AUTO,B,TEMP,FZMAT,X,CONEXP,
     &           NLAG,WORK)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XMAT(NRXMAT,NVAR*NLAG),PIVAL(NRXMAT,NS),YMAT(NS,NVAR)
      REAL*8 WVEC(NS),EXVAL(NVAR),SIGINV(NVAR,NVAR)
      REAL*8 AUTO(NVAR,NVAR*NLAG),B(NVAR),TEMP(NS),FZMAT(NS)
      REAL*8 X(NVAR*NLAG),CONEXP(NVAR),WORK(1)
      I1=1
      I2=I1+NS*NVAR
      CALL PHIMAT(YMAT,EXVAL,SIGINV,FZMAT,WORK(I1),WORK(I2),DETSIG,
     &            NS,NVAR)
      DO 10 I=1,NS
   10 TEMP(I)= WVEC(I)/FZMAT(I)
      DO 60 J=1,NRXMAT
         DO 20 K=1,NVAR*NLAG
   20    X(K)= XMAT(J,K)
         DO 40 L=1,NVAR
            CONEXP(L)=0
            DO 30 K=1,NVAR*NLAG
   30       CONEXP(L)=CONEXP(L) + X(K)*AUTO(L,K)
            CONEXP(L)=CONEXP(L)+B(L)
   40    CONTINUE
         CALL PHIMAT(YMAT,CONEXP,SIGINV,FZMAT,WORK(I1),WORK(I2),
     &               DETSIG,NS,NVAR)
         DO 50 I=1,NS
   50    PIVAL(J,I)=FZMAT(I)*TEMP(I)
   60 CONTINUE
      RETURN
      END
C.HR SUBROUTINE MARKST 
C@
      SUBROUTINE MARKST(PIMAT,PISTA,TEMP,IST,ISTOLD,JOLD,NSTATE,NS,NLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 IST(NLAG),ISTOLD(NS,NLAG),JOLD(NS)
      REAL*8 PIMAT(NSTATE,NS),PISTA(NSTATE),TEMP(NSTATE)
C--------------------------------------------------------------------
C SUBROUTINE TO CALCULATE THE STATIONARY DISTRIBUTION OF THE MARKOV 
C CHAIN  
C--------------------------------------------------------------------
C INTITIALIZE                                                        
C--------------------------------------------------------------------
      DATA EPSTA /.0001/
      DATA MAXIT /50/
      CRIT=EPSTA+1     
      NBASE = NS
      DO 10 J=1,NSTATE
   10 PISTA(J)=1/REAL(NSTATE)
C--------------------------------------------------------------------
C MAIN LOOP TO COMPUTE THE CHAIN'S STATIONARY DISTRIBUTION           
C--------------------------------------------------------------------
      DO 90 NUMIT=1,MAXIT
         IF (CRIT .GT. EPSTA) THEN
            IF (NLAG .EQ. 1) THEN
               DO 30 J=1,NSTATE
                  TEMP(J)=0
                  DO 20 K=1,NSTATE
   20             TEMP(J)=TEMP(J) + PISTA(K)*PIMAT(K,J)
   30          CONTINUE
            ELSE
               DO 70 J=1,NSTATE
                  CALL DECODE(J,IST,1,NBASE,NLAG)
                  DO 50 K=1,NS
                     DO 40 L=1,NLAG-1
   40                ISTOLD(K,L)=IST(L+1)
                     ISTOLD(K,NLAG)=K
   50             CONTINUE
                  CALL CODE(ISTOLD,JOLD,NS,NLAG,NBASE)
                  SUM=0
                  DO 60 K=1,NS
   60             SUM=SUM+PIMAT(JOLD(K),IST(1))*PISTA(JOLD(K))
                  TEMP(J)=SUM
   70          CONTINUE
            ENDIF
            CRIT=0
            DO 80 J=1,NSTATE
               DELTA=ABS(TEMP(J)-PISTA(J))
               IF (DELTA .GT. CRIT) CRIT=DELTA
               PISTA(J)=TEMP(J)
   80       CONTINUE
            WRITE(3,100) NUMIT,CRIT
         ENDIF
         IF ((NUMIT .EQ. MAXIT) .AND. (CRIT .GT. EPSTA)) THEN
            WRITE(3,110) NUMIT,DELTA
            CALL DGMPNT(PISTA,NSTATE,1)
         ENDIF
   90 CONTINUE
      RETURN
  100 FORMAT(1X,'NUMIT = ',I5,'  CRIT = ',F7.6)
  110 FORMAT(1X,'*** UNSUCCESSFUL ITERATIONS FOR STATIONARY ',
     &          'DISTRIBUTION',1X,'NUMIT = ',I5,'  DELTA = ',F10.6)
      END 

     !




C.hr DSWEEP
C@
C....*...1.........2.........3.........4.........5.........6.........7.*
C     DSWEEP   10/4/72
C
C     PURPOSE
C     INVERT A POSITIVE DEFINITE MATRIX IN PLACE BY A DIAGONAL SWEEP.
C
C     USAGE
C     CALL DSWEEP(A,N,EPS,IER)
C
C     ARGUMENTS
C     A   - SYMMETRIC POSITIVE DEFINITE N BY N MATRIX STORED COLUMNWISE
C           (STORAGE MODE OF 0).  ON RETURN CONTAINS THE INVERSE OF A
C           STORED COLUMNWISE.
C           REAL*8
C     N   - NUMBER OF ROWS AND COLUMNS OF A.
C           INTEGER
C     EPS - INPUT CONSTANT USED AS A RELATIVE TOLERANCE IN TESTING FOR
C           DEGENERATE RANK.  A REASONABLE VALUE FOR EPS IS 1.D-13.
C           REAL*8
C     IER - ERROR PARAMETER CODED AS FOLLOWS
C           IER=0  NO ERROR, RANK OF A IS N.
C           IER.GT.0  A IS SINGULAR, RANK OF A IS N-IER.
C           INTEGER
C
C     REMARK
C     IF IER.GT.0 THEN DSWEEP RETURNS A G-INVERSE.
C
C     REFERENCE
C     SCHATZOFF, M. ET AL.  EFFICIENT CALCULATION OF ALL POSSIBLE REG-
C     RESSIONS.  TECHNOMETRICS, 10. 769-79 (NOVEMBER 1968)
C
C     PROGRAMMER
C     DR. A. RONALD GALLANT
C     DEPARTMANT OF STATISTICS
C     NORTH CAROLINA STATE UNIVERSITY
C     RALEIGH, NORTH CAROLINA  27695-8203
C
C
C
      SUBROUTINE DSWEEP(A,N,EPS,IER)
      IMPLICIT REAL*8 (A-H,O-Z)
      save
      REAL*8 A(N*N)
      TOL=0.D0
      DO 5 I=1,N
      II=-N+N*I+I
      TEST=A(II)
5     IF(TEST.GT.TOL) TOL=TEST
      TOL=TOL*EPS
      IER=0
      DO 50 K=1,N
      KK=-N+N*K+K
      AKK=A(KK)
      IF(AKK.GT.TOL) GO TO 20
      DO 10 J=K,N
      KJ=-N+N*J+K
10    A(KJ)=0.D0
      IF(K.EQ.1) GO TO 16
      KLESS1=K-1
      DO 15 I=1,KLESS1
      IK=KK-I
15    A(IK)=0.D0
16    IER=IER+1
      GO TO 50
20    D=1.D0/AKK
      DO 25 I=1,N
      DO 25 J=I,N
      IF((I.EQ.K).OR.(J.EQ.K)) GO TO 25
      IJ=N*(J-1)+I
      IF(I.LT.K) AIK=A(-N+N*K+I)
      IF(I.GT.K) AIK=A(-N+N*I+K)
      IF(K.LT.J) AKJ=A(-N+N*J+K)
      IF(K.GT.J) AKJ=-A(-N+N*K+J)
      A(IJ)=A(IJ)-AIK*AKJ*D
25    CONTINUE
      DO 30 J=K,N
      KJ=-N+N*J+K
30    A(KJ)=A(KJ)*D
      IF(K.EQ.1) GO TO 36
      KLESS1=K-1
      DO 35 I=1,KLESS1
      IK=KK-I
35    A(IK)=-A(IK)*D
36    A(KK)=D
50    CONTINUE
      DO 55 I=1,N
      DO 55 J=I,N
      IF(I.EQ.J) GO TO 55
      IJ=-N+N*J+I
      JI=-N+N*I+J
      A(JI)=A(IJ)
55    CONTINUE
      RETURN
      END
      
C           C....*...1.........2.........3.........4.........5.........6.........7.*.......8
C
C     This is a conversion of Algorithm 358 from COMPLEX to REAL*8.
C
C     Algorithm 358
C     Singular Value Decomposition of a Complex Matrix
C     Peter A. Businger and Gene H. Golub
C     Communications of the ACM, Volume 12, Number 10,
C     October, 1969, 564-565.
C
C     SVD finds the singular values S(1) >= S(2) >= ... >= S(N) of an
C     M by N matrix A with M >= N.  The computed singular values are
C     stored in the vector S. (B, C, T are work vectors of dimension N
C     or more). SVD also finds the first NU columns of an M by N
C     orthogonal matrix U and the first NV columns of an N by N
C     orthogonal matrix V such that ||A - U*D*V'|| is negligible
C     relative to ||A||, where D = diag[S(1), S(2), ..., S(N)].  (The
C     only values permitted for NU are 0, N, or M; those for NV are 0
C     or N.) Moreover, the transformation T(U) is applied to the P
C     vectors given in columns N_1, N+2, N+P of the array A.  This
C     feature can be used as follows to find the least squares solution
C     of minimal Euclidean length (the pseudoinverse solution) of an
C     over determined system Ax ~= b.  Call SVD with NV = N and with
C     columns N+1, N+2, ..., N+P of A containing P right hand sides b.
C     From the computed singular values determine the rank r of
C     diag(S) and define R = diag[1/S(1), ..., 1/S(r), 0, ..., 0].
C     Now x = VRd, where d = U'b is furnished by SVD in place of each
C     right-hand side b.
C
C     SVD can also be used to solve a homogeneous system of linear
C     equations.  To find an orthonormal basis for all solutions of
C     the system Ax = 0, call SVD with NV = N.  The desired basis
C     consists of those columns of V which correspond to negligible
C     singular values.  Further applications are mentioned in the
C     references.
C
C     The constants used in the program for ETA and TOL are machine
C     dependent.  ETA is the relative machine precision, TOL the
C     smallest normalized positive number divided by ETA.  The
C     assignments made are valid for a GE635 computer (a two's
C     complement binary machine with a signed 27-bit mantissa and a
C     signed 7-bit exponent). For this machine, ETA = 2**(-26) = 1.5E-8
C     and TOL = 2**(-129)/2**(-26) = 1.E-31.
C
C     References
C
C     Golub, G. Least squares singular values and matrix approximations.
C     Aplickace Matematiky 13 (1968), 44-51.
C
C     Golub, G., and Kahan, W. Calculating the singular value and
C     pseudoinverse of a matrix.  J. SIAM Numer. Anal. (1965), 205-224.
C
C     Golub, G., and Reinsch, C. Singular value decomposition and least
C     squares solutions.  Numer. Math. (to appear).

      SUBROUTINE DSVD(A,MMAX,NMAX,M,N,P,NU,NV,S,U,V,B,C,T)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      save
      REAL*8 A(MMAX,1), U(MMAX,1), V(NMAX,1), S(1), B(1), C(1), T(1)
      INTEGER M,N,P,NU,NV
      DATA ETA,TOL/0.222D-15,0.4D-289/
      NP=N+P
      N1=N+1
C
C     HOUSEHOLDER REDUCTION
      C(1)=0.D0
      K=1
10    K1=K+1
C
C     ELIMINATION OF A(I,K), I=K+1....,M
      Z=0.D0
      DO 20 I=K,M
20    Z=Z+A(I,K)**2
      B(K)=0.D0
      IF(Z.LE.TOL)GO TO 70
      Z=DSQRT(Z)
      B(K)=Z
      W=DABS(A(K,K))
      Q=1.D0
      IF(W.NE.0.D0)Q=A(K,K)/W
      A(K,K)=Q*(Z+W)
      IF(K.EQ.NP)GO TO 70
      DO 50 J=K1,NP
      Q=0.D0
      DO 30 I=K,M
30    Q=Q+A(I,K)*A(I,J)
      Q=Q/(Z*(Z+W))
      DO 40 I=K,M
40    A(I,J)=A(I,J)-Q*A(I,K)
50    CONTINUE
C
C     PHASE TRANSFORMATION
      Q=-A(K,K)/DABS(A(K,K))
      DO 60 J=K1,NP
60    A(K,J)=Q*A(K,J)
C
C     ELIMINATION OF A(K,J), J=K+Z,...,N
70    IF(K.EQ.N)GO TO 140
      Z=0.D0
      DO 80 J=K1,N
80    Z=Z+A(K,J)**2
      C(K1)=0.D0
      IF(Z.LE.TOL)GO TO 130
      Z=DSQRT(Z)
      C(K1)=Z
      W=DABS(A(K,K1))
      Q=1.D0
      IF(W.NE.0.D0)Q=A(K,K1)/W
      A(K,K1)=Q*(Z+W)
      DO 110 I=K1,M
      Q=0.D0
      DO 90 J=K1,N
90    Q=Q+A(K,J)*A(I,J)
      Q=Q/(Z*(Z+W))
      DO 100 J=K1,N
100   A(I,J)=A(I,J)-Q*A(K,J)
110   CONTINUE
C
C     PHASE TRANSFORMATION
      Q=-A(K,K1)/DABS(A(K,K1))
      DO 120 I=K1,M
120   A(I,K1)=A(I,K1)*Q
130   K=K1
      GO TO 10
C
C     TOLERANCE FOR NEGLIGIBLE ELEMENTS
140   EPS=0.D0
      DO 150 K=1,N
      S(K)=B(K)
      T(K)=C(K)
150   EPS=DMAX1(EPS,S(K)+T(K))
      EPS=EPS*ETA
C
C     INITIALIZATION OF U AND V
      IF(NU.EQ.0)GO TO 180
      DO 170 J=1,NU
      DO 160 I=1,M
160   U(I,J)=0.D0
170   U(J,J)=1.D0
180   IF(NV.EQ.0)GO TO 210
      DO 200 J=1,NV
      DO 190 I=1,N
190   V(I,J)=0.D0
200   V(J,J)=1.D0
C
C     QR DIAGONALIZATION
210   DO 380 KK=1,N
      K=N1-KK
C
C     TEST FOR SPLIT
220   DO 230 LL=1,K
      L=K+1-LL
      IF(DABS(T(L)).LE.EPS)GO TO 290
      IF(DABS(S(L-1)).LE.EPS)GO TO 240
230   CONTINUE
C
C     CANCELLATION OF E(L)
240   CS=0.D0
      SN=1.D0
      L1=L-1
      DO 280 I=L,K
      F=SN*T(I)
      T(I)=CS*T(I)
      IF(DABS(F).LE.EPS)GO TO 290
      H=S(I)
      W=DSQRT(F*F+H*H)
      S(I)=W
      CS=H/W
      SN=-F/W
      IF(NU.EQ.0)GO TO 260
      DO 250 J=1,N
      X=U(J,L1)
      Y=U(J,I)
      U(J,L1)=X*CS+Y*SN
250   U(J,I)=Y*CS-X*SN
260   IF(NP.EQ.N)GO TO 280
      DO 270 J=N1,NP
      Q=A(L1,J)
      R=A(I,J)
      A(L1,J)=Q*CS+R*SN
270   A(I,J)=R*CS-Q*SN
280   CONTINUE
C
C     TEST FOR CONVERGENCE
290   W=S(K)
      IF(L.EQ.K)GO TO 360
C
C     ORIGIN SHIFT
      X=S(L)
      Y=S(K-1)
      G=T(K-1)
      H=T(K)
      F=((Y-W)*(Y+W)+(G-H)*(G+H))/(2.D0*H*Y)
      G=DSQRT(F*F+1.D0)
      IF(F.LT.0.D0)G=-G
      F=((X-W)*(X+W)+(Y/(F+G)-H)*H)/X
C
C     QR STEP
      CS=1.D0
      SN=1.D0
      L1=L+1
      DO 350 I=L1,K
      G=T(I)
      Y=S(I)
      H=SN*G
      G=CS*G
      W=DSQRT(H*H+F*F)
      T(I-1)=W
      CS=F/W
      SN=H/W
      F=X*CS+G*SN
      G=G*CS-X*SN
      H=Y*SN
      Y=Y*CS
      IF(NV.EQ.0)GO TO 310
      DO 300 J=1,N
      X=V(J,I-1)
      W=V(J,I)
      V(J,I-1)=X*CS+W*SN
300   V(J,I)=W*CS-X*SN
310   W=DSQRT(H*H+F*F)
      S(I-1)=W
      CS=F/W
      SN=H/W
      F=CS*G+SN*Y
      X=CS*Y-SN*G
      IF(NU.EQ.0)GO TO 330
      DO 320 J=1,N
      Y=U(J,I-1)
      W=U(J,I)
      U(J,I-1)=Y*CS+W*SN
320   U(J,I)=W*CS-Y*SN
330   IF(N.EQ.NP)GO TO 350
      DO 340 J=N1,NP
      Q=A(I-1,J)
      R=A(I,J)
      A(I-1,J)=Q*CS+R*SN
340   A(I,J)=R*CS-Q*SN
350   CONTINUE
      T(L)=0.D0
      T(K)=F
      S(K)=X
      GO TO 220
C
C     CONVERGENCE
360   IF(W.GE.0.D0)GO TO 380
      S(K)=-W
      IF(NV.EQ.0)GO TO 380
      DO 370 J=1,N
370   V(J,K)=-V(J,K)
380   CONTINUE
C
C     SORT SINGULAR VALUES
      DO 450 K=1,N
      G=-1.D0
      J=K
      DO 390 I=K,N
      IF(S(I).LE.G)GO TO 390
      G=S(I)
      J=I
390   CONTINUE
      IF(J.EQ.K)GO TO 450
      S(J)=S(K)
      S(K)=G
      IF(NV.EQ.0)GO TO 410
      DO 400 I=1,N
      Q=V(I,J)
      V(I,J)=V(I,K)
400   V(I,K)=Q
410   IF(NU.EQ.0)GO TO 430
      DO 420 I=1,N
      Q=U(I,J)
      U(I,J)=U(I,K)
420   U(I,K)=Q
430   IF(N.EQ.NP)GO TO 450
      DO 440 I=N1,NP
      Q=A(J,I)
      A(J,I)=A(K,I)
440   A(K,I)=Q
450   CONTINUE
C
C     BACK TRANSFORMATION
      IF(NU.EQ.0)GO TO 510
      DO 500 KK=1,N
      K=N1-KK
      IF(B(K).EQ.0.D0)GO TO 500
      Q=-A(K,K)/DABS(A(K,K))
      DO 460 J=1,NU
460   U(K,J)=Q*U(K,J)
      DO 490 J=1,NU
      Q=0.D0
      DO 470 I=K,M
470   Q=Q+A(I,K)*U(I,J)
      Q=Q/(DABS(A(K,K))*B(K))
      DO 480 I=K,M
480   U(I,J)=U(I,J)-Q*A(I,K)
490   CONTINUE
500   CONTINUE
510   IF(NV.EQ.0)GO TO 570
      IF(N.LT.2)GO TO 570
      DO 560 KK=2,N
      K=N1-KK
      K1=K+1
      IF(C(K1).EQ.0.D0)GO TO 560
      Q=-A(K,K1)/DABS(A(K,K1))
      DO 520 J=1,NV
520   V(K1,J)=Q*V(K1,J)
      DO 550 J=1,NV
      Q=0.D0
      DO 530 I=K1,N
530   Q=Q+A(K,I)*V(I,J)
      Q=Q/(DABS(A(K,K1))*C(K1))
      DO 540 I=K1,N
540   V(I,J)=V(I,J)-Q*A(K,I)
550   CONTINUE
560   CONTINUE
570   RETURN
      END
      
C     C....*...1.........2.........3.........4.........5.........6.........7.*.......8
C     DMPINV 8/22/73
C
C     PURPOSE
C     TO COMPUTE THE MOORE-PENROSE INVERSE OF A MATRIX.
C
C     USAGE
C     CALL DMPINV(A,N,M,AMP,IR,W)
C
C     SUBROUTINES CALLED
C     DAPLUS, DSVD, DGMTRA
C
C     ARGUMENTS
C     A   - AN N BY M MATRIX STORED COLUMNWISE (STORAGE MODE 0)
C           ELEMENTS ARE REAL*8
C     N   - NUMBER OF ROWS IN A (NUMBER OF COLUMNS IN AMP)
C           INTEGER
C     M   - NUMBER OF COLUMNS IN A (NUMBER OF ROWS IN AMP)
C           INTEGER
C     AMP - MOORE-PENROSE INVERSE OF A
C           AN M BY N MATRIX STORED COLUMNWISE (STORAGE MODE 0)
C           ELEMENTS ARE REAL*8
C     IR  - COMPUTED RANK OF A
C           INTEGER
C     W   - A WORK VECTOR LENGTH N*M+4*MIN(N,M)+MIN(N,M)**2
C           ELEMENTS ARE REAL*8
C
C     PROGRAMMER
C     DR. THOMAS M. GERIG
C     DEPARTMENT OF STATISTICS
C     NORTH CAROLINA STATE UNIVERSITY
C     RALEIGH, NORTH CAROLINA  27609
C
C
      SUBROUTINE DMPINV(A,N,M,AMP,IR,W)
      IMPLICIT REAL*8(A-H,O-Z)
      save
      REAL*8 A(1),AMP(1),W(1)
      MIN=MIN0(N,M)
      I2=01+N*M
      I3=I2+MIN
      I4=I3+MIN**2
      I5=I4+MIN
      I6=I5+MIN
      IF(N.GE.M) GO TO 10
      CALL DGMTRA(A,W,N,M)
      CALL DAPLUS(W,M,N,A,AMP,W(I2),W(I3),IR,1.D-13,W(I4),W(I5),W(I6))
      CALL DGMTRA(A,AMP,N,M)
      CALL DGMTRA(W,A,M,N)
      RETURN
 10   CALL DAPLUS(A,N,M,AMP,W,W(I2),W(I3),IR,1.D-13,W(I4),W(I5),W(I6))
      RETURN
      END
      
C         C.hr DMFSD
C@
C....*...1.........2.........3.........4.........5.........6.........7.*
C     SUBROUTINE DMFSD
C
C     PURPOSE
C     FACTOR A GIVEN SYMMETRIC POSITIVE DEFINITE MATRIX
C
C     USAGE
C     CALL DMFSD(A,N,EPS,IER)
C
C     ARGUMENTS
C     A   - UPPER TRIANGULAR PART OF GIVEN SYMMETRIC POSITIVE DEFINITE N
C           BY N MATRIX; ON RETURN, UPPER TRIANGULAR R, A=R'R.
C           REAL*8
C     N   - THE NUMBER OF ROWS (COLUMNS) IN GIVEN MATRIX
C           INTEGER*4
C     EPS - INPUT CONSTANT WHICH IS USED AS RELATIVE TOLERANCE FOR TEST
C           ON LOSS OF SIGNIFICANCE, A REASONABLE VALUE IS 1.D-13.
C           REAL*8
C     IER - ERROR PARAMETER, CODED AS FOLLOWS
C           IER=0  - NO ERROR
C           IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARAMETER N OR
C                    BECAUSE SOME RADICAND IS NON-POSITIVE (MATRIX A IS
C                    NOT POSITIVE DEFINITE, POSSIBLY DUE TO LOSS OF
C                    SIGNIFICANCE)
C           IER=K  - WARNING WHICH INDICATES LOSS OF SIGNIFICANCE.  THE
C                    RADICAND FORMED AT FACTORIZATION STEP K+1 WAS STILL
C                    POSITIVE BUT NO LONGER GREATER THAN
C                    ABS(EPS*A(K+1,K+1)
C
C     REMARK
C     THIS IS A COPY OF SSP ROUTINE DMFSD
C
      SUBROUTINE DMFSD(A,N,EPS,IER)
      save
      REAL*8 EPS,DPIV,DSUM
      REAL*8 A((N*N+N)/2)
      IF(N-1) 12,1,1
    1 IER=0
      KPIV=0
      DO 11 K=1,N
      KPIV=KPIV+K
      IND=KPIV
      LEND=K-1
      TOL=ABS(SNGL(EPS*A(KPIV)))
      DO 11 I=K,N
      DSUM=0.D0
      IF(LEND) 2,4,2
    2 DO 3 L=1,LEND
      LANF=KPIV-L
      LIND=IND-L
    3 DSUM=DSUM+A(LANF)*A(LIND)
    4 DSUM=A(IND)-DSUM
      IF(I-K) 10,5,10
    5 IF(SNGL(DSUM)-TOL) 6,6,9
    6 IF(DSUM) 12,12,7
    7 IF(IER) 8,8,9
    8 IER=K-1
    9 DPIV=DSQRT(DSUM)
      A(KPIV)=DPIV
      DPIV=1.D0/DPIV
      GO TO 11
   10 A(IND)=DSUM*DPIV
   11 IND=IND+I
      RETURN
   12 IER=-1
      RETURN
      END
C                           C.hr DGMTRA
C@
C....*...1.........2.........3.........4.........5.........6.........7.*
C     DGMTRA  8/23/73
C
C     PURPOSE
C     TRANSPOSE A MATRIX: R = A'
C
C     USAGE
C     CALL DGMTRA(A,R,N,M)
C
C     ARGUMENTS
C     A - INPUT N BY M MATRIX STORED COLUMNWISE (STORAGE MODE 0)
C         REAL*8
C     R - OUTPUT M BY N MATRIX STORED COLUMNWISE (STORAGE MODE 0)
C         REAL*8
C     N - NUMBER OF ROWS IN A AND COLUMNS IN R
C         INTEGER*4
C     M - NUMBER OF COLUMNS IN A AND ROWS IN R
C         INTEGER*4
C
C
      SUBROUTINE DGMTRA(A,R,N,M)
      implicit real*8 (a-h,o-z)
      save
      REAL*8 A(N*M),R(M*N)
      DO 10 I=1,N
      DO 10 J=1,M
 10   R((I-1)*M+J)=A((J-1)*N+I)
      RETURN
      END
      
C      C.hr dgmprd
C@
*....*...1.........2.........3.........4.........5.........6.........7.*
*     dgmprd  6/15/83, 12/27/94
*
*     purpose
*     multiply two matrices: R = AB
*
*     usage
*     call dgmprd(A,B,R,n,m,l)
*
*     arguments
*     A - input n by m matrix stored columwise, i.e. A must either be
*         dimensioned as A(n,m) in the calling program or be a vector
*         that is indexed as A(-n+n*j+i) where i is the row index and
*         j is the column index.
*         input, real*8
*     B - input m by l matrix stored columnwise.
*         input, real*8
*     R - output n by l matrix stored columnwise.
*         output, real*8
*     n - number of rows in A and R.
*         input, integer*4
*     m - number of columns in A and number of rows in B.
*         input, integer*4
*     l - number of columns in B and R.
*         input, integer*4

      subroutine dgmprd(A,B,R,n,m,l)
      implicit none

*     arguments

      real*8 A,B,R
      integer*4 n,m,l

      dimension A(n,m),B(m,l),R(n,l)

*     local variables

      integer*4 i,j,k

*     clear

      do j=1,l
      do i=1,n
        R(i,j)=0.d0
      end do
      end do

*     accumulate

      do j=1,l
      do k=1,m
      do i=1,n
        R(i,j)=R(i,j)+A(i,k)*B(k,j)
      end do
      end do
      end do

      return
      end


C     C.hr DGMPNT
C@
C....*...1.........2.........3.........4.........5.........6.........7.*
C     DGMPNT   7/16/87
C
C     PURPOSE
C     PRINT A MATRIX.
C
C     USAGE
C     CALL DGMPNT(A,M,N)
C
C     ARGUMENTS
C     A - AN M BY N MATRIX STORED COLUMNWISE (STORAGE MODE 0).
C         REAL*8
C     M - NUMBER OF ROWS IN A
C         INTEGER*4
C     N - NUMBER OF COLUMNS IN A
C         INTEGER*4
C
C     COMMENT
C     THE DEFAULT LINESIZE IS 133, 132 PLUS CARRIAGE CONTROL CHARACTER.
C     THE USAGE:
C     COMMON /ZLNSIZ/ LNSIZE
C     LNSIZE=80
C     WILL CHANGE THE LINESIZE TO 80.  LINESIZES BETWEEN 72 AND 133 ARE
C     PERMITTED.
C
C     PROGRAMMER
C     DR. A. RONALD GALLANT
C     DEPARTMENT OF STATISTICS
C     NORTH CAROLINA STATE UNIVERSITY
C     RALEIGH, NORTH CAROLINA  27695-8203
C
C
      SUBROUTINE DGMPNT(A,M,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      save
      INTEGER*4 START,STOP,NOUT
      REAL*8 A(M,N),F
      CHARACTER*1 DIGIT(10),DUMMY(8),COL(9)
      CHARACTER*8 TYPE(14),FMT(13)
      CHARACTER*104 CFMT
      COMMON /ZLNSIZ/ LNSIZE
      EQUIVALENCE (FMT(1),DUMMY(1)),(FMT(1),CFMT)
      DATA FMT   /'(   X,  ',11*'        ','       )'/
      DATA DIGIT /'0','1','2','3','4','5','6','7','8','9'/
      DATA COL   /6*' ','C','O','L'/
      DATA TYPE  /',0PF12.1',',1PD12.4',',0PF12.8',',0PF12.7',
     &            ',0PF12.6',',0PF12.5',',0PF12.4',',0PF12.3',
     &            ',0PF12.2',',0PF12.1',',0PF12.0',',9A1,I3 ',
     &            ' 6X     ','''ROW'',I3'/
      NOUT=3
      LNSIZ=LNSIZE
      IF((LNSIZE.LT.72).OR.(LNSIZE.GT.133)) LNSIZ=133
      MAXCOL=(LNSIZ-8)/12
      IF(N.LT.MAXCOL) MAXCOL=N
      IPAD=(LNSIZ-8-12*MAXCOL)/2+1
      IPAD10=IPAD/10
      IPAD1=IPAD-10*IPAD10
      DUMMY(3)=DIGIT(IPAD10+1)
      DUMMY(4)=DIGIT(IPAD1+1)
      START=1
11    STOP=START-1+MAXCOL
      IF(STOP.GT.N) STOP=N
      K=2
      DO 13 J=START,STOP
      K=K+1
13    FMT(K)=TYPE(12)
      FMT(2)=TYPE(13)
      WRITE(NOUT,3001)
      WRITE(NOUT,3001)
      WRITE(NOUT,CFMT) (COL,J,J=START,STOP)
      WRITE(NOUT,3001)
      FMT(2)=TYPE(14)
      DO 19 I=1,M
      K=2
      DO 18 J=START,STOP
      K=K+1
      FMT(K)=TYPE(2)
      F=DABS(A(I,J))
      IF(F.LT.1.D+8 ) FMT(K)=TYPE(11)
      IF(F.LT.1.D+5 ) FMT(K)=TYPE(10)
      IF(F.LT.1.D+4 ) FMT(K)=TYPE( 9)
      IF(F.LT.1.D+3 ) FMT(K)=TYPE( 8)
      IF(F.LT.1.D+2 ) FMT(K)=TYPE( 7)
      IF(F.LT.1.D+1 ) FMT(K)=TYPE( 6)
C     IF(F.LT.1.D+0 ) FMT(K)=TYPE( 6)
      IF(F.LT.1.D-1 ) FMT(K)=TYPE( 5)
      IF(F.LT.1.D-2 ) FMT(K)=TYPE( 3)
      IF(F.LT.1.D-4 ) FMT(K)=TYPE( 2)
      IF(F.LT.1.D-38) FMT(K)=TYPE( 1)
18    CONTINUE
19    WRITE(NOUT,*) I,(A(I,J),J=START,STOP)
      IF(STOP.EQ.N) RETURN
      START=STOP+1
      GO TO 11
3001  FORMAT(' ')
      END
      
C              C....*...1.........2.........3.........4.........5.........6.........7.*.......8
C     DAPLUS   1/10/72
C
C     PURPOSE
C     COMPUTE THE MOORE-PENROSE G-INVERSE OF A MATRIX - A+
C     OBTAIN THE SINGULAR VALUE DECOMPOSITION OF A MATRIX
C
C     USAGE
C     CALL DAPLUS(A,N,M,AMP,U,S,V,IR,EPS,D1,D2,D3)
C
C     SUBROUTINES CALLED
C     DSVD
C
C     ARGUMENTS
C     A    - AN N BY M MATRIX STORED COLUMNWISE (STORAGE MODE OF 0)
C            ELEMENTS OF A ARE REAL*8
C     AMP  - MOORE-PENROSE G-INVERSE OF A
C            AN M BY N MATRIX STORED COLUMNWISE (STORAGE MODE OF 0)
C            ELEMENTS ARE REAL*8
C     N    - NUMBER OF ROWS IN A  (ON OUTPUT THE NO. OF COLS. OF A+)
C     M    - NUMBER OF COLUMNS IN A (ON OUTPUT THE NO. OF ROWS OF A+)
C            N MUST BE GREATER THAN OR EQUAL TO M
C     U    - AN N BY M MATRIX STORED COLUMNWISE (STORAGE MODE OF 0)
C            ELEMENTS OF U ARE REAL*8
C     S    - A DIAGONAL MATRIX STORED AS AN M-VECTOR (STORAGE MODE OF 2)
C            ELEMENTS OF S ARE REAL*8
C     V    - AN M BY M MATRIX STORED COLUMNWISE (STORAGE MODE OF 0)
C            ELEMENTS OF V ARE REAL*8
C     IR   - COMPUTED RANK OF A
C            INTEGER
C     EPS  - VALUE USED TO TEST THE SINGULAR VALUES OF A AND DETERMINE
C            THE RANK OF A.  A REASONABLE VALUE IS EPS = 1.D-13
C            REAL*8
C     D1   - AN M VECTOR USED FOR WORKSPACE.
C     D2   - AN M VECTOR USED FOR WORKSPACE.
C     D3   - AN M VECTOR USED FOR WORKSPACE.
C            ELEMENTS OF D1,D2,D3 ARE REAL*8
C
C     REMARKS
C     A = U*S*(V-TRANSPOSE)
C     (U-TRANSPOSE)*U = (V-TRANSPOSE)*V = V*(V-TRANSPOSE) = I
C     A+ = V*(S+)*(U-TRANSPOSE)
C     S+(I) = 1/S(I) IF S(I).GT.0 AND S+(I)=0 OTHERWISE
C     IF A CAN BE DESTROYED THE FOLLOWING USAGE WILL CONSERVE CORE
C     CALL DAPLUS(A,N,M,AMP,A,S,V,IR,EPS,D1,D2,D3)
C     THE MATRIX A WILL CONTAIN U ON RETURN
C
C     REFERENCE
C     BUSINGER,P.A. AND GOLUB,G.H., SINGULAR VALUE DECOMPOSITION OF A
C     COMPLEX MATRIX. COMMUNICATIONS OF THE ACM 12: 564-565 (OCT. 1969)
C
      SUBROUTINE DAPLUS(A,N,M,AMP,U,S,V,IR,EPS,D1,D2,D3)
      IMPLICIT REAL*8 (A-H,O-Z)
      save
      REAL*8 A(N,M),AMP(M,N),U(N,M),S(M),V(M,M),D1(M),D2(M),D3(M)
      CALL Z1DAPLUS(A,N,M,AMP)
      CALL DSVD(AMP,N,M,N,M,0,M,M,S,U,V,D1,D2,D3)
C
C     RANK DETERMINATION
      IR=0
      TEST=S(1)*EPS
      DO 10 J=1,M
      IF(S(J).GT.TEST) IR=IR+1
      D1(J)=0.D0
10    IF(S(J).GT.TEST) D1(J)=1.D0/S(J)
C
      CALL Z2DAPLUS(AMP,M,N,V,D1,U,D2)
      RETURN
      END
      SUBROUTINE Z1DAPLUS(A,N,M,AMP)
      implicit real*8 (a-h,o-z)
      save
      REAL*8 A(1),AMP(1)
      L=M*N
      DO 10 I=1,L
10    AMP(I)=A(I)
      RETURN
      END
      SUBROUTINE Z2DAPLUS(A,M,N,V,D1,U,D2)
      implicit real*8 (a-h,o-z)
      save
      REAL*8 A(M,N), V(M,M),D1(M), U(N,M), D2(M)
      DO 30 I=1,M
      DO 10 K=1,M
10    D2(K)=V(I,K)*D1(K)
      DO 20 J=1,N
      A(I,J)=0.D0
      DO 20 K=1,M
20    A(I,J)=A(I,J)+D2(K)*U(J,K)
30    CONTINUE
      RETURN
      END
      
C      C.hr HEADER
C@
C....*...1.........2.........3.........4.........5.........6.........7.*
C     HEADER   7/16/87
C
C     PURPOSE
C     TITLE OUTPUT.
C
C     USAGE
C     CALL HEADER(TITLE)
C
C     ARGUMENTS
C     TITLE - A CHARACTER STRING ENCLOSED IN QUOTES.  THE LAST CHARACTER
C             MUST BE AN UNDERSCORE.  THE STRING IS MADE UP OF AN UN-
C             LIMITED NUMBER OF CONCATONATED LINES.  A LINE IS MADE UP
C             OF 0-68 CHARACTERS FOLLOWED BY A SLASH.
C
C     REMARK
C     AN EXAMPLE IS:
C     CALL HEADER('//LINE 1/LINE 2/LINE 3///_')
C
C     COMMENT
C     THE DEFAULT LINESIZE IS 133, 132 PLUS CARRIAGE CONTROL CHARACTER.
C     THE USAGE:
C     COMMON /ZLNSIZ/ LNSIZE
C     LNSIZE=80
C     WILL CHANGE THE LINESIZE TO 80.  LINESIZES BETWEEN 72 AND 133 ARE
C     PERMITTED.
C
C
      SUBROUTINE HEADER(TITLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      save
      COMMON /ZLNSIZ/ LNSIZE
      CHARACTER*1 DIGIT(10),LDUMMY(8)
      CHARACTER*1 TITLE(133),PRINT(68),SLASH,UNC,LTR,BLANK
      CHARACTER*4 DUMMY(18)
      CHARACTER*8 FMT1(11),FMT2(3),FMT3(3),RDUMMY
      CHARACTER*24 CFMT2,CFMT3
      CHARACTER*88 CFMT1
      EQUIVALENCE (FMT1(1),CFMT1),(FMT2(1),CFMT2),(FMT3,CFMT3)
      EQUIVALENCE (RDUMMY,LDUMMY(1)),(DUMMY(1),PRINT(1))
      DATA FMT1/'(1X,    ','X,70H***',8*'********','***  )  '/
      DATA FMT2/'(1X,    ','X,1H*,17','A4,1H*) '/
      DATA FMT3/'(1X,    ','X,1H*,68','X,1H*)  '/
      DATA DIGIT/'0','1','2','3','4','5','6','7','8','9'/
      DATA SLASH/'/'/,UNC/'_'/,BLANK/' '/
      NOUT=3
      LNSIZ=LNSIZE
      IF((LNSIZE.LT.72).OR.(LNSIZE.GT.133)) LNSIZ=133
      RDUMMY=FMT1(1)
      IPAD=(LNSIZ-72)/2+1
      IPAD10=IPAD/10
      IPAD1=IPAD-10*IPAD10
      LDUMMY(7)=DIGIT(IPAD10+1)
      LDUMMY(8)=DIGIT(IPAD1+1)
      FMT1(1)=RDUMMY
      FMT2(1)=RDUMMY
      FMT3(1)=RDUMMY
      WRITE(NOUT,3001)
      WRITE(NOUT,3001)
      WRITE(NOUT,CFMT1)
      M=0
      L=0
10    M=M+L
      L=0
      DO 20 I=1,69
      LTR=TITLE(M+I)
      IF((LTR.EQ.SLASH).OR.(LTR.EQ.UNC)) GO TO 25
20    L=I
25    IF((L.EQ.0).AND.(LTR.EQ.UNC)) GO TO 50
      IF(L.EQ.0) GO TO 60
      IPAD=(68-L)/2
      DO 30 I=1,68
30    PRINT(I)=BLANK
      DO 40 I=1,L
40    PRINT(I+IPAD)=TITLE(M+I)
      WRITE(NOUT,CFMT2) (DUMMY(I),I=1,17)
      L=L+1
      LTR=TITLE(M+L+1)
      IF(LTR.EQ.UNC) GO TO 50
      GO TO 10
50    WRITE(NOUT,CFMT1)
      RETURN
60    WRITE(NOUT,CFMT3)
      L=L+1
      GO TO 10
3001  FORMAT(' ')
      END


C     C.hr SPACER
C@
C....*...1.........2.........3.........4.........5.........6.........7.*
C     SPACER   5/10/73
C
C     PURPOSE
C     SPACE OUTPUT
C
C     USAGE
C     CALL SPACER(I)
C
C     ARGUMENT
C     I - NUMBER OF LINES TO BE SKIPPED.  IF I=0 THEN THE PRINTER IS
C         ADVACED TO A NEW PAGE.
C         INTEGER*4
C
C
      SUBROUTINE SPACER(I)
      implicit real*8 (a-h,o-z)
      save
      nout=3
      IF(I.EQ.0) WRITE(nout,301)
      IF(I/66.GT.0) WRITE(nout,301)
      IF(I.LE.0) RETURN
      M=MOD(I,66)
      DO 10 J=1,M
10    WRITE(nout,302)
      RETURN
301   FORMAT('1')
302   FORMAT(' ')
      END
