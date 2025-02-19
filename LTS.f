      PROGRAM LTS
c*******************  LTS.FOR  ************************************     
c---------------------------------------------------------------------- 
c PROGRAM   : L (ASER) T (RANSMISSION) S (PECTROSCOPY)                 
c                                                                      		
c
c             #12 LTSSO : CALCULATED SOLUTION                          	             
c             #13 LTSADPA : COEFFICIENTS AJ                             
c                                                                       
c VARIABLES : NMAX   : MAXIMUM NUMBER OF DATA POINTS                    
c             MMAX   : MAXIMUM NUMBER OF COEFFICIENTS AJ                
c             NSMAX  : MAXIMUM NUMBER OF POINTS WHERE THE SOL. IS CALCULATED  
c                                                                       
c             N      : NUMBER OF DATA POINTS                            
c             NE     : NUMBER OF DIFFERENT EQUATIONS                    
c             TI     : ABSZISSE OF THE DATA POINTS                      
c             GSTI   : ORDINATE OF THE DATA POINTS                      
c             GERRTI : ERROR OF THE DATA POINTS                         
c             ERROR  : SCALEFACTOR OF THE ERROR                         
c                                                                       
c             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
c             SMIN,                                                     
c             SMAX   : INTERVAL IN WHICH THE SOLUTION IS CALCULATED     
c             SJ     : ABSZISSES WHERE THE SOL. IS CALCULATED           
c             FSSJ   : VALUES OF THE CALCULATED SOLUTION                
c             FMSJ   : MIDPOINT FOR THE ERRORS OF THE CALCULATED SOL.   
c             FERRSJ : ERRORS OF THE CALCULATED SOLUTION                
c                                                                       
c             M      : NUMBER OF COEFFICIENTS AJ                        
c             ASJ    : VALUES OF THE CALCULATED COEFFICIENTS AJ         
c             AERRJ  : ERRORS OF THE CALCULATED COEFFICIENTS AJ         
c                                                                       
c             B      : ADDITIONAL LINEAR TERMS                          
c             K      : APPROXIMATED INTEGRAL OPERATOR                   
c                                                                       
c             NL     : FIRST DIMENSION OF THE REG. OPERATOR             
c             L      : APPROXIMATED REGULARIZATION OPERATOR             
c                                                                       
c             DISMOD : MODE FOR THE DISCRETIZATION OF THE INTEGRALS     
c             ERRMOD : MODE FOR THE ESTIMATION OF THE ERROR             
c             REGMOD : MODE FOR THE CALCULATION OF THE SOLUTION         
c             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
c                                                                       
c             LAMBDA : REGULARIZATION PARAMETER                         
c             LAMB?? : PARAMETER FOR THE CALCULATION OF LAMBDA          
c                                                                       
c             IFWK   : INTEGER WORKSPACE                                
c             DPFWK  : DOUBLE PRECISION WORKSPACE                       
c                                                                       
c COMMENT   : THE SECOND DIMENSION OF K MUST BE AT LEAST (M + NS)       
c                                                                       
c             THE SIZE OF THE INTEGER WORKSPACE MUST BE AT LEAST        
c             MAX(2(M + NS) + NEFF + NL,N)                              
c                                                                       
c             THE SIZE OF THE DOUBLE PRECISION WORKSPACE MUST BE AT     
c             LEAST (M + NS) (M + NS + NEFF + NL + 4) + NEFF + NL       
c                                                                       
c             (NEFF = MIN(N,M + NS))                                    
c---------------------------------------------------------------------- 
      IMPLICIT NONE 
                                                                  
      INTEGER NMAX,MMAX,NSMAX                                           
      PARAMETER(NMAX = 7000,MMAX = 5,NSMAX = 7000)                        
                                                                        
      INTEGER N,NE,J                                                      
      DOUBLE PRECISION ERROR                                            
      DOUBLE PRECISION TI(NMAX),GSTI(NMAX),GERRTI(NMAX)                 
      
      COMPLEX RI             
                                                           
      INTEGER NS                                                        
      DOUBLE PRECISION SMIN,SMAX                                        
      DOUBLE PRECISION SJ(NSMAX),FSSJ(NSMAX),FMSJ(NSMAX),FERRSJ(NSMAX)  
                                                                        
      INTEGER M                                                         
      DOUBLE PRECISION ASJ(MMAX),AERRJ(MMAX)                            
                                                                        
      DOUBLE PRECISION B(NMAX,MMAX)                                     
      DOUBLE PRECISION K(NMAX,MMAX + NSMAX)                             
                                                                        
      INTEGER NL                                                        
      DOUBLE PRECISION L(NSMAX + 2,NSMAX)                               
                                                                        
      INTEGER DISMOD,ERRMOD,REGMOD,INFMOD                               
                                                                        
      INTEGER LAMBIT                                                    
      DOUBLE PRECISION LAMBDA,LAMBST,LAMBSP,LAMBRA,LAMBPR               
                                                                        
      INTEGER IFWK(NMAX + 3 * (MMAX + NSMAX) + NSMAX + 2)               
      DOUBLE PRECISION DPFWK((2 * (MMAX + NSMAX) + NSMAX + 6) *         
     *                       (MMAX + NSMAX) + MMAX + 2 * NSMAX + 2)     
      
                                                                        
c-----START PROGRAM---------------------------------------------------- 
                                                                        
      WRITE(*,*)'--------------------------------------------------'    
      WRITE(*,*)'-----------------LTS STARTED----------------------'    
                                                                        
c-----READ INPUT PARAMETER AND DATA------------------------------------ 
                                                                        
      CALL RDPARA(NS,NSMAX,SMIN,SMAX,M,MMAX,N,NMAX,NE,ERROR,LAMBDA,     
     *            DISMOD,ERRMOD,REGMOD,INFMOD,LAMBST,LAMBSP,LAMBRA,     
     *            LAMBPR,LAMBIT)                                        
      CALL RDDATA(N,NE,TI,GSTI,GERRTI,ERRMOD,INFMOD)                    
                                                                        
c-----SETUP DISCRETE PROBLEM------------------------------------------- 
                                                                        
      CALL SEDIPR(N,NE,TI,GSTI,GERRTI,NS,SMIN,SMAX,SJ,M,B,NMAX,K,NMAX,  
     *            NL,L,NSMAX + 2,DISMOD,REGMOD,INFMOD)                  
                                                                        
c-----CALCULATE REGULARIZED SOLUTION----------------------------------- 
                                                                        
      IF (MOD(REGMOD / 10,10).EQ.1) THEN                                
       CALL CARSFR(N,GSTI,M,ASJ,AERRJ,NS,FSSJ,FMSJ,FERRSJ,B,NMAX,K,     
     *             NMAX,NL,L,NSMAX + 2,ERROR,LAMBDA,REGMOD / 100,       
     *             INFMOD,IFWK,DPFWK)                                   
      ELSE                                                              
       CALL CARSSC(N,GSTI,M,ASJ,AERRJ,NS,FSSJ,FMSJ,FERRSJ,B,NMAX,K,     
     *             NMAX,NL,L,NSMAX + 2,ERROR,LAMBDA,ERRMOD / 10,        
     *             REGMOD / 100,INFMOD,LAMBST,LAMBSP,LAMBRA,LAMBPR,     
     *             LAMBIT,IFWK,DPFWK)                                   
      ENDIF                                                             
                                                                        
c-----WRITE OUTPUT DATA------------------------------------------------ 
                                                                        
      CALL WERESO(NS,SJ,FSSJ,FMSJ,FERRSJ)                               
      CALL WECOAJ(M,ASJ,AERRJ)                                          
                                                                        
      WRITE(*,1001)'MAIN   > REGULARIZATION PARAMETER    :',LAMBDA      
      WRITE(*,1001)'MAIN   > ERROR                       :',ERROR       
      WRITE(*,*)'------------------LTS ENDED-------------------'    
      WRITE(*,*)'--------------------------------------------------'    
                                                                        
      STOP                                                              
                                                                        
1001  FORMAT(1X,A38,D18.8)                                              
                                                                        
      END                                                               
                                                                      
                                                                        
c---------------------------------------------------------------------- 
c FUNCTION  : KTISJ1                                                    
c                                                                       
c PURPOSE   : CALCULATE VALUE OF KERNELFUNCTION 1                       
c                                                                       
c VARIABLES : (ON INPUT)                                                
c             T,SFORW    : ARGUMENTS OF THE KERNELFUNCTION 1                
c                                                                       
c             (ON OUTPUT)                                               
c             KTISJ1 : CALCULATED VALUE OF KERNELFUNCTION 1             
c---------------------------------------------------------------------- 
                                                                        
      DOUBLE PRECISION FUNCTION KTISJ1(T,RESFORW)                             
                                                                        
      DOUBLE PRECISION T,RESFORW
                                                   
      DOUBLE PRECISION pi
      PARAMETER( pi=3.141593 )                                                                   
                                                   
                                                                        
      KTISJ1 =((T**2)/pi)*RESFORW                                
                                                                        
      RETURN                                                            
      END                                                               
      
c      DOUBLE PRECISION FUNCTION KTISJ1(T,S)                             
c                                                                        
c      DOUBLE PRECISION T,S                                              
c                                                                        
c      DOUBLE PRECISION DPWK                                             
c                                                                        
c      KTISJ1 = DEXP(-T / S)                                             
c                                                                        
c      RETURN                                                            
c      END 	                                                                  
                                                                        
c---------------------------------------------------------------------- 
c FUNCTION  : KTISJ2                                                    
c                                                                       
c PURPOSE   : CALCULATE VALUE OF KERNELFUNCTION 2                       
c                                                                       
c VARIABLES : (ON INPUT)                                                
c             T,S    : ARGUMENTS OF THE KERNELFUNCTION 2                
c                                                                       
c             (ON OUTPUT)                                               
c             KTISJ2 : CALCULATED VALUE OF KERNELFUNCTION 2             
c---------------------------------------------------------------------- 
                                                                        
       DOUBLE PRECISION FUNCTION KTISJ2(T,S)                             
                                                                        
      DOUBLE PRECISION T,S                                              
                                                                        
      DOUBLE PRECISION DPWK                                             
                                                                        
      DPWK = T * S                                                      
      KTISJ2 = DPWK / (1.D0 + DPWK * DPWK)                              
                                                                       
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
c---------------------------------------------------------------------- 
c FUNCTION  : BTIJ1                                                     
c                                                                       
c PURPOSE   : CALCULATE VALUE OF THE ADDITIONAL LINEAR TERM 1,J         
c                                                                       
c VARIABLES : (ON INPUT)                                                
c             T      : ARGUMENT OF THE ADDITIONAL LINEAR TERM           
c             J      : NUMBER OF THE ADDITIONAL LINEAR TERM             
c                                                                       
c             (ON OUTPUT)                                               
c             BTIJ1  : CALCULATED VALUE OF THE ADD. LINEAR TERM 1,J     
c---------------------------------------------------------------------- 
                                                                        
      DOUBLE PRECISION FUNCTION BTIJ1(T,J)                              
                                                                        
      DOUBLE PRECISION T                                                
                                                                        
      INTEGER J                                                         
                                                                        
      CHARACTER*60 ERRTXT                                               
                                                                        
      IF (J.EQ.1) THEN                                                 
       BTIJ1 = 0.D0                                                    
      ELSE IF (J.EQ.2) THEN                                            
       BTIJ1 = 0.D0                                                    
      ELSE                                                             
       ERRTXT = 'FUNCTION NOT DEFINED'                       
       CALL WEERRM(ERRTXT,.TRUE.)                                      
      ENDIF
                                                            
c Nel nostro caso il termine lineare è nullo       

c      BTIJ1 = 0.D0                                                      
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
c---------------------------------------------------------------------- 
c FUNCTION  : BTIJ2                                                     
c                                                                       
c PURPOSE   : CALCULATE VALUE OF THE ADDITIONAL LINEAR TERM 2,J         
c                                                                       
c VARIABLES : (ON INPUT)                                                
c             T      : ARGUMENT OF THE ADDITIONAL LINEAR TERM           
c             J      : NUMBER OF THE ADDITIONAL LINEAR TERM             
c                                                                       
c             (ON OUTPUT)                                               
c             BTIJ2  : CALCULATED VALUE OF THE ADD. LINEAR TERM 2,J     
c---------------------------------------------------------------------- 
                                                                        
      DOUBLE PRECISION FUNCTION BTIJ2(T,J)                              
                                                                        
      DOUBLE PRECISION T                                                
                                                                        
      INTEGER J                                                         
                                                                        
      CHARACTER*60 ERRTXT                                               
                                                                        
      IF (J.EQ.1) THEN                                                 
       BTIJ2 = 0.D0                                                    
      ELSE IF (J.EQ.2) THEN                                            
       BTIJ2 = 1.D0 / T                                                
      ELSE                                                             
       ERRTXT = ' BTIJ2  > FUNCTION NOT DEFINED'                       
       CALL WEERRM(ERRTXT,.TRUE.)                                      
      ENDIF
                                                            
c Nel nostro caso il termine lineare è nullo
        
C      BTIJ2 = 0.D0                                                      
      RETURN                                                            
      END                                                               


c----------------------------------------------------------------------
c SUBROUTINE: R(ERFACTIVE) I(NDEX) MA(TERIAL)
c----------------------------------------------------------------------
	 
        SUBROUTINE RIMA(TI,R,C)
	
        DOUBLE PRECISION TI,U,U1
        DOUBLE PRECISION R,C,R1,C1,R2,C2
c--------NANOMETER---->MICROMETER--------------------------------------
        
        U=TI/1000.d0

C-------PARTICLE MATERIALS---------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------


c--------POLYSTYRENE MA---------------------------------------------------
c        U1=U/1.5601

	R1=1.5725+0.0031080/U**2+0.00034779/U**4

c--------SULTANOVA
c	R1=DSQRT(1+(1.4435*U**2)/(U**2-0.020216))

c        R1= 1.6184-0.0017*U+3.2756E-5*U**2-4.1483E-7*U**3+2.9218E-9*U**4-8.5139E-12*U**5
	C1=0.0000001	
c---------------------------------------------------------------------- 
C--------SILVER--------------------------------------------------------
        R2=16.5934-152.0785*U+565.4346*U**2-1089.4660*U**3
     *     +1149.1601*U**4-630.4295*U**5 + 140.8213*U**6 
     
        C2=-55.6159+489.1477*U-1741.9084*U**2+3315.5199*U**3
     *    -3489.4885*U**4 + 1922.4835*U**5-432.9903*U**6
C----------------------------------------------------------------------
C-------SILICA---------------------------------------------------------
       R3=DSQRT(1+0.6961663/(1-(0.0684043/U)**2)
     *   + 0.4079426/(1-(0.1162414/U)**2)+0.8974794/(1-(9.896161/U)**2))
       C3=0.0000001


C-------RI-------------------------------------------------------------
 
	R=R1
	C=C1

	RETURN
	END
	
c----------------------------------------------------------------------
c SUBROUTINE: R(ERFACTIVE) I(NDEX) ME(DIUM)
c----------------------------------------------------------------------
	 
        SUBROUTINE RIME(TI,RE,CO)
	
        DOUBLE PRECISION TI,U
	DOUBLE PRECISION RE,CO
c--------NANOMETER---->MICROMETER--------------------------------------
        
        U=TI/1000.d0

C-------DISPERSANT MEDIUM----------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------

C-------WATER----------------------------------------------------------
c       RE=1.8412-4.7267*U+18.8292*U**2-41.5388*U**3+55.6176*U**4
c     *  - 46.9324*U**5+25.0486*U**6-8.1905*U**7+1.4958*U**8-0.1168*U**9
       RE=1.4430*exp(-U/0.0611) + 1.343 -0.017*U
       CO=1e-8



	RETURN
	END               
              
c----------------------------------------------------------------------  
c SUBROUTINE: SE(TUP) DI(SCRETE) PR(OBLEM)                              
c                                                                       
c PURPOSE   : SETUP VALUES FOR THE DISCRETE REGULARIZATION PROBLEM AS   
c             NEEDED BY THE SUBROUTINE CARSSC AND CARSFR                
c                                                                       
c VARIABLES : (ON INPUT)                                                
c             N      : NUMBER OF DATA POINTS                            
c             NE     : NUMBER OF DIFFERENT EQUATIONS                    
c             TI     : ABSZISSE OF THE DATA POINTS                      
c             GSTI   : ORDINATE OF THE DATA POINTS                     
c             GERRTI : ERROR OF THE DATA POINTS                         
c                                                                       
c             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
c             SMIN,                                                     
c             SMAX   : INTERVAL IN WHICH THE SOLUTION IS CALCULATED     
c                                                                       
c             M      : NUMBER OF COEFFICIENTS AJ                        
c             NBMAX  : LEADING DIMENSION OF THE FIELD B                 
c             NKMAX  : LEADING DIMENSION OF THE FIELD K                 
c             NLMAX  : LEADING DIMENSION OF THE FIELD L                 
c                                                                      
c             DISMOD : MODE FOR THE DISCRETIZATION OF THE INTEGRALS     
c             REGMOD : MODE FOR THE CALCULATION OF THE SOLUTION         
c              INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION       
c                                                                       
c           (ON OUTPUT)                                                 
c             GSTI   : CHANGED DATA POINTS FOR INPUT IN CARSSC AND      
c                      CARSFR                                           
c                                                                       
c             SJ     : ABSZISSES WHERE THE SOL. IS CALCULATED           
c                                                                       
c             B      : ADDITIONAL LINEAR TERMS                          
c                                                                       
c             K      : APPROXIMATED INTEGRAL OPERATOR                   
c                                                                       
c             NL     : FIRST DIMENSION OF THE REG. OPERATOR             
c             L      : APPROXIMATED REGULARIZATION OPERATOR             
c---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE SEDIPR(N,NE,TI,GSTI,GERRTI,NS,SMIN,SMAX,SJ,M,B,NBMAX,  
     *                  K,NKMAX,NL,L,NLMAX,DISMOD,REGMOD,INFMOD)        
                                                                      
      INTEGER N,NE                                                      
      DOUBLE PRECISION TI(N),GSTI(N),GERRTI(N)                          
                                                                        
      INTEGER NS                                                        
      DOUBLE PRECISION SMIN,SMAX                                        
      DOUBLE PRECISION SJ(NS)                                           
                                                                        
      INTEGER M                                                         
                                                                       
      INTEGER NBMAX                                                     
      DOUBLE PRECISION B(NBMAX,M)                                       
                                                                        
      INTEGER NKMAX                                                     
      DOUBLE PRECISION K(NKMAX,NS)                                      
                                                                        
      INTEGER NL,NLMAX                                                  
      DOUBLE PRECISION L(NLMAX,NS)                                      
                                                                                                                                                                           
      DOUBLE PRECISION KTISJ1,KTISJ2,BTIJ1,BTIJ2    

      INTEGER I,J                     
    
      COMPLEX RI,REFIN

      INTEGER DISMOD,REGMOD,INFMOD
                                                                                                                                                            
      DOUBLE PRECISION HS,T,X,U                                             
                                                                        
      INTEGER IWK
      DOUBLE PRECISION DPWK  
      
      DOUBLE PRECISION pi
      PARAMETER (pi=3.141593)
      
      DOUBLE PRECISION RESFORW
c ----------------------------------------------------------------------
c --------  I / O SPECIFICATIONS FOR SUBROUTINE MIEV0  -----------------
c ----------------------------------------------------------------------
      INTEGER  MAXANG, MOMDIM
      PARAMETER  ( MAXANG = 40, MOMDIM = 10 )
      LOGICAL  ANYANG, PERFCT, PRNT(2)
      INTEGER  IPOLZN, NUMANG, NMOM
      PARAMETER  ( IPOLZN=14 ) 
      
      REAL*8     GQSC, MIMCUT, PMOM( 0:MOMDIM, 4 ), QEXT, QSCA, SPIKE
      REAL*8     XMU(MAXANG), XX, REREFIN, IMAREFIN, R,C,RE,CO
      COMPLEX    CREFIN, SFORW, SBACK, S1( MAXANG ), S2( MAXANG )
      COMPLEX    TFORW( 2 ), TBACK( 2 )

                                          
                                                                        
c-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' SEDIPR > CALLED'                                     
      ENDIF                                                             
                                                                        
c-----CALCULATED CHANGED DATA POINTS FOR INPUT IN----------------------
c-----CARSSC AND CARSFR------------------------------------------------ 
                                                                        
      DO 10 I = 1,N                                                     
       GSTI(I) = GSTI(I) / GERRTI(I)                                    
10    CONTINUE                                                          
                                                                        
c-----CALCULATE ABSIZZES WHERE THE SOLUTION IS CALCULATED--------------  
                                                                 
      IF (DISMOD.EQ.1) THEN                                             
       HS = (SMAX - SMIN) / (NS - 1)                                    
       SJ(1) = SMIN                                                  
       DO 20 I = 2,NS                                                   
        SJ(I) = SJ(I - 1) + HS                                          
20     CONTINUE                                                         
      ELSE                                                              
       HS = (SMAX / SMIN)**(1.D0 / DBLE(NS - 1))                        
       SJ(1) = SMIN                                                     
       DO 30 I = 2,NS                                                   
        SJ(I) = SJ(I - 1) * HS                                          
30     CONTINUE                                                         
       HS = DLOG(HS)                                                    
      ENDIF 

                                                            
	                                                                       
c-----SETUP MATRIX FOR THE APPROXIMATED INTEGRAL OPERATOR-------------- 
	
C      WRITE(*,*) 'MIMCUT=MIN IMMAGINARY REFRACTIV INDEX (REAL POSITIVE)'	
C      READ(*,*) MIMCUT
        
        MIMCUT=0.000000000000001
      	PERFCT=.FALSE.
	ANYANG=.FALSE.
	NUMANG=0
	XMU(NUMANG)=1.0
	NMOM=0
	
      
        
      DO 40 I = 1,(N / NE) 
c----------- REFRACTIVE INDEX -----------------------------------------
	CALL RIMA(TI(I),R,C)
	CALL RIME(TI(I),RE,CO)

         REREFIN= R/RE
         IMAREFIN=C

         CREFIN=CMPLX(REREFIN,IMAREFIN) 
                       
       	DO 41 J = 1,NS  
                                                                                
        WRITE(16,*)HS      
    
c----------- MIE PARAMETER --------------------------------------------
         XX=(2*pi*SJ(J)*RE)/TI(I)
c         XX=(2*pi*SJ(J))/TI(I)


c----------- MIEV0 ----------------------------------------------------
      
         CALL MIEV0(XX,CREFIN,PERFCT,MIMCUT,ANYANG,NUMANG,XMU,NMOM,
     *              IPOLZN,MOMDIM,PRNT,QEXT,QSCA,GQSC,PMOM,SFORW,
     *              SBACK,S1,S2,TFORW,TBACK,SPIKE)   
         
	 WRITE(18,*) SFORW
C	 RESFORW=REAL(SFORW)                               
C         K(I,J) = KTISJ1(TI(I),RESFORW) * HS /  GERRTI(I)  
          K(I,J) = (pi*SJ(J)**2)*QEXT * HS / GERRTI(I)
C         WRITE(17,1001) KTISJ1(TI(I),RESFORW)/(pi*(SJ(J)/2)**2)    
41    CONTINUE  
      
      WRITE(19,1002)(K(I,J)/(pi*SJ(J)**2),J=1,NS)

40    CONTINUE

C------------------------------------------------------------
C---------------ATTENTION------------------------------------
C-----FORMAT MUST BE UPDATED WHEN CHANGING NS----------------
1002  FORMAT(2X,690D18.8) 

                                                         
      IF (NE.EQ.2) THEN                                                
       DO 50 I = 1,(N / 2)                                            
        IWK = I + (N / 2)                                              
        DO 50 J = 1,NS                                                 
         K(IWK,J) = KTISJ2(TI(I),SJ(J)) * HS / GERRTI(IWK)             
50      CONTINUE                                                       
      ENDIF                                                             
                                                                        
C-----SETUP MATRIX FOR THE ADDITIONAL LINEAR TERMS--------------------- 
                                                                        
      IF (M.GT.0) THEN                                                  
       DO 60 I = 1,(N / NE)                                             
        DO 60 J = 1,M                                                   
         B(I,J) = BTIJ1(TI(I),J) / GERRTI(I)                            
60     CONTINUE                                                         
       IF (NE.EQ.2) THEN                                                
        DO 70 I = 1,(N / 2)                                             
         IWK = I + (N / 2)                                              
         DO 70 J = 1,M                                                  
          B(IWK,J) = BTIJ2(TI(I),J) / GERRTI(IWK)                       
70      CONTINUE                                                        
       ENDIF                                                            
      ENDIF                            
                                                                       
c-----SETUP MATRIX FOR THE APPROXIMATED REGULARIZATION OPERATOR-------- 
                                                                        
      DO 80 I = 1,NLMAX                                                 
       DO 80 J = 1,NS                                                   
        L(I,J) = 0.D0                                                   
80    CONTINUE                                                          
                                                                        
      IF (MOD(REGMOD,10).EQ.1) THEN                                     
       NL = NS                                                          
       DPWK = DSQRT(HS)                                                 
       DO 90 I = 1,NL                                                   
        L(I,I) = DPWK                                                   
90     CONTINUE                                                         
      ELSE IF (MOD(REGMOD,10).EQ.2) THEN                                
       NL = NS - 2                                                      
       DPWK = DSQRT(HS) / HS**2                                         
       DO 100 I = 1,NL                                                  
        L(I,I) = DPWK                                                   
        L(I,I + 2) = DPWK                                               
        L(I,I + 1) = -2.D0 * DPWK                                       
100    CONTINUE                                                         
      ELSE                                                              
       NL = NS + 2                                                      
       DPWK = DSQRT(HS) / HS**2                                         
       DO 110 I = 3,NL - 2                                              
        L(I,I) = DPWK                                                   
        L(I,I - 2) = DPWK                                               
        L(I,I - 1) = -2.D0 * DPWK                                       
110    CONTINUE                                                         
       L(1,1) = DPWK                                                    
       L(2,1) = -2.D0 * DPWK                                            
       L(2,2) = DPWK                                                    
       L(NL - 1,NL - 3) = DPWK                                          
       L(NL - 1,NL - 2) = -2.D0 * DPWK                                  
       L(NL,NL - 2) = DPWK                                              
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
c---------------------------------------------------------------------- 
c SUBROUTINE: R(EA)D PARA(METER)                                        
c                                                                       
c PURPOSE   : READ PARAMETER FROM INPUT FILE #10                        
c                                                                       
c VARIABLES : (ON INPUT)                                                
c             NSMAX  : MAXIMUM NUMBER OF POINTS WHERE THE SOL. CAN BE  
c                        CALCULATED                                     
c             MMAX   : MAXIMUM NUMBER OF COEFFICIENTS AJ                
c             NMAX   : MAXIMUM NUMBER OF DATA POINTS                    
c                                                                       
c             (ON OUTPUT)                                               
c             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
c             SMIN,                                                     
c             SMAX   : INTERVAL IN WHICH THE SOLUTION IS CALCULATED     
c                                                                       
c             M      : NUMBER OF COEFFICIENTS AJ                        
c             N      : NUMBER OF DATA POINTS                            
c             NE     : NUMBER OF DIFFERENT EQUATIONS                    
c                                                                       
c             ERROR  : SCALEFACTOR OF THE ERROR                         
c             LAMBDA : REGULARIZATION PARAMETER                         
c                                                                       
c             DISMOD : MODE FOR THE DISCRETIZATION OF THE INTEGRALS     
c             ERRMOD : MODE FOR THE ESTIMATION OF THE ERROR             
c             REGMOD : MODE FOR THE CALCULATION OF THE SOLUTION         
c             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
c                                                                        
c            LAMB?? : PARAMETER FOR THE CALCULATION OF LAMBDA          
c---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE RDPARA(NS,NSMAX,SMIN,SMAX,M,MMAX,N,NMAX,NE,ERROR,      
     *                  LAMBDA,DISMOD,ERRMOD,REGMOD,INFMOD,LAMBST,      
     *                  LAMBSP,LAMBRA,LAMBPR,LAMBIT)                    
                                                                        
      INTEGER NS,NSMAX                                                  
      DOUBLE PRECISION SMIN,SMAX                                        
                                                                        
      INTEGER M,MMAX                                                    
      INTEGER N,NMAX                                                    
      INTEGER NE                                                        
                                                                        
      DOUBLE PRECISION ERROR,LAMBDA                                     
                                                                        
      INTEGER DISMOD,ERRMOD,REGMOD,INFMOD                               
                                                                        
      INTEGER LAMBIT                                                    
      DOUBLE PRECISION LAMBST,LAMBSP,LAMBRA,LAMBPR                      
                                                                        
      CHARACTER*60 ERRTXT                                               
                                                                        
c-----READ INPUT PARAMETER FROM FILE----------------------------------- 
                                                                        
      OPEN(UNIT=10,FILE='LTS.par',STATUS='OLD')                         
                                                                        
      READ(10,*)NS                                                      
      READ(10,*)SMIN                                                    
      READ(10,*)SMAX                                                    
      READ(10,*)DISMOD                                                  
      READ(10,*)M                                                       
      READ(10,*)N                                                       
      READ(10,*)NE                                                      
      READ(10,*)ERRMOD                                                  
      READ(10,*)ERROR                                                   
      READ(10,*)REGMOD                                                  
      READ(10,*)LAMBDA                                                  
      READ(10,*)INFMOD                                                  
      READ(10,*)LAMBST                                                  
      READ(10,*)LAMBSP                                                  
      READ(10,*)LAMBRA                                                  
      READ(10,*)LAMBPR                                                  
      READ(10,*)LAMBIT                                                  
                                                                        
      CLOSE(10)                                                         
                                                                        
c-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' RDPARA > CALLED'                                     
      ENDIF                                                             
                                                                        
c-----TEST INPUT PARAMETER FOR CORRECTNESS----------------------------- 
                                                                        
      IF (NS.GT.NSMAX) THEN                                             
       ERRTXT = ' RDPARA > NUMBER OF GRIDPOINTS TOO LARGE (NS)'         
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (NS.LT.2) THEN                                                 
       ERRTXT = ' RDPARA > ILLEGAL NUMBER OF GRIDPOINTS (NS)'           
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (SMIN.GE.SMAX) THEN                                            
       ERRTXT = ' RDPARA > NO INTERVAL DEFINED (SMIN,SMAX)'             
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF ((DISMOD.NE.1).AND.(DISMOD.NE.2)) THEN                         
       ERRTXT = ' RDPARA > ILLEGAL DISCRETIZATION MODE (DISMOD)'        
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (M.GT.MMAX) THEN                                               
       ERRTXT = ' RDPARA > NUMBER OF COEFFICIENTS AJ TOO LARGE (M)'     
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (M.LT.0) THEN                                                  
       ERRTXT = ' RDPARA > ILLEGAL NUMBER OF COEFFICIENTS AJ (M)'       
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (N.GT.NMAX) THEN                                               
       ERRTXT = ' RDPARA > NUMBER OF DATA POINTS TOO LARGE (N)'         
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (N.LE.(M + 2)) THEN                                            
       ERRTXT = ' RDPARA > ILLEGAL NUMBER OF DATA POINTS (N)'           
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF ((NE.LT.1).OR.(NE.GT.2)) THEN                                  
       ERRTXT = ' RDPARA > ILLEGAL NUMBER OF DIFFERENT EQUATIONS (NE)'  
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (MOD(N,NE).NE.0) THEN                                          
       ERRTXT = ' RDPARA > WRONG NUMBER OF DATA POINTS (N)'             
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF ((MOD(ERRMOD,10).LT.1).OR.(MOD(ERRMOD,10).GT.3).OR.            
     *    ((ERRMOD / 10).LT.1).OR.((ERRMOD / 10).GT.2)) THEN           
       ERRTXT = ' RDPARA > ILLEGAL ERROR MODE (ERRMOD)'                 
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (((ERRMOD / 10).EQ.1).AND.(ERROR.LE.0.D0)) THEN                
       ERRTXT = ' RDPARA > ILLEGAL ESTIMATE OF THE ERROR (ERROR)'       
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF ((MOD(REGMOD,10).LT.1).OR.(MOD(REGMOD,10).GT.3).OR.            
     *    (MOD(REGMOD / 10,10).LT.1).OR.(MOD(REGMOD / 10,10).GT.2).OR.  
     *    ((REGMOD / 100).LT.1).OR.((REGMOD / 100).GT.2)) THEN          
       ERRTXT = ' RDPARA > ILLEGAL REGULARIZATION MODE (REGMOD)'        
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF ((MOD(REGMOD / 10,10).EQ.1).AND.(LAMBDA.LT.0.D0)) THEN         
       ERRTXT = ' RDPARA > ILLEGAL VALUE FOR THE REG. PARA. (LAMBDA)'   
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF ((INFMOD.LT.1).OR.(INFMOD.GT.3)) THEN                          
       ERRTXT = ' RDPARA > ILLEGAL INFORMATION MODE (INFMOD)'           
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (LAMBSP.LE.0.D0) THEN                                          
       ERRTXT = ' RDPARA > ILLEGAL VALUE (LAMBSP)'                      
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (LAMBRA.LE.LAMBSP) THEN                                        
       ERRTXT = ' RDPARA > ILLEGAL VALUE (LAMBRA)'                      
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (LAMBPR.LE.0.D0) THEN                                          
       ERRTXT = ' RDPARA > ILLEGAL VALUE (LAMBPR)'                      
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      IF (LAMBIT.LE.0) THEN                                             
       ERRTXT = ' RDPARA > ILLEGAL VALUE (LAMBIT)'                      
       CALL WEERRM(ERRTXT,.TRUE.)                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
c---------------------------------------------------------------------- 
c SUBROUTINE: R(EA)D DATA                                               
c                                                                       
c PURPOSE   : READ DATA POINTS FROM INPUT FILE #11                      
c                                                                       
c VARIABLES : (ON INPUT)                                                
c             N      : NUMBER OF DATA POINTS                            
c             NE     : NUMBER OF DIFFERENT EQUATIONS                    
c                                                                       
c             ERRMOD : MODE FOR THE ESTIMATION OF THE ERROR             
c             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
c                                                                       
c             (ON OUTPUT)                                               
c             TI     : ABSZISSES OF THE DATA POINTS                    
c             GSTI   : ORDINATE OF THE DATA POINTS                     
c             GERRTI : ERROR OF THE DATA POINTS                         
c---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE RDDATA(N,NE,TI,GSTI,GERRTI,ERRMOD,INFMOD)              
                                                                        
      INTEGER N,NE                                                      
      DOUBLE PRECISION TI(N),GSTI(N),GERRTI(N)                          
                                                                        
      INTEGER ERRMOD,INFMOD                                             
                                                                        
      INTEGER I                                                         
                                                                        
      INTEGER IWK                                                       
                                                                        
c-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' RDDATA > CALLED'                                     
      ENDIF                                                             
                                                                        
c-----READ DATA POINTS FROM FILE--------------------------------------- 
                                                                        
      OPEN(UNIT=11,FILE='LTS.dat',STATUS='OLD')                         
                                                                        
      IF (NE.EQ.1) THEN                                                 
                                                                        
       DO 10 I = 1,N                                                    
        IF (MOD(ERRMOD,10).EQ.1) THEN                                   
         READ(11,*)TI(I),GSTI(I),GERRTI(I)                              
        ELSE IF (MOD(ERRMOD,10).EQ.2) THEN                              
         READ(11,*)TI(I),GSTI(I)                                        
         GERRTI(I) = 1.D0                                               
        ELSE IF (MOD(ERRMOD,10).EQ.3) THEN                              
         READ(11,*)TI(I),GSTI(I)                                        
         GERRTI(I) = GSTI(I)                                            
        ENDIF                                                           
10     CONTINUE                                                         
                                                                        
      ELSE                                                          
                                                                        
       DO 20 I = 1,(N / 2)                                              
        IWK = I + (N / 2)                                               
        IF (MOD(ERRMOD,10).EQ.1) THEN                                   
         READ(11,*)TI(I),GSTI(I),GERRTI(I),GSTI(IWK),GERRTI(IWK)        
        ELSE IF (MOD(ERRMOD,10).EQ.2) THEN                              
         READ(11,*)TI(I),GSTI(I),GSTI(IWK)                              
         GERRTI(I) = 1.D0                                               
         GERRTI(IWK) = 1.D0                                             
        ELSE IF (MOD(ERRMOD,10).EQ.3) THEN                              
         READ(11,*)TI(I),GSTI(I),GSTI(IWK)                              
         GERRTI(I) = GSTI(I)                                            
         GERRTI(IWK) = GSTI(IWK)                                        
        ENDIF                                                           
20     CONTINUE                                                         
                                                                        
      ENDIF                                                             
                                                                        
      CLOSE(11)                                                         
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
c---------------------------------------------------------------------- 
c SUBROUTINE: W(RIT)E RE(GULARIZED) SO(LUTION)                          
c                                                                       
c PURPOSE   : WRITE REGULARIZED SOLUTION ON FILE #12                    
c                                                                       
c VARIABLES : (ON INPUT)                                                
c             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
c             SJ     : ABSZISSES WHERE THE SOL. IS CALCULATED           
c             FSSJ   : VALUES OF THE CALCULATED SOLUTION                
c             FMSJ   : MIDPOINT FOR THE ERRORS OF THE CALCULATED SOL.   
c             FERRSJ : ERRORS OF THE CALCULATED SOLUTION                
c---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE WERESO(NS,SJ,FSSJ,FMSJ,FERRSJ)                         
                                                                        
      INTEGER NS                                                        
      DOUBLE PRECISION SJ(NS),FSSJ(NS),FMSJ(NS),FERRSJ(NS)              
                                                                        
      INTEGER I                                                         
                                                                        
      OPEN(UNIT=12,FILE='LTS.sol',STATUS='UNKNOWN')                     
                                                                        
      DO 10 I = 1,NS                                                    
       WRITE(12,1001)SJ(I),FSSJ(I),FMSJ(I),FERRSJ(I)                    
10    CONTINUE                                                          
                                                                        
      CLOSE(12)                                                         
                                                                        
      RETURN                                                            
                                                                        
1001  FORMAT(4D18.8)                                                    
                                                                        
      END                                                               
                                                                        
                                                                        
c---------------------------------------------------------------------- 
c SUBROUTINE: W(RIT)E CO(EFFICIENTS) AJ                                 
c                                                                       
c PURPOSE   : WRITE COEFFICIENTS AJ ON FILE #13                         
c                                                                       
c VARIABLES : (ON INPUT)                                                
c             M      : NUMBER OF COEFFICIENTS AJ                        
c             ASJ    : VALUES OF THE CALCULATED COEFFICIENTS AJ         
c             AERRJ  : ERRORS OF THE CALCULATED COEFFICIENTS AJ         
c---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE WECOAJ(M,ASJ,AERRJ)                                    
                                                                        
      INTEGER M                                                         
      DOUBLE PRECISION ASJ(M),AERRJ(M)                                  
                                                                        
      INTEGER I                                                         
                                                                        
      IF (M.GT.0) THEN                                                  
       OPEN(UNIT=13,FILE='LTS.adp',STATUS='UNKNOWN')                    
                                                                        
       DO 10 I = 1,M                                                    
        WRITE(13,1001)I,ASJ(I),AERRJ(I)                                 
10     CONTINUE                                                         
                                                                        
       CLOSE(13)                                                        
                                                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
                                                                        
1001  FORMAT(I5,2D18.8)                                                 
                                                                        
      END                                                               
                                                                        
                                                                        
c---------------------------------------------------------------------- 
c SUBROUTINE: CA(LCULATE) R(EGULARIZED) S(OLUTION WITH) SC(-METHOD)     
c                                                                       
c PURPOSE   : SOLVE DISCRETE REGULARIZATION PROBLEM                     
c                                                                       
c VARIABLES : (ON INPUT)                                                
c             N      : NUMBER OF DATA POINTS                            
c             GSTI   : ORDINATE OF THE DATA POINTS                      
c                                                                       
c             M      : NUMBER OF COEFFICIENTS AJ                        
c             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
c                                                                       
c             B      : ADDITIONAL LINEAR TERMS                          
c             NBMAX  : LEADING DIMENSION OF THE FIELD B                 
c                                                                       
c             K      : APPROXIMATED INTEGRAL OPERATOR                   
c             NKMAX  : LEADING DIMENSION OF THE FIELD K                 
c                                                                       
c             NL     : FIRST DIMENSION OF THE REG. OPERATOR             
c             L      : APPROXIMATED REGULARIZATION OPERATOR             
c             NLMAX  : LEADING DIMENSION OF THE FIELD L                 
c                                                                       
c             ERROR  : SCALEFACTOR OF THE ERROR                        
c                                                                       
c             SCLMOD : MODE FOR THE ESTIMATION OF THE ERROR             
c             POSMOD : MODE FOR THE CALCULATION OF A POSITIVE SOLUTION  
c             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
c                                                                       
c             LAMB?? : PARAMETER FOR THE CALCULATION OF LAMBDA          
c                                                                       
c             IFWK   : INTEGER WORKSPACE                                
c             DPFWK  : DOUBLE PRECISION WORKSPACE                       
c                                                                       
c             (ON OUTPUT)                                               
c             ASJ    : VALUES OF THE CALCULATED COEFFICIENTS AJ         
c             AERRJ  : ERRORS OF THE CALCULATED COEFFICIENTS AJ         
c                                                                       
c             FSSJ   : VALUES OF THE CALCULATED SOLUTION                
c             FMSJ   : MIDPOINT FOR THE ERRORS OF THE CALCULATED SOL.   
c             FERRSJ : ERRORS OF THE CALCULATED SOLUTION                
c                                                                       
c             ERROR  : SCALEFACTOR OF THE ERROR OBTAINED BY SC-METHOD   
c             LAMBDA : REGULARIZATION PARAMETER OBTAINED BY SC-METHOD   
c                                                                       
c COMMENT   : THE SECOND DIMENSION OF K MUST BE AT LEAST (M + NS)       
c                                                                       
c             THE SIZE OF THE INTEGER WORKSPACE MUST BE AT LEAST        
c             MAX(2(M + NS) + NEFF + NL,N)                              
c                                                                       
c             THE SIZE OF THE DOUBLE PRECISION WORKSPACE MUST BE AT     
c             LEAST (M + NS) (M + NS + NEFF + NL + 4) + NEFF + NL       
c                                                                       
c             (NEFF = MIN(N,M + NS))                                    
c---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE CARSSC(N,GSTI,M,ASJ,AERRJ,NS,FSSJ,FMSJ,FERRSJ,B,NBMAX, 
     *                  K,NKMAX,NL,L,NLMAX,ERROR,LAMBDA,SCLMOD,POSMOD,  
     *                  INFMOD,LAMBST,LAMBSP,LAMBRA,LAMBPR,LAMBIT,      
     *                  IFWK,DPFWK)                                     
                                                                        
      INTEGER N                                                         
      DOUBLE PRECISION GSTI(N)                                          
                                                                        
      INTEGER M                                                         
      DOUBLE PRECISION ASJ(M),AERRJ(M)                                  
                                                                        
      INTEGER NS                                                        
      DOUBLE PRECISION FSSJ(NS),FMSJ(NS),FERRSJ(NS)                     
                                                                        
      INTEGER NBMAX                                                     
      DOUBLE PRECISION B(NBMAX,M)                                       
                                                                        
      INTEGER NKMAX                                                     
      DOUBLE PRECISION K(NKMAX,M + NS)                                  
                                                                        
      INTEGER NL,NLMAX                                                  
      DOUBLE PRECISION L(NLMAX,NS)                                      
                                                                        
      DOUBLE PRECISION ERROR,LAMBDA                                     
                                                                        
      INTEGER SCLMOD,POSMOD,INFMOD                                      
                                                                        
      INTEGER LAMBIT                                                    
      DOUBLE PRECISION LAMBST,LAMBSP,LAMBRA,LAMBPR                      
                                                                        
      INTEGER IFWK(*)                                                   
      DOUBLE PRECISION DPFWK(*)                                         
                                                                        
      DOUBLE PRECISION CANORM,CALAMB                                    
                                                                        
      LOGICAL FLAG                                                      
                                                                        
      INTEGER NEFF                                                      
                                                                        
      INTEGER NX,MX                                                     
                                                                        
      DOUBLE PRECISION BSCALE,LSCALE                                    
                                                                        
      DOUBLE PRECISION GSTIRQ                                           
                                                                        
      INTEGER I                                                         
                                                                     
      DOUBLE PRECISION DPDUM                                            
      DOUBLE PRECISION DPFDUM(1)                                        
                                                                        
      INTEGER IWK1,IWK2,IWK3,IWK4,IWK5,IWK6,IWK7                        
      DOUBLE PRECISION DPWK                                             
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' CARSSC > CALLED'                                     
      ENDIF                                                             
                                                                        
*-----INITIALIZE SCALEFACTOR FOR THE REG.PARA-------------------------- 
                                                                        
      DPWK = CANORM(N,NS,K,NKMAX)                                       
      LSCALE = DPWK / CANORM(NL,NS,L,NLMAX)                             
                                                                        
*-----SETUP MATRIX COMPOSED OUT OF B AND K----------------------------- 
                                                                        
      IF (M.GT.0) THEN                                                  
       BSCALE = DPWK / CANORM(N,M,B,NBMAX)                              
       CALL SEMABK(N,M,B,NBMAX,BSCALE,NS,K,NKMAX,INFMOD)                
      ENDIF                                                             
                                                                        
*-----REDUCE EFFECTIVE NUMBER OF DATA POINTS--------------------------- 
                                                                        
      CALL REENDP(N,GSTI,M,NS,K,NKMAX,NEFF,IFWK,DPFWK,INFMOD)           
                                                                        
      GSTIRQ = 0.D0                                                     
      DO 10 I = (NEFF + 1),N                                            
       GSTIRQ = GSTIRQ + GSTI(I) * GSTI(I)                              
10    CONTINUE                                                          
                                                                        
*-----INITIALIZE PARAMETER FOR WORKSPACE MANAGEMENT-------------------- 
                                                                        
      IWK1 = (M + NS) * (M + NS) + 1                                    
      IWK2 = IWK1 + (M + NS)                                            
      IWK3 = IWK2 + (NEFF + NL)                                         
      IWK4 = IWK3 + (M + NS)                                            
      IWK5 = IWK4 + (M + NS)                                            
                                                                        
      IWK6 = 1 + (M + NS)                                               
      IWK7 = IWK6 + (M + NS)                                            
                                                                        
      DO 20 I = 1,(M + NS)                                              
       IFWK(I) = 0                                                      
       IFWK(I + IWK6 - 1) = I                                           
20    CONTINUE                                                          
                                                                        
*-----CAL. REG. SOL. NEGLECTING POSITIVITY CONSTRAINTS----------------- 
                                                                        
      MX = M + NS                                                       
      CALL INMADE(NEFF,GSTI,M,NS,K,NKMAX,LSCALE,NL,L,NLMAX,IFWK(IWK6),  
     *            NX,MX,DPFWK(1),M + NS,DPFWK(IWK1),DPFWK(IWK2),        
     *            .TRUE.,IFWK(IWK7),DPFWK(IWK5),INFMOD)                 
                                                                        
      LAMBDA = CALAMB(M,N,NX,MX,DPFWK(1),M + NS,DPFWK(IWK1),            
     *                DPFWK(IWK2),GSTIRQ,ERROR,SCLMOD,LAMBST,LAMBSP,    
     *                LAMBRA,LAMBPR,LAMBIT,DPFWK(IWK5),INFMOD)          
      LAMBDA = 10.D0**LAMBDA                                            
                                                                        
      CALL CAFSMD(LAMBDA,NX,MX,DPFWK(1),M + NS,DPFWK(IWK1),             
     *            DPFWK(IWK2),DPFWK(IWK3),.TRUE.,DPFWK(IWK5),           
     *            INFMOD)                                               
                                                                        
      FLAG = .FALSE.                                                    
                                                                        
      IF (POSMOD.EQ.2) THEN                                             
                                                                        
*-----CAL. REG. SOL. WITH POSITIVITY CONSTRAINTS----------------------- 
                                                                        
       CALL INPAQP(M,MX,DPFWK(IWK3),IFWK(1),FLAG,INFMOD)                
       IF (FLAG) THEN                                                   
                                                                        
        LSCALE = DSQRT(LAMBDA) * LSCALE                                 
        CALL CAFSQP(NEFF,GSTI,M,NS,K,NKMAX,LSCALE,NL,L,NLMAX,MX,        
     *              IFWK(1),IFWK(IWK6),DPFWK(1),M + NS,IFWK(IWK7),      
     *              DPFWK(IWK5),DPFWK(IWK2),DPFWK(IWK3),INFMOD)         
                                                                        
        CALL INMADE(NEFF,GSTI,M,NS,K,NKMAX,LSCALE,NL,L,NLMAX,           
     *              IFWK(IWK6),NX,MX,DPFWK(1),M + NS,DPFWK(IWK1),       
     *              DPFWK(IWK2),.FALSE.,IFWK(IWK7),DPFWK(IWK5),         
     *              INFMOD)                                             
                                                                        
        LAMBDA = CALAMB(M,N,NX,MX,DPFWK(1),M + NS,DPFWK(IWK1),          
     *                  DPFWK(IWK2),DPDUM,ERROR,1,0.D0,LAMBSP,          
     *                  LAMBRA,LAMBPR,LAMBIT,DPFWK(IWK5),INFMOD)        
        LAMBDA = 10.D0**LAMBDA                                          
                                                                        
        CALL CAFSMD(LAMBDA,NX,MX,DPFWK(1),M + NS,DPFWK(IWK1),           
     *              DPFWK(IWK2),DPFWK(IWK3),.TRUE.,DPFWK(IWK5),         
     *              INFMOD)                                             
                                                                        
        DPWK = DSQRT(LAMBDA) * LSCALE                                   
        CALL TEACCO(NEFF,GSTI,M,NS,K,NKMAX,DPWK,NL,L,NLMAX,             
     *              MX,DPFWK(IWK3),IFWK(IWK6),FLAG,DPFWK(IWK5),         
     *              DPFWK(IWK2),INFMOD)                                 
                                                                        
        IF (FLAG) THEN                                                  
         LSCALE = DSQRT(LAMBDA) * LSCALE                                
         LAMBDA = 1.D0                                                  
         CALL CAFSQP(NEFF,GSTI,M,NS,K,NKMAX,LSCALE,NL,L,NLMAX,MX,       
     *               IFWK(1),IFWK(IWK6),DPFWK(1),M + NS,IFWK(IWK7),     
     *               DPFWK(IWK5),DPFWK(IWK2),DPFWK(IWK3),INFMOD)        
                                                                        
         CALL INREIO(NEFF,M,NS,K,NKMAX,LSCALE,NL,L,NLMAX,IFWK(IWK6),    
     *               NX,MX,DPFWK(1),M + NS,.FALSE.,IFWK(IWK7),          
     *               DPFWK(IWK5),INFMOD)                                
                                                                        
         CALL CAFSMD(DPDUM,NX,MX,DPFWK(1),M + NS,DPFDUM,GSTI,           
     *               DPFWK(IWK3),.FALSE.,DPFDUM,INFMOD)                 
                                                                        
         CALL CAERMD(DPDUM,ERROR,NX,MX,DPFWK(1),M + NS,DPFDUM,          
     *               DPFWK(IWK4),.FALSE.,DPFDUM,INFMOD)                 
                                                                        
        ENDIF                                                           
                                                                        
       ENDIF                                                            
                                                                        
      ENDIF                                                             
                                                                        
      IF (.NOT.FLAG) THEN                                               
       CALL CAERMD(LAMBDA,ERROR,NX,MX,DPFWK(1),M + NS,DPFWK(IWK1),      
     *             DPFWK(IWK4),.TRUE.,DPFWK(IWK5),INFMOD)               
      ENDIF                                                             
                                                                        
      CALL RERESO(MX,IFWK(IWK6),DPFWK(IWK3),DPFWK(IWK4),BSCALE,M,ASJ,   
     *            AERRJ,NS,FSSJ,FMSJ,FERRSJ,POSMOD,INFMOD)              
                                                                        
      LAMBDA = LAMBDA * LSCALE**2                                       
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: CA(LCULATE) R(EGULARIZED) S(OLUTION FOR) F(IXED)          
*             R(EGULARIZATION PARAMETER)                                
*                                                                       
* PURPOSE   : SOLVE DISCRETE REGULARIZATION PROBLEM                     
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             N      : NUMBER OF DATA POINTS                            
*             GSTI   : ORDINATE OF THE DATA POINTS                      
*                                                                       
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
*                                                                       
*             B      : ADDITIONAL LINEAR TERMS                          
*             NBMAX  : LEADING DIMENSION OF THE FIELD B                 
*                                                                       
*             K      : APPROXIMATED INTEGRAL OPERATOR                   
*             NKMAX  : LEADING DIMENSION OF THE FIELD K                 
*                                                                       
*             NL     : FIRST DIMENSION OF THE REG. OPERATOR             
*             L      : APPROXIMATED REGULARIZATION OPERATOR             
*             NLMAX  : LEADING DIMENSION OF THE FIELD L                 
*                                                                       
*             ERROR  : SCALEFACTOR OF THE ERROR                         
*             LAMBDA : REGULARIZATION PARAMETER                        
*                                                                       
*             POSMOD : MODE FOR THE CALCULATION OF A POSITIVE SOLUTION  
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
*             IFWK   : INTEGER WORKSPACE                                
*             DPFWK  : DOUBLE PRECISION WORKSPACE                       
*                                                                       
*             (ON OUTPUT)                                               
*             ASJ    : VALUES OF THE CALCULATED COEFFICIENTS AJ         
*             AERRJ  : ERRORS OF THE CALCULATED COEFFICIENTS AJ         
*                                                                       
*             FSSJ   : VALUES OF THE CALCULATED SOLUTION                
*             FMSJ   : MIDPOINT FOR THE ERRORS OF THE CALCULATED SOL.   
*             FERRSJ : ERRORS OF THE CALCULATED SOLUTION                
*                                                                       
* COMMENT   : THE SECOND DIMENSION OF K MUST BE AT LEAST (M + NS)       
*                                                                       
*             THE SIZE OF THE DOUBLE PRECISION WORKSPACE MUST BE AT     
*             LEAST (M + NS) (M + NS + NEFF + NL + 2) + NEFF + NL       
*                                                                       
*             THE SIZE OF THE INTEGER WORKSPACE MUST BE AT LEAST        
*             MAX(2(M + NS) + NEFF + NL,N)                              
*                                                                       
*             (NEFF = MIN(N,M + NS))                                    
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE CARSFR(N,GSTI,M,ASJ,AERRJ,NS,FSSJ,FMSJ,FERRSJ,B,NBMAX, 
     *                  K,NKMAX,NL,L,NLMAX,ERROR,LAMBDA,POSMOD,INFMOD,  
     *                  IFWK,DPFWK)                                     
                                                                        
      INTEGER N                                                         
      DOUBLE PRECISION GSTI(N)                                          
                                                                        
      INTEGER M                                                         
      DOUBLE PRECISION ASJ(M),AERRJ(M)                                  
                                                                        
      INTEGER NS                                                        
      DOUBLE PRECISION FSSJ(NS),FMSJ(NS),FERRSJ(NS)                     
                                                                        
      INTEGER NBMAX                                                     
      DOUBLE PRECISION B(NBMAX,M)                                       
                                                                        
      INTEGER NKMAX                                                     
      DOUBLE PRECISION K(NKMAX,M + NS)                                  
                                                                        
      INTEGER NL,NLMAX                                                  
      DOUBLE PRECISION L(NLMAX,NS)                                      
                                                                        
      DOUBLE PRECISION ERROR,LAMBDA                                     
                                                                        
      INTEGER POSMOD,INFMOD                                             
                                                                        
      INTEGER IFWK(*)                                                   
      DOUBLE PRECISION DPFWK(*)                                         
                                                                        
      DOUBLE PRECISION CANORM                                           
                                                                        
      INTEGER NEFF                                                      
                                                                        
      INTEGER NX,MX                                                     
                                                                        
      DOUBLE PRECISION BSCALE,LSCALE                                    
                                                                        
      INTEGER I                                                         
                                                                        
      DOUBLE PRECISION DPDUM                                            
      DOUBLE PRECISION DPFDUM(1)                                        
                                                                        
      INTEGER IWK1,IWK2,IWK3,IWK4,IWK5,IWK6                             
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' CARSFR > CALLED'                                     
      ENDIF                                                             
                                                                        
*-----SETUP MATRIX COMPOSED OUT OF B AND K----------------------------- 
                                                                        
      IF (M.GT.0) THEN                                                  
       BSCALE = CANORM(N,NS,K,NKMAX) / CANORM(N,M,B,NBMAX)              
       CALL SEMABK(N,M,B,NBMAX,BSCALE,NS,K,NKMAX,INFMOD)                
      ENDIF                                                             
                                                                        
*-----REDUCE EFFECTIVE NUMBER OF DATA POINTS--------------------------- 
                                                                        
      CALL REENDP(N,GSTI,M,NS,K,NKMAX,NEFF,IFWK,DPFWK,INFMOD)           
                                                                        
*-----INITIALIZE PARAMETER FOR WORKSPACE MANAGEMENT-------------------- 
                                                                        
      IWK1 = (M + NS) * (M + NS) + 1                                    
      IWK2 = IWK1 + NEFF + NL                                           
      IWK3 = IWK2 + (M + NS)                                            
      IWK4 = IWK3 + (M + NS)                                            
                                                                        
      IWK5 = 1 + (M + NS)                                               
      IWK6 = IWK5 + (M + NS)                                            
                                                                        
      DO 10 I = 1,(M + NS)                                              
       IFWK(I) = 0                                                      
       IFWK(I + IWK5 - 1) = I                                           
10    CONTINUE                                                          
                                                                        
*-----CAL. REG. SOL. NEGLECTING POSITIVITY CONSTRAINTS----------------- 
                                                                        
      LSCALE = DSQRT(LAMBDA)                                            
                                                                        
      IF (POSMOD.EQ.1) THEN                                             
                                                                        
       MX = M + NS                                                      
                                                                        
       CALL INREIO(NEFF,M,NS,K,NKMAX,LSCALE,NL,L,NLMAX,IFWK(IWK5),      
     *             NX,MX,DPFWK(1),M + NS,.TRUE.,IFWK(IWK6),             
     *             DPFWK(IWK4),INFMOD)                                  
                                                                        
      ELSE                                                              
                                                                        
       CALL CAFSQP(NEFF,GSTI,M,NS,K,NKMAX,LSCALE,NL,L,NLMAX,MX,         
     *             IFWK(1),IFWK(IWK5),DPFWK(1),M + NS,IFWK(IWK6),       
     *             DPFWK(IWK4),DPFWK(IWK1),DPFWK(IWK2),INFMOD)          
                                                                        
       CALL INREIO(NEFF,M,NS,K,NKMAX,LSCALE,NL,L,NLMAX,IFWK(IWK5),      
     *             NX,MX,DPFWK(1),M + NS,.FALSE.,IFWK(IWK6),            
     *             DPFWK(IWK4),INFMOD)                                  
                                                                        
      ENDIF                                                             
                                                                        
      CALL CAFSMD(DPDUM,NX,MX,DPFWK(1),M + NS,DPFDUM,GSTI,              
     *            DPFWK(IWK2),.FALSE.,DPFDUM,INFMOD)                    
                                                                        
      CALL CAERMD(DPDUM,ERROR,NX,MX,DPFWK(1),M + NS,DPFDUM,             
     *            DPFWK(IWK3),.FALSE.,DPFDUM,INFMOD)                    
                                                                        
      CALL RERESO(MX,IFWK(IWK5),DPFWK(IWK2),DPFWK(IWK3),BSCALE,M,ASJ,   
     *            AERRJ,NS,FSSJ,FMSJ,FERRSJ,POSMOD,INFMOD)              
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: SE(T UP) MA(TRIX COMPOSED OUT OF) B (AND) K               
*                                                                       
* PURPOSE   : SETUP MATRIX (B,K)                                        
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             N      : NUMBER OF DATA POINTS                            
*                                                                       
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             B      : ADDITIONAL LINEAR TERMS                          
*             NBMAX  : LEADING DIMENSION OF THE FIELD B                 
*             BSCALE : SCALEFACTOR FOR THE ADDITIONAL LINEAR TERMS      
*                                                                       
*             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
*             K      : APPROXIMATED INTEGRAL OPERATOR                   
*             NKMAX  : LEADING DIMENSION OF THE FIELD K                 
*                                                                       
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
*             (ON OUTPUT)                                               
*             K      : COMPOSED MATRIX                                  
*---------------------------------------------------------------------- 
                                                                        
       SUBROUTINE SEMABK(N,M,B,NBMAX,BSCALE,NS,K,NKMAX,INFMOD)          
                                                                        
       INTEGER N                                                        
                                                                        
       INTEGER M,NBMAX                                                  
       DOUBLE PRECISION BSCALE                                          
       DOUBLE PRECISION B(NBMAX,M)                                      
                                                                        
       INTEGER NS,NKMAX                                                 
       DOUBLE PRECISION K(NKMAX,M + NS)                                 
                                                                        
       INTEGER INFMOD                                                   
                                                                        
       INTEGER I,J                                                      
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' SEMABK > CALLED'                                     
      ENDIF                                                             
                                                                        
*------SETUP MATRIX COMPOSED OUT OF B AND K---------------------------- 
                                                                        
       DO 10 I = NS,1,-1                                                
        DO 10 J = 1,N                                                   
         K(J,I + M) = K(J,I)                                            
10     CONTINUE                                                         
                                                                        
       DO 20 I = 1,M                                                    
        DO 20 J = 1,N                                                   
         K(J,I) = BSCALE * B(J,I)                                       
20     CONTINUE                                                         
                                                                        
       RETURN                                                           
       END                                                              
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: SE(TUP) MA(TRIX COMPOSED OUT OF B,) K (AND) L             
*                                                                       
* PURPOSE   : SETUP MATRIX COMPOSED OF K, B AND L                       
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             NEFF   : EFFECTIVE NUMBER OF DATA POINTS                  
*                                                                       
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
*             BK     : MATRIX COMPOSED OUT OF B AND K                   
*             NBKMAX : LEADING DIMENSION OF THE FIELD BK                
*                                                                       
*             LSCALE : SCALEFACTOR FOR THE APP. REG. OPERATOR           
*             NL     : FIRST DIMENSION OF THE REG. OPERATOR             
*             L      : APPROXIMATED REGULARIZATION OPERATOR             
*             NLMAX  : LEADING DIMENSION OF THE FIELD L                 
*                                                                       
*             MX     : FIRST DIMENSION OF THE MATRIX X                  
*                                                                       
*             ORDER  : ORDER OF THE PARAMETER                           
*                                                                       
*             NXMAX  : LEADING DIMENSION OF THE FIELD X                 
*                                                                       
*             (ON OUTPUT)                                               
*             X      : MATRIX COMPOSED OF K, B AND (OPTIONALLY) L       
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE SEMAKL(NEFF,M,NS,BK,NBKMAX,LSCALE,NL,L,NLMAX,MX,ORDER, 
     *                  X,NXMAX)                                        
                                                                        
      INTEGER NEFF                                                      
                                                                        
      INTEGER M,NS,NBKMAX                                               
      DOUBLE PRECISION BK(NBKMAX,M + NS)                                
                                                                        
      INTEGER NL,NLMAX                                                  
      DOUBLE PRECISION LSCALE                                           
      DOUBLE PRECISION L(NLMAX,NL)                                      
                                                                        
      INTEGER MX                                                        
      INTEGER ORDER(MX)                                                 
                                                                        
      INTEGER NXMAX                                                     
      DOUBLE PRECISION X(NXMAX,MX)                                      
                                                                        
      INTEGER I,J                                                       
                                                                        
      DO 10 I = 1,MX                                                    
       DO 10 J = 1,NEFF                                                 
        X(J,I) = BK(J,ORDER(I))                                         
10    CONTINUE                                                          
                                                                        
      DO 40 I = 1,MX                                                    
       IF (ORDER(I).LE.M) THEN                                          
        DO 20 J = 1,NL                                                  
         X(J + NEFF,I) = 0.D0                                           
20      CONTINUE                                                        
       ELSE                                                             
        DO 30 J = 1,NL                                                  
         X(J + NEFF,I) = LSCALE * L(J,ORDER(I) - M)                     
30      CONTINUE                                                        
       ENDIF                                                            
40    CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: RE(DUCE) E(FFECTIVE) N(UMBER OF) D(ATA) P(OINTS)          
*                                                                       
* PURPOSE   : REDUCE EFFECTIVE NUMBER OF DATA POINTS BY ORTHOGONAL      
*             TRANSFORMATION                                            
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             N      : NUMBER OF DATA POINTS                            
*             GSTI   : ORDINATE OF THE DATA POINTS                      
*                                                                       
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
*             BK     : MATRIX COMPOSED OUT OF B AND K                   
*             NBKMAX : LEADING DIMENSION OF THE FIELD BK                
*                                                                       
*             IFWK   : INTEGER WORKSPACE                                
*             DPFWK  : DOUBLE PRECISION WORKSPACE                       
*                                                                       
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
*             (ON OUTPUT)                                               
*             GSTI   : DATA POINTS GSTI MULTIPLIED WITH THE TRANSPOSE   
*                        OF AN ORTHOGONAL MATRIX                        
*                                                                       
*             BK     : MATRIX K MULTIPLIED WITH THE TRANSPOSE OF AN     
*                        ORTHOGONAL MATRIX                              
*                                                                       
*             NEFF   : EFFECTIVE NUMBER OF DATA POINTS                  
*                                                                       
* COMMENT   : THE SIZE OF THE INTEGER WORKSPACE MUST BE AT LEAST N      
*                                                                       
*             THE SIZE OF THE DOUBLE PRECISION WORKSPACE MUST BE AT     
*             LEAST (M + NS)                                            
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE REENDP(N,GSTI,M,NS,BK,NBKMAX,NEFF,IFWK,DPFWK,INFMOD)   
                                                                        
      INTEGER N                                                         
      DOUBLE PRECISION GSTI(N)                                          
                                                                        
      INTEGER M,NS,NBKMAX                                               
      DOUBLE PRECISION BK(NBKMAX,M + NS)                                
                                                                        
      INTEGER NEFF                                                      
                                                                        
      INTEGER IFWK(N)                                                   
      DOUBLE PRECISION DPFWK(M + NS)                                    
                                                                        
      INTEGER INFMOD                                                    
                                                                        
      INTEGER I,J                                                       
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' REENDP > CALLED'                                     
      ENDIF                                                             
                                                                        
*-----INITIALIZE VARIABLES--------------------------------------------- 
                                                                        
      NEFF = N                                                          
                                                                        
*-----REDUCE EFF. NUMBER OF DATA POINTS FOR N >= M + NS---------------- 
                                                                        
      IF (NEFF.GE.(M + NS)) THEN                                        
                                                                        
       CALL CAQRDE(N,M + NS,BK,NBKMAX,DPFWK,IFWK)                       
                                                                        
       CALL AHOTRV(N,M + NS,BK,NBKMAX,DPFWK,GSTI)                       
                                                                        
       DO 10 I = 1,(M + NS - 1)                                         
        DO 10 J = (I + 1),(M + NS)                                      
         BK(J,I) = 0.D0                                                 
10     CONTINUE                                                         
                                                                        
       NEFF = M + NS                                                    
                                                                        
      ENDIF                                                             
                                                                        
      IF (INFMOD.GE.3) THEN                                             
       WRITE(*,1001)'REENDP > RED. NUM. OF DATA POINTS    :',NEFF       
       DO 20 I = 1,NEFF                                                 
        WRITE(*,1002)'REENDP > TRANSFORMED DATA POINT      :',          
     *               I,GSTI(I)                                          
20     CONTINUE                                                         
      ENDIF                                                             
                                                                        
      RETURN                                                            
                                                                        
1001  FORMAT(1X,A38,I5)                                                 
1002  FORMAT(1X,A38,I5,D18.8)                                           
                                                                        
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: IN(ITIALIZE) MA(TRIX) DE(COMPOSITIONS)                    
*                                                                       
* PURPOSE   : INITIALIZE THE MATRIX DECOMPOSITION FOR THE CALCULATION   
*             OF THE REGULARIZATION PARAMETER / REGULARIZED SOLUTION    
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             NEFF   : EFFECTIVE NUMBER OF DATA POINTS                  
*             GSTI   : ORDINATE OF THE DATA POINTS                      
*                                                                       
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
*             BK     : MATRIX COMPOSED OUT OF B AND K                   
*             NBKMAX : LEADING DIMENSION OF THE FIELD BK                
*                                                                       
*             LSCALE : SCALEFACTOR FOR THE APP. REG. OPERATOR           
*             NL     : FIRST DIMENSION OF THE REG. OPERATOR             
*             L      : APPROXIMATED REGULARIZATION OPERATOR             
*             NLMAX  : LEADING DIMENSION OF THE FIELD L                 
*                                                                       
*             ORDER  : ORDER OF THE PARAMETER                           
*                                                                       
*             MX     : SECOND DIMENSION OF THE MATRIX X                 
*             NXMAX  : LEADING DIMENSION OF THE FIELD X                 
*                                                                       
*             MODUS  : MODUS OF INMADE                                  
*                                                                       
*             IFWK   : INTEGER WORKSPACE                                
*             DPFWK  : DOUBLE PRECISION WORKSPACE                       
*                                                                       
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
*             (ON OUTPUT)                                               
*             NX     : FIRST DIMENSION OF THE MATRIX X                  
*             X      : MATRIX FOR THE CALCULATION OF THE SOLUTION       
*                                                                       
*             W      : SINGULAR VALUES OF THE FIRST PART OF Q           
*             UTGS   : GSTI MULTIPLIED WITH THE LEFT SINGULAR VECTORS   
*                      OF THE FIRST PART OF Q                           
*                                                                       
* COMMENT   : THE LEADING DIMENSION NXMAX OF THE FIELD X MUST BE AT     
*             LEAST MAX(NEFF,MX)                                        
*                                                                       
*             THE SIZE OF THE INTEGER WORKSPACE MUST BE AT LEAST        
*             (NEFF + NL)                                               
*                                                                       
*             THE SIZE OF THE DOUBLE PRECISION WORKSPACE MUST BE AT     
*             LEAST (NEFF + NL) * MX + MIN(NEFF,MX)                     
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE INMADE(NEFF,GSTI,M,NS,BK,NBKMAX,LSCALE,NL,L,NLMAX,     
     *                  ORDER,NX,MX,X,NXMAX,W,UTGS,MODUS,IFWK,DPFWK,    
     *                  INFMOD)                                         
                                                                        
      INTEGER NEFF                                                      
      DOUBLE PRECISION GSTI(NEFF)                                       
                                                                        
      INTEGER M,NS,NBKMAX                                               
      DOUBLE PRECISION BK(NBKMAX,M + NS)                                
                                                                        
      INTEGER NL,NLMAX                                                  
      DOUBLE PRECISION LSCALE                                           
      DOUBLE PRECISION L(NLMAX,NS)                                      
                                                                        
      INTEGER ORDER(M + NS)                                             
                                                                        
      INTEGER NX,MX,NXMAX                                               
      DOUBLE PRECISION X(NXMAX,MX)                                      
                                                                        
*-----DOUBLE PRECISION W(MIN0(NEFF,MX))-------------------------------- 
      DOUBLE PRECISION W(MX)                                            
                                                                        
      DOUBLE PRECISION UTGS(NEFF)                                       
                                                                        
      LOGICAL MODUS                                                     
                                                                        
      INTEGER IFWK(NEFF + NL)                                           
*-----DOUBLE PRECISION DPFWK((NEFF + NL) * MX + MIN0(NEFF,MX))--------- 
      DOUBLE PRECISION DPFWK((NEFF + NL) * MX + MX)                     
                                                                        
      INTEGER INFMOD                                                    
                                                                        
      INTEGER I,J                                                       
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' INMADE > CALLED'                                     
      ENDIF                                                             
                                                                        
*-----CALCULATE PART OF QR-DECOMPOSITION------------------------------- 
                                                                        
      IF (MODUS) THEN                                                   
                                                                        
       CALL SEMAKL(NEFF,M,NS,BK,NBKMAX,LSCALE,NL,L,NLMAX,MX,ORDER,      
     *             DPFWK,NEFF + NL)                                     
                                                                        
       CALL CAQRDE(NEFF + NL,MX,DPFWK,NEFF + NL,X,IFWK)                 
                                                                        
       DO 10 I = 1,MX                                                   
        DO 10 J = 1,I                                                   
         X(J,I) = DPFWK(J + (I - 1) * (NEFF + NL))                      
10     CONTINUE                                                         
                                                                        
      ENDIF                                                             
                                                                        
      DO 20 I = 1,MX                                                    
       DO 20 J = 1,NEFF                                                 
        DPFWK(J + (I - 1) * NEFF) = BK(J,ORDER(I))                      
20    CONTINUE                                                          
                                                                        
      CALL SOLTME(MX,X,NXMAX,NEFF,DPFWK,NEFF)                           
                                                                        
*-----CALCULATE SVD-DECOMPOSITION OF Q--------------------------------- 
                                                                        
      DO 30 I = 1,NEFF                                                  
       UTGS(I) = GSTI(I)                                                
30    CONTINUE                                                          
        
                                                                                                                                
      CALL CASVDE(NEFF,MX,DPFWK,NEFF,W,UTGS,DPFWK(MX * NEFF + 1))       
                                                                        
      NX = MIN0(NEFF,MX)                                                
                                                                        
      IF (INFMOD.GE.3) THEN                                             
       DO 40 I = 1,NX                                                   
        WRITE(*,1001)'INMADE > SINGULAR VALUE              :',I,W(I)    
40     CONTINUE                                                         
       DO 50 I = 1,NX                                                   
        WRITE(*,1001)'INMADE > TRANSFORMED DATA POINT      :',          
     *               I,UTGS(I)                                          
50     CONTINUE                                                         
      ENDIF                                                             
                                                                        
*-----SOLVE EQUATION RX = V-------------------------------------------- 
                                                                        
      CALL SOUTME(MX,X,NXMAX,NX,DPFWK,NEFF)                             
                                                                        
      DO 60 I = 1,MX                                                    
       DO 60 J = 1,NX                                                   
        X(J,I) = DPFWK(J + (I - 1) * NEFF)                              
60    CONTINUE                                                          
                                                                        
      RETURN                                                            
                                                                        
1001  FORMAT(1X,A38,I5,D18.8)                                           
                                                                        
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: IN(VERT) RE(GULARIZED) I(NTEGRAL) O(PERATOR)              
*                                                                       
* PURPOSE   : CALCULATE THE INVERSE OF THE REGULARIZED INTEGRAL         
*             OPERATOR                                                  
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             NEFF   : EFFECTIVE NUMBER OF DATA POINTS                  
*                                                                       
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
*             BK     : MATRIX COMPOSED OUT OF B AND K                   
*             NBKMAX : LEADING DIMENSION OF THE FIELD BK                
*                                                                       
*             LSCALE : SCALEFACTOR FOR THE APP. REG. OPERATOR           
*             NL     : FIRST DIMENSION OF THE REG. OPERATOR             
*             L      : APPROXIMATED REGULARIZATION OPERATOR             
*             NLMAX  : LEADING DIMENSION OF THE FIELD L                 
*                                                                       
*             ORDER  : ORDER OF THE PARAMETER                           
*                                                                       
*             MX     : SECOND DIMENSION OF THE MATRIX X                 
*             NXMAX  : LEADING DIMENSION OF THE FIELD X                 
*                                                                       
*             MODUS  : MODUS OF INREIO                                 
*                                                                       
*             IFWK   : INTEGER WORKSPACE                                
*                                                                       
*             DPFWK  : DOUBLE PRECISION WORKSPACE                       
*                                                                       
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
*             (ON OUTPUT)                                               
*             NX     : FIRST DIMENSION OF THE MATRIX X                  
*             X      : MATRIX FOR THE CALCULATION OF THE SOLUTION       
*                                                                       
* COMMENT   : THE LEADING DIMENSION NXMAX OF THE FIELD X MUST BE AT     
*             LEAST MAX(NEFF,MX)                                        
*                                                                       
*             THE SIZE OF THE INTEGER WORKSPACE MUST BE AT LEAST        
*             (NEFF + NL)                                               
*                                                                       
*             THE SIZE OF THE DOUBLE PRECISION WORKSPACE MUST BE AT     
*             LEAST (NEFF + NL) * MX                                    
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE INREIO(NEFF,M,NS,BK,NBKMAX,LSCALE,NL,L,NLMAX,ORDER,    
     *                  NX,MX,X,NXMAX,MODUS,IFWK,DPFWK,INFMOD)          
                                                                        
      INTEGER NEFF                                                      
                                                                        
      INTEGER M,NS,NBKMAX                                               
      DOUBLE PRECISION BK(NBKMAX,M + NS)                                
                                                                        
      INTEGER NL,NLMAX                                                  
      DOUBLE PRECISION LSCALE                                           
      DOUBLE PRECISION L(NLMAX,NS)                                      
                                                                        
      INTEGER ORDER(MX)                                                 
                                                                        
      INTEGER NX,MX,NXMAX                                               
      DOUBLE PRECISION X(NXMAX,MX)                                      
                                                                        
      LOGICAL MODUS                                                     
                                                                        
      INTEGER IFWK(NEFF + NL)                                           
                                                                        
      DOUBLE PRECISION DPFWK(NEFF + NL,MX)                              
                                                                        
      INTEGER INFMOD                                                    
                                                                        
      INTEGER I,J                                                       
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' INREIO > CALLED'                                     
      ENDIF                                                             
                                                                        
*-----CALCULATE PART OF QR-DECOMPOSITION------------------------------- 
                                                                        
      IF (MODUS) THEN                                                   
                                                                        
       CALL SEMAKL(NEFF,M,NS,BK,NBKMAX,LSCALE,NL,L,NLMAX,MX,ORDER,      
     *             DPFWK,NEFF + NL)                                     
                                                                        
       CALL CAQRDE(NEFF + NL,MX,DPFWK,NEFF + NL,X,IFWK)                 
                                                                        
      ELSE                                                              
                                                                        
       DO 10 I = 1,MX                                                   
        DO 10 J = 1,I                                                   
         DPFWK(J,I) = X(J,I)                                            
10     CONTINUE                                                         
                                                                        
      ENDIF                                                             
                                                                        
      DO 20 I = 1,MX                                                    
       DO 20 J = 1,NEFF                                                 
        X(J,I) = BK(J,ORDER(I))                                         
20    CONTINUE                                                          
                                                                        
      CALL SOLTME(MX,DPFWK,NEFF + NL,NEFF,X,NXMAX)                      
                                                                        
      NX = NEFF                                                         
                                                                        
*-----SOLVE EQUATION RX = V-------------------------------------------- 
                                                                        
      CALL SOUTME(MX,DPFWK,NEFF + NL,NX,X,NXMAX)                        
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* FUNCTION  : CA(LCULATE) LAMB(DA)                                      
*                                                                       
* PURPOSE   : CALCULATE THE REGULARIZATION PARAMETER BY SEARCHING A     
*             ROOT IN THE FUNCTION CASCCR                               
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             N      : NUMBER OF DATA POINTS                            
*                                                                       
*             NX     : FIRST DIMENSION OF THE MATRIX X                  
*             MX     : SECOND DIMENSION OF THE MATRIX X                 
*             X      : MATRIX FOR THE CALCULATION OF THE SOLUTION       
*             NXMAX  : LEADING DIMENSION OF THE FIELD X                 
*                                                                       
*             W      : SINGULAR VALUES OF THE FIRST PART OF Q           
*                                                                       
*             UTGS   : GSTI MULTIPLIED WITH THE LEFT SINGULAR VECTORS   
*                      OF THE FIRST PART OF Q                           
*             GSTIRQ : QUANTITY FOR THE CALCULATION OF THE DISCREPANCY  
*                                                                       
*             ERROR  : SCALEFACTOR OF THE ERROR                         
*             SCLMOD : MODE FOR THE ESTIMATION OF THE ERROR             
*                                                                       
*             LAMB?? : PARAMETER FOR THE CALCULATION OF LAMBDA          
*                                                                       
*             DPFWK  : WORKSPACE                                        
*                                                                       
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
* COMMENT   : THE PARAMETER M,N,NX,MX,X,NXMAX,W,UTGS,ERROR,SCLMOD,DPFWK 
*             ARE ONLY USED BY THE FUNCTION CASCCR CALLED BY CALAMB     
*                                                                       
*             THE SIZE OF THE DOUBLE PRECISION WORKSPACE MUST BE AT     
*             LEAST 3 NX                                                
*---------------------------------------------------------------------- 
                                                                        
      DOUBLE PRECISION FUNCTION CALAMB(M,N,NX,MX,X,NXMAX,W,UTGS,GSTIRQ, 
     *                                 ERROR,SCLMOD,LAMBST,LAMBSP,      
     *                                 LAMBRA,LAMBPR,LAMBIT,DPFWK,      
     *                                 INFMOD)                          
                                                                        
      INTEGER M,N                                                       
                                                                        
      INTEGER NX,MX,NXMAX                                               
      DOUBLE PRECISION X(NXMAX,MX)                                      
                                                                        
      DOUBLE PRECISION W(NX)                                            
                                                                        
      DOUBLE PRECISION GSTIRQ                                           
      DOUBLE PRECISION UTGS(NX)                                         
                                                                        
      INTEGER SCLMOD                                                    
      DOUBLE PRECISION ERROR                                            
                                                                        
      INTEGER LAMBIT                                                    
      DOUBLE PRECISION LAMBST,LAMBSP,LAMBRA,LAMBPR                      
                                                                        
      DOUBLE PRECISION DPFWK(3 * NX)                                    
                                                                        
      INTEGER INFMOD                                                    
                                                                        
      DOUBLE PRECISION CASCCR                                           
                                                                        
      CHARACTER*60 ERRTXT                                               
                                                                        
      DOUBLE PRECISION LAMMIN,LAMMAX                                    
                                                                        
      INTEGER IWK1,IWK2                                                 
      DOUBLE PRECISION DPWK1,DPWK2,DPWK3                                
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' CALAMB > CALLED'                                     
      ENDIF                                                             
                                                                        
*-----GET INTERVAL FOR LAMBDA------------------------------------------ 
                                                                        
      IF (INFMOD.GE.3) THEN                                             
       WRITE(*,*)' CALAMB > GET INTERVAL FOR LAMBDA'                    
      ENDIF                                                             
                                                                        
      LAMMIN = LAMBST                                                   
      LAMMAX = LAMBST                                                   
                                                                        
      DPWK3 = 10.D0**LAMBST                                             
      DPWK1 = CASCCR(DPWK3,M,N,NX,MX,X,NXMAX,W,UTGS,GSTIRQ,ERROR,       
     *               SCLMOD,DPFWK)                                      
      IWK1 = 1                                                          
      IF (INFMOD.GE.3) THEN                                             
       WRITE(*,1001)'CALAMB > ',LAMBST,DPWK1                            
      ENDIF                                                             
                                                                        
      IF (DPWK1.GT.0.D0) THEN                                           
                                                                        
10     LAMMIN = LAMMIN - LAMBSP                                         
       IF (DABS(LAMBST - LAMMIN).GT.LAMBRA) THEN                        
        ERRTXT = ' CALAMB > ERROR DETECTED (1)'                         
        CALL WEERRM(ERRTXT,.TRUE.)                                      
       ENDIF                                                            
       DPWK2 = DPWK1                                                    
                                                                        
       DPWK3 = 10.D0**LAMMIN                                            
       DPWK1 = CASCCR(DPWK3,M,N,NX,MX,X,NXMAX,W,UTGS,GSTIRQ,ERROR,      
     *                SCLMOD,DPFWK)                                     
       IWK1 = IWK1 + 1                                                  
       IF (INFMOD.GE.3) THEN                                            
        WRITE(*,1001)'CALAMB > ',LAMMIN,DPWK1                           
       ENDIF                                                            
                                                                        
       IF (DPWK1.GT.0.D0) THEN                                          
        GOTO 10                                                         
       ENDIF                                                            
       LAMMAX = LAMMIN + LAMBSP                                         
                                                                        
      ELSE                                                              
                                                                        
       DPWK2 = DPWK1                                                    
20     LAMMAX = LAMMAX + LAMBSP                                         
       IF (DABS(LAMBST - LAMMAX).GT.LAMBRA) THEN                        
        ERRTXT = ' CALAMB > ERROR DETECTED (2)'                         
        CALL WEERRM(ERRTXT,.TRUE.)                                      
       ENDIF                                                            
       DPWK1 = DPWK2                                                    
                                                                        
       DPWK3 = 10.D0**LAMMAX                                            
       DPWK2 = CASCCR(DPWK3,M,N,NX,MX,X,NXMAX,W,UTGS,GSTIRQ,ERROR,      
     *                SCLMOD,DPFWK)                                     
       IWK1 = IWK1 + 1                                                  
       IF (INFMOD.GE.3) THEN                                            
        WRITE(*,1001)'CALAMB > ',LAMMAX,DPWK2                           
       ENDIF                                                            
                                                                        
       IF (DPWK2.LT.0.D0) THEN                                          
        GOTO 20                                                         
       ENDIF                                                            
       LAMMIN = LAMMAX - LAMBSP                                         
                                                                        
      ENDIF                                                             
                                                                        
*-----CALCULATE LAMBDA WITH SC-CRITERIUM------------------------------- 
                                                                        
      IF (INFMOD.GE.3) THEN                                             
       WRITE(*,*)' CALAMB > CALCULATE LAMBDA'                           
      ENDIF                                                             
                                                                        
      IWK2 = 1                                                          
                                                                        
30    CALAMB = (DPWK2 * LAMMIN - DPWK1 * LAMMAX) / (DPWK2 - DPWK1)      
      IF ((LAMMAX - CALAMB).LT.(LAMBPR / 2.D0)) THEN                    
       CALAMB = LAMMAX - LAMBPR / 2.D0                                  
      ELSE IF ((CALAMB - LAMMIN).LT.(LAMBPR / 2.D0)) THEN               
       CALAMB = LAMMIN + LAMBPR / 2.D0                                  
      ENDIF                                                             
      DPWK3 = 10.D0**CALAMB                                             
      DPWK3 = CASCCR(DPWK3,M,N,NX,MX,X,NXMAX,W,UTGS,GSTIRQ,ERROR,       
     *               SCLMOD,DPFWK)                                      
      IWK1 = IWK1 + 1                                                   
      IF (DPWK3.GT.0.D0) THEN                                           
       LAMMAX = (DPWK2 * CALAMB - DPWK3 * LAMMAX)                       
       IF ((DPWK3.GE.DPWK2).OR.                                         
     *     (LAMMAX.LE.(LAMMIN * (DPWK2 - DPWK3)))) THEN                 
        LAMMAX = (LAMMIN + CALAMB) / 2.D0                               
       ELSE                                                             
        LAMMAX = LAMMAX / (DPWK2 - DPWK3)                               
       ENDIF                                                            
       DPWK2 = 10.D0**LAMMAX                                            
       DPWK2 = CASCCR(DPWK2,M,N,NX,MX,X,NXMAX,W,UTGS,GSTIRQ,ERROR,      
     *                SCLMOD,DPFWK)                                     
       IWK1 = IWK1 + 1                                                  
       IF (DPWK2.LE.0.D0) THEN                                          
        LAMMIN = LAMMAX                                                 
        DPWK1 = DPWK2                                                   
        LAMMAX = CALAMB                                                 
        DPWK2 = DPWK3                                                   
       ENDIF                                                            
      ELSE                                                              
       LAMMIN = (DPWK1 * CALAMB - DPWK3 * LAMMIN)                       
       IF ((DPWK3.LE.DPWK1).OR.                                         
     *     (LAMMIN.LE.(LAMMAX * (DPWK1 - DPWK3)))) THEN                 
        LAMMIN = (LAMMAX + CALAMB) / 2.D0                               
       ELSE                                                             
        LAMMIN = LAMMIN / (DPWK1 - DPWK3)                               
       ENDIF                                                            
       DPWK1 = 10.D0**LAMMIN                                            
       DPWK1 = CASCCR(DPWK1,M,N,NX,MX,X,NXMAX,W,UTGS,GSTIRQ,ERROR,      
     *                SCLMOD,DPFWK)                                     
       IWK1 = IWK1 + 1                                                  
       IF (DPWK1.GT.0.D0) THEN                                          
        LAMMAX = LAMMIN                                                 
        DPWK2 = DPWK1                                                   
        LAMMIN = CALAMB                                                 
        DPWK1 = DPWK3                                                   
       ENDIF                                                            
      ENDIF                                                             
                                                                        
      IF (INFMOD.GE.3) THEN                                             
       WRITE(*,1002)'CALAMB > ',LAMMIN,LAMMAX,DPWK1,DPWK2               
      ENDIF                                                             
                                                                        
      IF ((LAMMAX - LAMMIN).GT.LAMBPR) THEN                             
       IWK2 = IWK2 + 1                                                  
       IF (IWK2.EQ.LAMBIT) THEN                                         
        ERRTXT = ' CALAMB > ERROR DETECTED (3)'                         
        CALL WEERRM(ERRTXT,.FALSE.)                                     
       ELSE                                                             
        GOTO 30                                                         
       ENDIF                                                            
      ENDIF                                                             
                                                                        
      IF (INFMOD.GE.3) THEN                                             
       WRITE(*,1003)'CALAMB > FUNCTION EVALUATIONS        :',IWK1       
       WRITE(*,1003)'CALAMB > ITERATIONS                  :',IWK2       
       DPWK1 = LAMMAX - LAMMIN                                          
       WRITE(*,1004)'CALAMB > PRECISION                   :',DPWK1      
       WRITE(*,1004)'CALAMB > PARAMETER                   :',CALAMB     
       WRITE(*,1004)'CALAMB > ESTIMATED ERROR             :',ERROR      
      ENDIF                                                             
                                                                        
      RETURN                                                            
                                                                        
1001  FORMAT(1X,A9,2D16.7)                                              
1002  FORMAT(1X,A9,4D16.7)                                              
1003  FORMAT(1X,A38,I5)                                                 
1004  FORMAT(1X,A38,D18.8)                                              
                                                                        
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* FUNCTION  : CA(LCULATE) SC CR(ITERIUM)                                
*                                                                       
* PURPOSE   : CALCULATE THE SC CRITERIUM FOR A GIVEN REGULARIZATION     
*             PARAMETER                                                 
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             LAMBDA : GIVEN REGULARIZATION PARAMETER                   
*                                                                       
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             N      : NUMBER OF DATA POINTS                            
*                                                                       
*             NX     : FIRST DIMENSION OF THE MATRIX X                  
*             MX     : SECOND DIMENSION OF THE MATRIX X                 
*             X      : MATRIX FOR THE CALCULATION OF THE SOLUTION       
*             NXMAX  : LEADING DIMENSION OF THE FIELD X                 
*                                                                       
*             W      : SINGULAR VALUES OF THE FIRST PART OF Q           
*                                                                       
*             UTGS   : GSTI MULTIPLIED WITH THE LEFT SINGULAR VECTORS   
*                      OF THE FIRST PART OF Q                           
*             GSTIRQ : QUANTITY FOR THE CALCULATION OF THE DISCREPANCY  
*                                                                       
*             ERROR  : SCALEFACTOR OF THE ERROR                         
*             SCLMOD : MODE FOR THE ESTIMATION OF THE ERROR             
*                                                                       
*             DPFWK  : WORKSPACE                                        
*                                                                       
* COMMENT   : THE SIZE OF THE WORKSPACE MUST BE AT LEAST 3 NX           
*---------------------------------------------------------------------- 
                                                                        
      DOUBLE PRECISION FUNCTION CASCCR(LAMBDA,M,N,NX,MX,X,NXMAX,W,UTGS, 
     *                                 GSTIRQ,ERROR,SCLMOD,DPFWK)       
                                                                        
      DOUBLE PRECISION LAMBDA                                           
                                                                        
      INTEGER M,N                                                       
                                                                        
      INTEGER NX,MX,NXMAX                                               
      DOUBLE PRECISION X(NXMAX,MX)                                      
                                                                        
      DOUBLE PRECISION W(NX)                                            
                                                                        
      DOUBLE PRECISION GSTIRQ                                           
      DOUBLE PRECISION UTGS(NX)                                         
                                                                        
      INTEGER SCLMOD                                                    
      DOUBLE PRECISION ERROR                                            
                                                                        
      DOUBLE PRECISION DPFWK(3 * NX)                                    
                                                                        
      DOUBLE PRECISION CASCER                                           
                                                                        
      INTEGER I,J                                                       
                                                                        
      INTEGER IWK1,IWK2,IWK3                                            
      DOUBLE PRECISION DPWK1,DPWK2,DPWK3,DPWK4,DPWK5,DPWK6              
                                                                        
      IWK1 = 0                                                          
      IWK2 = IWK1 + NX                                                  
      IWK3 = IWK2 + NX                                                  
      DO 10 I = 1,NX                                                    
       DPWK1 = W(I)                                                     
       DPWK2 = DPWK1 * DPWK1                                            
       DPWK3 = 1.D0 - DPWK2                                             
       DPWK4 = LAMBDA * DPWK3 + DPWK2                                   
       DPWK5 = DPWK4 * DPWK4                                            
       DPWK6 = DPWK5 * DPWK4                                            
       DPFWK(I + IWK1) = DPWK1 * DPWK3 / DPWK5 * UTGS(I)                
       DPFWK(I + IWK2) = DPWK1 * DPWK2 * DPWK3 / DPWK6 * UTGS(I)        
       DPFWK(I + IWK3) = DPWK2 * DPWK3 / DPWK6                          
10    CONTINUE                                                          
                                                                        
      DPWK1 = 0.D0                                                      
      DPWK2 = 0.D0                                                      
      DO 30 I = (M + 1),MX                                              
       DPWK3 = 0.D0                                                     
       DPWK4 = 0.D0                                                     
       DPWK5 = 0.D0                                                     
       DO 20 J = 1,NX                                                   
        DPWK3 = DPWK3 + X(J,I) * DPFWK(J + IWK1)                        
        DPWK4 = DPWK4 + X(J,I) * DPFWK(J + IWK2)                        
        DPWK5 = DPWK5 + X(J,I) * X(J,I) * DPFWK(J + IWK3)               
20     CONTINUE                                                         
       DPWK1 = DPWK1 + DPWK3 * DPWK4                                    
       DPWK2 = DPWK2 + DPWK5                                            
30    CONTINUE                                                          
                                                                        
      IF (SCLMOD.EQ.2) THEN                                             
       ERROR = CASCER(LAMBDA,N,NX,W,UTGS,GSTIRQ)                        
      ENDIF                                                             
                                                                        
      CASCCR = LAMBDA * DPWK1 / (ERROR * ERROR) / DPWK2 - 1.D0          
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* FUNCTION  : CA(LCULATE) SC ER(ROR)                                    
*                                                                       
* PURPOSE   : CALCULATE THE ERROR ESTIMATED WITH HELP OF THE SC-METHOD  
*             FOR A GIVEN REGULARIZATION PARAMETER                      
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             LAMBDA : GIVEN REGULARIZATION PARAMETER                   
*                                                                       
*             N      : NUMBER OF DATA POINTS                            
*             NX     : FIRST DIMENSION OF THE MATRIX X                  
*                                                                       
*             W      : SINGULAR VALUES OF THE FIRST PART OF Q           
*                                                                       
*             UTGS   : GSTI MULTIPLIED WITH THE LEFT SINGULAR VECTORS   
*                      OF THE FIRST PART OF Q                           
*             GSTIRQ : QUANTITY FOR THE CALCULATION OF THE DISCREPANCY  
*---------------------------------------------------------------------- 
                                                                        
      DOUBLE PRECISION FUNCTION CASCER(LAMBDA,N,NX,W,UTGS,GSTIRQ)       
                                                                        
      DOUBLE PRECISION LAMBDA                                           
                                                                        
      INTEGER N,NX                                                      
                                                                        
      DOUBLE PRECISION W(NX)                                            
                                                                        
      DOUBLE PRECISION GSTIRQ                                           
      DOUBLE PRECISION UTGS(NX)                                         
                                                                        
      INTEGER I                                                         
                                                                        
      DOUBLE PRECISION DPWK1,DPWK2,DPWK3,DPWK4,DPWK5,DPWK6,DPWK7        
                                                                        
      DPWK1 = GSTIRQ / (LAMBDA * LAMBDA)                                
      DPWK2 = 0.D0                                                      
      DPWK3 = (N - NX) / (LAMBDA * LAMBDA)                              
      DO 10 I = 1,NX                                                    
       DPWK4 = W(I) * W(I)                                              
       DPWK5 = 1.D0 - DPWK4                                             
       DPWK6 = LAMBDA * DPWK5 + DPWK4                                   
       DPWK7 = DPWK5 / DPWK6                                            
       DPWK3 = DPWK3 + DPWK7 * DPWK7                                    
       DPWK7 = DPWK7 * UTGS(I)                                          
       DPWK1 = DPWK1 + DPWK7 * DPWK7                                    
       DPWK7 = DPWK7 * DPWK4 / DPWK6                                    
       DPWK2 = DPWK2 + DPWK7 * DPWK7                                    
10    CONTINUE                                                          
                                                                        
      CASCER = DSQRT((DPWK1 - DPWK2) / DPWK3)                           
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: TE(ST) AC(TIVE) CO(NSTRAINTS)                             
*                                                                       
* PURPOSE   : TEST THE SET OF THE ACTIVE CONSTRAINTS BY CALCULATION OF  
*             THE LAGRANGE MULTIPLIERS                                  
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             NEFF   : EFFECTIVE NUMBER OF DATA POINTS                  
*             GSTI   : ORDINATE OF THE DATA POINTS                      
*                                                                       
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
*             BK     : MATRIX COMPOSED OUT OF B AND K                   
*             NBKMAX : LEADING DIMENSION OF THE FIELD BK                
*                                                                       
*             LSCALE : SCALEFACTOR FOR THE APP. REG. OPERATOR           
*             NL     : FIRST DIMENSION OF THE REG. OPERATOR             
*             L      : APPROXIMATED REGULARIZATION OPERATOR             
*             NLMAX  : LEADING DIMENSION OF THE FIELD L                 
*                                                                       
*             MX     : NUMBER OF FREE VARIABLES                         
*             REGSOL : REGULARIZED SOLUTION WITH MX FREE VARIABLES      
*             ORDER  : ORDER OF THE PARAMETER                           
*                                                                       
*             DPFWK1  : DOUBLE PRECISION WORKSPACE                      
*             DPFWK2  : DOUBLE PRECISION WORKSPACE                      
*                                                                       
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
*             (ON OUTPUT)                                               
*             FLAG   : FLAG FOR THE CORRECTNESS OF THE SET OF THE       
*                      ACTIVE CONSTRAINTS                               
*                                                                       
* COMMENT   : THE SIZE OF THE DOUBLE PRECISION WORKSPACE DPFWK1 MUST BE 
*             AT LEAST (NEFF + NL) * (M + NS)                           
*                                                                       
*           : THE SIZE OF THE DOUBLE PRECISION WORKSPACE DPFWK2 MUST BE 
*             AT LEAST (NEFF + NL)                                      
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE TEACCO(NEFF,GSTI,M,NS,BK,NBKMAX,LSCALE,NL,L,NLMAX,     
     *                  MX,REGSOL,ORDER,FLAG,DPFWK1,DPFWK2,INFMOD)      
                                                                        
      INTEGER NEFF                                                      
      DOUBLE PRECISION GSTI(NEFF)                                       
                                                                        
      INTEGER M,NS,NBKMAX                                               
      DOUBLE PRECISION BK(NBKMAX,M + NS)                                
                                                                        
      INTEGER NL,NLMAX                                                  
      DOUBLE PRECISION LSCALE                                           
      DOUBLE PRECISION L(NLMAX,NS)                                      
                                                                        
      INTEGER MX                                                        
                                                                        
      DOUBLE PRECISION REGSOL(MX)                                       
                                                                        
      INTEGER ORDER(M + NS)                                             
                                                                        
      LOGICAL FLAG                                                      
                                                                        
      DOUBLE PRECISION DPFWK1(NEFF + NL,M + NS)                         
      DOUBLE PRECISION DPFWK2(NEFF + NL)                                
                                                                        
      INTEGER INFMOD                                                    
                                                                        
      INTEGER I,J                                                       
                                                                        
      INTEGER IWK                                                       
      DOUBLE PRECISION DPWK                                             
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' TEACCO > CALLED'                                     
      ENDIF                                                             
                                                                        
*-----TEST INACTIVE CONSTRAINTS---------------------------------------- 
                                                                        
      FLAG = .FALSE.                                                    
      DO 10 I = 1,MX                                                    
       IF (INFMOD.GE.3) THEN                                            
        IWK = ORDER(I)                                                  
        DPWK = REGSOL(I)                                                
        WRITE(*,1001)'TEACCO > REGULARIZED SOLUTION        :',IWK,DPWK  
       ENDIF                                                            
       IF ((ORDER(I).GT.M).AND.(REGSOL(I).LT.0.D0)) THEN                
        FLAG = .TRUE.                                                   
        IF (INFMOD.LE.2) THEN                                           
         GOTO 50                                                        
        ENDIF                                                           
       ENDIF                                                            
10    CONTINUE                                                          
                                                                        
*-----TEST ACTIVE CONSTRAINTS WITH LAGRANGE MULTI.--------------------- 
                                                                        
      CALL SEMAKL(NEFF,M,NS,BK,NBKMAX,LSCALE,NL,L,NLMAX,M + NS,ORDER,   
     *            DPFWK1,NEFF + NL)                                     
                                                                        
      DO 20 I = 1,(NEFF + NL)                                           
       IF (I.LE.NEFF) THEN                                              
        DPFWK2(I) = -GSTI(I)                                            
       ELSE                                                             
        DPFWK2(I) = 0.D0                                                
       ENDIF                                                            
       DO 20 J = 1,MX                                                   
        DPFWK2(I) = DPFWK2(I) + DPFWK1(I,J) * REGSOL(J)                 
20    CONTINUE                                                          
                                                                        
      DO 40 I = (MX + 1),(M + NS)                                       
       DPWK = 0.D0                                                      
       DO 30 J = 1,(NEFF + NL)                                          
        DPWK = DPWK + DPFWK1(J,I) * DPFWK2(J)                           
30     CONTINUE                                                         
       IF (INFMOD.GE.3) THEN                                            
        IWK = ORDER(I)                                                  
        WRITE(*,1001)'TEACCO > LAGRANGE MULTIPLIER         :',IWK,DPWK  
       ENDIF                                                            
       IF (DPWK.LT.0.D0) THEN                                           
        FLAG = .TRUE.                                                   
        IF (INFMOD.LE.2) THEN                                           
         GOTO 50                                                        
        ENDIF                                                           
       ENDIF                                                            
40    CONTINUE                                                          
                                                                        
50    IF ((FLAG).AND.(INFMOD.GE.3)) THEN                                
       WRITE(*,*)' TEACCO > CHANGED SET OF ACTIVE CONSTRAINTS'          
      ENDIF                                                             
                                                                        
      RETURN                                                            
                                                                        
1001  FORMAT(1X,A38,I5,D18.8)                                           
                                                                        
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: IN(ITIALIZE) PA(RAMETERS) FOR Q(UADRATIC) P(ROGRAMMING)   
*                                                                       
* PURPOSE   : INITIALIZE THE SET OF ACTIVE CONSTRAINTS                  
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             MX     : TOTAL NUMBER OF PARAMETERS                       
*                                                                       
*             REGSOL : REGULARIZED SOLUTION (CALCULATED WITHOUT CONST.) 
*                                                                       
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
*             (ON OUTPUT)                                               
*             ACTCON : INITIAL GUESS OF THE SET OF ACTIVE CONSTRAINTS   
*                                                                       
*             FLAG   : FLAG FOR EXISTENCE OF ACTIVE CONSTRAINTS         
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE INPAQP(M,MX,REGSOL,ACTCON,FLAG,INFMOD)                 
                                                                        
      INTEGER M,MX                                                      
                                                                        
      DOUBLE PRECISION REGSOL(MX)                                       
                                                                        
      INTEGER ACTCON(MX)                                                
                                                                        
      LOGICAL FLAG                                                      
                                                                        
      INTEGER INFMOD                                                    
                                                                        
      INTEGER I                                                         
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' INPAQP > CALLED'                                     
      ENDIF                                                             
                                                                        
*-----CALCULATE INITIAL GUESS OF THE SET OF ACTIVE CONSTRAINTS--------- 
                                                                        
      DO 10 I = 1,M                                                     
       ACTCON(I) = 0                                                    
10    CONTINUE                                                          
                                                                        
      FLAG = .FALSE.                                                    
      DO 20 I = (M + 1),MX                                              
       ACTCON(I) = 0                                                    
       IF (REGSOL(I).LE.0.D0) THEN                                      
        ACTCON(I) = 1                                                   
        FLAG = .TRUE.                                                   
       ENDIF                                                            
20    CONTINUE                                                          
                                                                        
      IF ((.NOT.FLAG).AND.(INFMOD.GE.3)) THEN                           
       WRITE(*,*)' INPAQP > NO ACTIVE CONSTRAINTS'                      
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: CA(LCULATE) FS (WITH) Q(UADRATIC) P(ROGRAMMING)           
*                                                                       
* PURPOSE   : CALCULATE THE REGULARIZED SOLUTION WITH POSITIVITY        
*             CONSTRAINTS BY QUADRATIC PROGRAMMING                      
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             NEFF   : EFFECTIVE NUMBER OF DATA POINTS                  
*             GSTI   : ORDINATE OF THE DATA POINTS                      
*                                                                       
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
*             BK     : MATRIX COMPOSED OUT OF B AND K                   
*             NBKMAX : LEADING DIMENSION OF THE FIELD BK                
*                                                                       
*             LSCALE : SCALEFACTOR FOR THE APP. REG. OPERATOR           
*             NL     : FIRST DIMENSION OF THE REG. OPERATOR             
*             L      : APPROXIMATED REGULARIZATION OPERATOR             
*             NLMAX  : LEADING DIMENSION OF THE FIELD L                 
*                                                                       
*             MX     : NUMBER OF FREE VARIABLES                         
*                                                                       
*             NRMAX  : LEADING DIMENSION OF THE FIELD R                 
*                                                                       
*             IFWK   : INTEGER WORKSPACE                                
*                                                                       
*             DPFWK1 : DOUBLE PRECISION WORKSPACE                       
*             DPFWK2 : DOUBLE PRECISION WORKSPACE                       
*             DPFWK3 : DOUBLE PRECISION WORKSPACE                       
*                                                                       
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
*             (ON OUTPUT)                                               
*             ACTCON : FIELD WITH THE ACTIVE CONSTRAINTS                
*             ORDER  : ORDER OF THE PARAMETER                           
*                                                                       
*             R      : PART OF THE QR DECOMPOSITION OF THE COMPOSED     
*                      MATRIX                                           
*                                                                       
* COMMENT   : THE SIZE OF THE INTEGER WORKSPACE MUST BE AT LEAST        
*             (NEFF + NL)                                               
*                                                                       
*             THE SIZE OF THE DOUBLE PRECISION WORKSPACE DPFWK1 MUST BE 
*             AT LEAST (M + NS) * (NEFF + NL)                           
*                                                                       
*             THE SIZE OF THE DOUBLE PRECISION WORKSPACE DPFWK2 MUST BE 
*             AT LEAST (NEFF + NL)                                      
*                                                                       
*             THE SIZE OF THE DOUBLE PRECISION WORKSPACE DPFWK3 MUST BE 
*             AT LEAST (M + NS)                                         
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE CAFSQP(NEFF,GSTI,M,NS,BK,NBKMAX,LSCALE,NL,L,NLMAX,MX,  
     *                  ACTCON,ORDER,R,NRMAX,IFWK,DPFWK1,DPFWK2,DPFWK3, 
     *                  INFMOD)                                         
                                                                        
      INTEGER NEFF                                                      
      DOUBLE PRECISION GSTI(NEFF)                                       
                                                                        
      INTEGER M,NS,NBKMAX                                               
      DOUBLE PRECISION BK(NBKMAX,M + NS)                                
                                                                        
      INTEGER NL,NLMAX                                                  
      DOUBLE PRECISION LSCALE                                           
      DOUBLE PRECISION L(NLMAX,NS)                                      
                                                                        
      INTEGER MX                                                        
                                                                        
      INTEGER ACTCON(M + NS),ORDER(M + NS)                              
                                                                        
      INTEGER NRMAX                                                     
      DOUBLE PRECISION R(NRMAX,M + NS)                                  
                                                                        
      INTEGER IFWK(NEFF + NL)                                           
                                                                        
      DOUBLE PRECISION DPFWK1(NEFF + NL,M + NS)                         
      DOUBLE PRECISION DPFWK2(NEFF + NL)                                
      DOUBLE PRECISION DPFWK3(M + NS)                                   
                                                                        
      INTEGER INFMOD                                                    
                                                                        
      INTEGER I,J                                                       
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' CAFSQP > CALLED'                                     
      ENDIF                                                             
                                                                        
*-----SETUP VALUES FOR QUADRATIC MINIMIZATION-------------------------- 
                                                                        
      DO 10 I = 1,(M + NS)                                              
       ORDER(I) = I                                                     
10    CONTINUE                                                          
                                                                        
      CALL SEMAKL(NEFF,M,NS,BK,NBKMAX,LSCALE,NL,L,NLMAX,M + NS,ORDER,   
     *            DPFWK1,NEFF + NL)                                     
                                                                        
      DO 20 I = 1,NEFF                                                  
       DPFWK2(I) = GSTI(I)                                              
20    CONTINUE                                                          
                                                                        
      DO 30 I = (NEFF + 1),(NEFF + NL)                                  
       DPFWK2(I) = 0.D0                                                 
30    CONTINUE                                                          
                                                                        
*-----CALCULATE REGULARIZED SOLUTION BY QUADRATIC MINIMIZATION--------- 
                                                                        
      CALL MIQFWC(NEFF + NL,M + NS,NS,MX,DPFWK1,NEFF + NL,DPFWK2,       
     *            ACTCON,ORDER,IFWK,DPFWK3)                             
                                                                        
      IF (INFMOD.GE.3) THEN                                             
       WRITE(*,1001)'CAFSQP > NUMBER OF FREE PARAMETERS   :',MX         
       DO 40 I = (MX + 1),(M + NS)                                      
        WRITE(*,1001)'CAFSQP > ACTIVE CONSTRAINT           :',ORDER(I)  
40     CONTINUE                                                         
      ENDIF                                                             
                                                                        
      DO 50 I = 1,MX                                                    
       DO 50 J = 1,I                                                    
        R(J,I) = DPFWK1(J,I)                                            
50    CONTINUE                                                          
                                                                        
      RETURN                                                            
                                                                        
1001  FORMAT(1X,A38,I5)                                                 
                                                                        
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: CA(LCULATE) FS (WITH) M(ATRIX) D(ECOMPOSITION)            
*                                                                       
* PURPOSE   : CALCULATE FS WITH MATRIX DECOMPOSITION AS CALCULATED BY   
*             INMADE OR INREIO                                          
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             LAMBDA : REGULARIZATION PARAMETER                         
*                                                                       
*             NX     : FIRST DIMENSION OF THE MATRIX X                  
*             MX     : SECOND DIMENSION OF THE MATRIX X                 
*             X      : MATRIX FOR THE CALCULATION OF THE SOLUTION       
*             NXMAX  : LEADING DIMENSION OF THE FIELD X                 
*                                                                       
*             W      : SINGULAR VALUES OF THE FIRST PART OF Q           
*             UTGS   : GSTI MULTIPLIED WITH THE LEFT SINGULAR VECTORS   
*                      OF THE FIRST PART OF Q                           
*                                                                       
*             MODUS  : MODUS OF CAFSMD                                  
*                                                                       
*             DPFWK  : DOUBLE PRECISION WORKSPACE                       
*                                                                       
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
*             (ON OUTPUT)                                               
*             REGSOL : VALUES OF THE REGULARIZED SOLUTION               
*                                                                       
* COMMENT   : THE SIZE OF THE DOUBLE PRECISION  WORKSPACE MUST BE AT    
*             LEAST NX                                                  
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE CAFSMD(LAMBDA,NX,MX,X,NXMAX,W,UTGS,REGSOL,MODUS,DPFWK, 
     *                  INFMOD)                                         
                                                                        
      DOUBLE PRECISION LAMBDA                                           
                                                                        
      INTEGER NX,MX,NXMAX                                               
      DOUBLE PRECISION X(NXMAX,MX)                                      
                                                                        
      DOUBLE PRECISION W(NX)                                            
      DOUBLE PRECISION UTGS(NX)                                         
                                                                        
      DOUBLE PRECISION REGSOL(MX)                                       
                                                                        
      LOGICAL MODUS                                                     
                                                                        
      DOUBLE PRECISION DPFWK(NX)                                        
                                                                        
      INTEGER INFMOD                                                    
                                                                        
      INTEGER I,J                                                       
                                                                        
      DOUBLE PRECISION DPWK                                             
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' CAFSMD > CALLED'                                     
      ENDIF                                                             
                                                                        
      IF (MODUS) THEN                                                   
                                                                        
*-----CALCULATE REGULARIZED SOLUTION WITH MATRIX DECOMPOSITION--------- 
                                                                        
       DO 10 I = 1,NX                                                   
        DPFWK(I) = W(I) / (W(I)**2 + LAMBDA * (1.D0 - W(I)**2))         
        DPFWK(I) = DPFWK(I) * UTGS(I)                                   
10     CONTINUE                                                         
                                                                        
       DO 30 I = 1,MX                                                   
        DPWK = 0.D0                                                     
        DO 20 J = 1,NX                                                  
         DPWK = DPWK + X(J,I) * DPFWK(J)                                
20      CONTINUE                                                        
        REGSOL(I) = DPWK                                                
30     CONTINUE                                                         
                                                                        
      ELSE                                                              
                                                                        
*-----CALCULATE REGULARIZED SOLUTION WITH THE REG. INVERSE INT. OP.---- 
                                                                        
       DO 50 I = 1,MX                                                   
        DPWK = 0.D0                                                     
        DO 40 J = 1,NX                                                  
         DPWK = DPWK + X(J,I) * UTGS(J)                                 
40      CONTINUE                                                        
        REGSOL(I) = DPWK                                                
50     CONTINUE                                                         
                                                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: CA(LCULATE) ER(RORS WITH) M(ATRIX) D(ECOMPOSITION)        
*                                                                       
* PURPOSE   : CALCULATE ERRORS WITH MATRIX DECOMPOSITION AS CALCULATED  
*             BY INMADE OR INREIO                                       
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             LAMBDA : REGULARIZATION PARAMETER                         
*             ERROR  : SCALEFACTOR OF THE ERROR                         
*                                                                       
*             NX     : FIRST DIMENSION OF THE MATRIX X                  
*             MX     : SECOND DIMENSION OF THE MATRIX X                 
*             X      : MATRIX FOR THE CALCULATION OF THE SOLUTION       
*             NXMAX  : LEADING DIMENSION OF THE FIELD X                 
*                                                                       
*             W      : SINGULAR VALUES OF THE FIRST PART OF Q           
*                                                                       
*             MODUS  : MODUS OF CAERMD                                  
*                                                                       
*             DPFWK  : DOUBLE PRECISION WORKSPACE                       
*                                                                       
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
*             (ON OUTPUT)                                               
*             REGERR : ERRORS OF THE REGULARIZED SOLUTION               
*                                                                       
* COMMENT   : THE SIZE OF THE DOUBLE PRECISION WORKSPACE MUST BE AT     
*             LEAST NX                                                  
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE CAERMD(LAMBDA,ERROR,NX,MX,X,NXMAX,W,REGERR,MODUS,      
     *                  DPFWK,INFMOD)                                   
                                                                        
      DOUBLE PRECISION LAMBDA,ERROR                                     
                                                                        
      INTEGER NX,MX,NXMAX                                               
      DOUBLE PRECISION X(NXMAX,MX)                                      
                                                                        
      DOUBLE PRECISION W(NX)                                            
                                                                        
      DOUBLE PRECISION REGERR(MX)                                       
                                                                        
      LOGICAL MODUS                                                     
                                                                        
      DOUBLE PRECISION DPFWK(NX)                                        
                                                                        
      INTEGER INFMOD                                                    
                                                                        
      INTEGER I,J                                                       
                                                                        
      DOUBLE PRECISION DPWK                                             
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' CAERMD > CALLED'                                     
      ENDIF                                                             
                                                                        
      IF (MODUS) THEN                                                   
                                                                        
*-----CALCULATE ERRORS OF REG. SOL. WITH MATRIX DECOMPOSITION---------- 
                                                                        
       DO 10 I = 1,NX                                                   
        DPFWK(I) = W(I) / (W(I)**2 + LAMBDA * (1.D0 - W(I)**2))         
10     CONTINUE                                                         
                                                                        
       DO 30 I = 1,MX                                                   
        DPWK = 0.D0                                                     
        DO 20 J = 1,NX                                                  
         DPWK = DPWK + (X(J,I) * DPFWK(J))**2                           
20      CONTINUE                                                        
        REGERR(I) = ERROR * DSQRT(DPWK)                                 
30     CONTINUE                                                         
                                                                        
      ELSE                                                              
                                                                        
*-----CALCULATE ERRORS OF REG. SOL. WITH THE REG. INVERSE INT. OP.----- 
                                                                        
       DO 50 I = 1,MX                                                   
        DPWK = 0.D0                                                     
        DO 40 J = 1,NX                                                  
         DPWK = DPWK + (X(J,I))**2                                      
40      CONTINUE                                                        
        REGERR(I) = ERROR * DSQRT(DPWK)                                 
50     CONTINUE                                                         
                                                                        
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: RE(ORDER) RE(GULARIZED) SO(LUTION)                        
*                                                                       
* PURPOSE   : REORDER REGULARIZED SOLUTION AND COPY IT TO ASJ, AERRJ,   
*             FSSJ, FMSJ AND FERRSJ                                     
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             MX     : NUMBER OF FREE VARIABLES                         
*             ORDER  : ORDER OF THE VARIABLES                           
*             REGSOL : REGULARIZED SOLUTION WITH THE ORDER GIVEN BY     
*                      ORDER                                            
*             REGERR : ERRORS OF THE REGULARIZED SOLUTION WITH THE      
*                      ORDER GIVEN BY ORDER                             
*                                                                       
*             BSCALE : SCALEFACTOR FOR THE ADDITIONAL PARAMETERS        
*                                                                       
*             M      : NUMBER OF COEFFICIENTS AJ                        
*             NS     : NUMBER OF POINTS WHERE THE SOL. IS CALCULATED    
*                                                                       
*             POSMOD : MODE FOR THE CALCULATION OF A POSITIVE SOLUTION  
*                                                                       
*             INFMOD : MODE FOR OUTPUT OF ADDITIONAL INFORMATION        
*                                                                       
*             (ON OUTPUT)                                               
*             ASJ    : VALUES OF THE CALCULATED COEFFICIENTS AJ         
*             AERRJ  : ERRORS OF THE CALCULATED COEFFICIENTS AJ         
*                                                                       
*             FSSJ   : VALUES OF THE CALCULATED SOLUTION                
*             FMSJ   : MIDPOINT FOR THE ERRORS OF THE CALCULATED SOL.   
*             FERRSJ : ERRORS OF THE CALCULATED SOLUTION                
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE RERESO(MX,ORDER,REGSOL,REGERR,BSCALE,M,ASJ,AERRJ,NS,   
     *                  FSSJ,FMSJ,FERRSJ,POSMOD,INFMOD)                 
                                                                        
      INTEGER MX                                                        
      INTEGER ORDER(M + NS)                                             
      DOUBLE PRECISION REGSOL(MX),REGERR(MX)                            
                                                                        
      DOUBLE PRECISION BSCALE                                           
                                                                        
      INTEGER M                                                         
      DOUBLE PRECISION ASJ(M),AERRJ(M)                                  
                                                                        
      INTEGER NS                                                        
      DOUBLE PRECISION FSSJ(NS),FMSJ(NS),FERRSJ(NS)                     
                                                                        
      INTEGER POSMOD                                                    
                                                                        
      INTEGER INFMOD                                                    
                                                                        
      INTEGER I                                                         
                                                                        
*-----WRITE SUBROUTINE INFO-------------------------------------------- 
                                                                        
      IF (INFMOD.GE.2) THEN                                             
       WRITE(*,*)' RERESO > CALLED'                                     
      ENDIF                                                             
                                                                        
*-----REORDER REGULARIZED SOLUTION AND COPY IT TO OUTPUT VARIABLES----- 
                                                                        
      DO 10 I = 1,MX                                                    
       IF (ORDER(I).LE.M) THEN                                          
        ASJ(ORDER(I)) = REGSOL(I) * BSCALE                              
        AERRJ(ORDER(I)) = REGERR(I) * BSCALE                            
       ELSE                                                             
        FSSJ(ORDER(I) - M) = REGSOL(I)                                  
        FMSJ(ORDER(I) - M) = REGSOL(I)                                  
        FERRSJ(ORDER(I) - M) = REGERR(I)                                
       ENDIF                                                            
10    CONTINUE                                                          
                                                                        
*-----CORRECT SOLUTION IF POSITIVITY CONSTRAINTS ARE ACTIVE------------ 
                                                                        
      IF (POSMOD.NE.1) THEN                                             
       DO 20 I = (MX + 1),(M + NS)                                      
        FSSJ(ORDER(I) - M) = 0.D0                                       
        FMSJ(ORDER(I) - M) = 0.D0                                       
        FERRSJ(ORDER(I) - M) = 0.D0                                     
20     CONTINUE                                                         
       DO 30 I = 1,NS                                                   
        IF (FSSJ(I).LT.0.D0) THEN                                       
         IF (INFMOD.GE.3) THEN                                          
          WRITE(*,1001)'RERESO > NEGATIVE VALUE              :',        
     *                  I,FSSJ(I)                                       
         ENDIF                                                          
         FSSJ(I) = 0.D0                                                 
         FMSJ(I) = 0.D0                                                 
        ENDIF                                                           
        IF (FERRSJ(I).GT.FSSJ(I)) THEN                                  
         FMSJ(I) = (FERRSJ(I) + FSSJ(I)) / 2.D0                         
         FERRSJ(I) = FMSJ(I)                                            
        ENDIF                                                           
30     CONTINUE                                                         
      ENDIF                                                             
                                                                        
      RETURN                                                            
                                                                        
1001  FORMAT(1X,A38,I5,D18.8)                                           
                                                                        
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: SO(LVE) U(PPER) T(RIANGULAR) M(ATRIX) E(QUATION)          
*                                                                       
* PURPOSE   : SOLVE THE MATRIX EQUATION RXT = BT (R UPPER TRIANGULAR)   
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             NR     : FIRST / SECOND DIMENSION OF THE MATRIX R         
*             R      : UPPER TRIANGULAR MATRIX                          
*             NRMAX  : LEADING DIMENSION OF THE FIELD R                 
*                                                                       
*             NB     : FIRST DIMENSION OF THE MATRIX B                  
*             B      : MATRIX                                           
*             NBMAX  : LEADING DIMENSION OF THE FIELD B                 
*                                                                       
*             (ON OUTPUT)                                               
*             B      : SOLUTION X                                       
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE SOUTME(NR,R,NRMAX,NB,B,NBMAX)                          
                                                                        
      INTEGER NR,NRMAX                                                  
      DOUBLE PRECISION R(NRMAX,NR)                                      
                                                                        
      INTEGER NB,NBMAX                                                  
      DOUBLE PRECISION B(NBMAX,NR)                                      
                                                                        
      INTEGER I,J,K                                                     
                                                                        
*-----SOLVE UPPER TRIANGULAR MATRIX EQUATION--------------------------- 
                                                                        
      DO 20 I = NR,1,-1                                                 
       DO 10 J = 1,NB                                                   
        B(J,I) = B(J,I) / R(I,I)                                        
10     CONTINUE                                                         
       DO 20 J = (I - 1),1,-1                                           
        DO 20 K = 1,NB                                                  
         B(K,J) = B(K,J) - B(K,I) * R(J,I)                              
20    CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: SO(LVE) L(OWER) T(RIANGULAR) M(ATRIX) E(QUATION)          
*                                                                       
* PURPOSE   : SOLVE THE MATRIX EQUATION RTXT = BT (R UPPER TRIANGULAR)  
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             NR     : FIRST / SECOND DIMENSION OF THE MATRIX R         
*             R      : UPPER TRIANGULAR MATRIX                          
*             NRMAX  : LEADING DIMENSION OF THE FIELD R                 
*                                                                       
*             NB     : FIRST DIMENSION OF THE MATRIX B                  
*             B      : MATRIX                                           
*             NBMAX  : LEADING DIMENSION OF THE FIELD B                 
*                                                                       
*             (ON OUTPUT)                                               
*             B      : SOLUTION X                                       
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE SOLTME(NR,R,NRMAX,NB,B,NBMAX)                          
                                                                        
      INTEGER NR,NRMAX                                                  
      DOUBLE PRECISION R(NRMAX,NR)                                      
                                                                        
      INTEGER NB,NBMAX                                                  
      DOUBLE PRECISION B(NBMAX,NR)                                      
                                                                        
      INTEGER I,J,K                                                     
                                                                        
*-----SOLVE LOWER TRIANGULAR MATRIX EQUATION--------------------------- 
                                                                        
      DO 20 I = 1,NR                                                    
       DO 10 J = 1,NB                                                   
        B(J,I) = B(J,I) / R(I,I)                                        
10     CONTINUE                                                         
       DO 20 J = (I + 1),NR                                             
        DO 20 K = 1,NB                                                  
         B(K,J) = B(K,J) - B(K,I) * R(I,J)                              
20    CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: CA(LCULATE) QR DE(COMPOSITION)                            
*                                                                       
* PURPOSE   : CALCULATE R OF THE QR DECOMPOSITION OF THE MATRIX B       
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             NB     : FIRST DIMENSION OF THE MATRIX B                  
*             MB     : SECOND DIMENSION OF THE MATRIX B                 
*             B      : MATRIX                                           
*             NBMAX  : LEADING DIMENSION OF THE FIELD B                 
*                                                                       
*             IFWK   : INTEGER WORKSPACE                                
*                                                                       
*             (ON OUTPUT)                                               
*             B      : R (UPPER TRIANGULAR PART) AND INFORMATION ABOUT  
*                      THE ORTHOGONAL TRANSFORMATIONS (LOWER TRIANGULAR 
*                      PART)                                            
*             SCALE  : SCALING FACTOR FOR THE ORTHOGONAL TRANSFORM.     
*                                                                       
* COMMENT   : THE SIZE OF THE INTEGER WORKSPACE MUST BE AT LEAST NB     
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE CAQRDE(NB,MB,B,NBMAX,SCALE,IFWK)                       
                                                                        
      INTEGER NB,MB,NBMAX                                               
      DOUBLE PRECISION B(NBMAX,MB)                                      
                                                                        
      DOUBLE PRECISION SCALE(MB)                                        
                                                                        
      INTEGER IFWK(NB)                                                  
                                                                        
      INTEGER NBEFF                                                     
                                                                        
      INTEGER I,J,K                                                     
                                                                        
      DOUBLE PRECISION DPWK                                             
                                                                        
      DO 50 I = 1,MB                                                    
                                                                        
*-----CONSTRUCT HOUSEHOLDER TRANSFORMATION----------------------------- 
                                                                        
       DPWK = 0.D0                                                      
       DO 10 J = (I + 1),NB                                             
        DPWK = DPWK + B(J,I) * B(J,I)                                   
10     CONTINUE                                                         
       IF (DPWK.GT.0.D0) THEN                                           
        DPWK = -DSIGN(DSQRT(DPWK + B(I,I) * B(I,I)),B(I,I))             
        B(I,I) = B(I,I) - DPWK                                          
        SCALE(I) = DPWK                                                 
                                                                        
*------APPLY HOUSEHOLDER TRANSFORMATION-------------------------------- 
                                                                        
        NBEFF = I - 1                                                   
        DO 20 J = I,NB                                                  
         IF (B(J,I).NE.0.D0) THEN                                       
          NBEFF = NBEFF + 1                                             
          IFWK(NBEFF) = J                                               
         ENDIF                                                          
20      CONTINUE                                                        
                                                                        
        IF (NBEFF.GE.I) THEN                                            
         DO 40 J = (I + 1),MB                                           
          DPWK = 0.D0                                                   
          DO 30 K = I,NBEFF                                             
           DPWK = DPWK + B(IFWK(K),I) * B(IFWK(K),J)                    
30        CONTINUE                                                      
          DPWK = DPWK / B(I,I) / SCALE(I)                               
          DO 40 K = I,NBEFF                                             
           B(IFWK(K),J) = B(IFWK(K),J) + DPWK * B(IFWK(K),I)            
40       CONTINUE                                                       
        ENDIF                                                           
                                                                        
        DPWK = SCALE(I)                                                 
        SCALE(I) = B(I,I)                                               
        B(I,I) = DPWK                                                   
                                                                        
                                                                        
*------NO TRANSFORMATION NECESSARY------------------------------------- 
                                                                        
       ELSE                                                             
        SCALE(I) = 0.D0                                                 
       ENDIF                                                            
                                                                        
50    CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: A(PPLY) HO(USEHOLDER) TR(ANSFORMATIONS TO) V(ECTOR)       
*                                                                       
* PURPOSE   : MULTIPLY VECTOR WITH THE TRANSPOSE OF THE ORTHOGONAL      
*             MATRIX OBTAINED FROM A QR-FACTORIZATION                   
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             NB     : FIRST DIMENSION OF THE MATRIX B                  
*             MB     : SECOND DIMENSION OF THE MATRIX B                 
*             B      : MATRIX CONTAINING THE HOUSEHOLDER TRANSFORM.     
*                      BELOW THE DIAGONAL ELEMENTS                      
*             NBMAX  : LEADING DIMENSION OF THE FIELD B                 
*                                                                       
*             SCALE  : SCALEFACTORS FOR THE HOUSEHOLDER TRANSFORM.      
*                                                                       
*             X      : VECTOR                                           
*                                                                       
*             (ON OUTPUT)                                               
*             X      : VECTOR MULTIPLIED WITH QT                        
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE AHOTRV(NB,MB,B,NBMAX,SCALE,X)                          
                                                                        
      INTEGER NB,MB,NBMAX                                               
      DOUBLE PRECISION B(NBMAX,MB)                                      
                                                                        
      DOUBLE PRECISION SCALE(MB)                                        
                                                                        
      DOUBLE PRECISION X(NB)                                            
                                                                        
      INTEGER I,J                                                       
                                                                        
      DOUBLE PRECISION DPWK                                             
                                                                        
                                                                        
*-----APPLY HOUSEHOLDER TRANSFORMATIONS TO VECTOR---------------------- 
                                                                        
      DO 30 I = 1,MB                                                    
                                                                        
       IF (SCALE(I).NE.0.D0) THEN                                       
        DPWK = SCALE(I) * X(I)                                          
        DO 10 J = (I + 1),NB                                            
         DPWK = DPWK + B(J,I) * X(J)                                    
10      CONTINUE                                                        
        X(I) = X(I) + DPWK / B(I,I)                                     
        DPWK = DPWK / SCALE(I) / B(I,I)                                 
        DO 20 J = (I + 1),NB                                            
         X(J) = X(J) + DPWK * B(J,I)                                    
20      CONTINUE                                                        
       ENDIF                                                            
                                                                        
30    CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: CA(LCULATE) S(INGULAR) V(ALUE) DE(COMPOSITION)            
*                                                                       
* PURPOSE   : CALCULATE VT,W AND UTX AS GIVEN BY A SINGULAR VALUE       
*             DECOMPOSITION                                             
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             NB     : FIRST DIMENSION OF THE MATRIX B                  
*             MB     : SECOND DIMENSION OF THE MATRIX B                 
*             B      : MATRIX                                           
*             NBMAX  : LEADING DIMENSION OF THE FIELD B                 
*                                                                       
*             X      : VECTOR                                           
*                                                                       
*             DPFWK  : DOUBLE PRECISION WORKSPACE                       
*                                                                       
*             (ON OUTPUT)                                               
*             B      : RIGHT HAND SINGULAR VALUES                       
*                                                                       
*             W      : SINGULAR VALUES                                  
*                                                                       
*             X      : VECTOR MULTIPLIED WITH UT                        
*                                                                       
* COMMENT   : THE SIZE OF THE FIELD W MUST BE AT LEAST MIN(NB,MB)       
*                                                                       
*             THE SIZE OF THE DOUBLE PRECISION WORKSPACE MUST BE AT     
*             LEAST MIN(NB,MB)                                          
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE CASVDE(NB,MB,B,NBMAX,W,X,DPFWK)                        
                                                                        
      INTEGER NB,MB,NBMAX                                               
      DOUBLE PRECISION B(NBMAX,MB)                                      
                                                                        
*-----DOUBLE PRECISION W(MIN0(NB,MB))---------------------------------- 
      DOUBLE PRECISION W(MB)                                            
                                                                        
      DOUBLE PRECISION X(NB)                                            
                                                                        
*-----DOUBLE PRECISION DPFWK(MIN0(NB,MB))------------------------------ 
      DOUBLE PRECISION DPFWK(MB)                                        
                                                                        
      CHARACTER*60 ERRTXT                                               
                                                                        
      DOUBLE PRECISION NORMB                                            
                                                                        
      INTEGER I,J,K,L                                                   
                                                                        
      INTEGER IWK                                                       
      DOUBLE PRECISION DPWK1,DPWK2,DPWK3,DPWK4,DPWK5,DPWK6              
                                                                        
      NORMB = 0.D0                                                      
      DPFWK(1) = 0.D0                                                   
      IF (NB.GE.MB) THEN                                                
                                                                        
*-----REDUCE B TO UPPER BIDIAGONAL MATRIX------------------------------ 
                                                                        
       DO 70 I = 1,MB                                                   
                                                                        
*-----CONSTRUCT HOUSEHOLDER TRANSFORMATION (LEFT)---------------------- 
                                                                        
        DPWK1 = 0.D0                                                    
        DO 10 J = (I + 1),NB                                            
         DPWK1 = DPWK1 + B(J,I) * B(J,I)                                
10      CONTINUE                                                        
        IF (DPWK1.GT.0.D0) THEN                                         
         DPWK1 = -DSIGN(DSQRT(DPWK1 + B(I,I) * B(I,I)),B(I,I))          
         B(I,I) = B(I,I) - DPWK1                                        
         W(I) = DPWK1                                                   
                                                                        
*-----APPLY HOUSEHOLDER TRANSFORMATION (LEFT)-------------------------- 
                                                                        
         DO 30 J = (I + 1),MB                                           
          DPWK1 = 0.D0                                                  
          DO 20 K = I,NB                                                
           DPWK1 = DPWK1 + B(K,J) * B(K,I)                              
20        CONTINUE                                                      
          DPWK1 = DPWK1 / B(I,I) / W(I)                                 
          DO 30 K = I,NB                                                
           B(K,J) = B(K,J) + DPWK1 * B(K,I)                             
30       CONTINUE                                                       
                                                                        
*-----NO TRANSFORMATION NECESSARY-------------------------------------- 
                                                                        
        ELSE                                                            
         W(I) = B(I,I)                                                  
         B(I,I) = 0.D0                                                  
        ENDIF                                                           
                                                                        
        IF ((I + 1).LT.MB) THEN                                         
                                                                        
*-----CONSTRUCT HOUSEHOLDER TRANSFORMATION (RIGHT)--------------------- 
                                                                        
         DPWK1 = 0.D0                                                   
         DO 40 J = (I + 2),MB                                           
          DPWK1 = DPWK1 + B(I,J) * B(I,J)                               
40       CONTINUE                                                       
         IF (DPWK1.GT.0.D0) THEN                                        
          DPWK1 = -DSIGN(DSQRT(DPWK1 + B(I,I + 1) * B(I,I + 1)),        
     *                   B(I,I + 1))                                    
          B(I,I + 1) = B(I,I + 1) - DPWK1                               
          DPFWK(I + 1) = DPWK1                                          
                                                                        
*-----APPLY HOUSEHOLDER TRANSFORMATION (RIGHT)------------------------- 
                                                                        
          DO 60 J = (I + 1),NB                                          
           DPWK1 = 0.D0                                                 
           DO 50 K = (I + 1),MB                                         
            DPWK1 = DPWK1 + B(J,K) * B(I,K)                             
50         CONTINUE                                                     
           DPWK1 = DPWK1 / B(I,I + 1) / DPFWK(I + 1)                    
           DO 60 K = (I + 1),MB                                         
            B(J,K) = B(J,K) + DPWK1 * B(I,K)                            
60        CONTINUE                                                      
                                                                        
*-----NO TRANSFORMATION NECESSARY-------------------------------------- 
                                                                        
         ELSE                                                           
          DPFWK(I + 1) = B(I,I + 1)                                     
          B(I,I + 1) = 0.D0                                             
         ENDIF                                                          
                                                                        
        ELSE IF (I.LT.MB) THEN                                          
         DPFWK(MB) = B(I,I + 1)                                         
        ENDIF                                                           
                                                                        
*-----CALCULATE NORM OF B---------------------------------------------- 
                                                                        
        NORMB = DMAX1(NORMB,DABS(W(I)) + DABS(DPFWK(I)))                
                                                                        
70     CONTINUE                                                         
                                                                        
*-----APPLY TRANSFORMATIONS TO VECTOR---------------------------------- 
                                                                        
       DO 100 I = 1,MB                                                  
        IF (B(I,I).NE.0.D0) THEN                                        
         DPWK1 = 0.D0                                                   
         DO 80 J = I,NB                                                 
          DPWK1 = DPWK1 + B(J,I) * X(J)                                 
80       CONTINUE                                                       
         DPWK1 = DPWK1 / B(I,I) / W(I)                                  
         DO 90 J = I,NB                                                 
          X(J) = X(J) + DPWK1 * B(J,I)                                  
90       CONTINUE                                                       
        ENDIF                                                           
100    CONTINUE                                                         
                                                                        
*------CALCULATE THE TRANSPOSE OF V----------------------------------   
                                                                        
       B(MB,MB) = 1.D0                                                  
       DO 140 I = (MB - 2),1,-1                                         
        B(I + 1,I + 1) = 1.D0                                           
        DO 110 J = (I + 2),MB                                           
         B(I + 1,J) = 0.D0                                              
         B(J,I + 1) = 0.D0                                              
110     CONTINUE                                                        
        IF (B(I,I + 1).NE.0.D0) THEN                                    
         DO 130 J = (I + 1),MB                                          
          DPWK1 = 0.D0                                                  
          DO 120 K = (I + 1),MB                                         
           DPWK1 = DPWK1 + B(J,K) * B(I,K)                              
120       CONTINUE                                                      
          DPWK1 = DPWK1 / B(I,I + 1) / DPFWK(I + 1)                     
          DO 130 K = (I + 1),MB                                         
           B(J,K) = B(J,K) + DPWK1 * B(I,K)                             
130      CONTINUE                                                       
        ENDIF                                                           
140    CONTINUE                                                         
       B(1,1) = 1.D0                                                    
       DO 150 I = 2,MB                                                  
        B(I,1) = 0.D0                                                   
        B(1,I) = 0.D0                                                   
150    CONTINUE                                                         
                                                                        
*-----TRANSFORM UPPER BIDIAGONAL MATRIX TO DIAGONAL MATRIX------------- 
                                                                        
       DO 250 I = MB,1,-1  
	                                            
        DO 220 J = 1,500                                                
                                                                        
*-----SPLIT MATRIX---------------------------------------------------   
                                                                        
         DO 160 K = I,1,-1                                              
          IF ((DABS(DPFWK(K)) + NORMB).LE.NORMB) THEN                   
           IWK = K                                                      
           GOTO 190                                                     
          ELSE IF ((DABS(W(K - 1)) + NORMB).LE.NORMB) THEN              
           IWK = K                                                      
           GOTO 170                                                     
          ENDIF                                                         
160      CONTINUE                                                       
                                                                        
*-----APPLY GIVENS ROTATIONS TO DELETE DPFWK(K)----------------------   
                                                                        
170      DPWK1 = 0.D0                                                   
         DPWK2 = 1.D0                                                   
         DO 180 K = IWK,I                                               
          DPWK2 = DPWK2 * DPFWK(K)                                      
          DPFWK(K) = DPWK1 * DPFWK(K)                                   
          IF ((DABS(DPWK2) + NORMB).LE.NORMB) THEN                      
           GOTO 190                                                     
          ENDIF                                                         
          DPWK1 = W(K)                                                  
          W(K) = DSQRT(DPWK1 * DPWK1 + DPWK2 * DPWK2)                   
          DPWK1 = DPWK1 / W(K)                                          
          DPWK2 = -DPWK2 / W(K)                                         
          DPWK3 = X(IWK - 1)                                            
          DPWK4 = X(K)                                                  
          X(IWK - 1) = DPWK1 * DPWK3 + DPWK2 * DPWK4                    
          X(K) = -DPWK2 * DPWK3 + DPWK1 * DPWK4                         
180      CONTINUE                                                       
                                                                        
*-----TEST FOR CONVERGENCE-------------------------------------------   
                                                                        
190      IF (IWK.EQ.I) THEN                                             
          GOTO 230                                                      
         ENDIF                                                          
                                                                        
*-----CALCULATE PARAMETER FOR QR TRANSFORMATION----------------------   
                                                                        
         DPWK1 = (W(I - 1) - W(I)) * (W(I - 1) + W(I))                  
         DPWK2 = (DPFWK(I-1) - DPFWK(I)) * (DPFWK(I-1) + DPFWK(I))      
         DPWK3 = (W(IWK) - W(I)) * (W(IWK) + W(I))                      
         DPWK1 = (DPWK1 + DPWK2) / (2.D0 * DPFWK(I) * W(I - 1))         
         DPWK2 = DSQRT(DPWK1 * DPWK1 + 1.D0)                            
         DPWK2 = W(I - 1) / (DPWK1 + DSIGN(DPWK2,DPWK1)) - DPFWK(I)     
         DPWK3 = (DPWK3 + DPFWK(I) * DPWK2) / W(IWK)                    
         DPWK6 = W(IWK)                                                 
                                                                        
*-----PERFORM QR TRANSFORMATION--------------------------------------   
                                                                        
         DPWK1 = 1.D0                                                   
         DPWK2 = 1.D0                                                   
         DO 210 K = IWK,(I - 1)                                         
          DPWK5 = DPWK2 * DPFWK(K + 1)                                  
          DPWK4 = DPWK1 * DPFWK(K + 1)                                  
                                                                        
*-----CONSTRUCT GIVENS ROTATION (RIGHT)-------------------------------- 
                                                                        
          DPFWK(K) = DSQRT(DPWK3 * DPWK3 + DPWK5 * DPWK5)               
          IF (DPFWK(K).GT.0.D0) THEN                                    
           DPWK1 = DPWK3 / DPFWK(K)                                     
           DPWK2 = DPWK5 / DPFWK(K)                                     
          ELSE                                                          
           DPWK1 = 1.D0                                                 
           DPWK2 = 0.D0                                                 
          ENDIF                                                         
                                                                        
*-----APPLY GIVENS ROTATION (RIGHT) TO VT------------------------------ 
                                                                        
          DO 200 L = 1,MB                                               
           DPWK3 = B(K,L)                                               
           DPWK5 = B(K + 1,L)                                           
           B(K,L) = DPWK3 * DPWK1 + DPWK5 * DPWK2                       
           B(K + 1,L) = -DPWK3 * DPWK2 + DPWK5 * DPWK1                  
200       CONTINUE                                                      
                                                                        
*-----APPLY GIVENS ROTATION (RIGHT) TO BIDIAGONAL MATRIX--------------- 
                                                                        
          DPWK3 = DPWK6 * DPWK1 + DPWK4 * DPWK2                         
          DPWK4 = -DPWK6 * DPWK2 + DPWK4 * DPWK1                        
          DPWK5 = W(K + 1) * DPWK2                                      
          DPWK6 = W(K + 1) * DPWK1                                      
                                                                        
*-----CONSTRUCT GIVENS ROTATION (LEFT)--------------------------------- 
                                                                        
          W(K) = DSQRT(DPWK3 * DPWK3 + DPWK5 * DPWK5)                   
          IF (W(K).GT.0.D0) THEN                                        
           DPWK1 = DPWK3 / W(K)                                         
           DPWK2 = DPWK5 / W(K)                                         
          ELSE                                                          
           DPWK1 = 1.D0                                                 
           DPWK2 = 0.D0                                                 
          ENDIF                                                         
                                                                        
*-----APPLY GIVENS ROTATION (LEFT) TO X-------------------------------- 
                                                                        
          DPWK3 = X(K)                                                  
          DPWK5 = X(K + 1)                                              
          X(K) = DPWK1 * DPWK3 + DPWK2 * DPWK5                          
          X(K + 1) = -DPWK2 * DPWK3 + DPWK1 * DPWK5                     
                                                                        
*-----APPLY GIVENS ROTATION (LEFT) TO BIDIAGONAL MATRIX---------------- 
                                                                        
          DPWK3 = DPWK1 * DPWK4 + DPWK2 * DPWK6                         
          DPWK6 = -DPWK2 * DPWK4 + DPWK1 * DPWK6                        
210      CONTINUE                                                       
                                                                        
         DPFWK(IWK) = 0.D0                                              
         DPFWK(I) = DPWK3                                               
         W(I) = DPWK6                                                   
                                                                       
220     CONTINUE                                                        
        ERRTXT = ' CASVDE > ERROR DETECTED (1)'                         
        CALL WEERRM(ERRTXT,.TRUE.)                                      
                                                                        
*-----MAKE SINGULAR VALUE POSITIV-------------------------------------- 
                                                                        
230     IF (W(I).LT.0.D0) THEN                                          
         W(I) = -W(I)                                                   
         DO 240 J = 1,MB                                                
          B(I,J) = -B(I,J)                                              
240      CONTINUE                                                       
        ENDIF                                                           
250    CONTINUE                                                         
                                                                        
      ELSE                                                              
                                                                        
*-----REDUCE B TO LOWER BIDIAGONAL MATRIX------------------------------ 
                                                                        
       DO 320 I = 1,NB                                                  
                                                                        
*-----CONSTRUCT HOUSEHOLDER TRANSFORMATION (RIGHT)--------------------- 
                                                                        
        DPWK1 = 0.D0                                                    
        DO 260 J = (I + 1),MB                                           
         DPWK1 = DPWK1 + B(I,J) * B(I,J)                                
260     CONTINUE                                                        
        IF (DPWK1.GT.0.D0) THEN                                         
         DPWK1 = -DSIGN(DSQRT(DPWK1 + B(I,I) * B(I,I)),B(I,I))          
         B(I,I) = B(I,I) - DPWK1                                        
         W(I) = DPWK1                                                   
                                                                        
*-----APPLY HOUSEHOLDER TRANSFORMATION (RIGHT)------------------------- 
                                                                        
         DO 280 J = (I + 1),NB                                          
          DPWK1 = 0.D0                                                  
          DO 270 K = I,MB                                               
           DPWK1 = DPWK1 + B(J,K) * B(I,K)                              
270       CONTINUE                                                      
          DPWK1 = DPWK1 / B(I,I) / W(I)                                 
          DO 280 K = I,MB                                               
           B(J,K) = B(J,K) + DPWK1 * B(I,K)                             
280      CONTINUE                                                       
                                                                        
*-----NO TRANSFORMATION NECESSARY-------------------------------------- 
                                                                        
        ELSE                                                            
         W(I) = B(I,I)                                                  
         B(I,I) = 0.D0                                                  
        ENDIF                                                           
                                                                        
        IF ((I + 1).LT.NB) THEN                                         
                                                                        
*-----CONSTRUCT HOUSEHOLDER TRANSFORMATION (LEFT)---------------------- 
                                                                        
         DPWK1 = 0.D0                                                   
         DO 290 J = (I + 2),NB                                          
          DPWK1 = DPWK1 + B(J,I) * B(J,I)                               
290      CONTINUE                                                       
         IF (DPWK1.GT.0.D0) THEN                                        
          DPWK1 = -DSIGN(DSQRT(DPWK1 + B(I + 1,I) * B(I + 1,I)),        
     *                   B(I + 1,I))                                    
          B(I + 1,I) = B(I + 1,I) - DPWK1                               
          DPFWK(I + 1) = DPWK1                                          
                                                                        
*-----APPLY HOUSEHOLDER TRANSFORMATION (LEFT)-------------------------- 
                                                                        
          DO 310 J = (I + 1),MB                                         
           DPWK1 = 0.D0                                                 
           DO 300 K = (I + 1),NB                                        
            DPWK1 = DPWK1 + B(K,J) * B(K,I)                             
300        CONTINUE                                                     
           DPWK1 = DPWK1 / B(I + 1,I) / DPFWK(I + 1)                    
           DO 310 K = (I + 1),NB                                        
            B(K,J) = B(K,J) + DPWK1 * B(K,I)                            
310       CONTINUE                                                      
                                                                        
*-----NO TRANSFORMATION NECESSARY-------------------------------------- 
                                                                        
         ELSE                                                           
          DPFWK(I + 1) = B(I + 1,I)                                     
          B(I + 1,I) = 0.D0                                             
         ENDIF                                                          
                                                                        
        ELSE IF (I.LT.NB) THEN                                          
         DPFWK(NB) = B(I + 1,I)                                         
        ENDIF                                                           
                                                                        
*-----CALCULATE NORM OF B---------------------------------------------- 
                                                                        
        NORMB = DMAX1(NORMB,DABS(W(I)) + DABS(DPFWK(I)))                
                                                                        
320    CONTINUE                                                         
                                                                        
*-----APPLY TRANSFORMATIONS TO VECTOR---------------------------------- 
                                                                        
       DO 350 I = 1,(NB - 2)                                            
        IF (B(I + 1,I).NE.0.D0) THEN                                    
         DPWK1 = 0.D0                                                   
         DO 330 J = (I + 1),NB                                          
          DPWK1 = DPWK1 + B(J,I) * X(J)                                 
330      CONTINUE                                                       
         DPWK1 = DPWK1 / B(I + 1,I) / DPFWK(I + 1)                      
         DO 340 J = (I + 1),NB                                          
          X(J) = X(J) + DPWK1 * B(J,I)                                  
340      CONTINUE                                                       
        ENDIF                                                           
350    CONTINUE                                                         
                                                                        
*------CALCULATE THE TRANSPOSE OF V----------------------------------   
                                                                        
       DO 400 I = NB,1,-1                                               
        IF (B(I,I).NE.0.D0) THEN                                        
         DO 370 J = (I + 1),NB                                          
          DPWK1 = 0.D0                                                  
          DO 360 K = (I + 1),MB                                         
           DPWK1 = DPWK1 + B(J,K) * B(I,K)                              
360       CONTINUE                                                      
          DPWK1 = DPWK1 / W(I)                                          
          B(J,I) = DPWK1                                                
          DPWK1 = DPWK1 / B(I,I)                                        
          DO 370 K = (I + 1),MB                                         
           B(J,K) = B(J,K) + DPWK1 * B(I,K)                             
370      CONTINUE                                                       
         DO 380 K = (I + 1),MB                                          
          B(I,K) = B(I,K) / W(I)                                        
380      CONTINUE                                                       
         B(I,I) = 1.D0 + B(I,I) / W(I)                                  
        ELSE                                                            
         DO 390 J = (I + 1),NB                                          
          B(J,I) = 0.D0                                                 
390      CONTINUE                                                       
         B(I,I) = 1.D0                                                  
        ENDIF                                                           
400    CONTINUE                                                         
                                                                        
*-----TRANSFORM LOWER BIDIAGONAL MATRIX TO DIAGONAL MATRIX------------- 

      DO 500 I = NB,1,-1                                                
       DO 470 J = 1,500                                                 
                                                                        
*-----SPLIT MATRIX---------------------------------------------------   
                                                                        
        DO 410 K = I,1,-1                                               
         IF ((DABS(DPFWK(K)) + NORMB).LE.NORMB) THEN                    
          IWK = K                                                       
          GOTO 440                                                      
         ELSE IF ((DABS(W(K - 1)) + NORMB).LE.NORMB) THEN               
          IWK = K                                                       
          GOTO 420                                                      
         ENDIF                                                          
410     CONTINUE                                                        
                                                                        
*-----APPLY GIVENS ROTATIONS TO DELETE DPFWK(K)----------------------   
                                                                        
420     DPWK1 = 0.D0                                                    
        DPWK2 = 1.D0                                                    
        DO 430 K = IWK,I                                                
         DPWK2 = DPWK2 * DPFWK(K)                                       
         DPFWK(K) = DPWK1 * DPFWK(K)                                    
         IF ((DABS(DPWK2) + NORMB).LE.NORMB) THEN                       
          GOTO 440                                                      
         ENDIF                                                          
         DPWK1 = W(K)                                                   
         W(K) = DSQRT(DPWK1 * DPWK1 + DPWK2 * DPWK2)                    
         DPWK1 = DPWK1 / W(K)                                           
         DPWK2 = -DPWK2 / W(K)                                          
         DO 430 L = 1,MB                                                
          DPWK3 = B(IWK - 1,L)                                          
          DPWK4 = B(K,L)                                                
          B(IWK - 1,L) = DPWK1 * DPWK3 + DPWK2 * DPWK4                  
          B(K,L) = -DPWK2 * DPWK3 + DPWK1 * DPWK4                       
430     CONTINUE                                                        
                                                                        
*-----TEST FOR CONVERGENCE-------------------------------------------   
                                                                        
440     IF (IWK.EQ.I) THEN                                              
         GOTO 480                                                       
        ENDIF                                                           
                                                                        
*-----CALCULATE PARAMETER FOR QR TRANSFORMATION----------------------   
                                                                        
        DPWK1 = (W(I - 1) - W(I)) * (W(I - 1) + W(I))                   
        DPWK2 = (DPFWK(I - 1) - DPFWK(I)) * (DPFWK(I - 1) + DPFWK(I))   
        DPWK3 = (W(IWK) - W(I)) * (W(IWK) + W(I))                       
        DPWK1 = (DPWK1 + DPWK2) / (2.D0 * DPFWK(I) * W(I - 1))          
        DPWK2 = DSQRT(DPWK1 * DPWK1 + 1.D0)                             
        DPWK2 = W(I - 1) / (DPWK1 + DSIGN(DPWK2,DPWK1)) - DPFWK(I)      
        DPWK3 = (DPWK3 + DPFWK(I) * DPWK2) / W(IWK)                     
        DPWK6 = W(IWK)                                                  
                                                                        
*-----PERFORM QR TRANSFORMATION--------------------------------------   
                                                                        
        DPWK1 = 1.D0                                                    
        DPWK2 = 1.D0                                                    
        DO 460 K = IWK,(I - 1)                                          
         DPWK5 = DPWK2 * DPFWK(K + 1)                                   
         DPWK4 = DPWK1 * DPFWK(K + 1)                                   
                                                                        
*-----CONSTRUCT GIVENS ROTATION (LEFT)--------------------------------- 
                                                                        
         DPFWK(K) = DSQRT(DPWK3 * DPWK3 + DPWK5 * DPWK5)                
         IF (DPFWK(K).GT.0.D0) THEN                                     
          DPWK1 = DPWK3 / DPFWK(K)                                      
          DPWK2 = DPWK5 / DPFWK(K)                                      
         ELSE                                                           
          DPWK1 = 1.D0                                                  
          DPWK2 = 0.D0                                                  
         ENDIF                                                          
                                                                        
*-----APPLY GIVENS ROTATION (LEFT) TO X-------------------------------- 
                                                                        
         DPWK3 = X(K)                                                   
         DPWK5 = X(K + 1)                                               
         X(K) = DPWK1 * DPWK3 + DPWK2 * DPWK5                           
         X(K + 1) = -DPWK2 * DPWK3 + DPWK1 * DPWK5                      
                                                                        
*-----APPLY GIVENS ROTATION (LEFT) TO BIDIAGONAL MATRIX---------------- 
                                                                        
         DPWK3 = DPWK6 * DPWK1 + DPWK4 * DPWK2                          
         DPWK4 = -DPWK6 * DPWK2 + DPWK4 * DPWK1                         
         DPWK5 = W(K + 1) * DPWK2                                       
         DPWK6 = W(K + 1) * DPWK1                                       
                                                                        
*-----CONSTRUCT GIVENS ROTATION (RIGHT)-------------------------------- 
                                                                        
         W(K) = DSQRT(DPWK3 * DPWK3 + DPWK5 * DPWK5)                    
         IF (W(K).GT.0.D0) THEN                                         
          DPWK1 = DPWK3 / W(K)                                          
          DPWK2 = DPWK5 / W(K)                                          
         ELSE                                                           
          DPWK1 = 1.D0                                                  
          DPWK2 = 0.D0                                                  
         ENDIF                                                          
                                                                        
*-----APPLY GIVENS ROTATION (RIGHT) TO VT------------------------------ 
                                                                        
         DO 450 L = 1,MB                                                
          DPWK3 = B(K,L)                                                
          DPWK5 = B(K + 1,L)                                            
          B(K,L) = DPWK1 * DPWK3 + DPWK2 * DPWK5                        
          B(K + 1,L) = -DPWK2 * DPWK3 + DPWK1 * DPWK5                   
450      CONTINUE                                                       
                                                                        
*-----APPLY GIVENS ROTATION (RIGHT) TO BIDIAGONAL MATRIX--------------- 
                                                                        
         DPWK3 = DPWK1 * DPWK4 + DPWK2 * DPWK6                          
         DPWK6 = -DPWK2 * DPWK4 + DPWK1 * DPWK6                         
460     CONTINUE                                                        
                                                                        
        DPFWK(IWK) = 0.D0                                               
        DPFWK(I) = DPWK3                                                
        W(I) = DPWK6                                                    
                                                                        
470    CONTINUE                                                         
       ERRTXT = ' CASVDE > ERROR DETECTED (2)'                          
       CALL WEERRM(ERRTXT,.TRUE.)                                       
                                                                        
*-----MAKE SINGULAR VALUE POSITIV-------------------------------------- 
                                                                        
480    IF (W(I).LT.0.D0) THEN                                           
        W(I) = -W(I)                                                    
        DO 490 J = 1,MB                                                 
         B(I,J) = -B(I,J)                                               
490     CONTINUE                                                        
       ENDIF                                                            
                                                                        
500   CONTINUE                                                          
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
*---------------------------------------------------------------------- 
* SUBROUTINE: MI(NIMIZE) Q(UADRATIC) F(UNCTION) W(ITH) C(ONSTRAINTS)    
*                                                                       
* PURPOSE   : MINIMIZE A QUADRATIC FUNCTION WITH POSITIVITY CONSTRAINTS 
*             USING A SIMPLIFIED QUADRATIC PROGRAMMING ALGORITHM        
*                                                                       
* VARIABLES : (ON INPUT)                                                
*             NB     : FIRST DIMENSION OF THE MATRIX B                  
*             MB     : SECOND DIMENSION OF THE MATRIX B                 
*             MBCO   : NUMBER OF CONSTRAINT VARIABLES                   
*             B      : MATRIX                                           
*             NBMAX  : LEADING DIMENSION OF THE FIELD B                 
*                                                                       
*             A      : VECTOR CONTAINING THE DATA                       
*                                                                       
*             ACTCON : GUESS OF THE SET OF ACTIVE CONSTRAINTS           
*                                                                       
*             IFWK   : INTEGER WORKSPACE                                
*             DPFWK  : DOUBLE PRECISION WORKSPACE                       
*                                                                       
*             (ON OUTPUT)                                               
*             MBFR   : NUMBER OF THE FREE VARIABLES AT THE MINIMUM      
*                                                                       
*             B      : MATRIX B MULTIPLIED BY AN ORTHOGONAL MATRIX      
*                      THE MBFR X MBFR PART IS AN UPPER TRIANGULAR      
*                      MATRIX AND BELONGS TO THE FREE VARIABLES         
*                                                                       
*             A      : VECTOR A MULTIPLIED BY AN ORTHOGONAL MATRIX      
*                                                                       
*             ACTCON : SET OF ACTIVE CONSTRAINTS AT THE MINIMUM         
*                                                                       
*             ORDER  : ORDER OF THE COLUMNS OF THE MATRIX B             
*                                                                       
* COMMENT   : THE SIZE OF THE INTEGER WORKSPACE MUST BE AT LEAST NB     
*                                                                       
*             THE SIZE OF THE DOUBLE PRECISION WORKSPACE MUST BE AT     
*             LEAST MB                                                  
*---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE MIQFWC(NB,MB,MBCO,MBFR,B,NBMAX,A,ACTCON,ORDER,IFWK,    
     *                  DPFWK)                                          
                                                                        
      INTEGER NB,MB,MBCO,MBFR,NBMAX                                     
      DOUBLE PRECISION B(NBMAX,MB)                                      
                                                                        
      DOUBLE PRECISION A(NB)                                            
                                                                        
      INTEGER ACTCON(MB)                                                
      INTEGER ORDER(MB)                                                 
                                                                        
      INTEGER IFWK(NB)                                                  
      DOUBLE PRECISION DPFWK(MB)                                        
                                                                        
      INTEGER MBUNCO                                                    
                                                                        
      DOUBLE PRECISION MIN                                              
                                                                        
      INTEGER I,J                                                       
                                                                        
      INTEGER IWK1,IWK2                                                 
                                                                        
      DOUBLE PRECISION DPWK1,DPWK2,DPWK3,DPWK4                          
                                                                        
*-----INITIALIZE VARIABLES--------------------------------------------- 
                                                                        
      MBUNCO = MB - MBCO                                                
      MBFR = MB                                                         
      DO 10 I = (MBUNCO + 1),MB                                         
       IF (ACTCON(I).EQ.1) THEN                                         
        MBFR = MBFR - 1                                                 
       ENDIF                                                            
10    CONTINUE                                                          
                                                                        
      DO 20 I = 1,MB                                                    
       ORDER(I) = I                                                     
20    CONTINUE                                                          
                                                                        
      MIN = 0.D0                                                        
      DO 30 I = 1,NB                                                    
       MIN = MIN + A(I) * A(I)                                          
30    CONTINUE                                                          
                                                                        
*-----SORT B FOR VARIABLES WITH ACTIVE CONSTRAINTS--------------------- 
                                                                        
      IWK1 = MB                                                         
      DO 60 I = (MBUNCO + 1),MBFR                                       
       IF (ACTCON(I).EQ.1) THEN                                         
40      IF (ACTCON(IWK1).EQ.1) THEN                                     
         IWK1 = IWK1 - 1                                                
         GOTO 40                                                        
        ENDIF                                                           
        DO 50 J = 1,NB                                                  
         DPWK1 = B(J,I)                                                 
         B(J,I) = B(J,IWK1)                                             
         B(J,IWK1) = DPWK1                                              
50      CONTINUE                                                        
        ORDER(I) = IWK1                                                 
        ORDER(IWK1) = I                                                 
        IWK1 = IWK1 - 1                                                 
       ENDIF                                                            
60    CONTINUE                                                          
                                                                        
*-----INITIALIZE QR DECOMPOSITION-------------------------------------- 
                                                                        
      CALL CAQRDE(NB,MB,B,NBMAX,DPFWK,IFWK)                             
      CALL AHOTRV(NB,MB,B,NBMAX,DPFWK,A)                                
      DO 70 I = 1,(MB - 1)                                              
       DO 70 J = (I + 1),MB                                             
        B(J,I) = 0.D0                                                   
70    CONTINUE                                                          
                                                                        
*-----MINIMIZE QUADRATIC FUNCTION WITH CONSTRAINTS--------------------- 
                                                                        
*-----CALCULATE SOLUTION FOR THE GIVEN SET OF CONSTRAINTS-------------- 
                                                                        
80    DO 90 I = (MBUNCO + 1),MBFR                                       
       DPFWK(I) = A(I)                                                  
90    CONTINUE                                                          
      IWK1 = 0                                                          
      DPWK1 = 0.D0                                                      
      DO 100 I = MBFR,(MBUNCO + 1),-1                                   
       DPWK2 = DPFWK(I) / B(I,I)                                        
       IF (DPWK2.LT.DPWK1) THEN                                         
        DPWK1 = DPWK2                                                   
        IWK1 = I                                                        
       ENDIF                                                            
       DO 100 J = (MBUNCO + 1),(I - 1)                                  
        DPFWK(J) = DPFWK(J) - B(J,I) * DPWK2                            
100   CONTINUE                                                          
                                                                        
      IF (IWK1.NE.0) THEN                                               
                                                                        
*-----ADD CONSTRAINT TO ACTIVE SET------------------------------------- 
                                                                        
*-----REORDER MATRIX--------------------------------------------------- 
                                                                        
       MBFR = MBFR - 1                                                  
       IWK2 = ORDER(IWK1)                                               
       DO 110 I = 1,IWK1                                                
        DPFWK(I) = B(I,IWK1)                                            
110    CONTINUE                                                         
       DO 120 I = IWK1,MBFR                                             
        ORDER(I) = ORDER(I + 1)                                         
        DO 120 J = 1,(I + 1)                                            
         B(J,I) = B(J,I + 1)                                            
120    CONTINUE                                                         
       ORDER(MBFR + 1) = IWK2                                           
       DO 130 I = 1,IWK1                                                
        B(I,MBFR + 1) = DPFWK(I)                                        
130    CONTINUE                                                         
       DO 140 I = (IWK1 + 1),MB                                         
        B(I,MBFR + 1) = 0.D0                                            
140    CONTINUE                                                         
                                                                        
*-----APPLY GIVENS ROTATIONS TO DELETE ELEMENTS BELOW THE DIAGONAL----- 
                                                                        
       DO 160 I = IWK1,MBFR                                             
        DPWK1 = B(I,I)                                                  
        DPWK2 = B(I + 1,I)                                              
        IF (DPWK2.NE.0.D0) THEN                                         
                                                                        
*-----CONSTRUCT GIVENS ROTATION---------------------------------------- 
                                                                        
         B(I,I) = DSQRT(DPWK1 * DPWK1 + DPWK2 * DPWK2)                  
         B(I + 1,I) = 0.D0                                              
         DPWK1 = DPWK1 / B(I,I)                                         
         DPWK2 = DPWK2 / B(I,I)                                         
                                                                        
*-----APPLY GIVENS ROTATION TO MATRIX---------------------------------- 
                                                                        
         DO 150 J = (I + 1),MB                                          
          DPWK3 = B(I,J)                                                
          DPWK4 = B(I + 1,J)                                            
          B(I,J) = DPWK1 * DPWK3 + DPWK2 * DPWK4                        
          B(I + 1,J) = -DPWK2 * DPWK3 + DPWK1 * DPWK4                   
150      CONTINUE                                                       
                                                                        
*-----APPLY GIVENS ROTATION TO VECTOR---------------------------------- 
                                                                        
         DPWK3 = A(I)                                                   
         DPWK4 = A(I + 1)                                               
         A(I) = DPWK1 * DPWK3 + DPWK2 * DPWK4                           
         A(I + 1) = -DPWK2 * DPWK3 + DPWK1 * DPWK4                      
        ENDIF                                                           
                                                                        
160    CONTINUE                                                         
                                                                        
       GOTO 80                                                          
                                                                        
      ENDIF                                                             
                                                                        
C-----TEST IF LAST ITERATION IMPROVED THE RESULT----------------------- 
                                                                        
      DPWK1 = 0.D0                                                      
      DO 170 I = (MBFR + 1),MB                                          
       DPWK1 = DPWK1 + A(I) * A(I)                                      
170   CONTINUE                                                          
                                                                        
      IF (DPWK1.LT.MIN) THEN                                            
                                                                        
       MIN = DPWK1                                                      
                                                                        
C-----CALCULATE LAGRANGE MULTIPLIERS----------------------------------- 
                                                                        
       IWK1 = 0                                                         
       DPWK1 = 0.D0                                                     
       DO 190 I = (MBFR + 1),MB                                         
        DPWK2 = 0.D0                                                    
        DO 180 J = (MBFR + 1),MB                                        
         DPWK2 = DPWK2 - B(J,I) * A(J)                                  
180     CONTINUE                                                        
        IF (DPWK2.LT.DPWK1) THEN                                        
         DPWK1 = DPWK2                                                  
         IWK1 = I                                                       
        ENDIF                                                           
190    CONTINUE                                                         
                                                                        
       IF (IWK1.NE.0) THEN                                              
                                                                        
c-----REMOVE CONSTRAINT FROM ACTIVE SET-------------------------------- 
                                                                        
c-----REORDER MATRIX--------------------------------------------------- 
                                                                        
        MBFR = MBFR + 1                                                 
        IWK2 = ORDER(IWK1)                                              
        ORDER(IWK1) = ORDER(MBFR)                                       
        ORDER(MBFR) = IWK2                                              
        DO 200 I = 1,MB                                                 
         DPWK1 = B(I,IWK1)                                              
         B(I,IWK1) = B(I,MBFR)                                          
         B(I,MBFR) = DPWK1                                              
200     CONTINUE                                                        
                                                                        
c-----APPLY HOUSEHOLDER TRANSFORMATION TO DELETE ELEMENTS BELOW THE---- 
c-----DIAGONAL--------------------------------------------------------- 
                                                                        
c-----CONSTRUCT HOUSEHOLDER TRANSFORMATION----------------------------- 
                                                                        
        DO 210 I = (MBFR + 1),MB                                        
         DPWK1 = DPWK1 + B(I,MBFR) * B(I,MBFR)                          
210     CONTINUE                                                        
        IF (DPWK1.GT.0.D0) THEN                                         
         DPWK1 = -DSIGN(DSQRT(DPWK1 + B(MBFR,MBFR) * B(MBFR,MBFR)),     
     *                  B(MBFR,MBFR))                                   
         B(MBFR,MBFR) = B(MBFR,MBFR) - DPWK1                            
                                                                        
c-----APPLY HOUSEHOLDER TRANSFORMATION TO MATRIX----------------------- 
                                                                        
         DO 230 I = (MBFR + 1),MB                                       
          DPWK2 = 0.D0                                                  
          DO 220 J = MBFR,MB                                            
           DPWK2 = DPWK2 + B(J,I) * B(J,MBFR)                           
220       CONTINUE                                                      
          DPWK2 = DPWK2 / DPWK1 / B(MBFR,MBFR)                          
          DO 230 J = MBFR,MB                                            
           B(J,I) = B(J,I) + DPWK2 * B(J,MBFR)                          
230      CONTINUE                                                       
                                                                        
c-----APPLY HOUSEHOLDER TRANSFORMATION TO VECTOR----------------------- 
                                                                        
         DPWK2 = 0.D0                                                   
         DO 240 I = MBFR,MB                                             
          DPWK2 = DPWK2 + A(I) * B(I,MBFR)                              
240      CONTINUE                                                       
         DPWK2 = DPWK2 / DPWK1 / B(MBFR,MBFR)                           
         DO 250 I = MBFR,MB                                             
          A(I) = A(I) + DPWK2 * B(I,MBFR)                               
250      CONTINUE                                                       
                                                                        
         B(MBFR,MBFR) = DPWK1                                           
         DO 260 I = (MBFR + 1),MB                                       
          B(I,MBFR) = 0.D0                                              
260      CONTINUE                                                       
                                                                        
        ENDIF                                                           
                                                                        
        GOTO 80                                                         
                                                                        
       ENDIF                                                            
                                                                        
      ENDIF                                                             
                                                                        
c-----UPDATE THE FIELD OF THE ACTIVE CONSTRAINTS----------------------- 
                                                                        
      DO 270 I = 1,MBFR                                                 
       ACTCON(ORDER(I)) = 0                                             
270   CONTINUE                                                          
                                                                        
      DO 280 I = (MBFR + 1),MB                                          
       ACTCON(ORDER(I)) = 1                                             
280   CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
c---------------------------------------------------------------------- 
c FUNCTION  : CA(LCULATE) NORM                                          
c                                                                       
c PURPOSE   : CALCULATE THE NORM OF A MATRIX                            
c                                                                       
c VARIABLES : (ON INPUT)                                                
c             NB     : FIRST DIMENSION OF THE MATRIX B                  
c             MB     : SECOND DIMENSION OF THE MATRIX B                 
c             B      : MATRIX                                           
c             NBMAX  : LEADING DIMENSION OF THE FIELD B                 
c---------------------------------------------------------------------- 
                                                                        
      DOUBLE PRECISION FUNCTION CANORM(NB,MB,B,NBMAX)                   
                                                                        
      INTEGER NB,MB,NBMAX                                               
      DOUBLE PRECISION B(NBMAX,MB)                                      
                                                                        
      INTEGER I,J                                                       
      DOUBLE PRECISION DPWK                                             
                                                                        
      CANORM = 0.D0                                                     
      DO 20 I = 1,NB                                                    
       DPWK = 0.D0                                                      
       DO 10 J = 1,MB                                                   
        DPWK = DPWK + DABS(B(I,J))                                      
10     CONTINUE                                                         
       CANORM = DMAX1(CANORM,DPWK)                                      
20    CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
c---------------------------------------------------------------------- 
c SUBROUTINE: W(RIT)E ERR(OR) M(ESSAGE)                                 
c                                                                       
c PURPOSE   : WRITE ERROR MSG. TO OUTPUT AND STOP PROGRAM IF NECESSARY  
c                                                                       
c VARIABLES : ERRTXT : ERROR MESSAGE                                    
c             STOPFL : FLAG FOR PROGRAMSTOP                             
c---------------------------------------------------------------------- 
                                                                        
      SUBROUTINE WEERRM(ERRTXT,STOPFL)                                  
                                                                        
      CHARACTER*60 ERRTXT                                               
                                                                        
      LOGICAL STOPFL                                                    
                                                                        
      WRITE(*,*)'--------------------------------------------------'    
      WRITE(*,*)ERRTXT                                                  
                                                                        
      IF (STOPFL) THEN                                                  
       WRITE(*,*)'-----------------FTIKREG STOPPED------------------'   
       WRITE(*,*)'--------------------------------------------------'   
       STOP                                                             
      ENDIF                                                             
                                                                        
      WRITE(*,*)'----------------FTIKREG CONTINUES-----------------'    
      WRITE(*,*)'--------------------------------------------------'    
                                                                        
      RETURN                                                            
      END                                                               
c-------------------------------------------------------------------------------

      SUBROUTINE MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU,
     &                  NMOM, IPOLZN, MOMDIM, PRNT, QEXT, QSCA, GQSC,
     &                  PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK,
     &                  SPIKE )

c       NoPMOM version:  saves storage by deleting everything
c       having to do with PMOM calculation:

c       * deletes modules LPCoef, LPCO1T, LPCO2T

c       * NMOM,IPOLZN,PMOM,MOMDIM,CalcMo,NPQuan deleted from argument
c         lists and everywhere else they occur internally, although the
c         first four are retained in the MIEV0 argument list for
c         compatibility with complete version of MIEV0

c       * MIEV0:   shrinks size of LitA,LitB arrays to 2;  deletes call
c                  to LPCoef;  saving of An,Bn in LitA,LitB arrays
c                  deleted

c       * MiPrnt:  PMOM print deleted

c       * CkInMi:  NMOM test changed to allow only zero value (nonzero
c                  value would have no impact but indicates user
c                  conceptual error or negligence)

c       * TestMi:  saving and restoring of NMOM,IPOLZN deleted;
c                  30-loop deleted

c    Computes Mie scattering and extinction efficiencies; asymmetry
c    factor;  forward- and backscatter amplitude;  scattering
c    amplitudes vs. scattering angle for incident polarization parallel
c    and perpendicular to the plane of scattering;
c    some quantities needed in polarized radiative transfer;  and
c    information about whether or not a resonance has been hit.

c    Input and output arguments are described in file MIEV.doc;
c    internal variables are described in MIEV0.f file.
c    Many statements are accompanied by comments referring to 
c    references in MIEV.doc, notably the NCAR Mie report which is now
c    available electronically and which is referred to using the
c    shorthand (Rn), meaning Eq. (n) of the report.

c    CALLING TREE:

c        MIEV0
c            TESTMI
c                TSTBAD
c                MIPRNT
c                ERRMSG
c            CKINMI
c                WRTBAD
c                WRTDIM
c                ERRMSG
c            SMALL1
c            SMALL2
c            ERRMSG
c            BIGA
c                CONFRA
c                    ERRMSG
c            MIPRNT
c ----------------------------------------------------------------------


      IMPLICIT  NONE

c ----------------------------------------------------------------------
c --------  I / O SPECIFICATIONS FOR SUBROUTINE MIEV0  -----------------
c ----------------------------------------------------------------------
      LOGICAL  ANYANG, PERFCT, PRNT(*)
      INTEGER  IPOLZN, MOMDIM, NUMANG, NMOM
      REAL*8     GQSC, MIMCUT, PMOM( 0:MOMDIM, * ), QEXT, QSCA, SPIKE,
     &         XMU(*), XX
      COMPLEX  CREFIN, SFORW, SBACK, S1(*), S2(*), TFORW(*), TBACK(*)
c ----------------------------------------------------------------------

c                                  ** NOTE --  MaxTrm = 10100  is neces-
c                                  ** sary to do some of the test probs,
c                                  ** but 1100 is sufficient for most
c                                  ** conceivable applications
c     .. Parameters ..

      INTEGER   MAXANG, MXANG2
      PARAMETER ( MAXANG = 501, MXANG2 = MAXANG / 2 + 1 )
      INTEGER   MAXTRM
      PARAMETER ( MAXTRM = 10100 )
      REAL*8      ONETHR
      PARAMETER ( ONETHR = 1. / 3. )
c     ..
c     .. Local Scalars ..

      LOGICAL   NOABS, PASS1, YESANG
      INTEGER   I, J, N, NANGD2, NTRM
      REAL*8      CHIN, CHINM1, COEFF, DENAN, DENBN, FN, MIM, MM, MRE,
     &          NP1DN, PSIN, PSINM1, RATIO, RIORIV, RN, RTMP, TAUN,
     &          TCOEF, TWONP1, XINV
      COMPLEX   AN, ANM1, ANP, ANPM, BN, BNM1, BNP, BNPM, CDENAN,
     &          CDENBN, CIOR, CIORIV, CSUM1, CSUM2, CTMP, ZET, 
     &          ZETN, ZETNM1
c     ..
c     .. Local Arrays ..

      REAL*8      PIN( MAXANG ), PINM1( MAXANG ), RBIGA( MAXTRM )
      COMPLEX   CBIGA( MAXTRM ), LITA( 2 ), LITB( 2 ), SM( MAXANG ),
     &          SMS( MXANG2 ), SP( MAXANG ), SPS( MXANG2 )
c     ..
c     .. External Subroutines ..

      EXTERNAL  BIGA, CKINMI, ERRMSG, MIPRNT, SMALL1, SMALL2, TESTMI
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, CMPLX, CONJG, COS, MAX, MIN, REAL, SIN
c     ..
      SAVE      PASS1

c     .. Statement Functions ..

      REAL*8      SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = REAL( CTMP )**2 + AIMAG( CTMP )**2
c     ..
      DATA      PASS1 /.TRUE./


c                    ** Save some input variables and replace them
c                    ** with values needed to do the self-test

      IF( PASS1 ) CALL TESTMI( .FALSE., XX, CREFIN, MIMCUT, PERFCT,
     &                         ANYANG, NUMANG, XMU, QEXT, QSCA, GQSC,
     &                         SFORW, SBACK, S1, S2, TFORW, TBACK )

   10 CONTINUE
c                                        ** Check input and calculate
c                                        ** certain variables from input

      CALL CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, NMOM, ANYANG,
     &             XMU )


      IF( PERFCT .AND. XX.LE.0.1 ) THEN
c                                            ** Use totally-reflecting
c                                            ** small-particle limit

         CALL SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK,
     &                S1, S2, TFORW, TBACK, LITA, LITB )

         NTRM = 2
         GO TO  100

      END IF


      NOABS = .TRUE.

      IF( .NOT.PERFCT ) THEN

         CIOR = CREFIN
         IF( AIMAG( CIOR ).GT.0.0 ) CIOR   = CONJG( CIOR )

         MRE    =   REAL( CIOR )
         MIM    = - AIMAG( CIOR )
         NOABS  = MIM.LE.MIMCUT
         CIORIV = 1.0 / CIOR
         RIORIV = 1.0 / MRE


         IF( XX*MAX( 1.0, ABS(CIOR) ) .LE. 0.1 ) THEN

c                                    ** Use general-refractive-index
c                                    ** small-particle limit

            CALL SMALL2( XX, CIOR, MIM.GT.MIMCUT, NUMANG, XMU, QEXT,
     &                   QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW,
     &                   TBACK, LITA, LITB )

            NTRM = 2
            GO TO  100

         END IF

      END IF


      NANGD2 = ( NUMANG + 1 ) / 2
      YESANG = NUMANG.GT.0

c                              ** Number of terms in Mie series; Eq R50
      IF( XX.LE.8.0 ) THEN

         NTRM   = XX + 4.*XX**ONETHR + 1.

      ELSE IF( XX.LT.4200. ) THEN

         NTRM   = XX + 4.05*XX**ONETHR + 2.

      ELSE

         NTRM   = XX + 4.*XX**ONETHR + 2.
         
      END IF

      IF( NTRM + 1.GT.MAXTRM ) 
     &    CALL ERRMSG('MIEV0--PARAMETER MaxTrm TOO SMALL',.TRUE.)

c                            ** Calculate logarithmic derivatives of
c                            ** J-Bessel-fcn., A-sub-(1 to NTrm)

      IF( .NOT.PERFCT ) CALL BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA,
     &                             CBIGA )

c                            ** Initialize Ricatti-Bessel functions
c                            ** (psi,chi,zeta)-sub-(0,1) for upward
c                            ** recurrence ( Eq. R19 )
      XINV   = 1.0 / XX
      PSINM1 = SIN( XX )
      CHINM1 = COS( XX )
      PSIN   = PSINM1*XINV - CHINM1
      CHIN   = CHINM1*XINV + PSINM1
      ZETNM1 = CMPLX( PSINM1, CHINM1 )
      ZETN   = CMPLX( PSIN, CHIN )
c                                     ** Initialize previous coeffi-
c                                     ** cients for GQSC series
      ANM1 = ( 0.0, 0.0 )
      BNM1 = ( 0.0, 0.0 )
c                             ** Initialize angular function  pi
c                             ** and sums for S+, S- ( Ref. 2, p. 1507 )
      IF( ANYANG ) THEN

         DO 20 J = 1, NUMANG
c                             ** Eq. R39
            PINM1( J ) = 0.0
            PIN( J )   = 1.0
            SP( J )    = ( 0.0, 0.0 )
            SM( J )    = ( 0.0, 0.0 )
   20    CONTINUE

      ELSE

         DO 30 J = 1, NANGD2
c                             ** Eq. R39
            PINM1( J ) = 0.0
            PIN( J )   = 1.0
            SP( J )    = ( 0.0, 0.0 )
            SM( J )    = ( 0.0, 0.0 )
            SPS( J )   = ( 0.0, 0.0 )
            SMS( J )   = ( 0.0, 0.0 )
   30    CONTINUE

      END IF
      
c                       ** Initialize Mie sums for efficiencies, etc.
      QSCA  = 0.0
      GQSC  = 0.0
      SFORW = ( 0., 0.)
      SBACK = ( 0., 0.)
      CSUM1 = ( 0., 0. )
      CSUM2 = ( 0., 0. )


c ---------  LOOP TO SUM MIE SERIES  -----------------------------------

      MM     = +1.0
      SPIKE  = 1.0

      DO 60  N = 1, NTRM
c                           ** Compute various numerical coefficients
         FN     = N
         RN     = 1.0 / FN
         NP1DN  = 1.0 + RN
         TWONP1 = 2*N + 1
         COEFF  = TWONP1 / ( FN*( N + 1 ) )
         TCOEF  = TWONP1*( FN*( N + 1 ) )

c                              ** Calculate Mie series coefficients
         IF( PERFCT ) THEN
c                                 ** Totally-reflecting case; Eq R/A.1,2

            AN     = ( ( FN*XINV )*PSIN - PSINM1 ) /
     &               ( ( FN*XINV )*ZETN - ZETNM1 )
            BN     = PSIN / ZETN

         ELSE IF( NOABS ) THEN
c                                      ** No-absorption case; Eq (R16)

            CDENAN = ( RIORIV*RBIGA(N) + (FN*XINV) ) * ZETN - ZETNM1
            AN =   ( ( RIORIV*RBIGA(N) + (FN*XINV) ) * PSIN - PSINM1 )
     &             / CDENAN
            CDENBN = (  MRE * RBIGA(N) + (FN*XINV) ) * ZETN - ZETNM1
            BN =   ( (  MRE * RBIGA(N) + (FN*XINV) ) * PSIN - PSINM1 )
     &             / CDENBN

         ELSE
c                                       ** Absorptive case; Eq (R16)

            CDENAN = ( CIORIV * CBIGA(N) + (FN*XINV) ) * ZETN - ZETNM1
            CDENBN = (   CIOR * CBIGA(N) + (FN*XINV) ) * ZETN - ZETNM1
            AN =   ( ( CIORIV * CBIGA(N) + (FN*XINV) ) * PSIN - PSINM1 )
     &             / CDENAN
            BN =   ( (   CIOR * CBIGA(N) + (FN*XINV) ) * PSIN - PSINM1 )
     &             / CDENBN
c                                         ** Eq (R7)

            QSCA = QSCA + TWONP1 * ( SQ( AN ) + SQ( BN ) )

         END IF


         IF( .NOT.PERFCT .AND. N.GT.XX ) THEN
c                                               ** Flag resonance spikes
            DENAN  = ABS( CDENAN )
            DENBN  = ABS( CDENBN )
c                                                   ** Eq. R/B.9
            RATIO  = DENAN / DENBN
c                                                   ** Eq. R/B.10
            IF( RATIO.LE.0.2 .OR. RATIO.GE.5.0 ) 
     &          SPIKE  = MIN( SPIKE, DENAN, DENBN )

         END IF
c                                  ** Increment Mie sums for non-angle-
c                                  ** dependent quantities

c                                                   ** Eq. R/B.2
         SFORW = SFORW + TWONP1 * ( AN + BN )
c                                                   ** Eq. R/B.5,6
         CSUM1 = CSUM1 + TCOEF *( AN - BN )
c                                                   ** Eq. R/B.1
         SBACK = SBACK      + ( MM * TWONP1 ) * ( AN - BN )
c                                                   ** Eq. R/B.7,8
         CSUM2 = CSUM2 + ( MM*TCOEF ) *( AN + BN )

c                                         ** Eq (R8)

         GQSC  = GQSC + ( FN - RN ) * REAL( ANM1 * CONJG( AN )
     &                                    + BNM1 * CONJG( BN ) )
     &          + COEFF * REAL( AN * CONJG( BN ) )


         IF( YESANG ) THEN
c                                      ** Put Mie coefficients in form
c                                      ** needed for computing S+, S-
c                                      ** ( Eq R10 )
            ANP = COEFF * ( AN + BN )
            BNP = COEFF * ( AN - BN )
c                                      ** Increment Mie sums for S+, S-
c                                      ** while upward recursing
c                                      ** angular functions pi and tau
            IF( ANYANG ) THEN
c                                         ** Arbitrary angles

c                                              ** vectorizable loop
               DO 40 J = 1, NUMANG
c                                                 ** Eq. (R37b)

                  RTMP = ( XMU( J ) * PIN( J ) ) - PINM1( J )

c                                                 ** Eq. (R38b)
                  TAUN =  FN * RTMP - PINM1( J )

c                                                   ** Eq (R10)

                  SP( J )  = SP( J ) + ANP * ( PIN( J ) + TAUN )
                  SM( J )  = SM( J ) + BNP * ( PIN( J ) - TAUN )
                  PINM1( J ) = PIN( J )
c                                                 ** Eq. R37c

                  PIN( J ) = ( XMU( J ) * PIN( J ) ) + NP1DN * RTMP
   40          CONTINUE

            ELSE
c                                  ** Angles symmetric about 90 degrees
               ANPM   = MM * ANP
               BNPM   = MM * BNP
c                                          ** vectorizable loop
               DO 50 J = 1, NANGD2
c                                                 ** Eq. (R37b)

                  RTMP = ( XMU( J ) * PIN( J ) ) - PINM1( J )

c                                                 ** Eq. (R38b)
                  TAUN =  FN * RTMP - PINM1( J )

c                                                 ** Eq (R10,12)

                  SP ( J ) = SP ( J ) +  ANP * ( PIN( J ) + TAUN )
                  SMS( J ) = SMS( J ) + BNPM * ( PIN( J ) + TAUN )
                  SM ( J ) = SM ( J ) +  BNP * ( PIN( J ) - TAUN )
                  SPS( J ) = SPS( J ) + ANPM * ( PIN( J ) - TAUN )
                  PINM1( J ) = PIN( J )
c                                                 ** Eq. R37c

                  PIN( J ) = ( XMU( J ) * PIN( J ) ) + NP1DN * RTMP
   50          CONTINUE

            END IF

         END IF
c                          ** Update relevant quantities for next
c                          ** pass through loop
         MM     = -MM
         ANM1   = AN
         BNM1   = BN
c                           ** Upward recurrence for Ricatti-Bessel
c                           ** functions ( Eq. R17 )

         ZET    = ( TWONP1*XINV )*ZETN - ZETNM1
         ZETNM1 = ZETN
         ZETN   = ZET
         PSINM1 = PSIN
         PSIN   = REAL( ZETN )

   60 CONTINUE

c ---------- END LOOP TO SUM MIE SERIES --------------------------------


c                                         ** Eq (R6)
      QEXT   = 2./ XX**2 * REAL( SFORW )

      IF( PERFCT .OR. NOABS ) THEN

         QSCA = QEXT

      ELSE

         QSCA = 2./ XX**2 * QSCA
         
      END IF

      GQSC   = 4./ XX**2 * GQSC
      SFORW = 0.5 * SFORW
      SBACK = 0.5 * SBACK
      TFORW( 1 ) =  0.5*SFORW - 0.125*CSUM1
      TFORW( 2 ) =  0.5*SFORW + 0.125*CSUM1
      TBACK( 1 ) = -0.5*SBACK + 0.125*CSUM2
      TBACK( 2 ) =  0.5*SBACK + 0.125*CSUM2


      IF( YESANG ) THEN
c                                ** Recover scattering amplitudes
c                                ** from S+, S- ( Eq (R11) )

         IF( ANYANG ) THEN
c                                         ** vectorizable loop
            DO 70 J = 1, NUMANG
c                                                   ** Eq (R11)
               S1( J ) = 0.5*( SP( J ) + SM( J ) )
               S2( J ) = 0.5*( SP( J ) - SM( J ) )
   70       CONTINUE

         ELSE
c                                         ** vectorizable loop
            DO 80 J = 1, NANGD2
c                                                   ** Eq (R11)
               S1( J ) = 0.5*( SP( J ) + SM( J ) )
               S2( J ) = 0.5*( SP( J ) - SM( J ) )
   80       CONTINUE
c                                         ** vectorizable loop
            DO 90 J = 1, NANGD2
               S1( NUMANG + 1 - J ) = 0.5*( SPS( J ) + SMS( J ) )
               S2( NUMANG + 1 - J ) = 0.5*( SPS( J ) - SMS( J ) )
   90       CONTINUE

         END IF

      END IF


  100 CONTINUE
      IF( AIMAG( CREFIN ).GT.0.0 ) THEN
c                                         ** Take complex conjugates
c                                         ** of scattering amplitudes
         SFORW  = CONJG( SFORW )
         SBACK  = CONJG( SBACK )

         DO 110 I = 1, 2
            TFORW( I ) = CONJG( TFORW( I ) )
            TBACK( I ) = CONJG( TBACK( I ) )
  110    CONTINUE

         DO 120 J = 1, NUMANG
            S1( J ) = CONJG( S1( J ) )
            S2( J ) = CONJG( S2( J ) )
  120    CONTINUE

      END IF
 

      IF( PASS1 ) THEN
c                           ** Compare test case results with
c                           ** correct answers and abort if bad;
c                           ** otherwise restore user input and proceed

         CALL TESTMI( .TRUE., XX, CREFIN, MIMCUT, PERFCT, ANYANG,
     &                NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK, S1,
     &                S2, TFORW, TBACK )

         PASS1  = .FALSE.
         GO TO  10

      END IF


      IF ( PRNT(1) .OR. PRNT(2) )
     &   CALL  MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                 QSCA, GQSC, SFORW, SBACK, TFORW, TBACK, S1, S2 )

      RETURN

      END

      SUBROUTINE CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, NMOM,
     &                   ANYANG, XMU )

c        Check for bad input to MIEV0
c        (NoPMOM version)

c     Routines called :  ERRMSG, WRTBAD, WRTDIM


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   ANYANG, PERFCT
      INTEGER   MAXANG, NMOM, NUMANG
      REAL*8      XX
      COMPLEX   CREFIN
c     ..
c     .. Array Arguments ..

      REAL*8      XMU( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   INPERR
      INTEGER   I
c     ..
c     .. External Functions ..

      LOGICAL   WRTBAD, WRTDIM
      EXTERNAL  WRTBAD, WRTDIM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, REAL
c     ..


      INPERR = .FALSE.

      IF( NUMANG.GT.MAXANG ) INPERR = WRTDIM( 'MaxAng', NUMANG )
      IF( NUMANG.LT.0 ) INPERR = WRTBAD( 'NUMANG' )

      IF( XX.LT.0.) INPERR = WRTBAD( 'XX' )

      IF( .NOT.PERFCT .AND. REAL( CREFIN ).LE.0. )
     &    INPERR = WRTBAD( 'CREFIN' )


      IF( NMOM.NE.0 ) THEN

         INPERR = WRTBAD( 'NMOM' )
         WRITE( *, '(A)' )
     &      ' For nonzero NMOM, use complete version of MIEV0'
      END IF


      IF( ANYANG ) THEN
c                                ** Allow for slight imperfections in
c                                ** computation of cosine
         DO 10 I = 1, NUMANG

             IF ( XMU(I).LT.-1.00001 .OR. XMU(I).GT.1.00001 )
     &            INPERR = WRTBAD( 'XMU' )

   10    CONTINUE

      ELSE

         DO 20 I = 1, ( NUMANG + 1 ) / 2

             IF ( XMU(I).LT.-0.00001 .OR. XMU(I).GT.1.00001 )
     &            INPERR = WRTBAD( 'XMU' )

   20    CONTINUE

      END IF


      IF( INPERR ) CALL ERRMSG( 'MIEV0--INPUT ERROR(S).  Aborting...',
     &                          .TRUE.)

      IF( XX.GT.20000.0 .OR. REAL( CREFIN ).GT.10.0 .OR.
     &    ABS( AIMAG( CREFIN ) ).GT.10.0 )
     &    CALL ERRMSG( 'MIEV0--XX or CREFIN outside tested range',
     &    .FALSE. )

      RETURN
      END

      SUBROUTINE BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA, CBIGA )

c        Calculate logarithmic derivatives of J-Bessel-function

c     Input :  CIOR, XX, NTRM, NOABS, YESANG  (defined in MIEV0)

c    Output :  RBIGA or CBIGA  (defined in MIEV0)

c    Routines called :  CONFRA

c    INTERNAL VARIABLES :

c       CONFRA     Value of Lentz continued fraction for cBigA(NTrm),
c                     used to initialize downward recurrence

c       DOWN       = True, use down-recurrence.  False, do not.

c       F1,F2,F3   Arithmetic statement functions used in determining
c                     whether to use up-  or down-recurrence
c                     ( Eqs. R47a,b,c )

c       MRE        Real refractive index
c       MIM        Imaginary refractive index

c       REZINV     1 / ( MRE * XX ); temporary variable for recurrence
c       ZINV       1 / ( CIOR * XX ); temporary variable for recurrence


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   NOABS, YESANG
      INTEGER   NTRM
      REAL*8      XX
      COMPLEX   CIOR
c     ..
c     .. Array Arguments ..

      REAL*8      RBIGA( * )
      COMPLEX   CBIGA( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   DOWN
      INTEGER   N
      REAL*8      MIM, MRE, REZINV, RTMP
      COMPLEX   CTMP, ZINV
c     ..
c     .. External Functions ..

      COMPLEX   CONFRA
      EXTERNAL  CONFRA
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, COS, EXP, REAL, SIN
c     ..
c     .. Statement Functions ..

      REAL*8      F1, F2, F3
c     ..
c     .. Statement Function definitions ..

c                                                   ** Eq. R47c
      F1( MRE ) = -8.0 + MRE**2*( 26.22 +
     &            MRE*( -0.4474 + MRE**3*( 0.00204 - 0.000175*MRE ) ) )
     
c                                                   ** Eq. R47b
      F2( MRE ) = 3.9 + MRE*( -10.8 + 13.78*MRE )
c                                                   ** Eq. R47a      
      F3( MRE ) = -15.04 + MRE*( 8.42 + 16.35*MRE )
c     ..

c                                  ** Decide whether BigA can be
c                                  ** calculated by up-recurrence
      MRE = REAL( CIOR )
      MIM = ABS( AIMAG( CIOR ) )

      IF( MRE.LT.1.0 .OR. MRE.GT.10.0 .OR. MIM.GT.10.0 ) THEN

         DOWN = .TRUE.

      ELSE IF( YESANG ) THEN

         DOWN = .TRUE.
c                                                    ** Eq. R48
         IF( MIM*XX .LT. F2( MRE ) ) DOWN   = .FALSE.

      ELSE

         DOWN = .TRUE.
c                                                    ** Eq. R48
         IF( MIM*XX .LT. F1( MRE ) ) DOWN   = .FALSE.

      END IF


      ZINV   = 1.0 / ( CIOR*XX )
      REZINV = 1.0 / ( MRE*XX )

      IF( DOWN ) THEN
c                          ** Compute initial high-order BigA using
c                          ** Lentz method ( Ref. 1, pp. 17-20 )

         CTMP = CONFRA( NTRM, ZINV )

c                                   *** Downward recurrence for BigA
         IF( NOABS ) THEN
c                                        ** No-absorption case; Eq (R23)
            RBIGA( NTRM ) = REAL( CTMP )

            DO 10 N = NTRM, 2, -1
               RBIGA( N - 1 ) = ( N*REZINV ) -
     &                          1.0 / ( ( N*REZINV ) + RBIGA( N ) )
   10       CONTINUE

         ELSE
c                                         ** Absorptive case; Eq (R23)
            CBIGA( NTRM ) = CTMP

            DO 20 N = NTRM, 2, -1
               CBIGA( N-1 ) = (N*ZINV) - 1.0 / ( (N*ZINV) + CBIGA( N ) )
   20       CONTINUE

         END IF


      ELSE
c                              *** Upward recurrence for BigA
         IF( NOABS ) THEN
c                                  ** No-absorption case; Eq (R20,21)
            RTMP   = SIN( MRE*XX )
            RBIGA( 1 ) = -REZINV + RTMP /
     &                   ( RTMP*REZINV - COS( MRE*XX ) )

            DO 30 N = 2, NTRM
               RBIGA( N ) = -( N*REZINV ) +
     &                      1.0 / ( ( N*REZINV ) - RBIGA( N - 1 ) )
   30       CONTINUE

         ELSE
c                                     ** Absorptive case; Eq (R20,22)

            CTMP = EXP( - (0.,2.) * CIOR * XX )
            CBIGA( 1 ) = - ZINV + (1.-CTMP) /
     &                          ( ZINV * (1.-CTMP) - (0.,1.)*(1.+CTMP) )

            DO 40 N = 2, NTRM
               CBIGA( N ) = - (N*ZINV) + 1.0 / ((N*ZINV) - CBIGA( N-1 ))
   40       CONTINUE

         END IF

      END IF


      RETURN
      END

      COMPLEX FUNCTION CONFRA( N, ZINV )

c         Compute Bessel function ratio A-sub-N from its
c         continued fraction using Lentz method

c         ZINV = Reciprocal of argument of A


c    I N T E R N A L    V A R I A B L E S
c    ------------------------------------

c    CAK      Term in continued fraction expansion of A (Eq. R25)

c    CAPT     Factor used in Lentz iteration for A (Eq. R27)

c    CNUMER   Numerator   in capT  ( Eq. R28A )
c    CDENOM   Denominator in capT  ( Eq. R28B )

c    CDTD     Product of two successive denominators of capT factors
c                 ( Eq. R34C )
c    CNTN     Product of two successive numerators of capT factors
c                 ( Eq. R34B )

c    EPS1     Ill-conditioning criterion
c    EPS2     Convergence criterion

c    KK       Subscript k of cAk  ( Eq. R25B )

c    KOUNT    Iteration counter ( used to prevent infinite looping )

c    MAXIT    Max. allowed no. of iterations

c    MM       + 1  and - 1, alternately
c --------------------------------------------------------------------

      IMPLICIT  NONE

c     .. Scalar Arguments ..

      INTEGER   N
      COMPLEX   ZINV
c     ..
c     .. Local Scalars ..

      INTEGER   KK, KOUNT, MAXIT, MM
      REAL*8      EPS1, EPS2
      COMPLEX   CAK, CAPT, CDENOM, CDTD, CNTN, CNUMER
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, REAL
c     ..
      DATA      EPS1 / 1.E-2 / , EPS2 / 1.E-8 /
      DATA      MAXIT / 10000 /


c                                 ** Eq. R25a
      CONFRA = ( N + 1 ) * ZINV
      MM     = - 1
      KK     = 2*N + 3
c                                 ** Eq. R25b, k=2
      CAK    = ( MM*KK ) * ZINV
      CDENOM = CAK
      CNUMER = CDENOM + 1.0 / CONFRA
      KOUNT  = 1

   10 CONTINUE
      KOUNT = KOUNT + 1

      IF( KOUNT.GT.MAXIT )
     &    CALL ERRMSG('ConFra--Iteration failed to converge',.TRUE.)

      MM  = - MM
      KK  = KK + 2
c                                 ** Eq. R25b
      CAK = ( MM*KK ) * ZINV
c                                          ** Eq. R32
      IF( ABS( CNUMER / CAK ).LE.EPS1 .OR.
     &    ABS( CDENOM / CAK ).LE.EPS1 ) THEN

c                                  ** Ill-conditioned case -- stride
c                                  ** two terms instead of one

c                                       ** Eq. R34
         CNTN   = CAK * CNUMER + 1.0
         CDTD   = CAK * CDENOM + 1.0
c                                           ** Eq. R33
         CONFRA = ( CNTN / CDTD ) * CONFRA

         MM  = - MM
         KK  = KK + 2
c                                 ** Eq. R25b
         CAK = ( MM*KK ) * ZINV
c                                      ** Eq. R35
         CNUMER = CAK + CNUMER / CNTN
         CDENOM = CAK + CDENOM / CDTD
         KOUNT  = KOUNT + 1
         GO TO  10

      ELSE
c                           *** Well-conditioned case

c                                  ** Eq. R27
         CAPT   = CNUMER / CDENOM
c                                  ** Eq. R26
         CONFRA = CAPT * CONFRA
c                                  ** Check for convergence; Eq. R31

         IF (      ABS( REAL (CAPT) - 1.0 ).GE.EPS2
     &        .OR. ABS( AIMAG(CAPT) )      .GE.EPS2 )  THEN

c                                        ** Eq. R30
            CNUMER = CAK + 1.0 / CNUMER
            CDENOM = CAK + 1.0 / CDENOM
            GO TO  10

         END IF

      END IF


      RETURN

      END

      SUBROUTINE MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                   QSCA, GQSC, SFORW, SBACK, TFORW, TBACK, S1,
     &                   S2 )

c         Print scattering quantities of a single particle
c         (NoPMOM version)


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   PERFCT
      INTEGER   NUMANG
      REAL*8      GQSC, QEXT, QSCA, XX
      COMPLEX   CREFIN, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      LOGICAL   PRNT( * )
      REAL*8      XMU( * )
      COMPLEX   S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      INTEGER   I
      REAL*8      I1, I2
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC AIMAG, CONJG, REAL
c     ..

      IF( PERFCT ) WRITE( *, '(''1'',10X,A,1P,E11.4)' )
     &    'Perfectly Conducting Case, size parameter =', XX

c      IF( .NOT.PERFCT ) WRITE( *, '(''1'',10X,3(A,1P,E11.4))' )
c     &    'Refractive Index:  Real ', REAL( CREFIN ), '  Imag ',
c     &    AIMAG( CREFIN ), ',   Size Parameter =', XX


      IF( PRNT( 1 ) .AND. NUMANG.GT.0 ) THEN

         WRITE( *, '(/,A)' )
     &      '    cos(angle)  ------- S1 ---------  ------- S2 ---------'
     &      // '  --- S1*conjg(S2) ---   i1=S1**2   i2=S2**2  (i1+i2)/2'
     &      // '  DEG POLZN'

         DO 10 I = 1, NUMANG
            I1     = REAL( S1( I ) )**2 + AIMAG( S1( I ) )**2
            I2     = REAL( S2( I ) )**2 + AIMAG( S2( I ) )**2
            WRITE( *, '( I4, F10.6, 1P,10E11.3 )'   )
     &              I, XMU(I), S1(I), S2(I), S1(I)*CONJG(S2(I)),
     &              I1, I2, 0.5*(I1+I2), (I2-I1)/(I2+I1)
   10    CONTINUE

      END IF


c      IF( PRNT( 2 ) ) THEN

c         WRITE ( *, '(/,A,9X,A,17X,A,17X,A,/,(0P,F7.2, 1P,6E12.3) )' )
c     &           '  Angle', 'S-sub-1', 'T-sub-1', 'T-sub-2',
c     &               0.0,     SFORW,    TFORW(1),  TFORW(2),
c     &              180.,     SBACK,    TBACK(1),  TBACK(2)
c         WRITE ( *, '(/,4(A,1P,E11.4))' )
c     &           ' Efficiency Factors,  extinction:', QEXT,
c     &                              '   scattering:', QSCA,
c     &                              '   absorption:', QEXT-QSCA,
c     &                           '   rad. pressure:', QEXT-GQSC

c      END IF


      RETURN

      END

      SUBROUTINE SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW,
     &                   SBACK, S1, S2, TFORW, TBACK, A, B )

c       Small-particle limit of Mie quantities in totally reflecting
c       limit ( Mie series truncated after 2 terms )

c        A,B       First two Mie coefficients, with numerator and
c                  denominator expanded in powers of XX ( a factor
c                  of XX**3 is missing but is restored before return
c                  to calling program )  ( Ref. 2, p. 1508 )


      IMPLICIT  NONE

c     .. Parameters ..

      REAL*8      TWOTHR, FIVTHR, FIVNIN
      PARAMETER ( TWOTHR = 2. / 3., FIVTHR = 5. / 3., FIVNIN = 5. / 9. )
c     ..
c     .. Scalar Arguments ..

      INTEGER   NUMANG
      REAL*8      GQSC, QEXT, QSCA, XX
      COMPLEX   SBACK, SFORW
c     ..
c     .. Array Arguments ..

      REAL*8      XMU( * )
      COMPLEX   A( * ), B( * ), S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      INTEGER   J
      REAL*8      RTMP
      COMPLEX   CTMP
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC AIMAG, CMPLX, CONJG, REAL
c     ..
c     .. Statement Functions ..

      REAL*8      SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = REAL( CTMP )**2 + AIMAG( CTMP )**2
c     ..

c                                                       ** Eq. R/A.5
      A( 1 ) = CMPLX( 0., TWOTHR*( 1. - 0.2*XX**2 ) ) /
     &         CMPLX( 1. - 0.5*XX**2, TWOTHR*XX**3 )
c                                                       ** Eq. R/A.6
      B( 1 ) = CMPLX( 0., -( 1. - 0.1*XX**2 ) / 3. ) /
     &         CMPLX( 1. + 0.5*XX**2, -XX**3 / 3. )
c                                                       ** Eq. R/A.7,8
      A( 2 ) = CMPLX( 0.,   XX**2 / 30. )
      B( 2 ) = CMPLX( 0., - XX**2 / 45. )
c                                                       ** Eq. R/A.9
      QSCA = 6.* XX**4 *( SQ( A(1) ) + SQ( B(1) ) +
     &                      FIVTHR*( SQ( A(2) ) + SQ( B(2) ) ) )
      QEXT = QSCA
c                                                       ** Eq. R/A.10
      GQSC = 6.* XX**4 * REAL( A(1)*CONJG( A(2) + B(1) ) +
     &                        ( B(1) + FIVNIN*A(2) )*CONJG( B(2) ) )

      RTMP   = 1.5 * XX**3
      SFORW  = RTMP * ( A(1) + B(1) + FIVTHR * ( A(2) + B(2) ) )
      SBACK  = RTMP * ( A(1) - B(1) - FIVTHR * ( A(2) - B(2) ) )
      TFORW( 1 ) = RTMP*( B(1) + FIVTHR * ( 2.*B(2) - A(2) ) )
      TFORW( 2 ) = RTMP*( A(1) + FIVTHR * ( 2.*A(2) - B(2) ) )
      TBACK( 1 ) = RTMP*( B(1) - FIVTHR * ( 2.*B(2) + A(2) ) )
      TBACK( 2 ) = RTMP*( A(1) - FIVTHR * ( 2.*A(2) + B(2) ) )


      DO 10 J = 1, NUMANG
c                                                    ** Eq. R/A.11,12
         S1( J ) = RTMP*( A(1) + B(1)*XMU( J ) +
     &                    FIVTHR*( A(2)*XMU( J ) + 
     &                             B(2)*( 2.*XMU( J )**2 - 1.) ) )
         S2( J ) = RTMP*( B(1) + A(1)*XMU( J ) +
     &                    FIVTHR*( B(2)*XMU( J ) + 
     &                             A(2)*( 2.*XMU( J )**2 - 1.) ) )
   10 CONTINUE

c                                     ** Recover actual Mie coefficients
      A( 1 ) = XX**3 * A(1)
      A( 2 ) = XX**3 * A(2)
      B( 1 ) = XX**3 * B(1)
      B( 2 ) = XX**3 * B(2)

      RETURN
      END

      SUBROUTINE SMALL2( XX, CIOR, CALCQE, NUMANG, XMU, QEXT, QSCA,
     &                   GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK, 
     &                   A, B )

c       Small-particle limit of Mie quantities for general refractive
c       index ( Mie series truncated after 2 terms )

c        A,B       First two Mie coefficients, with numerator and
c                  denominator expanded in powers of XX ( a factor
c                  of XX**3 is missing but is restored before return
c                  to calling program )

c        CIORSQ    Square of refractive index


      IMPLICIT  NONE

c     .. Parameters ..

      REAL*8      TWOTHR, FIVTHR
      PARAMETER  ( TWOTHR = 2./3., FIVTHR = 5./3. )
c     ..
c     .. Scalar Arguments ..

      LOGICAL   CALCQE
      INTEGER   NUMANG
      REAL*8      GQSC, QEXT, QSCA, XX
      COMPLEX   CIOR, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      REAL*8      XMU( * )
      COMPLEX   A( * ), B( * ), S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      INTEGER   J
      REAL*8      RTMP
      COMPLEX   CIORSQ, CTMP
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC AIMAG, CMPLX, CONJG, REAL
c     ..
c     .. Statement Functions ..

      REAL*8      SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = REAL( CTMP )**2 + AIMAG( CTMP )**2
c     ..


      CIORSQ = CIOR**2
      CTMP   = CMPLX( 0., TWOTHR )*( CIORSQ - 1.0 )

c                                           ** Eq. R42a
      A( 1 ) = CTMP*( 1. - 0.1*XX**2 +
     &         ( CIORSQ / 350.+ 1./ 280.)*XX**4 ) /
     &         ( CIORSQ + 2.+ ( 1.- 0.7*CIORSQ )*XX**2 -
     &         ( CIORSQ**2 / 175.- 0.275*CIORSQ + 0.25 )*XX**4 +
     &         XX**3 * CTMP*( 1.- 0.1*XX**2 ) )

c                                           ** Eq. R42b
      B( 1 ) = ( XX**2 / 30.) * CTMP * ( 1.+
     &         ( CIORSQ / 35.- 1./ 14.)*XX**2 ) /
     &         ( 1.- ( CIORSQ / 15.- 1./ 6.)*XX**2 )

c                                           ** Eq. R42c

      A( 2 ) = ( 0.1*XX**2 )*CTMP*( 1.- XX**2 / 14.) /
     &         ( 2.*CIORSQ + 3.- ( CIORSQ / 7.- 0.5 )*XX**2 )

c                                           ** Eq. R40a

      QSCA = (6.*XX**4) * ( SQ( A(1) ) + SQ( B(1) ) +
     &                       FIVTHR*SQ( A(2) ) )

c                                           ** Eq. R40b
      QEXT = QSCA
      IF( CALCQE ) QEXT = 6.*XX * REAL( A(1) + B(1) + FIVTHR*A(2) )

c                                           ** Eq. R40c

      GQSC = (6.*XX**4) * REAL( A(1)*CONJG( A(2) + B(1) ) )

      RTMP   = 1.5 * XX**3
      SFORW  = RTMP*( A(1) + B(1) + FIVTHR*A(2) )
      SBACK  = RTMP*( A(1) - B(1) - FIVTHR*A(2) )
      TFORW( 1 ) = RTMP*( B(1) - FIVTHR*A(2) )
      TFORW( 2 ) = RTMP*( A(1) + 2.*FIVTHR*A(2) )
      TBACK( 1 ) = TFORW(1)
      TBACK( 2 ) = RTMP*( A(1) - 2.*FIVTHR*A(2) )


      DO 10 J = 1, NUMANG
c                                      ** Eq. R40d,e

         S1( J ) = RTMP*( A(1) + ( B(1) + FIVTHR*A(2) )*XMU( J ) )
         S2( J ) = RTMP*( B(1) + A(1)*XMU( J ) +
     &                    FIVTHR * A(2)*( 2.*XMU( J )**2 - 1.) )
   10 CONTINUE

c                                     ** Recover actual Mie coefficients
      A( 1 ) = XX**3 * A(1)
      A( 2 ) = XX**3 * A(2)
      B( 1 ) = XX**3 * B(1)
      B( 2 ) = ( 0., 0.)

      RETURN
      END

      SUBROUTINE TESTMI( COMPAR, XX, CREFIN, MIMCUT, PERFCT, ANYANG,
     &                   NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK,
     &                   S1, S2, TFORW, TBACK )

c         Set up to run test case when  COMPAR = False;  when  = True,
c         compare Mie code test case results with correct answers
c         and abort if even one result is inaccurate.
c         (NoPMOM version)

c         The test case is :  Mie size parameter = 10
c                             refractive index   = 1.5 - 0.1 i
c                             scattering angle = 140 degrees
c                             1 Sekera moment

c         Results for this case may be found among the test cases
c         at the end of reference (1).

c         *** NOTE *** When running on some computers, esp. in single
c         precision, the Accur criterion below may have to be relaxed.
c         However, if Accur must be set larger than 10**-3 for some
c         size parameters, your computer is probably not accurate
c         enough to do Mie computations for those size parameters.

c     Routines called :  ERRMSG, MIPRNT, TSTBAD


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   ANYANG, COMPAR, PERFCT
      INTEGER   NUMANG
      REAL*8      GQSC, MIMCUT, QEXT, QSCA, XX
      COMPLEX   CREFIN, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      REAL*8      XMU( * )
      COMPLEX   S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   ANYSAV, OK, PERSAV
      INTEGER   N, NUMSAV
      REAL*8      ACCUR,  MIMSAV, TESTGQ, TESTQE, TESTQS,
     &          XMUSAV, XXSAV
      REAL*8      CALC, EXACT
      COMPLEX   CRESAV, TESTS1, TESTS2, TESTSB, TESTSF
c     ..
c     .. Local Arrays ..

      LOGICAL   PRNT( 2 )
      COMPLEX   TESTTB( 2 ), TESTTF( 2 )
c     ..
c     .. External Functions ..

      LOGICAL   TSTBAD
      EXTERNAL  TSTBAD
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, MIPRNT
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, REAL
c     ..
c     .. Statement Functions ..

      LOGICAL   WRONG
c     ..
      SAVE      XXSAV, CRESAV, MIMSAV, PERSAV, ANYSAV, NUMSAV, XMUSAV

      DATA   TESTQE / 2.459791 /,
     &       TESTQS / 1.235144 /,
     &       TESTGQ / 1.139235 /,
     &       TESTSF / ( 61.49476, -3.177994 ) /,
     &       TESTSB / ( 1.493434, 0.2963657 ) /,
     &       TESTS1 / ( -0.1548380, -1.128972) /,
     &       TESTS2 / ( 0.05669755, 0.5425681) /,
     &       TESTTF / ( 12.95238, -136.6436 ), ( 48.54238, 133.4656 ) /,
     &       TESTTB / ( 41.88414, -15.57833 ), ( 43.37758, -15.28196 )/
     
      DATA      ACCUR / 1.E-1 /
c     ..
c     .. Statement Function definitions ..

      WRONG( CALC, EXACT ) = ABS( ( CALC - EXACT ) / EXACT ).GT.ACCUR
c     ..


      IF( .NOT.COMPAR ) THEN
c                                   ** Save certain user input values
         XXSAV  = XX
         CRESAV = CREFIN
         MIMSAV = MIMCUT
         PERSAV = PERFCT
         ANYSAV = ANYANG
         NUMSAV = NUMANG
         XMUSAV = XMU( 1 )
c                                   ** Reset input values for test case
         XX     = 10.0
         CREFIN = ( 1.5, -0.1 )
         MIMCUT = 0.0
         PERFCT = .FALSE.
         ANYANG = .TRUE.
         NUMANG = 1
         XMU( 1 ) = -0.7660444

      ELSE
c                                    ** Compare test case results with
c                                    ** correct answers and abort if bad
         OK     = .TRUE.

        IF ( WRONG( QEXT,TESTQE ) )
     &       OK =  TSTBAD( 'QEXT', DBLE(ABS((QEXT - TESTQE) / TESTQE)) )
     
        IF ( WRONG( QSCA,TESTQS ) )
     &       OK =  TSTBAD( 'QSCA', DBLE(ABS((QSCA - TESTQS) / TESTQS)) )
     
        IF ( WRONG( GQSC,TESTGQ ) )
     &       OK =  TSTBAD( 'GQSC', DBLE(ABS((GQSC - TESTGQ) / TESTGQ)) )

         IF( WRONG( DBLE(REAL( SFORW )),DBLE(REAL( TESTSF )) ) .OR.
     &       WRONG( DBLE(AIMAG( SFORW )),DBLE(AIMAG( TESTSF )) ) )
     &    OK = TSTBAD('SFORW', DBLE(ABS( ( SFORW - TESTSF ) / TESTSF )))

         IF( WRONG( DBLE(REAL( SBACK )),DBLE(REAL( TESTSB )) ) .OR.
     &       WRONG( DBLE(AIMAG( SBACK )),DBLE(AIMAG( TESTSB )) ) )
     &   OK = TSTBAD( 'SBACK', DBLE(ABS( ( SBACK - TESTSB ) / TESTSB )))

         IF( WRONG( DBLE(REAL( S1( 1 ) )),DBLE(REAL( TESTS1 )) ) .OR.
     &       WRONG( DBLE(AIMAG( S1( 1 ) )),DBLE(AIMAG( TESTS1 ) )) )
     &    OK = TSTBAD( 'S1', DBLE(ABS( ( S1( 1 ) - TESTS1 ) / TESTS1 )))

         IF( WRONG( DBLE(REAL( S2( 1 ) )),DBLE(REAL( TESTS2 )) ) .OR.
     &       WRONG( DBLE(AIMAG( S2( 1 ) )),DBLE(AIMAG( TESTS2 ) )) )
     &    OK = TSTBAD( 'S2', DBLE(ABS( ( S2( 1 ) - TESTS2 ) / TESTS2 )))


         DO 10 N = 1, 2

         IF ( WRONG(  DBLE(REAL(TFORW(N))), DBLE(REAL(TESTTF(N))) ) .OR.
     &        WRONG( DBLE(AIMAG(TFORW(N))), DBLE(AIMAG(TESTTF(N)) )) )
     &        OK =  TSTBAD( 'TFORW', DBLE(ABS( (TFORW(N) - TESTTF(N)) /
     &                                       TESTTF(N) )) )
         IF ( WRONG(  DBLE(REAL(TBACK(N))), DBLE(REAL(TESTTB(N))) ) .OR.
     &        WRONG( DBLE(AIMAG(TBACK(N))),DBLE(AIMAG(TESTTB(N))) ) )
     &        OK =  TSTBAD( 'TBACK', DBLE(ABS( (TBACK(N) - TESTTB(N)) /
     &                                        TESTTB(N) )) )

   10    CONTINUE


         IF( .NOT.OK ) THEN

            PRNT( 1 ) = .TRUE.
            PRNT( 2 ) = .TRUE.

            CALL  MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                    QSCA, GQSC, SFORW, SBACK, TFORW, TBACK, S1,S2)

            CALL ERRMSG( 'MIEV0 -- Self-test failed', .TRUE.)

         END IF
c                                       ** Restore user input values
         XX     = XXSAV
         CREFIN = CRESAV
         MIMCUT = MIMSAV
         PERFCT = PERSAV
         ANYANG = ANYSAV
         NUMANG = NUMSAV
         XMU( 1 ) = XMUSAV

      END IF

      RETURN
      END

C----------------------------------------------------------------------

      SUBROUTINE ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error
c        after making symbolic dump (machine-specific)

c     .. Scalar Arguments ..

      CHARACTER MESSAG*(*)
      LOGICAL   FATAL
c     ..
c     .. Local Scalars ..

      LOGICAL   MSGLIM
      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

ccccc EXTERNAL  SYMDUMP
c     ..
      SAVE      MAXMSG, NUMMSG, MSGLIM
      DATA      NUMMSG / 0 /,  MAXMSG / 100 /,  MSGLIM /.FALSE./


      IF( FATAL ) THEN

         WRITE( *, '(//,2A,//)' ) ' ****** ERROR *****  ', MESSAG

c                                 ** Example symbolic dump call for Cray
ccccc    CALL SYMDUMP( '-B -c3' )

         STOP

      END IF


      NUMMSG = NUMMSG + 1

      IF( MSGLIM ) RETURN

      IF( NUMMSG.LE.MAXMSG ) THEN

         WRITE( *, '(/,2A,/)' ) ' ****** WARNING *****  ', MESSAG

      ELSE

         WRITE( *, 9000 )
         MSGLIM = .True.
      END IF


      RETURN

 9000 FORMAT( / , / , ' ****** TOO MANY WARNING MESSAGES --  ',
     &      'They will no longer be printed *******', / , / )
      END

      LOGICAL FUNCTION WrtBad( VarNam )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

c     .. Scalar Arguments ..

      CHARACTER VarNam*(*)
c     ..
c     .. Local Scalars ..

      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

      EXTERNAL  ErrMsg
c     ..
      SAVE      NUMMSG, MAXMSG
      DATA      NUMMSG / 0 /, MAXMSG / 50 /


      WrtBad = .TRUE.
      NUMMSG = NUMMSG + 1
      WRITE( *, '(3A)' ) ' ****  Input variable  ', VarNam,
     &                   '  in error  ****'

      IF( NUMMSG.EQ.MAXMSG )
     &    CALL ErrMsg( 'Too many input errors.  Aborting...',.TRUE.)

      RETURN
      END

      LOGICAL FUNCTION WrtDim( DimNam, Minval )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

c     .. Scalar Arguments ..

      CHARACTER DimNam*(*)
      INTEGER   Minval
c     ..

      WRITE( *, '(3A,I7)' ) ' ****  Symbolic dimension  ', DimNam,
     &                      '  should be increased to at least ', Minval
      WrtDim = .TRUE.
      RETURN
      END

      LOGICAL FUNCTION TstBad( VarNam, RelErr )

c       Write name (VarNam) of variable failing self-test and its
c       percent error from the correct value;  return  'FALSE'.

c     .. Scalar Arguments ..

      CHARACTER VarNam*(*)
      REAL*8      RelErr
c     ..

      TstBad = .FALSE.
      WRITE( *, '(/,3A,1P,E11.2,A)' ) ' *** Output variable ', VarNam,
     &   ' differed by ', 100.*RelErr,
     &   ' per cent from correct value.  Self-test failed.'
      RETURN
      END
