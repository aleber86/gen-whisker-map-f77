      PROGRAM WM
      IMPLICIT NONE
      INTEGER DIMX, DIMY
      REAL*8 TANG(3), INITC(3), LOGFY
      REAL*8 MLCE(1024), VELEM, OM2, L2, L1, MEGNOF(1024)
      REAL*8 MAXWD, MINWD, DPI, ETA, LINVMA(1024), LINVMI(1024)
      REAL*8 MAXINV, H
      REAL*8 INISPR, MODUL, NORMA, DOTPRO
      REAL*8 TSTART, TEND
      REAL*8 SHANN, KSENTR
      INTEGER ANGT, ANGX, YWIDE, YOFFST, INDXER
      REAL*8 ANGS, YS, YMN, MEGNO, HW, INVL1, INFORA
c      INTEGER INDBC, INDBCX, INDBP, INDBPX
      INTEGER L1I, INICI, NIT, IIT,SECT, COUNTER, COUNTT
C     FILE ARGUMENTS*********************************
      INTEGER FILEI, NARGS, NUMARG, IFILE
      INTEGER INICD, FORR, FOCC
C     ***********************************************
      PARAMETER(DIMX = 2048, DIMY = 4096, NUMARG=32,
     &INISPR = 1.D-7, NIT=10**7, SECT = 32) 
C     FILE ARGUMENTS********************************
      REAL*8 L1V(NUMARG), L2V(NUMARG), OM2V(NUMARG)
      REAL*8 VV(NUMARG), MUV(NUMARG), HWV(NUMARG)
      REAL*8 ETAV(NUMARG)
C     ***********************************************
      INTEGER  CMAPT(DIMY, DIMX), CMAPX(DIMY, DIMX)
      INTEGER   INFOT(DIMX), INFOX(DIMX)
      REAL MAPOUT(DIMX*DIMY, 3), MAPOUX(DIMX*DIMY, 3)
      CHARACTER*80 CMDARG
      CHARACTER*3 NUMBFI
      REAL*8 OUTPUT(SECT, 10)  
      EXTERNAL MODUL, DOTPRO, INDXER, SHANN

      CALL CPU_TIME(TSTART)
      NARGS = IARGC()
      IF(NARGS .EQ. 0) THEN
         WRITE(*,*) 'Two args needed. file and ini_cond'
         STOP
      END IF
      
      CALL GETARG(1, CMDARG)
C     LOADS FILE DATA PRE-CATCH 
      OPEN(50, FILE=CMDARG, STATUS='OLD', ACTION='READ')
      DO IFILE = 1, NUMARG
         READ(50, *) L1V(IFILE), L2V(IFILE), OM2V(IFILE),
     &        MUV(IFILE), ETAV(IFILE), VV(IFILE), HWV(IFILE)
      END DO
      CLOSE(50)
C     ***************************************************
C     GET NUM. OF INITIAL CONDITIONS
      CALL GETARG(2, CMDARG)
      READ(CMDARG, '(I10)') INICD

      CALL GETARG(3, CMDARG)
      READ(CMDARG, '(I3)') FILEI
      WRITE(NUMBFI, '(I3)') FILEI
      
C     ***************************************************
C     2.PI DEFINITION
      DPI = 8.D0*DATAN(1.D0)
      ANGS = DPI/(1.D0*DIMX)
      YOFFST = DIMY/2
      INFORA = (1.D0*DIMX) / (2.D0 * NIT * INICD * DLOG(1.D0*DIMX))


      COUNTER =  0

      DO L1I = 1, NUMARG
C     FILE DATA LACE FOR EVERY LAYER ON LAMBDA 1
         L1 = L1V(L1I)
         L2 = L2V(L1I)
         VELEM = VV(L1I)
         OM2 = OM2V(L1I)
         ETA = ETAV(L1I)
         HW = HWV(L1I)
         MAXWD = 0.D0
         MINWD = 0.D0
         YS =2.D0 * HW/(1.D0*DIMY)
         INVL1 = 1.D0/L1
         CMAPX = 0
         CMAPT = 0
         INFOT = 0
         INFOX = 0
         MLCE = 0.D0
         KSENTR = 0.D0
         MAPOUT = 0.0
         MAPOUX = 0.0
         MEGNOF = 0.D0

         
         DO INICI = 1, INICD
C     INITIAL CONDITION LACE ITERATION OVER EVERY INI COND
            INITC(1) = RAND()*INISPR
            INITC(2) = RAND()*INISPR
            INITC(3) = RAND()*INISPR
C     TANGENT VECTOR DEFINITION
            TANG(1) = 0.D0
            TANG(2) = 1.D0
            TANG(3) = 0.D0
            MEGNO = 0.D0
            YMN = 0.D0
            LINVMA(INICI) = 0.D0
            LINVMI(INICI) = 0.D0

            DO IIT=1, NIT
C     MAP ITERATION
               INITC(3) = INITC(3) + INVL1*DSIN(INITC(2))
     &              - VELEM*INVL1*DCOS(INITC(1))
               LOGFY = DLOG(DABS(INITC(3)))
               INITC(2) = INITC(2) - L1*LOGFY + ETA
               INITC(1) = INITC(1) - L2*LOGFY + OM2*ETA

               INITC(1) = MODUL(INITC(1), DPI)
               INITC(2) = MODUL(INITC(2), DPI)
C     **************************************************
C     TANGENT VECTOR EVOLUTION
               CALL JACOB(INITC, L1, INVL1, L2, VELEM, TANG, DOTPRO)
               NORMA = DSQRT(TANG(1)**2 + TANG(2)**2 + TANG(3)**2)
               LOGFY = DLOG(NORMA)
               TANG(1) = TANG(1)/NORMA 
               TANG(2) = TANG(2)/NORMA 
               TANG(3) = TANG(3)/NORMA
               MLCE(INICI) = MLCE(INICI) + LOGFY
               MEGNO = MEGNO + 2.D0*LOGFY*(1.D0*IIT)
               YMN = YMN + MEGNO/(1.D0*IIT)
               
               ANGX = INDXER(0, ANGS, INITC(1), DIMX)
               ANGT = INDXER(0, ANGS, INITC(2), DIMX)
               YWIDE = INDXER(YOFFST, YS, INITC(3), DIMY)
               IF(MAXWD .LT. INITC(3)) MAXWD = INITC(3)
               IF(MINWD .GT. INITC(3)) MINWD = INITC(3)
               IF(LINVMA(INICI) .LT. INITC(3)) LINVMA(INICI) = INITC(3)
               IF(LINVMI(INICI) .GT. INITC(3)) LINVMI(INICI) = INITC(3)
               
          CMAPT(YWIDE, ANGT) = CMAPT(YWIDE, ANGT) + 1
          CMAPX(YWIDE, ANGX) = CMAPX(YWIDE, ANGX) + 1     
          INFOT(ANGT) = INFOT(ANGT) + 1
          INFOX(ANGX) = INFOX(ANGX) + 1

          END DO
           
          MEGNOF(INICI) = YMN / (1.D0* NIT)



         END DO
         MAXINV = 0.D0
         H = 0.D0
         DO FORR = 1, INICD
            LINVMA(FORR) = (LINVMA(FORR) - LINVMI(FORR))/2.D0
            H = H + LINVMA(FORR) * MLCE(FORR)
            IF(MAXINV .LT. LINVMA(FORR)) MAXINV = LINVMA(FORR)
         END DO
         H = H/(1.D0*INICD*NIT)

         COUNTER = COUNTER + 1
         OUTPUT(L1I,1) = L1
         OUTPUT(L1I, 2) = (MAXWD - MINWD)/2.D0
         OUTPUT(L1I, 3) = SUM(MLCE)/(1.D0*INICD*NIT)
         OUTPUT(L1I, 4) = 0.D0
         DO FORR = 1, DIMX
         DO FOCC = 1,DIMY
            IF(CMAPT(FOCC, FORR) .GE. 1) OUTPUT(L1I, 4) = OUTPUT(L1I,4) 
     &       + 1.D0
            IF(CMAPX(FOCC, FORR) .GE. 1) OUTPUT(L1I, 5) = OUTPUT(L1I, 5) 
     &       + 1.D0
         END DO
         END DO
         OUTPUT(L1I, 4) = OUTPUT(L1I,4)/(1.D0*DIMY*DIMX)
         OUTPUT(L1I, 5) = OUTPUT(L1I,5)/(1.D0*DIMY*DIMX)
         OUTPUT(L1I, 6) = (1.D0 - 1.D0/(DLOG(1.D0*DIMX))
     &    *SHANN(INICD*NIT, INFOT, DIMX)) / INFORA
         OUTPUT(L1I, 7) = (1.D0 - 1.D0/(DLOG(1.D0*DIMX))
     &    *SHANN(INICD*NIT, INFOX, DIMX)) / INFORA
c         OUTPUT(L1I, 8) = OUTPUT(L1I, 3) * OUTPUT(L1I,2) * L1
         OUTPUT(L1I, 8) = H 
         OUTPUT(L1I, 9) = 2.D0*DABS(SUM(MEGNOF)/(1.D0*INICD) - 2.D0)/NIT
         OUTPUT(L1I, 10) = MAXINV
         IF(MOD(COUNTER, SECT).EQ.0) THEN
         COUNTT = 1
C***************COPY ONE MAP*****************************************         
         DO FORR = 1, DIMX
         DO FOCC = 1,DIMY
            MAPOUT(COUNTT, 2) = YS * (FOCC-YOFFST) * L1
            MAPOUT(COUNTT, 1) = ANGS * FORR
            MAPOUX(COUNTT, 2) = MAPOUT(COUNTT, 2)
            MAPOUX(COUNTT, 1) = MAPOUT(COUNTT, 1)
            MAPOUT(COUNTT, 3) = CMAPT(FOCC, FORR)
            MAPOUX(COUNTT, 3) = CMAPX(FOCC, FORR)
            COUNTT = COUNTT + 1
        END DO
        END DO
C********************************************************************        
        OPEN(30, FILE='out'//NUMBFI//'.dat')
        WRITE(30,*)"#lambda", "y_hw", "<L>", "sigma_(t,y)", 
     &    "sigma_(x,y)", "I_t/I_r", "I_x/I_r", "h/sigma_(t_y)",
     &    "<L>(MEGNO)", "y_hw_lb"
        DO FORR = 1, SECT
            WRITE(30, *) (OUTPUT(FORR, FOCC), FOCC = 1,10)
        END DO
        CLOSE(30)
C***********LAST MAP COPY TO DISK************************************            
C            OPEN(31, FILE='map_t_'//NUMBFI//'.map')
C            COUNTT = 1
C            DO FORR = 1, DIMX
C            DO FOCC = 1, DIMY
C                WRITE(31,*) MAPOUT(COUNTT,1), MAPOUT(COUNTT, 2),
C     &           MAPOUT(COUNTT, 3)
C                COUNTT = COUNTT+1
C            END DO
C            END DO
C                
C            CLOSE(31)
C*********************************************************************            
            CALL CPU_TIME(TEND)
            WRITE(*,*) (TEND - TSTART)
         END IF
         


      END DO
      
      END


      REAL*8 FUNCTION MODUL(ARG, FMOD)
      IMPLICIT NONE
      REAL*8 ARG, FMOD, VALUE
      VALUE = DMOD(ARG, FMOD)
      IF(VALUE .LT. 0.) VALUE = VALUE + FMOD
      MODUL = VALUE
      RETURN 
      END

      SUBROUTINE JACOB(VECIN, L1, INVL1, L2, V,TANGEN, DOTPP)
      IMPLICIT NONE
      REAL*8 VECIN(3), L1, INVL1, L2, V, VECOUT(3)
      REAL*8 SINX, SINT, COSX, COST, ARGS, INVARG
      REAL*8 SIGNED, X1(3), X2(3), X3(3), TANGEN(3)
      REAL*8 DOTPP
      SINT = DSIN(VECIN(2))
      SINX = DSIN(VECIN(1))
      COST = DCOS(VECIN(2))
      COSX = DCOS(VECIN(1))

      ARGS = VECIN(3) - V*INVL1 * COSX + INVL1*SINT
      INVARG = 1.D0/DABS(ARGS)
      SIGNED = DSIGN(1.D0, ARGS)

      X1(1) = -L2*V*INVL1*SIGNED*INVARG*SINX + 1.D0
      X1(2) = -L2*INVL1*SIGNED*INVARG*COST
      X1(3) = -L2*SIGNED*INVARG

      X2(1) = -V*SIGNED*INVARG*SINX
      X2(2) = -SIGNED*INVARG*COST + 1.D0
      X2(3) = SIGNED*L1*INVARG

      X3(1) = V*INVL1*SINX
      X3(2) = INVL1*COST
      X3(3) = 1.D0

      VECOUT(1) = DOTPP(X1, TANGEN, 3)
      VECOUT(2) = DOTPP(X2, TANGEN, 3)
      VECOUT(3) = DOTPP(X3, TANGEN, 3)

      TANGEN(1) = VECOUT(1)
      TANGEN(2) = VECOUT(2)
      TANGEN(3) = VECOUT(3)

      RETURN
      END

      REAL*8 FUNCTION DOTPRO(VEC1, VEC2, N)
      IMPLICIT NONE
      INTEGER N, I
      REAL*8 VEC1(N), VEC2(N)
      DOTPRO = 0.D0

      DO I=1, N
         DOTPRO = DOTPRO + VEC1(I) * VEC2(I)
      END DO

      RETURN
      END

      INTEGER FUNCTION INDXER(OFFSET, SPREAD, COORD, ARRAYD)
      IMPLICIT NONE
      INTEGER OFFSET, ARRAYD, COUNT, INDEX
      REAL*8 SPREAD, COORD

      COUNT = INT(COORD/SPREAD)
      IF((COUNT + OFFSET) .LT. 1) THEN
         INDEX = 1

      ELSE IF((COUNT + OFFSET) .GE. OFFSET + ARRAYD) THEN
         INDEX = OFFSET + ARRAYD - 1

      ELSE
         INDEX = COUNT + OFFSET
      END IF
      INDXER = INDEX
      RETURN
      END
      
         
      REAL*8 FUNCTION SHANN(N, CCOUNT, CARDINA)
      IMPLICIT NONE
      INTEGER CARDINA
      INTEGER N, IT
      INTEGER CCOUNT(CARDINA)
      REAL*8 SUMM
      SUMM = 0.D0

      DO IT = 1, CARDINA
      IF (CCOUNT(IT) .GT. 0) THEN
        SUMM = SUMM + CCOUNT(IT)*DLOG(1.D0*CCOUNT(IT))
      END IF
      END DO
      SHANN = DLOG(1.D0*N) - 1.D0/(1.D0*N)*SUMM
      RETURN 
      END
