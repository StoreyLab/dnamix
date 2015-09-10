C Last modified August 14, 1998

      PROGRAM DNAMIX

***************************************************************
*                                                             *
* This program will perform calculations for the methods and  * 
* formulas presented in:                                      *
* Curran, J.M., Triggs, C.M., Buckleton, J., and B.S. Weir.   *
* 1998. Interpreting DNA mixtures in structured populations,  *
* preprint.                                                   *
*                                                             *
* Questions and/or comments ahould be sent to:                *
* storey@statgen.ncsu.edu  -or-                               *
* weir@stat.ncsu.edu                                          *
*                                                             *
***************************************************************


      INTEGER ALLELE, CONT_L, CONT_U, ANS, 
     +        SMRY, NBAS, CONT_L2, CONT_U2,
     +        K_CONT, V_CONT, K_ALE(2), T_SIZ,
     +        K_CONT2, V_CONT2, LN
      DOUBLE PRECISION CHECK
      CHARACTER*8, FIL, LOCUS
      CHARACTER*1, PAWS
      INTEGER I,J,K,L,A, HETER
      INTEGER, DIMENSION(:), ALLOCATABLE :: INDT, T, V
      CHARACTER*8, DIMENSION(:), ALLOCATABLE :: NAME,
     +                                          BASE
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::
     +                                THETA
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::
     +    FREQ, LKH_NM, LKH_DN


*********************************************************************
* Gets the data, i.e. allele names and frequencies present in the 
* evidence sample and the databases to be used to calculate
* likelihoods
*********************************************************************

      CALL BANNER()

      ANS = 0

 601  CONTINUE
      PRINT *
      WRITE(*, '(1X,A)', ADVANCE='YES') 
     +'Enter the number of alleles in the mixture: '
      READ *, ALLELE
      WRITE(*,*)
      IF (ALLELE .LT. 1) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 601
      END IF

      A = ALLELE

      ALLOCATE (INDT(A), T(A), V(A), NAME(A)) 

      DO K=1, ALLELE
         WRITE(*, '(1X, A, I2, A)', ADVANCE='YES') 
     +   'Enter the name of allele', K, ': '
         READ(*,100,ADVANCE='YES') NAME(K)
      END DO
 100  FORMAT(A8)
      WRITE(*,*)

 643  CONTINUE
      WRITE(*, '(1X, A)', ADVANCE='YES')
     +   'Enter the number of databases to be used: '
      READ *, NBAS
      IF (NBAS .LT. 1) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         PRINT *
         GO TO 643
      END IF

      PRINT *
      
      ALLOCATE (BASE(NBAS), FREQ(A, NBAS), THETA(NBAS))

      DO 250 J = 1, NBAS
         WRITE(*, '(1X,A,I2,A)', ADVANCE='YES') 
     +      'Enter the name for database number', J, ': '
         READ '(A)', BASE(J)
      DO 10 I= 1 , ALLELE
         CALL WORDLN(LN, NAME(I))
         WRITE(*, '(1X, 3A)', ADVANCE='YES') 
     +'Enter the frequency of allele ',NAME(I)(:LN), ': '
         READ *, FREQ(I, J)
 10   CONTINUE

      CHECK = 0.0
      DO 120 K=1, ALLELE
         CHECK = CHECK + FREQ(K, J)
 120  CONTINUE
      IF ( NINT(1000*CHECK) .GE. 1005) THEN
         PRINT *, 'WARNING: The sum of the allele ',
     +            'frequencies may be greater than 1.00!'
         PRINT *
      END IF

 793  CONTINUE
      WRITE(*,'(1X,A,I2,A)',ADVANCE='YES') 
     +'Enter the coancestry coefficient for database number', 
     +J, ': '
      READ *, THETA(J)
      IF((NINT(10000*THETA(J)) .GE. 10000) .OR.
     +   (NINT(10000*THETA(J)) .LT. 0)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 793
      END IF
         

      WRITE(*,*)
      
 250  CONTINUE


 650  CONTINUE
      WRITE(*, '(1X, A)', ADVANCE='YES')
     +     'Do you want a summary file? ( No=0, Yes=1 ): '
      READ *, SMRY
      IF ((SMRY .NE. 0) .AND. (SMRY .NE. 1)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 650
      END IF
      IF (SMRY .EQ. 1) THEN
         PRINT *
         WRITE(*, '(1X, A)', ADVANCE = 'YES')
     +   'Enter the name of the file: '
         READ '(A8)', FIL
         OPEN (UNIT=3, FILE=FIL, STATUS='UNKNOWN')
         PRINT *
         WRITE(*, '(1X, A)', ADVANCE='YES') 
     +   'Enter a title for the summary (such as a locus name): '
         READ '(A)', LOCUS
      END IF

      ANS = 0

      IF (SMRY .EQ. 1) THEN
      IF (ANS .EQ. 0) THEN
      WRITE (3,*)
      WRITE (3,*) '*** Output for DNAMIX version 2 ***'
      WRITE (3,*)
      WRITE (3,*) 'Locus: ', LOCUS
      DO 876 L=1, NBAS, 3
      WRITE (3,*)
      WRITE (3,*) ('Database: ', BASE(J),'      ', J=L,MIN(L+2,NBAS)) 
      WRITE (3,*) ('Allele        Freq      ', I=L,MIN(L+2,NBAS))
      WRITE (3,*) ('------------------      ', K=L,MIN(L+2,NBAS))
      DO 34 A = 1, ALLELE
      WRITE(3,521) (NAME(A), ' ',FREQ(A,I),'      ',I=L,MIN(L+2,NBAS))
 521  FORMAT (1X, 12(A8, A, F9.6,A))
 34   CONTINUE
      WRITE (3,*) ('------------------      ', K=L,MIN(L+2,NBAS))
      WRITE (3,47) ('THETA = ',THETA(K),'      ',K=L,MIN(L+2,NBAS))
 47   FORMAT (1X, 12(A8,F10.6,A))
 876  CONTINUE
      WRITE (3,'(//1X, A,A/)') '*************************************',
     +       '***************************************'
      END IF
      END IF


*********************************************************************
* Gets the alleles whose sources are un/known and the possible number 
* of unknown contributors for the numaerator of the likelihood ratio
*********************************************************************


 999  CONTINUE
      PRINT *
      PRINT *, '**********************************************',
     +         '**************************'
      PRINT *, '    THE FOLLOWING INPUT IS FOR THE NUMERATOR',
     +         ' OF THE LIKELIHOOD RATIO  '
      PRINT *, '**********************************************',
     +         '**************************'

 5    FORMAT (1X, A)
 6    FORMAT (1X, 2A)

 603  CONTINUE
      PRINT *
      WRITE(*, 5, ADVANCE ='YES')  
     + 'Enter the number of known contributors: '
      READ *, K_CONT

      IF(K_CONT .LT. 0) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 603
      END IF

      IF (SMRY .EQ. 1) THEN
         WRITE(3,*)
         WRITE(3,*) ' ----------------------------------------'
         WRITE(3,*) '| NUMERATOR HYPOTHESIS AND PROBABILITIES |'
         WRITE(3,*) ' ----------------------------------------'
         WRITE(3,*)
         WRITE(3,*) 'Number of known contributors: ', K_CONT
      IF (K_CONT .NE. 0) THEN
         WRITE(3,*) 'Individual   Genotype'
         WRITE(3,*) '-------------------------------'
      END IF
      END IF

      DO 66 I = 1, ALLELE
         INDT(I) = 0
         T(I) = 0
         V(I) = 0
 66   CONTINUE
      T_SIZ = 0
      HETER = 0

      DO K=1, K_CONT
      PRINT *
      PRINT '(1X,2A,I2,A)','Enter the numbers corresponding to the ',
     +         'genotype for known contributor number', K, '.'
      
      DO 12 I= 1, ALLELE
         PRINT 101, I, ' = ', NAME(I)
 12   CONTINUE
 101  FORMAT (1X, I2, 2A)

      DO 13 I= 1, 2
 604     CONTINUE
         WRITE(*, 5, ADVANCE='NO')
     +    'Enter the number corresponding to the'
         IF(I .EQ. 1) WRITE(*, 5, ADVANCE='YES')
     +    'first allele: '
         IF(I .EQ. 2) WRITE(*, 5, ADVANCE='YES')
     +    'second allele: '
         READ *, K_ALE(I)
         IF ((K_ALE(I) .GT. ALLELE) .OR. (K_ALE(I) .LT. 1)) THEN
            PRINT *, 'NOT A VALID ENTRY!'
            GO TO 604
         END IF
         INDT(K_ALE(I)) = 1
         T(K_ALE(I)) = T(K_ALE(I)) + 1
 13   CONTINUE

      IF (SMRY .EQ. 1) THEN
         CALL WORDLN(LN, NAME(K_ALE(1)))
         WRITE(3,744) K, ' ', NAME(K_ALE(1))(:LN), 
     +                   ', ', NAME(K_ALE(2))
 744  FORMAT (1X, I2, A11, 3A) 
      END IF

      IF (K_ALE(1) .NE. K_ALE(2)) THEN
         HETER = HETER + 1
      END IF

      END DO

      IF ((SMRY .EQ. 1) .AND. (K_CONT .NE. 0)) THEN
         WRITE(3,*) '-------------------------------'
         WRITE(3,*)
      END IF

      DO I=1, ALLELE
         T_SIZ = T_SIZ + INDT(I)
      END DO

 529  CONTINUE
      PRINT *
      WRITE(*, 6, ADVANCE='YES')
     + 'Enter the number of individuals known NOT to have',
     + ' contributed to the sample: '
      READ *, V_CONT

      IF(V_CONT .LT. 0) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 529
      END IF

      IF (SMRY .EQ. 1) THEN
         WRITE(3, *)
         WRITE(3, *) 'Number of individuals known NOT to have',
     +               ' contributed: ', V_CONT
      IF (V_CONT .NE. 0) THEN
         WRITE(3,*) 'Individual   Genotype'
         WRITE(3,*) '-------------------------------'
      END IF
      END IF

      DO K=1, V_CONT
      PRINT *
      PRINT '(1X,2A,I2,A)','Enter the numbers corresponding to the ',
     +         'genotype for known non-contributor number', K, '.'
      
      DO I= 1, ALLELE
         PRINT 101, I, ' = ', NAME(I)
      END DO
      
      DO 531 I= 1, 2
 530     CONTINUE
         WRITE(*, 5, ADVANCE='NO')
     +    'Enter the number corresponding to the'
         IF(I .EQ. 1) WRITE(*, 5, ADVANCE='YES')
     +    'first allele: '
         IF(I .EQ. 2) WRITE(*, 5, ADVANCE='YES')
     +    'second allele: '
         READ *, K_ALE(I)
         IF ((K_ALE(I) .GT. ALLELE) .OR. (K_ALE(I) .LT. 1)) THEN
            PRINT *, 'NOT A VALID ENTRY!'
            GO TO 530
         END IF
         V(K_ALE(I)) = V(K_ALE(I)) + 1
 531  CONTINUE

      IF (SMRY .EQ. 1) THEN
         CALL WORDLN(LN, NAME(K_ALE(1)))
         WRITE(3,744) K, ' ', NAME(K_ALE(1))(:LN), 
     +                   ', ', NAME(K_ALE(2))
      END IF

      IF (K_ALE(1) .NE. K_ALE(2)) THEN
         HETER = HETER + 1
      END IF

      END DO

      IF ((SMRY .EQ. 1) .AND. (V_CONT .NE. 0)) THEN
         WRITE(3,*) '-------------------------------'
         WRITE(3,*)
      END IF

 605  CONTINUE
      PRINT *
      WRITE(*, 6, ADVANCE ='YES')
     + 'Enter the lower bound on the number of ',
     + 'unknown contributors: '
      READ *, CONT_L

      IF ((CONT_L.LT.0).OR.((2*CONT_L-ALLELE+T_SIZ).LT.0)
     +    ) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 605
      END IF

 705  CONTINUE
      WRITE(*, 6, ADVANCE='YES') 
     + 'Enter the upper bound on the number of ',
     +         'unknown contributors: '
      READ *, CONT_U

      IF (CONT_U .LT. CONT_L) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 705
      END IF

      IF ((2*CONT_U + T_SIZ) .GT. 13) THEN
         PRINT *
         PRINT *, 'YOU MUST ENTER SMALLER BOUNDS ON THE NUMBER',
     + ' OF UNKNOWN CONTRIBUTORS DUE TO MEMORY CONSTRAINTS!'
         GO TO 605
      END IF

      ALLOCATE (LKH_NM(NBAS, CONT_U+1))

      DO 262 J=CONT_L , CONT_U 

      PRINT *
      PRINT '(1X, A, I2)', 'Unknown Contributors: ', J
      PRINT *, 'Database       Numerator Probability'
      PRINT *, '------------------------------------'

      IF(SMRY .EQ. 1) THEN
         WRITE(3,*)
         WRITE(3,'(1X, A, I2)') 'Unknown Contributors: ', J
         WRITE(3,*) 'Database       Numerator Probability'
         WRITE(3,*) '------------------------------------'
      END IF

      DO 252 I=1, NBAS

      CALL CALC (FREQ, ALLELE, INDT, K_CONT, J, 
     +           THETA(I), LKH_NM(I, J+1), I, T_SIZ, T, 
     +           V, V_CONT, NBAS, HETER)
     
      IF (NINT(1000*LKH_NM(I,J+1)) .GT. 1000) THEN
         LKH_NM(I,J+1) = 1.00000
      END IF

      WRITE(*,800) BASE(I), LKH_NM(I,J+1)
 800  FORMAT (1X, A8, ' ', E22.6)

      IF(SMRY .EQ. 1) THEN
         WRITE(3,800) BASE(I), LKH_NM(I,J+1)
      END IF

 252  CONTINUE
      PRINT *
      IF(SMRY .EQ. 1) THEN
         WRITE(3,*) '------------------------------------'
      END IF
 262  CONTINUE

      PRINT *
      PRINT *, 'Press RETURN to continue.'
      READ '(A)', PAWS


*********************************************************************
* Gets the alleles whose sources are un/known and the possible number 
* of unknown contributors for the denominator of the likelihood ratio
*********************************************************************


      PRINT *
      PRINT *, '**********************************************',
     +         '**************************'
      PRINT *, '   THE FOLLOWING INPUT IS FOR THE DENOMINATOR ',
     +         'OF THE LIKELIHOOD RATIO   '
      PRINT *, '**********************************************',
     +         '**************************'

 606  CONTINUE
      PRINT *
      WRITE(*, 5, ADVANCE = 'YES') 
     + 'Enter the number of known contributors: '
      READ *, K_CONT2

      IF (K_CONT2 .LT. 0) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GOTO 606
      END IF

      DO 67 I = 1, ALLELE
         INDT(I) = 0
         T(I) = 0
         V(I) = 0
 67   CONTINUE
      T_SIZ = 0
      HETER = 0

      IF (SMRY .EQ. 1) THEN
         WRITE(3,*)
         WRITE(3,*)
         WRITE(3,*) ' ------------------------------------------'
         WRITE(3,*) '| DENOMINATOR HYPOTHESIS AND PROBABILITIES |'
         WRITE(3,*) ' ------------------------------------------'
         WRITE(3,*)
         WRITE(3,*) 'Number of known contributors: ', K_CONT2
      IF (K_CONT2 .NE. 0) THEN
         WRITE(3,*) 'Individual   Genotype'
         WRITE(3,*) '-------------------------------'
      END IF
      END IF

      DO K=1, K_CONT2
      PRINT *
      PRINT '(1X,2A,I2,A)','Enter the numbers corresponding to the ',
     +         'genotype for known contributor number', K, '.'
      
      DO I= 1, ALLELE
         PRINT 101, I, ' = ', NAME(I)
      END DO

      DO 15 I= 1, 2
 607  CONTINUE
         WRITE(*, 5, ADVANCE='NO')
     +    'Enter the number corresponding to the'
         IF(I .EQ. 1) WRITE(*, 5, ADVANCE='YES')
     +    'first allele: '
         IF(I .EQ. 2) WRITE(*, 5, ADVANCE='YES')
     +    'second allele: '
         READ *, K_ALE(I)
         IF ((K_ALE(I) .GT. ALLELE) .OR. (K_ALE(I) .LT. 1)) THEN
            PRINT *, 'NOT A VALID ENTRY!'
            GO TO 607
         END IF
         INDT(K_ALE(I)) = 1
         T(K_ALE(I)) = T(K_ALE(I)) + 1
 15   CONTINUE

      IF (SMRY .EQ. 1) THEN
         CALL WORDLN(LN, NAME(K_ALE(1)))
         WRITE(3,744) K, ' ', NAME(K_ALE(1))(:LN), 
     +                   ', ', NAME(K_ALE(2))
      END IF

      IF (K_ALE(1) .NE. K_ALE(2)) THEN
         HETER = HETER + 1
      END IF

      END DO

      DO I=1, ALLELE
         T_SIZ = T_SIZ + INDT(I)
      END DO

      IF ((SMRY .EQ. 1) .AND. (K_CONT2 .NE. 0)) THEN
         WRITE(3,*) '-------------------------------'
         WRITE(3,*)
      END IF

 991  CONTINUE
      PRINT *
      WRITE(*, 6, ADVANCE='YES')
     + 'Enter the number of individuals known NOT to have',
     +         ' contributed to the sample: '
      READ *, V_CONT2

      IF(V_CONT2 .LT. 0) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 991
      END IF

      IF (SMRY .EQ. 1) THEN
         WRITE(3,*)
         WRITE(3, *) 'Number of individuals known NOT to have',
     +               ' contributed: ', V_CONT2
      IF (V_CONT2 .NE. 0) THEN
         WRITE(3,*) 'Individual   Genotype'
         WRITE(3,*) '-------------------------------'
      END IF
      END IF

      DO K=1, V_CONT2
      PRINT *
      PRINT '(1X,2A,I2,A)','Enter the numbers corresponding to the ',
     +         'genotype for known non-contributor number', K, '.'
      
      DO I= 1, ALLELE
         PRINT 101, I, ' = ', NAME(I)
      END DO
      
      DO 640 I= 1, 2
 645     CONTINUE
         WRITE(*, 5, ADVANCE='NO')
     +    'Enter the number corresponding to the'
         IF(I .EQ. 1) WRITE(*, 5, ADVANCE='YES')
     +    'first allele: '
         IF(I .EQ. 2) WRITE(*, 5, ADVANCE='YES')
     +    'second allele: '
         READ *, K_ALE(I)
         IF ((K_ALE(I) .GT. ALLELE) .OR. (K_ALE(I) .LT. 1)) THEN
            PRINT *, 'NOT A VALID ENTRY!'
            GO TO 645
         END IF
         V(K_ALE(I)) = V(K_ALE(I)) + 1
 640  CONTINUE

      IF (SMRY .EQ. 1) THEN
         CALL WORDLN(LN, NAME(K_ALE(1)))
         WRITE(3,744) K, ' ', NAME(K_ALE(1))(:LN), 
     +                   ', ', NAME(K_ALE(2))
      END IF

      IF (K_ALE(1) .NE. K_ALE(2)) THEN
         HETER = HETER + 1
      END IF

      END DO

      IF ((SMRY .EQ. 1) .AND. (V_CONT2 .NE. 0)) THEN
         WRITE(3,*) '-------------------------------'
         WRITE(3,*)
      END IF

 608  CONTINUE
      PRINT *
      WRITE(*, 6, ADVANCE='YES')
     + 'Enter the lower bound on the number of ',
     +         'unknown contributors: '
      READ *, CONT_L2

      IF ((CONT_L2.LT.0).OR.((2*CONT_L2-ALLELE+T_SIZ).LT.0)
     +    ) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 608
      END IF

 706  CONTINUE
      WRITE(*, 6, ADVANCE='YES') 
     + 'Enter the upper bound on the number of ',
     +         'unknown contributors: '
      READ *, CONT_U2

      IF (CONT_U2 .LT. CONT_L2) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 706
      END IF

      IF ((2*CONT_U2 + T_SIZ) .GT. 13) THEN
         PRINT *
         PRINT *, 'YOU MUST ENTER SMALLER BOUNDS ON THE NUMBER',
     + ' OF UNKNOWN CONTRIBUTORS DUE TO MEMORY CONSTRAINTS!'
         GO TO 608
      END IF

      ALLOCATE (LKH_DN(NBAS, CONT_U2+1))

      DO 263 J=CONT_L2 , CONT_U2 

      PRINT *
      PRINT '(1X, A, I2)', 'Unknown Contributors: ', J
      PRINT *, 'Database     Denominator Probability'
      PRINT *, '------------------------------------'

      IF(SMRY .EQ. 1) THEN
         WRITE(3,*)
         WRITE(3,'(1X, A, I2)') 'Unknown Contributors: ', J
         WRITE(3,*) 'Database       Denominator Probability'
         WRITE(3,*) '--------------------------------------'
      END IF

      DO 254 I=1, NBAS

      CALL CALC (FREQ, ALLELE, INDT, K_CONT2, J, THETA(I),
     +           LKH_DN(I,J+1), I, T_SIZ, T, V, 
     +           V_CONT2, NBAS, HETER)

      IF (NINT(1000*LKH_DN(I,J+1)) .GT. 1000) THEN
         LKH_DN(I,J+1) = 1.00000
      END IF

      WRITE(*,800) BASE(I), LKH_DN(I,J+1)

      IF(SMRY .EQ. 1) THEN
         WRITE(3,800) BASE(I), LKH_DN(I,J+1)
      END IF

 254  CONTINUE

      PRINT *
      IF (SMRY .EQ. 1) THEN
         WRITE(3,*) '--------------------------------------'
      END IF
 263  CONTINUE

      PRINT *
      PRINT *, 'Press RETURN to continue.'
      READ '(A)', PAWS


*********************************************************************
* Displays the likelihood ratios according to databases and number of 
* unknown contributors in the numerator and denominator
*********************************************************************

      PRINT *
      PRINT *, '**********************************************',
     +         '**************************'
      PRINT *, '       THE FOLLOWING ARE THE LIKELIHOOD RATIOS',
     +         ' FOR EACH DATABASE       '
      PRINT *, '**********************************************',
     +         '**************************'
      
      PRINT *
      PRINT *, 'Press RETURN to see the results for the next ',
     +         'database or to continue.'
      PRINT *
      PRINT *
      PRINT *, '           Numerator Unknown    Denominator Unknown'
      PRINT *, 'Database      Contributors         Contributors    ',
     +         '     Likelihood Ratio'
      PRINT *, '---------------------------------------------------',
     +         '---------------------'
      DO 751 I = 1, NBAS         
      DO 752 J = CONT_L, CONT_U
      DO 753 K = CONT_L2, CONT_U2
         PRINT '(1X,A8,I12,I21,F31.4)', BASE(I), J, K, 
     +         LKH_NM(I,J+1)/LKH_DN(I,K+1)
 753  CONTINUE
 752  CONTINUE
      READ '(A)', PAWS
 751  CONTINUE

*********************************************************************
* Writes a summary of all the input and output performed so far
*********************************************************************

      IF (SMRY .EQ. 1) THEN
      WRITE(3,*)
      WRITE(3,*)
      WRITE(3,*) ' -------------------'
      WRITE(3,*) '| LIKELIHOOD RATIOS |'
      WRITE(3,*) ' -------------------'
      WRITE(3,*)
      WRITE(3,*)'           Numerator Unknown    Denominator Unknown'
      WRITE(3,*)'Database      Contributors         Contributors    ',
     +         '     Likelihood Ratio'
      WRITE(3,*)'---------------------------------------------------',
     +         '---------------------'
      DO 851 I = 1, NBAS         
      DO 852 J = CONT_L, CONT_U
      DO 853 K = CONT_L2, CONT_U2
         WRITE(3,'(1X,A8,I12,I21,F31.4)') BASE(I), J, K, 
     +         LKH_NM(I,J+1)/LKH_DN(I,K+1)
 853  CONTINUE
 852  CONTINUE
      IF (I .NE. NBAS) THEN
      WRITE(3,*)
      ELSE
      WRITE(3,*)'---------------------------------------------------',
     +         '---------------------'
      END IF
 851  CONTINUE
      WRITE(3,*)
      WRITE (3,'(//1X, A,A/)') '*************************************',
     +       '***************************************'
      WRITE(3,*)
      END IF

*********************************************************************
* Prompts the user if he/she wants to calculate another likelihood
* ratio. Returns to the numerator input if 'yes'.
*********************************************************************

      DEALLOCATE (LKH_NM, LKH_DN)
      
 651  CONTINUE
      PRINT *
      WRITE(*, 6, ADVANCE='YES') 
     + 'Do you want to compute another likelihood',
     +         ' ratio? ( No=0, Yes=1 ): '
      READ *, ANS
      IF ((ANS .NE. 0) .AND. (ANS .NE. 1)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 651
      END IF
      IF (ANS .EQ. 1) GO TO 999

      PRINT *
      WRITE(*,431)'Program written by John Storey on August 1, 1998.'
      WRITE(*,431)'  Questions and/or comments should be sent to:   '
      WRITE(*,431)'             storey@statgen.ncsu.edu             ' 
 431  FORMAT (11(' '), A)
      PRINT *

      DEALLOCATE (INDT, T, V, NAME, 
     +            BASE, THETA, FREQ)

      END

*********************************************************************
*Calculates the likelihood ratio for the given data
*********************************************************************

      SUBROUTINE CALC (FQ, ALLELE, INDT, K_CONT, U_CONT, 
     +                 THETA, ANS, NDBAS, T_SIZ, T, 
     +                 V, V_CONT, NBAS, HETER)

      INTEGER ALLELE, NBAS, HETER
      INTEGER INDT(ALLELE), K_CONT, U_CONT, NDBAS, 
     +        T_SIZ, RM, T(ALLELE), 
     +        UVALS, I, J, K,
     +        FACTOR, V(ALLELE), V_CONT 
      DOUBLE PRECISION FQ(ALLELE, NBAS), THETA, ANS
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PSUM
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INDEX
      INTEGER, DIMENSION(:), ALLOCATABLE :: U
      INTEGER, DIMENSION(:), ALLOCATABLE :: R
      INTRINSIC DEXP, DLOG

      ANS = 0

      ALLOCATE (R(ALLELE), U(ALLELE))

      RM = 2*U_CONT - (ALLELE - T_SIZ)

      UVALS = (FACTOR(ALLELE+RM-1)/(FACTOR(ALLELE-1)*
     +         FACTOR(RM)))

      ALLOCATE (INDEX(UVALS,ALLELE), PSUM(UVALS))

      CALL LIST(RM, ALLELE, INDEX, UVALS)

      DO I=1, UVALS

      DO J=1, ALLELE
         R(J) = INDEX(I,J) 
         U(J) = R(J) + (1 - INDT(J))
      END DO
         
      PSUM(I) = 0
         
      PSUM(I) = DLOG(DBLE(FACTOR(2*U_CONT))) + 
     +          DLOG(DBLE(2**HETER))

      DO J=1, ALLELE
         DO K=0, (T(J)+U(J)+V(J)-1)
            PSUM(I) = PSUM(I) + DLOG(DBLE((1-THETA)*FQ(J,NDBAS) +
     +                              K*THETA))
         END DO
      END DO

      DO J=1, ALLELE
         PSUM(I) = PSUM(I) - DLOG(DBLE(FACTOR(U(J))))
      END DO

      DO J=0, (2*U_CONT+2*K_CONT+2*V_CONT-1)
         PSUM(I) = PSUM(I) - DLOG(DBLE(1-THETA+J*THETA))
      END DO

      PSUM(I) = DEXP(PSUM(I))

      END DO

      DO I=1, UVALS
         ANS = ANS + PSUM(I)
      END DO

      DEALLOCATE (INDEX, PSUM)
      DEALLOCATE (R, U)

      END

*********************************************************************
* Lists all possible ways to have k objects in n urns
*********************************************************************

      SUBROUTINE LIST(RM, AL, IN, UV)

      INTEGER AL, RM, UV, J, CURR, K
      INTEGER IN(UV,AL)
      INTEGER, DIMENSION(:), ALLOCATABLE :: INDX

      ALLOCATE (INDX(AL))


      DO J=1, AL
         INDX(J) = 0
      END DO 

      CURR = 1
      K = 1

      CALL MAKENEW(AL, INDX, RM, CURR, IN, UV, AL, K)

      DEALLOCATE (INDX)

      END


      RECURSIVE SUBROUTINE MAKENEW(DEPTH, INDX, UP, CURR,
     +                             IN, UV, AL, K)

      INTEGER DEPTH, UP, SUM, CURR, I
      INTEGER UV, AL, K
      INTEGER INDX(AL), IN(UV,AL)

      SUM = 0
 
      DO I=1, CURR-1
         SUM = SUM + INDX(I)
      END DO 
 
      IF(CURR .LT. DEPTH) THEN
         INDX(CURR) = 0
         DO WHILE(INDX(CURR) .LE. UP-SUM)
            CALL MAKENEW(DEPTH, INDX, UP, CURR+1,
     +                   IN, UV, AL, K)
            INDX(CURR) = INDX(CURR) + 1
         END DO
      ELSE IF(CURR .EQ. DEPTH) THEN
         INDX(CURR) = UP-SUM
         CALL MAKENEW(DEPTH, INDX, UP, CURR+1,
     +                IN, UV, AL, K)
      ELSE
         CALL SHOW(INDX, DEPTH, IN, UV, AL, K) 
      END IF

      END



      SUBROUTINE SHOW(INDX, DEPTH, IN, UV, AL, K)

      INTEGER DEPTH, J, UV, AL, K
      INTEGER INDX(AL), IN(UV, AL)

      DO J=1, DEPTH
         IN(K,J) = INDX(J)
      END DO
      
      K = K+1

      END

      
*********************************************************************
* Factorial function
*********************************************************************

      FUNCTION FACTOR (N)

      INTEGER FACTOR, N, I

      FACTOR = 1
      
      DO 10 I = 2, N
         FACTOR = FACTOR * I
 10   CONTINUE

      END

*********************************************************************
* Displays a banner giving a reference to the paper from which the 
* formulae come
*********************************************************************

      SUBROUTINE BANNER()

      CHARACTER*1, PAWS

      PRINT *
      PRINT *
      PRINT *, '**********************************************',
     +         '**************************'
      PRINT *, '*                                             ',
     +         '                         *'
      PRINT *, '*                           DNAMIX version 2  ',
     +         '                         *'
      PRINT *, '*                                             ',
     +         '                         *'
      PRINT *, '**********************************************',
     +         '**************************'
      PRINT *, '*                                             ',
     +         '                         *'
      PRINT *, '*    This program performs calculations for the m',
     +         'ethods and formulas   *'
      PRINT *, '*    presented in:                            ',
     +         '                         *'
      PRINT *, '*                                             ',
     +         '                         *'   
      PRINT *, '*    Curran, J.M., Triggs, C.M., Buckleton, J., a',
     +         'nd B.S. Weir. 1998.   *'
      PRINT *, '*    Interpreting DNA mixtures in structured popu',
     +         'lations, preprint.    *'
      PRINT *, '*                                             ',
     +         '                         *'
      PRINT *, '**********************************************',
     +         '**************************'
      PRINT *
      PRINT *, 'Press RETURN to continue.'
      READ '(A)', PAWS
      END

*********************************************************************
* Finds the length of the caharcter string up to the first blank 
* character.
*********************************************************************

      SUBROUTINE WORDLN(LN, CH)

      INTEGER LN
      CHARACTER *(*) CH

      IF(INDEX(CH, ' ') .EQ. 0) THEN
         LN = 8
      ELSE
         LN = MAX(1, INDEX(CH, ' ')-1)
      END IF

      END


*********************************************************************
* End of program
*********************************************************************
