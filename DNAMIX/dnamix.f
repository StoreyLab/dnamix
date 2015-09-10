      PROGRAM DNAMIX

**************************************************************
*                                                            *
* This program will perform calculations for the methods and * 
* formulas presented in:                                     *
* Weir BS, Triggs CM, Starling L, Stowell LI, Walsh KAJ,     *
* Buckleton J. Interpreting DNA Mixtures. J Forensic Sci     *
* 1997; 42(2):213-222.                                       *
*                                                            *
* Questions and/or comments ahould be sent to:               *
* storey@statgen.ncsu.edu                                    *
*                                                            *
**************************************************************

      INTEGER MAX, MAX2
      PARAMETER (MAX=12, MAX2=924)

      INTEGER ALLELE, FORM, UNK_AL, CONT_L, CONT_U, ANS, 
     +        US(MAX), SMRY, NBAS, CONT_L2, CONT_U2, US1(MAX),
     +        UNK_AL2
      DOUBLE PRECISION FREQ(MAX,MAX), FQ(MAX, MAX), 
     +        LKH_NM(MAX,MAX+1), LKH_DN(MAX,MAX+1), CHECK 
      CHARACTER*8, NAME(MAX), FIL, BASE(MAX), LOCUS
      CHARACTER*1, PAWS
      INTEGER I,J,K,L,A

*********************************************************************
* Gets the data, i.e. allele names and frequencies present in the 
* evidence sample and the databases to be used to calculate
* likelihoods
*********************************************************************

      CALL BANNER()

      ANS = 0

      PRINT *
 601  CONTINUE
      PRINT *
      PRINT *, 'Enter the number of alleles in the mixture ',
     +         '(maximum of 12):'
      READ *, ALLELE

      IF ((ALLELE .GT. MAX) .OR. (ALLELE .LT. 1)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 601
      END IF

      PRINT*
      PRINT *, 'Enter the names of the alleles ',
     +         '(one name per line):'

      DO 11 I=1, ALLELE
         READ 100, NAME(I)
 11   CONTINUE
 100  FORMAT(A)

 643  CONTINUE
      PRINT *
      PRINT *, 'Enter the number of databases to be used ',
     +         '(maximum of 12):'
      READ *, NBAS
      IF ((NBAS .GT. MAX) .OR. (NBAS .LT. 1)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 643
      END IF

      DO 250 J = 1, NBAS
         PRINT *
         PRINT *
         PRINT '(1X,2A,I2,A)', 'Enter the name for database',
     +         ' number', J, ':'
         READ '(A)', BASE(J)
      PRINT *
      DO 10 I= 1 , ALLELE
         PRINT *, 'Enter the frequency of ',NAME(I)
         READ *, FREQ(I, J)
 10   CONTINUE

      DO 61 L=1, ALLELE
      FQ(L, J)  = FREQ(L, J)
 61   CONTINUE

      CHECK = 0.0
      DO 120 K=1, ALLELE
         CHECK = CHECK + FREQ(K, J)
 120  CONTINUE
      IF ( NINT(1000*CHECK) .GE. 1005) THEN
         PRINT *
         PRINT *, 'WARNING: The sum of the allele ',
     +            'frequencies may be greater than 1.00!'
      END IF
 250  CONTINUE

      PRINT *

 602  CONTINUE
      PRINT *
      PRINT *, 'Do you want to use the "p^2" formulas or the',
     +         ' "2p" modification?'
      PRINT *, '("p^2" = 0, "2p" = 1)'
      READ *, FORM
      IF ((FORM .NE. 0) .AND. (FORM .NE. 1)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 602
      END IF

 650  CONTINUE
      PRINT *
      PRINT *, 'Do you want a summary file? ( No=0, Yes=1 ):'
      READ *, SMRY
      IF ((SMRY .NE. 0) .AND. (SMRY .NE. 1)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 650
      END IF
      IF (SMRY .EQ. 1) THEN
         PRINT *
         PRINT *, 'Enter the name of the file:'
         READ '(A8)', FIL
         OPEN (UNIT=3, FILE=FIL, STATUS='UNKNOWN')
         PRINT *
         PRINT *, 'Enter a title for the summary ',
     +            '(such as a locus name):'
         READ '(A)', LOCUS
      END IF


*********************************************************************
* Gets the alleles whose sources are unknown and the possible number 
* of unknown contributors for the numaerator of the likelihood ratio
*********************************************************************


 999  CONTINUE
      PRINT *
      PRINT *
      PRINT *, '**********************************************',
     +         '**************************'
      PRINT *, '    THE FOLLOWING IS INPUT FOR THE NUMERATOR',
     +         ' OF THE LIKELIHOOD RATIO  '
      PRINT *, '**********************************************',
     +         '**************************'

 603  CONTINUE
      PRINT *
      PRINT *, 'Enter the number of alleles whose sources ',
     +         'are unknown:'
      READ *, UNK_AL

      IF ((UNK_AL .GT. ALLELE) .OR. (UNK_AL .LT. 0)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 603
      END IF

      DO 66 I = 1, ALLELE
         US(I) = I
 66   CONTINUE

      IF ((UNK_AL .EQ. 0) .OR. (UNK_AL .EQ. ALLELE)) GO TO 990

      PRINT *
      PRINT *, 'Enter the number(s) corresponding to the ',
     +         'unknown allele(s):'
      
      DO 12 I= 1, ALLELE
         PRINT 101, I, ' = ', NAME(I)
 12   CONTINUE
 101  FORMAT (1X, I2, 2A)
      PRINT *, '(one entry per line):'

      DO 13 I= 1, UNK_AL
 604     CONTINUE
         READ *, US(I)
         IF ((US(I) .GT. ALLELE) .OR. (US(I) .LT. 1)) THEN
            PRINT *, 'NOT A VALID ENTRY!'
            GO TO 604
         END IF
         DO 69 J = 1, I-1
            IF (US(I) .EQ. US(J)) THEN
            PRINT *, 'THIS ALLELE HAS ALREADY BEEN ENTERED!'
            GO TO 604
            END IF
 69      CONTINUE
 13   CONTINUE

      DO 251 I=1, NBAS
      DO 62 L=1, ALLELE
      FREQ(L, I)  = FQ(L, I)
 62   CONTINUE
 251  CONTINUE

      CALL SORT(FREQ, US, UNK_AL, NBAS)

 990  CONTINUE

 605  CONTINUE
      PRINT *
      PRINT *, 'Enter the lower bound on the number of ',
     +         'unknown contributors:'
      PRINT *, '(maximum of 12)'
      READ *, CONT_L

      IF ((CONT_L.LT.0).OR.(UNK_AL.GT.CONT_L*2).OR.
     +    (CONT_L.GT.MAX)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 605
      END IF

      DO 378 I=1, UNK_AL
         US1(I) = US(I)
 378  CONTINUE


 705  CONTINUE
      PRINT *
      PRINT *, 'Enter the upper bound on the number of ',
     +         'unknown contributors:'
      PRINT *, '(maximum of 12)'
      READ *, CONT_U

      IF ((CONT_U .LT. CONT_L) .OR. (CONT_U.GT.MAX)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 705
      END IF


      DO 262 J=CONT_L , CONT_U 

      PRINT *
      PRINT '(1X, A, I2)', 'Unknown Contributors: ', J
      PRINT *, 'Database       Numerator Probability'
      PRINT *, '------------------------------------'

      DO 252 I=1, NBAS
     
      IF (FORM .EQ. 0) THEN
      CALL CALC (FREQ, ALLELE, UNK_AL, J, LKH_NM(I,J+1), I)
      ELSE IF (FORM .EQ. 1) THEN
      CALL CALC2P (FREQ, ALLELE, UNK_AL, J, LKH_NM(I,J+1), I)
      END IF

      IF (NINT(1000*LKH_NM(I,J+1)) .GT. 1000) THEN
         LKH_NM(I,J+1) = 1.00000
      END IF

      WRITE(*,800) BASE(I), LKH_NM(I,J+1)
 800  FORMAT (1X, A8, ' ', E22.6)

 252  CONTINUE
      PRINT *
 262  CONTINUE


*********************************************************************
* Gets the alleles whose sources are unknown and the possible number 
* of unknown contributors for the denominator of the likelihood ratio
*********************************************************************


      PRINT *
      PRINT *, '**********************************************',
     +         '**************************'
      PRINT *, '   THE FOLLOWING IS INPUT FOR THE DENOMINATOR ',
     +         'OF THE LIKELIHOOD RATIO   '
      PRINT *, '**********************************************',
     +         '**************************'

 606  CONTINUE
      PRINT *
      PRINT *, 'Enter the number of alleles whose sources are ',
     +         'unknown:'
      READ *, UNK_AL2

      IF ((UNK_AL2 .GT. ALLELE) .OR. (UNK_AL2 .LT. 0)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 606
      END IF

      DO 67 I=1, ALLELE
         US(I) = I
 67   CONTINUE

      IF ((UNK_AL2 .EQ. 0) .OR. (UNK_AL2 .EQ. ALLELE)) GO TO 991

      PRINT *
      PRINT *, 'Enter the number(s) corresponding to the ',
     +         'unknown allele(s):'
      
      DO 14 I= 1, ALLELE
         PRINT 101, I, ' = ', NAME(I)
 14   CONTINUE
      PRINT *, '(one entry per line):'

      DO 15 I= 1, UNK_AL2
 607     CONTINUE
         READ *, US(I)
         IF ((US(I) .GT. ALLELE) .OR. (US(I) .LT. 1)) THEN
            PRINT *, 'NOT A VALID ENTRY!'
            GO TO 607
         END IF
         DO 71 J = 1, I-1
            IF (US(I) .EQ. US(J)) THEN
            PRINT *, 'THIS ALLELE HAS ALREADY BEEN ENTERED!'
            GO TO 607
            END IF
 71      CONTINUE
 15   CONTINUE

      DO 253 I=1, NBAS
      DO 63 L=1, ALLELE
      FREQ(L, I)  = FQ(L, I)
 63   CONTINUE
 253  CONTINUE

      CALL SORT(FREQ, US, UNK_AL2, NBAS)

 991  CONTINUE

 608  CONTINUE

      PRINT *
      PRINT *, 'Enter the lower bound on the number of ',
     +         'unknown contributors:'
      PRINT *, '(maximum of 12)'
      READ *, CONT_L2

      IF ((CONT_L2 .LT. 0) .OR. (UNK_AL2 .GT. CONT_L2*2)
     +    .OR.(CONT_L2 .GT. MAX)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 608
      END IF

 706  CONTINUE
      PRINT *
      PRINT *, 'Enter the upper bound on the number of ',
     +         'unknown contributors:'
      PRINT *, '(maximum of 12)'
      READ *, CONT_U2

      IF ((CONT_U2 .LT. CONT_L2).OR.(CONT_U2.GT.MAX)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 706
      END IF

      DO 263 J=CONT_L2 , CONT_U2 

      PRINT *
      PRINT '(1X, A, I2)', 'Unknown Contributors: ', J
      PRINT *, 'Database     Denominator Probability'
      PRINT *, '------------------------------------'

      DO 254 I=1, NBAS

      IF (FORM .EQ. 0) THEN
      CALL CALC (FREQ, ALLELE, UNK_AL2, J, LKH_DN(I,J+1), I)
      ELSE IF (FORM .EQ. 1) THEN
      CALL CALC2P (FREQ, ALLELE, UNK_AL2, J, LKH_DN(I,J+1), I)
      END IF

      IF (NINT(1000*LKH_DN(I,J+1)) .GT. 1000) THEN
         LKH_DN(I,J+1) = 1.00000
      END IF

      WRITE(*,800) BASE(I), LKH_DN(I,J+1)
 
 254  CONTINUE
      PRINT *
 263  CONTINUE


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
      PRINT *, '(Press return after each database)'
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
         PRINT '(1X,A8,I12,I21,F31.2)', BASE(I), J, K, 
     +         LKH_NM(I,J+1)/LKH_DN(I,K+1)
 753  CONTINUE
 752  CONTINUE
      READ '(A)', PAWS
 751  CONTINUE

*********************************************************************
* Writes a summary of all the input and output performed so far
*********************************************************************

      IF (SMRY .EQ. 1) THEN
      IF (ANS .EQ. 0) THEN
      WRITE (3,*)
      WRITE (3,*) 'Title: ', LOCUS
      IF (FORM .EQ. 0) WRITE (3, '(1X,A)') 'Form : p^2'
      IF (FORM .EQ. 1) WRITE (3, '(1X,A)') 'Form : 2p'
      DO 876 L=1, NBAS, 3
      WRITE (3,*)
      WRITE (3,*) ('Database: ', BASE(J),'      ', J=L,MIN(L+2,NBAS)) 
      WRITE (3,*) ('Allele        Freq      ', I=L,MIN(L+2,NBAS))
      WRITE (3,*) ('------------------      ', K=L,MIN(L+2,NBAS))
      DO 34 A = 1, ALLELE
        WRITE(3,521) (NAME(A), ' ', FQ(A,I),'      ',I=L,MIN(L+2,NBAS))
 521  FORMAT (1X, 12(A8, A, F9.6,A))
 34   CONTINUE
      WRITE (3,*) ('------------------      ', K=L,MIN(L+2,NBAS))
 876  CONTINUE
      WRITE (3,'(//1X, A,A/)') '*************************************',
     +       '***************************************'
      END IF
      WRITE (3, *) 'Unknown Alleles in the Numerator:'
      DO 35 A = 1, UNK_AL
         WRITE (3, '(1X, A8)') NAME(US1(A))
 35   CONTINUE
      IF (UNK_AL .EQ. 0) WRITE(3,*) 'None'
      WRITE (3,*)
      DO 662 J=CONT_L , CONT_U 
      WRITE(3,*)
      WRITE(3,'(1X, A, I2)') 'Unknown Contributors: ', J
      WRITE(3,*)'Database       Numerator Probability'
      WRITE(3,*)'------------------------------------'
      DO 652 I=1, NBAS
      WRITE (3,800) BASE(I), LKH_NM(I,J+1)
 652  CONTINUE
      WRITE (3,*)
 662  CONTINUE

      WRITE (3, '(/1X, A)') 'Unknown Alleles in the Denominator:'
      DO 36 A = 1, UNK_AL2
         WRITE (3, '(1X, A8)') NAME(US(A))
 36   CONTINUE
      IF (UNK_AL2 .EQ. 0) WRITE(3,*) 'None'
      WRITE (3,*)
      DO 663 J=CONT_L2 , CONT_U2 
      WRITE(3,*)
      WRITE(3,'(1X, A, I2)') 'Unknown Contributors: ', J
      WRITE(3,*) 'Database     Denominator Probability'
      WRITE(3,*) '------------------------------------'
      DO 654 I=1, NBAS
      WRITE(3,800) BASE(I), LKH_DN(I,J+1)
 654  CONTINUE
      WRITE(3,*)
 663  CONTINUE

      WRITE(3,*)
      WRITE(3,*)
      WRITE(3,*)'           Numerator Unknown    Denominator Unknown'
      WRITE(3,*)'Database      Contributors         Contributors    ',
     +         '     Likelihood Ratio'
      WRITE(3,*)'---------------------------------------------------',
     +         '---------------------'
      DO 851 I = 1, NBAS         
      DO 852 J = CONT_L, CONT_U
      DO 853 K = CONT_L2, CONT_U2
         WRITE(3,'(1X,A8,I12,I21,F31.2)') BASE(I), J, K, 
     +         LKH_NM(I,J+1)/LKH_DN(I,K+1)
 853  CONTINUE
 852  CONTINUE
      WRITE(3,*)
 851  CONTINUE

      WRITE (3,'(/1X, A,A/)') '*************************************',
     +      '***************************************'
      END IF

*********************************************************************
* Prompts the user if he/she wants to calculate another likelihood
* ratio. Returns to the numerator input if 'yes'.
*********************************************************************

 651  CONTINUE
      PRINT *
      PRINT *, 'Do you want to compute another likelihood',
     +         ' ratio? ( No=0, Yes=1 ):'
      READ *, ANS
      IF ((ANS .NE. 0) .AND. (ANS .NE. 1)) THEN
         PRINT *, 'NOT A VALID ENTRY!'
         GO TO 651
      END IF
      IF (ANS .EQ. 1) GO TO 999

      PRINT *
      WRITE(*,431)'Program written by John Storey on July 7, 1997.'
      WRITE(*,431)' Questions and/or comments should be sent to:  '
      WRITE(*,431)'          storey@statgen.ncsu.edu              ' 
 431  FORMAT (12(' '), A)
      PRINT *

      END

*********************************************************************
* Displays a banner giving a reference to the paper from which the 
* formulae come
*********************************************************************

      SUBROUTINE BANNER()

      PRINT *
      PRINT *, '******************************************************'
      PRINT *, '* This program performs calculations for the methods *'
      PRINT *, '* and formulas presented in:                         *'
      PRINT *, '*                                                    *'   
      PRINT *, '* Weir BS, Triggs CM, Starling L, Stowell LI,        *'
      PRINT *, '* Walsh KAJ, Buckleton J. Interpreting DNA           *'
      PRINT *, '* Mixtures. J Forensic Sci 1997; 42(2):213-222.      *'
      PRINT *, '******************************************************'

      END

*********************************************************************
* Places the alleles whose sources are unknown in the front of the 
* array of allele frequencies
*********************************************************************


      SUBROUTINE SORT(FREQ, US, UNK_AL, NBAS)

      INTEGER MAX, MAX2
      PARAMETER (MAX=12, MAX2=924)

      INTEGER UNK_AL, US(MAX), Q, NBAS
      DOUBLE PRECISION FREQ(MAX, MAX), TMP
      INTEGER I, J, K

      DO 51 I = 1, UNK_AL
         DO 52 J = 1, UNK_AL-1
            IF(US(J) .GT. US(J+1)) THEN
               TMP = US(J)
               US(J) = US(J+1)
               US(J+1) = TMP
            END IF
 52      CONTINUE
 51   CONTINUE

      DO 55 I=1, NBAS
      DO 53 J=1, UNK_AL
         Q = US(J)
         TMP = FREQ(Q, I)
         DO 54 K = US(J)-1, J, -1
            FREQ(K+1, I) = FREQ(K, I)
 54      CONTINUE
         FREQ(J, I) = TMP
 53   CONTINUE
 55   CONTINUE

      END

*********************************************************************
* Performs the calculations for the 'p^2' method
*********************************************************************

      SUBROUTINE CALC(F, ALLELE, UNK_AL, X, SUM, NBAS)

      INTEGER MAX, MAX2
      PARAMETER (MAX=12, MAX2=924)

      DOUBLE PRECISION F(MAX, MAX), EXT, SUM, PAR
      INTEGER ALLELE, UNK_AL, X, CM(MAX2,MAX), C, N, NBAS
      INTEGER I,J,K,L,M

      EXT = 0.0
      SUM = 0.0
      PAR = 0.0

      DO 10 I = UNK_AL+1, ALLELE
         EXT = EXT + F(I, NBAS)
 10   CONTINUE

      DO 20 J = 1, ALLELE
         PAR = PAR + F(J, NBAS)
 20   CONTINUE

      SUM = SUM + (PAR)**(2*X)

      N = 1
      DO 30 K = UNK_AL-1, 1, -1
         CALL COMB(UNK_AL, K, CM, C)
         DO 40 L = 1, C-1
            PAR = 0.0
            DO 50 M = 1, K
               PAR = PAR + F(CM(L,M), NBAS)
 50         CONTINUE
            SUM = SUM + ((-1)**N)*((PAR+EXT)**(2*X))
 40      CONTINUE
         N = N+1
 30   CONTINUE

      IF (UNK_AL .NE. 0) THEN
      SUM = SUM + ((-1)**N)*(EXT**(2*X))
      END IF

      END

*********************************************************************
* Computes all possible combinations of unknown alleles to be
* included in the various partial sums
*********************************************************************

      SUBROUTINE COMB(N, K, CM, C)

      INTEGER MAX, MAX2
      PARAMETER (MAX=12, MAX2=924)

      INTEGER N, K, CM(MAX2,MAX), C, ARR(MAX)
      INTEGER M, J, L

      C = 1
      
      DO 10 M=1, MAX
         ARR(M) = M
 10   CONTINUE
      
 99   CONTINUE

      DO 20 J=1, K
         CM(C,J) = ARR(J)
 20   CONTINUE

      C = C + 1

      DO 30 J = K, 1, -1
         ARR(J) = ARR(J)+1
         IF (ARR(J) .LE. J+N-K) GO TO 101
 30   CONTINUE

 101  CONTINUE

      IF (ARR(1) .GT. (N-K+1)) GO TO 999

      DO 40 L=J, K
         ARR(L+1) = ARR(L)+1
 40   CONTINUE
      
      GO TO 99

 999  CONTINUE
      
      END

*********************************************************************
* Performs the calculations for the '2p' method
*********************************************************************
         
      SUBROUTINE CALC2P(F, ALLELE, UNK_AL, X, SUM, NBAS)

      INTEGER MAX, MAX2
      PARAMETER (MAX=12, MAX2=924)

      DOUBLE PRECISION F(MAX, MAX), PAR, SUM, G(MAX)
      INTEGER ALLELE, UNK_AL, X, CM(MAX2, MAX), C, N, NBAS
      INTEGER K, L, M, P, Q, R, S, T, U, V

      PAR = 0.0
      SUM = 0.0
      N = 0

      DO 10 K = UNK_AL, 1, -1
         CALL COMB(UNK_AL, K, CM, C)
         DO 20 L = 1, C-1
            PAR = 0.0
            DO 30 M = 1, K
               G(M) = F(CM(L, M), NBAS)
 30         CONTINUE
            DO 40 P = K+1, K+ALLELE-UNK_AL
               G(P) = F(P+UNK_AL-K, NBAS)
 40         CONTINUE
            DO 50 Q=1, K+ALLELE-UNK_AL
               PAR = PAR + G(Q)
 50         CONTINUE
            DO 60 R=1, K+ALLELE-UNK_AL-1
               DO 70 S = R+1, K+ALLELE-UNK_AL
                  PAR = PAR + G(R)*G(S)
 70            CONTINUE
 60         CONTINUE
            SUM = SUM + ((-1)**N)*(PAR**X)
 20      CONTINUE
         N = N + 1
 10   CONTINUE

      PAR = 0.0
      DO 80 T = UNK_AL+1, ALLELE
         PAR = PAR + F(T,NBAS)
 80   CONTINUE

      DO 90 U = UNK_AL+1, ALLELE-1
         DO 100 V = U+1, ALLELE
            PAR = PAR + F(U,NBAS)*F(V,NBAS)
 100     CONTINUE
 90   CONTINUE

      SUM = SUM + ((-1)**N)*(PAR**X)

      SUM = SUM * (2**X)

      END

*********************************************************************
* End of program
*********************************************************************
 

