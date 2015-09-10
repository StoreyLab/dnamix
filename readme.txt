7/15/97

DNAMIX version 1 is a FORTRAN 77 program that will perform calculations for the
methods and formulae presented in:

Weir BS, Triggs CM, Starling L, Stowell LI, Walsh KAJ, Buckleton J. Interpreting
DNA Mixtures. J Forensic Sci 1997; 42(2):213-222.
 
The file dnamix.f contains the source code which should compile under a FORTRAN 
77 or a FORTRAN 90 compiler. There is also a version written in C which is 
described below. This program was written by John Storey under the guidance of 
Professor Bruce Weir, Program in Statistical Genetics, Statistics Department, 
North Carolina State University.

The prompts and input should be fairly straightforward. There are many checks 
that should prevent illegal or nonsensical input by the user. If the allele 
frequencies of a database sum to a number greater than one a warning is 
displayed, but this is allowed. This may be done on purpose if one is doing 
calculations with a possible unseen allele with a frequency of one to be 
conservative. It is the case, though, that if a calculated probability ends up 
being greater than one (because of a situation as just described) it will be 
changed back to one for obvious reasons.

Calculations for up to 12 alleles can be performed with as many as 12 unknown 
contributors and 12 databases of the allele frequencies. If you need to do 
calculations for more than 12 alleles, then contact me and I can change the 
parameters very easily. This program allows calculations for both the 'p^2' 
method and the '2p' modification described in the paper. The double precision 
data type is used for all calculations, so hopefully rounding errors will not be
too much of a factor even though some of the probabilities may be very small. 
All input for character strings may be up to eight characters long and can be 
skipped by pressing return. This will result in a character string of all 
blanks.

There is an option to create a summary file that will record all input and 
output performed. I recommend using this to avoid cutting and pasting and having
to organize your results from the output to the screen.

Questions and/or comments should be sent to:

storey@statgen.ncsu.edu   or
weir@stat.ncsu.edu

For an example of the program see example.txt.

DISCLAIMER: I assume no responsibility for mistakes in the program or caused by
the program. Use DNAMIX at your own risk.

John Storey
storey@statgen.ncsu.edu

-----------------

8/4/97

I have just added a C version of the program. The file dnamix.c should compile
under a C compiler. Do not forget to include "-lm" when you compile since we
included the math library. The output and input are very similar to the
FORTRAN version, so read the rest of this document for information about the
program. 

John Storey
storey@statgen.ncsu.edu
