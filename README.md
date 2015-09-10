# Introduction

DNAMIX is a computer program to calcluate likelihood ratios as they pertain to mixed DNA samples encountered in forensic science. The latest version incorporates population structure into the calculations (see DNAMIX v. 2 below). If you have any problems downloading the program or if you have questions about the program,

DISCLAIMER: I assume no responsibility for mistakes in the program or caused by the program. Use DNAMIX at your own risk.

## DNAMIX v. 2

DNAMIX v. 2 is based on formulas and methods derived in the paper:

- Curran JM, Triggs CM, Buckleton J, Weir BS. 1999. Interpreting DNA mixtures in structured populations. Journal of Forensic Sciences, 44: 937-995.

DNAMIX v. 2 is written in FORTRAN 90; if you have a FORTRAN 90 compiler, then you may want to download the source file `DNAMIX2/dnamix2.f` and compile it yourself. The executable file is available for Unix and Windows. A user manual is also available.

### Windows

- `DNAMIX2/dnamix2_windows.f` (source code file)
- `DNAMIX2/dnamix2.exe` (executable file)

### Linux

- `DNAMIX2/dnamix2_unix.f` (source code file)

### Manual
- `DNAMIX2/manual.pdf`

## DNAMIX v. 1
DNAMIX v. 1 is written in both FORTRAN 77 and C. If you have a compiler for either language, then you may want to download the source code and compile it yourself. Currently, the executable file is only available for PC's. DNAMIX v. 1 is based on formulas derived in the paper:

- Weir BS, Triggs CM, Starling L, Stowell LI, Walsh KAJ, Buckleton J. 1997. Interpreting DNA Mixtures. Journal of Forensic Sciences, 42: 213-222.

### Relevant files

- `DNAMIX/dnamix.c` (source code for C version)
- `DNAMIX/dnamix.f` (source code for FORTRAN 77 version)
- `DNAMIX/example.txt` (shows an example run of the program)
- `DNAMIX/readme.txt` (gives directions on how to use the program)
- `DNAMIX/DNAMIX.EXE` (the executable file for Windows)

Please be sure to obtain the files `DNAMIX/example.txt` and `DNAMIX/readme.txt` as they explain how to use the program.

## Related References

- Curran JM, Triggs CM, Buckleton J, Weir BS. 1999. Interpreting DNA mixtures in structured populations. Journal of Forensic Sciences, 44: 937-995.
- Evett IW, Weir BS. 1998. Interpreting DNA evidence: Statistical genetics for forensic science. Sunderland, MA: Sinauer.
- Weir BS. 1998. The coancestry coefficient in forensic science. Proc 8th Int Symp Hum Identification. Madison, WI: Promega.
- Weir BS, Triggs CM, Starling L, Stowell LI, Walsh KAJ, Buckleton J. 1997. Interpreting DNA Mixtures. Journal of Forensic Sciences, 42: 213-222.
- National Research Council. 1996. The evaluation of forensic DNA evidence. Washington, DC: National Academy Press.
- Wright S. 1951. The genetical structure of populations. Annals of Eugenics, 15: 323-354.
