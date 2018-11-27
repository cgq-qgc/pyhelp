FORTRAN code very close to the original source code of version 3.07.

The code can be found here:
https://github.com/cgq-qgc/pyhelp/commit/a9da11c8107275a4c2b64b857f3301c744aaf64e

Only 3 modifications were made on the original source code of HELP 3.07:

1: The saturated conductivity (RC) of the soil was defined as a double-precision.
https://github.com/cgq-qgc/pyhelp/commit/a9da11c8107275a4c2b64b857f3301c744aaf64e

2: Remove reference to UTLTY for compilation
https://github.com/cgq-qgc/pyhelp/commit/ac585c1694447e1b8c0fbb67ae269324470ce96a

3: Fix GETDAT and GETTIM compatibility issue 
https://github.com/cgq-qgc/pyhelp/commit/f1b2bf20e746ba8b6f0b9622ca96ed40f785a641