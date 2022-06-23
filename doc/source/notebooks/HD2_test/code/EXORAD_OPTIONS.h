C  CPP options file for EXORAD package

#ifndef EXORAD_OPTIONS_H
#define EXORAD_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_EXORAD
C     Package-specific Options & Macros go here

C  allow for full radiative transfer (EXORAD_FULL = .TRUE. is still needed)
#define ALLOW_EXORAD_FULL

C  allow debugging output during runtime to stdout
C  Note: This may slow down the execution significantly!
C #define ALLOW_EXORAD_PRINT
#undef ALLOW_EXORAD_PRINT

C  allow for clouds (by Sven Kiefer)
#undef ALLOW_EXORAD_CLOUDS

C allow for NG acceleration
C Note: This may result in an unefficient code
#undef ALLOW_EXORAD_NG


#endif /* ALLOW_EXORAD */
#endif /* EXORAD_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
