C     Created 2022-04-02
C     This file is created automatically. This file sets the metadata input for the opacities used in exorad.
C     The content of this file may be known at compile time and may not change thereafter!
C     However, you can still change your opacities, if this does not require a change of this file.
C     Move this file to the code directory!
C     Copyright Aaron Schneider, 2021
C
C
      CHARACTER(len=19) EXORAD_OPACITY_ABS_FILE
      CHARACTER(len=20) EXORAD_OPACITY_SCAT_FILE
      CHARACTER(len=18) EXORAD_OPACITY_GWEIGHTS_FILE
      CHARACTER(len=15) EXORAD_OPACITY_PGRID_FILE
      CHARACTER(len=15) EXORAD_OPACITY_TGRID_FILE
      CHARACTER(len=15) EXORAD_OPACITY_FREQEDGES_FILE
      CHARACTER(len=12) EXORAD_OPACITY_I0_FILE
      CHARACTER(len=19) EXORAD_OPACITY_MUWEIGHTS_FILE
      CHARACTER(len=12) EXORAD_OPACITY_MU_FILE
      INTEGER EXORAD_OPACITY_NF
      INTEGER EXORAD_OPACITY_NP
      INTEGER EXORAD_OPACITY_NT
      INTEGER EXORAD_OPACITY_NG
      INTEGER EXORAD_OPACITY_NMU
      PARAMETER(
     &EXORAD_OPACITY_ABS_FILE = "opac_kappa_abs_R_S0",
     &EXORAD_OPACITY_SCAT_FILE = "opac_kappa_scat_R_S0",
     &EXORAD_OPACITY_GWEIGHTS_FILE = "opac_gweights_R_S0",
     &EXORAD_OPACITY_PGRID_FILE = "opac_pgrid_R_S0",
     &EXORAD_OPACITY_TGRID_FILE = "opac_tgrid_R_S0",
     &EXORAD_OPACITY_FREQEDGES_FILE = "opac_fgrid_R_S0",
     &EXORAD_OPACITY_I0_FILE = "opac_I0_R_S0",
     &EXORAD_OPACITY_MUWEIGHTS_FILE = "opac_muweights_R_S0",
     &EXORAD_OPACITY_MU_FILE = "opac_mu_R_S0",
     &EXORAD_OPACITY_NF = 7,
     &EXORAD_OPACITY_NP = 48,
     &EXORAD_OPACITY_NT = 1000,
     &EXORAD_OPACITY_NG = 16,
     &EXORAD_OPACITY_NMU = 3
     &)
CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***