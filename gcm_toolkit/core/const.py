"""
==============================================================
                      gcm_toolkit const
==============================================================
Consistent namings of variables

These variable names are currently in sync with the naming scheme
from the cubedsphere package.
The advantage of that is that we do not need to worry about
variable naming inconsistancies with mitgcm data.
==============================================================
"""
VARNAMES = dict(
    FACEDIM="face",  # index of the facedimension
    iter="iter",  # iter index
    time="time",
    j="j",  # Y index
    i="i",  # X index
    i_g="i_g",  # X index at interface
    j_g="j_g",  # Y index at interface
    k="k",  # Z index
    k_l="k_l",  # upper Z interface
    k_p1="k_p1",  # outer Z interface
    k_u="k_u",  # lower Z interface
    Z="Z",  # Z index
    Z_l="Z_l",  # lower Z interface
    Z_p1="Z_p1",  # outer Z interface
    Z_u="Z_u",  # upper Z interface
    Z_geo="Z_geo",  # geometrical height
    T="T",  # Temperature
    U="U",
    V="V",
    W="W",
    Ttave="Ttave",  # Temperature averaged
    wVeltave="wVeltave",  # vertical velocity timeaveraged
    drW="drW",
    drS="drS",
    HFacW="HFacW",
    HFacS="HFacS",
    HFacC="HFacC",
    drF="drF",
    drC="drC",
    dxC="dxC",
    dxG="dxG",
    dyC="dyC",
    dyG="dyG",
    rA="rA",
    rAz="rAz",
    rAs="rAs",
    rAw="rAw",
    lon="lon",
    lon_b="lon_b",
    lat_b="lat_b",
    lat="lat",
    AngleCS="AngleCS",
    AngleSN="AngleSN",
    dxF="dxF",
    dyU="dyU",
    dxV="dxV",
    dyF="dyF",
    P_rot="P_rot",
    P_orb="P_orb",
    g="g",
    R_p="R_p",
    Kappa="Kappa",
    cp="cp",
    R_s="R_s",
    p_ref="p_ref",
    p0="p0",
    dt="dt",
)

SUPPORTED_GCMS = ["MITgcm"]
