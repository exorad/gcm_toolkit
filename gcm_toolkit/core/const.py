"""
==============================================================
                      gcm_toolkit const
==============================================================
Consistent namings of variables

These variable names are currently synced with the naming scheme
from the cubedsphere package.
The advantage of that is that we do not need to worry about
variable naming inconsistancies with exorad data.
Warning: If at some points other standards are used here,
    we would need to make sure to rename variables
    in the readin routines of exorad data!
==============================================================
"""
import cubedsphere.const as cc

VARNAMES = dict(
    FACEDIM=cc.FACEDIM,  # index of the facedimension
    iter="iter",  # iter index
    time=cc.time,
    j=cc.j,  # Y index
    i=cc.i,  # X index
    i_g=cc.i_g,  # X index at interface
    j_g=cc.j_g,  # Y index at interface
    k=cc.k,  # Z index
    k_l=cc.k_l,  # upper Z interface
    k_p1=cc.k_p1,  # outer Z interface
    k_u=cc.k_u,  # lower Z interface
    Z=cc.Z,  # Z index
    Z_l=cc.Z_l,  # lower Z interface
    Z_p1=cc.Z_p1,  # outer Z interface
    Z_u=cc.Z_u,  # upper Z interface
    Z_geo=cc.Z_geo,  # geometrical height
    T=cc.T,  # Temperature
    U="U",
    V="V",
    W=cc.W,
    Ttave=cc.Ttave,  # Temperature averaged
    wVeltave=cc.wVeltave,  # vertical velocity timeaveraged
    drW=cc.drW,
    drS=cc.drS,
    HFacW=cc.HFacW,
    HFacS=cc.HFacS,
    HFacC=cc.HFacC,
    drF=cc.drF,
    drC=cc.drC,
    dxC=cc.dxC,
    dxG=cc.dxG,
    dyC=cc.dyC,
    dyG=cc.dyG,
    rA=cc.rA,
    rAz=cc.rAz,
    rAs=cc.rAs,
    rAw=cc.rAw,
    lon=cc.lon,
    lon_b=cc.lon_b,
    lat_b=cc.lat_b,
    lat=cc.lat,
    AngleCS=cc.AngleCS,
    AngleSN=cc.AngleSN,
    dxF=cc.dxF,
    dyU=cc.dyU,
    dxV=cc.dxV,
    dyF=cc.dyF,
    P_rot="P_rot",
    P_orb="P_orb",
    g="g",
    R_p="R_p",
)

SUPPORTED_GCMS = ["MITgcm"]
