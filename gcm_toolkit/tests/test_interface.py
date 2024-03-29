"""
General GCMT tests
"""
import os

import numpy as np
import pytest

from gcm_toolkit import GCMT
from gcm_toolkit.utils.interface import Interface
from gcm_toolkit.tests.test_gcmtools_common import (
    petitradtrans_testdata,
    all_raw_testdata,
    all_nc_testdata,
)


def test_prt_interface(petitradtrans_testdata, all_raw_testdata):
    """
    Test petitRADTRANS interface.
    Note: Currently requires raw readin, because we do some regridding of the data!
    """
    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    dirname_prt, expected_prt = petitradtrans_testdata
    os.environ["pRT_input_data_path"] = dirname_prt

    from petitRADTRANS import Radtrans

    pRT = Radtrans(
        line_species=expected_prt["line_species"],
        rayleigh_species=expected_prt["rayleigh_species"],
        continuum_opacities=expected_prt["continuum_opacities"],
        wlen_bords_micron=expected_prt["wlen_bords_micron"],
        do_scat_emis=True,
    )
    # Note: the pRT.setup_opa_structure is done by tools internally

    tools = GCMT(p_unit="bar", time_unit="day")  # create a GCMT object
    tools.read_raw(
        gcm=expected["gcm"], data_path=data_path, d_lat=15, d_lon=15
    )

    # add cloud properties
    dsi_clouds = tools.get_models("0")
    tmp = np.ones(
        (
            len(dsi_clouds["time"]),
            len(dsi_clouds["lat"]),
            len(dsi_clouds["lon"]),
            len(dsi_clouds["Z"]),
        )
    )
    dsi_clouds["ClAb"] = (("time", "lat", "lon", "Z"), tmp * 1e-5)
    dsi_clouds["ClDs"] = (("time", "lat", "lon", "Z"), tmp * 1e2)
    dsi_clouds["ClDr"] = (("time", "lat", "lon", "Z"), tmp * 1e3)
    dsi_clouds["ClVf_Test1[s]"] = (("time", "lat", "lon", "Z"), tmp * 0.5)
    dsi_clouds["ClVf_Test2[s]"] = (("time", "lat", "lon", "Z"), tmp * 0.5)
    tools._replace_model("0", dsi_clouds)

    phases = np.linspace(0, 1, 50)
    interface = tools.get_prt_interface(pRT)

    interface.set_data(time=expected["times"][-1])
    interface.chem_from_poorman("T", co_ratio=0.55, feh_ratio=0.0)

    spectra = interface.calc_phase_spectrum(
        mmw=expected["MMW"],
        Rstar=expected["Rstar"],
        Tstar=expected["Tstar"],
        semimajoraxis=expected["semimajoraxis"],
        normalize=True,
    )

    ph_without_filename = interface.phase_curve(phases, spectra=spectra)

    filename = "spectrum_test.nc"
    interface.calc_phase_spectrum(
        mmw=expected["MMW"],
        Rstar=expected["Rstar"],
        Tstar=expected["Tstar"],
        semimajoraxis=expected["semimajoraxis"],
        filename=filename,
        normalize=True,
    )

    ph_with_filename = interface.phase_curve(phases, filename=filename)

    with pytest.raises(ValueError):
        interface.phase_curve(
            phases
        )  # We need a spectrum or a filename with a spectrum

    assert (ph_with_filename == ph_without_filename).all()
    assert np.isclose(ph_with_filename.max(), expected["prt_max"])
    assert np.isclose(ph_with_filename.min(), expected["prt_min"])
    assert np.isclose(ph_with_filename.sum(), expected["prt_sum"])

    R_p = interface.dsi.attrs["R_p"]
    del interface.dsi.attrs["R_p"]
    with pytest.raises(ValueError):
        interface.calc_phase_spectrum(
            mmw=expected["MMW"],
            Rstar=expected["Rstar"],
            Tstar=expected["Tstar"],
            semimajoraxis=expected["semimajoraxis"],
            normalize=True,
        )
    interface.dsi.attrs["R_p"] = R_p

    with pytest.raises(ValueError):
        interface.calc_phase_spectrum(
            mmw=expected["MMW"],
            Rstar=expected["Rstar"],
            Tstar=None,
            semimajoraxis=expected["semimajoraxis"],
            normalize=True,
        )

    # test transit calculation set up
    interface.set_data(
        time=expected["times"][-1], terminator_avg=True, lon_resolution=60
    )
    assert (interface.dsi["lon"].values == [-90, 90]).all()
    assert interface.dsi["lat"].values == [0]

    # test transit calculation
    interface.chem_from_poorman("T", co_ratio=0.55, feh_ratio=0.0)
    wave, spectra = interface.calc_transit_spectrum(mmw=2.33)
    assert sum(spectra) == 45237213620.18512

    # test transit calculation
    interface.chem_from_poorman("T", co_ratio=0.55, feh_ratio=0.0)
    wave, spectra = interface.calc_transit_spectrum(mmw=2.33, clouds=True)
    assert sum(spectra) == 45331292591.13728

    # test transit calculation
    interface.chem_from_poorman("T", co_ratio=0.55, feh_ratio=0.0)
    wave, spectra = interface.calc_transit_spectrum(
        mmw=2.33, clouds=True, use_bruggemann=True
    )
    assert sum(spectra) == 45237213620.53514

    # Test if Pa works
    interface.dsi.attrs["p_unit"] = "Pa"
    interface.set_data(time=expected["times"][-1])

    # Test if something weird works
    interface.dsi.attrs["p_unit"] = "wrong"
    with pytest.raises(NotImplementedError):
        interface.set_data(time=expected["times"][-1])

    os.remove(filename)


def test_raise_error_on_large_data(all_nc_testdata, petitradtrans_testdata):
    dirname, expected = all_nc_testdata

    dirname_prt, expected_prt = petitradtrans_testdata
    os.environ["pRT_input_data_path"] = dirname_prt

    from petitRADTRANS import Radtrans

    pRT = Radtrans(
        line_species=expected_prt["line_species"],
        rayleigh_species=expected_prt["rayleigh_species"],
        continuum_opacities=expected_prt["continuum_opacities"],
        wlen_bords_micron=expected_prt["wlen_bords_micron"],
        do_scat_emis=True,
    )

    tools = GCMT(p_unit="bar", time_unit="day")  # create a GCMT object
    tools.read_reduced(data_path=dirname)
    interface = tools.get_prt_interface(pRT)

    # Now continue without regrid lowres
    interface.set_data(time=expected["times"][-1])
    interface.chem_from_poorman("T", co_ratio=0.55, feh_ratio=0.0)

    with pytest.raises(ValueError):
        # Resolution is too high, we should raise an error!
        interface.calc_phase_spectrum(
            mmw=expected["MMW"],
            Rstar=expected["Rstar"],
            Tstar=expected["Tstar"],
            semimajoraxis=expected["semimajoraxis"],
            normalize=True,
        )


def test_general_interface(all_nc_testdata, petitradtrans_testdata):
    """General test for interface class"""
    dirname, expected = all_nc_testdata

    dirname_prt, expected_prt = petitradtrans_testdata
    os.environ["pRT_input_data_path"] = dirname_prt

    tools = GCMT(p_unit="bar", time_unit="day")  # create a GCMT object
    tools.read_reduced(data_path=dirname)

    interface = Interface(tools)

    with pytest.raises(ValueError):
        # We need to set the data first
        interface.chemistry.chem_from_poorman(
            "T", co_ratio=0.55, feh_ratio=0.0
        )

    with pytest.raises(ValueError):
        # We need to set the data first
        interface.chem_from_poorman("T", co_ratio=0.55, feh_ratio=0.0)

    # Test that regrid_lowres works as well
    interface.set_data(time=expected["times"][-1], regrid_lowres=True)
    assert len(interface.dsi.lon) * len(interface.dsi.lat) == 288

    interface.chemistry.dsi = tools.get("0")
    assert "time" in interface.chemistry.dsi.dims
    with pytest.raises(ValueError):
        # We should not have time in the dimensions
        interface.chem_from_poorman("T", co_ratio=0.55, feh_ratio=0.0)

    # Testing units
    for unit in ["Pa", "wrong"]:
        dsi = tools.get("0").sel(time=expected["times"][-1])
        dsi.attrs["p_unit"] = unit
        interface.chemistry.dsi = dsi
        if unit == "Pa":
            interface.chem_from_poorman("T", co_ratio=0.55, feh_ratio=0.0)
        else:
            with pytest.raises(NotImplementedError):
                # other units should not work
                interface.chem_from_poorman("T", co_ratio=0.55, feh_ratio=0.0)

    for unit in ["Pa", "wrong"]:
        interface.chemistry.abunds.attrs.update({"p_unit": unit})
        if unit == "Pa":
            interface.chemistry.to_prt(["H2O"], np.array(0.5))
        else:
            with pytest.raises(NotImplementedError):
                interface.chemistry.to_prt(["H2O"], np.array(0.5))

    with pytest.raises(ValueError):
        interface.chemistry.to_prt(["wrong_species"], np.array(0.5))


def test_PAC_interface(all_nc_testdata):
    """General test for PAC interface class."""
    dirname, expected = all_nc_testdata

    tools = GCMT(p_unit="bar", time_unit="day")  # create a GCMT object
    tools.read_reduced(data_path=dirname)

    from gcm_toolkit.utils.interface import PACInterface

    pac1d = PACInterface(tools, 1)
    pac2d = PACInterface(tools, 2)
    with pytest.raises(ValueError):
        PACInterface(tools, 3)

    pac1d.set_data(time=expected["times"][-1])
    pac2d.set_data(time=expected["times"][-1])

    with pytest.raises(OSError):
        pac1d.write_inputfile("non/existent/path")
        pac2d.write_inputfile("non/existent/path")
        pac2d.generate_lptfile("non/existent/path")
        pac2d.generate_aptfiles("non/existent/path")
