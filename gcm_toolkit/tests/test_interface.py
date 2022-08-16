"""
General GCMT tests
"""
import os

import numpy as np
from gcm_toolkit import GCMT
from gcm_toolkit.tests.test_gcmtools_common import (
    petitradtrans_testdata,
    all_raw_testdata,
)
from petitRADTRANS import Radtrans


def test_prt_interface(petitradtrans_testdata, all_raw_testdata):
    dirname, expected = all_raw_testdata
    data_path = expected.get("rel_data_dir", "{}").format(dirname)

    dirname_prt, expected_prt = petitradtrans_testdata
    os.environ["pRT_input_data_path"] = dirname_prt

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

    phases = np.linspace(0, 1, 50)
    interface = tools.get_prt_interface(pRT)
    interface.set_data(time=expected["times"][-1])
    interface.chem_from_poorman(
        "T", co_ratio=0.55, feh_ratio=0.0
    )  # the temperature key sets the gcm key used for the temperature

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

    assert (ph_with_filename == ph_without_filename).all()
    assert np.isclose(ph_with_filename.max(), expected["prt_max"])
    assert np.isclose(ph_with_filename.min(), expected["prt_min"])
    assert np.isclose(ph_with_filename.mean(), expected["prt_mean"])
